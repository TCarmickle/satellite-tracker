"""RINEX 3.04 observation and navigation files (the exchange format every GNSS processor reads: RTKLIB's
rnx2rtkp, NRCan CSRS-PPP, OPUS...), and a reader for both so what is written can be checked by reading it back.

Only what this receiver produces is covered: one frequency (L1 / E1), observation types C, L, D, S per signal;
GPS (1019) and Galileo I/NAV (1046) broadcast orbits. Layout: RINEX 3.04 (IGS/RTCM, 2018), sections 5-6 and
tables A2-A4, A6-A8, A10. Times are GPS time: RTKLIB-style GPST seconds (counted like Unix time, no leap seconds),
so time.gmtime() of one gives its GPST calendar date.
"""
import calendar
import math
import os
import time

import rtcm_msm

VERSION = 3.04
LEAP_S = 18                       # GPST - UTC since 2017 (no leap second announced after it)


def _hl(text, label):
    """one header line: 60 columns of content, then the label"""
    return '%-60.60s%-20s\n' % (text, label)


def _cal(t):
    """GPST seconds -> (year, month, day, hour, minute, seconds as float)"""
    whole = math.floor(t + 1e-9)
    g = time.gmtime(whole)
    return g.tm_year, g.tm_mon, g.tm_mday, g.tm_hour, g.tm_min, g.tm_sec + (t - whole if t - whole > 1e-9 else 0.0)


def obs_types(epochs):
    """{system: ['C1C', 'L1C', 'D1C', 'S1C', ...]} for the signals present, in a fixed order"""
    codes = {}
    for _, obs in epochs:
        for sat, o in obs.items():
            codes.setdefault(sat[0], set()).add(o['code'])
    return {s: [k + c for c in sorted(cs) for k in 'CLDS'] for s, cs in sorted(codes.items())}


def _ssi(cn0):
    """RINEX signal strength indicator 1-9 from C/N0 (dB-Hz): 1 < 12, 2 12-17, ..., 9 >= 54"""
    return min(9, max(1, int(cn0 // 6)))


def write_obs(path, epochs, header):
    """epochs: [(t_gpst, {sat: {'code', 'P', 'L', 'D', 'cn0', 'lli'}})] (sat 'G05', 'E14', 'S31'); header: dict
    with 'marker', 'receiver' (number, type, version), 'antenna' (number, type), 'xyz' (approximate ECEF metres),
    optional 'delta_hen' (antenna height, east, north offsets m), 'program', 'run_by', 'comments' [str], 'interval'.
    Written to path + '.tmp' first, then moved over (a cut-off run never leaves half a file)."""
    types = obs_types(epochs)
    systems = ''.join(types)
    sysch = systems if len(systems) == 1 else 'M'
    first, last = epochs[0][0], epochs[-1][0]
    if 'interval' in header:
        interval = header['interval']
    else:
        gaps = sorted(b[0] - a[0] for a, b in zip(epochs, epochs[1:]) if b[0] > a[0])
        interval = gaps[len(gaps) // 2] if gaps else 1.0
    lines = [_hl('%9.2f%11s%-20s%-20s' % (VERSION, '', 'OBSERVATION DATA', sysch), 'RINEX VERSION / TYPE'),
             _hl('%-20.20s%-20.20s%-20.20s' % (header.get('program', 'trk2rinex'), header.get('run_by', ''),
                                                time.strftime('%Y%m%d %H%M%S UTC', time.gmtime())),
                 'PGM / RUN BY / DATE')]
    for c in header.get('comments', []):
        lines.append(_hl(c, 'COMMENT'))
    rec, ant = header.get('receiver', ('', '', '')), header.get('antenna', ('', ''))
    x, y, z = header.get('xyz', (0.0, 0.0, 0.0))
    h, e, n = header.get('delta_hen', (0.0, 0.0, 0.0))
    lines += [_hl(header.get('marker', 'UNKN'), 'MARKER NAME'),
              _hl('NON_GEODETIC', 'MARKER TYPE'),
              _hl('%-20.20s%-40.40s' % (header.get('observer', ''), header.get('agency', '')), 'OBSERVER / AGENCY'),
              _hl('%-20.20s%-20.20s%-20.20s' % tuple(rec), 'REC # / TYPE / VERS'),
              _hl('%-20.20s%-20.20s' % tuple(ant), 'ANT # / TYPE'),
              _hl('%14.4f%14.4f%14.4f' % (x, y, z), 'APPROX POSITION XYZ'),
              _hl('%14.4f%14.4f%14.4f' % (h, e, n), 'ANTENNA: DELTA H/E/N')]
    for s, ts in types.items():
        for i in range(0, len(ts), 13):
            head = '%-1s  %3d' % (s, len(ts)) if i == 0 else '      '
            lines.append(_hl(head + ''.join(' %-3s' % k for k in ts[i:i + 13]), 'SYS / # / OBS TYPES'))
    lines += [_hl('DBHZ', 'SIGNAL STRENGTH UNIT'),
              _hl('%10.3f' % interval, 'INTERVAL'),
              _hl('%6d%6d%6d%6d%6d%13.7f%5s%-3s' % (_cal(first) + ('', 'GPS')), 'TIME OF FIRST OBS'),
              _hl('%6d%6d%6d%6d%6d%13.7f%5s%-3s' % (_cal(last) + ('', 'GPS')), 'TIME OF LAST OBS'),
              _hl('%6d' % 0, 'RCV CLOCK OFFS APPL')]
    for s, ts in types.items():
        for k in ts:
            if k[0] == 'L':
                lines.append(_hl('%-1s %-3s %8.5f' % (s, k, 0.0), 'SYS / PHASE SHIFT'))
    lines += [_hl('%3d' % 0, 'GLONASS SLOT / FRQ #'),
              _hl(' C1C    0.000 C1P    0.000 C2C    0.000 C2P    0.000', 'GLONASS COD/PHS/BIS'),
              _hl('%6d' % LEAP_S, 'LEAP SECONDS'),
              _hl('', 'END OF HEADER')]
    for t, obs in epochs:
        if not obs:
            continue
        lines.append('> %4d %02d %02d %02d %02d%11.7f  %1d%3d\n' % (_cal(t) + (0, len(obs))))
        for sat in sorted(obs):
            o = obs[sat]
            row = sat
            for k in types[sat[0]]:
                v = None
                if k[1:] == o['code']:
                    v = {'C': o['P'], 'L': o['L'], 'D': o['D'], 'S': o['cn0']}[k[0]]
                if v is None:
                    row += ' ' * 16
                elif k[0] == 'L':
                    lli = o['lli'] & 3
                    row += '%14.3f%1s%1d' % (v, lli if lli else ' ', _ssi(o['cn0']))
                else:
                    row += '%14.3f  ' % v
            lines.append(row.rstrip() + '\n')
    _atomic(path, ''.join(lines))


def _d19(v):
    return '%19.12E' % v


def write_nav(path, ephs, header):
    """ephs: rtcm_msm.decode_ephemeris dicts with 'week' (full GPS week) and 'ttr' (GPST); one record per
    satellite and toe (the first copy received). header: 'program', 'run_by', 'comments'."""
    seen, recs = set(), []
    for e in sorted(ephs, key=lambda e: (e['sat'], e['week'], e['toe'], e['ttr'])):
        key = (e['sat'], e['week'], e['toe'], e['iode'])
        if key in seen:
            continue
        seen.add(key)
        recs.append(e)
    systems = ''.join(sorted({e['sys'] for e in recs})) or 'G'
    sysch = systems if len(systems) == 1 else 'M'
    lines = [_hl('%9.2f%11s%-20s%-20s' % (VERSION, '', 'N: GNSS NAV DATA', sysch), 'RINEX VERSION / TYPE'),
             _hl('%-20.20s%-20.20s%-20.20s' % (header.get('program', 'trk2rinex'), header.get('run_by', ''),
                                                time.strftime('%Y%m%d %H%M%S UTC', time.gmtime())),
                 'PGM / RUN BY / DATE')]
    for c in header.get('comments', []):
        lines.append(_hl(c, 'COMMENT'))
    lines += [_hl('%6d' % LEAP_S, 'LEAP SECONDS'), _hl('', 'END OF HEADER')]
    gps_epoch = 315964800
    for e in recs:
        t_toc = gps_epoch + e['week'] * 604800 + e['toc']
        t_toe = gps_epoch + e['week'] * 604800 + e['toe']
        if t_toc - t_toe > 302400:                  # toc in the week before / after the toe's
            t_toc -= 604800
        elif t_toe - t_toc > 302400:
            t_toc += 604800
        y, mo, d, h, mi, s = _cal(t_toc)
        # transmission time: seconds of the toe's week (may fall outside 0..604800 when it was sent the week before)
        ttr_sow = e['ttr'] - (gps_epoch + e['week'] * 604800)
        lines.append('%s %4d %02d %02d %02d %02d %02d%s%s%s\n' % (e['sat'], y, mo, d, h, mi, int(round(s)),
                                                                 _d19(e['af0']), _d19(e['af1']), _d19(e['af2'])))
        if e['sys'] == 'G':
            rows = [(e['iode'], e['crs'], e['deln'], e['M0']), (e['cuc'], e['e'], e['cus'], e['sqrtA']),
                    (e['toe'], e['cic'], e['OMG0'], e['cis']), (e['i0'], e['crc'], e['omg'], e['OMGd']),
                    (e['idot'], e['code'], e['week'], e['l2p']),
                    (rtcm_msm.ura_m(e['ura']), e['svh'], e['tgd'], e['iodc']), (ttr_sow, e['fit_h'], None, None)]
        else:
            # data sources 517 = I/NAV E1-B (bit 0) + E5b (bit 2) + clock for the E5b,E1 pair (bit 9)
            rows = [(e['iode'], e['crs'], e['deln'], e['M0']), (e['cuc'], e['e'], e['cus'], e['sqrtA']),
                    (e['toe'], e['cic'], e['OMG0'], e['cis']), (e['i0'], e['crc'], e['omg'], e['OMGd']),
                    (e['idot'], 517, e['week'], None),
                    (rtcm_msm.sisa_m(e['sisa']), e['svh'], e['bgd_a'], e['bgd_b']), (ttr_sow, None, None, None)]
        for r in rows:
            lines.append('    ' + ''.join(' ' * 19 if v is None else _d19(float(v)) for v in r).rstrip() + '\n')
    _atomic(path, ''.join(lines))


def _atomic(path, text):
    tmp = path + '.tmp'
    with open(tmp, 'w', newline='\n') as fh:
        fh.write(text)
    os.replace(tmp, path)


# -- reading back (for checks; also reads other receivers' single-frequency RINEX 3)
def _float(s):
    s = s.strip().replace('D', 'E').replace('d', 'e')
    return float(s) if s else None


def read_obs(path):
    """a RINEX 3 observation file -> (header {label: [content lines]}, types {sys: [codes]},
    [(t_gpst, {sat: {type: (value, lli, ssi)}})])"""
    header, types, epochs = {}, {}, []
    with open(path) as fh:
        last_sys = None
        for line in fh:
            label = line[60:].strip()
            header.setdefault(label, []).append(line[:60])
            if label == 'SYS / # / OBS TYPES':
                if line[0] != ' ':
                    last_sys = line[0]
                    types[last_sys] = []
                types[last_sys] += line[7:60].split()
            if label == 'END OF HEADER':
                break
        cur = None
        for line in fh:
            if line.startswith('>'):
                y, mo, d, h, mi = (int(line[2:6]), int(line[7:9]), int(line[10:12]), int(line[13:15]),
                                   int(line[16:18]))
                t = calendar.timegm((y, mo, d, h, mi, 0)) + float(line[18:29])
                cur = {}
                epochs.append((t, cur))
                continue
            if cur is None or len(line) < 3:
                continue
            sat = line[:3]
            vals = {}
            for i, k in enumerate(types.get(sat[0], [])):
                field = line[3 + 16 * i:3 + 16 * (i + 1)]
                v = _float(field[:14]) if len(field) >= 14 else _float(field)
                if v is None:
                    continue
                lli = field[14:15].strip()
                ssi = field[15:16].strip()
                vals[k] = (v, int(lli) if lli else 0, int(ssi) if ssi else 0)
            cur[sat] = vals
    return header, types, epochs


def read_nav(path):
    """a RINEX 3 navigation file -> [{'sat', 'toc' (GPST s), 'values': [af0, af1, af2, then 4 per line]}]"""
    out = []
    with open(path) as fh:
        for line in fh:
            if line[60:].strip() == 'END OF HEADER':
                break
        cur = None
        for line in fh:
            if line[0] != ' ':
                y, mo, d, h, mi, s = (int(line[4:8]), int(line[9:11]), int(line[12:14]), int(line[15:17]),
                                      int(line[18:20]), int(line[21:23]))
                cur = {'sat': line[:3], 'toc': calendar.timegm((y, mo, d, h, mi, s)),
                       'values': [_float(line[23 + 19 * i:42 + 19 * i]) for i in range(3)]}
                out.append(cur)
            elif cur is not None:
                cur['values'] += [_float(line[4 + 19 * i:23 + 19 * i]) for i in range(4)]
    return out

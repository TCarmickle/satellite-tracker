r"""pocket_trk's raw measurements of a GNSS capture -> RINEX 3.04 observation + navigation files, for working out
the antenna's position precisely offline (RTKLIB rnx2rtkp against a CORS station, or NRCan CSRS-PPP; README.txt
"PRECISE ANTENNA POSITION" has the recipe). Nothing is sent anywhere: the files are only written.

Usage:
  python trk2rinex.py gnss-20261106-205411            every run of that capture (data\captures\<name>.trkN.*)
  python trk2rinex.py C:\path\gnss-...trk1.rtcm         one run (its .log beside it is used if there)
options:
  --out DIR          where to write (default data\rinex\<capture>\)
  --source auto|log|rtcm   auto: the log's $OBS lines when the log is there, else the RTCM MSM7 (same numbers;
                     the log's carrier needs no re-anchoring repair)
  --galileo          include Galileo (E1-B) observations and orbits (few satellites decode an orbit here; leave
                     it out for a GPS-only CORS base)
  --no-sbas          leave the WAAS satellites out (CORS bases don't record them; RTKLIB ignores them without orbits)
  --marker NAME      the 4-character station name in the file names (default SDR1)
Only what pocket_trk measured goes in: GPS L1 C/A (+ WAAS L1, + Galileo E1-B), code, carrier, Doppler, C/N0.
"""
import argparse
import glob
import math
import os
import re
import statistics
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.join(os.path.dirname(HERE), 'app')
sys.path.insert(0, APP)
import rinex      # noqa: E402
import rtcm_msm   # noqa: E402

# the data folder, found the way paths.py finds it (paths itself is not imported: it moves files on import)
DATA = os.environ.get('IRIDIUM_TRACKER_DATA') or os.path.join(os.path.dirname(APP), 'data')
GPS_EPOCH = 315964800


def find_runs(arg):
    """a capture name or one run's file -> (capture name, [(run stem, run number)])"""
    m = re.match(r'(.*)\.trk(\d+)(\.rtcm|\.log)?$', arg)
    if m:
        stem = arg[:-len(m.group(3))] if m.group(3) else arg
        return os.path.basename(m.group(1)), [(stem, int(m.group(2)))]
    name = os.path.basename(arg)
    folder = os.path.dirname(arg) or os.path.join(DATA, 'captures')
    runs = set()
    for f in glob.glob(os.path.join(folder, glob.escape(name) + '.trk*.*')):
        r = re.search(r'\.trk(\d+)\.(rtcm|log)$', f)
        if r:
            runs.add((f[:-len(r.group(2)) - 1], int(r.group(1))))
    return name, sorted(runs, key=lambda r: r[1])


def capture_time(name):
    """the start time in a capture name (local time; only used to fix the GPS week, so hours off do not matter)"""
    m = re.search(r'(\d{8})-(\d{6})', name)
    if not m:
        return None
    return time.mktime(time.strptime(m.group(1) + m.group(2), '%Y%m%d%H%M%S')) + 18


def read_run(stem, source, name):
    """one pocket_trk run -> (ephemerides, info: 'raw' epochs, their code-carrier 'rate', 'source', 'check'...)"""
    rtcm_path, log_path = stem + '.rtcm', stem + '.log'
    have_rtcm, have_log = os.path.exists(rtcm_path), os.path.exists(log_path)
    info = {'run': os.path.basename(stem), 'pos': [], 'version': None}
    log_epochs, ref = [], None
    if have_log:
        log_epochs, ref, info['pos'] = rtcm_msm.read_log(log_path)
        with open(log_path, errors='replace') as fh:
            v = re.search(r'VER=(\S+)', fh.readline())
            info['version'] = v.group(1) if v else None
    if ref is None:
        ref = capture_time(name) or (os.path.getmtime(rtcm_path if have_rtcm else log_path) + 18)
    msm, ephs = rtcm_msm.read_rtcm(open(rtcm_path, 'rb').read(), ref) if have_rtcm else ([], [])
    use = source if source != 'auto' else ('log' if have_log else 'rtcm')
    if use == 'log' and not have_log or use == 'rtcm' and not have_rtcm:
        raise SystemExit('%s: no %s file' % (stem, use))
    info['source'] = use
    if have_log and have_rtcm:
        info['check'] = cross_check(msm, log_epochs)
    info['raw'] = log_epochs if use == 'log' else msm
    info['rate'] = rtcm_msm.estimate_rate(info['raw'])
    return ephs, info


def cross_check(msm, log_epochs):
    """the MSM7 numbers against the log's $OBS for the same epochs: largest differences (carrier: apart from the
    whole cycles the encoder re-anchors by)"""
    log = {round(t, 3): o for t, o in log_epochs}
    worst = {'P': 0.0, 'L': 0.0, 'D': 0.0, 'cn0': 0.0}
    n = both = 0
    for t, obs in msm:
        lo = log.get(round(t, 3))
        if lo is None:
            continue
        both += 1
        for s, o in obs.items():
            l = lo.get(s)
            if not l or None in (o['P'], o['L'], o['D']):
                continue
            n += 1
            d = o['L'] - l['L']
            worst['P'] = max(worst['P'], abs(o['P'] - l['P']))
            worst['L'] = max(worst['L'], abs(d - round(d)))
            worst['D'] = max(worst['D'], abs(o['D'] - l['D']))
            worst['cn0'] = max(worst['cn0'], abs(o['cn0'] - l['cn0']))
    return {'msm_epochs': len(msm), 'log_epochs': len(log_epochs), 'common_epochs': both, 'values': n,
            'worst': worst}


def convert(infos):
    """each run's raw epochs -> RINEX-ready epochs (rtcm_msm.consistent); a run too short to measure its
    code-carrier rate borrows the median of the others (the same tuner, so the same frequency error)"""
    known = [i['rate'] for i in infos if i['rate'] is not None]
    out = []
    for i in infos:
        rate = i['rate'] if i['rate'] is not None else (statistics.median(known) if known else 0.0)
        epochs, i['stats'] = rtcm_msm.consistent(i.pop('raw'), rate)
        out.append(epochs)
    return out


def gnss_eph(e):
    """a decode_ephemeris dict -> the form gnss.py computes positions from"""
    import gnss
    toe = GPS_EPOCH + e['week'] * 604800 + e['toe']
    toc = GPS_EPOCH + e['week'] * 604800 + e['toc']
    toc += -604800 if toc - toe > 302400 else 604800 if toe - toc > 302400 else 0
    d = {k: e[k] for k in ('iode', 'svh', 'e', 'i0', 'OMG0', 'omg', 'M0', 'deln', 'OMGd', 'idot', 'crc', 'crs',
                           'cuc', 'cus', 'cic', 'cis', 'af0', 'af1', 'af2')}
    d.update({'iodc': e.get('iodc', e['iode']), 'sva': e.get('ura', e.get('sisa', 0)), 'toe_gpst': int(toe),
              'toc_gpst': int(toc), 'ttr_gpst': int(e['ttr']), 'A': e['sqrtA'] ** 2, 'toes': e['toe'],
              'fit': e.get('fit_h', 0.0), 'tgd': e.get('tgd', e.get('bgd_b', 0.0)), 'code': e.get('code', 0),
              'flag': e.get('l2p', 0)})
    return gnss.finish_eph(d, e['sys'], e['prn'])


def spp_fixes(epochs, ephs, every=10):
    """plain single-point fixes (ECEF m, receiver clock m, residual rms m) from the GPS pseudoranges of every
    `every`-th epoch: the approximate position for the header, and a check that the ranges are sound"""
    import numpy as np
    import gnss
    by_sat = {}
    for e in ephs:
        if e['sys'] == 'G' and e['svh'] == 0:
            by_sat.setdefault(e['sat'], []).append(gnss_eph(e))
    omge, out, x = 7.2921151467e-5, [], np.zeros(4)
    for t, obs in epochs[::every]:
        rows = []
        for s, o in obs.items():
            cand = by_sat.get(s)
            if not cand:
                continue
            e = min(cand, key=lambda e: abs(e['toe_gpst'] - t))
            if abs(e['toe_gpst'] - t) > 4 * 3600:
                continue
            ttx = t - o['P'] / rtcm_msm.C
            for _ in range(2):
                dts = float(np.ravel(gnss.sat_clock(e, gnss.gpst_to_unix(ttx)))[0]) - e['tgd']
                ttx = t - o['P'] / rtcm_msm.C - dts
            rows.append((np.array(gnss.sat_ecef(e, gnss.gpst_to_unix(ttx))[0]) * 1000.0, dts, o['P']))
        if len(rows) < 5:                    # 4 would fit exactly and say nothing about the ranges
            continue
        xi = x.copy() if out else np.zeros(4)
        for _ in range(10):
            H, v = [], []
            for pos, dts, P in rows:
                th = omge * np.linalg.norm(pos - xi[:3]) / rtcm_msm.C       # Earth turns while the signal flies
                ps = np.array([pos[0] * math.cos(th) + pos[1] * math.sin(th),
                               -pos[0] * math.sin(th) + pos[1] * math.cos(th), pos[2]])
                r = np.linalg.norm(ps - xi[:3])
                H.append(list((xi[:3] - ps) / r) + [1.0])
                v.append(P - (r + xi[3] - rtcm_msm.C * dts))
            dx = np.linalg.lstsq(np.array(H), np.array(v), rcond=None)[0]
            xi += dx
            if np.linalg.norm(dx) < 1e-3:
                break
        res = np.array(v) - np.array(H) @ dx
        rms = float(np.sqrt(np.mean(res ** 2)))
        if rms < 30 and 6.3e6 < np.linalg.norm(xi[:3]) < 6.4e6:
            x = xi
            out.append((t, xi[:3].copy(), float(xi[3]), rms))
    return out


def ecef(lat, lon, h):
    a, f = 6378137.0, 1 / 298.257223563
    e2 = f * (2 - f)
    la, lo = math.radians(lat), math.radians(lon)
    n = a / math.sqrt(1 - e2 * math.sin(la) ** 2)
    return ((n + h) * math.cos(la) * math.cos(lo), (n + h) * math.cos(la) * math.sin(lo),
            (n * (1 - e2) + h) * math.sin(la))


def merge(runs):
    """several runs' epochs -> one time-ordered list (an epoch in two runs keeps the one with more satellites)"""
    by_t = {}
    for epochs in runs:
        for t, obs in epochs:
            k = round(t, 3)
            if k not in by_t or len(obs) > len(by_t[k][1]):
                by_t[k] = (t, obs)
    return [by_t[k] for k in sorted(by_t)]


def select(epochs, keep):
    """only the systems in `keep` ('GS', 'GSE'...) -> (epochs, number left out for having no GPS satellite).
    An epoch without a GPS satellite (pocket_trk's own fix lost, or not yet made) leaves nothing to tie the
    receiver clock to: its WAAS / Galileo ranges are of no use to a GPS solver, so it is left out - and every
    satellite starts a new arc after the hole (a loss of lock flagged in a left-out second must not get lost)."""
    kept, no_gps, hole = [], 0, False
    for t, obs in epochs:
        obs = {s: o for s, o in obs.items() if s[0] in keep}
        if not any(s[0] == 'G' for s in obs):
            no_gps += 1 if obs else 0
            hole = True
            continue
        if hole:
            for o in obs.values():
                o['lli'] |= 1
            hole = False
        kept.append((t, obs))
    return kept, no_gps


def rinex_name(marker, t0, t1, kind):
    """RINEX 3 long file name: SSSSMRCCC_R_YYYYDDDHHMM_PPU_FFU_XX.rnx (country unknown: XXX)"""
    g = time.gmtime(t0)
    minutes = max(1, int(math.ceil((t1 - t0 + 1) / 60.0)))
    if minutes < 60:
        dur = '%02dM' % minutes
    elif minutes < 24 * 60:
        dur = '%02dH' % int(math.ceil(minutes / 60.0))
    else:
        dur = '%02dD' % int(math.ceil(minutes / 1440.0))
    return '%-4.4s00XXX_R_%04d%03d%02d%02d_%s_01S_%s.rnx' % (marker.upper().ljust(4, '0'), g.tm_year, g.tm_yday,
                                                             g.tm_hour, g.tm_min, dur, kind)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('capture')
    ap.add_argument('--out')
    ap.add_argument('--source', choices=('auto', 'log', 'rtcm'), default='auto')
    ap.add_argument('--galileo', action='store_true')
    ap.add_argument('--no-sbas', action='store_true')
    ap.add_argument('--marker', default='SDR1')
    a = ap.parse_args(argv)
    name, runs = find_runs(a.capture)
    if not runs:
        raise SystemExit('no pocket_trk files (.trkN.rtcm / .trkN.log) found for %s' % a.capture)
    keep = 'G' + ('' if a.no_sbas else 'S') + ('E' if a.galileo else '')
    ephs, infos = [], []
    for stem, _ in runs:
        e, info = read_run(stem, a.source, name)
        ephs += e
        infos.append(info)
    all_epochs = convert(infos)
    epochs, no_gps = select(merge(all_epochs), keep)
    ephs = [e for e in ephs if e['sys'] in keep]
    if not epochs:
        raise SystemExit('%s: no usable observations' % name)

    # approximate position: our own single-point fixes from these very ranges, else pocket_trk's own fixes
    xyz, note = (0.0, 0.0, 0.0), 'none'
    try:
        fixes = spp_fixes(epochs, ephs)
    except Exception as exc:          # a broken orbit must not stop the conversion: the header value is optional
        fixes = []
        print('single-point fixes failed (%s); using pocket_trk\'s own fixes' % exc)
    pos = [p for i in infos for p in i['pos']]
    if fixes:
        xyz = tuple(statistics.median(f[1][k] for f in fixes) for k in range(3))
        note = 'median of %d own single-point fixes, residual rms median %.1f m' % (
            len(fixes), statistics.median(f[3] for f in fixes))
    elif pos:
        xyz = ecef(*(statistics.median(p[k] for p in pos) for k in range(3)))
        note = 'median of %d pocket_trk fixes' % len(pos)
    pos_check = None
    if fixes and pos:
        pk = ecef(*(statistics.median(p[k] for p in pos) for k in range(3)))
        pos_check = math.dist(xyz, pk)

    out = a.out or os.path.join(DATA, 'rinex', name)
    os.makedirs(out, exist_ok=True)
    t0, t1 = epochs[0][0], epochs[-1][0]
    rate = statistics.median(i['stats']['rate_m_s'] for i in infos)
    njumps = sum(len(i['stats']['jumps']) for i in infos)
    version = next((i['version'] for i in infos if i['version']), '?')
    sources = '/'.join(sorted({i['source'] for i in infos}))
    comments = ['pocket_trk %s raw measurements (%s), converted by' % (version, sources),
                'tools/trk2rinex.py of Satellite Tracker. Carrier phase and',
                'Doppler moved onto the code clock (code-carrier %.3f m/s,' % rate,
                '%d measurement-instant shifts). Loss of lock from carrier' % njumps,
                'continuity, not the MSM lock time. Antenna type unknown.',
                'Approx position: %s.' % note.split(',')[0]]
    head = {'marker': a.marker.upper(), 'receiver': ('', 'RTL-SDR', 'pocket_trk %s' % version),
            'antenna': ('', 'UNKNOWN         NONE'), 'xyz': xyz, 'program': 'trk2rinex', 'comments': comments}
    obs_path = os.path.join(out, rinex_name(a.marker, t0, t1, 'MO' if len({s[0] for _, o in epochs for s in o}) > 1
                                            else 'GO'))
    nav_path = os.path.join(out, rinex_name(a.marker, t0, t1, 'MN' if len({e['sys'] for e in ephs}) > 1 else 'GN'))
    rinex.write_obs(obs_path, epochs, head)
    rinex.write_nav(nav_path, ephs, {'program': 'trk2rinex', 'comments': ['broadcast orbits pocket_trk decoded']})

    print('%s: %d epochs, %.2f h, %d runs (%d epochs without GPS left out)' % (
        name, len(epochs), (t1 - t0) / 3600.0, len(runs), no_gps))
    for i in infos:
        s = i['stats']
        print('  %s from %s: %d epochs, %d observations, %d arcs, %d slips (+%d coarse), %d instant shifts %s, '
              'code-carrier %.3f m/s, %d MSM re-anchorings repaired (%d not), %d gaps'
              % (i['run'], i['source'], s['epochs'], s['obs'], s['arcs'], s['slips'], s['slips_doppler'],
                 len(s['jumps']), [round(j[1]) for j in s['jumps']], s['rate_m_s'], s['reanchored'],
                 s['reanchor_failed'], s['gaps']))
        if 'check' in i:
            c = i['check']
            print('    MSM7 vs $OBS: %d / %d epochs (%d common), largest difference P %.4f m, L %.4f cycle, '
                  'D %.4f Hz, C/N0 %.3f dB' % (c['msm_epochs'], c['log_epochs'], c['common_epochs'],
                                               c['worst']['P'], c['worst']['L'], c['worst']['D'], c['worst']['cn0']))
    print('  approximate position: %s%s' % (note, '' if pos_check is None else
                                             '; %.1f m from pocket_trk\'s own median' % pos_check))
    print('  wrote %s\n        %s' % (obs_path, nav_path))
    return obs_path, nav_path


if __name__ == '__main__':
    main()

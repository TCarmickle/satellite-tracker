"""The slow parts of the navigation messages: almanacs, UTC / leap seconds, Galileo's GPS time offset - decoded
from the raw message bits pocket_trk already logs ($NAV lines).

  $NAV,time,sat,sig,prn,nerr,size,hex       (PocketSDR sdr_nav.c out_log_nav)
    GPS L1 C/A    240 bits: one LNAV subframe, 10 words of 24 bits (parity removed)          IS-GPS-200 20.3.3
    Galileo E1B   128 bits: one I/NAV word (the even page's 112 data bits + the odd page's 16) Galileo OS SIS ICD 4.3.5
    SBAS L1       250 bits: one message, preamble first (CRC checked here)                    RTCA DO-229

An almanac is every satellite's rough orbit, sent by each satellite of the system for all of them (GPS: the whole
set every 12.5 min; Galileo: 3 satellites per minute). It is good to a few km for days - far too rough to draw a
satellite with, but enough to know which satellites are up: the 'Up but not heard' list, the search list and pass
predictions use it (labelled almanac). It is still our own reception: it came from the satellites to our antenna.

Pure functions plus one small stateful decoder (NavDecoder: Galileo's almanacs are split over several words, GPS's
almanac week comes on a separate page). Nothing here may raise on a bad line: a page that does not decode is skipped.
"""
import math

PI = 3.1415926535898            # the value the interface documents define semicircles with
GPS_EPOCH = 315964800           # 1980-01-06 00:00:00 UTC as Unix time (GPS week 0; RTKLIB-style GPST counts from it)
WEEK = 604800
GAL_SQRTA_NOM = math.sqrt(29600000.0)      # Galileo almanac semi-major axis is sent as a difference from this (m^0.5)
GAL_I_NOM = 56.0 / 180.0                    # ...and the inclination from 56 deg (semicircles)
GPS_I_NOM = 0.30                            # GPS almanac inclination: 0.30 semicircles (54 deg) + delta-i

# GALM line field order after the system tag (GPS / GAL): one almanac entry, enough for gnss._kepler
ALM_FIELDS = ['toa_gpst', 'svh', 'A', 'e', 'i0', 'OMG0', 'omg', 'M0', 'OMGd', 'af0', 'af1']
SBAS_ALM_FIELDS = ['t0', 'svh', 'x', 'y', 'z', 'vx', 'vy', 'vz']        # SBAS message type 17 (m, m/s; t0 unix)
UTC_FIELDS = ['A0', 'A1', 'tot', 'WNt', 'dtLS', 'WNlsf', 'DN', 'dtLSF']  # GPS SF4 p18 / Galileo word 6 (raw values)
GGTO_FIELDS = ['A0G', 'A1G', 't0G', 'WN0G']                             # Galileo word 10: GST - GPST


class Bits:
    """an MSB-first bit string from pocket_trk's hex dump"""

    def __init__(self, hexstr, n):
        self.n = int(n)
        nhex = len(hexstr)
        self.v = int(hexstr, 16) >> max(nhex * 4 - self.n, 0)
        if nhex * 4 < self.n:
            raise ValueError('fewer bits than the line says')

    def u(self, i, k):
        return (self.v >> (self.n - i - k)) & ((1 << k) - 1)

    def s(self, i, k):
        x = self.u(i, k)
        return x - (1 << k) if x & (1 << (k - 1)) else x


def nearest_week_time(sow, t_gpst, week_mod=None, week=None):
    """seconds of week -> the full RTKLIB-style GPST time nearest t_gpst; with a week number counted modulo
    week_mod, the instance of that week nearest t_gpst instead"""
    cur = int((t_gpst - GPS_EPOCH) // WEEK)
    if week_mod and week is not None:
        w = cur - ((cur - week) % week_mod)
        if cur - w > week_mod / 2:
            w += week_mod
        return GPS_EPOCH + w * WEEK + sow
    t = GPS_EPOCH + cur * WEEK + sow
    if t - t_gpst > WEEK / 2:
        t -= WEEK
    elif t_gpst - t > WEEK / 2:
        t += WEEK
    return t


def full_week(week, modulo, t_gpst):
    """a week number counted modulo `modulo` -> the full GPS week nearest t_gpst (RTKLIB-style GPST)"""
    cur = int((t_gpst - GPS_EPOCH) // WEEK)
    w = cur - ((cur - week) % modulo)
    if cur - w > modulo / 2:
        w += modulo
    return w


# ---------------------------------------------------------------- GPS LNAV (IS-GPS-200 20.3.3.5)
def gps_page(b):
    """(subframe id, SV id of the page or None) of a 240-bit LNAV subframe; None if it isn't one"""
    if b.n < 240 or b.u(0, 8) != 0x8B:
        return None
    sf = b.u(43, 3)
    return sf, (b.u(50, 6) if sf in (4, 5) else None)


def gps_almanac(b):
    """an almanac page (subframe 5 pages 1-24 for PRN 1-24, subframe 4 pages 2-5 / 7-10 for PRN 25-32) ->
    {prn, toas, svh, A, e, i0, OMG0, omg, M0, OMGd, af0, af1} (radians, metres, seconds), else None"""
    page = gps_page(b)
    if not page:
        return None
    sf, svid = page
    if not ((sf == 5 and 1 <= svid <= 24) or (sf == 4 and 25 <= svid <= 32)):
        return None
    w3, w4, w5, w6, w7, w8, w9, w10 = 48, 72, 96, 120, 144, 168, 192, 216
    sqrtA = b.u(w6, 24) * 2.0 ** -11
    if sqrtA < 4000:                    # a dummy page (no satellite in that slot)
        return None
    af0 = ((b.u(w10, 8) << 3) | b.u(w10 + 19, 3))            # 11 bits, split around af1
    af0 = (af0 - (1 << 11) if af0 & (1 << 10) else af0) * 2.0 ** -20
    return {'prn': svid, 'toas': b.u(w4, 8) * 2.0 ** 12, 'svh': b.u(w5 + 16, 8),
            'A': sqrtA * sqrtA, 'e': b.u(w3 + 8, 16) * 2.0 ** -21,
            'i0': (GPS_I_NOM + b.s(w4 + 8, 16) * 2.0 ** -19) * PI,
            'OMG0': b.s(w7, 24) * 2.0 ** -23 * PI, 'omg': b.s(w8, 24) * 2.0 ** -23 * PI,
            'M0': b.s(w9, 24) * 2.0 ** -23 * PI, 'OMGd': b.s(w5, 16) * 2.0 ** -38 * PI,
            'af0': af0, 'af1': b.s(w10 + 8, 11) * 2.0 ** -38}


def gps_health(b):
    """subframe 5 page 25 (SV id 51) -> {'toas', 'wna' (8 bits), 'health': {prn: 6-bit health}} for PRN 1-24;
    subframe 4 page 25 (SV id 63) -> {'health': {...}} for PRN 25-32; else None"""
    page = gps_page(b)
    if not page:
        return None
    sf, svid = page
    if sf == 5 and svid == 51:
        h = {}
        for k in range(24):
            h[k + 1] = b.u(72 + 24 * (k // 4) + 6 * (k % 4), 6)
        return {'toas': b.u(56, 8) * 2.0 ** 12, 'wna': b.u(64, 8), 'health': h}
    if sf == 4 and svid == 63:
        h = {25: b.u(168 + 18, 6)}
        for k in range(4):
            h[26 + k] = b.u(192 + 6 * k, 6)
        for k in range(3):
            h[30 + k] = b.u(216 + 6 * k, 6)
        return {'health': h}
    return None


def gps_utc(b):
    """subframe 4 page 18 (SV id 56), its second half -> the UTC parameters {A0 s, A1 s/s, tot s, WNt, dtLS s,
    WNlsf, DN, dtLSF s} (weeks counted modulo 256, as sent); the first half (Klobuchar) is spp.gps_iono's"""
    page = gps_page(b)
    if not page or page != (4, 56):
        return None
    w6, w7, w8, w9, w10 = 120, 144, 168, 192, 216
    a0 = (b.s(w7, 24) << 8) | b.u(w8, 8)                     # 32 bits: 24 MSBs in word 7, 8 LSBs in word 8
    return {'A0': a0 * 2.0 ** -30, 'A1': b.s(w6, 24) * 2.0 ** -50, 'tot': b.u(w8 + 8, 8) * 2.0 ** 12,
            'WNt': b.u(w8 + 16, 8), 'dtLS': b.s(w9, 8), 'WNlsf': b.u(w9 + 8, 8), 'DN': b.u(w9 + 16, 8),
            'dtLSF': b.s(w10, 8)}


# ---------------------------------------------------------------- Galileo I/NAV (OS SIS ICD 4.3.5)
def _gal_elements(b, i):
    """delta-sqrtA, e, omega, delta-i (the first part of an almanac entry) starting at bit i"""
    return {'dsqrtA': b.s(i, 13) * 2.0 ** -9, 'e': b.u(i + 13, 11) * 2.0 ** -16,
            'omg': b.s(i + 24, 16) * 2.0 ** -15 * PI, 'di': b.s(i + 40, 11) * 2.0 ** -14 * PI}


def gal_word(b):
    """one 128-bit I/NAV word -> (word type, its fields) for the types used here (5, 6, 7, 8, 9, 10), else
    (type, None)"""
    if b.n < 128:
        return None, None
    wt = b.u(0, 6)
    if wt == 5:         # ionosphere, both group delays, health, GST
        return wt, {'ai0': b.u(6, 11) * 2.0 ** -2, 'ai1': b.s(17, 11) * 2.0 ** -8, 'ai2': b.s(28, 14) * 2.0 ** -15,
                    'bgd_e5a': b.s(47, 10) * 2.0 ** -32, 'bgd_e5b': b.s(57, 10) * 2.0 ** -32,
                    'e5b_hs': b.u(67, 2), 'e1b_hs': b.u(69, 2), 'e5b_dvs': b.u(71, 1), 'e1b_dvs': b.u(72, 1),
                    'wn': b.u(73, 12), 'tow': b.u(85, 20)}
    if wt == 6:         # GST - UTC
        return wt, {'A0': b.s(6, 32) * 2.0 ** -30, 'A1': b.s(38, 24) * 2.0 ** -50, 'dtLS': b.s(62, 8),
                    'tot': b.u(70, 8) * 3600.0, 'WNt': b.u(78, 8), 'WNlsf': b.u(86, 8), 'DN': b.u(94, 3),
                    'dtLSF': b.s(97, 8)}
    if wt == 7:         # almanac: SVID1 (elements, its M0 included)
        d = {'ioda': b.u(6, 4), 'wna': b.u(10, 2), 't0a': b.u(12, 10) * 600.0, 'svid': b.u(22, 6)}
        d.update(_gal_elements(b, 28))
        d.update({'OMG0': b.s(79, 16) * 2.0 ** -15 * PI, 'OMGd': b.s(95, 11) * 2.0 ** -33 * PI,
                  'M0': b.s(106, 16) * 2.0 ** -15 * PI})
        return wt, d
    if wt == 8:         # SVID1's clock + health, then SVID2's elements (its M0 is in word 9)
        d = {'ioda': b.u(6, 4), 'af0': b.s(10, 16) * 2.0 ** -19, 'af1': b.s(26, 13) * 2.0 ** -38,
             'e5b_hs': b.u(39, 2), 'e1b_hs': b.u(41, 2), 'svid2': b.u(43, 6)}
        d.update({k + '2': v for k, v in _gal_elements(b, 49).items()})
        d.update({'OMG02': b.s(100, 16) * 2.0 ** -15 * PI, 'OMGd2': b.s(116, 11) * 2.0 ** -33 * PI})
        return wt, d
    if wt == 9:         # SVID2's M0, clock, health; SVID3's first elements
        d = {'ioda': b.u(6, 4), 'wna': b.u(10, 2), 't0a': b.u(12, 10) * 600.0, 'M0': b.s(22, 16) * 2.0 ** -15 * PI,
             'af0': b.s(38, 16) * 2.0 ** -19, 'af1': b.s(54, 13) * 2.0 ** -38, 'e5b_hs': b.u(67, 2),
             'e1b_hs': b.u(69, 2), 'svid3': b.u(71, 6)}
        d.update({k + '3': v for k, v in _gal_elements(b, 77).items()})
        return wt, d
    if wt == 10:        # SVID3's node, M0, clock, health; GST - GPST (GGTO)
        return wt, {'ioda': b.u(6, 4), 'OMG0': b.s(10, 16) * 2.0 ** -15 * PI, 'OMGd': b.s(26, 11) * 2.0 ** -33 * PI,
                    'M0': b.s(37, 16) * 2.0 ** -15 * PI, 'af0': b.s(53, 16) * 2.0 ** -19,
                    'af1': b.s(69, 13) * 2.0 ** -38, 'e5b_hs': b.u(82, 2), 'e1b_hs': b.u(84, 2),
                    'A0G': b.s(86, 16) * 2.0 ** -35, 'A1G': b.s(102, 12) * 2.0 ** -51, 't0G': b.u(114, 8) * 3600.0,
                    'WN0G': b.u(122, 6)}
    return wt, None


def _gal_entry(el, rest, wna, t0a, svid):
    """the pieces of one Galileo almanac entry -> the same dict gps_almanac gives (prn, toas, svh = E1-B health)"""
    sqrtA = GAL_SQRTA_NOM + el['dsqrtA']
    return {'prn': svid, 'toas': t0a, 'wna': wna, 'svh': rest['e1b_hs'], 'svh_e5b': rest['e5b_hs'],
            'A': sqrtA * sqrtA, 'e': el['e'], 'i0': GAL_I_NOM * PI + el['di'], 'OMG0': rest['OMG0'],
            'omg': el['omg'], 'M0': rest['M0'], 'OMGd': rest['OMGd'], 'af0': rest['af0'], 'af1': rest['af1']}


# ---------------------------------------------------------------- SBAS message type 17 (DO-229 A.4.4.12)
def sbas_almanac(data, t_unix):
    """message type 17's 212 data bits (list of 0/1, as gnss.sbas_message gives) -> [{prn, svh, x, y, z (m),
    vx, vy, vz (m/s), t0 (unix)}] for its up to three GEO satellites; t_unix resolves the day"""
    v = int(''.join(map(str, data[:212])), 2)
    u = lambda i, k: (v >> (212 - i - k)) & ((1 << k) - 1)
    s = lambda i, k: u(i, k) - (1 << k) if u(i, k) & (1 << (k - 1)) else u(i, k)
    t0_tod = u(201, 11) * 64.0
    out = []
    for k in range(3):
        i = 67 * k
        prn = u(i + 2, 8)
        if prn == 0:
            continue
        out.append({'prn': prn, 'svh': u(i + 10, 8),
                    'x': s(i + 18, 15) * 2600.0, 'y': s(i + 33, 15) * 2600.0, 'z': s(i + 48, 9) * 26000.0,
                    'vx': s(i + 57, 3) * 10.0, 'vy': s(i + 60, 3) * 10.0, 'vz': s(i + 63, 4) * 60.0})
    # t0 is a GPS time of day: the one nearest the moment it was received (as gnss.sbas_geo_nav does)
    day = math.floor(t_unix / 86400.0) * 86400.0
    t0 = day + t0_tod
    t0 += 86400.0 * round((t_unix - t0) / 86400.0)
    for a in out:
        a['t0'] = t0
    return out


# ---------------------------------------------------------------- formatting (the GALM line's fields)
def _g(x):
    return '%.12g' % x


def alm_text(sys_tag, a):
    """an almanac dict with 'toa_gpst' -> 'GPS,toa_gpst,svh,A,e,...' (GALM fields; sys_tag GPS / GAL / SBAS)"""
    keys = SBAS_ALM_FIELDS if sys_tag == 'SBAS' else ALM_FIELDS
    return ','.join([sys_tag] + [str(int(a[k])) if k in ('toa_gpst', 'svh') else
                                 '%.1f' % a[k] if k == 't0' else _g(a[k]) for k in keys])


def parse_alm(text):
    """GALM fields 'GPS,...' / 'GAL,...' / 'SBAS,...' -> dict with 'sys' ('G' / 'E' / 'S'), or None"""
    f = text.strip().split(',')
    keys = {'GPS': ALM_FIELDS, 'GAL': ALM_FIELDS, 'SBAS': SBAS_ALM_FIELDS}.get(f[0])
    if not keys or len(f) < 1 + len(keys):
        return None
    try:
        d = {k: float(v) for k, v in zip(keys, f[1:])}
    except ValueError:
        return None
    d['sys'] = {'GPS': 'G', 'GAL': 'E', 'SBAS': 'S'}[f[0]]
    d['svh'] = int(d['svh'])
    return d


def param_text(sys_tag, d, keys):
    return ','.join([sys_tag] + [_g(d[k]) for k in keys])


def parse_params(text, keys):
    f = text.strip().split(',')
    if len(f) < 1 + len(keys):
        return None
    try:
        d = {k: float(v) for k, v in zip(keys, f[1:])}
    except ValueError:
        return None
    d['sys'] = f[0]
    return d


# ---------------------------------------------------------------- the decoder
class NavDecoder:
    """$NAV lines in, new facts out: line(fields, t_unix, t_gpst) -> [(sat, text)] for the GALM lines -
    a new almanac entry (dedupe by satellite + its reference time), and the UTC / GGTO parameters when they change
    (sat 'UTCG' GPS, 'UTCE' Galileo, 'GGTO'). The newest of each is also kept: almanacs, utc, ggto."""
    GAL_PAIR_S = (20.0, 40.0)       # words 7+8 and 9+10 are 2 s apart, 8 then 9 28 s apart (measured on E1B)
    # Galileo sends a new almanac reference time every ~10 min (335 entries in 7 hours on the 11-06 recording): a
    # satellite's entry is passed on again only when it is this much newer - its error grows ~8 km a day, so 6 h
    # older costs ~2 km, plenty for what an almanac is used for
    GAL_STEP_S = 6 * 3600.0

    def __init__(self):
        self.almanacs = {}      # sat ('G05', 'E11', 'S31') -> (t_rx, text, toa) newest by reference time
        self.utc = {}           # 'UTCG' / 'UTCE' -> (t_rx, text)
        self.ggto = None        # (t_rx, text)
        self.gal5 = {}          # Galileo sat -> its newest word 5 (both group delays, health)
        self._gal = {}          # transmitting Galileo sat -> {word type: (t, fields)}
        self._wna = {}          # GPS almanac toa (s of week) -> WNa (8 bits), from subframe 5 page 25
        self.health = {}        # GPS prn -> 6-bit health from the page-25s

    def line(self, f, t_unix, t_gpst):
        """f: a $NAV line split on commas; t_unix / t_gpst: when it was received"""
        try:
            return self._line(f, t_unix, t_gpst)
        except (ValueError, IndexError, KeyError):
            return []

    def _line(self, f, t_unix, t_gpst):
        if len(f) < 8 or f[0] != '$NAV':
            return []
        sat, n = f[2], int(f[6])
        b = Bits(f[7], n)
        if sat[:1] == 'G' and n == 240:
            return self._gps(b, t_unix, t_gpst)
        if sat[:1] == 'E' and n == 128:
            return self._galileo(sat, b, t_unix, t_gpst)
        if sat[:1] == 'S' and n == 250:
            return self._sbas(b, t_unix)
        return []

    def _new_alm(self, sat, a, sys_tag, t_unix):
        key = a['t0'] if sys_tag == 'SBAS' else a['toa_gpst']
        old = self.almanacs.get(sat)
        if old and old[2] >= key:
            return []
        if old and sys_tag == 'GAL' and key - old[2] < self.GAL_STEP_S:
            return []
        if old and sys_tag == 'SBAS':
            # message type 17 comes every few minutes with a new time of day but (a geostationary satellite) the
            # same position to its 2.6 km steps: a new entry only once an hour, or when the position moved
            o = parse_alm(old[1])
            moved = o and math.dist((o['x'], o['y'], o['z']), (a['x'], a['y'], a['z'])) > 5000.0
            if key - old[2] < 3600.0 and not moved:
                return []
        text = alm_text(sys_tag, a)
        self.almanacs[sat] = (t_unix, text, key)
        return [(sat, text)]

    def _gps(self, b, t_unix, t_gpst):
        a = gps_almanac(b)
        if a:
            wna = self._wna.get(a['toas'])
            a['toa_gpst'] = (nearest_week_time(a['toas'], t_gpst, 256, wna) if wna is not None else
                             nearest_week_time(a['toas'], t_gpst))
            return self._new_alm('G%02d' % a['prn'], a, 'GPS', t_unix)
        h = gps_health(b)
        if h:
            if 'wna' in h:
                self._wna[h['toas']] = h['wna']
            self.health.update(h['health'])
            return []
        u = gps_utc(b)
        if u:
            return self._params('UTCG', param_text('GPS', u, UTC_FIELDS), t_unix)
        return []

    def _params(self, key, text, t_unix):
        """UTC / GGTO parameters -> a GALM fact when they are new: the leap-second fields changed, or a later
        reference time. (Each satellite sends its own upload's A0 / A1: comparing whole texts made the lines flip
        between satellites - 177 UTC lines in 7 hours on the 11-06 recording.)"""
        keys = GGTO_FIELDS if key == 'GGTO' else UTC_FIELDS
        new = parse_params(text, keys)
        old = self.ggto if key == 'GGTO' else self.utc.get(key)
        if old:
            o = parse_params(old[1], keys)
            if key == 'GGTO':
                ref, oref = (new['WN0G'], new['t0G']), (o['WN0G'], o['t0G'])
                later = ref != oref and ((new['WN0G'] - o['WN0G']) % 64 < 32)
                if not later:
                    return []
            else:
                leap = [new[k] for k in ('dtLS', 'WNlsf', 'DN', 'dtLSF')] != [o[k] for k in ('dtLS', 'WNlsf', 'DN',
                                                                                          'dtLSF')]
                ref, oref = (new['WNt'], new['tot']), (o['WNt'], o['tot'])
                later = ref != oref and (new['WNt'] - o['WNt']) % 256 < 128 and not (
                    new['WNt'] == o['WNt'] and new['tot'] < o['tot'])
                if not (leap or later):
                    return []
        if key == 'GGTO':
            self.ggto = (t_unix, text)
        else:
            self.utc[key] = (t_unix, text)
        return [(key, text)]

    def _galileo(self, sat, b, t_unix, t_gpst):
        wt, d = gal_word(b)
        if d is None:
            return []
        out = []
        if wt == 5:
            self.gal5[sat] = dict(d, t=t_unix)
            return []
        if wt == 6:
            return self._params('UTCE', param_text('GAL', d, UTC_FIELDS), t_unix)
        words = self._gal.setdefault(sat, {})
        words[wt] = (t_unix, d)
        near, far = self.GAL_PAIR_S

        def prev(k, lo, hi):
            got = words.get(k)
            return got[1] if got and lo < t_unix - got[0] <= hi and got[1]['ioda'] == d['ioda'] else None

        entry = None
        if wt == 8:
            w7 = prev(7, 0, near)
            if w7 and w7['svid']:
                entry = _gal_entry(w7, dict(w7, af0=d['af0'], af1=d['af1'], e1b_hs=d['e1b_hs'], e5b_hs=d['e5b_hs']),
                                   w7['wna'], w7['t0a'], w7['svid'])
        elif wt == 9:
            w8 = prev(8, 0, far)
            if w8 and w8['svid2']:
                el = {k: w8[k + '2'] for k in ('dsqrtA', 'e', 'omg', 'di')}
                entry = _gal_entry(el, dict(d, OMG0=w8['OMG02'], OMGd=w8['OMGd2']), d['wna'], d['t0a'], w8['svid2'])
        elif wt == 10:
            out += self._params('GGTO', param_text('GAL', d, GGTO_FIELDS), t_unix) if (d['A0G'] or d['A1G'] or
                                                                                        d['WN0G']) else []
            w9 = prev(9, 0, near)
            if w9 and w9['svid3']:
                el = {k: w9[k + '3'] for k in ('dsqrtA', 'e', 'omg', 'di')}
                entry = _gal_entry(el, d, w9['wna'], w9['t0a'], w9['svid3'])
        if entry and 1 <= entry['prn'] <= 36:
            # WNa is the GST week modulo 4 (GST week 0 = GPS week 1024, a multiple of 4: the same remainder)
            entry['toa_gpst'] = nearest_week_time(entry['toas'], t_gpst, 4, entry['wna'])
            out += self._new_alm('E%02d' % entry['prn'], entry, 'GAL', t_unix)
        return out

    def _sbas(self, b, t_unix):
        if b.u(8, 6) != 17:         # (only message type 17 is read here: skip the CRC of the other ~1 a second)
            return []
        bits = [(b.v >> (b.n - 1 - i)) & 1 for i in range(b.n)]
        import gnss
        msg = gnss.sbas_message(bits)
        if not msg or msg[0] != 17:
            return []
        out = []
        for a in sbas_almanac(msg[1], t_unix):
            if 120 <= a['prn'] <= 158:
                out += self._new_alm('S%02d' % (a['prn'] - 100), a, 'SBAS', t_unix)
        return out


def parse_utc_text(text):
    return parse_params(text, UTC_FIELDS)


def parse_ggto_text(text):
    return parse_params(text, GGTO_FIELDS)

"""GNSS satellite orbits from the broadcast ephemerides our own receiver decodes.

PocketSDR's pocket_trk logs every ephemeris it decodes as a line
  $EPH,time,sat,sig,IODE,IODC,SVA,SVH,Toe,Toc,Ttr,A,e,i0,OMEGA0,omega,M0,delta-n,OMEGAdot,Idot,
       Crc,Crs,Cuc,Cus,Cic,Cis,Toes,Fit,Af0,Af1,Af2,TGD,code,flag
(src/sdr_pvt.c out_log_eph, which prints RTKLIB's eph_t). Toe/Toc/Ttr are whole seconds of RTKLIB
gtime_t, i.e. GPS time counted like Unix time (GPST = UTC + leap seconds); BeiDou and Galileo times are
already converted to GPST there. Toes is seconds of the week in the satellite's own time scale (BDT
for BeiDou) and goes into the node angle. Angles are radians, A is metres.

The positions use the satellites' own interface documents: IS-GPS-200 (GPS, QZSS), Galileo OS SIS
ICD (same Keplerian form, Galileo GM), BeiDou ICD with CGCS2000 constants and the extra rotation of
the BeiDou GEO satellites (PRN 1-5, 59-63). Nothing here is downloaded: the satellite itself sent
these numbers to our antenna.
"""
import math

import numpy as np

import orbitfit

C_LIGHT = 299792458.0
GPS_EPOCH = 315964800          # 1980-01-06 00:00:00 UTC as Unix time (GPS week 0)
BDT_WEEK0 = 1356               # BeiDou week 0 starts at GPS week 1356 (2006-01-01), BDT = GPST - 14 s
GST_WEEK0 = 1024               # Galileo week 0 starts at GPS week 1024 (1999-08-22), GST = GPST
BDS_GEO_PRNS = set(range(1, 6)) | set(range(59, 64))

# (GM m^3/s^2, Earth rotation rad/s) per system's interface document
CONSTANTS = {'G': (3.9860050e14, 7.2921151467e-5),     # IS-GPS-200 (also QZSS)
             'J': (3.9860050e14, 7.2921151467e-5),
             'E': (3.986004418e14, 7.2921151467e-5),   # Galileo OS SIS ICD
             'C': (3.986004418e14, 7.292115e-5)}       # BeiDou ICD, CGCS2000

# GPST - UTC (s) from each Unix time on; no leap second has been announced after 2017. This table is the fallback:
# the satellites announce every leap second months ahead (GPS subframe 4 page 18, Galileo word 6 - nav_bits.py), and
# set_utc() adds a decoded one, so a future leap second is followed without a new version of this file.
LEAP_SECONDS = [(1341100800, 16), (1435708800, 17), (1483228800, 18)]
_LEAPS = list(LEAP_SECONDS)         # the table + what the satellites said (replaced whole, never changed in place)
UTC_SOURCE = ['table']              # 'table' or what set_utc was last given (for the log)
WEEK_S = 604800

FIELDS = ['iode', 'iodc', 'sva', 'svh', 'toe_gpst', 'toc_gpst', 'ttr_gpst', 'A', 'e', 'i0', 'OMG0', 'omg', 'M0',
          'deln', 'OMGd', 'idot', 'crc', 'crs', 'cuc', 'cus', 'cic', 'cis', 'toes', 'fit', 'af0', 'af1', 'af2',
          'tgd', 'code', 'flag']
INT_FIELDS = {'iode', 'iodc', 'sva', 'svh', 'toe_gpst', 'toc_gpst', 'ttr_gpst', 'code', 'flag'}

# position uncertainty model (km), see EphemerisOrbit.sigma_km
SIG_FIT_KM = 0.005        # inside the fit interval broadcast orbits are good to a few metres
SIG_HOURS = 2.0           # ...which runs this long either side of Toe
SIG_RATE_KM = 0.1         # beyond it the error grows by this per hour past it...
SIG_GROWTH_KM = 0.005     # ...plus this times (hours past it)^2, for manoeuvres and model drift


def gpst_minus_utc(t_unix):
    """GPS time minus UTC (s) at Unix time(s) t: the leap seconds since 1980"""
    t = np.asarray(t_unix, float)
    out = np.full(t.shape, 15.0)
    for start, leap in _LEAPS:
        out = np.where(t >= start, float(leap), out)
    return out


def leap_seconds(t_unix=None):
    """GPST - UTC (s, an int) at one Unix time (now by default) - the one place the receiver, the RTCM reader and
    the solver take it from"""
    import time
    return int(gpst_minus_utc(time.time() if t_unix is None else t_unix))


def set_utc(p, t_rx):
    """decoded UTC parameters (nav_bits.UTC_FIELDS: dtLS, WNlsf, DN, dtLSF as sent; GPS or Galileo - both count
    the same leap seconds) received at Unix time t_rx -> the leap-second table follows them. WNlsf is the week of
    the leap counted modulo 256 and within 127 weeks of now (IS-GPS-200 20.3.3.5.2.4), DN its day (1 = Sunday);
    the leap happens at the end of that day, UTC. Returns True if that changed anything. Values that make no sense
    (a damaged page) are ignored - the table stays."""
    global _LEAPS
    try:
        dtls, dtlsf, wn, dn = int(p['dtLS']), int(p['dtLSF']), int(p['WNlsf']), int(p['DN'])
    except (KeyError, TypeError, ValueError):
        return False
    if not (15 <= dtls <= 60 and 15 <= dtlsf <= 60 and abs(dtlsf - dtls) <= 1 and 1 <= dn <= 7):
        return False
    import nav_bits
    week = nav_bits.full_week(wn, 256, float(t_rx) + float(gpst_minus_utc(t_rx)))
    t_leap = GPS_EPOCH + week * WEEK_S + dn * 86400            # a UTC midnight as Unix time (GPS_EPOCH is one)
    if dtls == dtlsf:       # no leap pending: WNlsf is the last one's week, whose modulo-256 count can't be trusted
        t_leap = min(t_leap, float(t_rx))
    UTC_SOURCE[0] = 'decoded: GPST - UTC %d s%s' % (
        dtls, '' if dtlsf == dtls else ', %d s from %s' % (dtlsf, time_text(t_leap)))
    if int(gpst_minus_utc(t_leap)) == dtlsf and (dtls == dtlsf or int(gpst_minus_utc(t_leap - 1)) == dtls):
        return False                                            # (the table already has it)
    table = [(s, v) for s, v in _LEAPS if s < t_leap] + [(t_leap, dtlsf)]
    _LEAPS = table
    return True


def time_text(t_unix):
    import datetime
    return datetime.datetime.fromtimestamp(t_unix, datetime.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')


def gpst_to_unix(t_gpst):
    """RTKLIB-style GPST seconds -> Unix (UTC) seconds"""
    t = np.asarray(t_gpst, float)
    return t - gpst_minus_utc(t - gpst_minus_utc(t))       # (the inner step: GPST read as Unix time, ~18 s off)


def unix_to_gpst(t_unix):
    return np.asarray(t_unix, float) + gpst_minus_utc(t_unix)


def is_geo(eph):
    return eph['sys'] == 'C' and eph['prn'] in BDS_GEO_PRNS


URA_M = [2.4, 3.4, 4.85, 6.85, 9.65, 13.65, 24.0, 48.0, 96.0, 192.0, 384.0, 768.0, 1536.0, 3072.0, 6144.0]


def accuracy_m(sys, index):
    """the satellite's own forecast of its range accuracy (m): GPS URA index (IS-GPS-200 20.3.3.3.1.3), Galileo SISA
    index (OS SIS ICD 5.1.12); None when it says it has none (GPS 15, Galileo 255 'no accuracy prediction') or the
    system sends none we read (BeiDou URAI uses the GPS table)"""
    n = int(index)
    if sys == 'E':
        if n <= 49:
            return n * 0.01
        if n <= 74:
            return 0.5 + (n - 50) * 0.02
        if n <= 99:
            return 1.0 + (n - 75) * 0.04
        if n <= 125:
            return 2.0 + (n - 100) * 0.16
        return None
    return URA_M[n] if 0 <= n < len(URA_M) else None


def parse_eph(line):
    """One PocketSDR '$EPH,...' log line -> ephemeris dict, or None if it isn't one we can use
    (another log line, GLONASS, a malformed line)."""
    parts = line.strip().split(',')
    if len(parts) < 4 + len(FIELDS) or parts[0] != '$EPH':
        return None
    sat, sig = parts[2], parts[3]
    if not sat or sat[0] not in CONSTANTS:
        return None
    try:
        prn = int(sat[1:])
        eph = {'rcv_time': float(parts[1])}
        for name, val in zip(FIELDS, parts[4:4 + len(FIELDS)]):
            eph[name] = int(val) if name in INT_FIELDS else float(val)
    except ValueError:
        return None
    if eph['A'] <= 0:
        return None
    return finish_eph(eph, sat[0], prn, sig)


def finish_eph(eph, sys, prn, sig=''):
    """Adds the derived keys to a dict holding the FIELDS values (as $EPH prints them)."""
    if sys == 'E' and sig == 'E1B' and eph['svh'] & ~0x7:
        # an I/NAV ephemeris carries the health of both its signals (rtcm3 1046: E5b status << 3 | E1-B status);
        # we receive E1-B only, so only its own bits (health 2 bits << 1 | data validity) decide - an E5b outage
        # must not throw away a good E1-B satellite. The E5b bits are kept beside it.
        eph['svh_e5b'] = eph['svh'] >> 3
        eph['svh'] &= 0x7
    if sys in 'GJ' and not eph['fit']:
        eph['fit'] = 4.0        # fit interval flag 0 = 4 hours (lines written before 2026-11-07 stored it as 0)
    eph.update({'sys': sys, 'prn': prn, 'sat': '%s%02d' % (sys, prn), 'sig': sig,
                'sqrtA': math.sqrt(eph['A']), 'health': eph['svh'], 'accuracy_m': accuracy_m(sys, eph['sva']),
                'toe': float(gpst_to_unix(eph['toe_gpst'])), 'toc': float(gpst_to_unix(eph['toc_gpst'])),
                'ttr': float(gpst_to_unix(eph['ttr_gpst']))})
    t_sys = eph['toe_gpst'] - (14 if sys == 'C' else 0)          # in the satellite's own time scale
    week = int((t_sys - GPS_EPOCH) // 604800)
    eph['week'] = week - (BDT_WEEK0 if sys == 'C' else GST_WEEK0 if sys == 'E' else 0)
    eph['geo'] = is_geo(eph)
    return eph


def _kepler(eph, t_unix):
    """Satellite position (m, ECEF) and eccentric anomaly at Unix time(s) t."""
    mu, omge = CONSTANTS[eph['sys']]
    tk = unix_to_gpst(t_unix) - eph['toe_gpst']
    A = eph['A']
    n = math.sqrt(mu / A ** 3) + eph['deln']
    M = eph['M0'] + n * tk
    e = eph['e']
    E = M.copy()
    for _ in range(30):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E -= dE
        if np.all(np.abs(dE) < 1e-14):
            break
    nu = np.arctan2(math.sqrt(1 - e * e) * np.sin(E), np.cos(E) - e)
    phi = nu + eph['omg']
    s2, c2 = np.sin(2 * phi), np.cos(2 * phi)
    u = phi + eph['cus'] * s2 + eph['cuc'] * c2
    r = A * (1 - e * np.cos(E)) + eph['crs'] * s2 + eph['crc'] * c2
    i = eph['i0'] + eph['idot'] * tk + eph['cis'] * s2 + eph['cic'] * c2
    x, y = r * np.cos(u), r * np.sin(u)
    ci, si = np.cos(i), np.sin(i)
    if eph['geo']:
        # BeiDou GEO: the elements are in a frame tilted by 5 deg and not rotating with the Earth
        O = eph['OMG0'] + eph['OMGd'] * tk - omge * eph['toes']
        cO, sO = np.cos(O), np.sin(O)
        xg, yg, zg = x * cO - y * ci * sO, x * sO + y * ci * cO, y * si
        so, co = np.sin(omge * tk), np.cos(omge * tk)
        c5, s5 = math.cos(math.radians(-5)), math.sin(math.radians(-5))
        pos = np.stack([xg * co + yg * so * c5 + zg * so * s5,
                        -xg * so + yg * co * c5 + zg * co * s5,
                        -yg * s5 + zg * c5], -1)
    else:
        O = eph['OMG0'] + (eph['OMGd'] - omge) * tk - omge * eph['toes']
        cO, sO = np.cos(O), np.sin(O)
        pos = np.stack([x * cO - y * ci * sO, x * sO + y * ci * cO, y * si], -1)
    return pos, E


def sat_ecef(eph, t_unix):
    """Earth-fixed position (km) and velocity (km/s) at Unix time(s) t: (3,) each for a scalar t,
    (n, 3) for an array."""
    t = np.asarray(t_unix, float)
    tt = np.atleast_1d(t)
    pos = _kepler(eph, tt)[0]
    vel = (_kepler(eph, tt + 0.5)[0] - _kepler(eph, tt - 0.5)[0])
    if t.ndim == 0:
        pos, vel = pos[0], vel[0]
    return pos / 1000.0, vel / 1000.0


def sat_clock(eph, t_unix):
    """Satellite clock offset (s) at Unix time(s) t, relativistic term included (TGD not applied)."""
    mu = CONSTANTS[eph['sys']][0]
    dt = unix_to_gpst(t_unix) - eph['toc_gpst']
    E = _kepler(eph, np.atleast_1d(np.asarray(t_unix, float)))[1]
    rel = -2 * math.sqrt(mu * eph['A']) * eph['e'] * np.sin(E) / C_LIGHT ** 2
    out = eph['af0'] + eph['af1'] * dt + eph['af2'] * dt * dt + rel.reshape(np.shape(dt))
    return out


def _plane(r, v):
    """(node, inclination, argument of latitude) of the inertial state r, v (any units)"""
    h = np.cross(r, v)
    h /= np.linalg.norm(h, axis=-1, keepdims=True)
    inc = np.arccos(np.clip(h[..., 2], -1, 1))
    node = np.arctan2(h[..., 0], -h[..., 1])
    nvec = np.stack([np.cos(node), np.sin(node), np.zeros_like(node)], -1)
    u = np.arctan2(np.sum(np.cross(nvec, r) * h, -1), np.sum(nvec * r, -1))
    return node, inc, u


def _inertial(eph, t):
    """position (km) and velocity (km/s) in orbitfit's inertial frame (Earth turned back by GMST)"""
    t = np.atleast_1d(np.asarray(t, float))
    p, v = sat_ecef(eph, t)
    w = 2 * math.pi * 1.00273781191135448 / 86400      # the rate orbitfit.gmst turns at
    ri = orbitfit.ecef_to_eci(p, t)
    vi = orbitfit.ecef_to_eci(v, t) + w * np.stack([-ri[:, 1], ri[:, 0], np.zeros(len(t))], -1)
    return ri, vi


class EphemerisOrbit:
    """A satellite's broadcast ephemeris, wrapped to look like orbitfit.Orbit to the rest of the app.

    p = [node, inclination, u0, n, a] at t0 = Toe like orbitfit's circular model: node/inclination/
    argument of latitude of the inertial state at Toe (rad), n the mean angular rate (rad/s), a the
    semi-major axis (km). ecef() itself always uses the full broadcast model, not p."""

    def __init__(self, eph, received=None):
        self.eph = eph
        self.t0 = eph['toe']
        mu = CONSTANTS[eph['sys']][0]
        n = math.sqrt(mu / eph['A'] ** 3) + eph['deln']
        ri, vi = _inertial(eph, self.t0)
        node, inc, u0 = (float(x[0]) for x in _plane(ri, vi))
        self.p = np.array([node, inc, u0 % (2 * math.pi), n, eph['A'] / 1000.0])
        self.cov = None
        t_rx = eph['ttr'] if received is None else received
        self.info = {'fixes': 0, 'passes': 1, 'doppler_groups': 0, 'doppler_points': 0,
                     'first_obs': float(t_rx), 'last_obs': float(t_rx), 'fix_rms_km': 0.0,
                     'ephemeris': True, 'sat': eph['sat'], 'sig': eph['sig'], 'toe': eph['toe'],
                     'health': eph['health'], 'geo': eph['geo'],
                     # its own accuracy forecast (URA / SISA, m) and fit interval (h; GPS only, 0 = not sent)
                     'accuracy_m': eph.get('accuracy_m'), 'fit_h': eph.get('fit', 0.0)}

    def ecef(self, t):
        # position only (km): every caller (map, globe, sky view, Integrity's Doppler check by differencing two
        # positions) throws the velocity away, and sat_ecef's velocity costs two more Kepler solves - measured
        # 255 us per scalar call vs 89 us without it (measure_c.py), 57 GNSS orbits on every scene rebuild.
        # sat_ecef keeps the velocity for whoever needs it.
        return _kepler(self.eph, np.atleast_1d(np.asarray(t, float)))[0] / 1000.0

    def sigma_km(self, t):
        """1-sigma position uncertainty (km) at time(s) t. Model: SIG_FIT_KM within SIG_HOURS of Toe
        (the fit interval: GPS 4 h, Galileo and BeiDou are sent fresher still), then growing by
        SIG_RATE_KM per hour beyond it plus SIG_GROWTH_KM * hours^2. Checked against precise (SP3)
        orbits for 96 GPS/Galileo/BeiDou satellites over -1..+4 days: the error is 1-3 m in the fit
        interval, worst ~0.2 km at +-4 h, ~3 km at a day (BeiDou GEO/IGSO; 0.6 km median), ~15 km at
        4 days, all but 0.3 % inside 3 sigma. The model is 4.6 km at a day, 54 km at 4 days, 150 km
        (drawn faint) after about a week and 1500 km (hidden) after about three weeks."""
        x = np.maximum(np.abs(np.atleast_1d(np.asarray(t, float)) - self.t0) / 3600.0 - SIG_HOURS, 0.0)
        return SIG_FIT_KM + SIG_RATE_KM * x + SIG_GROWTH_KM * x * x

    def plane_at(self, t):
        """(node, inclination) in radians of the orbit plane at time t, in orbitfit's inertial frame"""
        ri, vi = _inertial(self.eph, t)
        node, inc, _ = _plane(ri, vi)
        if np.ndim(t) == 0:
            return float(node[0]), float(inc[0])
        return node, inc

    @property
    def plane_sigma_km(self):
        return SIG_FIT_KM

    @property
    def period_min(self):
        return 2 * math.pi / self.p[3] / 60


# ---------------------------------------------------------------- almanacs (nav_bits.py): rough orbits of every satellite
# Every GPS / Galileo satellite sends the almanac of the whole constellation: a Kepler orbit without the small
# corrections, good to a few km for days. It is used only to know where a satellite is roughly - which ones are up
# but not heard, which to search for, the next pass - never drawn as a satellite's orbit (that is its ephemeris).
# Measured 2026-11-07 (build\_tmp\later\orbits\alm_test.py: 6 sessions 11-03..11-07, against the same satellites'
# own ephemerides, the largest difference over toe-2h..toe+2h): GPS almanacs within 3 days of their reference time
# median 1.0-1.9 km, worst 4.7 km; 3-3.3 days 2.4 / 6.0 km; 5.3-6.3 days 8.2 / 12.8 km. Galileo's (a new set every
# ~10 min, t0a minutes before it is sent) 1.5 km fresh, then fast worse: 5.3 km median / 8.5 km worst at 0.5-0.8 days.
ALM_SIG_KM = 2.0            # 1-sigma near the reference time...
ALM_DAYS = {'G': 3.0, 'E': 0.0}         # ...for this long either side (days)...
ALM_RATE_KM = {'G': 2.0, 'E': 8.0}      # ...then growing this much per day (manoeuvres; nothing models drag)
ALM_KEEP_DAYS = 28.0        # older almanacs are not used at all (gnss-solver.json drops them too)


def almanac_eph(a, sat):
    """an almanac entry (nav_bits.parse_alm: GPS / Galileo) -> the dict _kepler / sat_ecef read"""
    toa = float(a['toa_gpst'])
    return {'sys': sat[0], 'prn': int(sat[1:]), 'sat': sat, 'A': a['A'], 'e': a['e'], 'i0': a['i0'],
            'OMG0': a['OMG0'], 'omg': a['omg'], 'M0': a['M0'], 'OMGd': a['OMGd'], 'deln': 0.0, 'idot': 0.0,
            'crc': 0.0, 'crs': 0.0, 'cuc': 0.0, 'cus': 0.0, 'cic': 0.0, 'cis': 0.0, 'toe_gpst': toa,
            'toes': (toa - GPS_EPOCH) % WEEK_S, 'toe': float(gpst_to_unix(toa)), 'geo': False,
            'af0': a.get('af0', 0.0), 'af1': a.get('af1', 0.0), 'af2': 0.0, 'toc_gpst': toa, 'svh': a['svh']}


class AlmanacOrbit:
    """a satellite's almanac with the orbit interface (ecef, sigma_km, plane_at) - deliberately NOT an
    EphemerisOrbit, so nothing that wants a real orbit (the map, the geometry, Integrity) mistakes it for one.
    SBAS almanacs (message type 17: a position and velocity to 2.6 km / 10 m/s) hold their position."""

    def __init__(self, text, sat, received):
        import nav_bits
        a = nav_bits.parse_alm(text)
        if a is None:
            raise ValueError('not an almanac: %r' % text[:40])
        self.alm, self.sat = a, sat
        self.sbas = a['sys'] == 'S'
        if self.sbas:
            self.t0 = a['t0']
            self.pos = np.array([a['x'], a['y'], a['z']]) / 1000.0
            self.eph = None
        else:
            self.eph = almanac_eph(a, sat)
            self.t0 = self.eph['toe']
        self.info = {'almanac': True, 'sat': sat, 'toa': self.t0, 'health': a['svh'], 'received': float(received)}

    def ecef(self, t):
        t = np.atleast_1d(np.asarray(t, float))
        if self.sbas:
            return np.repeat(self.pos[None, :], len(t), 0)
        return _kepler(self.eph, t)[0] / 1000.0

    def sigma_km(self, t):
        days = np.abs(np.atleast_1d(np.asarray(t, float)) - self.t0) / 86400.0
        if self.sbas:       # (2.6 km steps; a geostationary satellite's daily figure-eight is tens of km)
            return 5.0 + 0.0 * days
        k = self.sat[0]
        return ALM_SIG_KM + ALM_RATE_KM.get(k, 8.0) * np.maximum(days - ALM_DAYS.get(k, 0.0), 0.0)

    def usable(self, t):
        return abs(float(t) - self.t0) < ALM_KEEP_DAYS * 86400

    def plane_at(self, t):
        if self.sbas:
            return 0.0, 0.0
        ri, vi = _inertial(self.eph, t)
        node, inc, _ = _plane(ri, vi)
        if np.ndim(t) == 0:
            return float(node[0]), float(inc[0])
        return node, inc

    @property
    def period_min(self):
        if self.sbas:
            return GEO_PERIOD_MIN
        return 2 * math.pi / math.sqrt(CONSTANTS[self.sat[0]][0] / self.eph['A'] ** 3) / 60


# ---------------------------------------------------------------- the ephemeris archive
# Every ephemeris a satellite sent is kept (catalog 'eph_archive'; ~145 distinct ones over 3.4 days, ~18 kB a day),
# and each moment uses the one whose reference time is nearest: scrubbing the timeline to yesterday draws a
# satellite with yesterday's ephemeris, not today's extrapolated back (3 km wrong at a day, 15 km at 4 days).
class _Archive:
    """mixed in front of an orbit class: several decoded orbits of one satellite, each time sample served by the
    member whose reference time t0 is nearest. Looks like the newest member (eph / t0 / p / info come from it), so
    everything that reads those - and isinstance checks - keeps working. Members are built when first needed
    (a GLONASS one integrates its motion: ~30 ms)."""

    def _setup(self, members):
        """members: [(t0, received, factory)] with factory() -> an orbit of the base class"""
        best = {}
        for t0, rx, fac in members:
            if t0 not in best or rx > best[t0][0]:
                best[t0] = (rx, fac)                # (two at one reference time: the later reception)
        self._t0s = np.array(sorted(best))
        self._fac = [best[t][1] for t in self._t0s]
        self._built = {}
        newest = self._member(len(self._t0s) - 1)
        self.__dict__.update({k: v for k, v in newest.__dict__.items() if not k.startswith('_')})
        self.info = dict(newest.info)
        rxs = [best[t][0] for t in self._t0s]
        self.info.update({'archive': len(self._t0s), 'archive_first_toe': float(self._t0s[0]),
                          'archive_last_toe': float(self._t0s[-1]), 'archive_toes': [float(t) for t in self._t0s],
                          'first_obs': float(min(rxs)), 'last_obs': float(max(rxs))})

    def _member(self, k):
        o = self._built.get(k)
        if o is None:
            o = self._built[k] = self._fac[k]()
        return o

    def _pick(self, t):
        """index of the nearest-t0 member per sample"""
        i = np.searchsorted(self._t0s, t)
        i = np.clip(i, 1, len(self._t0s) - 1) if len(self._t0s) > 1 else np.zeros_like(i)
        if len(self._t0s) > 1:
            left = self._t0s[i - 1]
            i = np.where(np.abs(t - left) <= np.abs(self._t0s[i] - t), i - 1, i)
        return i

    def _each(self, t, fn):
        t = np.atleast_1d(np.asarray(t, float))
        idx = self._pick(t)
        ks = np.unique(idx)
        if len(ks) == 1:
            return fn(self._member(int(ks[0])), t)
        out = None
        for k in ks:
            m = idx == k
            got = np.asarray(fn(self._member(int(k)), t[m]))
            if out is None:
                out = np.zeros((len(t),) + got.shape[1:], got.dtype)
            out[m] = got
        return out

    def member_at(self, t):
        """the decoded orbit used at time t (scalar)"""
        return self._member(int(self._pick(np.atleast_1d(float(t)))[0]))

    def ecef(self, t):
        return self._each(t, lambda o, tt: o.ecef(tt))

    def sigma_km(self, t):
        return self._each(t, lambda o, tt: np.asarray(o.sigma_km(tt), float).reshape(-1))

    def plane_at(self, t):
        if np.ndim(t) == 0:
            return self.member_at(t).plane_at(float(t))
        got = self._each(t, lambda o, tt: np.stack(o.plane_at(tt), -1) if len(tt) > 1 else
                         np.array([o.plane_at(float(tt[0]))]))
        return got[:, 0], got[:, 1]


class MultiEphemerisOrbit(_Archive, EphemerisOrbit):
    def __init__(self, members):
        self._setup(members)


def archive_key(fields):
    """what makes a GEPH line's ephemeris a different one (catalog 'eph_archive' and the tracker dedupe by it):
    GPS / Galileo / BeiDou (IODE, toe); GLONASS its tb; WAAS (a new message type 9 every ~2 min, a geostationary
    satellite that hardly moves) one per hour"""
    f = fields.replace(' ', ',').split(',')
    if f[0] == 'SBAS':
        return 'SBAS %d' % (float(f[1]) // 3600)
    if f[0] == 'GLO':
        return 'GLO %.0f' % float(f[1])
    return '%s %s' % (f[1], f[5])          # (sig, IODE, IODC, SVA, SVH, toe_gpst, ...)


def archive_orbit(sat, entries):
    """a satellite's kept orbits [(received t, GEPH fields text)] -> one orbit object: the one orbit when there is
    one, else an archive wrapper choosing the nearest reference time. Entries that don't parse are left out;
    ValueError if none does."""
    members, kind = [], None
    for t, fields in entries:
        fields = fields.replace(' ', ',')
        try:
            if fields.startswith('SBAS'):
                k, t0 = 'S', float(fields.split(',')[1])
                fac = (lambda f=fields, t=t: SBASOrbit.from_fields(f, sat, t))
            elif fields.startswith('GLO'):
                k = 'R'
                t0 = glo_check_toe(float(fields.split(',')[1]), float(t))
                if t0 is None:
                    continue
                fac = (lambda f=fields, t=t: GlonassOrbit.from_fields(f, sat, t))
            else:
                k = 'K'
                eph = parse_eph('$EPH,%.3f,%s,%s' % (t, sat, fields))
                if eph is None:
                    continue
                t0 = eph['toe']
                fac = (lambda e=eph: EphemerisOrbit(e))
        except (ValueError, IndexError):
            continue
        members.append((k, t0, float(t), fac))
    if not members:
        raise ValueError('no usable ephemeris')
    kind = max(members, key=lambda m: m[2])[0]          # (one satellite sends one kind; the newest decides)
    members = [m[1:] for m in members if m[0] == kind]
    if len({m[0] for m in members}) == 1:
        return max(members, key=lambda m: m[1])[2]()
    cls = {'K': MultiEphemerisOrbit, 'S': MultiSBASOrbit, 'R': MultiGlonassOrbit}[kind]
    return cls(members)


# ---------------------------------------------------------------- SBAS (WAAS) geostationary satellites
# SBAS satellites send GPS-like signals on L1 and, every couple of minutes, their own position, velocity and
# acceleration (message type 9, RTCA DO-229) - so they are placed from what you receive, like the others.
SBAS_PREAMBLES = (0x53, 0x9A, 0xC6)
GEO_PERIOD_MIN = 1436.068          # one sidereal day


def _crc24q_bits(bits):
    crc = 0
    for b in bits:
        crc ^= b << 23
        crc = ((crc << 1) ^ 0x1864CFB) & 0xFFFFFF if crc & 0x800000 else (crc << 1) & 0xFFFFFF
    return crc


def sbas_message(bits):
    """a 250-bit SBAS message (list of 0/1, preamble first) -> (type, data bits) if its CRC checks, else None"""
    if len(bits) < 250:
        return None
    pre = int(''.join(map(str, bits[:8])), 2)
    if pre not in SBAS_PREAMBLES:
        return None
    if _crc24q_bits(bits[:226]) != int(''.join(map(str, bits[226:250])), 2):
        return None
    return int(''.join(map(str, bits[8:14])), 2), bits[14:226]


def _u(bits, i, n):
    return int(''.join(map(str, bits[i:i + n])), 2)


def _s(bits, i, n):
    x = _u(bits, i, n)
    return x - (1 << n) if x & (1 << (n - 1)) else x


def sbas_geo_nav(data, t_rx):
    """message type 9's 212 data bits -> {t0, x, y, z (m), vx..(m/s), ax..(m/s2)}; t_rx (unix) resolves the day"""
    i = 8                                   # (IODN / reserved)
    t0_tod = _u(data, i, 13) * 16.0
    i += 13 + 4                             # (and the URA index)
    x, y, z = _s(data, i, 30) * 0.08, _s(data, i + 30, 30) * 0.08, _s(data, i + 60, 25) * 0.4
    i += 85
    vx, vy, vz = _s(data, i, 17) * 0.000625, _s(data, i + 17, 17) * 0.000625, _s(data, i + 34, 18) * 0.004
    i += 52
    ax, ay, az = _s(data, i, 10) * 0.0000125, _s(data, i + 10, 10) * 0.0000125, _s(data, i + 20, 10) * 0.0000625
    # t0 is a GPS time of day: the one nearest the moment it was received
    gps = t_rx + gpst_minus_utc(t_rx)
    t0 = math.floor(gps / 86400) * 86400 + t0_tod
    t0 += 86400 * round((gps - t0) / 86400)
    return {'t0': t0 - gpst_minus_utc(t_rx), 'x': x, 'y': y, 'z': z, 'vx': vx, 'vy': vy, 'vz': vz,
            'ax': ax, 'ay': ay, 'az': az}


SBAS_FIELDS = ['t0', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax', 'ay', 'az']


class SBASOrbit:
    """a geostationary SBAS satellite from its own message type 9, looking like orbitfit.Orbit to the app.
    Its position is extrapolated with the broadcast velocity and acceleration for up to SBAS_HOLD_S either
    side of t0, then held (a geostationary satellite drifts only slowly in the Earth-fixed frame)."""
    SBAS_HOLD_S = 1800.0

    def __init__(self, nav, sat, received):
        self.nav, self.t0 = nav, nav['t0']
        r = np.array([nav['x'], nav['y'], nav['z']]) / 1000.0
        a = float(np.linalg.norm(r))
        self.p = np.array([0.0, 0.0, 0.0, 2 * math.pi / (GEO_PERIOD_MIN * 60), a])
        node, inc = self.plane_at(self.t0)
        self.p[0], self.p[1] = node, inc
        self.cov = None
        self.eph = {'e': 0.0}
        self.info = {'fixes': 0, 'passes': 1, 'doppler_groups': 0, 'doppler_points': 0,
                     'first_obs': float(received), 'last_obs': float(received), 'fix_rms_km': 0.0,
                     'ephemeris': True, 'sbas': True, 'sat': sat, 'sig': 'L1 SBAS', 'toe': self.t0,
                     'health': 0, 'geo': True}

    @classmethod
    def from_fields(cls, text, sat, received):
        """'SBAS,t0,x,y,z,vx,vy,vz,ax,ay,az' (a GEPH line's fields) -> SBASOrbit"""
        v = [float(x) for x in text.split(',')[1:]]
        return cls(dict(zip(SBAS_FIELDS, v)), sat, received)

    def ecef(self, t):
        t = np.atleast_1d(np.asarray(t, float))
        dt = np.clip(t - self.t0, -self.SBAS_HOLD_S, self.SBAS_HOLD_S)[:, None]
        n = self.nav
        p = np.array([n['x'], n['y'], n['z']])
        v = np.array([n['vx'], n['vy'], n['vz']])
        a = np.array([n['ax'], n['ay'], n['az']])
        return (p + v * dt + 0.5 * a * dt * dt) / 1000.0

    def sigma_km(self, t):
        """1-sigma (km): metres within 15 min of t0, then ~2 km per hour of holding it still (the slow daily
        figure-eight a geostationary satellite traces)"""
        h = np.maximum(np.abs(np.atleast_1d(np.asarray(t, float)) - self.t0) / 3600.0 - 0.25, 0.0)
        return 0.01 + 2.0 * np.minimum(h, 12.0)

    def plane_at(self, t):
        t = np.atleast_1d(np.asarray(t, float))
        p = self.ecef(t)
        n = self.nav
        v = np.repeat([[n['vx'] / 1000.0, n['vy'] / 1000.0, n['vz'] / 1000.0]], len(t), 0)
        w = 2 * math.pi * 1.00273781191135448 / 86400
        ri = orbitfit.ecef_to_eci(p, t)
        vi = orbitfit.ecef_to_eci(v, t) + w * np.stack([-ri[:, 1], ri[:, 0], np.zeros(len(t))], -1)
        node, inc, _ = _plane(ri, vi)
        if len(t) == 1:
            return float(node[0]), float(inc[0])
        return node, inc

    @property
    def plane_sigma_km(self):
        return 0.01

    @property
    def period_min(self):
        return GEO_PERIOD_MIN


# ---------------------------------------------------------------- GLONASS
# GLONASS broadcasts, every 30 minutes, its position, velocity and the Moon / Sun acceleration at a reference
# time tb (in the Earth-fixed PZ-90 frame, practically WGS-84); receivers integrate the motion from there under
# the Earth's gravity with its J2 term (GLONASS ICD, A.3.1.2; as in RTKLIB's glorbit).
GLO_MU, GLO_J2, GLO_RE, GLO_OMGE = 3.9860044e14, 1.0826257e-3, 6378136.0, 7.292115e-5
GLO_FIELDS = ['toe', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'ax', 'ay', 'az', 'taun', 'gamn', 'fcn', 'slot']
GLO_STEP_S = 30.0
GLO_TB_SLIP_S = 32 * 900.0        # see glo_check_toe
GLO_TB_TOL_S = 5400.0


def glo_check_toe(toe, received):
    """the reference time tb of a GLONASS broadcast, checked against when it was received: a broadcast is for the
    next quarter / half hour, so tb is within minutes of reception. Some come out of the decoder exactly 8 h late
    (tb's 32 bit set: seen for every broadcast sent 15:00-24:00 Moscow time, while position and velocity are
    right - checked against the public orbits to 1-3 km); those are put back. Anything else that doesn't fit
    -> None (not used)"""
    d = (toe - received + 43200.0) % 86400.0 - 43200.0
    if abs(d) <= GLO_TB_TOL_S:
        return received + d
    if abs(d - GLO_TB_SLIP_S) <= GLO_TB_TOL_S:
        return received + d - GLO_TB_SLIP_S
    return None
GLO_SPAN_S = 8 * 3600.0           # integrated this far either side of tb (a little over half an orbit... and more)


def _glo_deq(x, acc):
    r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2]
    r3 = r2 * math.sqrt(r2)
    a = 1.5 * GLO_J2 * GLO_MU * GLO_RE * GLO_RE / r2 / r3
    b = 5.0 * x[2] * x[2] / r2
    c = -GLO_MU / r3 - a * (1.0 - b)
    w2 = GLO_OMGE * GLO_OMGE
    return np.array([x[3], x[4], x[5],
                     (c + w2) * x[0] + 2 * GLO_OMGE * x[4] + acc[0],
                     (c + w2) * x[1] - 2 * GLO_OMGE * x[3] + acc[1],
                     (c - 2 * a) * x[2] + acc[2]])


def _glo_rk4(x, acc, h):
    k1 = _glo_deq(x, acc)
    k2 = _glo_deq(x + k1 * h / 2, acc)
    k3 = _glo_deq(x + k2 * h / 2, acc)
    k4 = _glo_deq(x + k3 * h, acc)
    return x + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6


class GlonassOrbit:
    """a GLONASS satellite from its own broadcast ephemeris, looking like orbitfit.Orbit to the app. The motion is
    integrated once (RK4, 30 s steps, GLO_SPAN_S either side of tb) and interpolated (Hermite, with the
    velocities); the broadcast Moon / Sun acceleration is only meant for +-15 min, so it is applied only there."""

    def __init__(self, g, sat, received):
        self.g, self.t0 = g, g['toe']
        x0 = np.array([g['x'], g['y'], g['z'], g['vx'], g['vy'], g['vz']])
        acc = np.array([g['ax'], g['ay'], g['az']])
        n = int(GLO_SPAN_S / GLO_STEP_S)
        fwd, back = [x0], [x0]
        for k in range(n):
            a = acc if (k + 1) * GLO_STEP_S <= 900 else np.zeros(3)
            fwd.append(_glo_rk4(fwd[-1], a, GLO_STEP_S))
            back.append(_glo_rk4(back[-1], a, -GLO_STEP_S))
        self.states = np.array(back[::-1] + fwd[1:])           # (2n + 1, 6), m and m/s
        self.ts = self.t0 + GLO_STEP_S * np.arange(-n, n + 1)
        r = float(np.linalg.norm(x0[:3])) / 1000.0
        vi = (x0[3:] + np.cross([0.0, 0.0, GLO_OMGE], x0[:3])) / 1000.0       # inertial velocity, km/s
        a_km = 1.0 / (2.0 / r - float(vi @ vi) / (GLO_MU / 1e9))              # (vis-viva)
        self.p = np.array([0.0, 0.0, 0.0, math.sqrt(GLO_MU / 1e9 / a_km ** 3), a_km])
        node, inc = self.plane_at(self.t0)
        self.p[0], self.p[1] = node, inc
        self.cov = None
        self.eph = {'e': 0.0}
        self.info = {'fixes': 0, 'passes': 1, 'doppler_groups': 0, 'doppler_points': 0,
                     'first_obs': float(received), 'last_obs': float(received), 'fix_rms_km': 0.0,
                     'ephemeris': True, 'glonass': True, 'sat': sat, 'sig': 'L1OF (channel %+d)' % int(g['fcn']),
                     'toe': self.t0, 'health': 0, 'geo': False}

    @classmethod
    def from_fields(cls, text, sat, received):
        """'GLO,toe,x,y,z,vx,vy,vz,ax,ay,az,taun,gamn,fcn,slot' (a GEPH line's fields) -> GlonassOrbit"""
        g = dict(zip(GLO_FIELDS, [float(x) for x in text.split(',')[1:]]))
        toe = glo_check_toe(g['toe'], float(received))     # (also repairs ones saved before the check existed)
        if toe is None:
            raise ValueError('its broadcast time does not match when it was received')
        g['toe'] = toe
        return cls(g, sat, received)

    def _state(self, t):
        t = np.clip(np.atleast_1d(np.asarray(t, float)), self.ts[0], self.ts[-1])
        i = np.clip(((t - self.ts[0]) // GLO_STEP_S).astype(int), 0, len(self.ts) - 2)
        h = GLO_STEP_S
        s = ((t - self.ts[i]) / h)[:, None]
        p0, p1 = self.states[i, :3], self.states[i + 1, :3]
        v0, v1 = self.states[i, 3:], self.states[i + 1, 3:]
        h00, h10, h01, h11 = 2 * s ** 3 - 3 * s ** 2 + 1, s ** 3 - 2 * s ** 2 + s, -2 * s ** 3 + 3 * s ** 2, s ** 3 - s ** 2
        pos = h00 * p0 + h10 * h * v0 + h01 * p1 + h11 * h * v1
        vel = v0 + (v1 - v0) * s
        return pos, vel

    def ecef(self, t):
        return self._state(t)[0] / 1000.0

    def sigma_km(self, t):
        """1-sigma (km): metres within 15 min of tb (what the broadcast is made for), then growing - the Moon and
        Sun pull is left out beyond that: ~1 km after 2 h, ~5 km after 6 h; past the integrated span, hidden"""
        dt = np.abs(np.atleast_1d(np.asarray(t, float)) - self.t0)
        h = np.maximum(dt / 3600.0 - 0.25, 0.0)
        return np.where(dt > GLO_SPAN_S, 5000.0, 0.005 + 0.3 * h + 0.1 * h * h)

    def plane_at(self, t):
        t = np.atleast_1d(np.asarray(t, float))
        p, v = self._state(t)
        w = 2 * math.pi * 1.00273781191135448 / 86400
        ri = orbitfit.ecef_to_eci(p / 1000.0, t)
        vi = orbitfit.ecef_to_eci(v / 1000.0, t) + w * np.stack([-ri[:, 1], ri[:, 0], np.zeros(len(t))], -1)
        node, inc, _ = _plane(ri, vi)
        if len(t) == 1:
            return float(node[0]), float(inc[0])
        return node, inc

    @property
    def plane_sigma_km(self):
        return 0.005

    @property
    def period_min(self):
        return 2 * math.pi / self.p[3] / 60


class MultiSBASOrbit(_Archive, SBASOrbit):
    def __init__(self, members):
        self._setup(members)


class MultiGlonassOrbit(_Archive, GlonassOrbit):
    def __init__(self, members):
        self._setup(members)

"""The GNSS receiver's line formats in one place (gnss_rx writes them; the catalog, the tracker, the console and
Integrity read them).

Every line is 'KIND: unix_t field field ...'. catalog.read_gnss gives (kind, t, fields after t); the functions here
take those fields. Field ORDER is the contract: a new field only ever goes on the end (older lines simply lack it).

    GCH:  t sat sig cn0 doppler_hz lock_s [sync nav err lol]
    GEPH: t sat <$EPH fields | GLO,... | SBAS,...>
    GPOS: t lat lon height_m nsat                                    (pocket_trk's own fix)
    GSPP: t lat lon height_m nsat hdop n_waas_iono n_fast n_longterm n_smoothed prior
                                                                     (this program's solver, spp.py; 'prior' 1 = the
                                                                      static antenna's height row held it; lines
                                                                      before 2026-11-07 have no 'prior')
    GGAP: t gap_ms source                                            (a sample-stream hole gnss_rx zero-filled)
    GLOG: t EVENT [sat] [values...]                                  (FOUND LOST FALSELOCK EPOCH SYNC START RESTART)
"""

# a GSPP fix with fewer satellites smoothed this long is a 'quick' one (spp.SETTLE_S): shown, but kept out of the
# accuracy figures - the code loops are still pulling in (11 m and -31 m errors beside the Garmin without it)
QUICK_BELOW_SMOOTHED = 4

GSPP_FIELDS = ('lat', 'lon', 'height_m', 'nsat', 'hdop', 'n_waas_iono', 'n_fast', 'n_longterm', 'n_smoothed', 'prior')


def gspp_line(t, s, n_iono, n_fast, n_lt, n_smoothed, prior):
    """the GSPP line for one spp.Solver solution s (gnss_rx)"""
    return 'GSPP: %.3f %.9f %.9f %.3f %d %.2f %d %d %d %d %d' % (
        t, s['lat'], s['lon'], s['h'], s['nsat'], s['hdop'], n_iono, n_fast, n_lt, n_smoothed, 1 if prior else 0)


def gspp_fix(t, f):
    """a GSPP line's fields (after t) -> the fix dict the catalog and the tracker keep:
    {t, lat, lon, height_m, nsat, hdop[, quick][, prior][, waas]}. Raises ValueError/IndexError on a damaged line."""
    fx = {'t': t, 'lat': float(f[0]), 'lon': float(f[1]), 'height_m': float(f[2]), 'nsat': int(f[3]),
          'hdop': float(f[4])}
    if len(f) >= 9 and int(f[8]) < QUICK_BELOW_SMOOTHED:
        fx['quick'] = True          # before carrier smoothing settled
    if len(f) >= 10 and f[9] == '1':
        fx['prior'] = True          # the height row held it: its height is the antenna's, not a measurement
    if len(f) >= 7 and (int(f[5]) > 0 or int(f[6]) > 0):
        fx['waas'] = True           # WAAS ionosphere or fast corrections were applied
    return fx


def gpos_fix(t, f):
    """a GPOS line's fields (after t) -> {t, lat, lon, height_m, nsat}"""
    return {'t': t, 'lat': float(f[0]), 'lon': float(f[1]), 'height_m': float(f[2]), 'nsat': int(f[3])}


def gch(f):
    """a GCH line's fields (after t) -> (sat, sig, cn0, doppler_hz, lock_s or None)"""
    lock = float(f[4]) if len(f) >= 5 else None
    return f[0], f[1], float(f[2]), float(f[3]), lock


def ggap(f):
    """a GGAP line's fields (after t) -> (ms, source)"""
    return float(f[0]), (f[1] if len(f) > 1 else '?')

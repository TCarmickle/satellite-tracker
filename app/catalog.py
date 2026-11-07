"""Running catalog of every satellite heard across all captures (Iridium and Orbcomm).

Usage: python catalog.py            (adds any new captures, then writes catalog.txt)
       python catalog.py --redo     (re-analyzes every capture from scratch)

Each capture is analyzed once and stored in catalog-data.json, so results don't change
later when the satellite orbit data (TLEs) gets older or is refreshed.
Satellites are named two ways:
  position  - a ring alert reported its position; the satellite overhead there is named
  Doppler   - the frequency curve of a plane/slot's ITL frames matches one satellite's orbit
Orbcomm satellites send their own position (ephemeris packets); the satellite at that spot is named,
and that name sticks to the satellite ID in its other packets.
"""
import collections
import datetime
import glob
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gnss_lines  # noqa: E402  (the GNSS receiver's line formats)
import paths  # noqa: E402  (where the data lives)
import storage  # noqa: E402
import systems  # noqa: E402
# iridium-toolkit: inside this folder (a packaged copy) or next to it
TOOLKIT = next((d for d in (os.path.join(HERE, 'iridium-toolkit'), os.path.join(HERE, '..', 'iridium-toolkit'))
                if os.path.exists(os.path.join(d, 'iridium-parser.py'))), os.path.join(HERE, 'iridium-toolkit'))
PARSER = os.path.join(TOOLKIT, 'iridium-parser.py')
# public orbit list, used only to put names on satellites (positions always come from our own fixes);
# kept fresh from CelesTrak, with SatDump's copy as a fallback
TLE_FILE = os.path.join(paths.DATA, 'iridium-names.tle')
TLE_URL = 'https://celestrak.org/NORAD/elements/gp.php?GROUP=%s&FORMAT=tle'


def tle_file(system='iridium'):
    """the name list of a system (both ham bands share the ham one)"""
    s = systems.get(system)
    return os.path.join(paths.DATA, '%s-names.tle' % s.get('system', s['key']))


SATDUMP_TLES = os.path.join(os.environ.get('APPDATA', ''), 'satdump', 'satdump_tles.txt')
# one file per capture in data\catalog, and an index of them (its time changes with every save: the tracker reloads)
CAT_DIR = os.path.join(paths.DATA, 'catalog')
DATA = os.path.join(CAT_DIR, 'index.json')
LEGACY = os.path.join(paths.DATA, 'catalog-data.json')                # the old single file (read until migrated)
PREVIOUS = os.path.join(paths.DATA, 'catalog-data.previous.json')     # ...and its last good copy
REPORT = os.path.join(paths.DATA, 'catalog.txt')
SETTINGS = os.path.join(paths.DATA, 'settings.json')
# one catalog run at a time: the tracker's Stop, the console's Stop and a --redo by hand can all start one, and two
# runs saving at once used to delete each other's new entries. The lock holds the PID of the run that owns it.
LOCK = os.path.join(CAT_DIR, '.lock')
LOCK_WAIT_S = 60
# the capture names the index listed when this process loaded the catalog - the only files save_data may delete
_LOADED = set()


def load_home():
    """(lat, lon, altitude km) of the antenna from settings.json, or None if not set up yet"""
    try:
        lat, lon, alt = json.load(open(SETTINGS))['home']
        return float(lat), float(lon), float(alt)
    except (OSError, ValueError, KeyError, TypeError):
        return None


def save_home(lat, lon, alt_km):
    try:
        settings = json.load(open(SETTINGS))
    except (OSError, ValueError):
        settings = {}
    settings['home'] = [lat, lon, alt_km]
    json.dump(settings, open(SETTINGS, 'w'), indent=1)


HOME = load_home()
UTC = datetime.timezone.utc

# ---------- orbits (TLE mean elements + J2 drift; good to a few km over days) ----------
MU, RE, J2, D2R, C = 398600.4418, 6378.137, 1.08262668e-3, math.pi / 180, 299792.458


class Sat:
    def __init__(self, name, l1, l2):
        self.name = name
        yy, doy = int(l1[18:20]), float(l1[20:32])
        self.epoch = datetime.datetime(2000 + yy, 1, 1, tzinfo=UTC) + datetime.timedelta(days=doy - 1)
        self.i, self.raan = float(l2[8:16]) * D2R, float(l2[17:25]) * D2R
        self.e = float('0.' + l2[26:33].strip())
        self.w, self.m = float(l2[34:42]) * D2R, float(l2[43:51]) * D2R
        self.n = float(l2[52:63]) * 2 * math.pi / 86400
        self.a = (MU / self.n ** 2) ** (1 / 3)
        p = self.a * (1 - self.e ** 2)
        k = 1.5 * self.n * J2 * (RE / p) ** 2
        self.raan_dot = -k * math.cos(self.i)
        self.w_dot = 0.5 * k * (5 * math.cos(self.i) ** 2 - 1)

    def ecef(self, t):
        dt = (t - self.epoch).total_seconds()
        M, raan, w = self.m + self.n * dt, self.raan + self.raan_dot * dt, self.w + self.w_dot * dt
        E = M
        for _ in range(8):
            E -= (E - self.e * math.sin(E) - M) / (1 - self.e * math.cos(E))
        nu = 2 * math.atan2(math.sqrt(1 + self.e) * math.sin(E / 2), math.sqrt(1 - self.e) * math.cos(E / 2))
        r, u = self.a * (1 - self.e * math.cos(E)), w + nu
        cO, sO, ci, si, cu, su = math.cos(raan), math.sin(raan), math.cos(self.i), math.sin(self.i), math.cos(u), math.sin(u)
        x, y, z = r * (cO * cu - sO * su * ci), r * (sO * cu + cO * su * ci), r * su * si
        jd = 2440587.5 + t.timestamp() / 86400
        th = (280.46061837 + 360.98564736629 * (jd - 2451545.0)) % 360 * D2R
        return x * math.cos(th) + y * math.sin(th), -x * math.sin(th) + y * math.cos(th), z


def site(lat, lon, h):
    f = 1 / 298.257223563
    e2, la, lo = f * (2 - f), lat * D2R, lon * D2R
    N = RE / math.sqrt(1 - e2 * math.sin(la) ** 2)
    return ((N + h) * math.cos(la) * math.cos(lo), (N + h) * math.cos(la) * math.sin(lo), (N * (1 - e2) + h) * math.sin(la))


def look(sat, t, lat, lon, h=0.0):
    """elevation (deg), azimuth (deg), range (km), range rate (km/s) from a ground point"""
    o, p = site(lat, lon, h), sat.ecef(t)
    p2 = sat.ecef(t + datetime.timedelta(seconds=1))
    d = [p[k] - o[k] for k in range(3)]
    rng = math.sqrt(sum(v * v for v in d))
    rng2 = math.sqrt(sum((p2[k] - o[k]) ** 2 for k in range(3)))
    la, lo = lat * D2R, lon * D2R
    S = math.sin(la) * math.cos(lo) * d[0] + math.sin(la) * math.sin(lo) * d[1] - math.cos(la) * d[2]
    Ea = -math.sin(lo) * d[0] + math.cos(lo) * d[1]
    Z = math.cos(la) * math.cos(lo) * d[0] + math.cos(la) * math.sin(lo) * d[1] + math.sin(la) * d[2]
    return math.degrees(math.asin(Z / rng)), (math.degrees(math.atan2(Ea, -S)) + 360) % 360, rng, rng2 - rng


def refresh_tles(max_age_days=2, system='iridium'):
    """download a fresh name list now and then (CelesTrak asks for no more than one fetch every couple of hours)"""
    s = systems.get(system)
    path = tle_file(system)
    groups = s.get('tle_groups') or (s['tle_group'],)
    try:            # (the groups the list was made from: a list missing a newly added group is fetched again)
        same = open(path + '.groups').read().split() == list(groups)
    except OSError:
        same = len(groups) == 1
    if same and os.path.exists(path) and time.time() - os.path.getmtime(path) < max_age_days * 86400:
        return
    try:
        import requests
        text = ''
        for group in groups:
            r = requests.get(TLE_URL % group, timeout=15, headers={'User-Agent': 'IridiumTracker/1.0'})
            if not (r.ok and '\n1 ' in r.text and (not s['prefix'] or s['prefix'].strip() in r.text)):
                return          # a group failed: keep the old list rather than half of a new one
            text += r.text.replace('\r\n', '\n').rstrip('\n') + '\n'
        open(path, 'w', newline='\n').write(text)
        open(path + '.groups', 'w').write(' '.join(groups))
    except Exception:
        pass        # offline: keep using what we have


def load_tles(system='iridium'):
    """the public catalog's satellites of one system (for names only)"""
    refresh_tles(system=system)
    prefix = (systems.get(system)['prefix'] or '').strip()
    for path in (tle_file(system), SATDUMP_TLES):
        if os.path.exists(path):
            lines = open(path).read().splitlines()
            sats = [Sat(lines[i].strip(), lines[i + 1], lines[i + 2])
                    for i in range(len(lines) - 2) if lines[i].startswith(prefix) and lines[i + 1].startswith('1 ')
                    and lines[i + 2].startswith('2 ') and not lines[i].startswith(('1 ', '2 '))]
            if sats:
                return sats
    return []


# ---------- reading a capture ----------
def capture_start(stem):
    """(UTC start of the burst clock, how it was known)"""
    try:
        text = open(stem + '.start').read().strip()
        t = datetime.datetime.strptime(text.rstrip('Z'), '%Y-%m-%d %H:%M:%S.%f')
        if text.endswith('Z'):
            return t.replace(tzinfo=UTC), 'measured'
        return t.astimezone(UTC), 'measured'
    except (OSError, ValueError):
        t = datetime.datetime.strptime(re.search(r'(\d{8}-\d{6})', stem).group(1), '%Y%m%d-%H%M%S')
        # older captures: the decoder started ~10-45 s after the file was named
        return t.astimezone(UTC) + datetime.timedelta(seconds=30), 'estimated'


def read_frames(lines):
    frames = []
    for line in lines:
        f = line.split()
        if len(f) < 6 or not f[0].endswith(':') or f[0] == 'ERR:':
            continue
        try:
            snr = float(f[5].split('|')[2]) if '|' in f[5] else 0.0
            frames.append((f[0][:-1], float(f[2]), int(f[3]), snr, line))
        except ValueError:
            continue
    return frames


# ---------- naming satellites ----------
def match_position(sats, t, lat, lon):
    """the satellite closest to straight above (lat, lon) within +-90 s of t"""
    best = None
    for s in sats:
        if look(s, t, lat, lon)[0] < 15:
            continue
        for dt in range(-90, 91, 2):
            el = look(s, t + datetime.timedelta(seconds=dt), lat, lon)[0]
            if best is None or el > best[1]:
                best = (s.name, el, dt)
    return best


def match_doppler(sats, pts, shifts):
    """fit the measured frequencies of one plane/slot pass to each satellite's Doppler curve"""
    fits = []
    for s in sats:
        best = None
        for sh in shifts:
            looks = [look(s, t + datetime.timedelta(seconds=sh), *HOME) for t, _ in pts]
            if min(l[0] for l in looks) < 0:
                continue
            diff = [hz + l[3] / C * 1626.1e6 for (_, hz), l in zip(pts, looks)]
            off = sum(diff) / len(diff)
            rms = math.sqrt(sum((d - off) ** 2 for d in diff) / len(diff))
            if best is None or rms < best:
                best = rms
        if best is not None:
            fits.append((best, s.name))
    fits.sort()
    return fits


def match_ecef(sats, t, xyz_km, max_km=300.0):
    """(name, km off) of the catalog satellite nearest an Earth-fixed position the satellite reported
    itself, or (None, km) if none is within max_km"""
    best = None
    for s in sats:
        p = s.ecef(t)
        d = math.sqrt(sum((a - b) ** 2 for a, b in zip(p, xyz_km)))
        if best is None or d < best[1]:
            best = (s.name, d)
    if best is None or best[1] > max_km:
        return None, best[1] if best else None
    return best


def ecef_to_lla(x, y, z):
    """Earth-fixed km -> latitude, longitude (degrees), height (km) on WGS84"""
    f = 1 / 298.257223563
    e2 = f * (2 - f)
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1 - e2))
    for _ in range(5):
        N = RE / math.sqrt(1 - e2 * math.sin(lat) ** 2)
        h = p / math.cos(lat) - N
        lat = math.atan2(z, p * (1 - e2 * N / (N + h)))
    N = RE / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    return math.degrees(lat), math.degrees(lon), p / math.cos(lat) - N


def read_orbcomm(lines):
    """Orbcomm receiver lines -> [(kind, unix t, Hz, sat id or None, snr, rest)]
       OEP: t Hz satid snr x y z vx vy vz   ephemeris: the satellite's own position (km) and velocity (km/s)
       OSY: t Hz satid snr                  sync packet
       OPK: t Hz satid snr type             any other packet
       OBAD: t Hz - snr                     heard, but failed its checksum"""
    out = []
    for line in lines:
        f = line.split()
        if len(f) < 5 or not f[0].endswith(':') or not f[0].startswith('O'):
            continue
        try:
            sid = None if f[3] == '-' else int(f[3], 16)
            out.append((f[0][:-1], float(f[1]), float(f[2]), sid, float(f[4]), f[5:]))
        except ValueError:
            continue
    return out


def analyze_orbcomm(stem, sats):
    """an Orbcomm capture: every position fix (named by matching it to the catalog) and, per satellite and
    pass, the measured carrier frequency of all its packets (its Doppler curve)"""
    start, how = capture_start(stem)
    src = frames_source(stem)
    if src.endswith(storage.REDUCED):
        meta, lines = storage.read_reduced(src)
    else:
        meta, lines = {}, open(src, errors='replace').read().splitlines()
    pk = read_orbcomm(lines)
    iso = lambda t: datetime.datetime.fromtimestamp(t, UTC).isoformat(timespec='seconds')
    heard = [p[1] for p in pk if p[0] != 'OEP'] or [p[1] for p in pk]     # reception times (ephemeris times are the satellite's)
    span = (max(heard) - min(heard)) / 3600 if heard else 0
    cap = {'capture': os.path.basename(stem), 'system': 'orbcomm', 'start_utc': iso(start.timestamp()), 'clock': how,
           'hours': meta.get('hours', round(span, 3)), 'counts': meta.get('counts', dict(collections.Counter(p[0] for p in pk))),
           'slots': [], 'positions': [], 'passes': []}
    ids = collections.defaultdict(collections.Counter)      # satellite ID -> names its positions matched
    for kind, t, hz, sid, snr, rest in pk:
        if kind != 'OEP' or len(rest) < 6:
            continue
        xyz = [float(v) for v in rest[:3]]
        vel = [float(v) for v in rest[3:6]]
        name, off = match_ecef(sats, datetime.datetime.fromtimestamp(t, UTC), xyz)
        if name:
            ids[sid][name] += 1
        lat, lon, alt = ecef_to_lla(*xyz)
        cap['positions'].append({'time': iso(t), 't': round(t, 3), 'net_id': sid, 'beam': None, 'lat': round(lat, 4),
                                 'lon': round(lon, 4), 'alt_km': round(alt), 'snr': snr, 'name': name, 'hz': hz,
                                 'catalog_off_km': round(off, 1) if off is not None else None,
                                 'xyz_km': [round(v, 3) for v in xyz], 'vel_kms': [round(v, 4) for v in vel]})
    id_name = {sid: c.most_common(1)[0][0] for sid, c in ids.items()}
    for p in cap['positions']:          # a position that matched nothing still belongs to its satellite ID
        p['name'] = p['name'] or id_name.get(p['net_id'])
    by_sat = collections.defaultdict(list)
    for kind, t, hz, sid, snr, rest in pk:
        if kind != 'OBAD' and sid in id_name:
            by_sat[id_name[sid], sid].append((t, hz, snr))
    for (name, sid), hits in sorted(by_sat.items()):
        hits.sort()
        groups, cur = [], [hits[0]]
        for h in hits[1:]:      # separate passes at gaps over 15 minutes
            if h[0] - cur[-1][0] > 900:
                groups.append(cur)
                cur = []
            cur.append(h)
        groups.append(cur)
        for g in groups:
            cap['passes'].append({'name': name, 'sat_id': sid, 'first': iso(g[0][0]), 'last': iso(g[-1][0]),
                                  'frames': len(g), 'best_snr': max(h[2] for h in g),
                                  'doppler': [[round(h[0], 3), h[1]] for h in g]})
    return cap


def norad_names(system='ham'):
    """NORAD number -> name in the public catalog's list (ham satellites are named the catalog's way,
    e.g. gr-satellites' 'AO-73' -> 'FUNCUBE-1 (AO-73)', so the info window can find them)"""
    out = {}
    try:
        lines = open(tle_file(system)).read().splitlines()
    except OSError:
        return out
    for i in range(len(lines) - 2):
        if lines[i + 1].startswith('1 ') and lines[i + 2].startswith('2 ') and not lines[i].startswith(('1 ', '2 ')):
            try:
                out[int(lines[i + 1][2:7])] = lines[i].strip()
            except ValueError:
                pass
    return out


def read_ham(lines):
    """ham receiver lines -> [(kind, unix t, Hz, norad, name or None, snr, callsign, candidates)]
       HFR: t Hz norad name snr bytes callsign candidates   a decoded frame (name '?' = which satellite is unknown)
       HSIG: t Hz - snr seconds                             a signal heard but not decoded"""
    out = []
    for line in lines:
        f = line.split()
        if len(f) < 5 or f[0] not in ('HFR:', 'HSIG:'):
            continue
        try:
            if f[0] == 'HFR:' and len(f) >= 8:
                name = None if f[4] == '?' else f[4].replace('_', ' ')
                out.append(('HFR', float(f[1]), float(f[2]), int(f[3]), name, float(f[5]), f[7],
                            f[8].replace('_', ' ').split('|') if len(f) > 8 else []))
            elif f[0] == 'HSIG:':
                out.append(('HSIG', float(f[1]), float(f[2]), 0, None, float(f[4]), '-', []))
        except ValueError:
            continue
    return out


def ham_name(norad, name, names):
    """the catalog's name for a decoded satellite (falls back to gr-satellites' name)"""
    return names.get(norad) or name


def analyze_ham(stem, names):
    """a 2 m / 70 cm capture: per identified satellite and pass, its frames' measured frequencies (its Doppler
    curve) and callsigns; frames that can't be tied to one satellite are only counted"""
    start, how = capture_start(stem)
    src = frames_source(stem)
    if src.endswith(storage.REDUCED):
        meta, lines = storage.read_reduced(src)
    else:
        meta, lines = {}, open(src, errors='replace').read().splitlines()
    fr = read_ham(lines)
    iso = lambda t: datetime.datetime.fromtimestamp(t, UTC).isoformat(timespec='seconds')
    span = (max(f[1] for f in fr) - min(f[1] for f in fr)) / 3600 if fr else 0
    key = 'ham70' if os.path.basename(stem).startswith('ham70') else 'ham2m'
    cap = {'capture': os.path.basename(stem), 'system': key, 'start_utc': iso(start.timestamp()), 'clock': how,
           'hours': meta.get('hours', round(span, 3)),
           'counts': meta.get('counts', dict(collections.Counter(f[0] for f in fr))),
           'slots': [], 'positions': [], 'passes': [],
           'unidentified': sum(1 for f in fr if f[0] == 'HFR' and not f[4])}
    by_sat = collections.defaultdict(list)
    for kind, t, hz, norad, name, snr, call, _ in fr:
        if kind == 'HFR' and name:
            by_sat[ham_name(norad, name, names), norad].append((t, hz, snr, call))
    for (name, norad), hits in sorted(by_sat.items()):
        hits.sort()
        groups, cur = [], [hits[0]]
        for h in hits[1:]:      # separate passes at gaps over 15 minutes
            if h[0] - cur[-1][0] > 900:
                groups.append(cur)
                cur = []
            cur.append(h)
        groups.append(cur)
        for g in groups:
            cap['passes'].append({'name': name, 'norad': norad, 'first': iso(g[0][0]), 'last': iso(g[-1][0]),
                                  'frames': len(g), 'best_snr': max(h[2] for h in g),
                                  'callsigns': sorted({h[3] for h in g if h[3] != '-'}),
                                  'doppler': [[round(h[0], 3), h[1]] for h in g]})
    return cap


GNSS_NAMES = {'G': 'GPS', 'E': 'Galileo', 'C': 'BeiDou', 'J': 'QZSS', 'R': 'GLONASS'}


def gnss_name(sat):
    """'G05' -> 'GPS G05', 'E11' -> 'Galileo E11', 'C23' -> 'BeiDou C23', 'S31' -> 'WAAS 131'"""
    if sat[:1] == 'S':
        return 'WAAS %d' % (100 + int(sat[1:]))
    return '%s %s' % (GNSS_NAMES.get(sat[:1], 'GNSS'), sat)


def gnss_catalog_names():
    """'G22' / 'C06' -> the public catalog's name ('GPS BIIR-5  (PRN 22)', 'BEIDOU-2 IGSO-1 (C06)'), for the info
    window. (Galileo's catalog names don't carry the PRN; those are matched by position instead.)"""
    out = {}
    try:
        lines = open(tle_file('gnss')).read().splitlines()
    except OSError:
        return out
    for i in range(len(lines) - 2):
        n = lines[i].strip()
        if not lines[i + 1].startswith('1 ') or '(' not in n:
            continue
        tag = n[n.rindex('(') + 1:].rstrip(')').strip()
        if tag.startswith('PRN ') and n.startswith('GPS '):
            out['G%02d' % int(tag[4:])] = n
        elif '/PRN ' in tag and 120 <= int(tag.split('/PRN ')[1]) <= 158:     # SBAS: 'SES-15 (WAAS/PRN 133)' -> S33
            out['S%02d' % (int(tag.split('/PRN ')[1]) - 100)] = n
        elif tag[:1] == 'C' and tag[1:].isdigit():
            out['C%02d' % int(tag[1:])] = n
    return out


def read_gnss(lines):
    """GNSS receiver lines -> [(kind, unix t, fields)]
       GCH: t sat sig cn0 doppler_hz lock_s    GEPH: t sat <$EPH fields>    GPOS: t lat lon height nsat
       GSPP: t lat lon height nsat hdop ...    (our own WAAS-corrected, carrier-smoothed fix: spp.py)
       GGAP: t gap_ms source                   (a hole in the sample stream the receiver zero-filled: usb, queue)
       GLOG: t EVENT [sat] [values...]         (receiver events: FOUND LOST FALSELOCK EPOCH RESTART START)"""
    out = []
    for line in lines:
        f = line.split(None, 3)
        if len(f) < 3 or f[0] not in ('GCH:', 'GEPH:', 'GPOS:', 'GSPP:', 'GGAP:', 'GLOG:', 'GALM:', 'GHEARD:'):
            continue
        try:
            out.append((f[0][:-1], float(f[1]), line.split()[2:]))
        except ValueError:
            continue
    return out


def _eph_key(text):
    import gnss
    try:
        return gnss.archive_key(text)
    except (ValueError, IndexError):
        return text                 # (an odd line: kept, keyed by itself)


def analyze_gnss(stem):
    """a GNSS capture: per satellite and pass its Doppler curve and signal level, the newest broadcast
    ephemeris of each satellite ('ephemerides') and every different one ('eph_archive'), the newest almanac of each
    satellite ('almanacs') and the UTC / GGTO parameters ('nav_params'), and the position fixes"""
    start, how = capture_start(stem)
    src = frames_source(stem)
    if src.endswith(storage.REDUCED):
        meta, lines = storage.read_reduced(src)
    else:
        meta, lines = {}, open(src, errors='replace').read().splitlines()
    rows = read_gnss(lines)
    iso = lambda t: datetime.datetime.fromtimestamp(t, UTC).isoformat(timespec='seconds')
    ts = [t for _, t, _ in rows]
    cap = {'capture': os.path.basename(stem), 'system': systems.of_capture(os.path.basename(stem))['key'],
           'start_utc': iso(start.timestamp()), 'clock': how,
           'hours': meta.get('hours', round((max(ts) - min(ts)) / 3600, 3) if ts else 0),
           'counts': meta.get('counts', dict(collections.Counter(k for k, _, _ in rows))),
           'slots': [], 'positions': [], 'passes': [], 'ephemerides': {}, 'gnss_fixes': []}
    by_sat = collections.defaultdict(list)
    bad = restarts = 0
    gaps = []
    heard = collections.Counter()   # satellite -> GHEARD lines (one per 60 s search that found it)
    archive, alm, params = {}, {}, {}
    for kind, t, f in rows:
        try:        # one line cut short by a crash (a field missing or half written) must not lose the capture
            if kind == 'GCH' and len(f) >= 4:
                by_sat[f[0]].append((t, GNSS_L1 + float(f[3]), float(f[2])))
            elif kind == 'GEPH' and len(f) >= 2:
                text = ' '.join(f[1:])
                cap['ephemerides'][f[0]] = [round(t, 3), text]     # the newest of each (what older readers use)
                # ...and every different one (by IODE + reference time), so the tracker can draw a satellite at an
                # earlier time with the orbit it was sending then (~18 kB a day; the later copy of one is kept)
                archive.setdefault(f[0], {})[_eph_key(text)] = [round(t, 3), text]
            elif kind == 'GALM' and len(f) >= 2:
                # an almanac entry (the receiver sends each satellite's only when its reference time is newer) or the
                # UTC / GGTO parameters: the last of each per capture
                (params if f[0] in ('UTCG', 'UTCE', 'GGTO') else alm)[f[0]] = [round(t, 3), ' '.join(f[1:])]
            elif kind == 'GPOS' and len(f) >= 4:
                cap['gnss_fixes'].append(gnss_lines.gpos_fix(round(t, 3), f))
            elif kind == 'GSPP' and len(f) >= 5:
                fx = gnss_lines.gspp_fix(round(t, 3), f)        # (quick / prior flags: see gnss_lines)
                fx.pop('waas', None)                            # (not kept: the session stats have their own)
                cap.setdefault('gnss_spp', []).append(fx)
            elif kind == 'GGAP' and len(f) >= 1:
                gaps.append([round(t, 3), float(f[0])])
            elif kind == 'GLOG' and len(f) >= 1 and f[0] == 'RESTART':
                restarts += 1
            elif kind == 'GHEARD' and len(f) >= 1:     # heard by our own search while pocket_trk was not tracking it
                heard[f[0]] += 1
        except (ValueError, IndexError):
            bad += 1
    if bad:
        cap['counts']['bad'] = bad
    if gaps:
        cap['gaps'] = gaps
    if restarts:
        cap['restarts'] = restarts
    if heard:
        cap['heard_not_tracked'] = dict(heard)
    if archive:
        cap['eph_archive'] = {sat: sorted(v.values()) for sat, v in sorted(archive.items())}
    if alm:
        cap['almanacs'] = dict(sorted(alm.items()))
    if params:
        cap['nav_params'] = params
    cap['_cn0'] = [h[2] for v in by_sat.values() for h in v]
    for sat, hits in sorted(by_sat.items()):
        hits.sort()
        groups, cur = [], [hits[0]]
        for h in hits[1:]:
            if h[0] - cur[-1][0] > 900:
                groups.append(cur)
                cur = []
            cur.append(h)
        groups.append(cur)
        for g in groups:
            ps = {'name': gnss_name(sat), 'sat': sat, 'first': iso(g[0][0]), 'last': iso(g[-1][0]),
                  'frames': len(g), 'best_snr': max(h[2] for h in g)}
            ps.update(slim_pass(g))
            cap['passes'].append(ps)
    cn0 = cap.pop('_cn0')
    try:
        cap['session'] = session_stats(cap, cn0)
    except Exception as e:      # the statistics are a summary; losing them must not lose the capture
        cap['session'], cap['session_error'] = None, str(e)
        print('%s: session statistics failed (%s)' % (cap['capture'], e))
    garmin = analyze_garmin(stem, rows, cap)    # (a Garmin handheld recorded beside it: receiver.py's sidecar)
    if garmin:
        cap['garmin'] = garmin
    return cap


def analyze_garmin(stem, rows, cap):
    """the Garmin handheld recorded beside a GNSS capture (<capture>.garmin, written by garmin_usb.py while it
    received) against our receiver, second by second and satellite by satellite (garmin_usb.compare); None without
    one. Like the other raw files the .garmin goes after a week (storage.RAW_EXT) - a --redo then keeps the
    comparison made while it was there. A comparison that fails is noted, never the capture lost."""
    path = stem + '.garmin'
    if not os.path.exists(path):
        try:
            with open(_cap_file(cap['capture'])) as f:
                return json.load(f).get('garmin')
        except (OSError, ValueError, AttributeError):
            return None
    try:
        import garmin_usb
        g, o = garmin_usb.LogReader(), garmin_usb.Ours()
        with open(path, errors='replace') as f:
            g.feed_text(f.read(), final=True)
        o.feed_rows(rows)
        return garmin_usb.compare(g, o, cap.get('system'), HOME)
    except Exception as e:
        log_error('comparing %s with the Garmin' % cap['capture'], e)
        return {'error': '%s: %s' % (type(e).__name__, e)}


def session_stats(cap, cn0=None):
    """a GNSS capture in a few numbers, kept with it (like a handheld's track statistics): when and how long, how
    long the first fixes took (pocket_trk's and our own solver's), satellites and their signal, and the session's
    middle position with its accuracy figures (posquality)"""
    import posquality
    start = datetime.datetime.fromisoformat(cap['start_utc']).timestamp()
    fx, sp = cap.get('gnss_fixes', []), cap.get('gnss_spp', [])
    good = [f for f in sp if not f.get('quick')]
    s = {'start': start, 'hours': cap.get('hours', 0),
         'first_fix_s': round(min(f['t'] for f in fx) - start, 1) if fx else None,
         'first_waas_s': round(min(f['t'] for f in sp) - start, 1) if sp else None,
         'fixes': len(fx), 'waas_fixes': len(good),
         'satellites': len({p['sat'] for p in cap.get('passes', [])})}
    if cn0:
        s['cn0_median'] = round(float(np.median(cn0)), 1)
    if cap.get('gaps'):             # holes in the sample stream (dropped USB buffers) and receiver restarts
        s['gaps'] = len(cap['gaps'])
    if cap.get('restarts'):
        s['restarts'] = cap['restarts']
    for key, pts in (('pocket', fx), ('waas', good)):
        q = posquality.quality(pts, home=HOME) if pts else None
        if q and q.get('enough'):
            s[key] = {k: round(q[k], 8 if k in ('lat', 'lon') else 2) for k in
                      ('lat', 'lon', 'height_m', 'r95_m', 'cep50_m', 'avg95_m', 'span_s', 'from_home_m') if k in q}
    if cap.get('heard_not_tracked'):
        s['heard_not_tracked'] = cap['heard_not_tracked']
    try:        # the sky's signal map and lock table (the sky view): a summary - failing must not lose the rest
        sky = sky_stats(cap)
        if sky:
            s.update(sky)
    except Exception as e:
        s['sky_error'] = '%s: %s' % (type(e).__name__, e)
    return s


# ---------- GNSS passes kept compact, and what they say about the sky (2026-11-07)
# A GNSS pass has a GCH line every second for hours. The catalog kept each second's Doppler - 21 of the 26 MB of the
# GNSS catalog on 2026-11-07 - although no GNSS orbit is fitted from it (their orbits are the satellites' own broadcast
# ones), and threw each second's C/N0 away. Now a pass keeps its Doppler every DOPPLER_STEP_S ('doppler_n': how many
# seconds of tracking each point stands for, so counts of 'seconds heard' stay right) and its C/N0 every second in
# CN0_STEP steps, packed (pack_cn0 / cn0_series). Passes catalogued before keep their 1 Hz Doppler and have no C/N0
# (one 'catalog.py --redo' converts them while the capture's frames are still there).
DOPPLER_STEP_S = 10
CN0_STEP = 0.5                  # dB-Hz (pocket_trk prints 0.1, but its own C/N0 wanders by ~1 dB second to second)
SCINT_WINDOW_S = 60             # the S4 index is defined over 60 s
SCINT_MIN_S = 50                # ...a window needs this many seconds with C/N0
SKY_CELL_DEG = 5                # the signal map's cells (elevation and azimuth)
LOCK_BANDS = (10, 20, 30, 45, 60, 90)   # elevation bands of the lock table (deg; below 10 the solver ignores them)
LOCK_STEP_S = 30                # the lock table looks at every satellite this often
LOCK_NEAR_S = 15                # ...and counts it locked with a GCH line within this of the moment
ORBIT_OK_KM = 300.0             # an orbit placed worse than this (1 sigma) is not used: ~1 deg seen from here


def slim_pass(g):
    """a GNSS pass's hits [(t, Hz, C/N0)] (sorted) -> its compact fields: 'doppler' [[t, Hz]] every DOPPLER_STEP_S (the
    first line of each step), 'doppler_n' (seconds each point stands for), 'cn0' (every second, pack_cn0) and 'scint'
    (scint_index, when the pass is long enough)"""
    dop, dn, step = [], [], None
    for t, hz, _ in g:
        k = int(t // DOPPLER_STEP_S)
        if k != step:
            step = k
            dop.append([round(t, 3), round(hz, 1)])
            dn.append(0)
        dn[-1] += 1
    out = {'doppler': dop, 'doppler_n': dn, 'cn0': pack_cn0([(h[0], h[2]) for h in g])}
    sc = scint_index(np.array([h[0] for h in g], float), np.array([h[2] for h in g], float))
    if sc:
        out['scint'] = sc
    return out


def pack_cn0(hits):
    """[(t, C/N0)] sorted by time -> {'t0': first t, 'v': [...]}: one entry per second after t0, round(C/N0 / CN0_STEP),
    and a run of k seconds without a line as -k (a pass may have holes of up to 15 minutes). 3-4 characters a second
    in the catalog, where each second's Doppler took 22"""
    if not hits:
        return {'t0': 0, 'v': []}
    t0 = hits[0][0]
    v, last = [], -1
    for t, c in hits:
        k = int(round(t - t0))
        if k <= last:
            continue                # (two lines in one second: the first counts)
        if k > last + 1:
            v.append(-(k - last - 1))
        v.append(max(1, int(round(c / CN0_STEP))))
        last = k
    return {'t0': round(t0, 3), 'v': v}


def cn0_series(ps):
    """a GNSS pass's C/N0 every second as numpy arrays (unix times, dB-Hz); empty for a pass catalogued before
    C/N0 was kept"""
    c = ps.get('cn0') if isinstance(ps, dict) else None
    if not c or not c.get('v'):
        return np.zeros(0), np.zeros(0)
    v = np.asarray(c['v'], float)
    gap = v < 0
    step = np.where(gap, -v, 1.0)
    pos = np.cumsum(step) - step        # seconds after t0 each entry starts at
    return float(c['t0']) + pos[~gap], v[~gap] * CN0_STEP


def scint_index(ts, cn0):
    """an S4-like scintillation index from 1 Hz C/N0: per 60 s window (with 50+ seconds of C/N0), the standard deviation
    of the signal power (C/N0 as linear intensity) after dividing out its straight-line trend (the slow rise and fall
    with elevation and the antenna pattern). The true S4 uses 50 Hz power; pocket_trk's C/N0 is already averaged
    over about a second, so fast scintillation is smoothed out, and its own estimation noise sets a floor that grows
    as the signal gets weaker (see the measured values in the satellite info window's notes). Returns
    {'s4': median over the windows, 's4_max', 'n': windows} or None"""
    if len(ts) < SCINT_MIN_S:
        return None
    w = ((ts - ts[0]) // SCINT_WINDOW_S).astype(int)
    edges = np.flatnonzero(np.diff(w)) + 1
    vals = []
    for idx in np.split(np.arange(len(ts)), edges):
        if len(idx) < SCINT_MIN_S:
            continue
        tt = ts[idx] - ts[idx].mean()
        ii = 10.0 ** (cn0[idx] / 10.0)
        a, b = np.polyfit(tt, ii, 1)
        trend = a * tt + b
        if trend.min() <= 0:
            continue
        vals.append(float(np.std(ii / trend)))
    if not vals:
        return None
    return {'s4': round(float(np.median(vals)), 3), 's4_max': round(max(vals), 3), 'n': len(vals)}


def _sat_orbits(cap):
    """satellite -> an orbit object (ecef(t) in km, sigma_km(t)) from the broadcast ephemerides OUR receiver decoded:
    this capture's newest of each satellite, else the one the solver kept (gnss-solver.json, any session). Never a
    public orbit: the sky map and lock table place satellites from what they sent us"""
    import gnss
    pool = dict(cap.get('ephemerides', {}))
    try:
        with open(os.path.join(paths.DATA, 'gnss-solver.json')) as f:
            for sat, (t, text) in json.load(f).get('ephs', {}).items():
                pool.setdefault(sat, [t, text])
    except (OSError, ValueError, TypeError, AttributeError):
        pass
    out = {}
    for sat, (t, fields) in pool.items():
        try:
            fields = fields.replace(' ', ',')
            if fields.startswith('SBAS'):
                o = gnss.SBASOrbit.from_fields(fields, sat, t)
            elif fields.startswith('GLO'):
                o = gnss.GlonassOrbit.from_fields(fields, sat, t)
            else:
                e = gnss.parse_eph('$EPH,%.3f,%s,%s' % (t, sat, fields))
                if not e or e.get('svh', 0):        # (an unhealthy satellite may not be sending at all)
                    continue
                o = gnss.EphemerisOrbit(e)
            out[sat] = o
        except Exception:
            continue
    return out


def _el_az(xyz_km):
    """(n, 3) Earth-fixed positions (km) -> elevation, azimuth (deg) from home"""
    lat, lon = HOME[0] * D2R, HOME[1] * D2R
    d = np.asarray(xyz_km, float).reshape(-1, 3) - np.array(site(HOME[0], HOME[1], HOME[2]))
    e = -math.sin(lon) * d[:, 0] + math.cos(lon) * d[:, 1]
    n = -math.sin(lat) * math.cos(lon) * d[:, 0] - math.sin(lat) * math.sin(lon) * d[:, 1] + math.cos(lat) * d[:, 2]
    u = math.cos(lat) * math.cos(lon) * d[:, 0] + math.cos(lat) * math.sin(lon) * d[:, 1] + math.sin(lat) * d[:, 2]
    return np.degrees(np.arctan2(u, np.hypot(e, n))), np.degrees(np.arctan2(e, n)) % 360


def sky_stats(cap):
    """what a GNSS session says about the sky here, each satellite placed by its own broadcast orbit:
      skymap  {'cell_deg', 'cells': [[el_row, az_col, seconds, median C/N0 / CN0_STEP]]} - how strong the signals
              were in each 5 x 5 deg patch of sky (which directions are weak: antenna placement)
      lock    {'step_s', 'bands', 'up', 'trk', 'no_orbit'} - GPS satellites only (whose orbits we nearly all have):
              per elevation band x 8 compass sectors (N, NE, ...), how many LOCK_STEP_S moments a satellite was up
              there and at how many of them it was tracked (the review's sky_cover.py table, per session)"""
    if HOME is None:
        return None
    orbits = _sat_orbits(cap)
    series = collections.defaultdict(list)
    for ps in cap.get('passes', []):
        ts, c = cn0_series(ps)
        if len(ts):
            series[ps['sat']].append((ts, c))
    if not series:
        return None
    n_az = 360 // SKY_CELL_DEG
    cells, vals = [], []
    for sat, parts in series.items():
        o = orbits.get(sat)
        if o is None:
            continue
        ts = np.concatenate([p[0] for p in parts])
        c = np.concatenate([p[1] for p in parts])
        good = np.asarray(o.sigma_km(ts)) < ORBIT_OK_KM
        if not good.any():
            continue
        ts, c = ts[good], c[good]
        el, az = _el_az(o.ecef(ts))
        up = el >= 0
        cells.append((np.minimum(el[up], 89.999) // SKY_CELL_DEG).astype(int) * n_az +
                     (az[up] // SKY_CELL_DEG).astype(int) % n_az)
        vals.append(c[up])
    out = {}
    if cells:
        cell, val = np.concatenate(cells), np.concatenate(vals)
        order = np.lexsort((val, cell))
        cell, val = cell[order], val[order]
        cut = np.flatnonzero(np.diff(cell)) + 1
        rows = []
        for a, b in zip(np.r_[0, cut], np.r_[cut, len(cell)]):
            if b <= a:
                continue
            k = int(cell[a])
            rows.append([k // n_az, k % n_az, int(b - a), int(round(float(np.median(val[a:b])) / CN0_STEP))])
        out['skymap'] = {'cell_deg': SKY_CELL_DEG, 'cells': rows}
    # the lock table: the GPS satellites we have an orbit of our own for, every LOCK_STEP_S over the session
    t_lo = min(float(p[0][0]) for parts in series.values() for p in parts)
    t_hi = max(float(p[0][-1]) for parts in series.values() for p in parts)
    if cap.get('system') == 'gnss1' and t_hi - t_lo >= 300:
        grid = np.arange(t_lo, t_hi, float(LOCK_STEP_S))
        nb = len(LOCK_BANDS) - 1
        up_n, trk_n = np.zeros((nb, 8), int), np.zeros((nb, 8), int)
        no_orbit = 0
        for k in range(1, 33):
            sat = 'G%02d' % k
            o = orbits.get(sat)
            good = np.asarray(o.sigma_km(grid)) < ORBIT_OK_KM if o is not None else np.zeros(len(grid), bool)
            if not good.any():
                no_orbit += 1
                continue
            g = grid[good]
            el, az = _el_az(o.ecef(g))
            band = np.searchsorted(LOCK_BANDS, el, side='right') - 1
            ok = (band >= 0) & (band < nb)
            if not ok.any():
                continue
            sector = (((az + 22.5) % 360) // 45).astype(int) % 8
            heard_ts = np.sort(np.concatenate([p[0] for p in series[sat]])) if sat in series else np.zeros(0)
            if len(heard_ts):
                i = np.minimum(np.searchsorted(heard_ts, g - LOCK_NEAR_S), len(heard_ts) - 1)
                near = (heard_ts[i] >= g - LOCK_NEAR_S) & (heard_ts[i] < g + LOCK_NEAR_S)
            else:
                near = np.zeros(len(g), bool)
            np.add.at(up_n, (band[ok], sector[ok]), 1)
            np.add.at(trk_n, (band[ok & near], sector[ok & near]), 1)
        out['lock'] = {'step_s': LOCK_STEP_S, 'bands': list(LOCK_BANDS), 'up': up_n.tolist(), 'trk': trk_n.tolist(),
                       'no_orbit': no_orbit}
    return out


GNSS_L1 = 1575.42e6


def frames_source(stem):
    """the decoded frames of a capture: full (.parsed / .live) or the reduced .frames.gz"""
    for ext in ('.parsed', '.live', storage.REDUCED):
        if os.path.exists(stem + ext) and os.path.getsize(stem + ext) > 0:
            return stem + ext
    return None


def analyze(stem, sats):
    start, how = capture_start(stem)
    src = frames_source(stem)
    if src.endswith(storage.REDUCED):       # only ring alerts and time & location frames were kept
        meta, lines = storage.read_reduced(src)
    else:
        meta, lines = {}, open(src, errors='replace')
    frames = read_frames(lines)
    at = lambda ms: start + datetime.timedelta(milliseconds=ms)
    iso = lambda t: t.isoformat(timespec='seconds')
    cap = {'capture': os.path.basename(stem), 'system': 'iridium', 'start_utc': iso(start), 'clock': how,
           'hours': meta.get('hours', round(max([f[1] for f in frames] or [0]) / 3.6e6, 2)),
           'counts': meta.get('counts', dict(collections.Counter(f[0] for f in frames))), 'slots': [], 'positions': []}
    clock_errs = []

    for kind, ms, hz, snr, line in frames:
        m = re.search(r'sat:(\d+) beam:(\d+).*?pos=\(([-+\d.]+)/([-+\d.]+)\) alt=(-?\d+)', line)
        if kind != 'IRA' or not m or int(m.group(5)) < 100:
            continue
        lat, lon = float(m.group(3)), float(m.group(4))
        hit = match_position(sats, at(ms), lat, lon)
        name = hit[0] if hit and hit[1] > 85 else None
        if name:
            clock_errs.append(hit[2])
        x = re.search(r'xyz=\(([-+\d]+),([-+\d]+),([-+\d]+)\)', line)
        cap['positions'].append({'time': iso(at(ms)), 'ms': ms, 'net_id': int(m.group(1)), 'beam': int(m.group(2)),
                                 'lat': lat, 'lon': lon, 'alt_km': int(m.group(5)), 'snr': snr, 'name': name,
                                 'overhead_deg': round(hit[1], 1) if hit else None, 'hz': hz,
                                 # Earth-fixed position as sent by the satellite (4 km units -> km)
                                 'xyz_km': [4 * int(v) for v in x.groups()] if x else None})

    # a capture's measured clock error (from position fixes) also corrects its Doppler fits
    shift = sorted(clock_errs)[len(clock_errs) // 2] if clock_errs else 0
    shifts = [shift] if clock_errs or how == 'measured' else range(-20, 21, 5)

    passes = collections.defaultdict(list)
    for kind, ms, hz, snr, line in frames:
        m = re.search(r' P(\d) S(\d+) ', line)
        if kind == 'ITL' and m:
            passes[(int(m.group(1)), int(m.group(2)))].append((at(ms), hz, snr))
    for (p, s), hits in sorted(passes.items()):
        hits.sort()
        groups, cur = [], [hits[0]]
        for h in hits[1:]:   # split into separate passes at gaps over 15 minutes
            if (h[0] - cur[-1][0]).total_seconds() > 900:
                groups.append(cur)
                cur = []
            cur.append(h)
        groups.append(cur)
        for g in groups:
            slot = {'plane': p, 'slot': s, 'first': iso(g[0][0]), 'last': iso(g[-1][0]), 'frames': len(g),
                    'best_snr': max(h[2] for h in g), 'name': None, 'confidence': None, 'rms_hz': None,
                    # measured frequency of every frame: (unix time, Hz) - the Doppler curve used for orbit fitting
                    'doppler': [[round(h[0].timestamp(), 3), h[1]] for h in g]}
            if len(g) >= 3 and (g[-1][0] - g[0][0]).total_seconds() >= 20:
                fits = match_doppler(sats, [(h[0], h[1]) for h in g], shifts)
                if fits:
                    rms, name = fits[0]
                    second = fits[1][0] if len(fits) > 1 else 1e9
                    slot['name'], slot['rms_hz'] = name, round(rms)
                    slot['confidence'] = 'confident' if rms < 300 and second > 2.5 * rms else 'possible'
            cap['slots'].append(slot)
    return cap


# ---------- report ----------
def local(iso_utc):
    return datetime.datetime.fromisoformat(iso_utc).astimezone()


def fmt(iso_utc):
    t = local(iso_utc)
    return t.strftime('%b %d %I:%M %p').replace(' 0', ' ')


def build_report(caps):
    out = []
    w = out.append
    decoded = sum(n for c in caps for k, n in c['counts'].items() if k != 'RAW')
    raw = sum(c['counts'].get('RAW', 0) for c in caps)
    hours = sum(c['hours'] for c in caps)

    slots = collections.defaultdict(list)          # (plane, slot) -> sightings
    for c in caps:
        for s in c['slots']:
            slots[(s['plane'], s['slot'])].append(s)
    names = collections.defaultdict(lambda: {'positions': [], 'slots': set(), 'net_ids': set(), 'seen': []})
    for c in caps:
        for p in c['positions']:
            if p['name']:
                n = names[p['name']]
                n['positions'].append(p)
                n['net_ids'].add(p['net_id'])
                n['seen'].append((p['time'], p['snr']))
        for s in c['slots']:
            if s['name'] and s['confidence'] == 'confident':
                n = names[s['name']]
                n['slots'].add((s['plane'], s['slot']))
                n['seen'].append((s['first'], s['best_snr']))
    slot_name = {}
    for nm, n in names.items():
        for ps in n['slots']:
            slot_name.setdefault(ps, set()).add(nm)

    w('IRIDIUM CATALOG  -  updated %s' % datetime.datetime.now().strftime('%b %d %Y %I:%M %p'))
    w('=' * 78)
    w('Captures: %d    listening time: %.1f h    decoded: %d    damaged (RAW): %d'
      % (len(caps), hours, decoded, raw))
    w('Plane/slot positions heard: %d of 66    satellites named: %d    position fixes: %d'
      % (len(slots), len(names), sum(len(c['positions']) for c in caps)))
    w('')
    w('PLANES AND SLOTS HEARD  (x = heard, N = heard and named by Doppler)')
    w('            ' + ''.join('S%-3d' % s for s in range(1, 12)))
    for p in range(1, 7):
        row = ''
        for s in range(1, 12):
            row += ('N   ' if (p, s) in slot_name else 'x   ') if (p, s) in slots else '.   '
        w('  Plane %d   %s  %d/11' % (p, row, sum(1 for s in range(1, 12) if (p, s) in slots)))
    w('')
    w('EVERY PLANE/SLOT HEARD')
    w('  %-12s %-6s %-7s %-17s %-17s %-8s %s' % ('plane/slot', 'passes', 'frames', 'first heard', 'last heard', 'best dB', 'name'))
    for (p, s), lst in sorted(slots.items()):
        nm = ', '.join(sorted(slot_name.get((p, s), []))) or '-'
        w('  P%d S%-9s %-6d %-7d %-17s %-17s %-8.1f %s'
          % (p, '%02d' % s, len(lst), sum(x['frames'] for x in lst), fmt(min(x['first'] for x in lst)),
             fmt(max(x['last'] for x in lst)), max(x['best_snr'] for x in lst), nm))
    w('')
    w('SATELLITES IDENTIFIED BY NAME')
    w('  %-13s %-10s %-11s %-9s %-17s %-17s %s' % ('name', 'network ID', 'plane/slot', 'pos fixes', 'first seen', 'last seen', 'best dB'))
    for nm in sorted(names, key=lambda x: int(re.sub(r'\D', '', x) or 0)):
        n = names[nm]
        w('  %-13s %-10s %-11s %-9d %-17s %-17s %.1f'
          % (nm, ','.join('%03d' % i for i in sorted(n['net_ids'])) or '-',
             ','.join('P%d S%02d' % ps for ps in sorted(n['slots'])) or '-', len(n['positions']),
             fmt(min(t for t, _ in n['seen'])), fmt(max(t for t, _ in n['seen'])), max(x for _, x in n['seen'])))
    w('  (network ID = the number the satellite uses for itself in ring alerts)')
    w('')
    w('POSITION FIXES  (ring alerts that reported where the satellite was)')
    for c in caps:
        for p in c['positions']:
            w('  %-17s net ID %03d  lat %+6.2f lon %+7.2f  %d km  %4.1f dB  -> %s'
              % (fmt(p['time']), p['net_id'], p['lat'], p['lon'], p['alt_km'], p['snr'],
                 p['name'] or 'no orbit match'))
    w('')
    w('DOPPLER MATCHES  (plane/slot passes with 3+ frames; confident = clear best fit)')
    for c in caps:
        for s in c['slots']:
            if s['name']:
                w('  %-17s P%d S%02d  %3d frames  -> %-13s %-9s fit %d Hz'
                  % (fmt(s['first']), s['plane'], s['slot'], s['frames'], s['name'], s['confidence'], s['rms_hz']))
    w('')
    w('CAPTURES')
    for c in caps:
        dec = sum(n for k, n in c['counts'].items() if k != 'RAW')
        w('  %-24s %-17s %5.2f h  %5d decoded  %5d RAW  %2d fixes  clock %s'
          % (c['capture'], fmt(c['start_utc']), c['hours'], dec, c['counts'].get('RAW', 0), len(c['positions']), c['clock']))
    return '\n'.join(out) + '\n'


def _load_legacy():
    """the old single-file catalog (catalog-data.json, or its previous copy if cut off)"""
    for path in (LEGACY, PREVIOUS):
        if not os.path.exists(path):
            continue
        try:
            with open(path) as f:
                data = json.load(f)
            if path != LEGACY:
                print('catalog-data.json could not be read - using the previous copy (%s)' % path)
            return data
        except ValueError:
            continue
    if os.path.exists(LEGACY):  # unreadable and no good copy: stop rather than start the catalog from nothing
        sys.exit('catalog-data.json can not be read and there is no good previous copy - see README (recovery)')
    return {}


def _cap_file(name):
    return os.path.join(CAT_DIR, name + '.json')


def _alive(pid):
    """is a process with this ID still running? (the same check receiver.py makes on a capture's .running marker)"""
    if os.name != 'nt':
        try:
            os.kill(pid, 0)
            return True
        except OSError:
            return False
    import ctypes
    h = ctypes.windll.kernel32.OpenProcess(0x1000, False, pid)        # PROCESS_QUERY_LIMITED_INFORMATION
    if not h:
        return False
    try:
        code = ctypes.c_ulong()
        return bool(ctypes.windll.kernel32.GetExitCodeProcess(h, ctypes.byref(code))) and code.value == 259  # STILL_ACTIVE
    finally:
        ctypes.windll.kernel32.CloseHandle(h)


def _lock_owner():
    """PID written in the lock file, or 0 if it can't be read"""
    try:
        return int(open(LOCK).read().split()[0])
    except (OSError, ValueError, IndexError):
        return 0


def acquire_lock(wait_s=LOCK_WAIT_S):
    """take the run lock: True when this process owns it. Another live run keeps us waiting (polling every 2 s, up
    to wait_s); a lock left by a run that died (its PID is gone) is removed. Creating the file with O_EXCL is what
    makes the hand-over safe: of two runs trying at the same moment only one succeeds"""
    os.makedirs(CAT_DIR, exist_ok=True)
    deadline = time.time() + wait_s
    told = False
    while True:
        try:
            fd = os.open(LOCK, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            pid = _lock_owner()
            if pid and _alive(pid) and pid != os.getpid():
                if time.time() >= deadline:
                    return False
                if not told:
                    print('another catalog run (PID %d) is busy - waiting for it' % pid)
                    told = True
                time.sleep(2)
                continue
            try:                # stale (its run is gone, or the file is unreadable): clear it and try again
                if _lock_owner() == pid:        # (unless a fresh run took it in the meantime)
                    os.remove(LOCK)
            except OSError:
                pass
            time.sleep(0.1)
            continue
        with os.fdopen(fd, 'w') as f:
            f.write('%d %s\n' % (os.getpid(), time.strftime('%Y-%m-%d %H:%M:%S')))
        return True


def release_lock():
    if _lock_owner() == os.getpid():
        try:
            os.remove(LOCK)
        except OSError:
            pass


def load_data():
    """every capture's results: one file each in data\\catalog (the index lists them). Before the first save in this
    layout, the old single catalog-data.json. A capture file that can't be read is left out (with a message) - the
    rest load; it is analysed again from its capture files at the next catalog run if they are still there."""
    if not os.path.exists(DATA):
        data = _load_legacy()
        _LOADED.clear()
        _LOADED.update(data)
        return data
    names = _index_names()
    data = {}
    for name in names:
        try:
            with open(_cap_file(name)) as f:
                data[name] = json.load(f)
        except (OSError, ValueError) as e:
            print('catalog: %s could not be read (%s)' % (name, e))
    _LOADED.clear()
    _LOADED.update(data)
    return data


def _index_names():
    """the capture names index.json lists now (the folder's files if the index can't be read)"""
    try:
        with open(DATA) as f:
            return list(json.load(f)['captures'])
    except (OSError, ValueError, KeyError, TypeError):
        return [os.path.basename(p)[:-5] for p in glob.glob(os.path.join(CAT_DIR, '*.json'))
                if os.path.basename(p) != 'index.json']


def _write_atomic(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, separators=(',', ':'))
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def save_data(data):
    """one small file per capture, each written to a temporary file and swapped in, only when it changed; then the
    index. An interrupted save can at most leave one capture's previous version (never a cut-off catalog), and a
    save rewrites only what is new instead of the whole history.
    Only captures this process loaded and then dropped are deleted - never "every file not in the dict": a second
    run that finished in between may have written entries we never saw, and the index written here keeps them"""
    os.makedirs(CAT_DIR, exist_ok=True)
    for name, cap in data.items():
        text = json.dumps(cap, separators=(',', ':'))
        path = _cap_file(name)
        try:
            with open(path) as f:
                if f.read() == text:
                    continue
        except OSError:
            pass
        _write_atomic(path, cap)
    dropped = _LOADED - set(data)           # captures the caller removed from the catalog
    for name in dropped:
        try:
            os.remove(_cap_file(name))
        except OSError:
            pass
    others = {n for n in _index_names() if n not in dropped and os.path.exists(_cap_file(n))} if os.path.exists(DATA) else set()
    _write_atomic(DATA, {'captures': sorted(set(data) | others), 'saved': time.time()})
    _LOADED.clear()
    _LOADED.update(data)


def log_error(what, exc):
    """a traceback into the tracker's error log (this program's own stderr is thrown away when the tracker runs it)"""
    text = ''.join(traceback.format_exception(type(exc), exc, exc.__traceback__))
    try:
        with open(paths.ERROR_LOG, 'a') as f:
            f.write('%s  catalog: %s\n%s\n' % (time.ctime(), what, text))
    except OSError:
        pass


def error_stub(name, stem, exc):
    """what the catalog keeps for a capture whose analysis failed: every key the tracker reads, no frames (so it
    is skipped like an empty run), and the error. --redo analyses it again; its raw files are kept meanwhile"""
    try:
        start, how = capture_start(stem)
        start_utc = start.isoformat(timespec='seconds')
    except Exception:
        start_utc, how = None, None
    return {'capture': name, 'system': systems.of_capture(name)['key'], 'start_utc': start_utc, 'clock': how,
            'hours': 0, 'counts': {}, 'slots': [], 'positions': [], 'passes': [], 'ephemerides': {}, 'gnss_fixes': [],
            'error': str(exc) or type(exc).__name__}


def main():
    if HOME is None:
        sys.exit('No antenna location yet - start the tracker once and enter it (it is saved in settings.json).')
    if not acquire_lock():
        sys.exit('another catalog run (PID %d) has been busy for over %d s - giving up; the next run catches up'
                 % (_lock_owner(), LOCK_WAIT_S))
    try:
        _main()
    finally:
        release_lock()


def _main():
    redo = '--redo' in sys.argv
    ready = {a[8:] for a in sys.argv if a.startswith('--ready=')}   # captures known to be finished
    data = load_data()      # (with --redo, captures whose files are gone keep their old results)
    sats = {s['key']: None for s in systems.ORDER}     # loaded when first needed
    recording = set()
    for name, files in sorted(storage.captures().items()):
        stem = os.path.join(paths.CAPTURES, name)
        bits = stem + '.bits'
        if os.path.exists(bits) and name not in ready and time.time() - os.path.getmtime(bits) < 20:
            recording.add(name)         # still being recorded
            continue
        if name in data and not redo:
            continue
        if systems.of_capture(name).get('system') in ('meteor', 'radio'):
            continue        # pictures, not frames: images.py makes them (receiver.py runs it after Stop)
        live = stem + '.live'
        if systems.of_capture(name) is not systems.IRIDIUM and name not in ready and os.path.exists(live) \
                and time.time() - os.path.getmtime(live) < 20:
            recording.add(name)
            continue
        if not frames_source(stem) and os.path.exists(bits) and os.path.getsize(bits) >= 200:
            subprocess.run([sys.executable, PARSER, '--no-stats', bits], cwd=paths.CAPTURES, capture_output=True)
        src = frames_source(stem)
        if not src:
            continue
        print('analyzing %s ...' % name)
        key = systems.of_capture(name)['key']
        try:        # one capture that can't be analysed must not stop the rest (or the save) - it used to
            if systems.get(key).get('system') == 'ham':
                load_tles(key)                  # refresh the name list
                cap = analyze_ham(stem, norad_names())
            elif systems.get(key).get('system') == 'gnss':
                cap = analyze_gnss(stem)        # (named from the satellite codes: no name list needed)
            else:
                if sats[key] is None:
                    sats[key] = load_tles(key)
                cap = analyze_orbcomm(stem, sats[key]) if key == 'orbcomm' else analyze(stem, sats[key])
            if not src.endswith(storage.REDUCED):
                storage.write_reduced(stem, src, cap['counts'], cap['hours'])
        except Exception as e:
            log_error('analysing %s' % name, e)
            print('%s could not be analysed (%s: %s) - kept for a later --redo; details in %s'
                  % (name, type(e).__name__, e, os.path.basename(paths.ERROR_LOG)))
            cap = error_stub(name, stem, e)
        data[name] = cap
    # session statistics for GNSS captures analysed before they existed (from what the catalog holds)
    for name, cap in data.items():
        if cap.get('system') == 'gnss1' and 'session' not in cap and cap.get('start_utc') and not cap.get('error'):
            try:
                cap['session'] = session_stats(cap)
            except Exception as e:
                log_error('session statistics of %s' % name, e)
    save_data(data)
    caps = [c for _, c in sorted(data.items()) if sum(c['counts'].values()) and c.get('system', 'iridium') == 'iridium']
    open(REPORT, 'w').write(build_report(caps))     # the text report covers Iridium
    print('catalog written to %s' % REPORT)
    # (a capture whose analysis failed keeps its raw files: they are all it has until a --redo succeeds)
    n, freed = storage.enforce({n for n, c in data.items() if not c.get('error')}, protect=recording)
    if n:
        print('deleted %d raw capture files (%s) - see storage settings' % (n, storage.fmt_bytes(freed)))


if __name__ == '__main__':
    main()

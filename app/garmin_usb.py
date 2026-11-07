"""Live data from a Garmin handheld GPS over USB (Garmin's own USB protocol, through Garmin's USB driver).

With the handheld's Interface set to Garmin Serial or Garmin Spanner, and "No" to mass storage, Windows shows it as
"Garmin USB GPS". This talks to it the way Garmin's protocol document (Garmin Device Interface Specification) and
GPSBabel do, and asks for its position once a second:
  PVT (packet 51, D800)   position, height (above the WGS84 ellipsoid and above sea level), Garmin's own error
                          estimates (EPE / EPH / EPV), fix type (2D / 3D, 'diff' = with WAAS corrections), GPS time
  satellites (packet 114) each channel's satellite, signal level, elevation, azimuth and status (has ephemeris,
                          WAAS-corrected, used in the fix) - undocumented, as GPSBabel reads it

  python garmin_usb.py [seconds]        prints what it receives
  python garmin_usb.py log FILE [--parent-pid PID]
                                        records it to FILE until stopped, reconnecting as it comes and goes (and
                                        ending when the program PID ends - receiver.py runs it this way beside a GNSS
                                        capture, writing <capture>.garmin)
  GarminUSB().open(); .product(); .start_pvt(); .read() -> ('pvt', dict) / ('sats', list) / (other, None)

The side-by-side comparison (catalog.py after Stop, the tracker's position panel while it records):
  LogReader / Ours collect the two receivers' lines, compare() sums them up, panel_text() says it in a line,
  Beside does it live for the position panel.
"""
import ctypes
import math
import os
import struct
import sys
import threading
import time
from ctypes import wintypes

GUID_GARMIN = (0x2C9C45C2, 0x8E7D, 0x4C08, (0xA1, 0x2D, 0x81, 0x6B, 0xBA, 0xE7, 0x22, 0xC0))
IOCTL_ASYNC_IN = 0x00222140
IOCTL_USB_PACKET_SIZE = 0x00222144
IOCTL_API_VERSION = 0x00222000
INTERRUPT_SIZE = 64
FIX_TYPES = {0: 'none', 1: 'invalid', 2: '2D', 3: '3D', 4: '2D+WAAS', 5: '3D+WAAS'}
GARMIN_EPOCH = 631065600            # 1989-12-31 00:00 UTC (Garmin's day count starts here)

# packet types / ids
T_USB, T_APP = 0, 20
P_DATA_AVAILABLE, P_START_SESSION, P_SESSION_STARTED = 2, 5, 6
P_COMMAND, P_PROTOCOLS, P_PRODUCT_RQST, P_PRODUCT_DATA = 10, 253, 254, 255
P_PVT, P_SATS = 51, 114
CMD_START_PVT, CMD_STOP_PVT = 49, 50


class GUID(ctypes.Structure):
    _fields_ = [('Data1', wintypes.DWORD), ('Data2', wintypes.WORD), ('Data3', wintypes.WORD),
                ('Data4', ctypes.c_ubyte * 8)]


class SP_DEVICE_INTERFACE_DATA(ctypes.Structure):
    _fields_ = [('cbSize', wintypes.DWORD), ('InterfaceClassGuid', GUID), ('Flags', wintypes.DWORD),
                ('Reserved', ctypes.c_void_p)]


def _guid():
    g = GUID()
    g.Data1, g.Data2, g.Data3 = GUID_GARMIN[0], GUID_GARMIN[1], GUID_GARMIN[2]
    for i, b in enumerate(GUID_GARMIN[3]):
        g.Data4[i] = b
    return g


def device_paths():
    """the Garmin USB devices present (their Windows device paths)"""
    sa = ctypes.WinDLL('setupapi', use_last_error=True)
    sa.SetupDiGetClassDevsW.restype = ctypes.c_void_p
    sa.SetupDiGetClassDevsW.argtypes = [ctypes.POINTER(GUID), ctypes.c_wchar_p, ctypes.c_void_p, wintypes.DWORD]
    sa.SetupDiEnumDeviceInterfaces.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.POINTER(GUID), wintypes.DWORD,
                                               ctypes.POINTER(SP_DEVICE_INTERFACE_DATA)]
    sa.SetupDiGetDeviceInterfaceDetailW.argtypes = [ctypes.c_void_p, ctypes.POINTER(SP_DEVICE_INTERFACE_DATA),
                                                    ctypes.c_void_p, wintypes.DWORD, ctypes.POINTER(wintypes.DWORD),
                                                    ctypes.c_void_p]
    sa.SetupDiDestroyDeviceInfoList.argtypes = [ctypes.c_void_p]
    g = _guid()
    h = sa.SetupDiGetClassDevsW(ctypes.byref(g), None, None, 0x12)      # DIGCF_PRESENT | DIGCF_DEVICEINTERFACE
    if not h or h == ctypes.c_void_p(-1).value:
        return []
    out = []
    try:
        i = 0
        while True:
            d = SP_DEVICE_INTERFACE_DATA()
            d.cbSize = ctypes.sizeof(d)
            if not sa.SetupDiEnumDeviceInterfaces(h, None, ctypes.byref(g), i, ctypes.byref(d)):
                break
            need = wintypes.DWORD()
            sa.SetupDiGetDeviceInterfaceDetailW(h, ctypes.byref(d), None, 0, ctypes.byref(need), None)
            buf = ctypes.create_string_buffer(need.value)
            # SP_DEVICE_INTERFACE_DETAIL_DATA_W: cbSize (8 on 64-bit Windows) then the path
            struct.pack_into('I', buf, 0, 8 if ctypes.sizeof(ctypes.c_void_p) == 8 else 6)
            if sa.SetupDiGetDeviceInterfaceDetailW(h, ctypes.byref(d), buf, need, None, None):
                out.append(ctypes.wstring_at(ctypes.addressof(buf) + 4))
            i += 1
    finally:
        sa.SetupDiDestroyDeviceInfoList(h)
    return out


class GarminUSB:
    def __init__(self, path=None):
        self.path = path
        self.h = None
        self.k32 = ctypes.WinDLL('kernel32', use_last_error=True)
        self.k32.CreateFileW.restype = wintypes.HANDLE
        self.k32.CreateFileW.argtypes = [ctypes.c_wchar_p, wintypes.DWORD, wintypes.DWORD, ctypes.c_void_p,
                                         wintypes.DWORD, wintypes.DWORD, wintypes.HANDLE]
        self.k32.DeviceIoControl.argtypes = [wintypes.HANDLE, wintypes.DWORD, ctypes.c_void_p, wintypes.DWORD,
                                             ctypes.c_void_p, wintypes.DWORD, ctypes.POINTER(wintypes.DWORD),
                                             ctypes.c_void_p]
        self.bulk = False
        self.packet_size = 64
        self.unit_id = None

    def open(self):
        if not self.path:
            paths = device_paths()
            if not paths:
                raise OSError('no Garmin USB GPS found (driver installed? interface Garmin Serial/Spanner, and "No" '
                              'to mass storage?)')
            self.path = paths[0]
        self.h = self.k32.CreateFileW(self.path, 0xC0000000, 0, None, 3, 0x80, None)
        if self.h == wintypes.HANDLE(-1).value:
            raise OSError('could not open the Garmin (%d) - is another program using it?' % ctypes.get_last_error())
        n = wintypes.DWORD()
        size = ctypes.c_uint32()
        if self.k32.DeviceIoControl(self.h, IOCTL_USB_PACKET_SIZE, None, 0, ctypes.byref(size), 4, ctypes.byref(n), None):
            self.packet_size = size.value or 64
        # start a session: the unit answers with its ID
        for _ in range(3):
            self.send(T_USB, P_START_SESSION)
            for _ in range(10):
                t, pid, data = self.recv()
                if t == T_USB and pid == P_SESSION_STARTED:
                    self.unit_id = struct.unpack('<I', data[:4])[0] if len(data) >= 4 else None
                    return self
        raise OSError('the Garmin did not answer the session start')

    def close(self):
        if self.h:
            try:
                self.send(T_APP, P_COMMAND, struct.pack('<H', CMD_STOP_PVT))
            except OSError:
                pass
            self.k32.CloseHandle(self.h)
            self.h = None

    def send(self, ptype, pid, data=b''):
        pkt = struct.pack('<B3xH2xI', ptype, pid, len(data)) + data
        n = wintypes.DWORD()
        if not self.k32.WriteFile(self.h, pkt, len(pkt), ctypes.byref(n), None):
            raise OSError('write to the Garmin failed (%d)' % ctypes.get_last_error())
        if len(pkt) % self.packet_size == 0:
            self.k32.WriteFile(self.h, b'', 0, ctypes.byref(n), None)

    def recv(self):
        """one packet -> (type, id, data); (None, None, b'') when nothing came"""
        buf = ctypes.create_string_buffer(4096)
        n = wintypes.DWORD()
        if self.bulk:
            ok = self.k32.ReadFile(self.h, buf, 4096, ctypes.byref(n), None)
            if not ok or n.value == 0:
                self.bulk = False
                return None, None, b''
        else:
            ok = self.k32.DeviceIoControl(self.h, IOCTL_ASYNC_IN, None, 0, buf, INTERRUPT_SIZE, ctypes.byref(n), None)
            if not ok or n.value == 0:
                return None, None, b''
        raw = buf.raw[:n.value]
        if len(raw) < 12:
            return None, None, b''
        ptype, pid, size = struct.unpack('<B3xH2xI', raw[:12])
        if ptype == T_USB and pid == P_DATA_AVAILABLE:
            self.bulk = True        # the rest comes on the bulk pipe
        return ptype, pid, raw[12:12 + size]

    def product(self):
        """(product id, software version, description)"""
        self.send(T_APP, P_PRODUCT_RQST)
        for _ in range(30):
            t, pid, d = self.recv()
            if pid == P_PRODUCT_DATA and len(d) >= 4:
                prod, sw = struct.unpack('<Hh', d[:4])
                return prod, sw / 100.0, d[4:].split(b'\0')[0].decode(errors='replace')
        return None

    def start_pvt(self):
        self.send(T_APP, P_COMMAND, struct.pack('<H', CMD_START_PVT))

    def read(self):
        t, pid, d = self.recv()
        if t == T_APP and pid == P_PVT and len(d) >= 64:
            return 'pvt', parse_pvt(d)
        if t == T_APP and pid == P_SATS:
            return 'sats', parse_sats(d)
        return pid, None


def parse_pvt(d):
    """D800 position / velocity / time"""
    alt, epe, eph, epv, fix, tow, lat, lon, ve, vn, vu, msl, leap, wn_days = struct.unpack('<ffffhdddffffhi', d[:64])
    t = GARMIN_EPOCH + wn_days * 86400 + tow - leap
    return {'t': t, 'lat': math.degrees(lat), 'lon': math.degrees(lon), 'alt_ell': alt, 'msl': alt + msl,
            'geoid_sep': -msl, 'epe': epe, 'eph': eph, 'epv': epv, 'fix': FIX_TYPES.get(fix, str(fix)),
            'vel': (ve, vn, vu)}


def parse_sats(d):
    """the satellite records (7 bytes each, as GPSBabel reads them): prn, snr, elevation, azimuth, status bits
    (1 has ephemeris, 2 WAAS-corrected, 4 used in the fix, 8 rising)"""
    out = []
    for i in range(0, len(d) - 6, 7):
        svid, snr, el, az, st = struct.unpack('<BHBHB', d[i:i + 7])
        if svid == 255:
            continue
        out.append({'prn': svid, 'snr': snr, 'el': el, 'az': az, 'status': st,
                    'eph': bool(st & 1), 'diff': bool(st & 2), 'used': bool(st & 4)})
    return out


# ---- side by side: the Garmin's record (<capture>.garmin) against our receiver's lines of the same capture
# A Garmin second counts as the truth only with a WAAS-corrected 3D fix and its own horizontal error estimate (EPH)
# within 5 m: without WAAS it was 9 m (median) from home indoors with an EPH of 43 m during the 2026-11-06 review -
# no better than what it would be judging. Its heights are barometric (they drift 12-50 m), so only horizontal
# distances are compared, and its "SNR" is its own scale with no relation to our C/N0 (fitted slope 0.05, review) -
# levels are never compared, only who holds which satellite.
TRUTH_FIX = '3D+WAAS'
TRUTH_EPH_M = 5.0
MIN_USED_S = 10         # 'used by the Garmin, never by us' needs this long in its fix (a blip is not a satellite)
MOSTLY_SHARE = 0.1      # ... and 'mostly not by us': we held it under a tenth of the time it used it
# what each of our bands could have received that the Garmin also does (it has GPS, SBAS and GLONASS; no Galileo or
# BeiDou) - a GLONASS satellite it uses is no gap of the GPS band's
BAND_LETTERS = {'gnss1': 'GS', 'glo1': 'R', 'bds1': ''}
SETUP_HINT = ('On the Garmin: Setup > System > GPS + GLONASS and WAAS/EGNOS On; Interface Garmin Spanner and "No" '
              'to its mass-storage question; put it 1-2 m from the GPS puck, never on it (it blocks the puck).')


def sat_name(prn):
    """the Garmin's satellite number -> our name: 1-32 GPS (G05), 33-51 SBAS = PRN 120-138 (S31 = PRN 131, WAAS),
    65-96 GLONASS slot + 64 (R03)"""
    if 1 <= prn <= 32:
        return 'G%02d' % prn
    if 33 <= prn <= 51:
        return 'S%02d' % (prn + 87 - 100)
    if 65 <= prn <= 96:
        return 'R%02d' % (prn - 64)
    return '?%d' % prn


def dist_m(lat1, lon1, lat2, lon2):
    """horizontal distance in metres (flat over the few metres compared here)"""
    return math.hypot(math.radians(lon2 - lon1) * 6371000 * math.cos(math.radians(lat1)),
                      math.radians(lat2 - lat1) * 6371000)


def compass(az):
    return ('N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW')[int((az % 360) / 45 + 0.5) % 8]


class LogReader:
    """the Garmin's lines (GSTATE / GPVT / GSATS, see log_forever), fed whole files or a growing file's new text"""

    def __init__(self):
        self.pvt = []           # {'pc', 't', 'lat', 'lon', 'epe', 'eph', 'fix'}, in order
        self.sats = {}          # PC second -> [(our name, snr, el, az, status)]
        self.unit = None
        self.connected = None   # the newest GSTATE: True connected, False lost, None none yet
        self.connects = 0
        self.bad = 0
        self.rest = ''

    def feed_text(self, text, final=False):
        """new text of the file; a last line without its newline waits for the rest (or, final, is dropped - the
        recorder is killed on Stop and may have been half way through a line)"""
        lines = (self.rest + text).split('\n')
        self.rest = '' if final else lines.pop()
        if final and lines and not text.endswith('\n'):
            lines.pop()
        self.feed(lines)

    def feed(self, lines):
        for line in lines:
            f = line.split()
            if not f:
                continue
            try:
                if f[0] == 'GPVT:' and len(f) >= 11:
                    self.pvt.append({'pc': float(f[1]), 't': float(f[2]), 'lat': float(f[3]), 'lon': float(f[4]),
                                     'epe': float(f[7]), 'eph': float(f[8]), 'fix': f[10]})
                elif f[0] == 'GSATS:' and len(f) >= 2:
                    recs = []
                    for e in f[2:]:
                        prn, snr, el, az, st = map(int, e.split(':'))
                        recs.append((sat_name(prn), snr, el, az, st))
                    self.sats[int(round(float(f[1])))] = recs
                elif f[0] == 'GSTATE:' and len(f) >= 3:
                    float(f[1])
                    self.connected = f[2] == 'connected'
                    if self.connected:
                        self.connects += 1
                        self.unit = ' '.join(f[3:]) or self.unit
            except (ValueError, IndexError):
                self.bad += 1


class Ours:
    """our receiver's lines of the same capture (catalog.read_gnss rows): per satellite the seconds it was locked
    (GCH lines are stamped with this PC's clock, as the Garmin's satellite lists are), and the position fixes
    (stamped with GPS time in UTC, as the Garmin's PVT is)"""

    def __init__(self):
        self.gch = {}           # sat -> {PC second: C/N0}
        self.gpos = {}          # UTC second -> (lat, lon, ok)       pocket_trk's fixes
        self.gspp = {}          # UTC second -> (lat, lon, ok)       our own solver's; ok = not a quick fix

    def feed_rows(self, rows):
        for kind, t, f in rows:
            try:
                if kind == 'GCH' and len(f) >= 3:
                    self.gch.setdefault(f[0], {})[int(round(t))] = float(f[2])
                elif kind == 'GPOS' and len(f) >= 4:
                    self.gpos[int(round(t))] = (float(f[0]), float(f[1]), True)
                elif kind == 'GSPP' and len(f) >= 5:
                    # (a quick fix before carrier smoothing settled was 11-31 m off beside the Garmin: not compared)
                    self.gspp[int(round(t))] = (float(f[0]), float(f[1]), not (len(f) >= 9 and int(f[8]) < 4))
            except (ValueError, IndexError):
                pass


def compare(g, o, band=None, home=None):
    """a LogReader and an Ours of the same capture -> the summary the catalog keeps (cap['garmin']) and the position
    panel shows: the Garmin's fixes (types, EPH), per-second horizontal distances to our two solutions over the truth
    window (or, without one, its plain 3D fixes - 'basis' says which), both receivers' middles over it, its GPS-minus-PC
    clock, and per satellite who held it ({garmin_s, garmin_used_s, our_locked_s, both_s, el, az, our_cn0}) with
    the ones the Garmin used that we never locked. home (lat, lon[, h]) adds the Garmin middle's distance from it."""
    import collections
    import numpy as np
    med = lambda v: round(float(np.median(v)), 2)
    s = {'unit': g.unit, 'connects': g.connects, 'pvt_s': len(g.pvt)}
    pcs = [p['pc'] for p in g.pvt] + list(g.sats)
    if not pcs:
        s['note'] = 'nothing from the Garmin (is its Interface Garmin Spanner, and was "No" the answer to mass storage?)'
        return s
    t0, t1 = min(pcs), max(pcs)
    s['first'], s['last'] = round(t0, 1), round(t1, 1)
    s['fix'] = dict(collections.Counter(p['fix'] for p in g.pvt))
    fixes = [p for p in g.pvt if p['fix'][:2] in ('2D', '3D')]
    if fixes:
        s['eph_m'] = med([p['eph'] for p in fixes])
        s['epe_m'] = med([p['epe'] for p in fixes])
        # its GPS time against this PC's clock (a cross-check of gnss_rx's GTIME: the packet arrives a little after the
        # fix's moment, so the true offset is slightly larger than this)
        s['clock_s'] = med([p['t'] - p['pc'] for p in fixes])
    s['waas_seen'] = (any(p['fix'].endswith('WAAS') for p in g.pvt) or
                      any(r[4] & 2 for recs in g.sats.values() for r in recs))
    truth = [p for p in fixes if p['fix'] == TRUTH_FIX and p['eph'] <= TRUTH_EPH_M]
    basis = truth or [p for p in fixes if p['fix'].startswith('3D')]
    s['basis'] = 'waas' if truth else '3d' if basis else None
    s['basis_s'] = len(basis)
    ours = (('own', o.gspp), ('pocket', o.gpos))
    if basis:
        s['basis_eph_m'] = med([p['eph'] for p in basis])
        for key, pts in ours:          # the same second: the Garmin's fix time and ours are both GPS time
            d = []
            for p in basis:
                v = pts.get(int(round(p['t'])))
                if v and v[2]:
                    d.append(dist_m(p['lat'], p['lon'], v[0], v[1]))
            if d:
                s['vs_' + key] = {'n': len(d), 'median_m': med(d), 'p95_m': round(float(np.percentile(d, 95)), 2)}
        glat, glon = float(np.median([p['lat'] for p in basis])), float(np.median([p['lon'] for p in basis]))
        mid = {'garmin_n': len(basis)}
        if home:
            mid['garmin_from_home_m'] = round(dist_m(home[0], home[1], glat, glon), 2)
        b0, b1 = min(p['t'] for p in basis) - 0.5, max(p['t'] for p in basis) + 0.5
        for key, pts in ours:          # each of ours over the same stretch of time
            v = [x for sec, x in pts.items() if b0 <= sec <= b1 and x[2]]
            if v:
                mid[key + '_n'] = len(v)
                mid[key + '_m'] = round(dist_m(glat, glon, float(np.median([x[0] for x in v])),
                                               float(np.median([x[1] for x in v]))), 2)
        s['middle'] = mid
    # per satellite, over the time the Garmin was recorded (seconds of this PC's clock)
    acc = {}
    for sec, recs in g.sats.items():
        for name, snr, el, az, st in recs:
            r = acc.setdefault(name, [0, 0, [], [], set()])
            if snr > 0:
                r[0] += 1
                r[4].add(sec)
            r[1] += bool(st & 4)
            if el > 0 or az > 0:        # (0/0: it doesn't know where the satellite is yet)
                r[2].append(el)
                r[3].append(az)
    lo, hi = t0 - 0.5, t1 + 0.5
    letters = BAND_LETTERS[band] if band in BAND_LETTERS else ''.join(sorted({n[0] for n in o.gch}))
    heard_g = {name for name, r in acc.items() if r[0]}
    table = {}
    for name in sorted(set(acc) | {n for n in o.gch if n[0] in 'GSR'}):
        r = acc.get(name, [0, 0, [], [], set()])
        mine = {sec: c for sec, c in o.gch.get(name, {}).items() if lo <= sec <= hi}
        if not (r[0] or r[1] or mine):
            continue
        e = {'garmin_s': r[0], 'garmin_used_s': r[1], 'our_locked_s': len(mine), 'both_s': len(r[4] & set(mine))}
        if r[2]:
            e['el'] = int(round(float(np.median(r[2]))))
            # (the azimuth's mean direction, not its median: a satellite in the north crosses 0/360)
            a = np.radians(r[3])
            e['az'] = int(round(math.degrees(math.atan2(np.sin(a).mean(), np.cos(a).mean())))) % 360
        if mine:
            e['our_cn0'] = round(float(np.median(list(mine.values()))), 1)
        table[name] = e
    s['sats'] = table
    used = [(n, e) for n, e in table.items() if n[:1] in letters and n[:1] != '?' and e['garmin_used_s'] >= MIN_USED_S]
    s['garmin_only'] = [n for n, e in sorted(used, key=lambda x: -x[1]['garmin_used_s']) if not e['our_locked_s']]
    s['garmin_mostly'] = [n for n, e in sorted(used, key=lambda x: -x[1]['garmin_used_s'])
                          if 0 < e['our_locked_s'] < MOSTLY_SHARE * e['garmin_used_s']]
    s['heard'] = dict(collections.Counter(n[0] for n in heard_g))
    return s


def _where(names, table):
    """'W, 20-60°' for a list of satellites: their compass directions and elevation range"""
    dirs, els = [], []
    for n in names:
        e = table.get(n, {})
        if 'az' in e:
            c = compass(e['az'])
            if c not in dirs:
                dirs.append(c)
            els.append(e['el'])
    if not els:
        return ''
    el = '%d°' % els[0] if min(els) == max(els) else '%d-%d°' % (min(els), max(els))
    return ' (%s, %s)' % ('/'.join(dirs), el)


def panel_text(s, live=True, gtime=None):
    """a compare() summary -> (one line for the position panel, its tooltip). live: it is being recorded now (else
    it is the newest catalogued capture with one, dated); gtime: our receiver's GPS-minus-PC clock offset right now
    (gnss_rx's GTIME) for the clock cross-check"""
    if not s:
        return '', ''
    when = time.strftime('%d %b %H:%M', time.localtime(s['first'])) if s.get('first') else ''
    head = 'Garmin beside you' + ('' if live else ' (last time, %s)' % when)
    tip = ['Garmin beside the antenna%s - %s' % ('' if live else ' (the last capture with one)', s.get('unit') or
                                                  'Garmin handheld')]
    if not s.get('pvt_s'):
        line = head + (': no position from it yet' if live else ': it sent no position')
        if s.get('note'):
            line += ' - ' + s['note']
        return line, '\n'.join(tip + ['', SETUP_HINT])
    fix = '3D+WAAS' if s.get('basis') == 'waas' else ('3D without WAAS (not a truth reference: turn WAAS on)'
                                                       if s.get('basis') else 'no 3D fix')
    parts = [fix]
    if s.get('basis_eph_m') is not None:
        parts.append('EPH %.0f m' % s['basis_eph_m'])
    mid = s.get('middle', {})
    if 'own_m' in mid:
        parts.append('%.1f m from our average' % mid['own_m'])
    elif 'pocket_m' in mid:
        parts.append('%.1f m from pocket_trk\'s average' % mid['pocket_m'])
    elif s.get('basis'):
        parts.append('no fix of ours to compare yet' if live else 'no fix of ours to compare')
    line = head + ': ' + ', '.join(parts)
    table = s.get('sats', {})
    if s.get('garmin_only'):
        line += '; it holds %s%s that we do not' % (' '.join(s['garmin_only']), _where(s['garmin_only'], table))
    # the tooltip: the numbers behind the line
    span = (s['last'] - s['first']) if s.get('first') else 0
    tip.append('   %d s of positions over %.0f min (%s); its fixes: %s; EPH median %s m, EPE %s m' % (
        s['pvt_s'], span / 60, when, ', '.join('%s %d s' % kv for kv in sorted(s.get('fix', {}).items())),
        '%.0f' % s['eph_m'] if 'eph_m' in s else '-', '%.0f' % s['epe_m'] if 'epe_m' in s else '-'))
    if s.get('basis') == 'waas':
        tip.append('   truth window: %d s with 3D+WAAS and EPH <= %.0f m' % (s['basis_s'], TRUTH_EPH_M))
    else:
        tip.append('   no truth window: it had no WAAS-corrected fix with EPH <= %.0f m%s - compared with its plain 3D '
                   'fixes (%d s), which are no better than ours' % (TRUTH_EPH_M, '' if s.get('waas_seen') else
                                                                    ' (WAAS was off or not heard)', s['basis_s']))
    for key, label in (('own', 'our solver'), ('pocket', 'pocket_trk')):
        v = s.get('vs_' + key)
        tip.append('   same second, horizontally, %s: %s' % (label, 'median %.1f m, 95 %% within %.1f m (%d s)' % (
            v['median_m'], v['p95_m'], v['n']) if v else 'no fixes at the same seconds'))
    if mid:
        m = ['our solver %.1f m' % mid['own_m'] if 'own_m' in mid else None,
             'pocket_trk %.1f m' % mid['pocket_m'] if 'pocket_m' in mid else None,
             'the Garmin\'s from home %.1f m' % mid['garmin_from_home_m'] if 'garmin_from_home_m' in mid else None]
        tip.append('   middles over the same time, from the Garmin\'s: ' + '; '.join(x for x in m if x))
    if s.get('clock_s') is not None:
        tip.append('   PC clock: the Garmin\'s GPS time is %+.2f s from this PC\'s%s (its packets arrive a moment '
                   'after the fix, so the true figure is slightly larger)' % (
                       s['clock_s'], '; our receiver says %+.2f s' % gtime if gtime is not None else ''))
    if table:
        tip.append('')
        tip.append('Satellites - seconds the Garmin heard / used in its fix / we locked (and both heard):')
        order = sorted(table.items(), key=lambda kv: (-kv[1]['garmin_used_s'], -kv[1]['our_locked_s'], kv[0]))
        for n, e in order[:16]:
            note = (' - used by it, never locked by us' if n in s.get('garmin_only', []) else
                    ' - used by it, hardly locked by us' if n in s.get('garmin_mostly', []) else '')
            tip.append('   %-4s %5d / %5d / %5d (%d)  %s%s%s' % (
                n, e['garmin_s'], e['garmin_used_s'], e['our_locked_s'], e['both_s'],
                '%2d° up, az %3d° %-2s' % (e['el'], e['az'], compass(e['az'])) if 'az' in e else 'position ?',
                '  ours %.0f dB-Hz' % e['our_cn0'] if 'our_cn0' in e else '', note))
        if len(order) > 16:
            tip.append('   ... and %d more' % (len(order) - 16))
        if s.get('heard') is not None and not s['heard'].get('R'):
            tip.append('   it heard no GLONASS satellite (Setup > System > GPS + GLONASS)')
    tip += ['', 'The Garmin\'s "SNR" is its own scale with no relation to our C/N0, so levels are not compared; its '
            'heights are barometric, so only horizontal distances are.', SETUP_HINT]
    return line, '\n'.join(tip)


class Beside:
    """the position panel's Garmin line, asked every 2 s by the tracker: live from the recorder's file and the
    capture's own live file while a GNSS capture records one; else the newest catalogued capture that had one.
    Reads only new text each time (at most CHUNK a go, so a long file catching up never freezes the window) and sums
    up at most every SUM_S."""
    CHUNK = 1 << 20
    SUM_S = 10.0

    def __init__(self):
        self.stem = None
        self.cached = ('', '')

    def text(self, rx, caps=None, gtime=None):
        """rx: the tracker's Receiver; caps: the catalog's captures (name -> capture, or name -> (key, capture)).
        -> (line, tooltip), ('', '') when there is nothing to say. Never raises."""
        try:
            return self._text(rx, caps, gtime)
        except Exception as e:
            print('Garmin beside: %r' % e)
            return '', ''

    def _text(self, rx, caps, gtime):
        path = rx.garmin_file() if hasattr(rx, 'garmin_file') else None
        if path:
            if path != self.stem:
                self.stem, self.g, self.o = path, LogReader(), Ours()
                self.cached = ('', '')
                self.pos = {'g': 0, 'o': 0}
                self.orest, self.done_t, self.behind = '', 0.0, False
                self.band = getattr(rx, 'system', {}).get('key')
            got = (self._tail('g', path, self.g.feed_text),
                   self._tail('o', path[:-len('.garmin')] + '.live', self._feed_ours))
            if 2 in got:            # still catching up with a long file: sum up once it has all been read
                self.behind = True
                return self.cached if self.cached[0] else ('Garmin beside you: reading its record...', '')
            due = any(got) and time.time() - self.done_t >= self.SUM_S
            if due or self.behind or not self.cached[0]:
                self.behind = False
                self.done_t = time.time()
                import catalog
                self.cached = panel_text(compare(self.g, self.o, self.band, catalog.HOME), True, gtime)
            return self.cached
        self.stem = None
        last = None
        for name in sorted(caps or {}, reverse=True):
            c = caps[name]
            c = c[1] if isinstance(c, tuple) else c
            if isinstance(c, dict) and c.get('garmin') and not c['garmin'].get('error'):
                last = c['garmin']
                break
        return panel_text(last, False, gtime) if last else ('', '')

    def _tail(self, key, path, feed):
        """the file's new text into feed -> 0 nothing new, 1 new text, 2 new text and more still to read"""
        try:
            with open(path, 'rb') as f:
                f.seek(self.pos[key])
                data = f.read(self.CHUNK)
        except OSError:
            return 0
        if not data:
            return 0
        self.pos[key] += len(data)
        feed(data.decode('utf-8', 'replace').replace('\r\n', '\n'))
        return 2 if len(data) == self.CHUNK else 1

    def _feed_ours(self, text):
        import catalog
        lines = (self.orest + text).split('\n')
        self.orest = lines.pop()
        self.o.feed_rows(catalog.read_gnss(lines))


def parent_alive(pid):
    """is the program that started us still running? (as sdr_taps.parent_alive, without loading GNU Radio)"""
    if not pid:
        return True
    k = ctypes.windll.kernel32
    h = k.OpenProcess(0x1000, False, pid)          # PROCESS_QUERY_LIMITED_INFORMATION
    if not h:
        return False
    try:
        code = ctypes.c_ulong()
        return bool(k.GetExitCodeProcess(h, ctypes.byref(code))) and code.value == 259     # STILL_ACTIVE
    finally:
        k.CloseHandle(h)


def _watch_parent(pid):
    """end this process when the tracker that started it is gone: a recorder left behind would hold the Garmin (only
    one program can open it) and keep writing into a finished capture. Its file needs no tidying - every line is
    appended and closed on its own."""
    def run():
        while parent_alive(pid):
            time.sleep(2.0)
        os._exit(0)
    threading.Thread(target=run, daemon=True).start()


def log_forever(path, parent_pid=0):
    """record everything the Garmin sends to a file, reconnecting whenever it is unplugged or powered off:
       GSTATE: pc_time text              connected / lost, the product
       GPVT: pc_time gps_time lat lon h_ellipsoid h_msl epe eph epv fix
       GSATS: pc_time prn:snr:el:az:status ...
    pc_time is this PC's clock; gps_time the Garmin's own (both Unix seconds)"""
    def w(line):
        with open(path, 'a') as f:
            f.write(line + '\n')

    if parent_pid:
        _watch_parent(parent_pid)
    while True:
        try:
            g = GarminUSB().open()
        except OSError:
            time.sleep(1.0)
            continue
        try:
            prod = g.product()
            w('GSTATE: %.3f connected %s' % (time.time(), prod[2] if prod else '?'))
            g.start_pvt()
            geoid_done = False
            last = time.time()
            while time.time() - last < 5:       # nothing for 5 s: unplugged or switched off
                k, v = g.read()
                now = time.time()
                if k is None:
                    time.sleep(0.05)
                    continue
                last = now
                if k == 'pvt':
                    if v['fix'].startswith('3D') and not geoid_done:
                        geoid_done = True       # the local sea-level offset, for the tracker's heights
                        try:
                            import refpoints
                            refpoints.set_geoid_sep(v['alt_ell'] - v['msl'], prod[2] if prod else 'Garmin')
                        except Exception:
                            pass
                    w('GPVT: %.3f %.3f %.8f %.8f %.2f %.2f %.2f %.2f %.2f %s' % (
                        now, v['t'], v['lat'], v['lon'], v['alt_ell'], v['msl'], v['epe'], v['eph'], v['epv'], v['fix']))
                elif k == 'sats':
                    w('GSATS: %.3f %s' % (now, ' '.join('%d:%d:%d:%d:%d' % (s['prn'], s['snr'], s['el'], s['az'],
                                                                             s['status']) for s in v)))
        except OSError:
            pass
        w('GSTATE: %.3f lost' % time.time())
        try:
            g.close()
        except OSError:
            pass
        time.sleep(1.0)


def main():
    if len(sys.argv) > 2 and sys.argv[1] == 'log':
        pid = int(sys.argv[sys.argv.index('--parent-pid') + 1]) if '--parent-pid' in sys.argv[3:-1] else 0
        log_forever(sys.argv[2], pid)
        return
    secs = float(sys.argv[1]) if len(sys.argv) > 1 else 20
    print('Garmin USB devices:', device_paths())
    g = GarminUSB().open()
    print('session started, unit id', g.unit_id, '- packet size', g.packet_size)
    print('product:', g.product())
    g.start_pvt()
    end = time.time() + secs
    kinds = {}
    while time.time() < end:
        k, v = g.read()
        if k is None:
            time.sleep(0.05)
            continue
        kinds[k] = kinds.get(k, 0) + 1
        if k == 'pvt':
            print('PVT  %s  fix %-7s  %.6f %.6f  h(ell) %.1f  msl %.1f  EPE %.1f EPH %.1f EPV %.1f m' % (
                time.strftime('%H:%M:%S', time.gmtime(v['t'])), v['fix'], v['lat'], v['lon'], v['alt_ell'], v['msl'],
                v['epe'], v['eph'], v['epv']))
        elif k == 'sats':
            print('SATS ' + '  '.join('%d:%d/%d%s' % (s['prn'], s['snr'], s['el'], '*' if s['used'] else '')
                                      for s in v))
    print('packets:', kinds)
    g.close()


if __name__ == '__main__':
    main()

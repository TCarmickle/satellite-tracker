"""Sky plot: where in your sky the satellites are now, and where they were when your antenna heard them.

Centre = straight overhead, edge = the horizon, north at the top, east to the right (like a map).
  satellites  the ones above your horizon now (from our own orbits), named; a green ring = heard in the last
              minute, blue brackets = the one selected, whose pass across your sky is drawn as a path
  big dots    positions the satellites reported (Iridium ring alerts, Orbcomm ephemeris), coloured by signal
  shading     everything else heard from identified satellites (frames, packets, GNSS tracking), placed with
              our own orbit for that satellite and counted in 3-degree patches of sky (brighter = more heard)
  signal      (GNSS, the 'Signal' button) the middle C/N0 of everything tracked in each 5-degree patch of sky, each
              satellite placed by its own broadcast orbit: which directions your antenna hears weakly
  dashed ring a satellite our own search heard although the receiver is not tracking it (GHEARD)
Empty parts of the circle are directions or heights your antenna isn't hearing - blocked, or in its
pattern's blind spot (a vertical whip, for instance, hears little straight up).
"""
import math
import time

import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QPointF, Qt

TEXT, MUTED, GRID = '#e6e9ef', '#8b93a1', '#2f3541'
SNR_LO, SNR_HI = 12.0, 30.0
SECTORS = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW']
# 'heard' receptions are counted in cells of the sky this many degrees on a side (elevation and azimuth): a GNSS pass
# has a point per second for hours (760 000 points in the catalog after a few days); drawing a dot for each froze the
# window for seconds, and the cap of 200 dots per pass that replaced it (73fabfc) thinned long passes out
CELL_DEG = 3
N_EL, N_AZ = 90 // CELL_DEG, 360 // CELL_DEG
# the signal layer's cells: the catalog's per-session sky map (catalog.SKY_CELL_DEG), 5 x 5 degrees; C/N0 is counted in
# CN0_STEP bins so the middle value over any period (sessions + the run going on) is one cumulative sum per cell
SIG_DEG = 5
SIG_EL, SIG_AZ = 90 // SIG_DEG, 360 // SIG_DEG
CN0_STEP, CN0_BINS = 0.5, 120              # 0-60 dB-Hz
CN0_LO, CN0_HI = 25.0, 45.0                # the colour scale: GPS here runs ~28 (low, shaded) to ~45 dB-Hz
LOCK_LABELS = ['10-20', '20-30', '30-45', '45-60', '60-90']
SIG_MIN_S = 5                              # a signal cell is drawn from this many seconds on


def snr_colour(snr):
    """weak = blue, medium = green/yellow, strong = red"""
    x = min(1.0, max(0.0, (snr - SNR_LO) / (SNR_HI - SNR_LO)))
    stops = [(0.0, (60, 120, 255)), (0.4, (40, 210, 140)), (0.7, (250, 220, 60)), (1.0, (255, 70, 60))]
    for (x0, c0), (x1, c1) in zip(stops, stops[1:]):
        if x <= x1:
            k = (x - x0) / (x1 - x0)
            return QtGui.QColor(*(int(a + (b - a) * k) for a, b in zip(c0, c1)))
    return QtGui.QColor(*stops[-1][1])


def cn0_colour(cn0):
    """the signal layer's colour: the same scale as the reported positions (weak blue ... strong red), over C/N0"""
    return snr_colour(SNR_LO + (cn0 - CN0_LO) / (CN0_HI - CN0_LO) * (SNR_HI - SNR_LO))


def sky_fixes(data, from_home, since=0.0, shows=None):
    """[(el, az, snr, t)]: the positions satellites reported (Iridium ring alerts, Orbcomm) since `since`;
    shows(name, capture) -> bool limits it to the systems being shown"""
    shows = shows or (lambda name, cap=None: True)
    fx = [f for f in data.all_fixes() if f.get('xyz_km') and f['t'] >= since and shows(f['name'], f['cap'])]
    if not fx:
        return []
    X, Y, Z = np.array([f['xyz_km'] for f in fx], float).T
    el, az, _ = from_home(X, Y, Z)
    return [(float(e), float(a), f['snr'], f['t']) for e, a, f in zip(el, az, fx) if e > -2]


class HeardCells:
    """everything heard from identified satellites (frames, packets, GNSS tracking), each reception placed with our
    own orbit for that satellite and counted in a CELL_DEG sky cell. Worked out incrementally: each Doppler curve's
    receptions are placed once and kept (cell number + time); a curve still growing live only has its new points
    placed; a satellite whose orbit was refitted is placed again; a moved home starts over.
    Measured on the 11-06 catalog + an 11 h live session (week-2 sky_load.py): 'Everything' took 190 ms to gather
    and 360 ms to paint every 5 s even with the 200-dot cap; now placing all 783 000 receptions takes 400 ms once
    (when the sky is first shown), each 5 s refresh 5-9 ms and its paint 3 ms (25 ms when the cells changed)."""

    def __init__(self):
        self.curves = {}        # id(curve) -> [curve, orbit, points done, cells (int16), times, name, all-time counts]
        self.home = None

    def update(self, data, from_home, shows):
        probe = tuple(float(np.asarray(v).reshape(-1)[0]) for v in from_home(np.array([7000.0]), np.array([0.0]),
                                                                            np.array([0.0])))
        if probe != self.home:          # (home moved: every place in the sky is different)
            self.curves, self.home = {}, probe
        seen = set()
        weights = getattr(data, 'curve_weights', lambda g: None)
        for name, groups in data.doppler.items():
            orbit = data.orbits.get(name)
            if orbit is None:
                continue
            visible = shows(name)
            for g in groups:
                k = id(g)
                e = self.curves.get(k)
                if e is not None and e[0] is g and not visible:
                    seen.add(k)             # (hidden for now: kept as it was, placed when shown)
                    continue
                if not visible:
                    continue
                if e is None or e[0] is not g or e[1] is not orbit or len(g) < e[2]:
                    e = self.curves[k] = [g, orbit, 0, np.zeros(0, np.int16), np.zeros(0), name, None, np.zeros(0)]
                seen.add(k)
                if len(g) > e[2]:
                    w = weights(g)
                    self._place(e, g[e[2]:], from_home, None if w is None else w[e[2]:len(g)])
        for k in [k for k in self.curves if k not in seen]:
            del self.curves[k]

    @staticmethod
    def _place(e, new, from_home, w=None):
        """place new points of a curve; w = how many seconds each stands for (catalog GNSS passes keep a point every
        10 s, catalog.slim_pass) - so the counts stay 'seconds heard', as when every second was kept"""
        orbit = e[1]
        e[2] += len(new)
        ts = np.fromiter((t for t, _ in new), float, len(new))
        wt = np.ones(len(ts)) if w is None or len(w) != len(ts) else np.asarray(w, float)
        good = orbit.sigma_km(ts) < 300             # only where our orbit for it is good
        ts, wt = ts[good], wt[good]
        if not len(ts):
            return
        X, Y, Z = orbit.ecef(ts).T
        el, az, _ = from_home(X, Y, Z)
        keep = el > -2                              # (just below the horizon counts as on it)
        el, az, ts, wt = el[keep], az[keep], ts[keep], wt[keep]
        cell = (np.clip(el, 0, 89.999) // CELL_DEG).astype(np.int16) * N_AZ + ((az % 360) // CELL_DEG).astype(np.int16)
        e[3] = np.concatenate([e[3], cell])
        e[4] = np.concatenate([e[4], ts])
        e[7] = np.concatenate([e[7], wt])
        e[6] = None

    def counts(self, since, shows):
        """(N_EL, N_AZ) receptions per cell since `since`, of the satellites shown"""
        acc = np.zeros(N_EL * N_AZ)
        for e in self.curves.values():
            ts = e[4]
            if not len(ts) or not shows(e[5]):
                continue
            if since <= ts.min():
                if e[6] is None:
                    e[6] = np.bincount(e[3], weights=e[7], minlength=N_EL * N_AZ)
                acc += e[6]
            elif since <= ts.max():
                m = ts >= since
                acc += np.bincount(e[3][m], weights=e[7][m], minlength=N_EL * N_AZ)
        return np.rint(acc).astype(np.int64).reshape(N_EL, N_AZ)


class SignalCells:
    """the signal layer: C/N0 by 5 x 5 degree sky cell, each GNSS satellite placed by its own broadcast orbit.
    Finished sessions bring their sky map from the catalog (catalog.sky_stats: seconds and middle C/N0 per cell);
    the run going on (and any capture not catalogued yet) is placed here from the tracker's live C/N0, incrementally
    like HeardCells. The middle value over a period is a weighted median: each session cell counts as its seconds at
    its own middle value, each live second as itself."""

    def __init__(self):
        self.curves = {}        # (capture, name) -> [list, orbit, points done, cells (int16), times, C/N0 bins, name]
        self.home = None

    def update(self, data, from_home, shows):
        probe = tuple(float(np.asarray(v).reshape(-1)[0]) for v in from_home(np.array([7000.0]), np.array([0.0]),
                                                                            np.array([0.0])))
        if probe != self.home:
            self.curves, self.home = {}, probe
        live = getattr(data, 'live_cn0', {})
        for key in [k for k in self.curves if k not in live]:
            del self.curves[key]
        for key, g in list(live.items()):
            orbit = data.orbits.get(key[1])
            if orbit is None or not shows(key[1]):
                continue
            e = self.curves.get(key)
            if e is None or e[0] is not g or e[1] is not orbit or len(g) < e[2]:
                e = self.curves[key] = [g, orbit, 0, np.zeros(0, np.int16), np.zeros(0), np.zeros(0, np.int16), key[1]]
            if len(g) > e[2]:
                new = np.asarray(g[e[2]:], float).reshape(-1, 2)
                e[2] = len(g)
                ts, c = new[:, 0], new[:, 1]
                good = orbit.sigma_km(ts) < 300
                ts, c = ts[good], c[good]
                if not len(ts):
                    continue
                el, az, _ = from_home(*orbit.ecef(ts).T)
                up = el >= 0
                cell = ((np.minimum(el[up], 89.999) // SIG_DEG).astype(np.int16) * SIG_AZ +
                        ((az[up] % 360) // SIG_DEG).astype(np.int16))
                e[3] = np.concatenate([e[3], cell])
                e[4] = np.concatenate([e[4], ts[up]])
                e[5] = np.concatenate([e[5], np.clip(np.rint(c[up] / CN0_STEP), 0, CN0_BINS - 1).astype(np.int16)])

    def histogram(self, data, since, shows, now=None):
        """(SIG_EL * SIG_AZ, CN0_BINS) seconds per cell and C/N0 bin since `since`: live seconds + the sky maps of
        the catalogued sessions that overlap the period (whole sessions: they are summarised per session)"""
        now = time.time() if now is None else now
        hist = np.zeros(SIG_EL * SIG_AZ * CN0_BINS)
        for e in self.curves.values():
            if not len(e[4]) or not shows(e[6]):
                continue
            m = e[4] >= since
            if m.any():
                hist += np.bincount(e[3][m].astype(np.int64) * CN0_BINS + e[5][m], minlength=len(hist))
        hist = hist.reshape(SIG_EL * SIG_AZ, CN0_BINS)
        for start, end, st, cap in getattr(data, 'gnss_sky', []):
            sm = st.get('skymap')
            if not sm or sm.get('cell_deg') != SIG_DEG or end < since or start > now or not shows('', cap):
                continue
            rows = np.asarray(sm.get('cells', []), float).reshape(-1, 4)
            if not len(rows):
                continue
            cell = rows[:, 0].astype(int) * SIG_AZ + rows[:, 1].astype(int)
            b = np.clip(np.rint(rows[:, 3] * 0.5 / CN0_STEP), 0, CN0_BINS - 1).astype(int)   # (stored / 2 = dB-Hz)
            np.add.at(hist, (cell, b), rows[:, 2])
        return hist


def signal_cells(hist):
    """histogram -> (seconds, median C/N0) per cell, (SIG_EL, SIG_AZ) each (median NaN where nothing)"""
    n = hist.sum(axis=1)
    cum = np.cumsum(hist, axis=1)
    idx = np.argmax(cum >= (n / 2)[:, None], axis=1)
    med = np.where(n > 0, idx * CN0_STEP, np.nan)
    return n.reshape(SIG_EL, SIG_AZ), med.reshape(SIG_EL, SIG_AZ)


def lock_table(data, since, shows, now=None):
    """(up, tracked, sessions) summed over the catalogued GPS sessions overlapping the period: LOCK_STEP_S moments
    a GPS satellite was up per elevation band x 8 sectors, and how many of them it was tracked (catalog.sky_stats)"""
    now = time.time() if now is None else now
    up = trk = None
    n = 0
    for start, end, st, cap in getattr(data, 'gnss_sky', []):
        lk = st.get('lock')
        if not lk or end < since or start > now or not shows('', cap):
            continue
        u, t = np.asarray(lk['up'], float), np.asarray(lk['trk'], float)
        if u.shape != (len(LOCK_LABELS), 8):
            continue
        up, trk = (u, t) if up is None else (up + u, trk + t)
        n += 1
    return up, trk, n


def signal_summary(sig_n, sig_med, lock, heard_nt):
    """the summary lines about the signal layer, the lock table and satellites heard but not tracked"""
    lines = []
    if sig_n is not None and sig_n.sum() > 0:
        el_c, az_c = np.meshgrid((np.arange(SIG_EL) + 0.5) * SIG_DEG, (np.arange(SIG_AZ) + 0.5) * SIG_DEG, indexing='ij')
        ok = sig_n > 0
        n, med, els, azs = sig_n[ok], sig_med[ok], el_c[ok], az_c[ok]
        sector = (((azs + 22.5) % 360) // 45).astype(int)

        def wmed(mask):
            if not mask.any():
                return None
            o = np.argsort(med[mask])
            c = np.cumsum(n[mask][o])
            return float(med[mask][o][np.searchsorted(c, c[-1] / 2)])
        by_sec = [(SECTORS[k], wmed(sector == k)) for k in range(8)]
        got = [(s, v) for s, v in by_sec if v is not None]
        bands = [('0-15', els < 15), ('15-30', (els >= 15) & (els < 30)), ('30-60', (els >= 30) & (els < 60)),
                 ('60-90', els >= 60)]
        lines.append('<b>Signal</b> (middle C/N0 of %s tracked, by direction, each satellite placed by its own orbit): '
                     % _secs(n.sum()) + ', '.join('%s %.0f' % sv for sv in got) + ' dB-Hz; by height ' +
                     ', '.join('%s\u00b0 %.0f' % (b, v) for b, v in ((b, wmed(m)) for b, m in bands) if v is not None)
                     + ' dB-Hz.')
        if len(got) >= 4:
            worst = min(got, key=lambda sv: sv[1])
            best = max(got, key=lambda sv: sv[1])
            if best[1] - worst[1] >= 3:
                lines.append('<span style="color:#ffb454">Weakest towards %s (%.0f dB-Hz, %.0f dB below %s)</span> - '
                             'something there shades the antenna, or its pattern is weak that way.'
                             % (worst[0], worst[1], best[1] - worst[1], best[0]))
    up, trk, n_sess = lock
    if up is not None and up.sum() > 0:
        pct = lambda a, b: '%d%%' % round(100 * a / b) if b else '-'
        lines.append('<b>Locked when up</b> (GPS, own orbits, %d finished session%s): ' % (n_sess, '' if n_sess == 1 else 's')
                     + ', '.join('%s\u00b0 %s' % (lab, pct(trk[i].sum(), up[i].sum())) for i, lab in enumerate(LOCK_LABELS))
                     + '.')
        hi_up, hi_trk = up[1:].sum(axis=0), trk[1:].sum(axis=0)     # above 20 deg, by sector
        secs = [(SECTORS[k], hi_trk[k] / hi_up[k]) for k in range(8) if hi_up[k] >= 20]
        if len(secs) >= 4:
            secs.sort(key=lambda sv: sv[1])
            lines.append('Above 20\u00b0 by direction: ' + ', '.join('%s %d%%' % (s_, round(100 * v)) for s_, v in
                                                                    sorted(secs, key=lambda sv: SECTORS.index(sv[0]))) +
                         '; worst %s.' % ', '.join(s_ for s_, v in secs[:2]))
    if heard_nt:
        lines.append('<b>Heard, not tracked</b> (our own search found them; the receiver has not locked them): ' +
                     ', '.join('%s ~%.0f dB-Hz' % (name, h['cn0']) for name, h in heard_nt) + '.')
    return '<br>'.join(lines)


def _secs(s):
    s = float(s)
    return '%.0f s' % s if s < 120 else '%.0f min' % (s / 60) if s < 7200 else '%.1f h' % (s / 3600)


class SkyPlot(QtWidgets.QWidget):
    satelliteClicked = QtCore.pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.setMinimumSize(300, 240)
        self.fixes = []
        self.cells = np.zeros((N_EL, N_AZ), np.int64)       # receptions per sky cell (HeardCells.counts)
        self.now = []           # satellites up now: dicts with el, az, name, selected, live, el colour, system
        self.path = []          # the selected satellite's pass: [(el, az, t)]
        self.show_heard = self.show_now = self.show_path = True
        self.show_signal = False
        self.signal = None      # (seconds, median C/N0) per SIG cell, or None
        self.hits = []
        # the grid and everything heard are drawn once into a picture and reused: the satellites up now change every
        # second, what was heard only every 5 s (and often not at all) - so a second's repaint draws only the overlay
        self._layer, self._layer_key, self._layer_version = None, None, 0

    @property
    def n_heard(self):
        return int(self.cells.sum())

    def set_points(self, fixes, cells):
        if fixes != self.fixes or not np.array_equal(cells, self.cells):
            self.fixes, self.cells = fixes, cells
            self._layer_version += 1
            self.update()

    def set_signal(self, signal):
        same = (signal is None and self.signal is None) or (
            signal is not None and self.signal is not None and np.array_equal(signal[0], self.signal[0]) and
            np.array_equal(signal[1], self.signal[1], equal_nan=True))
        if not same:
            self.signal = signal
            self._layer_version += 1
            self.update()

    def set_now(self, sats, path):
        self.now, self.path = sats, path
        self.update()

    def mouseReleaseEvent(self, e):
        best = min(((abs(x - e.x()) + abs(y - e.y()), n) for x, y, n in self.hits), default=None)
        self.satelliteClicked.emit(best[1] if best and best[0] < 18 else '')     # ('': empty space - deselect)

    def _xy(self, el, az, cx, cy, R):
        r = (90 - max(0.0, el)) / 90 * R
        a = math.radians(az)
        return QPointF(cx + r * math.sin(a), cy - r * math.cos(a))

    def _geometry(self):
        w, h = self.width(), self.height()
        return w, h, min(w, h) / 2 - 34, w / 2, h / 2

    def _heard_layer(self):
        """the background picture: grid, labels, heard cells, reported positions and their colour scale"""
        dpr = self.devicePixelRatioF()
        key = (self.width(), self.height(), dpr, self.show_heard, self.show_signal, self._layer_version)
        if self._layer is not None and self._layer_key == key:
            return self._layer
        pix = QtGui.QPixmap(max(1, int(self.width() * dpr)), max(1, int(self.height() * dpr)))
        pix.setDevicePixelRatio(dpr)
        p = QtGui.QPainter(pix)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        p.setFont(self.font())
        self._paint_heard(p)
        p.end()
        self._layer, self._layer_key = pix, key
        return pix

    def _paint_heard(self, p):
        w, h, R, cx, cy = self._geometry()
        p.fillRect(QtCore.QRectF(0, 0, w, h), QtGui.QColor('#14171c'))
        font = p.font()
        font.setPointSize(8)
        p.setFont(font)
        # grid: horizon, 30 and 60 degrees up, compass lines
        p.setBrush(QtGui.QColor('#1b1f26'))
        p.setPen(QtGui.QPen(QtGui.QColor(GRID), 1))
        p.drawEllipse(QPointF(cx, cy), R, R)
        # the signal layer (when chosen) instead of the heard shading: each 5-degree cell coloured by its middle C/N0
        # (cells resting on under SIG_MIN_S seconds are left out: a glimpse says little about that direction)
        sig = self.signal if self.show_signal else None
        if sig is not None and sig[0].any():
            p.setPen(Qt.NoPen)
            for i, j in zip(*np.nonzero(sig[0] >= SIG_MIN_S)):
                col = cn0_colour(float(sig[1][i, j]))
                col.setAlpha(220)
                p.setBrush(col)
                p.drawPath(self._cell_path(int(i), int(j), cx, cy, R, SIG_DEG))
            p.setPen(QtGui.QPen(QtGui.QColor(GRID), 1))
        # what was heard, under the grid lines: one wedge per sky cell, brighter where more was heard
        elif self.show_heard and self.cells.any():
            p.setPen(Qt.NoPen)
            top = math.log1p(float(self.cells.max()))
            for i, j in zip(*np.nonzero(self.cells)):
                alpha = int(45 + 150 * math.log1p(float(self.cells[i, j])) / top)
                p.setBrush(QtGui.QColor(150, 170, 200, alpha))
                p.drawPath(self._cell_path(int(i), int(j), cx, cy, R))
            p.setPen(QtGui.QPen(QtGui.QColor(GRID), 1))
        p.setBrush(Qt.NoBrush)
        for el in (30, 60):
            r = (90 - el) / 90 * R
            p.drawEllipse(QPointF(cx, cy), r, r)
        for az in range(0, 360, 45):
            p.drawLine(QPointF(cx, cy), self._xy(0, az, cx, cy, R))
        p.setPen(QtGui.QColor(MUTED))
        for el in (0, 30, 60):
            r = (90 - el) / 90 * R
            p.drawText(QPointF(cx + 4, cy - r + 12), '%d\u00b0' % el)
        p.drawText(QPointF(cx + 4, cy + 12), '90\u00b0 (overhead)')
        for name, az in (('N', 0), ('E', 90), ('S', 180), ('W', 270)):
            pt = self._xy(-8, az, cx, cy, R + 14)
            p.setPen(QtGui.QColor(TEXT))
            f2 = p.font()
            f2.setBold(True)
            p.setFont(f2)
            p.drawText(QtCore.QRectF(pt.x() - 10, pt.y() - 8, 20, 16), Qt.AlignCenter, name)
            p.setFont(font)
        if not self.show_heard:
            return
        # the position fixes on top
        for el, az, snr, _ in sorted(self.fixes, key=lambda f: f[2]):
            p.setPen(QtGui.QPen(QtGui.QColor(0, 0, 0, 140), 1))
            p.setBrush(snr_colour(snr))
            p.drawEllipse(self._xy(el, az, cx, cy, R), 4.5, 4.5)
        if sig is not None:                 # the signal layer's colour scale (bottom right; reported positions' left)
            x0, y0 = w - 112, h - 22
            for i in range(100):
                p.fillRect(QtCore.QRectF(x0 + i, y0, 1.2, 8), cn0_colour(CN0_LO + (CN0_HI - CN0_LO) * i / 99))
            p.setPen(QtGui.QColor(MUTED))
            p.drawText(QtCore.QRectF(w - 300, y0 - 18, 288, 14), Qt.AlignRight | Qt.AlignVCenter,
                       'signal: middle C/N0, %.0f (weak) - %.0f dB-Hz' % (CN0_LO, CN0_HI))
            if not sig[0].any():
                p.drawText(QtCore.QRectF(0, 4, w, 16), Qt.AlignCenter, 'No GNSS signal in this period yet')
        # colour scale (for the reported positions, when there are any)
        if self.fixes:
            x0, y0 = 12, h - 22
            for i in range(100):
                p.fillRect(QtCore.QRectF(x0 + i, y0, 1.2, 8), snr_colour(SNR_LO + (SNR_HI - SNR_LO) * i / 99))
            p.setPen(QtGui.QColor(MUTED))
            p.drawText(QPointF(x0, y0 - 4), 'reported position signal  %.0f - %.0f dB' % (SNR_LO, SNR_HI))

    @staticmethod
    def _cell_path(i, j, cx, cy, R, deg=CELL_DEG):
        """the wedge of sky cell (elevation row i, azimuth column j) of cells `deg` on a side"""
        r_out = (90 - i * deg) / 90 * R
        r_in = (90 - (i + 1) * deg) / 90 * R
        a0 = 90 - j * deg               # (Qt's angles: anticlockwise from east; azimuth: clockwise from north)
        outer = QtCore.QRectF(cx - r_out, cy - r_out, 2 * r_out, 2 * r_out)
        path = QtGui.QPainterPath()
        path.arcMoveTo(outer, a0)
        path.arcTo(outer, a0, -deg)
        if r_in > 0.5:
            path.arcTo(QtCore.QRectF(cx - r_in, cy - r_in, 2 * r_in, 2 * r_in), a0 - deg, deg)
        else:
            path.lineTo(cx, cy)
        path.closeSubpath()
        return path

    def paintEvent(self, _):
        p = QtGui.QPainter(self)
        p.drawPixmap(0, 0, self._heard_layer())
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        w, h, R, cx, cy = self._geometry()
        font = p.font()
        font.setPointSize(8)
        p.setFont(font)
        import iridium_tracker as it        # (the map's satellite drawing and colours)
        # the selected satellite's pass: rises, highest, sets
        if self.show_path and len(self.path) > 1:
            path = QtGui.QPainterPath(self._xy(self.path[0][0], self.path[0][1], cx, cy, R))
            for el, az, _ in self.path[1:]:
                path.lineTo(self._xy(el, az, cx, cy, R))
            p.setPen(QtGui.QPen(it.C_SELECT, 2, Qt.DashLine))
            p.setBrush(Qt.NoBrush)
            p.drawPath(path)
            p.setPen(QtGui.QColor(TEXT))
            for el, az, t in (self.path[0], self.path[-1]):
                q = self._xy(el, az, cx, cy, R)
                p.drawText(QPointF(q.x() + 6, q.y() + 14), time.strftime('%I:%M %p', time.localtime(t)).lstrip('0'))
        # the satellites up now
        self.hits = []
        if self.show_now:
            for s in sorted(self.now, key=lambda s: s['selected']):
                q = self._xy(s['el'], s['az'], cx, cy, R)
                ring = it.heard_colour(s.get('heard_age'), s.get('sig_q'))
                if ring is not None:
                    p.setPen(QtGui.QPen(ring, 2))
                    p.setBrush(Qt.NoBrush)
                    p.drawEllipse(q, 14, 14)
                elif s.get('heard_nt'):         # heard by our own search, not tracked: a hollow dashed ring
                    p.setPen(QtGui.QPen(QtGui.QColor(225, 230, 240), 1.4, Qt.DashLine))
                    p.setBrush(Qt.NoBrush)
                    p.drawEllipse(q, 14, 14)
                it.draw_satellite(p, q.x(), q.y(), QtGui.QColor(40, 190, 90), 1.1 if s['selected'] else 0.85, s['system'])
                if s['selected']:
                    it.draw_selection(p, q.x(), q.y(), 16.0)
                p.setPen(QtGui.QColor(TEXT))
                p.drawText(QPointF(q.x() + 13, q.y() - 7), s['label'])
                self.hits.append((q.x(), q.y(), s['name']))


def summary(fixes, cells):
    """a few lines on what the plot shows (the receptions counted by sky cell, at each cell's middle)"""
    el_c, az_c = np.meshgrid((np.arange(N_EL) + 0.5) * CELL_DEG, (np.arange(N_AZ) + 0.5) * CELL_DEG, indexing='ij')
    els = np.concatenate([np.array([f[0] for f in fixes], float), el_c.ravel()])
    azs = np.concatenate([np.array([f[1] for f in fixes], float), az_c.ravel()])
    wts = np.concatenate([np.ones(len(fixes)), np.asarray(cells, float).ravel()])
    total = wts.sum()
    if not total:
        return 'Nothing heard in this period yet.'
    sector = (((azs + 22.5) % 360) // 45).astype(int)
    counts = {s: wts[sector == k].sum() for k, s in enumerate(SECTORS)}
    dead = [s for s in SECTORS if counts[s] < 0.02 * total]
    bands = [('overhead (60-90\u00b0)', els >= 60), ('high (30-60\u00b0)', (els >= 30) & (els < 60)),
             ('low (0-30\u00b0)', els < 30)]
    share = {n: wts[m].sum() / total for n, m in bands}
    lines = ['<b>%d reported positions</b> and <b>%d other receptions</b> from identified satellites.'
             % (len(fixes), int(np.asarray(cells).sum())),
             'By height: ' + ', '.join('%s %d%%' % (n, round(100 * share[n])) for n, _ in bands) + '.',
             'By direction: ' + ', '.join('%s %d%%' % (s, round(100 * counts[s] / total)) for s in SECTORS) + '.']
    if dead:
        lines.append('<span style="color:#ff6b6b">Little or nothing heard from: %s.</span>' % ', '.join(dead))
    if share[bands[0][0]] < 0.05:
        lines.append('<span style="color:#ff6b6b">Almost nothing from overhead</span> - typical of a vertical whip '
                     '(blind straight up), or something above the antenna. An upward-facing patch antenna hears '
                     'overhead satellites best - and those are the closest and strongest.')
    return '<br>'.join(lines)


class SkyDialog(QtWidgets.QDialog):
    """the Sky window: the plot (satellites up now, the selected one's pass, and where you have heard them),
    what to show, the satellites up now, and a summary; refreshes while open"""
    PERIODS = [('This run', None), ('Last hour', 3600), ('Last 24 hours', 86400), ('Everything', 0)]

    def __init__(self, parent, data, from_home, current_run_start, shows=None):
        super().__init__(None)      # its own window: taskbar button, minimize/maximize (closed with the main window)
        self.setWindowFlags(Qt.Window | Qt.WindowMinMaxButtonsHint | Qt.WindowCloseButtonHint)
        self.setWindowTitle('Sky plot - where your antenna hears satellites')
        self.data, self.from_home, self.run_start, self.shows = data, from_home, current_run_start, shows
        self.tracker = parent
        self.plot = SkyPlot()
        self.plot.satelliteClicked.connect(parent.select)
        self.heard = HeardCells()
        self.signal = SignalCells()
        self.period = QtWidgets.QComboBox()
        for label, _ in self.PERIODS:
            self.period.addItem(label)
        self.period.setCurrentIndex(0)          # this run: 'Everything' is a choice, not the default it was (slow)
        self.period.currentIndexChanged.connect(self.refresh)
        self.text = QtWidgets.QLabel()
        self.text.setWordWrap(True)
        self.text.setTextFormat(Qt.RichText)
        top = QtWidgets.QHBoxLayout()
        self.layers = []
        for label, attr, tip in (('Satellites now', 'show_now', 'the satellites above your horizon now'),
                                 ('Selected pass', 'show_path', "the selected satellite's path across your sky"),
                                 ('Heard', 'show_heard', 'where you heard satellites (fixes and frames)'),
                                 ('Signal', 'show_signal', 'how strong GNSS signals are in each part of your sky '
                                                           '(middle C/N0 per 5\u00b0 patch, from each satellite\'s own '
                                                           'orbit) - in place of the Heard shading')):
            b = QtWidgets.QPushButton(label)
            b.setCheckable(True)
            b.setChecked(attr != 'show_signal')
            bold = b.font()
            bold.setBold(True)                  # (checked buttons are drawn bold: leave room for that)
            b.setMinimumWidth(QtGui.QFontMetrics(bold).width(label) + 34)
            b.setToolTip('Show ' + tip)
            b.toggled.connect(lambda on, a=attr: (setattr(self.plot, a, on), self.plot.update(),
                                                  self.refresh() if a == 'show_signal' else None))
            top.addWidget(b)
            self.layers.append(b)
        top.addSpacing(10)
        top.addWidget(QtWidgets.QLabel('Heard in:'))
        top.addWidget(self.period)
        top.addStretch(1)
        note = QtWidgets.QLabel('<span style="color:%s">Centre = straight overhead, edge = horizon, north up, east right. '
                                'Satellites: up now (click one to select it; ring = signal strength, red weak to '
                                'green strong, fading as it goes quiet). '
                                'Big dots: positions the satellites reported (colour = signal). Shading: everything '
                                'else heard from identified satellites. Signal: how strong GNSS signals are in each '
                                'part of the sky (weak blue to strong red). Dashed ring: heard by our own search but '
                                'not tracked. Move the antenna and watch the empty areas fill in.</span>' % MUTED)
        note.setWordWrap(True)
        self.up = QtWidgets.QTableWidget(0, 4)          # the satellites up now
        self.up.setHorizontalHeaderLabels(['Up now', 'Height', 'Direction', 'Heard'])
        self.up.verticalHeader().hide()
        self.up.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.up.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.up.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        self.up.horizontalHeader().setStretchLastSection(True)
        self.up.setMinimumWidth(270)
        self.up.cellClicked.connect(lambda r, c: parent.select(self.up.item(r, 0).data(Qt.UserRole)))
        self.up.setSelectionMode(QtWidgets.QAbstractItemView.NoSelection)    # (the selected satellite is coloured)
        self.up.setFocusPolicy(Qt.NoFocus)
        self._up_rows = None        # what the table shows now (update_now changes only what differs)
        mid = QtWidgets.QHBoxLayout()
        mid.addWidget(self.plot, 1)
        mid.addWidget(self.up)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addLayout(top)
        lay.addLayout(mid, 1)
        lay.addWidget(self.text)
        lay.addWidget(note)
        self.resize(900, 760)
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._tick)
        self.timer.start(5000)
        self.refresh()

    def _tick(self):
        """the 5 s refresh - only while the sky is actually on screen (as a view or its own window): once opened,
        it used to keep recomputing every 5 s behind the map, freezing the whole tracker each time"""
        if self.isVisible():
            self.refresh()

    def showEvent(self, event):
        super().showEvent(event)
        self.refresh()

    def refresh(self):
        span = self.PERIODS[self.period.currentIndex()][1]
        if span is None:
            since = self.run_start() or time.time()
        else:
            since = time.time() - span if span else 0.0
        shows = self.shows or (lambda name, cap=None: True)
        fixes = sky_fixes(self.data, self.from_home, since, shows)
        self.heard.update(self.data, self.from_home, shows)
        cells = self.heard.counts(since, shows)
        self.plot.set_points(fixes, cells)
        extra = ''
        try:        # the signal layer and lock table: an addition - never a reason for the sky view to fail
            sig = None
            if self.plot.show_signal:
                self.signal.update(self.data, self.from_home, shows)
                sig = signal_cells(self.signal.histogram(self.data, since, shows))
            self.plot.set_signal(sig)
            extra = signal_summary(sig[0] if sig else None, sig[1] if sig else None,
                                   lock_table(self.data, since, shows), self._heard_nt())
        except Exception as e:
            print('sky view: signal layer: %r' % e)
        self.text.setText(summary(fixes, cells) + ('<br>' + extra if extra else ''))
        self.update_now()

    def _heard_nt(self):
        """[(name, GHEARD info)] heard by gnss_rx's own search and not tracked now (newest first)"""
        got = getattr(self.data, 'heard_not_tracked', None)
        if got is None:
            return []
        import catalog
        return sorted(((catalog.gnss_name(s), h) for s, h in got(time.time()).items()), key=lambda x: -x[1]['t'])

    def update_now(self):
        """the satellites up now and the selected one's pass, from the tracker's current picture (every second)"""
        tr = self.tracker
        sc = getattr(tr.map, 'scene', None) or {}
        import systems
        import iridium_tracker as it
        hnt = dict(self._heard_nt())
        sats = [dict(el=s['el'], az=s['az'], name=s['name'], selected=s['selected'], live=s['live'],
                     heard_age=s.get('heard_age'), sig_q=s.get('sig_q'), heard_nt=hnt.get(s['name']),
                     system=s['system'], label=systems.short_name(s['name']))
                for s in sc.get('sats', []) if s['el'] > 0]
        path = []
        sel = tr.selected
        orbit = tr.data.orbits.get(sel) if sel else None
        np_ = tr.next_pass(sel) if orbit is not None else None
        if np_:
            ts = np.linspace(np_[0], np_[1], 80)
            el, az, _ = self.from_home(*orbit.ecef(ts).T)
            path = [(float(e), float(a), float(t)) for e, a, t in zip(el, az, ts) if e >= 0]
        self.plot.set_now(sats, path)
        sats.sort(key=lambda s: -s['el'])
        rows = []
        for s in sats:
            heard = 'now' if s['live'] else ''
            info = tr.data.info.get(s['name'], {}).get('now')
            if info and s['live'] and 'cn0' in info:
                heard = '%.0f dB-Hz' % info['cn0']
            elif s.get('heard_nt') and not s['live']:
                heard = 'heard, not tracked (~%.0f)' % s['heard_nt']['cn0']
            rows.append((s['name'], s['selected'], [s['label'], '%.0f\u00b0' % s['el'], '%s %.0f\u00b0' % (
                SECTORS[int(((s['az'] + 22.5) % 360) // 45)], s['az']), heard]))
        if rows == self._up_rows:
            return
        # only the cells that changed are touched: rebuilding the table with new items took ~1 ms a cell (its
        # columns are sized to their contents, measured again at every one) - 98 ms of each second's update
        # with 25 satellites up (week-2 sky_load.py)
        self._up_rows = rows
        self.up.setUpdatesEnabled(False)
        try:
            self.up.setRowCount(len(rows))
            for r, (name, selected, texts) in enumerate(rows):
                for c, v in enumerate(texts):
                    item = self.up.item(r, c)
                    if item is None:
                        item = QtWidgets.QTableWidgetItem(v)
                        self.up.setItem(r, c, item)
                    elif item.text() != v:
                        item.setText(v)
                    if c == 0 and item.data(Qt.UserRole) != name:
                        item.setData(Qt.UserRole, name)
                    if item.font().bold() != selected:
                        f = item.font()
                        f.setBold(selected)
                        item.setFont(f)
                    if selected:            # the main list's selection colour
                        item.setBackground(it.C_SELECTED_NAME)
                        item.setForeground(QtGui.QColor('#ffffff'))
                    elif item.background().style() != Qt.NoBrush:
                        item.setBackground(QtGui.QBrush())
                        item.setForeground(QtGui.QBrush())
        finally:
            self.up.setUpdatesEnabled(True)

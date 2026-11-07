"""'Heard but not tracked': our own search for GPS satellites that pocket_trk is not tracking.

pocket_trk only locks a satellite whose signal clears its acquisition threshold (~35 dB-Hz with its search settings),
so satellites at 25-33 dB-Hz - low ones, or ones behind a roof edge - can be there for an hour without ever being
tracked (the review of 2026-11-06 found G27 and G10 at 28-33 dB-Hz in one second of a recording, never tracked).
This module looks for exactly those, cheaply, because it knows where to look: gnss_rx gives it, once a minute, one
second of raw samples and the GPS satellites that are up by their own broadcast orbits (the saved ephemerides) but
not tracked, each with its predicted Doppler plus this SDR's own frequency offset (measured from the tracked
satellites). It searches only +-400 Hz around that, much longer than pocket_trk does (10 ms coherent x 100
non-coherent = 1 s, two 5 ms-apart bit phases), so a satellite ~10 dB weaker than pocket_trk can lock still stands
out. gnss_rx writes a GHEARD line for each one found:
  GHEARD: t sat cn0_est doppler_hz elevation_deg peak_ratio
It is a statement that the satellite's signal is reaching the antenna - never a position, never used for one.

Measured offline (build/_tmp/later/sky/acq_test.py, hook_test.py; four 40 s recordings of 2026-11-06):
every tracked GPS satellite searched the same way was found (113 of 113; C/N0 within 0.6 dB sd of pocket_trk's after
CN0_CAL_DB, Doppler within 25 Hz); untracked satellites found at 24-33 dB-Hz (G27/G10/G24 of the review recording,
G19/G13/G29/G06/G25 later), and during pocket_trk's lock-loss events the satellites it had just dropped were still
there at 35-39 dB-Hz; GPS satellites below the horizon: no detection in 453 searches (highest ratio 2.36, see
RATIO_MIN). 0.15 s of CPU per satellite searched: ~1 s a minute with 5-7 candidates, ~0.02 core.

It runs as its own low-priority process (python acqdetect.py --serve, started by gnss_rx), so neither the radio's
thread nor pocket_trk's reader can ever wait for it. Requests on stdin: a JSON line
  {"id": n, "n": samples, "rate": Hz, "cands": [[sat, prn, doppler_hz], ...], "refs": [[sat, doppler_hz, cn0], ...]}
followed by 2 * samples bytes of 8-bit I, Q as the SDR gives them (gr-osmosdr's convention: a satellite's carrier
sits at +Doppler as pocket_trk reports it - pocket_trk itself is fed Q negated and reports the mirror of its input,
checked on a recording: G12 at +180 Hz by pocket_trk was found at +180 Hz here, at -200 Hz in the Q-negated samples).
Answers on stdout:
  HIT sat cn0 doppler_hz ratio code_sample      (one per satellite found)
  DONE id cpu_s candidates
"""
import json
import math
import sys
import time

import numpy as np

L1_HZ = 1575.42e6
CHIP_HZ = 1.023e6
COH_MS = 10                 # coherent integration: a GPS data bit is 20 ms, so 10 ms blocks at two phases 5 ms apart
NONCOH = 100                # ...summed in power 100 times: one second
SNAP_MS = COH_MS * NONCOH + 5
DOP_HALF_HZ = 400.0         # searched around the prediction (orbit + measured SDR offset: both good to ~10 Hz)
DOP_STEP_HZ = 50.0          # (a 10 ms block's frequency response is 100 Hz wide: 50 Hz steps lose at most 0.9 dB)
# The detection: the strongest cell of the search (2 bit phases x 17 frequencies x 2048 code positions) against the
# median of its frequency's sum. White noise alone would stay under ~1.5 (each cell a sum of 100 exponentials, sd 0.1),
# but the real band is not white (the receiver's own spurs, the strong satellites' and WAAS codes leaking in): searches
# of GPS satellites 5+ deg BELOW the horizon on four recordings (453 searches, 2026-11-06, acq_test.py) peaked at up to
# 2.36 (99 % under 2.10). 3.0 leaves a clear margin; it is ~24 dB-Hz.
RATIO_MIN = 3.0
CN0_CAL_DB = 1.4            # added to 10 log10((ratio - 1) / 10 ms): tracked satellites came out 1.4 dB under pocket_trk (sd 0.6)
# A strong satellite's code leaks into another's search (C/A cross-correlation ~ -21 to -24 dB, at Doppler
# differences near a multiple of 1 kHz): a 'find' that weak, that close to such a line, is not trusted
XCORR_DB = 18.0
XCORR_HZ = 60.0
FREQ_CHUNK = 16             # frequencies folded and transformed together (~26 MB at a time)

# C/A code phase selectors (G2 taps) of PRN 1-32: IS-GPS-200 table 3-Ia
G2_TAPS = {1: (2, 6), 2: (3, 7), 3: (4, 8), 4: (5, 9), 5: (1, 9), 6: (2, 10), 7: (1, 8), 8: (2, 9), 9: (3, 10),
           10: (2, 3), 11: (3, 4), 12: (5, 6), 13: (6, 7), 14: (7, 8), 15: (8, 9), 16: (9, 10), 17: (1, 4),
           18: (2, 5), 19: (3, 6), 20: (4, 7), 21: (5, 8), 22: (6, 9), 23: (1, 3), 24: (4, 6), 25: (5, 7),
           26: (6, 8), 27: (7, 9), 28: (8, 10), 29: (1, 6), 30: (2, 7), 31: (3, 8), 32: (4, 9)}


def ca_code(prn):
    """the 1023-chip C/A code of GPS PRN 1-32 as +-1"""
    a, b = G2_TAPS[prn]
    g1, g2 = [1] * 10, [1] * 10
    out = np.empty(1023, np.int8)
    for i in range(1023):
        out[i] = g1[9] ^ g2[a - 1] ^ g2[b - 1]
        g1 = [g1[2] ^ g1[9]] + g1[:9]
        g2 = [g2[1] ^ g2[2] ^ g2[5] ^ g2[7] ^ g2[8] ^ g2[9]] + g2[:9]
    return 1 - 2 * out.astype(np.int16)


class Searcher:
    """the search itself, for one sample rate (a whole number of samples per millisecond)"""

    def __init__(self, rate):
        self.rate = float(rate)
        self.n = int(round(rate / 1000.0))
        if abs(self.n * 1000.0 - rate) > 1e-6:
            raise ValueError('the search needs a whole number of samples per ms (rate %.0f)' % rate)
        self._codes = {}
        self.k = np.arange(self.n)
        try:        # (imported here, not by gnss_rx, which only needs SNAP_MS and predict from this module)
            from scipy import fft       # single precision kept, ~20 % quicker than numpy's here (0.9 vs 1.3 s / 5 sats)
            self.fft = fft
        except ImportError:
            self.fft = np.fft

    def code_fft(self, prn):
        if prn not in self._codes:
            idx = (np.arange(self.n) * CHIP_HZ / self.rate).astype(int) % 1023
            self._codes[prn] = np.conj(self.fft.fft(ca_code(prn)[idx].astype(np.complex64))).astype(np.complex64)
        return self._codes[prn]

    def search(self, x, cands, refs=()):
        """x: complex samples (gr-osmosdr's convention), at least SNAP_MS ms; cands: [(sat, prn, doppler_hz)] -> found:
        [{'sat', 'cn0', 'dop', 'ratio', 'code'}]. Every candidate's frequencies are searched together, so a frequency
        two candidates share is folded and transformed once."""
        n = self.n
        need = SNAP_MS * n
        if len(x) < need:
            raise ValueError('need %d samples, got %d' % (need, len(x)))
        steps = np.arange(-DOP_HALF_HZ, DOP_HALF_HZ + 1e-6, DOP_STEP_HZ)
        want = {}                       # frequency (Hz, on the 50 Hz grid) -> candidates searching it
        for i, (sat, prn, f0) in enumerate(cands):
            if prn not in G2_TAPS:
                continue
            for f in np.round((f0 + steps) / DOP_STEP_HZ) * DOP_STEP_HZ:
                want.setdefault(float(f), []).append(i)
        best = {}                       # candidate -> (ratio, f, code sample)
        m_ph = np.arange(COH_MS) * 1e-3
        g_t = np.arange(NONCOH) * COH_MS * 1e-3
        freqs = sorted(want)
        for phase in (0, COH_MS // 2):
            seg = np.asarray(x[phase * n:phase * n + NONCOH * COH_MS * n], np.complex64).reshape(NONCOH, COH_MS, n)
            # (ms within the block, block x sample): folding a 10 ms block into one ms is then one matrix product
            seg_t = np.ascontiguousarray(seg.transpose(1, 0, 2)).reshape(COH_MS, NONCOH * n)
            for c0 in range(0, len(freqs), FREQ_CHUNK):
                fs_ = np.array(freqs[c0:c0 + FREQ_CHUNK])
                # carrier off: the per-ms phase steps fold each 10 ms block into one ms (the code repeats every ms;
                # within 10 ms the code's own Doppler moves it 0.03 samples), then the phase within the ms
                W = np.exp(-2j * np.pi * fs_[:, None] * m_ph[None, :]).astype(np.complex64)
                folds = (W @ seg_t).reshape(len(fs_), NONCOH, n)
                folds *= np.exp(-2j * np.pi * fs_[:, None] * self.k[None, :] / self.rate).astype(np.complex64)[:, None, :]
                Fs = self.fft.fft(folds, axis=2)
                for j, f in enumerate(fs_):
                    # over the second the code itself slides with the Doppler (and the SDR's clock error, which shows
                    # in the same measured frequency): ~4 samples at 3 kHz - realigned block by block before summing
                    d = np.round(-f * (g_t + phase * 1e-3) * self.rate / L1_HZ).astype(int)
                    cut = np.r_[0, np.flatnonzero(np.diff(d)) + 1, NONCOH]
                    for i in want[float(f)]:
                        z = self.fft.ifft(Fs[j] * self.code_fft(cands[i][1]), axis=1)
                        P = z.real * z.real + z.imag * z.imag
                        acc = np.zeros(n, np.float64)
                        for a, b in zip(cut[:-1], cut[1:]):
                            acc += np.roll(P[a:b].sum(axis=0), -int(d[a]))
                        k = int(np.argmax(acc))
                        med = float(np.median(acc))
                        r = float(acc[k] / med) if med > 0 else 0.0
                        if i not in best or r > best[i][0]:
                            best[i] = (r, float(f), k)
        found = []
        for i, (r, f, k) in best.items():
            if r < RATIO_MIN:
                continue
            cn0 = 10 * math.log10(max(r - 1, 1e-6) / (COH_MS * 1e-3)) + CN0_CAL_DB
            leak = [s for s, fr, c in refs if c - cn0 >= XCORR_DB and
                    abs((f - fr + 500) % 1000 - 500) <= XCORR_HZ]
            if leak:
                continue
            found.append({'sat': cands[i][0], 'cn0': round(cn0, 1), 'dop': f, 'ratio': round(r, 2), 'code': k})
        return found

    def all_ratios(self, x, cands):
        """for tests: every candidate's best (ratio, doppler, code sample), detected or not"""
        keep = globals()["RATIO_MIN"]
        try:
            globals()['RATIO_MIN'] = -1.0
            out = self.search(x, cands)
        finally:
            globals()['RATIO_MIN'] = keep
        return out


def predict(ephs, home, t):
    """GPS satellites' elevation (deg) and Doppler (Hz, approaching = positive, as pocket_trk reports it) at home
    (lat, lon deg, height km) and unix time t, from broadcast ephemerides (gnss.parse_eph dicts): {sat: (el, dop)}"""
    import gnss
    lat, lon = math.radians(home[0]), math.radians(home[1])
    f = 1 / 298.257223563
    e2 = f * (2 - f)
    N = 6378137.0 / math.sqrt(1 - e2 * math.sin(lat) ** 2)
    h = home[2] * 1000.0
    rx = np.array([(N + h) * math.cos(lat) * math.cos(lon), (N + h) * math.cos(lat) * math.sin(lon),
                   (N * (1 - e2) + h) * math.sin(lat)])
    up = np.array([math.cos(lat) * math.cos(lon), math.cos(lat) * math.sin(lon), math.sin(lat)])
    out = {}
    for e in ephs:
        try:
            p, v = gnss.sat_ecef(e, float(t))
        except Exception:
            continue
        d = np.asarray(p, float) * 1000.0 - rx
        r = float(np.linalg.norm(d))
        el = math.degrees(math.asin(float(np.dot(d, up)) / r))
        rr = float(np.dot(np.asarray(v, float) * 1000.0, d / r))
        out[e['sat']] = (el, -rr * L1_HZ / 299792458.0)
    return out


def serve():
    """the low-priority helper process: requests on stdin, answers on stdout, until stdin closes (gnss_rx gone)"""
    inp, out = sys.stdin.buffer, sys.stdout
    searchers = {}
    while True:
        head = inp.readline()
        if not head:
            return
        try:
            req = json.loads(head)
            n = int(req['n'])
            raw = inp.read(2 * n)
            if len(raw) < 2 * n:
                return
            c0 = time.process_time()
            iq = np.frombuffer(raw, np.int8).astype(np.float32)
            x = (iq[0::2] + 1j * iq[1::2]).astype(np.complex64)
            rate = float(req['rate'])
            s = searchers.get(rate) or searchers.setdefault(rate, Searcher(rate))
            found = s.search(x, [tuple(c) for c in req.get('cands', [])], [tuple(r) for r in req.get('refs', [])])
            for h in found:
                out.write('HIT %s %.1f %.1f %.2f %d\n' % (h['sat'], h['cn0'], h['dop'], h['ratio'], h['code']))
            out.write('DONE %s %.3f %d\n' % (req.get('id', 0), time.process_time() - c0, len(req.get('cands', []))))
            out.flush()
        except (ValueError, KeyError, TypeError) as e:
            out.write('ERROR %s\n' % str(e).replace('\n', ' ')[:200])
            out.write('DONE %s 0 0\n' % (req.get('id', 0) if isinstance(req, dict) else 0))
            out.flush()


if __name__ == '__main__':
    if '--serve' in sys.argv:
        try:
            serve()
        except (BrokenPipeError, OSError, KeyboardInterrupt):
            pass

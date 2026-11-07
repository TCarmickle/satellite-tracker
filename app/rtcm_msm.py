"""Raw observations (pseudorange, carrier phase, Doppler, C/N0) of pocket_trk, made fit for a RINEX file so the
antenna position can be worked out precisely offline (tools/trk2rinex.py writes the files, rinex.py the format).

Two sources carry the same numbers (checked: they agree to the last printed digit, 1 mm / 0.001 cycle / 0.001 Hz):
  * the RTCM 3 stream (`<capture>.trkN.rtcm`): MSM7 messages 1077 (GPS), 1097 (Galileo), 1107 (SBAS / WAAS), one
    set per second while pocket_trk's own position (PVT) is valid, plus the ephemerides 1019 (GPS) / 1046
    (Galileo). Field layout RTCM 10403.3 as RTKLIB's rtcm3e.c writes it (PocketSDR uses it); the frame splitter and
    bit reader are rtcm3.py's.
  * the `$OBS` lines of the log (`<capture>.trkN.log`), the same epochs.

What has to be done to them before a precise solver can use them (measured on the 2026-11-06 evening captures,
scratch scripts in build\\_tmp\\later\\rinex\\):
  1. ONE receiver clock for code and carrier. pocket_trk's pseudoranges carry the clock its time tags are on;
     its carrier phase runs on the tuner's own frequency, which is off by ~88 Hz (the E4000's synthesiser cannot
     hit 1575.42 MHz exactly), so code minus carrier grows a steady +16.73 m every second, the same for every
     satellite (16.71-16.75 m/s in every 10 minutes of three runs). On top of that pocket_trk now and then shifts
     its measurement instant by 20 or 40 ms ('PVT EPOCH ADJUSTED' in the log): every pseudorange jumps by exactly
     that much light time, the carrier does not. `consistent()` measures the rate (median over the run) and the
     jumps (median over satellites of code minus Doppler-predicted range change) and moves the CARRIER and the
     Doppler onto the code's clock. The code is left as it is because the time tags belong to its clock: moving
     the code instead would leave tags and ranges up to 40 ms apart (up to ~30 m of range error).
  2. Loss of lock from the carrier itself. The MSM lock-time indicator restarts every ~72 s for every satellite,
     because the encoder re-anchors the MSM phase to the pseudorange by a whole number of cycles whenever code
     minus carrier passes the field's 1171 m limit - so it is not a loss-of-lock flag, and those re-anchoring jumps
     are taken back out here (the jump is a whole number of cycles, found to ~0.1 cycle). A slip is a phase that
     leaves the smooth track of its own last three seconds (third difference, minus the part common to all
     satellites, which is the receiver's oscillator): 0.07 cycle typical, 99 % under 0.18, slip threshold
     SLIP_CYC. The first three seconds of an arc are checked against the Doppler instead (its noise is ~0.6 cycle
     a second, so only a coarse check). A missing second, an epoch shift (above) and a slip start a new arc:
     RINEX loss-of-lock indicator bit 0.
"""
import calendar
import statistics

import rtcm3

C = 299792458.0
F_L1 = 1575.42e6                 # GPS L1, Galileo E1 and SBAS L1 share the frequency
LAM_L1 = C / F_L1
RANGE_MS = C * 1e-3              # one millisecond of light, metres

MSM7 = {1077: 'G', 1097: 'E', 1107: 'S'}
EPHEMERIS = (1019, 1046)

# MSM signal id -> RINEX 3 observation code (RTCM 10403.3 tables 3.5-91, -100 and -102; pocket_trk sends only the
# L1 ones, the rest keep a correct name so nothing is mislabelled if it ever sends more)
SIGNALS = {'G': {2: '1C', 3: '1P', 4: '1W', 8: '2C', 9: '2P', 10: '2W', 15: '2S', 16: '2L', 17: '2X',
                 22: '5I', 23: '5Q', 24: '5X', 30: '1S', 31: '1L', 32: '1X'},
           'E': {2: '1C', 3: '1A', 4: '1B', 5: '1X', 6: '1Z', 14: '7I', 15: '7Q', 16: '7X', 22: '5I', 23: '5Q',
                 24: '5X'},
           'S': {2: '1C', 22: '5I', 23: '5Q', 24: '5X'}}

SLIP_CYC = 0.5          # third-difference residual (cycles) beyond which the phase slipped (a full slip shows 1.0)
DOPPLER_SLIP_CYC = 5.0  # phase vs Doppler over one second in an arc's first seconds (Doppler noise reaches ~3.7)
REANCHOR_TOL = 0.25     # an MSM re-anchoring jump is a whole number of cycles to within this
MAX_GAP_S = 1.5         # epochs further apart than this (s) break every arc
JUMP_M = 1000.0         # a common code jump bigger than this (m) is a shift of the measurement instant
MIN_P = 1e7             # shorter pseudoranges are not real (spp.py uses the same test)


def message_type(payload):
    return (payload[0] << 4) | (payload[1] >> 4) if len(payload) >= 2 else 0


def decode_msm7(payload):
    """one MSM7 payload -> {'msg', 'sys', 'tow_ms', 'more', 'steer', 'obs': [...]} or None.
    Each observation: {'sat', 'code' (RINEX, e.g. '1C'), 'P' m, 'L' cycles, 'D' Hz, 'cn0' dB-Hz, 'lock' (the MSM
    lock-time indicator, raw), 'half' (1 = half-cycle ambiguity not resolved)}; P/L/D are None where the message
    marks them invalid."""
    msg = message_type(payload)
    sys = MSM7.get(msg)
    if sys is None:
        return None
    b = rtcm3._Bits(payload)
    b.u(12 + 12)                                    # message number, station id
    tow_ms = b.u(30)                                # GPS / Galileo time of week, ms
    more = b.u(1)                                   # 1: more MSM messages for this epoch follow
    b.u(3 + 7)                                      # IODS, reserved
    steer = b.u(2)                                  # clock steering indicator (0 none, 1 applied, 2 unknown)
    b.u(2 + 1 + 3)                                  # external clock, smoothing, smoothing interval
    satmask, sigmask = b.u(64), b.u(32)
    sats = [i + 1 for i in range(64) if satmask >> (63 - i) & 1]
    sigs = [i + 1 for i in range(32) if sigmask >> (31 - i) & 1]
    ncell = len(sats) * len(sigs)
    if ncell > 64:                                  # not a valid MSM (the cell mask is at most 64 bits)
        return None
    cellmask = b.u(ncell) if ncell else 0
    cells = [(sats[i // len(sigs)], sigs[i % len(sigs)]) for i in range(ncell) if cellmask >> (ncell - 1 - i) & 1]
    n, m = len(sats), len(cells)
    if len(payload) * 8 < b.p + n * 36 + m * 80:    # truncated message
        return None
    rough = [b.u(8) for _ in range(n)]              # whole milliseconds of range (255 = invalid)
    [b.u(4) for _ in range(n)]                      # extended satellite info (GLONASS channel; unused)
    rough_mod = [b.u(10) for _ in range(n)]         # range modulo 1 ms, 2^-10 ms
    rough_rate = [b.s(14) for _ in range(n)]        # phase range rate, m/s (-8192 = invalid)
    fine_p = [b.s(20) for _ in range(m)]            # 2^-29 ms
    fine_l = [b.s(24) for _ in range(m)]            # 2^-31 ms
    lock = [b.u(10) for _ in range(m)]
    half = [b.u(1) for _ in range(m)]
    cnr = [b.u(10) for _ in range(m)]               # 2^-4 dB-Hz
    fine_rate = [b.s(15) for _ in range(m)]         # 0.0001 m/s (-16384 = invalid)
    obs = []
    for j, (sid, sig) in enumerate(cells):
        i = sats.index(sid)
        code = SIGNALS[sys].get(sig)
        if code is None:
            continue
        # MSM SBAS satellite id 1 is PRN 120; RINEX names SBAS satellites PRN - 100 ('S31' = PRN 131)
        sat = 'S%02d' % (sid + 19) if sys == 'S' else '%s%02d' % (sys, sid)
        P = L = D = None
        if rough[i] != 255:
            r = (rough[i] + rough_mod[i] / 1024.0) * RANGE_MS
            if fine_p[j] != -(1 << 19):
                P = r + fine_p[j] * 2.0 ** -29 * RANGE_MS
            if fine_l[j] != -(1 << 23):
                L = (r + fine_l[j] * 2.0 ** -31 * RANGE_MS) / LAM_L1
        if rough_rate[i] != -8192 and fine_rate[j] != -16384:
            D = -(rough_rate[i] + fine_rate[j] * 1e-4) / LAM_L1      # range growing = Doppler negative
        obs.append({'sat': sat, 'code': code, 'P': P, 'L': L, 'D': D, 'cn0': cnr[j] * 2.0 ** -4,
                    'lock': lock[j], 'half': half[j]})
    return {'msg': msg, 'sys': sys, 'tow_ms': tow_ms, 'more': more, 'steer': steer, 'obs': obs}


def decode_ephemeris(payload):
    """1019 (GPS) / 1046 (Galileo I/NAV) -> the raw broadcast values a RINEX navigation record needs, or None.
    rtcm3.decode turns the same messages into the solver's form (GPST times, A instead of sqrt(A), one group delay);
    RINEX wants them as sent: week, seconds of week, both Galileo group delays and the health bits."""
    msg = message_type(payload)
    if msg not in EPHEMERIS or len(payload) * 8 < (488 if msg == 1019 else 504):
        return None
    b = rtcm3._Bits(payload)
    b.u(12)
    e = {'msg': msg}
    if msg == 1019:
        e['sys'] = 'G'
        e['prn'], e['week10'], e['ura'], e['code'] = b.u(6), b.u(10), b.u(4), b.u(2)
        e['idot'] = b.s(14) * 2.0 ** -43 * rtcm3.SC2RAD
        e['iode'], e['toc'] = b.u(8), b.u(16) * 16.0
        e['af2'], e['af1'], e['af0'] = b.s(8) * 2.0 ** -55, b.s(16) * 2.0 ** -43, b.s(22) * 2.0 ** -31
        e['iodc'] = b.u(10)
    else:
        e['sys'] = 'E'
        e['prn'], e['week'], e['iode'], e['sisa'] = b.u(6), b.u(12) + 1024, b.u(10), b.u(8)   # GST week -> GPS
        e['idot'] = b.s(14) * 2.0 ** -43 * rtcm3.SC2RAD
        e['toc'] = b.u(14) * 60.0
        e['af2'], e['af1'], e['af0'] = b.s(6) * 2.0 ** -59, b.s(21) * 2.0 ** -46, b.s(31) * 2.0 ** -34
    e['crs'], e['deln'] = b.s(16) * 2.0 ** -5, b.s(16) * 2.0 ** -43 * rtcm3.SC2RAD
    e['M0'] = b.s(32) * 2.0 ** -31 * rtcm3.SC2RAD
    e['cuc'], e['e'], e['cus'] = b.s(16) * 2.0 ** -29, b.u(32) * 2.0 ** -33, b.s(16) * 2.0 ** -29
    e['sqrtA'] = b.u(32) * 2.0 ** -19
    e['toe'] = b.u(16) * 16.0 if msg == 1019 else b.u(14) * 60.0
    e['cic'], e['OMG0'] = b.s(16) * 2.0 ** -29, b.s(32) * 2.0 ** -31 * rtcm3.SC2RAD
    e['cis'], e['i0'] = b.s(16) * 2.0 ** -29, b.s(32) * 2.0 ** -31 * rtcm3.SC2RAD
    e['crc'], e['omg'] = b.s(16) * 2.0 ** -5, b.s(32) * 2.0 ** -31 * rtcm3.SC2RAD
    e['OMGd'] = b.s(24) * 2.0 ** -43 * rtcm3.SC2RAD
    if msg == 1019:
        e['tgd'] = b.s(8) * 2.0 ** -31
        e['svh'], e['l2p'], fit = b.u(6), b.u(1), b.u(1)
        e['fit_h'] = 4.0 if fit == 0 else 6.0       # IS-GPS-200 20.3.4.4: flag 1 = more than 4 hours
    else:
        e['bgd_a'], e['bgd_b'] = b.s(10) * 2.0 ** -32, b.s(10) * 2.0 ** -32      # BGD(E1,E5a), BGD(E1,E5b)
        e5b_hs, e5b_dvs, e1b_hs, e1b_dvs = b.u(2), b.u(1), b.u(2), b.u(1)
        # RINEX 3.04 Galileo health: bit 0 E1-B DVS, 1-2 E1-B HS, 6 E5b DVS, 7-8 E5b HS
        e['svh'] = e1b_dvs | (e1b_hs << 1) | (e5b_dvs << 6) | (e5b_hs << 7)
    e['sat'] = '%s%02d' % (e['sys'], e['prn'])
    return e


def ura_m(index):
    """GPS URA index -> metres (IS-GPS-200 20.3.3.3.1.3, the value RINEX carries)"""
    table = [2.4, 3.4, 4.85, 6.85, 9.65, 13.65, 24.0, 48.0, 96.0, 192.0, 384.0, 768.0, 1536.0, 3072.0, 6144.0]
    return table[index] if 0 <= index < len(table) else 6144.0


def sisa_m(index):
    """Galileo SISA index -> metres (OS SIS ICD 5.1.12); -1 = no accuracy prediction available"""
    if index < 50:
        return index * 0.01
    if index < 75:
        return 0.5 + (index - 50) * 0.02
    if index < 100:
        return 1.0 + (index - 75) * 0.04
    if index < 126:
        return 2.0 + (index - 100) * 0.16
    return -1.0


def full_week(week_mod, modulo, ref_week):
    """a week number counted modulo `modulo` -> the full week nearest ref_week"""
    return ref_week + ((week_mod - ref_week + modulo // 2) % modulo) - modulo // 2


def read_rtcm(buf, ref_gpst):
    """an RTCM 3 byte stream -> (epochs, ephemerides).
    epochs: [(t_gpst, {sat: obs})], the MSM7 of one epoch merged (GPS, Galileo and SBAS share the time of week);
    t_gpst is RTKLIB-style GPST seconds (counted like Unix time, no leap seconds). ref_gpst: any GPST within half a
    week of the first epoch (the log's first $TIME, or the capture's start time) to fix the week number.
    ephemerides: decode_ephemeris dicts with 'week' (full) and 'ttr' (GPST of the epoch before it arrived - the
    nearest thing to its transmission time the stream holds)."""
    payloads, _ = rtcm3.split(buf)
    raw, ephs = [], []
    for p in payloads:
        mt = message_type(p)
        if mt in MSM7:
            m = decode_msm7(p)
            if m is None:
                continue
            if not raw or raw[-1][0] != m['tow_ms']:
                raw.append((m['tow_ms'], {}))
            for o in m['obs']:
                raw[-1][1][o['sat']] = o
        elif mt in EPHEMERIS:
            e = decode_ephemeris(p)
            if e:
                e['_epoch'] = len(raw) - 1
                ephs.append(e)
    epochs = []
    if raw:
        week = round((ref_gpst - rtcm3.GPS_EPOCH - raw[0][0] / 1000.0) / 604800.0)
        prev = None
        for tow_ms, obs in raw:
            if prev is not None and tow_ms < prev - 302400000:      # Saturday midnight: the week rolled over
                week += 1
            prev = tow_ms
            epochs.append((rtcm3.GPS_EPOCH + week * 604800 + tow_ms / 1000.0, obs))
    ref_week = int((ref_gpst - rtcm3.GPS_EPOCH) // 604800)
    for e in ephs:
        e['ttr'] = epochs[max(e.pop('_epoch'), 0)][0] if epochs else ref_gpst
        if e['sys'] == 'G':
            e['week'] = full_week(e['week10'], 1024, ref_week)
        # the toe can sit in the week before / after the reception (an orbit for 02:00 sent before midnight)
        t_toe = rtcm3.GPS_EPOCH + e['week'] * 604800 + e['toe']
        if t_toe - e['ttr'] > 302400:
            e['week'] -= 1
        elif e['ttr'] - t_toe > 302400:
            e['week'] += 1
    return epochs, ephs


def read_log(path):
    """pocket_trk's log -> (epochs from its $OBS lines in read_rtcm's form, first $TIME as GPST or None,
    [(lat, lon, h)] of its own fixes ($POS))"""
    obs, ref, pos = {}, None, []
    with open(path, 'r', errors='replace') as fh:
        for line in fh:
            # pocket_trk writes from several threads: two records can share a line
            for rec in line.replace('$', '\n$').split('\n'):
                f = rec.strip().split(',')
                try:
                    if f[0] == '$OBS' and len(f) >= 15:
                        t = calendar.timegm((int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), 0)) + float(f[7])
                        obs.setdefault(round(t, 3), {})[f[8]] = {
                            'sat': f[8], 'code': f[9], 'cn0': float(f[10]), 'P': float(f[11]), 'L': float(f[12]),
                            'D': float(f[13]), 'lli': int(f[14]), 'lock': None, 'half': 0}
                    elif f[0] == '$TIME' and ref is None and len(f) >= 8:
                        ref = calendar.timegm((int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6]), 0)) + \
                            float(f[7])
                    elif f[0] == '$POS' and len(f) >= 11:
                        pos.append((float(f[8]), float(f[9]), float(f[10])))
                except (ValueError, IndexError):
                    continue
    return sorted(obs.items()), ref, pos


def _usable(o):
    return o['P'] is not None and o['L'] is not None and o['D'] is not None and o['P'] > MIN_P and o['L'] != 0


def code_carrier_rate(epochs, max_gap_s=MAX_GAP_S):
    """the steady growth of code minus carrier (m/s; the tuner's frequency error), or None if the run is too short.
    Each satellite's code minus carrier is cut into pieces wherever it jumps (a gap, a shifted instant, an MSM
    re-anchoring - all far over 100 m in a second) and a straight line fitted to every piece of a minute or more;
    the answer is the length-weighted median of the slopes. A line over minutes is needed: second-to-second
    differences carry the ~1 m code noise. What remains is the ionosphere's own change (code and carrier part by
    twice it, a few mm/s, different for every satellite), which no single-frequency receiver can separate."""
    pieces, open_ = [], {}
    prev_t = None

    def close(s):
        seg = open_.pop(s)
        if len(seg) >= 30 and seg[-1][0] - seg[0][0] >= 60:
            n = len(seg)
            mt = sum(x[0] for x in seg) / n
            my = sum(x[1] for x in seg) / n
            sxx = sum((x[0] - mt) ** 2 for x in seg)
            pieces.append((sum((x[0] - mt) * (x[1] - my) for x in seg) / sxx, seg[-1][0] - seg[0][0]))

    for t, obs in epochs:
        joined = prev_t is not None and 0 < t - prev_t <= max_gap_s
        for s in list(open_):
            if not joined or s not in obs or not _usable(obs[s]):
                close(s)
        for s, o in obs.items():
            if not _usable(o):
                continue
            y = o['P'] - LAM_L1 * o['L']
            seg = open_.get(s)
            lock = o.get('lock')
            if seg and (abs(y - seg[-1][1]) > 100.0 or o.get('lli', 0) & 1 or
                        (lock is not None and seg[-1][2] is not None and lock < seg[-1][2])):
                close(s)
                seg = None
            if seg is None:
                seg = open_[s] = []
            seg.append((t - epochs[0][0], y, lock))
        prev_t = t
    for s in list(open_):
        close(s)
    if not pieces:
        return None
    pieces.sort()
    half, acc = sum(w for _, w in pieces) / 2.0, 0.0
    for slope, w in pieces:
        acc += w
        if acc >= half:
            return slope
    return pieces[-1][0]


def consistent(epochs, rate=None, slip_cyc=SLIP_CYC, max_gap_s=MAX_GAP_S):
    """[(t, {sat: obs})] from read_rtcm or read_log -> ([(t, {sat: {'code', 'P', 'L', 'D', 'cn0', 'lli'}})], stats)
    with code and carrier on one receiver clock and loss of lock found from the carrier (see the module notes).
    'lli': RINEX loss-of-lock indicator, bit 0 = new arc (slip, gap, first sight), bit 1 = half-cycle ambiguity.
    Each arc's carrier is offset by whole cycles so it starts within a cycle of its pseudorange, as a receiver's
    own RINEX would. Observations missing P, L or D are dropped (a precise solver can't use them)."""
    if rate is None:
        rate = estimate_rate(epochs, max_gap_s) or 0.0
    return _consistent(epochs, rate, slip_cyc, max_gap_s)


def estimate_rate(epochs, max_gap_s=MAX_GAP_S):
    """code_carrier_rate for one run, refined: a first pass repairs the arcs with the rate from the raw pieces,
    then the rate is measured again on the repaired arcs, which are the same whether the numbers came from the log
    or the RTCM (whose re-anchoring cuts its raw pieces short; the two then agree to 0.002 cycle). A run without a
    minute-long arc gets the coarse median of second-to-second changes; None if even that has under 10 seconds -
    the caller should then use another run of the same tuner (16.73-16.74 m/s in every run of 2026-11-06)."""
    rate = code_carrier_rate(epochs, max_gap_s)
    if rate is None:
        per_epoch, prev_t, prev = [], None, {}
        for t, obs in epochs:
            if prev_t is not None and 0 < t - prev_t <= max_gap_s:
                d = [((o['P'] - prev[s]['P']) - LAM_L1 * (o['L'] - prev[s]['L'])) / (t - prev_t)
                     for s, o in obs.items() if s in prev and _usable(o) and _usable(prev[s])]
                if d:
                    per_epoch.append(statistics.median(d))
            prev_t, prev = t, obs
        return statistics.median(per_epoch) if len(per_epoch) >= 10 else None
    first, _ = _consistent(epochs, rate, SLIP_CYC, max_gap_s)
    return rate + (code_carrier_rate(first, max_gap_s) or 0.0)


def _consistent(epochs, rate, slip_cyc, max_gap_s):
    st = {'rate_m_s': rate, 'epochs': 0, 'obs': 0, 'arcs': 0, 'slips': 0, 'slips_doppler': 0, 'half': 0,
          'jumps': [], 'reanchored': 0, 'reanchor_failed': 0, 'lock_resets': 0, 'res3': [], 'gaps': 0}
    corr = 0.0             # metres added to lambda*L to put the carrier on the code's clock
    arc = {}               # sat -> {'off': whole cycles taken off, 'hist': [(t, L unwrapped)], 'lock', 'D'}
    out, prev_t = [], None
    for t, obs in epochs:
        cur = {s: o for s, o in obs.items() if _usable(o)}
        dt = None if prev_t is None else t - prev_t
        joined = dt is not None and 0 < dt <= max_gap_s
        if prev_t is not None and not joined:
            st['gaps'] += 1
        # a shift of pocket_trk's measurement instant: every code range jumps by the same light time while the
        # carrier and Doppler carry on (code minus Doppler-predicted change, median across satellites)
        jump = 0.0
        if joined:
            q = [(o['P'] - arc[s]['P']) + LAM_L1 * 0.5 * (o['D'] + arc[s]['D']) * dt - rate * dt
                 for s, o in cur.items() if s in arc]
            if q and abs(statistics.median(q)) > JUMP_M:
                jump = round(statistics.median(q) / RANGE_MS) * RANGE_MS
                st['jumps'].append((t, jump / RANGE_MS))
            corr += rate * dt + jump
        if not joined or jump:
            arc = {}       # nothing guarantees continuity across a gap or a shifted instant: every arc restarts
        # carrier prediction from each satellite's last three seconds; the part common to all is the oscillator
        res, reset = {}, set()
        for s, o in cur.items():
            a = arc.get(s)
            if a is None:
                continue
            Lu = o['L'] - a['off']
            if o['lock'] is not None and a['lock'] is not None and o['lock'] < a['lock']:
                reset.add(s)
                st['lock_resets'] += 1
            h = a['hist']
            if len(h) == 3 and abs((h[2][0] - h[1][0]) - dt) < 1e-3 and abs((h[1][0] - h[0][0]) - dt) < 1e-3:
                res[s] = ('3rd', Lu - (3 * h[2][1] - 3 * h[1][1] + h[0][1]))
            else:
                res[s] = ('dop', (Lu - h[-1][1]) + 0.5 * (o['D'] + a['D']) * dt)
        third = [r for s, (k, r) in res.items() if k == '3rd' and s not in reset] or \
            [r for s, (k, r) in res.items() if k == '3rd']
        common = statistics.median(third) if len(third) >= 2 else 0.0
        rec = {}
        for s, o in cur.items():
            a = arc.get(s)
            lli = 0
            if a is not None:
                kind, r = res[s]
                e = r - common if kind == '3rd' else r
                if s in reset and abs(e) > slip_cyc:
                    # the MSM encoder re-anchored this satellite's phase by whole cycles: take them back out
                    k = round(e)
                    if kind == '3rd' and abs(e - k) < REANCHOR_TOL:
                        a['off'] += k
                        e -= k
                        st['reanchored'] += 1
                    else:
                        st['reanchor_failed'] += 1
                if kind == '3rd' and abs(e) <= slip_cyc:
                    st['res3'].append(e)
                if abs(e) > (slip_cyc if kind == '3rd' else DOPPLER_SLIP_CYC):
                    st['slips' if kind == '3rd' else 'slips_doppler'] += 1
                    a = None
            if a is None:
                lli = 1
                st['arcs'] += 1
                # whole cycles chosen so this arc's carrier starts within a cycle of its pseudorange
                a = {'off': round(o['L'] + corr / LAM_L1 - o['P'] / LAM_L1), 'hist': []}
                arc[s] = a
            Lu = o['L'] - a['off']
            a['hist'] = (a['hist'] + [(t, Lu)])[-3:]
            a['lock'], a['D'], a['P'] = o['lock'], o['D'], o['P']
            if o['half']:
                lli |= 2
                st['half'] += 1
            rec[s] = {'code': o['code'], 'P': o['P'], 'L': Lu + corr / LAM_L1, 'D': o['D'] - rate / LAM_L1,
                      'cn0': o['cn0'], 'lli': lli}
        for s in [s for s in arc if s not in cur]:
            del arc[s]     # not seen this second: its next sighting is a new arc
        st['epochs'] += 1
        st['obs'] += len(rec)
        out.append((t, rec))
        prev_t = t
    return out, st

"""GNSS bench: run pocket_trk several times at real time on one raw recording and print one table of what it did -
acquisitions, losses, lock lengths, satellites, pocket_trk's own fixes and this program's solver (spp.py) on the
same measurements. For comparing pocket_trk settings (an options file) or what a damaged stream does to it.

    python build/tools/gnss_bench.py RECORDING.cs8 [--runs 3] [--opt my.opt] [--perturb cut:20:120]

RECORDING: signed 8-bit I,Q pairs (.cs8) as Integrity's snapshots and the iq-trigger write them (gr-osmosdr's IQ
convention, so Q is negated here for pocket_trk - pass --q-negated for a file already converted). The .json beside
it gives the rate (else --rate).
--perturb (one or more):
    cut:T:MS     remove MS milliseconds of samples at T seconds (an RTL-SDR/USB loss left as a cut)
    zero:T:MS    replace them with zeros (what gnss_rx's gap zero-fill does with a loss)
    spur:KHZ:DB  add a steady carrier KHZ kHz from the centre, DB dB above the noise (a narrowband spur/jammer)
Runs are one at a time and real time (-tscale 1): faster replays starve pocket_trk's one-channel-at-a-time search,
and its search timing varies between runs, so compare settings over several runs (the median row).
Positions: only spreads about each run's own middle, never coordinates.
Work files go to build/_tmp/bench/<recording>-<label>/ (git-ignored)."""
import argparse
import collections
import json
import math
import os
import re
import subprocess
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
APP = os.path.join(ROOT, 'app')
sys.path.insert(0, APP)
import gnss  # noqa: E402
import rtcm3  # noqa: E402
import spp  # noqa: E402

PT = os.path.join(APP, 'pocketsdr', 'pocket_trk.exe')
# the same satellites the live receiver searches (gnss_rx.SIGNALS)
SIGNALS = ['-sig', 'L1CA', '-prn', '1-32,131,133,135', '-sig', 'E1B', '-prn', '1-36']
# pocket_trk writes from several threads: now and then two records share a line (as gnss_rx._RECORD)
RECORD = re.compile(r'(?<!^)(?=[$](?:CH|OBS|NAV|EPH|POS|TIME|SAT|LOG),)')


def records(path):
    with open(path, errors='replace') as f:
        for line in f:
            line = line.rstrip()
            yield from (RECORD.split(line) if line.count('$') > 1 else (line,))


def prepare(src, out, rate, q_negated, perturb):
    """the recording as pocket_trk reads it (INT8X2, Q negated), with the perturbations applied"""
    m = np.fromfile(src, np.int8)
    if not q_negated:
        m[1::2] = -np.maximum(m[1::2], -127)        # (-128 would overflow)
    notes = []
    for p in perturb:
        kind, a, b = p.split(':')
        a, b = float(a), float(b)
        if kind in ('cut', 'zero'):
            i0, n = int(a * rate) * 2, int(b / 1000.0 * rate) * 2
            if kind == 'cut':
                m = np.concatenate([m[:i0], m[i0 + n:]])
            else:
                m[i0:i0 + n] = 0
            notes.append('%s %.0f ms at %.1f s' % (kind, b, a))
        elif kind == 'spur':
            x = m[0::2].astype(np.float32) + 1j * m[1::2].astype(np.float32)
            rms = float(np.sqrt(np.mean(np.abs(x[:int(rate)]) ** 2)))
            amp = rms * 10 ** (b / 20.0)
            ph = 2 * np.pi * a * 1000.0 / rate * np.arange(len(x))
            x = x + amp * np.exp(1j * ph)
            m = np.empty(2 * len(x), np.int8)
            m[0::2] = np.clip(np.round(x.real), -127, 127)
            m[1::2] = np.clip(np.round(x.imag), -127, 127)
            notes.append('spur %+.0f kHz %+.0f dB' % (a, b))
        else:
            raise SystemExit('unknown perturbation %r' % p)
    m.tofile(out)
    return len(m) // 2 / rate, notes


def run_once(folder, data, rate, opt):
    env = dict(os.environ)
    lib = os.path.join(ROOT, 'python', 'Library', 'bin')      # (SoapySDR.dll: else exit -1073741515)
    env['PATH'] = lib + os.pathsep + env.get('PATH', '')
    cmd = [PT] + SIGNALS + ['-fmt', 'INT8X2', '-f', '%.6f' % (rate / 1e6), '-tscale', '1',
                            '-log', 'trk.log', '-rtcm', 'trk.rtcm']
    if opt:
        cmd += ['-opt', os.path.abspath(opt)]
    cmd.append(data)
    t0 = time.time()
    p = subprocess.run(cmd, cwd=folder, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return time.time() - t0, p.returncode


def analyse(folder, t_ref):
    """one run's numbers from its log and RTCM"""
    found, lost, lock = {}, [], collections.defaultdict(float)
    early_lost = 0
    obs_sats = collections.defaultdict(set)
    pos_t = []
    for rec in records(os.path.join(folder, 'trk.log')):
        f = rec.split(',')
        try:
            if rec.startswith('$LOG') and len(f) > 3:
                t, sat = float(f[1]), f[2][:4] + f[3]
                if 'SIGNAL FOUND' in rec:
                    found.setdefault(sat, []).append(t)
                elif 'SIGNAL LOST' in rec:
                    lost.append((t, sat))
                    if found.get(sat) and t - found[sat][-1] <= 2.5:
                        early_lost += 1
            elif rec.startswith('$CH') and len(f) > 10:
                lock[f[4]] = max(lock[f[4]], float(f[7]))
            elif rec.startswith('$OBS') and len(f) > 8:
                obs_sats[f[1]].add(f[8])
            elif rec.startswith('$POS') and len(f) > 2:
                pos_t.append(float(f[1]))
        except ValueError:
            continue
    # the own solver on the same measurements; every orbit pocket_trk decoded is given at the start (a warm start:
    # the RTCM file carries no arrival times), so its first fix is 'as soon as 4 satellites are measured'
    s = spp.Solver(iono='sbas', sbas_orbit=True, smooth=True)
    try:
        payloads, _ = rtcm3.split(open(os.path.join(folder, 'trk.rtcm'), 'rb').read())
    except OSError:
        payloads = []
    for pl in payloads:
        e = rtcm3.decode(pl, t_ref)
        if e and not e[2].get('glo'):
            sat, sig, fields = e
            e = gnss.parse_eph('$EPH,0,%s,%s' % (sat, ','.join([sig] + [str(fields[k]) for k in gnss.FIELDS])))
            if e:
                s.eph(e)
    own = []
    for rec in records(os.path.join(folder, 'trk.log')):
        try:
            if rec.startswith('$NAV'):
                s.nav(rec)
            elif rec.startswith('$OBS'):
                r = s.obs(rec)
                if r and r['ok']:
                    own.append(r)
        except Exception:
            continue
    n_acq = sum(len(v) for v in found.values())
    gps = [v for k, v in lock.items() if k.startswith('G')]
    gal = [v for k, v in lock.items() if k.startswith('E')]
    row = {'acq': n_acq, 'lost': len(lost), 'lost_2s_pct': 100.0 * early_lost / n_acq if n_acq else None,
           'gps_lock_med_s': float(np.median(gps)) if gps else None,
           'gal_lock_med_s': float(np.median(gal)) if gal else None,
           'max_sats': max((len(v) for v in obs_sats.values()), default=0),
           'pt_first_fix_s': min(pos_t) if pos_t else None, 'pt_fixes': len(pos_t),
           'own_first_fix_s': own[0]['t'] if own else None, 'own_fixes': len(own), 'own_r95_m': None}
    if own and len(own) >= 30:
        lat0, lon0 = math.radians(np.median([r['lat'] for r in own])), math.radians(np.median([r['lon'] for r in own]))
        e = []
        for r in own:
            dn = (math.radians(r['lat']) - lat0) * 6371000.0
            de = (math.radians(r['lon']) - lon0) * 6371000.0 * math.cos(lat0)
            e.append(math.hypot(dn, de))
        row['own_r95_m'] = float(np.percentile(e, 95))
    # the own solver's first-fix time is a receiver time (seconds since the file start), like pocket_trk's
    return row


COLS = [('acq', 'acq', '%4d'), ('lost', 'lost', '%4d'), ('lost_2s_pct', 'lost<2.5s%', '%9.0f'),
        ('gps_lock_med_s', 'GPS lock med s', '%14.0f'), ('gal_lock_med_s', 'Gal lock med s', '%14.0f'),
        ('max_sats', 'max sats', '%8d'), ('pt_first_fix_s', 'pt 1st fix s', '%12.1f'), ('pt_fixes', 'pt fixes', '%8d'),
        ('own_first_fix_s', 'own 1st fix s', '%13.1f'), ('own_fixes', 'own fixes', '%9d'),
        ('own_r95_m', 'own R95 m', '%9.1f')]


def fmt(row):
    out = []
    for k, head, f in COLS:
        v = row.get(k)
        w = len(head)
        out.append(('%' + str(w) + 's') % '-' if v is None else f % v)
    return '  '.join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('recording')
    ap.add_argument('--runs', type=int, default=3)
    ap.add_argument('--opt', default='', help='a pocket_trk options file (lines name = value)')
    ap.add_argument('--perturb', action='append', default=[])
    ap.add_argument('--rate', type=float, default=0.0)
    ap.add_argument('--q-negated', action='store_true', help='the file is already in pocket_trk\'s IQ convention')
    ap.add_argument('--label', default='')
    ap.add_argument('--reuse', action='store_true', help='only analyse the runs already in the work folder')
    a = ap.parse_args()
    meta = {}
    js = os.path.splitext(a.recording)[0] + '.json'
    if os.path.exists(js):
        meta = json.load(open(js))
    rate = a.rate or float(meta.get('sample_rate', 2048000))
    t_ref = float(meta.get('first_sample_unix', os.path.getmtime(a.recording)))
    label = a.label or ('-'.join(p.replace(':', '_') for p in a.perturb) or 'plain') + (
        '-' + os.path.splitext(os.path.basename(a.opt))[0] if a.opt else '')
    work = os.path.join(ROOT, 'build', '_tmp', 'bench', '%s-%s' % (os.path.splitext(os.path.basename(a.recording))[0],
                                                                   label))
    os.makedirs(work, exist_ok=True)
    data = os.path.join(work, 'input.int8x2')
    if not a.reuse:
        secs, notes = prepare(a.recording, data, rate, a.q_negated, a.perturb)
        print('recording %.1f s at %.3f MS/s%s; options %s' % (secs, rate / 1e6, ''.join('; ' + n for n in notes),
                                                               a.opt or 'pocket_trk defaults'))
    rows = []
    for i in range(a.runs):
        d = os.path.join(work, 'run%d' % (i + 1))
        os.makedirs(d, exist_ok=True)
        if not a.reuse:
            wall, code = run_once(d, data, rate, a.opt)
            print('run %d: %.0f s, exit %d' % (i + 1, wall, code))
        if os.path.exists(os.path.join(d, 'trk.log')):
            rows.append(analyse(d, t_ref))
    print()
    print('run    ' + '  '.join(h for _, h, _ in COLS))
    for i, r in enumerate(rows):
        print('%-5d  %s' % (i + 1, fmt(r)))
    if len(rows) > 1:
        med = {}
        for k, _, _ in COLS:
            v = [r[k] for r in rows if r.get(k) is not None]
            med[k] = float(np.median(v)) if v else None
        print('median %s' % fmt(med))
    try:
        os.remove(data)         # (the converted copy: as big as the recording)
    except OSError:
        pass


if __name__ == '__main__':
    main()

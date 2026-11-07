"""Regression replay of this program's GNSS position solver (app/spp.py) on recorded sessions.

Every GNSS capture keeps pocket_trk's log and RTCM beside it (data/captures/<capture>.trk<N>.log / .rtcm). This
feeds them through spp.Solver exactly as gnss_rx does - all of the session's decoded orbits given at the start (a
warm start), then the navigation bits and measurements in order - and prints one line per session: how many epochs
gave a good / quick / no fix, and the good fixes' 95 % radius, 50 % radius and middle as distances from the home
setting in metres (never coordinates).

    python build/tests/spp_replay.py                     every capture with pocket_trk files
    python build/tests/spp_replay.py gnss-20261106-203728 [more captures]
    python build/tests/spp_replay.py --git HEAD          also replay git HEAD's spp.py (before / after a change)
    python build/tests/spp_replay.py --prior other       with the static height prior (the median of the OTHER
                                                         sessions' WAAS-corrected middles, never the session's own)
    --max-s 3600   stop each session after this many receiver seconds

Acceptance figures of week 2: R95 <= 5 m per session with the prior, session middles within ~1 m of each other.
Reads the data folder (IRIDIUM_TRACKER_DATA, else data/ beside app/); writes nothing there."""
import argparse
import collections
import glob
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
import tempfile

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, 'app'))
import gnss  # noqa: E402
import rtcm3  # noqa: E402
import spp  # noqa: E402

DATA = os.environ.get('IRIDIUM_TRACKER_DATA') or os.path.join(ROOT, 'data')
RECORD = re.compile(r'(?<!^)(?=[$](?:CH|OBS|NAV|EPH|POS|TIME|SAT|LOG),)')     # (as gnss_rx._RECORD)


def load_git(rev):
    src = subprocess.run(['git', '-C', ROOT, 'show', '%s:app/spp.py' % rev], capture_output=True, text=True,
                         encoding='utf-8')
    if src.returncode:
        raise SystemExit('git show %s:app/spp.py failed: %s' % (rev, src.stderr.strip()))
    path = os.path.join(tempfile.mkdtemp(prefix='spp_git_'), 'spp_git.py')
    with open(path, 'w', encoding='utf-8') as f:
        f.write(src.stdout)
    spec = importlib.util.spec_from_file_location('spp_git', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sessions(names):
    out = []
    for log in sorted(glob.glob(os.path.join(DATA, 'captures', 'gnss-*.trk*.log'))):
        rtcm = log[:-4] + '.rtcm'
        cap = os.path.basename(log).split('.trk')[0]
        if os.path.exists(rtcm) and (not names or cap in names):
            out.append((cap, log, rtcm))
    return out


def orbits(rtcm_path, t_ref):
    out = []
    payloads, _ = rtcm3.split(open(rtcm_path, 'rb').read())
    for pl in payloads:
        e = rtcm3.decode(pl, t_ref)
        if not e or e[2].get('glo'):
            continue
        sat, sig, fields = e
        e = gnss.parse_eph('$EPH,0,%s,%s' % (sat, ','.join([sig] + [str(fields[k]) for k in gnss.FIELDS])))
        if e:
            out.append(e)
    return out


def cap_time(cap):
    """the capture's start from its name (local time), for resolving the orbits' GPS week"""
    import datetime
    return datetime.datetime.strptime(cap.split('-', 1)[1], '%Y%m%d-%H%M%S').timestamp()


def prior_from_others(cap):
    mids = []
    for p in glob.glob(os.path.join(DATA, 'catalog', 'gnss-*.json')):
        if os.path.basename(p).startswith(cap):
            continue
        try:
            c = json.load(open(p))
        except (OSError, ValueError):
            continue
        w = (c.get('session') or {}).get('waas')
        if w and w.get('r95_m', 99) < 10 and w.get('span_s', 0) >= 300:
            mids.append((w['height_m'], w['lat'], w['lon']))
    if len(mids) < 2:
        return None
    return tuple(float(np.median([m[i] for m in mids])) for i in range(3))


def replay(mod, log, ephs, prior, max_s):
    s = mod.Solver(iono='sbas', sbas_orbit=True, smooth=True)
    if prior:
        s.set_prior(*prior)
    for e in ephs:
        s.eph(e)
    sols = []
    with open(log, errors='replace') as f:
        for line in f:
            line = line.rstrip()
            for rec in (RECORD.split(line) if line.count('$') > 1 else (line,)):
                p = rec.split(',', 2)
                if len(p) < 2:
                    continue
                try:
                    if float(p[1]) > max_s:
                        return sols
                except ValueError:
                    continue
                try:
                    if rec.startswith('$NAV'):
                        s.nav(rec)
                    elif rec.startswith('$OBS'):
                        r = s.obs(rec)
                        if r:
                            sols.append(r)
                except Exception:       # (a damaged record: skipped, as gnss_rx does)
                    continue
    return sols


def kind(x):
    """what gnss_rx would emit: a good fix, a quick one (good geometry, ranges agree) or none"""
    if x['ok']:
        return 'ok'
    need = 3 if x.get('prior') else 4
    if x['nsat'] >= need and x['hdop'] <= spp.HDOP_MAX and x.get('res_rms') is not None and x['res_rms'] < 15.0:
        return 'quick'
    return None


def summary(sols, home):
    k = collections.Counter(kind(x) for x in sols)
    row = {'epochs': len(sols), 'ok': k['ok'], 'quick': k['quick'], 'none': k[None]}
    ok = [x for x in sols if x['ok']]
    if len(ok) >= 30 and home:
        lat0, lon0, h0 = math.radians(home[0]), math.radians(home[1]), home[2] * 1000.0
        x0 = spp.lla_to_ecef(lat0, lon0, h0)
        rot = np.array([[-math.sin(lon0), math.cos(lon0), 0],
                        [-math.sin(lat0) * math.cos(lon0), -math.sin(lat0) * math.sin(lon0), math.cos(lat0)],
                        [math.cos(lat0) * math.cos(lon0), math.cos(lat0) * math.sin(lon0), math.sin(lat0)]])
        e = np.array([rot @ (spp.lla_to_ecef(math.radians(x['lat']), math.radians(x['lon']), x['h']) - x0)
                      for x in ok])
        mid = np.median(e, axis=0)
        rad = np.hypot(e[:, 0] - mid[0], e[:, 1] - mid[1])
        row.update(r95=float(np.percentile(rad, 95)), cep50=float(np.median(rad)),
                   mid_m=float(math.hypot(mid[0], mid[1])), dh=float(mid[2]), mid_en=(float(mid[0]), float(mid[1])))
    return row


def line(name, r):
    s = '%-28s epochs %6d  ok %6d  quick %5d  none %5d' % (name, r['epochs'], r['ok'], r['quick'], r['none'])
    if 'r95' in r:
        s += '  | R95 %5.1f  CEP50 %4.1f  middle %5.1f m from home  height %+6.1f m' % (
            r['r95'], r['cep50'], r['mid_m'], r['dh'])
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('captures', nargs='*')
    ap.add_argument('--git', default='', help='also replay this git revision\'s app/spp.py')
    ap.add_argument('--prior', choices=('none', 'other'), default='none')
    ap.add_argument('--max-s', type=float, default=1e12)
    a = ap.parse_args()
    try:
        home = json.load(open(os.path.join(DATA, 'settings.json')))['home']
    except (OSError, ValueError, KeyError):
        home = None
    mods = [('now', spp)] + ([(a.git, load_git(a.git))] if a.git else [])
    found = sessions(set(a.captures))
    if not found:
        raise SystemExit('no capture with pocket_trk files (<capture>.trkN.log + .rtcm) in %s' % DATA)
    mids = collections.defaultdict(list)
    for cap, log, rtcm in found:
        ephs = orbits(rtcm, cap_time(cap))
        prior = prior_from_others(cap) if a.prior == 'other' else None
        if a.prior == 'other' and not prior:
            print('%s: no prior (fewer than 2 other sessions with a WAAS middle)' % cap)
        for label, mod in mods:
            r = summary(replay(mod, log, ephs, prior, a.max_s), home)
            print(line('%s %s' % (os.path.basename(log)[:-4], label), r))
            if 'mid_en' in r:
                mids[label].append(r['mid_en'])
    for label, m in mids.items():
        if len(m) >= 2:
            m = np.array(m)
            c = np.median(m, axis=0)
            print('session middles (%s): %d sessions, rms %.1f m about their median' % (
                label, len(m), float(np.sqrt(np.mean(np.sum((m - c) ** 2, axis=1))))))


if __name__ == '__main__':
    main()

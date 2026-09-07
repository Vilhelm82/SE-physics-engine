"""Decoder surface L_d(h, p): logical failure of the rotated surface code (Codex's stack, stack noise 1e-3) for a synthetic
per-gate channel with herald h and unheralded depolarising Pauli rate p, on a log grid. Run once; it is a property of the code
and the stack, not of any gate. Every gate configuration is then a point on this surface.
Run from the repo root: ~/.cache/seated-root-qec-venv/bin/python experiments/2026-09-06/sweep_decoder_surface.py
"""
import sys, argparse, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
from rlq.receipt import Receipt
from rlq.sweep import grid, progress
from rlq.decoder import surface_circuit, run_case
from rlq import figures as F

ap = argparse.ArgumentParser(); ap.add_argument('--shots', type=int, default=200000); ap.add_argument('--shots-d7', type=int, default=50000)
ap.add_argument('--distances', type=int, nargs='+', default=[3, 5, 7]); ap.add_argument('--noise', type=float, default=1e-3); args = ap.parse_args()
H = [0.] + list(np.logspace(-5, np.log10(3e-2), 8)); P = [0.] + list(np.logspace(-6, -2, 9))
r = Receipt('sweep_decoder_surface', f'L_d(h, p) on a {len(H)}x{len(P)} log grid, d={args.distances}, stack noise {args.noise:g}', driver=__file__)
r.held_out('bicycle code', 'stack noise other than the one value', 'gate duration / idle cost (the stack does not price time)',
           'non-depolarising Pauli mixes (X, Y, Z taken equal)', 'single seeds per point')
codes = {d: [surface_circuit(d, b, args.noise) for b in ('X', 'Z')] for d in args.distances}
pts = grid(h=H, p=P); rows = []; t0 = time.time(); k = 0
for i, pt in enumerate(pts, 1):
    ch = dict(herald=float(pt['h']), pauli=[1 - pt['p'], pt['p']/3, pt['p']/3, pt['p']/3], total_time=0.)
    for d in args.distances:
        shots = args.shots_d7 if d >= 7 else args.shots; fails = 0
        for code in codes[d]:
            res = run_case(code, ch, shots, 20260906 + k*1009); k += 1; fails += res['flags_visible']['failures']
        rows.append(dict(h=pt['h'], p=pt['p'], d=d, failures=fails, shots=2*shots, rate=fails/(2*shots)))
    progress(i, len(pts), t0, every=10)
r.record('surface', rows); r.record('axes', dict(h=H, p=P, distances=args.distances))
base = {d: next(x['rate'] for x in rows if x['h'] == 0 and x['p'] == 0 and x['d'] == d) for d in args.distances}
r.record('baseline_rates', base)
for d in args.distances:
    r.check(f'd={d}: L increases from (0,0) to (3e-2, 1e-2)', -(next(x['rate'] for x in rows if x['d'] == d and x['h'] == H[-1] and x['p'] == P[-1]) - base[d]), 0.)
    r.fork(f'd={d} baseline (h=p=0)', f'{base[d]:.2e}', 'stack noise alone')
# figure: contours of L_d over (h, p) with FD's operating point and RH-3's r-line for reference
fig, ax = F.grid(1, len(args.distances), figsize=(4.5*len(args.distances), 4.2))
ax = np.atleast_1d(ax)
for a, d in zip(ax, args.distances):
    T = np.array([[next(x['rate'] for x in rows if x['d'] == d and x['h'] == h and x['p'] == p) for p in P] for h in H])
    hh = np.array(H) + 1e-6; pp = np.array(P) + 1e-7
    cs = a.contourf(pp, hh, np.log10(np.maximum(T, 1e-6)), levels=14, cmap='viridis'); a.set_xscale('log'); a.set_yscale('log')
    a.plot(4.45e-6 + 1e-7, 2.65e-3 + 1e-6, '*', color='white', ms=10); a.text(6e-6, 3.2e-3, 'FD', color='white', fontsize=8)
    a.plot(4.44e-6 + 1e-7, 1e-6, 'o', color='white', ms=6); a.text(6e-6, 1.2e-6, 'accumulator, ideal', color='white', fontsize=8)
    a.set(xlabel='unheralded Pauli rate p (+1e-7)', ylabel='herald h (+1e-6)', title=f'log10 logical failure, d={d}')
    fig.colorbar(cs, ax=a, shrink=.8)
F.finish(fig, ax, f'Decoder surface L_d(h, p): rotated surface code, stack noise {args.noise:g}'); r.figure(fig, 'decoder_surface')
sys.exit(r.close())

"""Instrument sweep: (h, p, T, decomposition) for every configuration of the FD instrument with and without the register, with no
decoding in the loop. Axes: n, gain, delta_d, background gamma, register (gamma_r, gamma_phi, monitored), capture e_c; 'fd' = capture off.
Run from the repo root: ~/.cache/seated-root-qec-venv/bin/python experiments/2026-09-06/sweep_instrument.py
Join with the decoder surface via sweep_readout.py.
"""
import sys, argparse, time, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
from rlq.receipt import Receipt
from rlq.sweep import grid, progress
from rlq.instruments import Model, channel, calibrate, N
from rlq.words import loop_endpoints_with_effects
from rlq.channels import native_channel, operating_point
from rlq import ROOT

ap = argparse.ArgumentParser(); ap.add_argument('--t-c', type=float, default=1.0); ap.add_argument('--quick', action='store_true'); args = ap.parse_args()
AX = dict(n=[1, 2, 3], gain=[1e-3, 3e-3, 1e-2, 3e-2], delta_d=[0., 1e-5, 1e-4, 1e-3], gamma=[0., 1e-9, 1e-7, 1e-6, 1e-5, 1e-4])
REG = [dict(name='fd', capture=False)]
for e_c in (0., 0.05):
    for gr in (0., 1e-3, 1e-2):
        for gphi in (0., 1e-3, 1e-2):
            for mon in (True, False):
                if gr == 0. and not mon: continue
                REG.append(dict(name=f"acc e_c={e_c:g} gr={gr:g} gphi={gphi:g} {'mon' if mon else 'unmon'}", capture=True, e_c=e_c, gamma_r=gr, gamma_phi=gphi, register_heralded=mon))
if args.quick: AX = dict(n=[1], gain=[1e-2], delta_d=[1e-4], gamma=[1e-6]); REG = REG[:3]
r = Receipt('sweep_instrument', f"{np.prod([len(v) for v in AX.values()])} physics points x {len(REG)} register/capture configs, t_c={args.t_c}", driver=__file__,
            modules=[ROOT/'rlq/reflection_loop_finite_dump.py', ROOT/'rlq/qec_distance_stack.py'])
r.held_out('kappa/a and kappa t_d fixed at FD values (10, 40)', 'delta_r = 0', 'register modes with their own detuning', 'capture pulse bandwidth',
           'occupancy readout and 10->2 recombination as physical operations', 'the trine word (five-loop only)')
row0 = operating_point(background=1e-6); fd0 = native_channel(row0)
rows = []; t0 = time.time(); phys = grid(**{k: AX[k] for k in ('n', 'gain', 'delta_d', 'gamma')}); done = 0
for n in AX['n']:
    cases = [(g, dd, 0.) for g in AX['gain'] for dd in AX['delta_d']]
    batch_cases, batch_rates = [], []
    for c in cases:
        for gm in AX['gamma']: batch_cases.append(c); batch_rates.append(gm)
    loops, effects = loop_endpoints_with_effects(batch_cases, batch_rates, n=n)
    for i, (c, gm) in enumerate(zip(batch_cases, batch_rates)):
        U, cal = calibrate(args.t_c, case=c, gamma=gm, n=n, loops=loops[i], effects=effects[i])
        for cfg in REG:
            kw = {k: v for k, v in cfg.items() if k != 'name'}
            m = Model(t_c=args.t_c, case=c, gamma=gm, n=n, loops=loops[i], effects=effects[i], **kw)
            ch = channel(m, U if cfg['capture'] else np.zeros((2*N, 2), dtype=complex))
            rows.append(dict(n=n, gain=c[0], delta_d=c[1], gamma=gm, config=cfg['name'], capture=cfg['capture'], e_c=cfg.get('e_c', 0.), gamma_r=cfg.get('gamma_r', 0.),
                             gamma_phi=cfg.get('gamma_phi', 0.), monitored=cfg.get('register_heralded', True), captured=cal['captured_fraction'], DB_stack=cal['DB_stack'],
                             h=ch['herald'], p=ch['unheralded_pauli_rate'], q_raw=ch['unheralded_residual'], q_repaired=ch['repaired_residual'], T=ch['total_time'],
                             unobserved=ch['unobserved_absorption'], surviving=ch['surviving_error'], recovered=ch['recovered_error']))
        done += 1; progress(done, len(phys), t0, every=8)
        if done % 24 == 0: r.record('rows', rows); (r.dir/'partial.json').write_text(json.dumps(dict(rows=rows), default=float))
r.record('axes', AX); r.record('configs', [c['name'] for c in REG]); r.record('rows', rows)
ref = next(x for x in rows if x['n'] == 1 and x['gain'] == 1e-2 and x['delta_d'] == 1e-4 and x['gamma'] == 1e-6 and x['config'] == 'fd')
r.check("fd config at the operating point reproduces native_channel herald", ref['h'] - fd0['herald'], 1e-9)
r.check("fd config at the operating point reproduces native_channel unheralded rate", ref['p'] - sum(fd0['pauli'][1:])*(1 - fd0['herald']), 1e-11)
r.check("fd config at the operating point reproduces FD receipt q_u", ref['q_raw'] - row0['unheralded_residual'], 1e-12)
ideal = [x for x in rows if x['capture'] and x['e_c'] == 0 and x['gamma_r'] == 0 and x['gamma_phi'] == 0]
r.check("ideal capture: herald below 1e-4 x captured fraction at every physics point", max(x['h']/max(x['captured'], 1e-30) for x in ideal), 1e-4)
r.fork('ideal capture unheralded vs FD background across the map', 'reported', f"max ratio p_acc/p_fd = {max(x['p']/max(next(y['p'] for y in rows if y['config']=='fd' and (y['n'],y['gain'],y['delta_d'],y['gamma'])==(x['n'],x['gain'],x['delta_d'],x['gamma'])),1e-30) for x in ideal):.3f}")
sys.exit(r.close())

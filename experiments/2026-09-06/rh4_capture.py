"""RH-4 on the harness: the capture model at the absorbing interface, as a one-page driver.
Acceptance test for rlq phase 1: every grid row must reproduce docs/rh4_capture-checks.json (the pre-harness run) to 1e-15.
Run from the repository root:  ~/.cache/seated-root-qec-venv/bin/python experiments/2026-09-06/rh4_capture.py [--qec]
"""
import sys, json, argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
from rlq.receipt import Receipt
from rlq.instruments import Model, channel, calibrate, N, GAMMA
from rlq.channels import native_channel, operating_point
from rlq import figures as F, ROOT

ap = argparse.ArgumentParser(); ap.add_argument('--t-c', type=float, default=1.0); ap.add_argument('--qec', action='store_true'); ap.add_argument('--shots', type=int, default=200000)
args = ap.parse_args()
r = Receipt('rh4_capture', 'capture model at the absorbing interface (harness port; acceptance against the pre-harness receipt)', driver=__file__,
            modules=[ROOT/'reflection_loop_finite_dump.py', ROOT/'qec_distance_stack.py'])
r.held_out('occupancy readout and 10->2 recombination as physical operations (folded into gamma_phi)', 'register modes with their own detuning',
           'finite-bandwidth capture pulses', 'd >= 7 and the bicycle code at this interface', 'single seeds', 'hardware')

row = operating_point(background=GAMMA); fd = native_channel(row)
off = channel(Model(capture=False), np.zeros((2*N, 2), dtype=complex))
r.check('capture off: herald == FD', off['herald'] - fd['herald'], 1e-9)
r.check('capture off: conditional Pauli vector == FD', float(np.abs(np.array(off['pauli']) - np.array(fd['pauli'])).max()), 1e-9)
U, cal = calibrate(args.t_c); r.record('calibration', cal)
r.check('capture takes what FD absorbed', cal['captured_fraction'] - fd['herald'], 2e-6)

base_unh = sum(fd['pauli'][1:])*(1 - fd['herald']); r_eff = lambda unh: (unh - base_unh)/(0.75*fd['herald'])
grid = []
for e_c in (0., 0.01, 0.05):
    for gr in (0., 1e-4, 1e-3, 1e-2):
        for gphi in (0., 1e-4, 1e-3, 1e-2):
            for rh in (True, False):
                if gr == 0. and not rh: continue
                if (e_c and (gr or gphi)) and not (gr in (0., 1e-3) and gphi in (0., 1e-3)): continue
                ch = channel(Model(e_c=e_c, t_c=args.t_c, gamma_r=gr, gamma_phi=gphi, register_heralded=rh), U)
                ch.update(e_c=e_c, gamma_r=gr, gamma_phi=gphi, register_heralded=rh, r_eff=r_eff(ch['unheralded_pauli_rate'])); grid.append(ch)
r.record('grid', grid)
# acceptance: bit-for-bit against the pre-harness receipt
old = json.load(open(ROOT/'docs/rh4_capture-checks.json'))['grid']
key = lambda g: (g['e_c'], g['gamma_r'], g['gamma_phi'], g['register_heralded'])
olds = {key(g): g for g in old}; worst = 0.
for g in grid:
    o = olds[key(g)]; worst = max(worst, abs(g['herald'] - o['herald']), abs(g['unheralded_pauli_rate'] - o['unheralded_pauli_rate']))
r.check(f'acceptance: {len(grid)} grid rows reproduce the pre-harness receipt (herald, unheralded)', worst, 1e-15)
ideal = next(g for g in grid if not (g['e_c'] or g['gamma_r'] or g['gamma_phi']))
r.fork('ideal register', 'a', f"herald {ideal['herald']:.1e}, unheralded {ideal['unheralded_pauli_rate']:.3e} = FD background")
for g in grid:
    if g['gamma_r'] == 1e-2 and g['gamma_phi'] == 0 and g['e_c'] == 0:
        r.fork(f"loss 1e-2, {'monitored' if g['register_heralded'] else 'unmonitored'}", 'ties FD' if g['register_heralded'] else 'loses to FD', f"r_eff {g['r_eff']:.3f}")
# figure: r_eff and herald against register rates
fig, ax = F.grid(1, 2, figsize=(12.5, 4.5))
for sel, color, ls, name in [((lambda g: g['gamma_phi'] == 0 and not g['register_heralded']), F.BR, '-', 'loss, unmonitored'),
                             ((lambda g: g['gamma_phi'] == 0 and g['register_heralded'] and g['gamma_r'] > 0), F.BR, '--', 'loss, monitored'),
                             ((lambda g: g['gamma_r'] == 0 and g['gamma_phi'] > 0), F.TE, '-', 'dephasing')]:
    pts = sorted([(max(g['gamma_r'], g['gamma_phi']), g) for g in grid if g['e_c'] == 0 and sel(g)], key=lambda t: t[0])
    ax[0].semilogx([p[0] for p in pts], [p[1]['r_eff'] for p in pts], 'o' + ls, color=color, label=name)
    ax[1].loglog([p[0] for p in pts], [max(p[1]['herald'], 1e-12) for p in pts], 'o' + ls, color=color, label=name)
ax[0].axhspan(0.1, 0.3, color=F.GY, alpha=.08); ax[0].set(xlabel='register rate', ylabel=r'$r_{\rm eff}$', title='RH-4: register requirement')
ax[1].axhline(fd['herald'], color=F.GY, ls='--', lw=.8); ax[1].set(xlabel='register rate', ylabel='herald per gate', title='what goes back to FD')
F.finish(fig, ax, 'RH-4 on the harness'); r.figure(fig, 'register_requirement')
if args.qec:
    from rlq.decoder import surface_circuit, run_case
    codes = [surface_circuit(d, b, 1e-3) for d in (3, 5) for b in ('X', 'Z')]; runs = []
    for name, chn in [('fd', fd), ('acc_ideal', ideal)]:
        for code in codes:
            res = run_case(code, chn, args.shots, 20260906 + len(runs)*1009); runs.append(dict(scenario=name, code=code.name, basis=code.basis, failures=res['flags_visible']['failures']))
    r.record('qec', runs)
sys.exit(r.close())

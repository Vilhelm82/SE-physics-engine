"""Readout: join the instrument sweep with the decoder surface and read the optimum off the chart.
For every instrument row, L_d = surface(d, h, p). Then: best config per physics point, the FD-vs-accumulator margin as a function of
every axis, and the register requirement as a heatmap. Run from the repo root after both sweeps.
"""
import sys, argparse, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import numpy as np
from rlq.receipt import Receipt
from rlq.sweep import DecoderSurface
from rlq import figures as F, ROOT

ap = argparse.ArgumentParser(); ap.add_argument('--surface', default='runs/2026-09-06/sweep_decoder_surface/receipt.json')
ap.add_argument('--instrument', default='runs/2026-09-06/sweep_instrument/receipt.json'); ap.add_argument('--d', type=int, default=5); args = ap.parse_args()
S = DecoderSurface.from_receipt(ROOT/args.surface); inst = json.load(open(ROOT/args.instrument)); rows = inst['records']['rows']; AX = inst['records']['axes']
d = args.d
r = Receipt('sweep_readout', f'instrument sweep x decoder surface at d={d}: optimum per physics point, margins, register requirement', driver=__file__)
r.held_out('surface interpolated in log-log between measured (h, p) nodes; single seeds per node', 'gate duration not priced (the stack does not idle-cost)',
           'd = 7 surface at 50k shots only', 'everything the two sweeps held out')
for x in rows:
    for dd in S.d: x[f'L{dd}'] = float(S(dd, x['h'], x['p'])[0])
key = lambda x: (x['n'], x['gain'], x['delta_d'], x['gamma'])
phys = sorted({key(x) for x in rows}); by = {}
for x in rows: by.setdefault(key(x), []).append(x)
best = []
for k in phys:
    fd = next(x for x in by[k] if x['config'] == 'fd'); acc = [x for x in by[k] if x['capture']]
    ideal = next(x for x in acc if x['e_c'] == 0 and x['gamma_r'] == 0 and x['gamma_phi'] == 0)
    b = min(acc, key=lambda x: x[f'L{d}'])
    best.append(dict(n=k[0], gain=k[1], delta_d=k[2], gamma=k[3], fd=fd[f'L{d}'], ideal=ideal[f'L{d}'], best_acc=b[f'L{d}'], best_config=b['config'],
                     winner='accumulator' if b[f'L{d}'] < fd[f'L{d}'] else 'fd', margin=fd[f'L{d}']/max(b[f'L{d}'], 1e-12),
                     worst_acc_that_still_beats_fd=max([x[f'L{d}'] for x in acc if x[f'L{d}'] < fd[f'L{d}']] or [np.nan])))
r.record('best_per_physics_point', best); r.record('rows_with_L', rows)
wins = sum(b['winner'] == 'accumulator' for b in best)
r.fork(f'd={d}: physics points where SOME accumulator config beats FD', f'{wins}/{len(best)}', f"ideal register beats FD at {sum(b['ideal'] < b['fd'] for b in best)}/{len(best)}")
op = next(b for b in best if (b['n'], b['gain'], b['delta_d'], b['gamma']) == (1, 1e-2, 1e-4, 1e-6))
r.fork('operating point (n=1, 1%, 1e-4, 1e-6)', op['winner'], f"FD L={op['fd']:.2e}, ideal L={op['ideal']:.2e}, best acc L={op['best_acc']:.2e} ({op['best_config']})")
r.check('surface interpolation reproduces the RH-3 FD count at d=3 within 20%', S(3, 2.6505e-3, 4.451e-6)[0]/(445/400000) - 1 if 3 in S.d else 0., 0.2)

# ---- figures
Ld = lambda x: x[f'L{d}']
fig, ax = F.grid(3, 2)
# (0,0) the surface at d with every config at the operating physics point on it
opk = (1, 1e-2, 1e-4, 1e-6); pts = by[opk]; hs = np.array(S.hs) + S.h_floor; ps = np.array(S.ps) + S.p_floor
T = np.array([[float(S(d, h, p)[0]) for p in S.ps] for h in S.hs])
cs = ax[0, 0].contourf(ps, hs, np.log10(T), levels=14, cmap='viridis'); ax[0, 0].set_xscale('log'); ax[0, 0].set_yscale('log')
for x in pts:
    c = F.GY if x['config'] == 'fd' else (F.TE if (x['gamma_r'] == 0 and x['gamma_phi'] == 0) else F.AC if not x['monitored'] else F.BL)
    ax[0, 0].plot(x['p'] + S.p_floor, x['h'] + S.h_floor, '*' if x['config'] == 'fd' else 'o', color=c, ms=10 if x['config'] == 'fd' else 5)
ax[0, 0].set(xlabel='unheralded p (+floor)', ylabel='herald h (+floor)', title=f'Every config at (n=1, 1%, 1e-4, 1e-6) on the d={d} surface\nstar FD; teal ideal register; blue monitored; red unmonitored')
fig.colorbar(cs, ax=ax[0, 0], shrink=.8)
# (0,1) L vs gamma at n=1, 1%, 1e-4: FD, ideal, best monitored 1e-3, unmonitored 1e-2
for name, sel, color, ls in [('FD', lambda x: x['config'] == 'fd', F.GY, '-'), ('ideal register', lambda x: x['capture'] and x['gamma_r'] == 0 and x['gamma_phi'] == 0 and x['e_c'] == 0, F.TE, '-'),
                             ('reg 1e-3/1e-3 monitored', lambda x: x['capture'] and x['gamma_r'] == 1e-3 and x['gamma_phi'] == 1e-3 and x['monitored'] and x['e_c'] == 0, F.BL, '--'),
                             ('reg 1e-2 loss unmonitored', lambda x: x['capture'] and x['gamma_r'] == 1e-2 and x['gamma_phi'] == 0 and not x['monitored'] and x['e_c'] == 0, F.AC, ':')]:
    sub = sorted([x for x in rows if x['n'] == 1 and x['gain'] == 1e-2 and x['delta_d'] == 1e-4 and sel(x)], key=lambda x: x['gamma'])
    ax[0, 1].loglog([x['gamma'] + 1e-10 for x in sub], [Ld(x) for x in sub], 'o' + ls, color=color, label=name)
ax[0, 1].set(xlabel=r'background $\gamma/a$ (+1e-10)', ylabel=f'logical failure L{d}', title='n=1, 1% gain, delta_d=1e-4')
# (1,0) L vs gain at n=1, delta 1e-4, gamma 1e-6
for name, sel, color in [('FD', lambda x: x['config'] == 'fd', F.GY), ('ideal register', lambda x: x['capture'] and x['gamma_r'] == 0 and x['gamma_phi'] == 0 and x['e_c'] == 0, F.TE),
                         ('reg 1e-3/1e-3 monitored', lambda x: x['capture'] and x['gamma_r'] == 1e-3 and x['gamma_phi'] == 1e-3 and x['monitored'] and x['e_c'] == 0, F.BL)]:
    sub = sorted([x for x in rows if x['n'] == 1 and x['delta_d'] == 1e-4 and x['gamma'] == 1e-6 and sel(x)], key=lambda x: x['gain'])
    ax[1, 0].loglog([x['gain'] for x in sub], [Ld(x) for x in sub], 'o-', color=color, label=name)
ax[1, 0].set(xlabel=r'gain error $\epsilon$', ylabel=f'L{d}', title='n=1, delta_d=1e-4, gamma=1e-6')
# (1,1) register requirement heatmap: L(fd)/L(acc) over (gamma_r, gamma_phi), monitored, at the operating point
grs = [0., 1e-3, 1e-2]; M = np.zeros((3, 3))
fdL = next(x for x in pts if x['config'] == 'fd')[f'L{d}']
for i, gr in enumerate(grs):
    for j, gp in enumerate(grs):
        x = next(y for y in pts if y['capture'] and y['gamma_r'] == gr and y['gamma_phi'] == gp and y['monitored'] and y['e_c'] == 0); M[i, j] = fdL/Ld(x)
im = ax[1, 1].imshow(np.log10(M), cmap='RdBu', vmin=-1, vmax=1, origin='lower'); ax[1, 1].set_xticks(range(3)); ax[1, 1].set_yticks(range(3))
ax[1, 1].set_xticklabels([f'{g:g}' for g in grs]); ax[1, 1].set_yticklabels([f'{g:g}' for g in grs])
for i in range(3):
    for j in range(3): ax[1, 1].text(j, i, f'{M[i, j]:.2f}x', ha='center', va='center', fontsize=9)
ax[1, 1].set(xlabel=r'register dephasing $\gamma_\phi$', ylabel=r'register loss $\gamma_r$ (monitored)', title=f'L{d}(FD) / L{d}(accumulator) at the operating point; >1 = accumulator wins'); fig.colorbar(im, ax=ax[1, 1], shrink=.8)
# (2,0) L vs n
for name, sel, color in [('FD', lambda x: x['config'] == 'fd', F.GY), ('ideal register', lambda x: x['capture'] and x['gamma_r'] == 0 and x['gamma_phi'] == 0 and x['e_c'] == 0, F.TE),
                         ('reg 1e-3/1e-3 monitored', lambda x: x['capture'] and x['gamma_r'] == 1e-3 and x['gamma_phi'] == 1e-3 and x['monitored'] and x['e_c'] == 0, F.BL)]:
    sub = sorted([x for x in rows if x['gain'] == 1e-2 and x['delta_d'] == 1e-4 and x['gamma'] == 1e-6 and sel(x)], key=lambda x: x['n'])
    ax[2, 0].semilogy([x['n'] for x in sub], [Ld(x) for x in sub], 'o-', color=color, label=name)
ax[2, 0].set(xlabel='return index n', ylabel=f'L{d}', title='1% gain, delta_d=1e-4, gamma=1e-6', xticks=[1, 2, 3])
# (2,1) who wins over (gain, gamma) at n=1, delta 1e-4: margin of the best monitored-register config over FD
G, GM = AX['gain'], AX['gamma']; W = np.zeros((len(GM), len(G)))
for i, gm in enumerate(GM):
    for j, g in enumerate(G):
        b = next(x for x in best if (x['n'], x['gain'], x['delta_d'], x['gamma']) == (1, g, 1e-4, gm)); W[i, j] = np.log10(b['fd']/max(b['best_acc'], 1e-12))
im2 = ax[2, 1].imshow(W, cmap='RdBu', vmin=-1, vmax=1, origin='lower', aspect='auto'); ax[2, 1].set_xticks(range(len(G))); ax[2, 1].set_xticklabels([f'{g:g}' for g in G])
ax[2, 1].set_yticks(range(len(GM))); ax[2, 1].set_yticklabels([f'{g:g}' for g in GM])
ax[2, 1].set(xlabel=r'gain $\epsilon$', ylabel=r'background $\gamma$', title=f'log10 L{d}(FD)/L{d}(best accumulator config), n=1, delta_d=1e-4'); fig.colorbar(im2, ax=ax[2, 1], shrink=.8)
F.finish(fig, ax, f'Reading the optimum off the chart: instrument sweep x decoder surface, d={d}'); r.figure(fig, 'readout')
# best-config table
lines = ['| n | gain | delta_d | gamma | L(FD) | L(ideal) | L(best acc) | best config | winner | margin |', '|---|---|---|---|---|---|---|---|---|---|']
for b in best: lines.append(f"| {b['n']} | {b['gain']:g} | {b['delta_d']:g} | {b['gamma']:g} | {b['fd']:.2e} | {b['ideal']:.2e} | {b['best_acc']:.2e} | {b['best_config']} | {b['winner']} | {b['margin']:.2f}x |")
(r.dir/'best_config_table.md').write_text('\n'.join(lines) + '\n'); print(f"  table {r.dir/'best_config_table.md'}")
sys.exit(r.close())

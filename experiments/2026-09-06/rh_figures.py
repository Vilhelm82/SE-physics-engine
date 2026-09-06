"""Figure for RH-2/2P/2A/3/4 from their receipts (no recomputation). Style follows reflection_loop_finite_dump.plot."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))  # repo root on sys.path (reorg 2026-09-06)
import json, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]/'docs'/'receipts'; FIG = ROOT.parent/'figures'
BR, TE, BL, GY = '#976b3b', '#147d72', '#546da8', '#252525'
x_of = lambda e: 15*np.pi/32*e

r2 = json.load(open(ROOT/'rh2-checks.json')); r2p = json.load(open(ROOT/'rh2_precision-checks.json'))
r5 = json.load(open(ROOT/'rh2p_five-checks.json')); r3 = json.load(open(ROOT/'rh3_qec-checks.json')); r4 = json.load(open(ROOT/'rh4_capture-checks.json'))

fig, ax = plt.subplots(3, 2, figsize=(12.5, 12), layout='constrained')
# (0,0) DB_stack vs eps: 50-digit points, float32 points, closed forms
ee = np.array(r2p['eps']); ef = np.logspace(np.log10(ee.min()), np.log10(ee.max()), 100)
for rec, color, name, a0 in [(r2p, BR, 'trine', 3375*np.sqrt(5)*np.pi**3/32768), (r5, TE, 'five', 2025*np.sqrt(5)*np.pi**3/32768)]:
    mp50 = [float(rec['table'][str(e)]['mp50']['k2']) for e in ee]; f32 = [rec['table'][str(e)]['f32']['k2'] for e in ee]
    ax[0, 0].loglog(ee, mp50, 'o', color=color, label=f'{name}: 50 digits')
    ax[0, 0].loglog(ee, f32, 'x', color=color, alpha=.45, label=f'{name}: float32')
    ax[0, 0].loglog(ef, a0*ef**3*(1 - 55*ef/32), '-', color=color, alpha=.8, label=f'{name}: $a_0\\epsilon^3(1-55\\epsilon/32)$')
ax[0, 0].set(xlabel=r'Gain error $\epsilon$', ylabel=r'$DB_{\rm stack}$', title=r'Tight-frame deviation: $\sqrt{5}\,x^3$ and $(3\sqrt{5}/5)\,x^3$, $x=15\pi\epsilon/32$')
# (0,1) Theorem 2 conditioning: float64 relative error by path (trine)
for k, color, name in [('k1', TE, 'svd'), ('k2', BL, 'stable 2x2'), ('k3', BR, 'tr/det'), ('k4', GY, r'$I-K^\dagger K$')]:
    rel = [r2p['theorem2'][str(e)]['paths'][k]['rel_err_64'] for e in ee]
    rel = [v if (v is not None and v > 0) else np.nan for v in rel]
    ax[0, 1].loglog(ee, rel, 'o-', color=color, label=name)
ax[0, 1].axhline(1, color=GY, lw=.6, ls=':')
ax[0, 1].set(xlabel=r'Gain error $\epsilon$', ylabel='float64 relative error of $DB_{\\rm stack}$', title='Cella U-0784: same $a$, path-dependent conditioning $b_k$ (trine)')
# (1,0) accumulator floor vs eps (RH-2 rows) with closed forms
for word, color in [('trine', BR), ('five', TE)]:
    rows = [r for r in r2['words'][word]['rows'] if r['case'][1] == 0 and r['case'][2] == 0]
    e = np.array([r['case'][0] for r in rows])
    ax[1, 0].loglog(e, [r['discard_unconditional'] for r in rows], 's--', color=color, alpha=.6, label=f'{word}: discard, unconditional')
    ax[1, 0].loglog(e, [r['per_dump']['infidelity'] for r in rows], 'v--', color=color, alpha=.6, label=f'{word}: per-dump fold-back (RH-1b)')
    ax[1, 0].loglog(e, [r['acc_self']['infidelity'] for r in rows], 'o-', color=color, label=f'{word}: accumulator, unconditional')
    ax[1, 0].loglog(e, [r['conditional'] for r in rows], '.', color=color, alpha=.5, label=f'{word}: no-click conditional')
    coef = 25/4 if word == 'trine' else 15/4
    ax[1, 0].loglog(ef, coef*x_of(ef)**8, ':', color=color, label=f'{word}: $({int(coef*4)}/4)\\,x^8$')
ax[1, 0].set(xlabel=r'Gain error $\epsilon$', ylabel='Average gate infidelity', title='RH-2: the return bus is eighth order (ideal interface)')
# (1,1) RH-3 decoder
runs = r3['runs']
def fails(scn, d): return sum(r['flags_visible']['failures'] for r in runs if r['scenario'] == scn and r['distance'] == d)
rs = [0, 0.01, 0.1, 0.3, 1.0]
for d, color in [(3, TE), (5, BL), (7, BR)]:
    ax[1, 1].semilogx([max(r, 3e-3) for r in rs], [fails(f'acc_r{r:g}', d) for r in rs], 'o-', color=color, label=f'd={d}: accumulator(r)')
    ax[1, 1].axhline(fails('fd', d), color=color, ls='--', alpha=.7); ax[1, 1].axhline(fails('baseline', d), color=color, ls=':', alpha=.7)
ax[1, 1].set_yscale('log'); ax[1, 1].set(xlabel='Register replacement error $r$ (0 plotted at $3\\times10^{-3}$)', ylabel='Logical failures per 400k (X+Z)',
                                          title='RH-3: decoder at stack noise $10^{-3}$ (dashed FD, dotted perfect gate)')
ax[1, 1].axvspan(0.1, 0.3, color=GY, alpha=.08)
# (2,0)/(2,1) RH-4 sweep: r_eff and herald vs rate
g = r4['grid']; fd_h = r4['fd_channel']['herald']
def curve(sel, key):
    pts = sorted([(x, row[key]) for row, x in sel], key=lambda t: t[0]); return [p[0] for p in pts], [p[1] for p in pts]
series = [([(row, row['gamma_r']) for row in g if row['e_c'] == 0 and row['gamma_phi'] == 0 and not row['register_heralded'] and row['gamma_r'] > 0], BR, '-', 'loss, unmonitored'),
          ([(row, row['gamma_r']) for row in g if row['e_c'] == 0 and row['gamma_phi'] == 0 and row['register_heralded'] and row['gamma_r'] > 0], BR, '--', 'loss, monitored'),
          ([(row, row['gamma_phi']) for row in g if row['e_c'] == 0 and row['gamma_r'] == 0 and row['gamma_phi'] > 0], TE, '-', 'dephasing'),
          ([(row, row['gamma_r']) for row in g if row['e_c'] == 0 and row['gamma_r'] == row['gamma_phi'] and row['gamma_r'] > 0 and not row['register_heralded']], BL, '-', 'both, unmonitored'),
          ([(row, row['gamma_r']) for row in g if row['e_c'] == 0 and row['gamma_r'] == row['gamma_phi'] and row['gamma_r'] > 0 and row['register_heralded']], BL, '--', 'both, monitored')]
for sel, color, ls, name in series:
    xs, ys = curve(sel, 'r_eff'); ax[2, 0].semilogx(xs, ys, 'o' + ls, color=color, label=name)
    xs, ys = curve(sel, 'herald'); ax[2, 1].loglog(xs, ys, 'o' + ls, color=color, label=name)
ax[2, 0].axhspan(0.1, 0.3, color=GY, alpha=.08); ax[2, 0].text(1.2e-4, 0.2, r'$r^*$ from RH-3', fontsize=8, color=GY)
ax[2, 0].set(xlabel=r'Register rate $\gamma_r$ or $\gamma_\phi$ (per $1/a$; word $\approx 130/a$)', ylabel=r'$r_{\rm eff}$ (captured shots effectively replaced)',
             title='RH-4: register requirement at the absorbing interface')
ax[2, 1].axhline(fd_h, color=GY, ls='--', lw=.8); ax[2, 1].text(1.2e-4, fd_h*0.55, 'FD herald 0.265%', fontsize=8, color=GY); ax[2, 1].set_ylim(top=fd_h*3)
ax[2, 1].set(xlabel=r'Register rate $\gamma_r$ or $\gamma_\phi$', ylabel='Herald per gate', title='RH-4: what the register hands back to FD')
for a in ax.flat: a.grid(alpha=.2); a.legend(frameon=False, fontsize=7)
ax[0, 0].legend(frameon=False, fontsize=7, loc='upper left'); ax[1, 0].legend(frameon=False, fontsize=6.5, loc='lower right', ncol=2)
fig.suptitle('Accumulator recovery, 2026-09-06: closed forms (RH-2A), precision separation (RH-2P), decoder (RH-3), register (RH-4)')
out = FIG/'rh-accumulator-tradeoffs.png'; fig.savefig(out, dpi=170); plt.close(fig); print('wrote', out)

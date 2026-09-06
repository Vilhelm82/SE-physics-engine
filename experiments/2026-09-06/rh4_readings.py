"""RH-4R -- FD's six readings (reflection_loop_finite_dump.plot) with the accumulator on the same axes.
FD curves come from FD's receipt; the accumulator is computed by rh4_capture at the same (n, gain, delta, gamma) points,
in two register variants: ideal (gamma_r = gamma_phi = 0) and realistic (gamma_r = gamma_phi = 1e-3, monitored).
Written 2026-09-06 by Claire."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))  # repo root on sys.path (reorg 2026-09-06)
import json, time, hashlib, subprocess
from pathlib import Path
import numpy as np
import rlq.reflection_loop_finite_dump as FD
import rlq.rh4_capture as C
ROOT = Path(__file__).resolve().parents[2]
REG = {'ideal': dict(gamma_r=0., gamma_phi=0., register_heralded=True), 'realistic': dict(gamma_r=1e-3, gamma_phi=1e-3, register_heralded=True)}
GAMMAS = [0., 1e-15, 1e-12, 1e-9, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
T_C = 1.0


def grid_point(n, case, gammas):
    loops, effects = FD.loop_endpoints_with_effects([case]*len(gammas), gammas, n=n)
    out = []
    for i, g in enumerate(gammas):
        U, cal = C.calibrate(T_C, case=case, gamma=g, n=n, loops=loops[i], effects=effects[i])
        rec = dict(n=n, gain=case[0], delta_d=case[1], background_rate=g, captured=cal['captured_fraction'])
        for name, kw in REG.items():
            m = C.Model(t_c=T_C, case=case, gamma=g, n=n, loops=loops[i], effects=effects[i], **kw)
            ch = C.channel(m, U); rec[name] = {k: ch[k] for k in ('herald', 'unheralded_residual', 'repaired_residual', 'unobserved_absorption',
                                                                  'surviving_error', 'recovered_error', 'total_time', 'unheralded_pauli_rate')}
            rec[name]['total_absorption'] = rec[name]['herald'] + rec[name]['unobserved_absorption']*(1 - rec[name]['herald'])
        out.append(rec)
    return out


def main():
    t0 = time.time(); fd = json.load(open(ROOT/'docs/receipts/reflection-loop-finite-dump-checks.json'))
    rows = []
    for n in (1, 2, 3):
        for case in ((1e-3, 1e-4, 0.), (1e-2, 1e-4, 0.), (1e-2, 0., 0.)):
            rows += grid_point(n, case, GAMMAS if case[1] > 0 else [0.]); print(f"  n={n} case={case} done ({time.time()-t0:.0f}s)", flush=True)
    # register sweep at n=1, 1% gain, gamma=1e-6 (raw residual and herald vs rate)
    loops, effects = FD.loop_endpoints_with_effects([(1e-2, 1e-4, 0.)], [1e-6], n=1); U, _ = C.calibrate(T_C, loops=loops[0], effects=effects[0])
    sweep = []
    for rate in (1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2):
        for label, kw in (('loss, unmonitored', dict(gamma_r=rate, register_heralded=False)), ('loss, monitored', dict(gamma_r=rate, register_heralded=True)),
                          ('dephasing', dict(gamma_phi=rate)), ('both, monitored', dict(gamma_r=rate, gamma_phi=rate, register_heralded=True))):
            ch = C.channel(C.Model(t_c=T_C, loops=loops[0], effects=effects[0], **kw), U)
            sweep.append(dict(rate=rate, variant=label, herald=ch['herald'], unheralded_residual=ch['unheralded_residual'], repaired_residual=ch['repaired_residual']))
    res = dict(scope='FD readings with the accumulator overlaid; FD from its receipt, accumulator from rh4_capture at the same points', t_c=T_C, register_variants=REG,
               gammas=GAMMAS, rows=rows, register_sweep=sweep,
               provenance=dict(date='2026-09-06', runtime_s=round(time.time()-t0, 1), head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                               sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ('experiments/2026-09-06/rh4_readings.py', 'rlq/rh4_capture.py', 'rlq/reflection_loop_finite_dump.py')}))
    (ROOT/'docs/receipts/rh4_readings-checks.json').write_text(json.dumps(res, indent=2) + '\n')
    plot(res, fd, ROOT/'docs/figures/rh4-readings-vs-fd.png'); print(f"done {time.time()-t0:.0f}s", flush=True)


def plot(res, fd, path):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    BR, TE, BL, GY, AC = '#976b3b', '#147d72', '#546da8', '#252525', '#b03a2e'
    fig, ax = plt.subplots(3, 2, figsize=(12.5, 12), layout='constrained')
    rows = res['rows']; ns = np.array([1, 2, 3])
    pick = lambda **kw: [r for r in rows if all(np.isclose(r[k], v) for k, v in kw.items())]
    # (0,0) endpoint flags at 1% gain, Delta = 0
    for kind, color in [('trine', BR), ('five', TE)]:
        e = [next(x for x in fd['cases'][f'{kind}_n{n}']['errors'] if x['gain'] == .01 and x['delta_d'] == 0)['erasure'] for n in ns]
        ax[0, 0].plot(ns, [100*v for v in e], 'o-', color=color, label=f'FD {kind}')
    acc = [pick(n=n, gain=1e-2, delta_d=0.)[0] for n in ns]
    ax[0, 0].plot(ns, [100*r['realistic']['herald'] for r in acc], 's-', color=AC, label='accumulator, realistic register')
    ax[0, 0].plot(ns, [100*r['ideal']['herald'] + 1e-14 for r in acc], 'v--', color=AC, alpha=.6, label='accumulator, ideal register ($\\approx 0$)')
    ax[0, 0].set(xlabel='Return index n', ylabel='Flag / herald probability (%)', title=r'Endpoint flags at 1% gain, $\Delta=0$', xticks=ns)
    # (0,1) conditional infidelity vs duration at eps=1e-3, delta_d=1e-4, zero background
    for kind, color in [('trine', BR), ('five', TE)]:
        pts = [(fd['cases'][f'{kind}_n{n}']['total_time'], next(x for x in fd['cases'][f'{kind}_n{n}']['errors'] if x['gain'] == .001 and x['delta_d'] == 1e-4)['conditional_infidelity']) for n in ns]
        ax[0, 1].semilogy(*zip(*pts), 'o-', color=color, label=f'FD {kind}')
    acc = [pick(n=n, gain=1e-3, delta_d=1e-4, background_rate=0.)[0] for n in ns]
    for name, ls in (('ideal', '-'), ('realistic', '--')):
        ax[0, 1].semilogy([r[name]['total_time'] for r in acc], [r[name]['unheralded_residual'] for r in acc], 's' + ls, color=AC, label=f'accumulator, {name}: accepted residual')
    ax[0, 1].set(xlabel=r'Total duration $aT$', ylabel='Conditional / accepted infidelity', title=r'$\epsilon=10^{-3},\ \Delta_d/a=10^{-4}$, zero background')
    # (1,0) register sweep (the accumulator's dump-action analogue)
    fdq = next(r for r in fd['operating_curves'] if r['n'] == 1 and r['gain'] == .01 and np.isclose(r['background_rate'], 1e-6))
    for label, color, ls in (('loss, unmonitored', BR, '-'), ('loss, monitored', BR, '--'), ('dephasing', TE, '-'), ('both, monitored', BL, '--')):
        pts = [s for s in res['register_sweep'] if s['variant'] == label]
        ax[1, 0].loglog([s['rate'] for s in pts], [s['unheralded_residual'] for s in pts], 'o' + ls, color=color, label=label)
    ax[1, 0].axhline(fdq['unheralded_residual'], color=GY, ls=':', lw=1); ax[1, 0].text(1.2e-5, fdq['unheralded_residual']*1.3, 'FD: $5.93\\times10^{-6}$', fontsize=8, color=GY)
    ax[1, 0].set(xlabel=r'Register loss/dephasing rate (per $1/a$); n=1, 1% gain, $\gamma/a=10^{-6}$', ylabel='Unheralded residual per accepted run',
                 title='Register sweep: the accumulator\'s dump-action tradeoff')
    # (1,1) total absorption vs gamma, 1% gain, delta_d=1e-4
    for n, color in zip(ns, [TE, BL, BR]):
        fr = [r for r in fd['operating_curves'] if r['n'] == n and r['gain'] == .01 and r['background_rate'] > 0]
        ax[1, 1].loglog([r['background_rate'] for r in fr], [r['total_absorption'] for r in fr], color=color, label=f'FD n={n}, total')
        ax[1, 1].loglog([r['background_rate'] for r in fr], [r['heralded_fraction'] for r in fr], '--', color=color, alpha=.55)
        ar = [r for r in pick(n=n, gain=1e-2, delta_d=1e-4) if r['background_rate'] > 0]
        ax[1, 1].loglog([r['background_rate'] for r in ar], [r['ideal']['total_absorption'] for r in ar], 's-', color=color, ms=3, alpha=.8, label=f'accumulator n={n}, ideal: total = unobserved')
        ax[1, 1].loglog([r['background_rate'] for r in ar], [r['realistic']['herald'] for r in ar], 's:', color=color, ms=3, alpha=.8)
    ax[1, 1].set(xlabel=r'Unmonitored background rate $\gamma/a$', ylabel='Total absorption probability', title=r'1% gain, $\Delta_d/a=10^{-4}$; dashed = FD heralded, dotted = realistic-register herald')
    # (2,0) paired operating points
    for n, color in zip(ns, [TE, BL, BR]):
        for eps, style in ((.001, '-'), (.01, '--')):
            fr = [r for r in fd['operating_curves'] if r['n'] == n and r['gain'] == eps]
            ax[2, 0].loglog([100*r['heralded_fraction'] for r in fr], [r['unheralded_residual'] for r in fr], style, color=color, label=rf'FD n={n}, $\epsilon={eps:g}$')
            ar = pick(n=n, gain=eps, delta_d=1e-4)
            ax[2, 0].loglog([100*r['realistic']['herald'] for r in ar], [r['realistic']['unheralded_residual'] for r in ar], style, color=color, alpha=.5, lw=2.5)
    p_fd = fdq; p_ac = pick(n=1, gain=1e-2, delta_d=1e-4, background_rate=1e-6)[0]['realistic']
    ax[2, 0].plot(100*p_fd['heralded_fraction'], p_fd['unheralded_residual'], '*', color=GY, ms=10)
    ax[2, 0].plot(100*p_ac['herald'], p_ac['unheralded_residual'], '*', color=AC, ms=10)
    ax[2, 0].annotate(f"FD n=1, $\\gamma/a=10^{{-6}}$\n({100*p_fd['heralded_fraction']:.3f}%, {p_fd['unheralded_residual']:.2e})", xy=(100*p_fd['heralded_fraction'], p_fd['unheralded_residual']),
                      xytext=(.03, 3e-8), arrowprops=dict(arrowstyle='->', color='#444444'), fontsize=8)
    ax[2, 0].annotate(f"accumulator, realistic register\n({100*p_ac['herald']:.3f}%, {p_ac['unheralded_residual']:.2e})\nideal register: herald $\\approx 0$", xy=(100*p_ac['herald'], p_ac['unheralded_residual']),
                      xytext=(.0025, 3e-10), arrowprops=dict(arrowstyle='->', color=AC), fontsize=8, color=AC)
    ax[2, 0].set(xlabel='Heralded fraction (%)', ylabel='Unheralded residual per accepted run', title=r'Paired operating points, $\Delta_d/a=10^{-4}$; thick faint = accumulator (realistic)')
    # (2,1) residual decomposition, n=1, eps=1e-3
    fr = [r for r in fd['operating_curves'] if r['n'] == 1 and r['gain'] == .001 and r['background_rate'] > 0]
    ax[2, 1].loglog([r['background_rate'] for r in fr], [r['unheralded_residual'] for r in fr], color=GY, alpha=.5, label='FD accepted residual')
    ar = [r for r in pick(n=1, gain=1e-3, delta_d=1e-4) if r['background_rate'] > 0]
    for key, label, color, ls in (('unheralded_residual', 'accumulator accepted residual', TE, '-'), ('unobserved_absorption', 'unobserved absorption', BR, '--'),
                                  ('surviving_error', 'surviving-state error', BL, ':'), ('recovered_error', 'recovered-branch error', AC, '-.')):
        ax[2, 1].loglog([r['background_rate'] for r in ar], [max(r['ideal'][key], 1e-18) for r in ar], ls, color=color, label=label)
    ax[2, 1].set(xlabel=r'Unmonitored background rate $\gamma/a$', ylabel='Conditional average infidelity', title=r'Residual decomposition: n=1, $\epsilon=10^{-3}$, ideal register')
    for a in ax.flat: a.grid(alpha=.2); a.legend(frameon=False, fontsize=7)
    ax[2, 0].legend(frameon=False, fontsize=7, loc='upper center', ncol=2)
    fig.suptitle(r'FD readings with the accumulator overlaid: $\kappa/a=10$, $\kappa t_d=40$, capture $t_c=1/a$; realistic register $\gamma_r=\gamma_\phi=10^{-3}$, monitored')
    fig.savefig(path, dpi=170); plt.close(fig)


if __name__ == '__main__':
    main()

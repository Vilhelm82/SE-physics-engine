"""RH-2 -- coherent accumulator recovery on the STACKED leak.  Spec: docs/results/2026-09-06/2026-09-06-RH-2-spec.md.
Builds on rh1_common.py (unmodified; Codex's loop matrices and metrics).
Interface (declared model): instantaneous coherent P/Q separation at every dump; the Q amplitude is SENT
to a fresh, frozen 2-mode register slot with no measurement; ONE 'register occupied' readout after the
last loop; ONE fixed isometry U_cal^dag (2N -> 2) back to P, then G.  Register amplitude outside
range(U_cal) is heralded leftover and replaced by I/2 (RH-1b's convention).  Nothing else is inserted.
Written 2026-09-06 by Claire after the spec and predictions were committed (ffa01f3).
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))  # repo root on sys.path (reorg 2026-09-06)
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import polar
from rlq.rh1_common import (word_loops, compose_discard, channel_metrics, replacement_kraus, foldback,
                        pattern_kraus, db, check, norm_json, G, I2)

EPS0 = 1e-2
CASES_GAIN = [(5e-3, 0., 0.), (1e-2, 0., 0.), (2e-2, 0., 0.)]
CASES_DET = [(0., 1e-4, 0.), (1e-2, 1e-4, 0.)]


def stacked(loops):
    k, leaks, out = compose_discard(loops)          # k = S_N..S_1, leaks = [F_1..F_N]
    F = np.vstack(leaks)                            # (2N) x 2 stacked leak
    ftf = F.conj().T@F
    u, h = polar(F)                                 # F = u h, u: 2N x 2 isometry, h >= 0 on P
    return dict(k=k, F=F, ftf=ftf, u=u, h=h, out=out, DB=db(F)['DB'],
                c=float(np.trace(ftf).real/2),
                identity_defect=float(np.linalg.norm(ftf + k.conj().T@k - I2)),
                isometry_defect=float(np.linalg.norm(u.conj().T@u - I2)))


def accumulator(st, u_cal):
    n = u_cal.shape[0]
    leftover = (np.eye(n) - u_cal@u_cal.conj().T)@st['F']
    kraus = [st['k'], G@u_cal.conj().T@st['F']] + replacement_kraus(leftover)
    m = channel_metrics(kraus)
    m['leftover'] = float(np.linalg.norm(leftover)**2/2)
    return m


def per_dump(loops, cal_loops):
    # RH-1b's per-dump fold-back, recomputed here as the comparison line (calibrated at cal_loops).
    branches = []
    for L, Lc in zip(loops, cal_loops):
        s, e = L[2:, 2:], L[:2, 2:]
        branches.append((s, foldback(Lc[2:, 2:], Lc[:2, 2:])@e))
    return channel_metrics(pattern_kraus(branches))


def run():
    t0 = time.time()
    res = dict(spec='docs/results/2026-09-06/2026-09-06-RH-2-spec.md', predictions_sealed='docs/prereg/RH-2/DESIGNER-PREDICTIONS.md',
               interface='coherent send into frozen register at every dump; one readout; one isometry back', words={}, checks=[])
    ch = res['checks']
    cases = CASES_GAIN + CASES_DET
    for kind in ('trine', 'five'):
        print(f"\n=== {kind} ===", flush=True)
        loops_all = word_loops(kind, cases)
        cal = stacked(loops_all[1])                                # U_cal at eps_0 = 1e-2, gain-only
        N = loops_all.shape[1]
        rows = []
        hdr = f"{'case':>22} {'DB_stack':>10} {'c':>10} {'cond':>10} {'discard':>10} {'per-dump':>10} {'acc self':>10} {'acc fixed':>10} {'floor':>10}"
        print(hdr, flush=True)
        for case, loops in zip(cases, loops_all):
            st = stacked(loops)
            floor = st['c']*st['DB']**2/6 if st['DB'] is not None else float('nan')
            row = dict(case=list(case), N=int(N), DB_stack=st['DB'], c=st['c'], singular_values=db(st['F'])['singular_values'],
                       identity_defect=st['identity_defect'], isometry_defect=st['isometry_defect'],
                       conditional=st['out']['conditional'], discard_unconditional=st['out']['unconditional'],
                       erasure=st['out']['erasure'], per_dump=per_dump(loops, loops_all[1]),
                       acc_self=accumulator(st, st['u']), acc_fixed=accumulator(st, cal['u']), floor=floor)
            rows.append(row)
            print(f"{str(case):>22} {row['DB_stack']:10.3e} {row['c']:10.3e} {row['conditional']:10.3e} {row['discard_unconditional']:10.3e} "
                  f"{row['per_dump']['infidelity']:10.3e} {row['acc_self']['infidelity']:10.3e} {row['acc_fixed']['infidelity']:10.3e} {floor:10.3e}", flush=True)

        # --- checks (numerical identities) and forks (outcomes, recorded as text) ---
        for row in rows:
            tag = f"{kind}{tuple(row['case'])}"
            check(ch, f"{tag} F^dagF + K^dagK = I", row['identity_defect'], 1e-11)
            check(ch, f"{tag} U isometry", row['isometry_defect'], 1e-11)
            for name in ('acc_self', 'acc_fixed', 'per_dump'):
                check(ch, f"{tag} {name} trace preserving", row[name]['tp_defect'], 1e-11)
                check(ch, f"{tag} {name} direct vs stable infidelity agree", row[name]['infidelity']-row[name]['direct_infidelity'], 1e-12)
            check(ch, f"{tag} self-calibrated leftover is zero", row['acc_self']['leftover'], 1e-24)
        g = {tuple(r['case']): r for r in rows}
        r0 = g[(EPS0, 0., 0.)]
        ratio_pd = r0['per_dump']['infidelity']/r0['discard_unconditional']
        check(ch, f"{kind} per-dump/discard ratio reproduces RH-1b (0.6668)", ratio_pd-0.6668, 2e-3)
        forks = {}
        for r in rows:
            fl = max(r['conditional'], r['floor'])
            ra, rb = r['acc_self']['infidelity']/fl, r['acc_self']['infidelity']/r['discard_unconditional']
            forks[str(tuple(r['case']))] = dict(self_over_floor=ra, self_over_discard=rb,
                                                fork='a' if ra <= 10 else ('b' if rb >= 0.1 else 'between'))
        if kind == 'trine':
            check(ch, "KILL 1: trine gain-only DB_stack <= 1e-2", max(g[c]['DB_stack'] for c in CASES_GAIN), 1e-2)
            check(ch, "KILL 2: trine self-cal accumulator NOT within 10x of discard (self/discard <= 0.1)",
                  max(forks[str(c)]['self_over_discard'] for c in CASES_GAIN), 0.1)
        ex5 = g[CASES_GAIN[0]]['acc_fixed']['infidelity']-g[CASES_GAIN[0]]['acc_self']['infidelity']
        ex20 = g[CASES_GAIN[2]]['acc_fixed']['infidelity']-g[CASES_GAIN[2]]['acc_self']['infidelity']
        expo = float(np.log(ex20/ex5)/np.log(2.)) if ex5 > 0 and ex20 > 0 else None
        det_change = g[(EPS0, 1e-4, 0.)]['acc_fixed']['infidelity']/g[(EPS0, 1e-4, 0.)]['discard_unconditional']
        res['words'][kind] = dict(rows=rows, U_cal=cal['u'], forks=forks, per_dump_over_discard_at_eps0=ratio_pd,
                                  fixed_cal_excess=dict(at_5e3=ex5, at_2e2=ex20, exponent_in_eps_minus_eps0=expo),
                                  detuning_case_acc_fixed_over_discard=det_change)
        print(f"forks: {json.dumps(forks, indent=1)}\nper-dump/discard at eps0: {ratio_pd:.6f}   fixed-cal excess: {ex5:.3e}, {ex20:.3e}, exponent {expo}", flush=True)

    files = ['experiments/2026-09-06/rh2.py', 'rlq/rh1_common.py', 'rlq/reflection_loop_composite.py', 'rlq/reflection_loop_dynamics.py',
             'rlq/reflection_loop_finite_dump.py', 'rlq/reflection_loop_reference_echo.py',
             'docs/results/2026-09-06/2026-09-06-RH-2-spec.md', 'docs/prereg/RH-2/DESIGNER-PREDICTIONS.md']
    res['provenance'] = dict(test='rh2', date='2026-09-06', runtime_s=round(time.time()-t0, 1),
        head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
        blinding='designer-run: spec and predictions written and sealed by the same agent (Claire) that wrote and ran this runner; '
                 'sealed BEFORE the runner existed (commit ffa01f3) but NOT independent.',
        sha256={p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in files})
    res['check_count'] = len(ch); res['checks_passed'] = int(sum(c['passed'] for c in ch))
    res['summary'] = {k: dict(forks=v['forks'], fixed_cal_excess=v['fixed_cal_excess']) for k, v in res['words'].items()}
    ap = argparse.ArgumentParser(); ap.add_argument('--json', nargs='?', const='docs/receipts/rh2-checks.json', default=None)
    a = ap.parse_args()
    if a.json:
        Path(a.json).write_text(json.dumps(res, default=norm_json, indent=2, allow_nan=False)+'\n')
    print(f"\nRESULT: {res['checks_passed']}/{res['check_count']} checks passed.  runtime {res['provenance']['runtime_s']}s", flush=True)
    failed = [c['name'] for c in ch if not c['passed']]
    if failed:
        print("FAILED:", *failed, sep='\n  ', flush=True); sys.exit(1)


if __name__ == '__main__':
    run()

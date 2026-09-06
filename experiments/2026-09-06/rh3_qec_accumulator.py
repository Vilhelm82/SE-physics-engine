"""RH-3 -- the accumulator channel inside Codex's QEC distance stack: the decoder-level comparison against FD.

The FD operating row (n = 1, gain 1%, background gamma/a = 1e-6, kappa/a = 10, kappa t_d = 40) supplies the exact
operators K (survival map), E_h (flagged/absorbed effect on P) and E_u (background effect).  FD's channel (Codex,
`native_channel`) is  Phi_FD(rho) = A rho A^dag + tr((E_u + B^dag B) rho) I/2  with herald h = tr(E_h)/2,
A = G K_PP, B = K_QP.  The ideal accumulator captures the absorbed amplitude coherently instead of flagging it;
the recovered branch's Kraus operator is G U^dag F = G sqrt(E_h) (the isometry cancels), so relative to the target
it is sqrt(E_h).  With register error r modelled as I/2 replacement of the recovered branch (r = 0 perfect register;
r = 1 register useless AND unheralded, which is strictly worse than FD's heralded discard):

    Phi_acc(rho) = A rho A^dag + (1 - r) sqrt(E_h) rho sqrt(E_h) + tr((E_u + B^dag B + r E_h) rho) I/2,  herald = 0.

Exact Pauli twirl: w_P = |tr(P A)/2|^2 + (1 - r) |tr(P sqrt(E_h))/2|^2 + tr(E_u + B^dag B + r E_h)/8.
Everything else in the stack (codes, decoders, stack noise p = 1e-3, seeds) is Codex's, unchanged.
DECLARED: ideal self-calibrated capture; register error = depolarising replacement of the recovered branch; no
register duration; background loss unchanged (the accumulator never sees it).  Written 2026-09-06 by Claire.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))  # repo root on sys.path (reorg 2026-09-06)
import argparse, hashlib, json, subprocess, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import sqrtm
import rlq.qec_distance_stack as Q


def accumulator_channel(row, r):
    k, eh, eu = map(Q.unpack, (row['survival_code_map'], row['dump_effect'], row['background_effect']))
    a, b = Q.G@k[2:], k[:2]
    eh = (eh + eh.conj().T)/2; seh = sqrtm(eh); seh = (seh + seh.conj().T)/2
    lost = eu + b.conj().T@b + r*eh
    rl = float(np.trace(lost).real/2)
    w = (np.abs(np.einsum('pij,ji->p', Q.PAULIS, a)/2)**2 + (1 - r)*np.abs(np.einsum('pij,ji->p', Q.PAULIS, seh)/2)**2 + rl/4)
    assert abs(w.sum() - 1) < 1e-10, w.sum()
    assert np.linalg.eigvalsh(lost).min() > -1e-13
    # exact-twirl verification on the Pauli basis, as verify_twirl does for FD
    err = 0.
    for rho in Q.PAULIS:
        direct = sum(p@(a@(p@rho@p)@a.conj().T + (1 - r)*seh@(p@rho@p)@seh + np.trace(lost@p@rho@p)*np.eye(2)/2)@p for p in Q.PAULIS)/4
        expected = sum(wp*p@rho@p for wp, p in zip(w, Q.PAULIS))
        err = max(err, float(np.linalg.norm(direct - expected)))
    w[0] = 1 - w[1:].sum()
    return dict(herald=0., pauli=w.tolist(), register_error=float(r), recovered_fraction=float(np.trace(eh).real/2),
                replaced_unheralded=rl, unheralded_pauli_rate=float(w[1:].sum()),
                recovered_branch_floor=float((np.abs(np.einsum('pij,ji->p', Q.PAULIS, seh)/2)**2)[1:].sum()),
                twirl_error=err, total_time=row['total_time'], n=row['n'], gain=row['gain'], background_rate=row['background_rate'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shots', type=int, default=50000); ap.add_argument('--bb-shots', type=int, default=0)
    ap.add_argument('--distances', type=int, nargs='+', default=[3, 5, 7]); ap.add_argument('--basis', nargs='+', default=['X', 'Z'])
    ap.add_argument('--r', type=float, nargs='+', default=[0.0, 0.1, 0.3, 1.0]); ap.add_argument('--noise', type=float, default=1e-3)
    ap.add_argument('--seed', type=int, default=20260906); ap.add_argument('--json', type=Path, default=Q.ROOT/'docs/receipts/rh3_qec-checks.json')
    ap.add_argument('--skip-fd', action='store_true')
    args = ap.parse_args(); t0 = time.time()
    row = Q.operating_point(background=1e-6)
    fd = Q.native_channel(row)
    chans = [('fd', fd)] if not args.skip_fd else []
    chans += [(f'acc_r{r:g}', accumulator_channel(row, r)) for r in args.r]
    chans.append(('baseline', dict(herald=0., pauli=[1., 0., 0., 0.], total_time=row['total_time'])))
    print("== per-gate channels (n=1, gain 1%, gamma/a=1e-6) ==", flush=True)
    for name, ch in chans:
        px = sum(ch['pauli'][1:])
        print(f"  {name:10s} herald {ch['herald']:.6e}  unheralded Pauli {px:.6e}  (X,Y,Z)={['%.3e' % v for v in ch['pauli'][1:]]}"
              + (f"  twirl_err {ch['twirl_error']:.1e}  floor {ch['recovered_branch_floor']:.2e}" if 'twirl_error' in ch else ''), flush=True)
    report = dict(scope='RH-3: accumulator channel in the FD QEC stack; ideal coherent capture of E_h with register error r; '
                        'background loss and stack noise unchanged; codes/decoders/seeds are qec_distance_stack.py',
                  stack_noise=args.noise, shots=args.shots, channels={n: c for n, c in chans}, runs=[],
                  provenance=dict(date='2026-09-06', head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                                  sha256={p: hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in ('experiments/2026-09-06/rh3_qec_accumulator.py', 'rlq/qec_distance_stack.py', 'docs/receipts/reflection-loop-finite-dump-checks.json')}))
    codes = [Q.surface_circuit(d, b, args.noise) for d in args.distances for b in args.basis]
    if args.bb_shots: codes += [Q.bicycle_circuit(b, args.noise) for b in args.basis]
    for name, ch in chans:
        for code in codes:
            shots = args.shots if code.name.startswith('rotated') else args.bb_shots
            seed = args.seed + len(report['runs'])*1009
            res = Q.run_case(code, ch, shots, seed); res['scenario'] = name; report['runs'].append(res)
            print(f"  >> {name:10s} {code.name} {code.basis}: failures visible {res['flags_visible']['failures']} / hidden {res['flags_hidden']['failures']} of {shots}  ({res['elapsed_seconds']:.0f}s)", flush=True)
            args.json.write_text(json.dumps(report, indent=2) + '\n')
    print("\n== summary: failures per shots (flags visible) ==", flush=True)
    for name, _ in chans:
        line = [f"{r['code'][-2:]}{r['basis']}={r['flags_visible']['failures']}" for r in report['runs'] if r['scenario'] == name]
        print(f"  {name:10s} " + '  '.join(line), flush=True)
    print(f"\nDONE runtime {time.time()-t0:.0f}s", flush=True)


if __name__ == '__main__':
    main()

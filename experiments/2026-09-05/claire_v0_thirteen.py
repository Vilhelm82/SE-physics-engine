#!/usr/bin/env python3
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))  # repo root on sys.path (reorg 2026-09-06)
# =============================================================================
# CLAIRE-V0-13 -- a thirteen-block reflection word with the eps^2 leakage column removed (v = 0), Will's eight
#   first-order conditions, and his second-order logical detuning conditions, in his own parameterisation.
#   (Claire, 2026-09-05.  Uses Will's evaluators -- reflection_loop_detuning_order, _reference_echo, _compression --
#   as the judge; the design is the input.  No existing doc or runner is edited.)
# IDEA: CM-3 eq. (10): the composite word's eps^2 bright column is l_n v with v = sum_j sigma_j (cos beta_j, sin beta_j).
#   The eleven-block has v = (0.6887, 0), so its gain infidelity is (15/2) b_n^2 theta^2 |v|^2 eps^4 = 8.46 eps^4.
#   Impose v = 0 as one more real equation (v_y = 0 is automatic for antisymmetric angles).  The eleven-block's system
#   is square, so this needs the two extra blocks of the thirteen; there they compete with Will's second-order
#   detuning conditions (logical_second = 0), and the two-parameter slack is spent on total stretch W.
# RESULTS (n = 1, nu = 10 a^2, Will's evaluators):
#   v0-only thirteen (his first-order set + v = 0, 10 unknowns):  W = 15.511, gain ~ eps^8 (2.2e-20 at 1e-3),
#     six generators cancelled, but A_Delta = 1.13e7 -- 8x worse than the eleven.  The two blocks bought gain, not Delta.
#   v0 + logical_second thirteen (his 13-unknown parameterisation):  gain ~ eps^8 (1.9e-20 at 1e-3), six generators
#     cancelled, order-4 detuning with A_Delta = 1.34e5, JOINT (eps = 1e-3, Delta/a = 1e-4) = 1.31e-11 at duration 1278/a
#     and exposure 72/a -- beats the nested echo (4.73e-11, 1386/a, 78.3/a) on all three axes with no wrapper; 1.6x the
#     raised wrapper's error at 1.8x less duration and exposure.  Residual error is entirely detuning.
# WHAT IT SAYS ABOUT HIS THIRTEEN: logical_second = 0 does not fix A_Delta by itself (this word has it and A_Delta = 1.3e5
#   vs his 4391), so his coefficient also relies on a small second-order Delta LEAKAGE column, which no condition
#   constrains.  Killing that column is two more real conditions: a fifteen-block target.
# TIER: [DERIVED | Will's evaluators as stated]; the existence certificate (interval Newton) is NOT done here -- the
#   centre is a double-precision root with residual ~1e-11.  KILL: any of the printed checks failing on re-run.
# =============================================================================
import numpy as np, math, json, time, sys
from scipy.optimize import least_squares
import rlq.reflection_loop_detuning_order as D
from rlq.reflection_loop_reference_echo import metrics, response_integrals
from rlq.reflection_loop_compression import noise_generators, bad_response

par_of = lambda x: tuple(D.expand_palindrome(x))
def full_resid(x):
    r = list(D.reduced_conditions(x))
    L = D.logical_second(par_of(x), 1); r += [L[0, 0].real, L[0, 1].real, L[0, 1].imag]
    b, a, w, sig = D.expand_palindrome(x); r.append(sig@np.cos(b))
    return np.array(r)
LB = np.array([-np.inf]*6 + [1.0]*7); UB = np.array([np.inf]*13)

def search(n_starts=30, seed=11, log=print):
    rng = np.random.default_rng(seed); sols = []; t0 = time.time()
    for k in range(n_starts):
        x0 = np.concatenate([rng.uniform(-math.pi, math.pi, 6), rng.uniform(1.0, 1.7, 7)])
        g = lambda x: np.concatenate([full_resid(x), [0.01*(D.expand_palindrome(x)[2].sum() - 12)]])
        s1 = least_squares(g, x0, bounds=(LB, UB), xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=1200)
        s2 = least_squares(full_resid, s1.x, bounds=(LB, UB), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=800)
        r = abs(full_resid(s2.x)).max(); w = D.expand_palindrome(s2.x)[2]
        log(f"start {k:2d}: residual {r:.1e}  W {w.sum():.4f}  min w {w.min():.4f}  [{time.time()-t0:.0f}s]")
        if r < 1e-9: sols.append((float(w.sum()), s2.x, float(r)))
    sols.sort(key=lambda t: t[0]); return sols

def evaluate(x, log=print):
    par = par_of(x); W = D.expand_palindrome(x)[2].sum(); out = dict(W=W, centre=[float(v) for v in x])
    out['residual'] = float(abs(full_resid(x)).max())
    for eps in (1e-2, 1e-3): out[f'gain_only_{eps:g}'] = metrics(D.gain_word(eps, 1, par))['average_infidelity']
    out['gain_order'] = math.log10(out['gain_only_0.01']/out['gain_only_0.001'])
    st = D.stages_from_parameters(par, 1); mom = response_integrals(st, noise_generators())
    out['six_generators_max'] = [float(np.max(np.abs(bad_response(m)))) for m in mom['moments']]
    out['endpoint_vs_composite'] = float(abs(mom['endpoint'] - D.composite()).max())
    deltas = [1e-4, 2e-4, 4e-4]
    U = D.finite_words(deltas, error=0.0, parameter_tuple=par, n=1, steps=120); inf = [metrics(u)['average_infidelity'] for u in U]
    out['detuning_only'] = inf; out['detuning_order'] = math.log(inf[2]/inf[0])/math.log(4); out['A_Delta'] = inf[0]/1e-16
    U = D.finite_words([1e-4], error=1e-3, parameter_tuple=par, n=1, steps=120); out['joint_1e-3_1e-4'] = metrics(U[0])['average_infidelity']
    out['duration_est'] = W*57.7302; out['exposure_est'] = W*3.26165
    for k, v in out.items():
        if k != 'centre': log(f"  {k}: {v}")
    return out

if __name__ == "__main__":
    n_starts = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    sols = search(n_starts)
    print(f"{len(sols)} admissible solutions (residual < 1e-9, all stretches >= 1)")
    if not sols: sys.exit(1)
    W, x, r = sols[0]; print(f"BEST W = {W:.6f}"); out = evaluate(x)
    out['all_admissible_W'] = [s[0] for s in sols]
    json.dump(out, open('docs/receipts/claire-v0-thirteen-checks.json', 'w'), indent=1)
    print("reference: his thirteen joint 3.739e-10 (A_Delta 4391, 948/a, 53.6/a); eleven 1.531e-10 (686/a, 38.7/a); nested echo 4.735e-11 (1386/a, 78.3/a); raised 8.41e-12 (2341/a, 132/a)")

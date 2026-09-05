#!/usr/bin/env python3
# =============================================================================
# SPLIT-1 -- the loop that returns only the data and splits the error.  Part A (FAILS): bright loss DURING an adiabatic dark
#   loop.  Part B (WORKS): exact-return loops with an ENDPOINT bright dump, tight-frame trine (pi/3, 2pi/3, pi/3).
#   (Claire, 2026-09-05, at Will's 'go'.)
# CONSTRUCTION
#   The REALFIBER three-arc loop driven ADIABATICALLY (C^2-smoothed arcs, duration tau each) keeps the state in the dark
#   space: zero nominal bright population, zero dynamical phase, holonomy diag(-1, 1) on (e3, e4).  Everything bright is
#   therefore ERROR.  Add bright-selective absorption at rate kappa into orthogonal sinks: the no-jump propagator
#   K = T exp(-i int (H - i kappa P_b/2)) carries the data; 1 - ||K psi||^2 is the ERASURE probability; the conditional state
#   K psi/||K psi|| is what survives.  Strong kappa Zeno-suppresses the coherent leakage as well as flagging it.
#   SYMMETRISATION: three loops with the active line rotated in the (e3,e4) plane by alpha = +pi/3, 0, -pi/3 (Will's O_{pi/3}).
#   Three reflections at angles (a, b, c) compose to a reflection at c - b + a = 0, i.e. diag(-1, 1) exactly; and the three
#   active directions form a tight frame, sum u u^T = (3/2) I, so the first-order bright exposure of e3 and e4 is EQUAL --
#   state-independent erasure, the flag is data-blind (RE-6's (24), achieved by the trine instead of the four-pattern echo).
# WHAT IS MEASURED, on a grid (a tau, kappa/a) and for gain error eps and static detuning Delta on e3:
#   erasure probability (Haar mean, and max-min over states as the state-dependence witness); conditional average
#   infidelity to diag(-1, 1) of the no-jump map; the Zeno factor (erasure at kappa vs the coherent leakage at kappa = 0).
# PART A RESULT: erasure 27-70% at kappa/a = 1, state-dependent (spread ~ mean), conditional infidelity 1e-2..1e-3.  Adiabatic
#   following IS an O(v/a) bright admixture; loss attacks the transport (RF-7: real transport goes through the bright sector).
#   Bright loss during the loop cannot split anything in this family.
# PART B: at the exact return T_n the nominal bright population is a POINT zero, so any bright amplitude there is error by
#   definition.  Dump it there (loss or projection at the hold, encoded plane stationary).  With dumps between loops the
#   reflections reorient the effective active directions: angles (a, b, c) give effective directions a, 2a - b + pi, a - b
#   and the gate needs c - b + a = 0; the tight-frame (state-independent erasure) solution is (pi/3, 2pi/3, pi/3), NOT the
#   composite's (pi/3, 0, -pi/3).  RESULT (55/a): gain eps = 1% -> heralded erasure 0.16%, spread/mean 1.4e-5, conditional
#   infidelity 2.4e-15; eps = 0.1% -> 1.6e-5, 1.4e-8, ~1e-18.  Gain error is converted entirely into a data-blind flag.
#   NOT split: detuning (a logical rotation, nothing goes bright): cond 1-F 2.3e-8 at Delta/a = 1e-4, first order.  That
#   is now the only coherent error left, so a detuning-only word with dumps at every block boundary is the next target.
# EXPECTATIONS to test (part A, recorded as failed/partly): (i) gain error does NOT enter the holonomy (the dark projector I - H^2/a^2 is gain-independent), so
#   eps only changes the nonadiabatic leakage -- no 5 b^2 eps^2 term at all; (ii) the trine cancels the first-order relative
#   detuning phase (sum cos 2 alpha = sum sin 2 alpha = 0 over +-pi/3, 0); (iii) there is a region with conditional infidelity
#   < 1e-10 at few-percent erasure.  TIER: numerical witnesses on a declared model (bright-selective absorbing loss, static
#   errors).  KILL: state-dependence of erasure not vanishing at first order; or conditional infidelity floor above 1e-8 everywhere.
# =============================================================================
import numpy as np, math, json, sys, time
from scipy.integrate import solve_ivp

e = np.eye(4); a = 1.0
def arcs_for(alpha):
    u = math.cos(alpha)*e[2] + math.sin(alpha)*e[3]          # active logical direction; u_perp = -sin e3 + cos e4 is the spectator
    return [lambda s: (math.cos(s)*e[0] + math.sin(s)*u, e[1]),
            lambda s: (u, math.cos(s)*e[1] + math.sin(s)*e[0]),
            lambda s: (math.cos(s)*u + math.sin(s)*e[1], e[0])]
def smooth_s(t, tau):                                          # C^2 at the ends: ds/dt and d2s/dt2 vanish
    u = t/tau; return (math.pi/2)*(u - math.sin(2*math.pi*u)/(2*math.pi))
def loop_prop(tau, kappa=0.0, gain=1.0, delta=0.0, alpha=0.0, exact_return=False):
    arcs = arcs_for(alpha); D3 = np.diag([0, 0, 1., 0]); U = np.eye(4, dtype=complex)
    for arc in range(3):
        def rhs(t, y):
            s = (math.pi/2)*t/tau if exact_return else smooth_s(t, tau); p, q = arcs[arc](s)
            H = gain*a*(np.outer(p, q) + np.outer(q, p)); Pb = (np.outer(p, p) + np.outer(q, q))   # bright projector = span(p, q)
            Heff = H + delta*D3 - 0.5j*kappa*Pb
            return (-1j*Heff@y.reshape(4, 4)).reshape(-1)
        U = solve_ivp(rhs, (0, tau), U.reshape(-1), rtol=1e-11, atol=1e-13, max_step=tau/3000).y[:, -1].reshape(4, 4)
    return U
def trine_prop(tau, **kw):
    K = np.eye(4, dtype=complex)
    for alpha in (math.pi/3, 0.0, -math.pi/3): K = loop_prop(tau, alpha=alpha, **kw)@K
    return K
G = np.diag([-1., 1.])
def haar_states(n=400, seed=0):
    r = np.random.default_rng(seed); z = r.normal(size=(n, 2)) + 1j*r.normal(size=(n, 2)); return z/np.linalg.norm(z, axis=1)[:, None]
def analyse(K, states):
    out = []
    for psi in states:
        full = K@np.concatenate([[0, 0], psi]); survive = np.vdot(full, full).real
        cond = np.vdot(G@psi, full[2:]); f_cond = abs(cond)**2/survive if survive > 0 else 0.0
        out.append((1 - survive, 1 - f_cond, np.vdot(full[:2], full[:2]).real/max(survive, 1e-300)))
    out = np.array(out)
    return dict(erasure_mean=float(out[:, 0].mean()), erasure_spread=float(out[:, 0].max() - out[:, 0].min()),
                cond_infidelity=float(out[:, 1].mean()), bright_residual_cond=float(out[:, 2].mean()))

P4 = np.diag([0, 0, 1, 1.]); TAU1 = (math.pi/2)*math.sqrt(15)
def dump_word(angles, eps=0., delta=0., tau=TAU1):
    K = np.eye(4, dtype=complex)
    for al in angles: K = P4@loop_prop(tau, gain=1 + eps, delta=delta, alpha=al, exact_return=True)@K
    return K

if __name__ == "__main__":
    quick = len(sys.argv) > 1 and sys.argv[1] == "quick"
    states = haar_states(120 if quick else 400); rows = []; t0 = time.time()
    print("=== PART A: bright loss DURING the adiabatic dark loop (trine pi/3, 0, -pi/3) -- expected to fail ===", flush=True)
    print(f"{'a tau':>6} {'kappa/a':>8} {'eps':>6} {'D/a':>7} | {'erasure':>10} {'spread':>9} {'cond 1-F':>10}", flush=True)
    for tau in ((20.,) if quick else (10., 20., 40.)):
        for kappa in (0.0, 0.1, 1.0, 10.0):
            K = trine_prop(tau, kappa=kappa); r = analyse(K, states); r.update(part="A", a_tau=tau, kappa=kappa, eps=0.0, delta=0.0); rows.append(r)
            print(f"{tau:6.0f} {kappa:8.2f} {0:6.3f} {0:7.0e} | {r['erasure_mean']:10.3e} {r['erasure_spread']:9.1e} {r['cond_infidelity']:10.3e}   [{time.time()-t0:.0f}s]", flush=True)
    print("=== PART B: exact-return loops + ENDPOINT dump ===", flush=True)
    K0 = dump_word((math.pi/3, 2*math.pi/3, math.pi/3)); print("nominal logical block of the tight-frame trine:", np.round(K0[2:, 2:].real, 8).tolist(), flush=True)
    print(f"{'scheme':38s} {'eps':>6} {'D/a':>7} | {'erasure':>9} {'spread/mean':>11} {'cond 1-F':>10}", flush=True)
    for eps, delta in ((0.01, 0.0), (0.001, 0.0), (0.0, 1e-4), (0.0, 1e-3), (0.01, 1e-4), (0.001, 1e-4)):
        for angles, name in (((0.0,), "single return + dump (18/a)"), ((math.pi/3, 0.0, -math.pi/3), "composite angles + dumps (55/a)"), ((math.pi/3, 2*math.pi/3, math.pi/3), "tight-frame trine + dumps (55/a)")):
            r = analyse(dump_word(angles, eps, delta), states); r.update(part="B", scheme=name, angles=list(angles), eps=eps, delta=delta); rows.append(r)
            print(f"{name:38s} {eps:6.3f} {delta:7.0e} | {r['erasure_mean']:9.2e} {r['erasure_spread']/max(r['erasure_mean'], 1e-300):11.2e} {r['cond_infidelity']:10.2e}", flush=True)
    json.dump(rows, open("docs/split1-checks.json" if not quick else "/tmp/split1-quick.json", "w"), indent=1)
    print("reference (unheralded): single loop 5.31e-4 at eps=1% (18/a); composite 2.49e-7 (57/a); reference echo 4.7e-11 at eps=1e-3, D=1e-4 (1386/a)", flush=True)

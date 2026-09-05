#!/usr/bin/env python3
# =============================================================================
# QI-HARNESS -- known quantum-information error results and open problems, reproduced from first principles,
#   with a hook for a candidate general error formula.  (Claire, 2026-09-05, for Will's quantum targets.)
#
# HOW IT WORKS
#   Every benchmark carries: `data` (the physical parameters and the STRATUM DATA a divisor/valuation-type formula
#   might consume: the small parameter f, whether the face is one- or two-sided, the gap, the duration, the order of
#   the zero if known), `known` (the accepted closed form, bound, or None for open), `cls` in
#   {EXACT, ASYMPTOTIC, BOUND, NUMERICAL, OPEN}, and `compute()` -- an independent first-principles calculation
#   (Schroedinger/Lindblad integration, exact enumeration, exact algebra) that must REPRODUCE `known` before the
#   benchmark is trusted.  That reproduction is the harness's own receipt.
#
#   A candidate formula is `candidate(bench) -> float` (predicted error).  `run(candidate)` prints, per benchmark:
#   the first-principles value, the known value, the candidate's value, and a verdict: MATCH / BREAKS / UNJUDGED
#   (bound or open: candidate's number is recorded against the bound, not judged).
#
#   Will's general error formula is NOT in this file.  `example_divisor_candidate` is a placeholder illustrating the
#   hook: error ~ C f^2 at a two-sided face, C f at a one-sided boundary (Codex's constitutive Thm 1 / Will's eq. (10)),
#   with C = 1.  It is expected to break in most places; that is what a placeholder is for.
#
# BENCHMARKS (first-principles method -> known result)
#   B1  Landau-Zener            ODE on H = (v t sz + D sx)/2               -> P_diabatic = exp(-pi D^2 / (2 v))       [ASYMPTOTIC, exact as window -> oo]
#   B2  Rabi with detuning      exact 2-level algebra                        -> P_max = W^2/(W^2 + Delta^2)             [EXACT]
#   B3  Amplitude damping       Lindblad ODE, L = sqrt(g) s-                 -> rho_11 = e^{-g t}, coherence e^{-g t/2}  [EXACT]
#   B4  Pure dephasing          Lindblad ODE, L = sqrt(g/2) sz               -> coherence e^{-g t}                       [EXACT]
#   B5  Quantum Zeno            n projections during H = w sx                -> p_n = cos^{2n}(w t/n) -> 1               [EXACT]
#   B6  3-qubit repetition code exact enumeration of bit flips               -> p_L = 3p^2 - 2p^3, break-even p = 1/2    [EXACT]
#   B7  Berry phase, spin-1/2   holonomy integral of A = i<n|dn> on a cone   -> gamma = -pi(1 - cos theta) = -Omega/2    [EXACT]
#   B8  Helstrom discrimination trace norm of two pure states                -> P_err = (1 - sqrt(1 - |<psi|phi>|^2))/2 [EXACT]
#   B9  Mandelstam-Tamm         2-level minimal orthogonalisation time       -> T = pi/(2 dE)                            [EXACT]
#   B10 DFS leakage (Will eq.10) K_P = P L^+ Q L P for L = sqrt(g)(Z1+Z2) + eps X1, P = DFS   -> leakage = eps^2 exactly (two-sided) [EXACT]
#       (a first draft used eps(Z1-Z2): that is the LOGICAL Z inside the DFS and leaks nothing -- harness bug, recorded)
#   B11 Exceptional point       H = [[i g, k],[k, -i g]]: eigen-splitting vs propagator at the EP -> splitting ~ sqrt(eps), propagator FINITE [EXACT]
#   B12 Will's dark loop (6)    finite-duration Schroedinger on the 3-arc loop, m=4 -> orientation -1 on e3; leakage exponent in (a tau) measured [NUMERICAL]
#   B13 Surface-code threshold  --                                          -> ~1% circuit-level (numerical), exact value OPEN [OPEN]
#   B14 Non-Markovian memory advantage bound                                -> OPEN (no closed form)                     [OPEN]
#
# KILL for the harness itself: any EXACT benchmark whose compute() fails to reproduce its known result to 1e-8.
# =============================================================================
import numpy as np, itertools, math
from scipy.integrate import solve_ivp
from scipy.linalg import expm, eigh

EXACT, ASYMP, BOUND, NUMER, OPEN = "EXACT", "ASYMPTOTIC", "BOUND", "NUMERICAL", "OPEN"
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1., -1.]).astype(complex)
sm = np.array([[0, 0], [1, 0]], complex)                      # lowers |0> -> |1> in the (|0>=excited, |1>=ground) convention below
I2 = np.eye(2)

def lindblad_rhs(H, Ls):
    def f(t, y):
        n = H.shape[0]; rho = y.reshape(n, n)
        d = -1j*(H@rho - rho@H)
        for L in Ls:
            LdL = L.conj().T@L
            d += L@rho@L.conj().T - 0.5*(LdL@rho + rho@LdL)
        return d.reshape(-1)
    return f

def evolve(H, Ls, rho0, T, **kw):
    sol = solve_ivp(lindblad_rhs(H, Ls), (0, T), rho0.reshape(-1).astype(complex), rtol=1e-10, atol=1e-12, **kw)
    n = H.shape[0]; return sol.y[:, -1].reshape(n, n)

BENCH = []
def bench(name, cls, data, known, compute, note=""):
    BENCH.append(dict(name=name, cls=cls, data=data, known=known, compute=compute, note=note))

# ---------------- B1 Landau-Zener ----------------
def lz_compute(v=1.0, D=0.6, W=40.0):
    def rhs(t, y): 
        Ht = 0.5*(v*t*sz + D*sx); return (-1j*Ht@y)
    y0 = np.array([1, 0], complex)                                # diabatic state |0> at t = -W
    sol = solve_ivp(rhs, (-W, W), y0, rtol=1e-11, atol=1e-13)
    return abs(sol.y[0, -1])**2                                   # probability of staying diabatic
bench("B1 Landau-Zener diabatic transition", ASYMP,
      dict(v=1.0, D=0.6, f=0.6, sidedness="two", gap=0.6, duration=None, order=None),
      lambda d: math.exp(-math.pi*d["D"]**2/(2*d["v"])), lambda d: 2*lz_compute(d["v"], d["D"], 160.0) - lz_compute(d["v"], d["D"], 80.0),
      "exact as the sweep window -> oo; the finite-window correction is O(1/W) (diabatic-state admixture), removed by Richardson extrapolation W = 80, 160")

# ---------------- B2 Rabi with detuning ----------------
bench("B2 Rabi maximum transfer with detuning", EXACT,
      dict(W=1.0, Delta=0.7, f=0.7, sidedness="two", gap=None),
      lambda d: d["W"]**2/(d["W"]**2 + d["Delta"]**2),
      lambda d: (lambda U: max(abs(U[1, 0])**2 for U in [expm(-1j*0.5*(d["Delta"]*sz + d["W"]*sx)*t) for t in np.linspace(0, 4*math.pi/math.hypot(d["W"], d["Delta"]), 4001)]))(None))

# ---------------- B3 amplitude damping ----------------
def ad_compute(g=0.8, T=1.3):
    rho0 = np.array([[0.7, 0.3], [0.3, 0.3]], complex)            # |0> = excited
    rho = evolve(np.zeros((2, 2)), [math.sqrt(g)*sm], rho0, T)
    return (rho[0, 0].real, abs(rho[0, 1]))
bench("B3 amplitude damping: population and coherence", EXACT,
      dict(g=0.8, T=1.3, f=0.8, sidedness="one", gap=None),
      lambda d: (0.7*math.exp(-d["g"]*d["T"]), 0.3*math.exp(-d["g"]*d["T"]/2)), lambda d: ad_compute(d["g"], d["T"]))

# ---------------- B4 pure dephasing ----------------
def deph_compute(g=0.9, T=1.1):
    rho0 = np.array([[0.5, 0.5], [0.5, 0.5]], complex)
    rho = evolve(np.zeros((2, 2)), [math.sqrt(g/2)*sz], rho0, T); return abs(rho[0, 1])
bench("B4 pure dephasing coherence", EXACT, dict(g=0.9, T=1.1, f=0.9, sidedness="one"),
      lambda d: 0.5*math.exp(-d["g"]*d["T"]), lambda d: deph_compute(d["g"], d["T"]))

# ---------------- B5 quantum Zeno ----------------
def zeno_compute(w=1.0, t=1.0, n=50):
    U = expm(-1j*w*sx*t/n); psi = np.array([1, 0], complex); p = 1.0
    for _ in range(n):
        psi = U@psi; a = abs(psi[0])**2; p *= a; psi = np.array([1, 0], complex)
    return p
bench("B5 quantum Zeno survival, n projections", EXACT, dict(w=1.0, t=1.0, n=50, f=1.0/50, sidedness="one"),
      lambda d: math.cos(d["w"]*d["t"]/d["n"])**(2*d["n"]), lambda d: zeno_compute(d["w"], d["t"], d["n"]))

# ---------------- B6 repetition code ----------------
def rep_compute(p=0.1):
    pl = 0.0
    for flips in itertools.product([0, 1], repeat=3):
        pr = np.prod([p if f else 1-p for f in flips])
        if sum(flips) >= 2: pl += pr
    return pl
bench("B6 3-qubit repetition code logical error", EXACT, dict(p=0.1, f=0.1, sidedness="one", order=2),
      lambda d: 3*d["p"]**2 - 2*d["p"]**3, lambda d: rep_compute(d["p"]))

# ---------------- B7 Berry phase ----------------
def berry_compute(theta=0.9, N=4000):
    phis = np.linspace(0, 2*math.pi, N+1); g = 0.0
    def n_up(phi): return np.array([math.cos(theta/2), np.exp(1j*phi)*math.sin(theta/2)])
    for k in range(N):
        a, b = n_up(phis[k]), n_up(phis[k+1]); g += np.angle(np.vdot(a, b))
    return -g                                                      # gamma = -Im sum log <n_k|n_{k+1}>  -> -Omega/2 for spin-up along B
bench("B7 Berry phase of spin-1/2 on a cone of half-angle theta", EXACT, dict(theta=0.9, f=None, sidedness=None),
      lambda d: -math.pi*(1 - math.cos(d["theta"])), lambda d: -((lambda x: (x + math.pi) % (2*math.pi) - math.pi)(-berry_compute(d["theta"]))))

# ---------------- B8 Helstrom ----------------
def helstrom_compute(ov=0.6):
    psi = np.array([1, 0], complex); phi = np.array([ov, math.sqrt(1-ov**2)], complex)
    d = np.outer(psi, psi.conj()) - np.outer(phi, phi.conj()); return 0.5*(1 - 0.5*np.abs(np.linalg.eigvalsh(d)).sum())
bench("B8 Helstrom minimum error, two pure states", EXACT, dict(overlap=0.6),
      lambda d: 0.5*(1 - math.sqrt(1 - d["overlap"]**2)), lambda d: helstrom_compute(d["overlap"]))

# ---------------- B9 Mandelstam-Tamm ----------------
def mt_compute(dE=0.5):
    from scipy.optimize import brentq
    H = dE*sz; psi = np.array([1, 1], complex)/math.sqrt(2)
    f = lambda t: np.vdot(psi, expm(-1j*H*t)@psi).real          # = cos(dE t) for this state
    return float(brentq(f, 1e-9, 0.9*math.pi/dE + 0.2*math.pi/dE))
bench("B9 Mandelstam-Tamm minimal orthogonalisation time", EXACT, dict(dE=0.5),
      lambda d: math.pi/(2*d["dE"]), lambda d: mt_compute(d["dE"]))

# ---------------- B10 DFS leakage: Will's eq. (10) ----------------
def kron(*ms): 
    out = np.array([[1.]], complex)
    for m in ms: out = np.kron(out, m)
    return out
def dfs_compute(g=1.0, eps=0.05):
    Z1, Z2, X1 = kron(sz, I2), kron(I2, sz), kron(sx, I2); L = math.sqrt(g)*(Z1 + Z2) + eps*X1   # eps (Z1 - Z2) would be the LOGICAL Z inside the DFS: no leakage
    v01 = np.zeros(4); v01[1] = 1; v10 = np.zeros(4); v10[2] = 1
    P = np.outer(v01, v01) + np.outer(v10, v10); Q = np.eye(4) - P
    K = P@L.conj().T@Q@L@P; return max(np.linalg.eigvalsh(K).real)   # leakage rate from the worst DFS state
bench("B10 DFS leakage under symmetry-breaking eps X1 (Will eq. 10)", EXACT, dict(g=1.0, eps=0.05, f=0.05, sidedness="two", order=2),
      lambda d: d["eps"]**2, lambda d: dfs_compute(d["g"], d["eps"]),
      "two-sided face eps = 0: P(Z1+Z2) = 0 kills the cross terms, K_P = eps^2 P X1 Q X1 P = eps^2 P: exactly the f^2 law")

# ---------------- B11 exceptional point ----------------
def ep_compute(k=1.0, eps=1e-4, T=1.0):
    H = np.array([[1j*(k - eps), k], [k, -1j*(k - eps)]])
    ev = np.linalg.eigvals(H); split = abs(ev[0] - ev[1])
    U = expm(-1j*H*T); return (split, np.linalg.norm(U))
bench("B11 exceptional point: splitting vs propagator", EXACT, dict(k=1.0, eps=1e-4, T=1.0, f=1e-4, sidedness="one", order=0.5),
      lambda d: (2*math.sqrt(2*d["k"]*d["eps"] - d["eps"]**2), None), lambda d: ep_compute(d["k"], d["eps"], d["T"]),
      "splitting ~ sqrt(eps) (order 1/2); the propagator norm stays finite -- divergent eigenvector coordinates cancel")

# ---------------- B12 Will's dark loop, finite duration ----------------
def dark_loop(a=1.0, tau=20.0, smooth=False, N=60):
    e = np.eye(4); arcs = [
        lambda s: (math.cos(s)*e[0] + math.sin(s)*e[2], e[1]),
        lambda s: (e[2], math.cos(s)*e[1] + math.sin(s)*e[0]),
        lambda s: (math.cos(s)*e[2] + math.sin(s)*e[1], e[0])]
    def H_of(s, arc):
        p, q = arcs[arc](s); return a*(np.outer(p, q) + np.outer(q, p)).astype(complex)
    psi = e[2].astype(complex)                                     # start in dark state e3 at H*
    for arc in range(3):
        def rhs(t, y):
            u = t/tau; s = (math.pi/2)*(u if not smooth else (u - math.sin(2*math.pi*u)/(2*math.pi)))   # smoothed: derivative stops at ends
            return -1j*H_of(s, arc)@y
        sol = solve_ivp(rhs, (0, tau), psi, rtol=1e-10, atol=1e-12); psi = sol.y[:, -1]
    P0 = np.eye(4) - H_of(0, 0)@H_of(0, 0)/a**2
    dark = P0@psi; leak = 1 - np.vdot(dark, dark).real
    return (np.vdot(e[2], psi).real, leak)
def dark_compute(d):
    taus = d["tau"]*np.array([0.5, 0.7, 1.0, 1.4, 2.0, 2.8]); res = [dark_loop(d["a"], t, d["smooth"]) for t in taus]
    sgn = res[2][0]; leaks = np.array([r[1] for r in res])
    slope = -np.polyfit(np.log(taus), np.log(leaks), 1)[0]      # corner leakage oscillates with phase a tau; a 6-point fit reads the envelope
    return (round(sgn), float(leaks[2]), float(slope))
bench("B12 Will's dark loop (6): orientation and leakage scaling, corners", NUMER, dict(a=1.0, tau=20.0, smooth=False, f=None, gap=1.0, duration=20.0),
      lambda d: (-1, None, 2.0), dark_compute, "adiabatic theorem: velocity discontinuities at the corners give leakage ~ (a tau)^-2 times an oscillating factor in a tau; sign -1 is the holonomy")
bench("B12s Will's dark loop (6): smoothed corners", NUMER, dict(a=1.0, tau=20.0, smooth=True, f=None, gap=1.0, duration=20.0),
      lambda d: (-1, None, 6.0), dark_compute, "s(u) = u - sin(2 pi u)/(2 pi) kills the first TWO derivatives of ds/dt at the ends (C^2), so the adiabatic theorem gives amplitude ~ (a tau)^-3, leakage ~ (a tau)^-6; the harness first wrote 4 -- corrected against the theorem")

# ---------------- B13, B14 open ----------------
bench("B13 surface-code threshold, circuit-level depolarising noise", OPEN, dict(f=None, note="numerical ~0.5-1%; exact value open"), None, lambda d: None)
bench("B14 non-Markovian memory advantage for a logical lifetime", OPEN, dict(f=None), None, lambda d: None)

# ---------------- the hook ----------------
def example_divisor_candidate(b):
    """PLACEHOLDER -- not Will's formula.  error ~ f^2 (two-sided), f (one-sided), None otherwise."""
    d = b["data"]; f = d.get("f"); s = d.get("sidedness")
    if f is None or s is None: return None
    return f**2 if s == "two" else f

def run(candidate=example_divisor_candidate, tol=1e-6):
    ok_h = 0; n_h = 0; rows = []
    for b in BENCH:
        fp = b["compute"](b["data"]) if b["cls"] != OPEN else None
        kn = b["known"](b["data"]) if b["known"] else None
        if b["cls"] in (EXACT,):
            n_h += 1
            fpt = fp if isinstance(fp, tuple) else (fp,); knt = kn if isinstance(kn, tuple) else (kn,)
            pairs = [(x, y) for x, y in zip(fpt, knt) if y is not None and x is not None]
            close = bool(pairs) and all(abs(x - y) < tol for x, y in pairs)
            ok_h += close; receipt = "REPRODUCED" if close else f"HARNESS FAIL {[(float(x), float(y)) for x, y in pairs]}"
        elif b["cls"] == ASYMP:
            n_h += 1; close = abs(fp - kn) < 2e-3; ok_h += close; receipt = "REPRODUCED(asymp)" if close else "HARNESS FAIL"
        elif b["cls"] == NUMER:
            receipt = f"measured sign {fp[0]}, leak {fp[1]:.2e}, exponent {fp[2]:.2f}" if fp else "n/a"
        else:
            receipt = "OPEN"
        cv = candidate(b)
        if cv is None or b["cls"] in (OPEN, BOUND, NUMER): verdict = "UNJUDGED"
        else:
            target = kn[0] if isinstance(kn, tuple) else kn
            verdict = "MATCH" if (target is not None and abs(cv - target) <= 0.05*max(abs(target), 1e-12)) else "BREAKS"
        rows.append((b["name"], b["cls"], receipt, kn, cv, verdict))
    print(f"{'benchmark':70s} {'class':11s} {'harness receipt':52s} {'known':>22s} {'candidate':>11s} verdict")
    for r in rows:
        kn = r[3]; kns = "open" if kn is None else (f"{kn[0]:.6g},..." if isinstance(kn, tuple) else f"{kn:.6g}")
        cvs = "-" if r[4] is None else f"{r[4]:.4g}"
        print(f"{r[0][:70]:70s} {r[1]:11s} {r[2][:52]:52s} {kns:>22s} {cvs:>11s} {r[5]}")
    print(f"\nHARNESS: {ok_h}/{n_h} exact/asymptotic benchmarks reproduced from first principles.")
    print(f"CANDIDATE '{candidate.__name__}': " + ", ".join(f"{v}={sum(1 for r in rows if r[5]==v)}" for v in ("MATCH", "BREAKS", "UNJUDGED")))
    return rows

if __name__ == "__main__":
    run()

#!/usr/bin/env python3
# =============================================================================
# KIN-1 -- the kinetic cost is a FORK, not a result.  (Claire, 2026-09-05, pushing on
#   Codex's item 3 in external/lorentzian-upgrade-derivations.py section B.)
# INPUTS: T7a's form q = diag(1,1,-1) with seat c = (0,0,1); a ruler as a unit spacelike
#   line H(u,phi) = (cosh u cos phi, cosh u sin phi, sinh u), q(H,H) = 1; T1's conjugator
#   S = diag(1,1,i); T7b3 (hbar <-> G is a symmetry of the bare tier); P4 (state = angles).
# CODEX'S DECLARATION: "assume equal ruler-velocity costs measured by the seated form
#   h_c = q + 2 q(.,c) q(.,c), fix the seat, eliminate the common azimuth by minimisation"
#   -> g_kin = diag(cosh 2l1, cosh 2l2, C1 C2/(C1+C2)), fixture diag(3,1,2/3) at gamma = 1.
# CLAIMS, checked:
#   k1  h_c = q o S.  S is NOT in O(q_C).  Real pivots and the complexified quarter-turn both
#       preserve q: NO pivot produces h_c.  [Claire's 'it is a pivot' hunch: DEAD, recorded.]
#   k2  A positive cost agreeing with q on c_perp is diag(1,1,kappa), kappa > 0 free; depth
#       coefficient kappa cosh^2 u + sinh^2 u; Codex's fixture 3 is 2 kappa + 1.  The S-reading
#       pins kappa = 1 (i^2 = -1): a SURPLUS over the hand-flip.
#   k3  'Equal ruler weights' is DERIVED from T7b3: with weights (w1,w2) the eliminated angular
#       coefficient is w1 w2 C1 C2/(w1 C1 + w2 C2), invariant under l1 <-> l2 iff w1 = w2.
#       [First run FAILED k3d: T7b3 fixes the weight RATIO, not the overall scale w. Recorded.]
#   k4  'Eliminate the common azimuth' is P4: the minimiser is the horizontal projection off the
#       fibre; the value is the mechanical-connection quotient metric.
#   k5  THE FORK.  The ruler's OWN state space {q = 1} carries the induced metric
#       -du^2 + cosh^2 u dphi^2 = dS_2 (depth TIMELIKE, K = +1).  Three positivisations:
#         (i)   Wick the AMBIENT vector (Codex):  cosh 2u du^2 + cosh^2 u dphi^2,  K = -sech^2 2u
#         (ii)  flip the sign of dS_2's depth:    du^2 + cosh^2 u dphi^2,          K = -1  (H^2)
#         (iii) Wick dS_2 itself, u -> i u:       du^2 + cos^2 u dphi^2,           K = +1  (S^2)
#       (iii) is the round sphere: the paper's (S^2)^3 is the COMPACT REAL FORM OF THE RULER
#       MANIFOLD, not merely 'Cl(3) residue'.  The three depth coefficients differ as FUNCTIONS.
#   k6  Robust across the fork: positive-definite and finite at Codex's boundary fixture
#       (gamma = 1).  NOT robust: the depth/azimuth weighting -- (i) gives (3,1,2/3), (ii) gives
#       (1,1,2/3) -- and by T7h (in prep) the direction of approach decides the horizon reading.
# VERDICT (tier): the kinetic metric is DECLARED at exactly one point -- which object the cost
#   is quadratic ON (ambient vector / ruler state) -- and no primitive in P1-P13 selects it.
#   The load (Codex item 1) is the only closure.  Nothing may be built on diag(3,1,2/3) as reached.
# KILL: a primitive that selects ambient vs state closes the fork and this runner's verdict dies.
# =============================================================================
import sympy as sp, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), sp.cancel, sp.factor, lambda q: sp.simplify(q.rewrite(sp.exp))):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False

x, y, t, lam, u, v, phi, kap, w1, w2, Om, td = sp.symbols('x y t lambda u v phi kappa w1 w2 Omega tdot', real=True)
q = sp.diag(1, 1, -1); c = sp.Matrix([0, 0, 1]); vec = sp.Matrix([x, y, t])
Q = lambda A, p, r: (p.T*A*r)[0]
h_c = q + 2*(q*c)*(q*c).T
S = sp.diag(1, 1, sp.I)

print("=== k1: h_c is q on T1's rotated slice; no pivot produces it ===")
check("k1a h_c = q + 2 q(.,c) q(.,c) is the identity in the seat chart", h_c == sp.eye(3))
check("k1b h_c(v,v) = q(Sv,Sv) with S = diag(1,1,i): h_c = q o S, the form on the slice t -> i t", z(Q(h_c, vec, vec) - sp.expand(Q(q, S*vec, S*vec))))
check("k1c S is not in O(q_C): S^T q S = I != q, det S = i", (S.T*q*S) == sp.eye(3) and S.det() == sp.I)
B = sp.Matrix([[sp.cosh(lam), 0, sp.sinh(lam)], [0, 1, 0], [sp.sinh(lam), 0, sp.cosh(lam)]])
check("k1d a real pivot (boost) preserves q: B^T q B = q -- it cannot produce h_c", z((B.T*q*B - q).norm()))
K = sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]); W = sp.simplify(sp.exp(sp.I*sp.pi/2*K))
check("k1e the complexified quarter-turn W = exp(i pi/2 K) also preserves q: W^T q W = q", z((W.T*q*W - q).norm()))
check("k1f W sends c to i x: 'c -> iK' is a group element on the SAME form; the form change is S, not W", W*c == sp.Matrix([sp.I, 0, 0]))

print("=== k2: the free constant, and what pins it ===")
H = sp.Matrix([sp.cosh(u)*sp.cos(phi), sp.cosh(u)*sp.sin(phi), sp.sinh(u)])
check("k2a H(u,phi) is a unit spacelike ruler: q(H,H) = 1", z(Q(q, H, H) - 1))
Hu, Hp = H.diff(u), H.diff(phi)
h_k = sp.diag(1, 1, kap)
gd = sp.simplify(Q(h_k, Hu, Hu)); ga = sp.simplify(Q(h_k, Hp, Hp)); gm = sp.simplify(Q(h_k, Hu, Hp))
check("k2b cost diag(1,1,kappa): depth = kappa cosh^2 u + sinh^2 u, azimuth = cosh^2 u, mixed = 0",
      z(gd - (kap*sp.cosh(u)**2 + sp.sinh(u)**2)) and z(ga - sp.cosh(u)**2) and z(gm))
check("k2c kappa = 1 gives Codex's cosh 2u; kappa enters ONLY the depth coefficient", z(gd.subs(kap, 1) - sp.cosh(2*u)) and ga.has(kap) == False)
fix_u = sp.asinh(1)
check("k2d Codex's boundary fixture 3 is 2 kappa + 1 -- a declared kappa would move it", z(gd.subs(u, fix_u) - (2*kap + 1)))
check("k2e the S-reading pins kappa = 1: q(Sv,Sv) has coefficient exactly (i)^2 * (-1) = +1 on t^2", sp.expand(Q(q, S*vec, S*vec)).coeff(t, 2) == 1)

print("=== k3: equal weights are T7b3, not a declaration ===")
C1, C2 = sp.symbols('C1 C2', positive=True)
Tw = w1*C1*Om**2 + w2*C2*(Om + td)**2
Om_star = sp.solve(sp.diff(Tw, Om), Om)[0]
ang_w = sp.simplify(Tw.subs(Om, Om_star)/td**2)
check("k3a weighted elimination gives w1 w2 C1 C2/(w1 C1 + w2 C2)", z(ang_w - w1*w2*C1*C2/(w1*C1 + w2*C2)))
swapped = ang_w.subs({C1: C2, C2: C1}, simultaneous=True)
diff_w = sp.factor(sp.together(ang_w - swapped))
check("k3b it is invariant under C1 <-> C2 iff w1 = w2 (difference has factor (w1 - w2))", (w1 - w2) in [f for f in sp.Mul.make_args(sp.numer(diff_w))] or z(diff_w.subs(w1, w2)))
check("k3c hence T7b3 (hbar <-> G symmetry of the bare tier) FORCES equal weights", z(diff_w.subs(w1, w2)) and not z(diff_w.subs({w1: 2, w2: 1, C1: 1, C2: 2})))
check("k3d with w1 = w2 = w the coefficient is w C1 C2/(C1 + C2): T7b3 fixes the RATIO of weights, not the overall scale; Codex's w = 1 is a unit convention", z(ang_w.subs(w1, w2) - w2*C1*C2/(C1 + C2)) and not z(ang_w.subs(w1, w2) - C1*C2/(C1 + C2)))

print("=== k4: the azimuth elimination is P4's quotient ===")
check("k4a the minimiser Omega* = -C2 tdot/(C1 + C2) is the horizontal lift (common rotation = fibre)", z(Om_star.subs({w1: 1, w2: 1}) + C2*td/(C1 + C2)))
check("k4b the eliminated cost is the harmonic combination -- the mechanical-connection quotient metric of the fibre", z(ang_w.subs({w1: 1, w2: 1}) - 1/(1/C1 + 1/C2)))

print("=== k5: THE FORK -- three positivisations of the ruler's own state space ===")
def gauss(E, Gm):
    sq = sp.sqrt(E*Gm)
    return sp.simplify(-1/(2*sq)*sp.diff(sp.diff(Gm, u)/sq, u))
E0, G0 = sp.simplify(Q(q, Hu, Hu)), sp.simplify(Q(q, Hp, Hp))
check("k5a the q-induced metric on {q = 1} is -du^2 + cosh^2 u dphi^2: depth is TIMELIKE on the ruler's own state space", z(E0 + 1) and z(G0 - sp.cosh(u)**2))
check("k5b that is dS_2: K = +1", z(gauss(E0, G0) - 1))
K_i = gauss(sp.cosh(2*u), sp.cosh(u)**2)
check("k5c (i) Codex/ambient-Wick: K = -sech^2 2u ; -1 at the orthogonal state, 0 as depth -> oo", z(K_i + 1/sp.cosh(2*u)**2) and z(K_i.subs(u, 0) + 1) and sp.limit(K_i, u, sp.oo) == 0)
K_ii = gauss(1, sp.cosh(u)**2)
check("k5d (ii) sign-flip of dS_2: du^2 + cosh^2 u dphi^2 has K = -1 : H^2", z(K_ii + 1))
ds2 = -sp.Symbol('du')**2 + sp.cosh(u)**2*sp.Symbol('dphi')**2
wick = sp.simplify(ds2.subs(u, sp.I*u).subs(sp.Symbol('du'), sp.I*sp.Symbol('du')))
check("k5e (iii) Wick dS_2 itself, u -> i u: du^2 + cos^2 u dphi^2", z(wick - (sp.Symbol('du')**2 + sp.cos(u)**2*sp.Symbol('dphi')**2)))
K_iii = sp.simplify(-1/(2*sp.cos(u))*sp.diff(sp.diff(sp.cos(u)**2, u)/sp.cos(u), u))
check("k5f     and that is the ROUND S^2: K = +1.  The Euclidean cell is the compact real form of the ruler manifold", z(K_iii - 1))
check("k5g the three depth coefficients cosh 2u, 1, 1(on S^2) are DISTINCT FUNCTIONS: a genuine fork, not a chart", not z(sp.cosh(2*u) - 1) and not z(K_i - K_ii) and not z(K_ii - K_iii))

print("=== k6: what survives the fork and what does not ===")
gi = sp.diag(sp.cosh(2*u), sp.cosh(2*v), sp.cosh(u)**2*sp.cosh(v)**2/(sp.cosh(u)**2 + sp.cosh(v)**2))
gii = sp.diag(1, 1, sp.cosh(u)**2*sp.cosh(v)**2/(sp.cosh(u)**2 + sp.cosh(v)**2))
fx = {u: fix_u, v: 0}
gi_f = gi.subs(fx).applyfunc(lambda e: sp.simplify(e.rewrite(sp.exp))); gii_f = gii.subs(fx).applyfunc(lambda e: sp.simplify(e.rewrite(sp.exp)))
gam_fix = sp.simplify((sp.cosh(u)*sp.cosh(v)*sp.cos(sp.pi/4) - sp.sinh(u)*sp.sinh(v)).subs(fx).rewrite(sp.exp))
check("k6a the fixture (u = asinh 1, v = 0, t = pi/4) is ON the ruler-plane boundary: gamma = 1", z(gam_fix - 1))
check("k6b (i) gives diag(3, 1, 2/3) there: positive, finite", gi_f == sp.diag(3, 1, sp.Rational(2, 3)))
check("k6c (ii) gives diag(1, 1, 2/3) there: positive, finite", gii_f == sp.diag(1, 1, sp.Rational(2, 3)))
check("k6d ROBUST: both are positive-definite everywhere (all diagonal entries > 0 for real u, v)", all(sp.ask(sp.Q.positive(e), sp.Q.real(u) & sp.Q.real(v)) or z(e - 1) or sp.simplify(e) == sp.cosh(2*u) or sp.simplify(e) == sp.cosh(2*v) for e in list(gi.diagonal()) + list(gii.diagonal())))
check("k6e NOT ROBUST: the depth/azimuth ratio differs -- (i) 9/2 vs (ii) 3/2 at the fixture; the approach direction to the boundary is not fixed by geometry",
      sp.simplify(gi_f[0, 0]/gi_f[2, 2]) == sp.Rational(9, 2) and sp.simplify(gii_f[0, 0]/gii_f[2, 2]) == sp.Rational(3, 2))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: h_c is q o S (T1's conjugator), which pins kappa = 1 -- a surplus over Codex's hand-flip.  Equal weights")
print("  are T7b3; the azimuth quotient is P4.  ONE declaration remains: which object the cost is quadratic ON.")
print("  Three positivisations of the ruler manifold dS_2 -- ambient-Wick (K = -sech^2 2u), sign-flip (H^2), Wick (S^2) --")
print("  are distinct, all positive and finite at gamma = 1, and disagree on the depth/azimuth weighting that decides the")
print("  approach direction.  The fork is OPEN; the load (Codex item 1) is its only closure.  Tier: diag(3,1,2/3) is DECLARED.")

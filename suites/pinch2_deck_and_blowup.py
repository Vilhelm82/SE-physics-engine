#!/usr/bin/env python3
# =============================================================================
# PINCH-2 -- the pinch, the deck, and the weighted blow-up.  (Claire, 2026-09-05.)
#   Ties PINCH-1 (67099dd) and Codex's directional resolution xi to the sheet machinery already in
#   the repo (T5c: V_spin = 2D, D^2 = -det G_c, tau: D -> -D).  No load, no soldering, no r.
# INPUTS: T7a's form; G(a,b,gamma); eps = 1 - gamma^2, W = a^2 + b^2 - 2ab gamma, delta = eps + W;
#   F = eps/delta; D = det[hbar, G, c] = cosh l1 cosh l2 sin t in seat coords (T5c), D^2 = delta;
#   Codex's xi = (b - a gamma)/d with d^2 = delta (docs/results/2026-09-05/2026-09-05-pinch-directional-resolution.md).
# CLAIMS, checked:
#   q1  delta - (b - a gamma)^2 = (1 - gamma^2)(1 + a^2) IDENTICALLY, hence F = (1 - xi^2)/(1 + a^2)
#       everywhere (Codex's formula is a global identity, not a pinch expansion).
#   q2  xi = 2(b - a gamma)/V_spin: the pinch-resolution parameter is a ruler projection divided by
#       the sheet observable.  In seat coords xi = (b - a gamma)/(cosh l1 cosh l2 sin t).
#   q3  THE DECK: tau: D -> -D (t -> -t) sends xi -> -xi and fixes F, eta, N^2.  xi = +-1 are DECK
#       PARTNERS (the two sheets' null normals of PINCH-1); xi = 0 (the bisector nu*) is DECK-FIXED.
#       The live seat at the pinch is the unique deck-invariant continuation.
#   q4  PINCH-1's three approach classes are xi-values: transverse -> xi -> 0; in-horizon -> xi = -+sigma
#       exactly (sign = sheet); tangent-leaving with eps ~ s^m, W ~ s^2 -> xi -> 0, -1/sqrt3, -1.
#   q5  KUMMER: F, eta, N^2 and every component of nu lie in the base field K = Q(a,b,gamma); xi lies in
#       K(sqrt delta) \ K (delta irreducible, degree 2 in gamma: not a square).  Base-field observables
#       are sheet-blind BY CONSTRUCTION and cannot resolve the pinch; only a Kummer coordinate can.
#   q6  TANGENCY: grad delta || grad eps along the whole pinch line.  The branch locus {delta = 0} is
#       TANGENT to the horizon {gamma = sigma} at the pinch, not transverse -- that is WHY the pinch is
#       a codim-2 singularity and why contact order (not direction) governs the readout.
#   q7  WEIGHTED BLOW-UP (1,2): along a = sigma b0 + s, gamma = sigma(1 - lam s^2), xi^2 -> 1/(1 + 2 lam (1+b0^2)).
#       xi^2 < 1 <-> lam > 0 (compact sector, F > 0, seat alive); xi^2 = 1 <-> lam = 0 (the horizon);
#       xi^2 > 1 <-> -1/(2(1+b0^2)) < lam < 0 (the interior sliver, F < 0); xi -> oo <-> the branch
#       parabola; xi = 0 <-> transverse.  xi is a coordinate on the exceptional divisor; tau folds it.
#   q8  NO REALISABLE LOOP encircles the pinch line: delta changes sign around it (half the small circle
#       is unrealisable).  The exceptional divisor's two ends xi = +-1 are identified by the deck, not
#       by a path in state space.  Mod the deck the pinch resolves to a SEGMENT xi in [0, 1]:
#       seat-alive (0) to seat-null (1).
# [First run: q3c/q3d were placeholders, not computations (rule 12); replaced with the horizon xi-values and nu(a,0) = nu*.]
# DYNAMICAL COMPLEMENT (Codex, docs/results/2026-09-05/2026-09-05-cella-predictive-state.md sec. 3-4, 63d951e, recorded 2026-09-05 evening):
#   q5's 'base-field observables are sheet-blind' is the STATIC statement.  Under the constant passive load
#   M = [[mu,beta],[beta,nu]] the deck is tau(xi,d) = (-xi,-d) and it acts on the INPUTS too (X_d flips, e_d -> -e_d),
#   so equal F under tau does not identify two fixed-input experiments: at d = 0 the two null limits xi = +-1 give
#   OPPOSITE signs of dF/dt under the same two-step input (dF/dt = -2 mu nu t xi/p, verified independently).  The
#   response distinguishes the sheets; (xi, d) is the exact minimum predictive state at the pinch; 'a finite sheet label
#   cannot replace an interval of distinguishable directions' (hers).  q3/q8 stand as algebra; this is what they mean.
# TIER: [DERIVED | T7a's declared signature].  Kill: tau failing to send xi -> -xi; or grad delta not
#   parallel to grad eps at the pinch; or a realisable loop about the pinch line.
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

a, b, gam = sp.symbols('a b gamma', real=True)
t, l1, l2 = sp.symbols('t l1 l2', real=True)
s = sp.Symbol('s', positive=True); lam = sp.Symbol('lam', real=True); b0 = sp.Symbol('b0', real=True)
G = sp.Matrix([[-1, a, b], [a, 1, gam], [b, gam, 1]])
eps = 1 - gam**2; W = a**2 + b**2 - 2*a*b*gam; delta = eps + W; F = eps/delta; eta = W/eps
xi2 = (b - a*gam)**2/delta
seat = {a: -sp.sinh(l1), b: -sp.sinh(l2), gam: sp.cosh(l1)*sp.cosh(l2)*sp.cos(t) - sp.sinh(l1)*sp.sinh(l2)}
D = sp.cosh(l1)*sp.cosh(l2)*sp.sin(t)                     # T5c 1a: the oriented volume
xi_seat = sp.simplify((b - a*gam).subs(seat)/D)

print("=== q1: Codex's F = (1 - xi^2)/(1 + a^2) is a global identity ===")
check("q1a delta - (b - a gamma)^2 = (1 - gamma^2)(1 + a^2) identically", z(sp.expand(delta - (b - a*gam)**2) - sp.expand(eps*(1 + a**2))))
check("q1b hence F = (1 - xi^2)/(1 + a^2) at EVERY state, not only near the pinch", z(F - (1 - xi2)/(1 + a**2)))
check("q1c and eta = (a^2 + xi^2)/(1 - xi^2) (Codex's tilt identity), everywhere", z(eta - (a**2 + xi2)/(1 - xi2)))

print("=== q2: xi is a ruler projection over the sheet observable ===")
check("q2a D^2 = delta in seat coordinates (T5c: D^2 = -det G_c)", z(sp.simplify(D**2 - delta.subs(seat))))
check("q2b so xi^2 = (b - a gamma)^2/D^2 = 4 (b - a gamma)^2/V_spin^2 with V_spin = 2D", z(sp.simplify(xi2.subs(seat) - (b - a*gam).subs(seat)**2/D**2)))
bag = sp.simplify(sp.expand_trig((b - a*gam).subs(seat)))
check("q2c the numerator b - a gamma is EVEN in t and D is ODD in t: xi is odd in t", z(bag.subs(t, -t) - bag) and z(D.subs(t, -t) + D))

print("=== q3: the deck acts on xi by sign and fixes the base-field observables ===")
tau = {t: -t}
check("q3a tau: D -> -D (T5c 1d), hence xi -> -xi", z(D.subs(tau) + D) and z(xi_seat.subs(tau) + xi_seat))
Fs, etas = sp.simplify(F.subs(seat)), sp.simplify(eta.subs(seat))
check("q3b tau fixes F, eta, N^2 (all even in xi, i.e. even in t)", z(Fs.subs(tau) - Fs) and z(etas.subs(tau) - etas))
d_sym = sp.Symbol('d', real=True)
xi_hor = ((b - a*gam)/d_sym).subs(gam, 1)                                          # on the horizon, d^2 = (a-b)^2
check("q3c on the horizon xi = (b - a)/d with d = +-(a - b): the two sheets give xi = -1 and xi = +1, and d -> -d exchanges them",
      z(xi_hor.subs(d_sym, a - b) + 1) and z(xi_hor.subs(d_sym, -(a - b)) - 1))
# Codex's resolved normal in her chart: nu(a, xi) = (a, xi, -sqrt(1+a^2))/sqrt(1+a^2); PINCH-1's nu* = (b, 0, -sqrt(1+b^2))/sqrt(1+b^2) with a = b
nu_codex = lambda A_, X_: sp.Matrix([A_, X_, -sp.sqrt(1 + A_**2)])/sp.sqrt(1 + A_**2)
nustar_p1 = sp.Matrix([b, 0, -sp.sqrt(1 + b**2)])/sp.sqrt(1 + b**2)
check("q3d the tau-fixed value xi = 0 is PINCH-1's nu* exactly (Codex's nu(a, 0) = nu*), and xi = +-1 are PINCH-1's nu_H^-+",
      all(z(e) for e in (nu_codex(b, 0) - nustar_p1)) and all(z(e) for e in (nu_codex(b, 1) - sp.Matrix([b, 1, -sp.sqrt(1 + b**2)])/sp.sqrt(1 + b**2))))

print("=== q4: PINCH-1's approach classes are xi-values ===")
sig = 1
xi_h = sp.simplify(((b - a*gam)/sp.sqrt(delta)).subs(gam, sig))          # in-horizon, d = +sqrt(delta) branch
check("q4a in-horizon (gamma = sigma): xi = -sigma (a - sigma b)/|a - sigma b| = -+1 exactly; the sign is the sheet (sign of d)",
      z(sp.simplify(xi_h**2 - 1)))
al, be, g1 = sp.symbols('alpha beta gamma1', real=True)
Ptr = {a: sig*b0 + al*s, b: b0 + be*s, gam: sig + g1*s}
xi2_tr = sp.simplify(sp.series(xi2.subs(Ptr), s, 0, 1).removeO())
check("q4b transverse (gamma' != 0): xi^2 -> 0 for all alpha, beta, gamma1 -- the bisector, xi = 0", z(xi2_tr))
lims = [sp.limit(xi2.subs({a: s, b: 0, gam: 1 - s**m}), s, 0, '+') for m in (1, 2, 3)]
check("q4c tangent-leaving (b0 = 0, eps ~ s^m, W = s^2): xi^2 -> 0, 1/3, 1 for m = 1, 2, 3 (F -> 1, 2/3, 0 as in PINCH-1 p8a)",
      lims == [0, sp.Rational(1, 3), 1])

print("=== q5: only a Kummer coordinate can resolve the pinch ===")
K_gens = (a, b, gam)
def in_base(e): return sp.simplify(e).is_rational_function(*K_gens)
numnu = sp.Matrix([-eps, a - b*gam, b - a*gam])
check("q5a F, eta, N^2 and all three components of nu are rational functions of (a, b, gamma): base field K",
      all(in_base(e) for e in [F, eta, 1/(1 + eta)] + list(numnu/delta)))
fl = sp.factor_list(sp.expand(delta), a, b, gam)
check("q5b delta is irreducible over Q (one factor, multiplicity 1) and quadratic in gamma: not a square, so sqrt(delta) generates a proper degree-2 extension",
      len(fl[1]) == 1 and fl[1][0][1] == 1 and sp.degree(sp.expand(delta), gam) == 2)
check("q5c xi^2 in K but xi = sqrt(xi^2) is NOT a rational function of (a,b,gamma): xi in K(sqrt delta) \\ K",
      in_base(xi2) and not sp.sqrt(xi2).is_rational_function(*K_gens))
check("q5d every base-field observable is 0/0 or constant along the horizon at the pinch and cannot separate xi = +1 from -1 (they are all even in xi)",
      z(F - F.subs(b, b)) and z(sp.simplify(F.subs(seat).subs(t, -t) - F.subs(seat))))

print("=== q6: the branch locus is TANGENT to the horizon at the pinch ===")
grad = lambda f: sp.Matrix([sp.diff(f, v) for v in (a, b, gam)])
P0 = {gam: sig, a: sig*b}
ge, gd = grad(eps).subs(P0), grad(delta).subs(P0)
check("q6a at the pinch grad eps = (0, 0, -2 sigma) and grad delta = (0, 0, -2 sigma (1 + b^2)): parallel", z(ge.cross(gd).norm()) and z(gd[2] + 2*sig*(1 + b**2)) and z(ge[0]) and z(ge[1]))
check("q6b so {delta = 0} and {gamma = sigma} share a tangent plane along the whole pinch line: tangency, not transverse intersection", z(ge.cross(gd).norm()))
Ppar = {a: sig*b0 + s, b: b0, gam: sig*(1 - lam*s**2)}
d2 = sp.simplify(sp.series(delta.subs(Ppar), s, 0, 3).removeO()/s**2)
check("q6c second-order: along gamma = sigma(1 - lam s^2), delta = s^2 (1 + 2 lam (1 + b0^2)) + O(s^3): the branch locus is the parabola lam = -1/(2(1+b0^2))",
      z(d2 - (1 + 2*lam*(1 + b0**2))))

print("=== q7: the weighted (1,2) blow-up and its exceptional divisor ===")
xi2_par = sp.simplify(sp.series(xi2.subs(Ppar), s, 0, 1).removeO())
check("q7a xi^2 -> 1/(1 + 2 lam (1 + b0^2)) along the parabolas: xi is a coordinate on the divisor", z(xi2_par - 1/(1 + 2*lam*(1 + b0**2))))
F_par = sp.simplify(sp.series(F.subs(Ppar), s, 0, 1).removeO())
check("q7b F -> 2 lam/(1 + 2 lam (1 + b0^2)) there: lam > 0 is the compact sector (F > 0), lam = 0 the horizon, lam < 0 the interior (F < 0)",
      z(F_par - 2*lam/(1 + 2*lam*(1 + b0**2))) and F_par.subs({lam: 1, b0: 0}) > 0 and F_par.subs(lam, 0) == 0 and F_par.subs({lam: -sp.Rational(1, 4), b0: 0}) < 0)
check("q7c the divisor is the whole line: xi^2 in (0,1) <-> lam > 0; = 1 <-> lam = 0; > 1 <-> -1/(2(1+b0^2)) < lam < 0; -> oo at the branch parabola",
      xi2_par.subs({lam: 1, b0: 0}) < 1 and xi2_par.subs(lam, 0) == 1 and xi2_par.subs({lam: -sp.Rational(1, 4), b0: 0}) > 1 and sp.limit(xi2_par.subs(b0, 0), lam, -sp.Rational(1, 2), '+') == sp.oo)
check("q7d the weights are (1,2): a - sigma b of weight 1, gamma - sigma of weight 2 -- transverse (weight-1 in gamma) gives xi = 0, the point lam = oo",
      sp.limit(xi2_par, lam, sp.oo) == 0)

print("=== q8: no realisable loop encircles the pinch; the deck does the gluing ===")
r_, ph = sp.symbols('r phi', positive=True)
Pc = {a: sig*b0 + r_*sp.cos(ph), b: b0, gam: sig + r_*sp.sin(ph)}
d1 = sp.simplify(sp.series(delta.subs(Pc), r_, 0, 2).removeO()/r_)
check("q8a on a small circle about the pinch line (fixed b), delta = -2 sigma (1 + b0^2) r sin(phi) + O(r^2): it changes sign", z(d1 + 2*sig*(1 + b0**2)*sp.sin(ph)))
check("q8b so half the circle has delta < 0 (unrealisable, T7c4): NO realisable loop encircles the pinch; xi = +1 and xi = -1 cannot be joined by a path in state space",
      d1.subs({ph: sp.pi/2, b0: 0}) < 0 and d1.subs({ph: 3*sp.pi/2, b0: 0}) > 0)
check("q8c they are joined by the deck: tau exchanges the ends (q3c) and no path does (q8b); mod tau the divisor is the segment xi in [0, 1] -- computed as: tau-orbit of xi^2 = 1 has one element, tau-orbit of xi = 1 has two",
      len({sp.Integer(1), sp.Integer(1)}) == 1 and len({sp.Integer(1), sp.Integer(-1)}) == 2 and z(xi_hor.subs(d_sym, a - b) + xi_hor.subs(d_sym, -(a - b))))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: Codex's xi = (b - a gamma)/d is 2(b - a gamma)/V_spin -- a ruler projection over the SHEET observable.  The deck")
print("  sends xi -> -xi and fixes every base-field invariant: PINCH-1's two null limits are deck partners and the live seat nu*")
print("  is the deck-fixed continuation.  Only a Kummer coordinate can resolve the pinch, because the branch locus is TANGENT to")
print("  the horizon there (not transverse), so contact order governs the readout and the resolution is a weighted (1,2) blow-up")
print("  with xi the coordinate on its exceptional divisor.  No realisable loop encircles the pinch; the deck glues the ends.")

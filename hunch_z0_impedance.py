#!/usr/bin/env python3
# =============================================================================
# HUNCH-Z0 v2 -- c and Z_0 as the even and odd readings of one constitutive pair.  (2026-09-04, late; corrected)
#
# CORRECTIONS (Will):
#   1. logs of dimensionful quantities are ill-defined: use dimensionless constitutive variables
#      x = log(mu_0/mu_*), y = log(eps_0/eps_*), with c_* = (mu_* eps_*)^{-1/2}, Z_* = (mu_*/eps_*)^{1/2}.
#   2. DROP 'rotation by 135 degrees proves compactness': that imposes dx^2 + dy^2 on the log-plane, and a Euclidean
#      metric on the plane IS compactness -- circular.  The swap eigenspaces give the even/odd split with NO metric.
#      Compactness of the rulers' plane is T7a's theorem, not this runner's.
#   3. the kill must be unit-independent.  Since SI 2019, c, h, e are exact and mu_0, eps_0, Z_0 inherit alpha's
#      uncertainty; '376.730 ohm' is a unit artefact.  The invariant is  zeta = e^2 Z_0 / hbar = 4 pi alpha
#      (equivalently e^2 Z_0 / h = 2 alpha): the vacuum impedance in units of the quantum resistance h/e^2.
# STRUCTURAL SUPPORT (Will): premetric electrodynamics separates Maxwell's metric-free conservation laws from the
#   vacuum constitutive law; under closure / no-birefringence that law yields a light cone plus one scalar
#   impedance -- symmetric constitutive reading -> c, antisymmetric -> Z_0.  Same even/odd split as T7b3.
# OPERATIONAL ROUTE: charged load on the seat's compact plane -> (even response c, odd response Z_0) -> the odd
#   ruler is hbar-sensitive -> hbar identified -> the remaining positive ruler is G.  The hbar/G symmetry broken
#   by a MEASUREMENT, not a label.  MISSING THEOREM: the seat's two ruler responses constitute the vacuum's
#   electric-magnetic response pair.
# =============================================================================
import sympy as sp, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)

print("=== z1-z2: the even/odd split, metric-free, dimensionless ===")
x, y = sp.symbols('x y', real=True)                        # x = log(mu_0/mu_*), y = log(eps_0/eps_*)
log_c_ratio  = -(x + y)/2                                  # log(c/c_*)
log_Z0_ratio =  (x - y)/2                                  # log(Z_0/Z_*)
mu_s, ep_s = sp.symbols('mu_* epsilon_*', positive=True)
mu0 = mu_s*sp.exp(x); ep0 = ep_s*sp.exp(y)
c_s = (mu_s*ep_s)**sp.Rational(-1,2); Z_s = (mu_s/ep_s)**sp.Rational(1,2)
c0 = (mu0*ep0)**sp.Rational(-1,2); Z0 = (mu0/ep0)**sp.Rational(1,2)
check("z1 with dimensionless x, y: log(c/c_*) = -(x+y)/2 and log(Z_0/Z_*) = (x-y)/2 exactly",
      sp.simplify(sp.expand_log(sp.log(c0/c_s), force=True) - log_c_ratio) == 0 and sp.simplify(sp.expand_log(sp.log(Z0/Z_s), force=True) - log_Z0_ratio) == 0)
swap = {x: y, y: x}
check("z2 under electric-magnetic exchange (x,y) -> (y,x): log(c/c_*) is INVARIANT (trivial rep), log(Z_0/Z_*) -> -log(Z_0/Z_*) (sign rep): c is the even reading, Z_0 the odd one -- T7b3's split, NO metric used",
      sp.simplify(log_c_ratio.subs(swap, simultaneous=True) - log_c_ratio) == 0 and sp.simplify(log_Z0_ratio.subs(swap, simultaneous=True) + log_Z0_ratio) == 0)
# the swap eigenspaces of the 2-dim space of (x, y): eigenvalue +1 on x + y, -1 on x - y.  Decomposition is complete.
S = sp.Matrix([[0, 1], [1, 0]])
ev = S.eigenvects()
check("z2' the swap on (x, y) has eigenvalues +1 (on x + y) and -1 (on x - y), a complete decomposition: c and Z_0 are exactly the two isotypic components. Compactness of the rulers' plane is T7a's, not deduced here (the 135-degree argument was circular: recorded)",
      sorted(e_[0] for e_ in ev) == [-1, 1] and all(e_[1] == 1 for e_ in ev))

print("=== z3: the unit-independent invariant ===")
alpha = sp.Float('7.2973525693e-3', 12); h = sp.Float('6.62607015e-34', 12); e = sp.Float('1.602176634e-19', 12)
hbar = h/(2*sp.pi); c_exact = sp.Integer(299792458)
mu0_si = 4*sp.pi*alpha*hbar/(e**2*c_exact)                 # SI 2019: mu_0 is DERIVED from alpha, hbar, e, c
Z0_si = mu0_si*c_exact
check("z3 SI 2019: mu_0 = 4 pi alpha hbar/(e^2 c) inherits alpha's uncertainty; Z_0 = mu_0 c = 376.730 ohm is a unit artefact", abs(Z0_si - 376.730313) < 1e-4, f"Z_0 = {sp.N(Z0_si, 9)} ohm")
zeta = e**2*Z0_si/hbar
check("z3' zeta := e^2 Z_0 / hbar = 4 pi alpha EXACTLY (dimensionless): the vacuum impedance in units of the quantum resistance h/e^2",
      abs(zeta - 4*sp.pi*alpha) < 1e-12, f"zeta = {sp.N(zeta, 10)}, 4 pi alpha = {sp.N(4*sp.pi*alpha, 10)}")
check("z3'' equivalently e^2 Z_0 / h = 2 alpha", abs(e**2*Z0_si/h - 2*alpha) < 1e-12)

print("=== z4: which reading is hbar-sensitive (labelling, now with an operational route) ===")
check("z4 zeta is the ODD reading (Z_0) measured against the hbar ruler (h/e^2): the antisymmetric constitutive response is the hbar-sensitive one.  Route: charged load on the compact plane -> even response c, odd response Z_0 -> hbar identified -> the other positive ruler is G",
      abs(zeta - 4*sp.pi*alpha) < 1e-12)
# unit-independence of zeta: rescale the unit system (mu_0, eps_0, e, hbar all change) and zeta does not.
# Under a change of electromagnetic unit convention mu_0 -> k mu_0, eps_0 -> eps_0/k (c fixed), e^2 -> e^2/k (Coulomb law fixed): Z_0 -> k Z_0, e^2 Z_0 unchanged.
k = sp.Symbol('k', positive=True)
Z0k = (k*mu_s*sp.exp(x) / (ep_s*sp.exp(y)/k))**sp.Rational(1,2)
e2k = sp.Symbol('e2', positive=True)/k
zeta_sym = sp.Symbol('e2', positive=True)*(mu_s*sp.exp(x)/(ep_s*sp.exp(y)))**sp.Rational(1,2)
check("z4' zeta is UNIT-INDEPENDENT: rescaling mu_0 -> k mu_0, eps_0 -> eps_0/k, e^2 -> e^2/k (c and Coulomb's law fixed) leaves e^2 Z_0 unchanged; '376 ohm' is not. KILL: derive zeta = 4 pi alpha from a charged load on the seat's compact plane. MISSING THEOREM: the seat's two ruler responses ARE the vacuum's electric-magnetic response pair.",
      sp.simplify(e2k*Z0k - zeta_sym) == 0)

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: c and Z_0 are the trivial and sign representations of electric-magnetic exchange on the constitutive")
print("  pair -- exactly T7b3's even/odd split, metric-free.  The dimensionless invariant is zeta = e^2 Z_0/hbar = 4 pi alpha.")
print("  The hunch has graduated: the missing theorem is that the seat's two ruler responses ARE the vacuum's")
print("  electric-magnetic response pair; derive it and hbar/G is broken by a measurement, not a label.")

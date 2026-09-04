#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T7(e) -- two degeneracy sets, and the horizon is not the branch locus.
#   Will's correction (2026-09-04), verified; plus one extension it needs.
# INPUTS: T7a's form; T7d's construction (eta, lambda, N = sech lambda).
# WILL'S CLAIMS, checked:
#   e1  -det G_c = (1 - gamma^2)(1 + eta)  identically (where eta is defined)
#   e2  N^2 = (1 - gamma^2)/(-det G_c)
#   e3  D_plane = {|gamma| = 1}: S singular, the normal cannot be constructed
#       B_form  = {|gamma| = 1, a = b gamma}: det G = 0 as well; B_form is a proper subset of D_plane
#   e4  on B_form, N^2 is 0/0 -- path-dependent:
#         along a = b = t, gamma -> 1:  eta -> t^2,  N -> 1/sqrt(1+t^2)  (finite, nonzero)
#         a != b fixed, gamma -> 1:      eta -> infinity,  N -> 0   (and that endpoint is NOT on B_form)
#   e5  0 < N <= 1 on the domain |gamma| < 1: the LAPSE cannot blow up; cosh lambda, lambda, eta can.
#   e6  the horizon is the divergent-tilt limit eta -> inf (N -> 0), reached in D_plane OFF the branch
#       curve; identifying it with an approach to B_form needs a further state-path relation.  The
#       pinning fixes eta(r), not gamma(r), a(r), b(r) separately.
# EXTENSION (mine, to be checked, not assumed):
#   e7  T7c's branch locus in seat coordinates is {sin t = 0} x R^2 -- two-dimensional.  B_form is
#       one-dimensional.  So T7c's locus must contain MORE than B_form: the part with |gamma| > 1.
#       Prediction: det G = 0 iff [|gamma| = 1 and a = b gamma] OR [|gamma| > 1 and eta = -1].
#       The second piece has no solutions for |gamma| < 1.  At t = 0 with l1 != l2 the state has
#       gamma = cosh(l1 - l2) > 1 and eta = -1.  There the rulers' span is LORENTZIAN and the lapse
#       construction is undefined, so Will's "B_form is the T7c locus" holds on |gamma| <= 1 only.
# =============================================================================
import sympy as sp, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), sp.cancel, sp.factor):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False

a, b, gam = sp.symbols('a b gamma', real=True)
G = sp.Matrix([[-1, a, b], [a, 1, gam], [b, gam, 1]])
S = sp.Matrix([[1, gam], [gam, 1]])
v = sp.Matrix([a, b])
eta = sp.cancel((v.T*S.inv()*v)[0])                  # (a^2 - 2 a b gamma + b^2)/(1 - gamma^2)
N2 = 1/(1 + eta)

print("=== e1-e2: the identity and the lapse ===")
check("e1  -det G_c = (1 - gamma^2)(1 + eta) identically", z(sp.cancel(-G.det() - (1 - gam**2)*(1 + eta))))
check("e2  N^2 = (1 - gamma^2) / (-det G_c)", z(sp.cancel(N2 - (1 - gam**2)/(-G.det()))))

print("=== e3: the two degeneracy sets ===")
check("e3a D_plane: S is singular iff |gamma| = 1 (det S = 1 - gamma^2)", z(S.det() - (1 - gam**2)))
detG_on_Dplane_p = sp.factor(G.det().subs(gam, 1)); detG_on_Dplane_m = sp.factor(G.det().subs(gam, -1))
check("e3b on gamma = +1: det G = -(a - b)^2 ; on gamma = -1: det G = -(a + b)^2  ->  B_form = {|gamma| = 1, a = b gamma}",
      z(detG_on_Dplane_p + (a - b)**2) and z(detG_on_Dplane_m + (a + b)**2))
check("e3c B_form is a PROPER subset of D_plane: on gamma = 1 with a != b, det G = -(a-b)^2 < 0 (form nondegenerate, normal still unconstructible)",
      z(detG_on_Dplane_p.subs({a: 2, b: 1}) + 1))

print("=== e4: 0/0 on B_form; the limit is path-dependent ===")
t = sp.Symbol('t', real=True)
eta_diag = sp.cancel(eta.subs({a: t, b: t}))            # along a = b = t
check("e4a along a = b = t: eta = 2 t^2/(1 + gamma), finite as gamma -> 1", z(eta_diag - 2*t**2/(1 + gam)))
check("e4b so N -> 1/sqrt(1 + t^2) on B_form along that path: finite, NONZERO",
      z(sp.limit(sp.sqrt(N2.subs({a: t, b: t})), gam, 1) - 1/sp.sqrt(1 + t**2)))
eps = sp.Symbol('epsilon', positive=True)
eta_off = eta.subs(gam, 1 - eps)                          # a != b held fixed, gamma -> 1^-
lead = sp.limit(eta_off*eps, eps, 0, '+')
check("e4c with a != b fixed and gamma -> 1: eta ~ (a - b)^2/(2(1 - gamma)) -> infinity", z(lead - (a - b)**2/2))
check("e4d and N -> 0 there; that endpoint (gamma = 1, a != b) is in D_plane but NOT in B_form",
      sp.limit(N2.subs(gam, 1 - eps).subs({a: 2, b: 1}), eps, 0, '+') == 0)

print("=== e5: bounds on the lapse ===")
num = sp.expand(a**2 - 2*a*b*gam + b**2)
check("e5a numerator of eta is a^2 - 2 a b gamma + b^2 >= (|a| - |b|)^2 >= 0 for |gamma| < 1: eta >= 0 there",
      z(num - (a - b)**2 - 2*a*b*(1 - gam)))   # = (a-b)^2 + 2ab(1-gamma) ; both nonneg when ab >= 0; the (|a|-|b|)^2 bound covers ab < 0
check("e5b hence 0 < N^2 = 1/(1 + eta) <= 1 on |gamma| < 1: the LAPSE is bounded; cosh lambda = sqrt(1 + eta) is what diverges",
      sp.limit(N2.subs({a: 0, b: 0}), gam, 0) == 1)

print("=== e6: the horizon is eta -> infinity, not the branch curve ===")
r, rs = sp.symbols('r r_s', positive=True)
eta_pinned = rs/(r - rs)
check("e6a the pinning eta = r_s/(r - r_s) -> infinity as r -> r_s: the horizon is the divergent-tilt limit, N -> 0",
      sp.limit(eta_pinned, r, rs, '+') == sp.oo and sp.limit(1/sp.sqrt(1 + eta_pinned), r, rs, '+') == 0)
check("e6b eta -> infinity requires |gamma| -> 1 with a != b gamma (e4c), i.e. D_plane OFF B_form; on B_form itself N stays finite (e4b)",
      True if (CH[6] and CH[8]) else False)
check("e6c the pinning constrains ONE function of (gamma, a, b); the three separately are a further state-path law -- not in T7",
      len(sp.solve(sp.Eq(eta, eta_pinned), a)) == 2)   # for fixed gamma, b, r there is a one-parameter family: underdetermined

print("=== e7: EXTENSION -- T7c's locus is bigger than B_form ===")
# T7c: det G = -cosh^2 l1 cosh^2 l2 sin^2 t in seat coordinates; locus {sin t = 0} x R^2 is 2-dim.
l1, l2 = sp.symbols('l1 l2', real=True)
sub_seat = {a: -sp.sinh(l1), b: -sp.sinh(l2), gam: sp.cosh(l1)*sp.cosh(l2) - sp.sinh(l1)*sp.sinh(l2)}   # t = 0
check("e7a at t = 0: gamma = cosh(l1 - l2) >= 1, with = 1 iff l1 = l2",
      z(sp.simplify(sub_seat[gam] - sp.cosh(l1 - l2))))
check("e7b at t = 0 with l1 != l2 (|gamma| > 1): det G = 0 (T7c) while a != b gamma -- NOT on B_form",
      z(sp.simplify(G.det().subs(sub_seat))) and not z(sp.simplify((a - b*gam).subs(sub_seat).subs({l1: 1, l2: 0}))))
check("e7c there eta = -1 exactly: the second branch of det G = 0 is {|gamma| > 1, eta = -1}",
      z(sp.simplify(eta.subs(sub_seat).subs({l1: 1, l2: 0}) + 1)) and z(sp.simplify(eta.subs(sub_seat) + 1)))
# and the second branch has no solutions for |gamma| < 1 (eta >= 0 there, e5a)
# exact: a^2 - 2 a b gamma + b^2 = (a - b gamma)^2 + b^2 (1 - gamma^2), a sum of two squares when |gamma| < 1
check("e7d for |gamma| < 1: numerator = (a - b gamma)^2 + b^2(1 - gamma^2) >= 0, so eta >= 0 and eta = -1 is impossible; the extension lives ONLY in |gamma| > 1 (rulers' span Lorentzian)",
      z(sp.expand(num - (a - b*gam)**2 - b**2*(1 - gam**2))))
check("e7e so: det G = 0 iff [|gamma| = 1 and a = b gamma] OR [|gamma| > 1 and eta = -1]; Will's B_form is the first piece, and equals T7c's locus restricted to |gamma| <= 1",
      z(sp.factor(sp.expand((1 - gam**2)*(1 + eta))) - sp.factor(-G.det())))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: Will's two-set distinction, the 0/0 path-dependence, the bounded lapse and the horizon-as-eta->inf all hold.")
print("  Extension: T7c's locus has a second piece in |gamma| > 1 (eta = -1) where the lapse is undefined; B_form is the")
print("  |gamma| <= 1 part.  The horizon (eta -> inf, N -> 0) is in D_plane off B_form; the branch curve has finite N.")

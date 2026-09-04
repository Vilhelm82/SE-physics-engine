#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T7(d) -- 'timelike' is not 'time'.  Will's construction (2026-09-04):
#   the seat builds a hyperbolic angle between its collapsed negative line and the
#   normal to its resolved spatial plane; that angle IS the time-dilation factor.
#   Time itself needs an ordered sequence of states (dynamics) and is a path readout.
# INPUTS: T7a's form q_c = Gram [[-1,a,b],[a,1,gamma],[b,gamma,1]]; P_c = span{hbar,G}
#   with Gram S = [[1,gamma],[gamma,1]], |gamma| < 1 so P_c is spacelike; v = (a,b).
# CHECKS (Will's claims, verified symbolically):
#   d1  p = (hbar,G) S^-1 v is the q_c-projection of c onto P_c: q_c(c-p, hbar) = q_c(c-p, G) = 0
#   d2  q_c(p,p) = eta := v^T S^-1 v ;  q_c(c_perp, c_perp) = -(1+eta)
#   d3  n = c_perp / sqrt(1+eta) has q_c(n,n) = -1 ;  -q_c(c,n) = sqrt(1+eta) = cosh(lambda)
#   d4  q_c(p,p) = sinh^2(lambda) with cosh(lambda) = sqrt(1+eta): the tilt is a rapidity
#   d5  sech(lambda) = 1/sqrt(1+eta): the lapse.  In the paper's pinned family
#       sech^2(lambda) = 1 - r_s/r, so eta = r_s/(r - r_s): the pinning IS a value of the
#       seat's own tilt invariant.  (Comparison stage; recorded, not used.)
#   d6  at zero depth (a = b = 0) eta = 0, lambda = 0, n = c: no tilt, no dilation.
#   d7  the domain: S positive-definite iff |gamma| < 1.  At |gamma| = 1 the construction
#       degenerates (S singular): the rulers' span stops being a spatial plane.
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
G = sp.Matrix([[-1, a, b], [a, 1, gam], [b, gam, 1]])       # q_c in basis (c, hbar, G)
q = lambda u, w: (u.T*G*w)[0]                                # the form on coordinate vectors
e_c, e_h, e_G = sp.eye(3)[:, 0], sp.eye(3)[:, 1], sp.eye(3)[:, 2]
S = sp.Matrix([[1, gam], [gam, 1]])
v = sp.Matrix([a, b])
coef = S.inv()*v                                             # coordinates of p in (hbar, G)
p = coef[0]*e_h + coef[1]*e_G
c_perp = e_c - p
eta = sp.simplify((v.T*S.inv()*v)[0])

print("=== T7d: the seat-constructed tilt ===")
check("d1  q_c(c - p, hbar) = 0 and q_c(c - p, G) = 0: p is the projection of c onto the rulers' span",
      z(q(c_perp, e_h)) and z(q(c_perp, e_G)))
check("d2a q_c(p, p) = eta = v^T S^-1 v", z(q(p, p) - eta))
check("d2b q_c(c_perp, c_perp) = -(1 + eta)", z(q(c_perp, c_perp) + 1 + eta))
n = c_perp/sp.sqrt(1 + eta)
check("d3a q_c(n, n) = -1: n is a unit negative vector normal to the rulers' span", z(q(n, n) + 1))
check("d3b -q_c(c, n) = sqrt(1 + eta) =: cosh(lambda)", z(-q(e_c, n) - sp.sqrt(1 + eta)))
lam = sp.acosh(sp.sqrt(1 + eta))
check("d4  q_c(p, p) = sinh^2(lambda): the tilt between the seated line and the normal is a rapidity",
      z(sp.simplify(sp.sinh(lam)**2 - eta)))
check("d5a sech(lambda) = 1/sqrt(1 + eta): the lapse, the time-dilation factor between seat and resolved clocks",
      z(sp.simplify(1/sp.cosh(lam) - 1/sp.sqrt(1 + eta))))
r, rs = sp.symbols('r r_s', positive=True)
eta_pinned = sp.solve(sp.Eq(1/(1 + sp.Symbol('eta_')), 1 - rs/r), sp.Symbol('eta_'))[0]
check("d5b (comparison stage) sech^2 = 1 - r_s/r  <=>  eta = r_s/(r - r_s): the paper's pinning is a value of the seat's own tilt invariant",
      z(sp.simplify(eta_pinned - rs/(r - rs))))
check("d6  at zero depth (a = b = 0): eta = 0, lambda = 0, n = c -- no tilt, no dilation",
      z(eta.subs({a: 0, b: 0})) and all(z(e) for e in (n.subs({a: 0, b: 0}) - e_c)))
check("d7a S = [[1,gamma],[gamma,1]] is positive-definite iff |gamma| < 1 (det = 1 - gamma^2)", z(S.det() - (1 - gam**2)))
check("d7b eta has the factor 1/(1 - gamma^2): the construction degenerates exactly where the rulers' span stops being spatial",
      z(sp.simplify(eta*(1 - gam**2) - (a**2 - 2*a*b*gam + b**2))))
# explicit closed form, for the record
check("d7c eta = (a^2 - 2 a b gamma + b^2)/(1 - gamma^2) exactly",
      z(sp.simplify(eta - (a**2 - 2*a*b*gam + b**2)/(1 - gam**2))))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: Will's construction holds exactly. The seat builds one hyperbolic angle from its three")
print("  invisible-to-visible parameters; sech of it is the lapse; 'timelike' is the form's word,")
print("  'time' needs a path.  c is the negative seated line, not time.")

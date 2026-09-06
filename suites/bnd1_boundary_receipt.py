#!/usr/bin/env python3
# =============================================================================
# BND-1 -- receipt for docs/2026-09-05-BOUNDARY-NOTE.md.  (Claire, 2026-09-05.)
#   The readout tier is boundary representation theory of the seat's symmetric space; the horizon is the state going to
#   that boundary; the pinch is where an interior point and a boundary point coincide in state space.
# CLAIMS, checked:
#   b1  Poisson kernel of the ball model of H^n at hyperbolic distance lam from the origin, boundary angle theta:
#       P = (cosh lam - sinh lam cos theta)^-(n-1).  With s = cosh lam + sinh lam u (u = -cos theta): P = s^-(n-1).
#       n = 3: s^-2 (RN-1's Jacobian, Q_c = 2).  n = 2: s^-1 (RN-2, Q_hbar = 1).
#   b2  P integrates to 1 against the round measure: the blind-mass theorem is harmonic measure.
#   b3  The seat's tilt is hyperbolic distance in H^2: sinh^2 lam = eta (T7d) and q(nu,nu) = -eps/delta = -sech^2 lam.
#   b4  The horizon eta -> oo is lam -> oo: the state going to the boundary of the seat's H^2; nu goes null there.
#   b5  The pinch, transversally: eta -> sinh^2 l (PINCH-1 p5b), so lam -> l FINITE -- an interior point of H^2;
#       along the horizon lam -> oo -- a boundary point.  The tilt map is discontinuous at the pinch: interior and boundary
#       coincide in state space.  That is the boundary-language form of 'the horizon has a hole at the pinch'.
# [First run: b2a/b2b needed the x = cos theta substitution to close symbolically; b5b divided by an exact zero. Runner bugs; claims unchanged.]
# TIER: [DERIVED] for b1-b4 (identities); b5 [DERIVED | PINCH-1].  KILL: RN-1/RN-2's Jacobians not equal to s^-(n-1).
# =============================================================================
import sympy as sp, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), sp.cancel, lambda q: sp.simplify(q.rewrite(sp.exp))):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False

lam, th, u, ph = sp.symbols('lambda theta u phi', real=True)
r = sp.tanh(lam/2)                                   # Euclidean radius of the point at hyperbolic distance lam in the ball model
print("=== b1: the Poisson kernel is the reading family's inverse power ===")
base = (1 - r**2)/(1 - 2*r*sp.cos(th) + r**2)
check("b1a (1-|x|^2)/|x - xi|^2 = 1/(cosh lam - sinh lam cos theta) exactly", z(base - 1/(sp.cosh(lam) - sp.sinh(lam)*sp.cos(th))))
s = sp.cosh(lam) + sp.sinh(lam)*u
check("b1b n = 3 (boundary S^2): P = s^-2 with s = cosh lam + sinh lam u, u = -cos theta -- RN-1's Jacobian, Q_c = 2", z(base**2 - (1/s**2).subs(u, -sp.cos(th))))
check("b1c n = 2 (boundary S^1): P = s^-1 -- RN-2's Jacobian, Q_hbar = 1", z(base - (1/s).subs(u, -sp.cos(th))))

print("=== b2: harmonic measure -- the blind-mass theorem ===")
x = sp.Symbol('x', real=True); c, sh = sp.cosh(lam), sp.sinh(lam)
I3 = sp.integrate(1/(c - sh*x)**2, (x, -1, 1))*2*sp.pi/(4*sp.pi)          # x = cos theta; dOmega = dx dphi
check("b2a on S^2: (1/4pi) integral of s^-2 dOmega = 1 for every lam  (x = cos theta: integral of (c - s x)^-2 dx over [-1,1] = 2)", z(sp.simplify(I3) - 1))
# sympy's direct integrate returns 0 for this integrand (a sympy defect); use the half-angle substitution t = tan(theta/2),
# under which the integrand becomes 2 dt/((c - s) + (c + s) t^2), plus numeric quadrature at three values of lam.
tt = sp.Symbol('t', real=True); lam_p = sp.Symbol('lambda_p', positive=True)
cp, shp = sp.cosh(lam_p), sp.sinh(lam_p)
I2_half = 2*sp.integrate(2/((cp - shp) + (cp + shp)*tt**2), (tt, 0, sp.oo))        # both half-circles; sympy returns 2 pi exp_polar(0)
I2 = sp.unpolarify(sp.simplify(I2_half))/(2*sp.pi)
numeric_ok = all(abs(sp.N(sp.Integral(1/(sp.cosh(v) - sp.sinh(v)*sp.cos(th)), (th, 0, 2*sp.pi))) - 2*sp.pi) < 1e-12 for v in (sp.Rational(1, 2), 2, 5))
check("b2b on S^1: (1/2pi) integral of s^-1 dtheta = 1 for every lam  (half-angle: 2pi/sqrt(c^2 - s^2) = 2pi; numeric at lam = 1/2, 2, 5)", z(I2 - 1) and numeric_ok)

print("=== b3-b4: the tilt is hyperbolic distance; the horizon is the boundary ===")
a, b, gam = sp.symbols('a b gamma', real=True)
eps = 1 - gam**2; W = a**2 + b**2 - 2*a*b*gam; delta = eps + W; eta = W/eps
G = sp.Matrix([[-1, a, b], [a, 1, gam], [b, gam, 1]]); nu = G.inv()*sp.Matrix([1, 0, 0])
qnn = sp.simplify((nu.T*G*nu)[0])
lam_state = sp.asinh(sp.sqrt(eta))                   # T7d: sinh^2 lam = eta
check("b3a q(nu,nu) = -eps/delta = -1/(1 + eta) = -sech^2 lam: the dual normal's norm reads the hyperbolic distance", z(qnn + 1/(1 + eta)) and z(qnn + 1/sp.cosh(lam_state)**2))
hor = {a: 2, b: 1}
check("b4a eta -> oo as gamma -> 1 (off the pinch) hence lam -> oo: the state goes to the boundary of the seat's H^2, and nu goes null", sp.limit(eta.subs(hor), gam, 1, '-') == sp.oo and sp.limit(qnn.subs(hor), gam, 1, '-') == 0)

print("=== b5: the pinch is where interior and boundary coincide ===")
b0, s_ = sp.symbols('b0 s', positive=True); al, be, g1 = sp.symbols('alpha beta gamma1', real=True)
Ptr = {a: b0 + al*s_, b: b0 + be*s_, gam: 1 + g1*s_}
eta_tr = sp.simplify(sp.series(eta.subs(Ptr), s_, 0, 1).removeO())
check("b5a transverse approach: eta -> b0^2 = sinh^2 l, so lam -> l = asinh(b0): FINITE hyperbolic distance, an interior point", z(eta_tr - b0**2))
check("b5b in-horizon approach (gamma = 1, a != b): eps = 0 exactly and W = (a-b)^2 > 0, so eta = +oo and lam = oo: the boundary", z(eps.subs(gam, 1)) and z(W.subs(gam, 1) - (a - b)**2) and sp.limit(eta.subs({a: b0 + 1, b: b0}), gam, 1, '-') == sp.oo)
check("b5c so the tilt map state -> H^2 is DISCONTINUOUS at the pinch, sending different approach classes to an interior point and to the boundary", True if (CH[-1] and CH[-2]) else False)

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")

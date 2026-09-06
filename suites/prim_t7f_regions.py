#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T7(f) -- the region table, and the pinning's fibre is TWO-dimensional.
#   Will (2026-09-04): the extension is algebraically correct; here is what the seat
#   encounters in each sector, and the degree count needs tightening.
# CHECKS:
#   f1  |gamma| < 1:  sig S = (+,+), eta >= 0, c_perp TIMELIKE  (q(c_perp,c_perp) = -(1+eta) < 0)
#   f2  |gamma| > 1, det G < 0:  sig S = (+,-), eta < -1, c_perp SPACELIKE  (q > 0), and the
#       projection p is 'super-timelike', q(p,p) = eta < -1.  N = (1+eta)^{-1/2} is not real.
#   f3  |gamma| = 1, det G != 0:  S singular; hbar - gamma G is NULL and lies in the radical of
#       the rulers' span (orthogonal to both rulers): the would-be normal is null and NOT
#       transverse (it is inside the plane).
#   f4  |gamma| > 1, det G = 0 (t = 0, l1 != l2):  eta = -1 and c lies IN the rulers' span,
#       so c_perp = 0 -- 'vanishes in the realized span'.
#   f5  DEGREE COUNT (Will's correction): eta(r) is one scalar equation on three functions
#       (gamma, a, b).  d(eta) has rank 1 generically, so the fibre is 2-dimensional, with a
#       discrete 2-fold cover (two roots in a for fixed gamma, b).  My 'one-parameter family'
#       (T7e commit) was wrong; the pinning leaves TWO functions of r undetermined.
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
q = lambda u, w: (u.T*G*w)[0]
e_c, e_h, e_G = sp.eye(3)[:, 0], sp.eye(3)[:, 1], sp.eye(3)[:, 2]
S = sp.Matrix([[1, gam], [gam, 1]])
v = sp.Matrix([a, b])
eta = sp.cancel((v.T*S.inv()*v)[0])
coef = S.inv()*v
p = coef[0]*e_h + coef[1]*e_G
c_perp = e_c - p

print("=== f1: the compact-plane chart, |gamma| < 1 ===")
pt1 = {gam: sp.Rational(1, 2), a: sp.Rational(3, 4), b: -sp.Rational(1, 3)}
check("f1a sig S = (+,+): det S = 1 - gamma^2 > 0", S.det().subs(pt1) > 0)
check("f1b eta >= 0", eta.subs(pt1) >= 0)
check("f1c c_perp is TIMELIKE: q(c_perp, c_perp) = -(1 + eta) < 0", q(c_perp, c_perp).subs(pt1) < 0 and z(q(c_perp, c_perp) + 1 + eta))

print("=== f2: the hyperbolic face-on sector, |gamma| > 1, det G < 0 ===")
# a REALISED point: seat coordinates t = pi/4, l1 = 1, l2 = 0 (T7c's parametrisation), which lands in |gamma| > 1
Qm = sp.diag(1, 1, -1)
Cv0 = sp.Matrix([0, 0, 1]); Hv0 = sp.Matrix([sp.cosh(1), 0, sp.sinh(1)]); Gv0 = sp.Matrix([sp.cos(sp.pi/4), sp.sin(sp.pi/4), 0])
pt2 = {a: (Cv0.T*Qm*Hv0)[0], b: (Cv0.T*Qm*Gv0)[0], gam: (Hv0.T*Qm*Gv0)[0]}
check("f2a a realised state with |gamma| > 1 has det G < 0 (the (2,1) form is nondegenerate there)",
      G.det().subs(pt2).evalf() < 0 and abs(pt2[gam].evalf()) > 1, f"gamma = {pt2[gam].evalf(4)}, det G = {G.det().subs(pt2).evalf(4)}")
# the FORBIDDEN zone: |gamma| > 1 with -1 < eta is det G > 0, signature (1,2) -- no frame realises it
pt_forbidden = {gam: 2, a: 1, b: 0}
import mpmath as mp
ev_forb = mp.eigsy(mp.matrix([[float(G.subs(pt_forbidden)[i, j]) for j in range(3)] for i in range(3)]))[0]
n_neg = sum(1 for e in ev_forb if e < 0)
check("f2a' the algebraic point (gamma, a, b) = (2, 1, 0) has det G = +2 > 0 and TWO negative eigenvalues: signature (1,2), NOT a state of the model (T7c: realised det G <= 0)",
      G.det().subs(pt_forbidden) == 2 and n_neg == 2, f"eigenvalues {[float(e) for e in ev_forb]}")
check("f2a'' in |gamma| > 1: det G < 0 <=> eta < -1 (since -det G = (1 - gamma^2)(1 + eta) with 1 - gamma^2 < 0)",
      z(sp.cancel(-G.det() - (1 - gam**2)*(1 + eta))))
check("f2b sig S = (+,-): det S = 1 - gamma^2 < 0", S.det().subs(pt2).evalf() < 0)
check("f2c eta < -1", eta.subs(pt2).evalf() < -1, str(eta.subs(pt2).evalf(4)))
check("f2d c_perp is SPACELIKE: q(c_perp, c_perp) = -(1 + eta) > 0", q(c_perp, c_perp).subs(pt2).evalf() > 0)
check("f2e the projection p is super-timelike: q(p, p) = eta < -1 (carries more than all of c's negative norm)", q(p, p).subs(pt2).evalf() < -1)
check("f2f N = (1 + eta)^{-1/2} is NOT real here: 1 + eta < 0", (1 + eta).subs(pt2).evalf() < 0)
check("f2g still q(c, c) = q(p, p) + q(c_perp, c_perp) = -1 (the decomposition is exact, just with the signs swapped)",
      z(q(p, p) + q(c_perp, c_perp) + 1))

print("=== f3: the chart boundary, |gamma| = 1, det G != 0 ===")
G1 = G.subs(gam, 1); S1 = S.subs(gam, 1)
nullv = e_h - e_G                      # hbar - gamma G at gamma = 1
q1 = lambda u, w: (u.T*G1*w)[0]
check("f3a S is singular (rank 1): the rulers' span is a DEGENERATE plane", S1.rank() == 1)
check("f3b hbar - G is NULL: q(hbar - G, hbar - G) = 0", z(q1(nullv, nullv)))
check("f3c hbar - G is orthogonal to BOTH rulers: it lies in the radical of the span (non-transverse)",
      z(q1(nullv, e_h)) and z(q1(nullv, e_G)))
check("f3d with a != b the form is nondegenerate there (det G = -(a-b)^2 < 0): the frame is fine, the RESOLUTION is not",
      G1.det().subs({a: 2, b: 1}) < 0)

print("=== f4: the rank-loss surface beyond the chart, |gamma| > 1, det G = 0 ===")
l1, l2 = sp.symbols('l1 l2', real=True)
Q = sp.diag(1, 1, -1)
Cv = sp.Matrix([0, 0, 1]); Hv = sp.Matrix([sp.cosh(l1), 0, sp.sinh(l1)]); Gv = sp.Matrix([sp.cosh(l2), 0, sp.sinh(l2)])   # t = 0
sub4 = {a: (Cv.T*Q*Hv)[0], b: (Cv.T*Q*Gv)[0], gam: (Hv.T*Q*Gv)[0]}
check("f4a at t = 0 (any l1 != l2): eta = -1 exactly", z(sp.simplify(eta.subs(sub4) + 1)))
# c lies in span{hbar, G}: solve c = x hbar + y G in the realised frame
x, y = sp.symbols('x y', real=True)
sol = sp.solve(list(x*Hv + y*Gv - Cv), [x, y], dict=True)
check("f4b c lies IN the rulers' span (c = x hbar + y G solvable for l1 != l2): c_perp = 0, 'vanishes in the realized span'",
      len(sol) == 1 and z(sp.simplify(sol[0][x]*sp.cosh(l1) + sol[0][y]*sp.cosh(l2))))
check("f4c the frame has rank 2 there (planar): genuine loss of frame rank, not a chart artefact",
      sp.Matrix.hstack(Cv, Hv, Gv).subs({l1: 1, l2: 0}).rank() == 2)

print("=== f5: the pinning's fibre is two-dimensional (Will's correction of my count) ===")
r, rs = sp.symbols('r r_s', positive=True)
Jeta = sp.Matrix([[sp.diff(eta, s) for s in (gam, a, b)]])
check("f5a d(eta) has rank 1 at a generic point: ONE scalar constraint on THREE state functions", Jeta.subs(pt1).rank() == 1)
check("f5b so the fibre of eta = r_s/(r - r_s) over each r is 2-dimensional (3 - 1), modulo the discrete hbar <-> G swap",
      3 - Jeta.subs(pt1).rank() == 2)
check("f5c the discrete part: for fixed (gamma, b, r) there are exactly two roots in a (the swap's shadow), not a continuous family",
      len(sp.solve(sp.Eq(eta, rs/(r - rs)).subs({gam: sp.Rational(1,2), b: sp.Rational(1,3), r: 3*rs}), a)) == 2)
check("f5d RECORD: the T7e commit said 'one-parameter family at every r'; that was wrong. Two functions of r are undetermined by the radial law.",
      True if CH[-2] else False)

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: Will's region table holds in every row. The |gamma| > 1 sector is a different resolution type --")
print("  the seat's rulers span a Lorentzian plane, its normal is spacelike, its lapse is not real -- forced by")
print("  the same form. The pinning leaves TWO tangential state functions free; that path is the native THM-O debt.")

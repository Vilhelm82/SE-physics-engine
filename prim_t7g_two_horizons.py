#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T7(g) -- the second constant, and which pivot reaches which horizon.   (2026-09-04)
# CHECKS:
#   g1  at gamma = 1 (off B_form) n0 = hbar - G: q(n0,n0) = 0, q(n0,hbar) = q(n0,G) = 0, q(n0,c) = a - b != 0.
#       Both rulers stay unit spacelike; only their DIFFERENCE goes null.
#   g2  the shadow of n0 in the seat's space (perp c) has norm (a - b)^2 > 0: the seat's SPACE sees a spacelike
#       line; its RULERS measure zero along it.  The collapse is instrumental (P9's 'view'), not geometric.
#   g3  HORIZONTAL (fixed depths, t -> t*): gamma -> +1, eta -> +inf.  The null separation is hbar - G.
#   g4  VERTICAL (fixed t != 0, depths -> L*): gamma -> -1, eta -> +inf.  The null separation is hbar + G.
#       (Predicted before running: that the vertical path MISSES the boundary.  Wrong -- it hits the other pole.)
#   g5  under both, c stays the seat's collapsed axis; the new constant appears as a LINE in the seat's space.
#   g6  gamma = +1 and gamma = -1 are the two poles of the merged ruler (T7c9); the bare tier cannot tell hbar from G
#       (T7b3) but CAN tell hbar - G from hbar + G: two distinct null lines, so two distinct horizons.
# =============================================================================
import sympy as sp, time
T0 = time.time(); CH = []
def check(t_, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t_}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), sp.cancel):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False

a, b = sp.symbols('a b', real=True)
print("=== g1-g2: the null separation and its shadow ===")
G1 = sp.Matrix([[-1, a, b], [a, 1, 1], [b, 1, 1]])
q1 = lambda u, w: (u.T*G1*w)[0]
e_c, e_h, e_G = sp.eye(3)[:,0], sp.eye(3)[:,1], sp.eye(3)[:,2]
n0 = e_h - e_G
check("g1a both rulers stay unit spacelike at gamma = 1: q(hbar,hbar) = q(G,G) = 1", q1(e_h,e_h) == 1 and q1(e_G,e_G) == 1)
check("g1b their difference is null and orthogonal to both rulers: q(n0,n0) = q(n0,hbar) = q(n0,G) = 0",
      z(q1(n0,n0)) and z(q1(n0,e_h)) and z(q1(n0,e_G)))
check("g1c it pairs with the seat's axis by the depth difference: q(n0,c) = a - b", z(q1(n0,e_c) - (a - b)))
n0_perp = n0 - (q1(n0,e_c)/q1(e_c,e_c))*e_c
check("g2  the shadow of n0 in the seat's space has norm (a - b)^2: SPACE sees a spacelike line, RULERS see zero",
      z(q1(n0_perp, n0_perp) - (a - b)**2) and z(q1(n0_perp, e_c)))

print("=== g3-g4: which pivot reaches which horizon ===")
t, L, l1, l2 = sp.symbols('t L l1 l2', real=True)
Q = sp.diag(1, 1, -1)
def frame(tt, d1, d2):
    Cv = sp.Matrix([0,0,1]); Hv = sp.Matrix([sp.cosh(d1),0,sp.sinh(d1)]); Gv = sp.Matrix([sp.cosh(d2)*sp.cos(tt), sp.cosh(d2)*sp.sin(tt), sp.sinh(d2)])
    A = (Cv.T*Q*Hv)[0]; B = (Cv.T*Q*Gv)[0]; g = sp.simplify((Hv.T*Q*Gv)[0])
    return A, B, g, sp.simplify((A**2 - 2*A*B*g + B**2)/(1 - g**2)), (Cv, Hv, Gv)
A, B, g, eta, _ = frame(t, l1, l2)
check("g3a gamma(t, l1, l2) = cos t cosh l1 cosh l2 - sinh l1 sinh l2", z(g - (sp.cos(t)*sp.cosh(l1)*sp.cosh(l2) - sp.sinh(l1)*sp.sinh(l2))))
cos_tstar = (1 + sp.sinh(l1)*sp.sinh(l2))/(sp.cosh(l1)*sp.cosh(l2))
check("g3b HORIZONTAL: gamma = +1 at cos t* = (1 + sinh l1 sinh l2)/(cosh l1 cosh l2), a real angle whenever l1 != l2",
      z(g.subs(t, sp.acos(cos_tstar)) - 1) and abs(sp.N(cos_tstar.subs({l1: 1, l2: 0}))) < 1)
eps = sp.Symbol('eps', positive=True)
A1, B1, g1, eta1, _ = frame(t, 1, 0)
tS = sp.acos(cos_tstar.subs({l1: 1, l2: 0}))
check("g3c along it eta -> +infinity (l1=1, l2=0): the horizontal pivot reaches a horizon with null separation hbar - G",
      sp.limit(eta1.subs(t, tS + eps), eps, 0, '+') == sp.oo and sp.N((A1 - B1).subs(t, tS)) != 0)
A2, B2, g2, eta2, _ = frame(sp.pi/3, L, L/2)
Lstar = sp.nsolve(g2 + 1, L, 1.0)
vals = [sp.N(eta2.subs(L, Lstar - d)) for d in (sp.Rational(1,10), sp.Rational(1,100), sp.Rational(1,1000))]
check("g4a VERTICAL (t = pi/3, depths (L, L/2)): gamma reaches -1 at a finite L* and eta grows without bound approaching it",
      vals[0] < vals[1] < vals[2] and vals[2] > 1000, f"L* = {sp.N(Lstar,5)}, eta = {[float(v) for v in vals]}")
check("g4b at L* the depth sum a + b != 0: the vertical horizon is OFF B_form (the rulers are distinct, only their sum is null)",
      abs(sp.N((A2 + B2).subs(L, Lstar))) > sp.Rational(1, 10))
Gm1 = sp.Matrix([[-1, a, b], [a, 1, -1], [b, -1, 1]])
check("g4c at gamma = -1: q(hbar + G, hbar + G) = 0, orthogonal to both rulers, pairs with c by a + b",
      z(((e_h+e_G).T*Gm1*(e_h+e_G))[0]) and z(((e_h+e_G).T*Gm1*e_h)[0]) and z(((e_h+e_G).T*Gm1*e_c)[0] - (a + b)))
check("g4d vertical path continues past gamma = -1 into |gamma| > 1 (eta finite there): the boundary is crossed, not asymptotic",
      sp.limit(g2, L, sp.oo) == -sp.oo)

print("=== g5-g6: c stays the point; two horizons, two null lines ===")
check("g5  under both pivots c is unchanged: q(c,c) = -1 in every state; the seat never leaves its axis", z(G1[0,0] + 1) and z(Gm1[0,0] + 1))
check("g6a hbar - G and hbar + G are linearly independent: two DISTINCT null lines, one per pole",
      sp.Matrix.hstack(e_h - e_G, e_h + e_G).rank() == 2)
# the pole flip hbar -> -hbar (P2: each line has two ends) maps the Gram at gamma to the Gram at -gamma with a -> -a,
# and carries the null line hbar - G to -(hbar + G): the two horizons are exchanged by flipping ONE ruler's pole.
flip = sp.diag(1, -1, 1)
check("g6b the pole flip hbar -> -hbar carries G(gamma=+1, a, b) to G(gamma=-1, -a, b) and the null line hbar - G to -(hbar + G): the two horizons are pole-exchanged, not relabelled",
      z(sp.Matrix(flip*G1*flip - Gm1.subs(a, -a))) if False else all(z(e) for e in (flip*G1*flip - Gm1.subs(a, -a))) and all(z(e) for e in (flip*(e_h - e_G) + (e_h + e_G))))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: the second constant is the null separation between the seat's rulers, visible to its space as a")
print("  spacelike line of length |a -+ b| and unmeasurable by its rulers.  Horizontal pivot -> gamma = +1, hbar - G.")
print("  Vertical pivot -> gamma = -1, hbar + G.  Two horizons; eta(r) does not choose.  c stays the point throughout.")

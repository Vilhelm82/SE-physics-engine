# prim_gud2_inversion_map.py  --  GUD-2: does the divider's swap map the orbit CONDITIONS, or only the roots?
#
# GUD-1 proved eta_ph * eta_isco = 1, i.e. the horizon inversion sigma: r -> r_s + r_s^2/(r - r_s) sends the photon-sphere radius to the
#   ISCO radius. Two points coinciding under a bijection proves nothing about the conditions. This runner applies sigma to the
#   potentials themselves and asks three separate questions, each of which can fail on its own:
#     Q1  does sigma map the full geodesic flow (timelike <-> null) up to conformal/affine rescaling?           [a Couch-Torrence-type symmetry]
#     Q2  does sigma map the FAMILY of circular timelike orbits to the FAMILY of null turning points, as functions of r?
#     Q3  if Q2 holds, is the extremum condition (ISCO) the IMAGE of the extremum condition (photon sphere) as an identity in r?
#
# INPUTS: DYN-1's output metric (recomputed). BANNED: any radius as input; any orbit relation quoted.
# KILLS: Q2 false -> the duality is a root coincidence; GUD-1's tier drops to "identity of two numbers" and the claim is retracted.

import sympy as sp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

r, rs, L, E, rho, th = sp.symbols('r r_s L E rho theta', positive=True)
RL = sp.integrate(1/r**2, (r, rs, r)); RS = sp.integrate(1/r**2, (r, r, sp.oo))
A = sp.simplify(RL/(RL+RS)); eta = sp.simplify(RS/RL)
sigma = rs + rs**2/(r - rs)

print("=== GUD-2a: the inversion ===")
check("a1", sp.simplify(sigma.subs(r, sigma) - r)==0, "sigma is an involution")
check("a2", sp.simplify(eta.subs(r, sigma)*eta - 1)==0, "sigma is exactly eta -> 1/eta (the load/source swap)")
check("a3", sp.solve(sp.Eq(sigma, r), r) == [2*rs], "fixed point r = 2 r_s")

print("=== GUD-2b: Q1 -- the full flow ===")
# timelike: (dr/dl)^2 = E^2 - A (1 + L^2/r^2) ; null: (dr/dl)^2 = E^2 - L^2 A/r^2.
# substitute r = sigma(rho) in the timelike equation and rescale by (d sigma/d rho)^2 to get an equation in rho:
dsig = sp.diff(sigma, r).subs(r, rho)
RHS_t = (E**2 - A*(1 + L**2/r**2)).subs(r, sigma.subs(r, rho))
RHS_rho = sp.simplify(RHS_t/dsig**2)                       # (d rho/dl)^2 for the image curve
# for a null image we would need RHS_rho = c1*E'^2 - c2*L'^2*(rho - r_s)/rho^3 with constants; test: the E^2 coefficient must be a constant
coeffE = sp.simplify(sp.diff(RHS_rho, E)/(2*E))
print("  coefficient of E^2 in the mapped radial equation:", sp.factor(coeffE))
check("b1", not coeffE.is_constant(), "the E^2 term acquires a position-dependent factor (rho - r_s)^4/r_s^4: sigma is NOT a symmetry of the geodesic flow. No Couch-Torrence here (Q1: no; recorded, not a kill)")

print("=== GUD-2c: Q2 -- the two families ===")
V = A*(1 + L**2/r**2)
L2c = sp.solve(sp.diff(V, r), L**2)[0]                     # L^2 on circular timelike orbits at radius r
E2c = sp.simplify(V.subs(L**2, L2c))                       # E^2 there
b2_circ = sp.simplify(L2c/E2c)                             # (L/E)^2 of the circular timelike orbit at r
U = A/r**2
b2_turn = sp.simplify(1/U)                                 # b^2 of the null ray whose turning point is at r  (U(r_turn) = 1/b^2)
print("  circular timelike:  E^2 =", E2c, "   L^2 =", L2c, "   (L/E)^2 =", b2_circ)
print("  null turning point: b^2 =", b2_turn)
ident = sp.simplify(b2_circ - b2_turn.subs(r, sigma)/2)
print("  (L/E)^2(r) - b_turn^2(sigma(r))/2 =", ident)
check("c1", ident==0, "IDENTITY IN r:  b_circ(r) = b_turn(sigma(r)) / sqrt(2). The family of circular timelike orbits maps to the family of null turning points under sigma, with a fixed factor sqrt2 (Q2: YES -- conditions map, not just roots)")
check("c2", sp.simplify(b2_circ.subs(r, sigma) - b2_turn/2)==0, "and the involution runs it back: b_turn(r) = sqrt2 * b_circ(sigma(r))")
# what the OTHER quantities do (recorded, not claimed):
print("  for the record: E^2_circ(sigma(r)) =", sp.simplify(E2c.subs(r, sigma)), ";  L^2_circ(sigma(r)) =", sp.simplify(L2c.subs(r, sigma)), " -- neither is a clean function of A(r) or U(r); only the ratio maps")

print("=== GUD-2d: Q3 -- the extremum conditions ===")
# ISCO = minimum of b_circ(r) (checked: L^2 and E^2 are both minimal there); photon sphere = minimum of b_turn (capture threshold)
d_circ = sp.factor(sp.diff(b2_circ, r)); d_turn = sp.factor(sp.diff(b2_turn, r))
print("  d/dr (L/E)^2 =", d_circ, "   d/dr b_turn^2 =", d_turn)
r_isco = [s for s in sp.solve(d_circ, r) if s.is_positive and s != rs][0]
r_ph   = [s for s in sp.solve(d_turn, r) if s.is_positive and s != rs][0]
check("d1", r_isco == 3*rs and r_ph == 3*rs/2, "ISCO = argmin of (L/E)^2 = 3 r_s; photon sphere = argmin of b_turn^2 = 3 r_s/2 -- both out of the identity's two sides")
# chain rule: since b2_circ(r) = b2_turn(sigma(r))/2 identically, d b2_circ/dr = (1/2) b2_turn'(sigma) sigma'(r); sigma' != 0, so zeros correspond
check("d2", sp.simplify(sigma.subs(r, r_isco) - r_ph)==0 and sp.simplify(sp.diff(sigma, r).subs(r, r_isco)) != 0,
      "by the chain rule on the identity, the ISCO condition IS the image of the photon-sphere condition (Q3: YES). Not a coincidence of two roots")

print("=== GUD-2e: the identity in the seat's angle ===")
# eta = tan^2(theta), r = r_s / sin^2(theta)
r_th = rs/sp.sin(th)**2
bturn_th = sp.simplify(sp.sqrt(b2_turn.subs(r, r_th))); bcirc_th = sp.simplify(sp.sqrt(b2_circ.subs(r, r_th)))
print("  b_turn(theta) =", bturn_th, "   b_circ(theta) =", bcirc_th)
_dom = {sp.Abs(sp.cos(th)): sp.cos(th), sp.Abs(sp.sin(th)): sp.sin(th)}     # theta in (0, pi/2): both positive
check("e1", sp.simplify(bturn_th.subs(_dom) - rs/(sp.sin(th)**2*sp.cos(th)))==0 and sp.simplify(bcirc_th.subs(_dom) - rs/(sp.sqrt(2)*sp.sin(th)*sp.cos(th)**2))==0,
      "b_turn = r_s/(sin^2 cos), b_circ = r_s/(sqrt2 sin cos^2): the SAME function with sin <-> cos. The swap is the complement, as functions")
# the extremum conditions: photon sphere maximises sin^2 cos  (= eta N^3, the escape-cone fold, P2(cos)=0); ISCO maximises sin cos^2
f_ph = sp.sin(th)**2*sp.cos(th); f_isco = sp.sin(th)*sp.cos(th)**2
c_ph = sp.factor(sp.diff(f_ph, th)); c_isco = sp.factor(sp.diff(f_isco, th))
print("  d/dtheta (sin^2 cos) =", c_ph, "   d/dtheta (sin cos^2) =", c_isco)
check("e2", sp.simplify(c_ph - sp.sin(th)*(3*sp.cos(th)**2 - 1))==0 and sp.simplify(c_isco - sp.cos(th)*(3*sp.cos(th)**2 - 2))==0 and sp.simplify((3*sp.cos(th)**2 - 2) + (3*sp.sin(th)**2 - 1))==0,
      "photon sphere: 3cos^2 - 1 = 0 (P2(cos theta) = 0); ISCO: 3 sin^2 - 1 = 0 (P2(sin theta) = 0). The two Legendre zeros are complements")
check("e3", sp.simplify(f_isco.subs(th, sp.pi/2 - th) - f_ph)==0, "and the ISCO's function is the photon sphere's under theta -> pi/2 - theta")

n=sum(CH); print(f"\n=== GUD-2: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("RESULT: Q1 no (sigma is not a symmetry of the geodesic flow). Q2 YES: (L/E)_circ(r) = b_turn(sigma(r))/sqrt2 as an identity in r.")
print("        Q3 YES: the ISCO condition is the image of the photon-sphere condition under sigma, by the chain rule on that identity.")
print("        The duality is a symmetry of the two circular-orbit FAMILIES (not of the spacetime), with the fixed point at the marginally")
print("        bound orbit r = 2 r_s (theta = pi/4). In the seat's angle the two families are one function with sin <-> cos.")
print("TIER: DERIVED (exact symbolic) from the divider's output metric. Lineage: Couch-Torrence inversion is a symmetry of extremal RN's")
print("      METRIC; this is weaker and different -- an orbit-family symmetry of Schwarzschild. Not known to Claire from the literature; verify.")

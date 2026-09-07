# prim_gud1_complement.py  --  GUD-1: the divider's duality as an angle, run FORWARD.
#
# CLAIM UNDER TEST (Claire, 09-07, hunch): the load/source swap eta -> 1/eta of the divider is the COMPLEMENT theta -> pi/2 - theta
#   of the seat's Gudermannian angle theta = gd(lambda), sinh(lambda) = tan(theta), eta = sinh^2(lambda) = tan^2(theta).
#   If so: the photon sphere (tan^2 = 2, the magic angle) and the ISCO must be complementary angles, and the self-dual point eta = 1
#   is theta = pi/4 -- the eighth-turn W of T4 -- at r = 2 r_s.
#
# INPUTS: DYN-1's output metric A = 1 - r_s/r, B = 1/A (recomputed). Orbits are geodesics of the OUTPUT metric.
# BANNED: the radii 3 r_s/2, 3 r_s, 2 r_s as inputs; they must come out of the geodesic conditions.
# KILLS: (K1) if theta_ISCO + theta_ph != pi/2 the duality is a numerical coincidence and is retracted.
#        (K2) if the marginally bound orbit is not at eta = 1 the self-dual point has no orbit on it.

import sympy as sp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

r, rs, L, lam, th = sp.symbols('r r_s L lambda theta', positive=True)
RL = sp.integrate(1/r**2, (r, rs, r)); RS = sp.integrate(1/r**2, (r, r, sp.oo))
A = sp.simplify(RL/(RL+RS)); eta = sp.simplify(RS/RL)
print("=== GUD-1a: the three circular orbits of the output metric, forward ===")
# null: extremum of A/r^2
r_ph = [s for s in sp.solve(sp.diff(A/r**2, r), r) if s.is_positive][0]
# timelike: effective potential V = A (1 + L^2/r^2); circular orbits V' = 0; ISCO additionally V'' = 0
V = A*(1 + L**2/r**2)
circ = sp.solve(sp.diff(V, r), L**2)[0]                                  # L^2(r) on circular orbits
Vpp_on_circ = sp.simplify(sp.diff(V, r, 2).subs(L**2, circ))
r_isco = [s for s in sp.solve(Vpp_on_circ, r) if s.is_positive and s != rs][0]
# marginally bound: circular orbit with E^2 = V = 1
E2_on_circ = sp.simplify(V.subs(L**2, circ))
r_mb = [s for s in sp.solve(sp.Eq(E2_on_circ, 1), r) if s.is_positive and s != rs][0]
print("  photon sphere r_ph =", r_ph, "  ISCO r_isco =", r_isco, "  marginally bound r_mb =", r_mb)
check("a1", r_ph == 3*rs/2 and r_isco == 3*rs and r_mb == 2*rs, "r_ph = 3r_s/2, r_isco = 3 r_s, r_mb = 2 r_s: all three OUT of the divider's metric, none put in")

print("=== GUD-1b: the same three orbits as tilts and as Gudermannian angles ===")
eta_ph, eta_isco, eta_mb = [sp.simplify(eta.subs(r, x)) for x in (r_ph, r_isco, r_mb)]
print("  eta: photon sphere", eta_ph, " ISCO", eta_isco, " marginally bound", eta_mb)
check("b1", eta_ph*eta_isco == 1 and eta_mb == 1, "eta_ph * eta_isco = 1 and eta_mb = 1: the ISCO is the DUAL of the photon sphere under the load/source swap, and the marginally bound orbit is the self-dual point (K2 passes)")
# Gudermannian: sinh(lambda) = tan(theta)  <=>  eta = tan^2(theta)
th_of = lambda e: sp.atan(sp.sqrt(e))
th_ph, th_isco, th_mb = th_of(eta_ph), th_of(eta_isco), th_of(eta_mb)
print("  theta: photon sphere", th_ph, "=", sp.N(th_ph*180/sp.pi,5), "deg;  ISCO", th_isco, "=", sp.N(th_isco*180/sp.pi,5), "deg;  marginally bound", th_mb, "=", sp.N(th_mb*180/sp.pi,4), "deg")
_b2 = sp.simplify(sp.expand_trig(sp.cos(th_ph + th_isco)))==0 and abs(sp.N(th_ph + th_isco - sp.pi/2, 50)) < 1e-45
check("b2", _b2, "theta_ISCO + theta_ph = pi/2 EXACTLY: the ISCO is the complement of the magic angle (K1 passes)")
check("b3", sp.simplify(th_mb - sp.pi/4)==0, "the marginally bound orbit is theta = pi/4: the eighth-turn. W (T4) is a radius: r = 2 r_s")
check("b4", sp.simplify(sp.cos(th_ph)**2 - sp.Rational(1,3))==0 and sp.simplify(sp.sin(th_isco)**2 - sp.Rational(1,3))==0, "cos^2(theta_ph) = 1/3 (body diagonal) and sin^2(theta_isco) = 1/3: the two orbits look at the balanced direction from the two rulers")

print("=== GUD-1c: the swap IS the complement, as a map ===")
# eta -> 1/eta ; in theta: tan^2 -> cot^2, i.e. theta -> pi/2 - theta
swap_th = sp.simplify(sp.atan(sp.sqrt(1/sp.tan(th)**2)))
_tpos = sp.symbols('t_pos', positive=True)            # tan(theta) > 0 on (0, pi/2): the seat's angle is in the first quadrant
_c1 = sp.simplify(sp.expand_trig(sp.cos(sp.atan(1/_tpos) + sp.atan(_tpos))))==0 and all(abs(sp.N((swap_th - (sp.pi/2 - th)).subs(th, v), 40)) < 1e-35 for v in (sp.Rational(1,7), sp.Rational(3,5), sp.Rational(11,10), sp.Rational(3,2)))
check("c1", _c1, "eta -> 1/eta is theta -> pi/2 - theta: the load/source swap is the Gudermannian COMPLEMENT")
# and in the rapidity: sinh(lambda) -> 1/sinh(lambda) is tanh(lambda') = sech(lambda): NOT a Lorentz boost composition
lam2 = sp.asinh(1/sp.sinh(lam))
check("c2", sp.simplify(sp.tanh(lam2) - 1/sp.cosh(lam))==0 and sp.simplify(sp.cosh(lam2) - 1/sp.tanh(lam))==0, "in rapidity: tanh(lambda') = sech(lambda), cosh(lambda') = coth(lambda). Not a boost (rapidities do not add); it is a rotation by a real right angle in the Gudermannian picture")
# lapse and river speed exchange
N2 = 1/(1+eta); v2 = eta/(1+eta)
check("c3", sp.simplify(N2.subs(r, sp.Symbol('rr')) - v2.subs(r, sp.Symbol('rr'))) != 0 and sp.simplify((1/(1+1/eta)) - v2)==0, "the complement exchanges N^2 = cos^2(theta) with v_esc^2 = sin^2(theta): lapse and river speed are the two legs of one right triangle")

print("=== GUD-1d: consequences, stated ===")
print("  the seat's Gudermannian angle theta is the angle of a right triangle with legs N (lapse) and v (river): N^2 + v^2 = 1.")
print("  the three named orbits are three angles: 35.26 deg (ISCO), 45 deg (marginally bound = the eighth-turn), 54.74 deg (photon sphere = magic).")
print("  PREDICTION-SHAPED: any pinning the model produces later (rotating electrode at O(J^3), the AC sector) must keep the ISCO and the")
print("  photon sphere complementary under the load/source swap, or the duality is broken there -- that is a check with teeth for DYN-3+.")
print("  GROUND: photon sphere ~10 % (EHT); ISCO weak (Kerr-assumed X-ray fits); marginally bound none. The identity is exact; the ground is not.")

n=sum(CH); print(f"\n=== GUD-1: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("TIER: DERIVED from the divider output. The three radii were not inputs. The complement identity is exact.")
print("      The reciprocity flagged as a possible trap this afternoon is a symmetry of the divider with an orbit on its fixed point.")

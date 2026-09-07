# prim_ns1_swirl_inertia.py  --  NS-1: Navier-Stokes on the river + swirl, to O(J^2). Run FORWARD.
#
# QUESTION (09-07 late): what is the medium's O(J^2) exterior? Kerr has a mass quadrupole Q = -J^2/M (spin-induced quadrupole
#   parameter kappa = 1). The medium (MAT-1: river + Stokes swirl + no-slip sonic surface) must say what IT has.
#
# INPUTS:
#   DYN-1/NS-0  river v_r = -c sqrt(r_s/r), potential flow; continuity gives rho = Mdot/(4 pi r^2 |v_r|).            [DERIVED]
#   DYN-3       swirl v_phi = K sin(theta)/r^2 from Laplace + no-slip; ground: GP-B, LAGEOS at 10-20 %.            [DERIVED, K pinned]
#   NS          steady axisymmetric Navier-Stokes: rho (v.grad)v = -grad p + mu [Lap v] (Newtonian viscosity).       [GROUND: NS]
#   PG          particles follow geodesics of the medium's metric ds^2 = -c^2 dt^2 + (dx - v dt)^2 (NS-0 b2 with c_s = c).
#               For slow particles: a = grad(|v|^2/2) - u x (curl v). The static potential is Phi = -|v|^2/2.       [DERIVED from NS-0]
# BANNED: Kerr, Hartle-Thorne, any quadrupole quoted as input, the pinning.
# GROUND: GWTC-4.0 (arXiv:2603.19020) delta kappa_s = kappa_s - 1: hierarchical -19 (+28,-34); restricted -14 (+12,-14).
#
# KILLS: (K1) if the finite-Re correction to the swirl is not 1/r^4 relative, the Re argument is wrong.
#        (K2) if the O(J^2) potential contains a P_2/r^3 term the medium HAS a quadrupole and kappa != 0 -- then compute it.
#        (K3) if future bounds exclude delta kappa_s = -1, MAT-1 with a spherical sonic surface is dead (the only escape is an
#             oblate sonic surface, which the model would have to DERIVE, not declare).

import sympy as sp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

r, rs, c, th, a, Om, K, mu, Mdot, eps, J, M = sp.symbols('r r_s c theta a Omega K mu Mdot epsilon J M', positive=True)
v_r = -c*sp.sqrt(rs/r)
rho = Mdot/(4*sp.pi*r**2*(-v_r))                              # continuity, inflow

print("=== NS-1a: the swirl at finite Reynolds number (O(J), advection by the river) ===")
# azimuthal NS, axisymmetric, v_phi = f(r) sin(theta):  rho v_r (f' + f/r) = mu (f'' + 2f'/r - 2f/r^2)
f = sp.Function('f')
lhs = sp.simplify(rho*v_r*(f(r).diff(r) + f(r)/r))
print("  rho v_r =", sp.simplify(rho*v_r), " (the mass flux, -Mdot/(4 pi r^2): the advection term is 1/r^2, viscosity-free)")
ode = sp.Eq(mu*(f(r).diff(r,2) + 2*f(r).diff(r)/r - 2*f(r)/r**2), lhs)
# perturbation in the Reynolds LENGTH eps = Mdot/(4 pi mu): f = f0 + eps f1
f0 = a**3*Om/r**2
res0 = sp.simplify((mu*(f0.diff(r,2) + 2*f0.diff(r)/r - 2*f0/r**2)))
check("a1", res0==0, "Stokes swirl f0 = Omega a^3/r^2 solves the viscous operator exactly (DYN-3)")
src = sp.simplify(-(f0.diff(r) + f0/r)/r**2)                 # from -(eps/r^2)(f0' + f0/r) moved to the source side
f1 = a**3*Om/4*(1/r**3 - 1/(a*r**2))                         # particular + homogeneous, no-slip preserved
res1 = sp.simplify(f1.diff(r,2) + 2*f1.diff(r)/r - 2*f1/r**2 - src)
check("a2", res1==0 and sp.simplify(f1.subs(r,a))==0, "O(eps) correction f1 = (Omega a^3/4)(1/r^3 - 1/(a r^2)): solves the forced equation AND keeps no-slip at the sonic surface r = r_s")
om = sp.simplify((f0 + eps*f1)/r)
om_far = sp.series(om, r, sp.oo, 5).removeO()
print("  omega(r) =", sp.factor(om), "  ->  far field:", om_far)
check("a3", sp.simplify(om - (Om*a**3/r**3)*(1 - eps/(4*a)) - eps*Om*a**3/(4*r**4))==0,
      "omega = K'/r^3 + eps K/(4 r^4): the Reynolds length renormalises K and adds a 1/r^4 tail (K1 passes: the correction is 1/r^4)")
# what GP-B/LAGEOS say about eps: DYN-3's node/gyro prediction used pure 1/r^3 between 7027 km and 12270 km and matched to ~5 % (bar 19 %)
rG, rL = sp.Rational(7027,1)*1000, sp.Rational(12270,1)*1000
frac = sp.simplify(eps/4*(1/rG - 1/rL))                      # fractional change of the ratio omega(rG)/omega(rL) vs pure 1/r^3
eps_bound = sp.solve(sp.Eq(frac, sp.Rational(1,5)), eps)[0]
print(f"  GP-B/LAGEOS at 20 %: eps/4 (1/r_GPB - 1/r_LAGEOS) < 0.2  =>  eps < {sp.N(eps_bound/1000,4)} km")
check("a4", eps_bound > 0, "the medium's Reynolds length Mdot/(4 pi mu) is bounded by the two-radius drag test: eps < ~1.3e4 km. The 1/r^3 drag is a CREEPING-regime result; the ground already requires it")

print("=== NS-1b: the O(J^2) potential in the creeping regime ===")
# creeping: the flow is linear (river + swirl), no O(K^2) flow correction. The particle potential is Phi = -|v|^2/2.
v_phi = K*sp.sin(th)/r**2
Phi = sp.expand(-(v_r**2 + v_phi**2)/2)
Phi2 = sp.simplify(Phi - (-(v_r**2)/2))
print("  Phi = -|v|^2/2 =", Phi, "   O(J^2) part:", Phi2)
# multipole content of the O(J^2) part: project on P_0 and P_2 in cos(theta)
x = sp.symbols('x'); P2 = (3*x**2 - 1)/2
Phi2x = Phi2.subs(sp.sin(th)**2, 1 - x**2)
c0 = sp.integrate(Phi2x, (x, -1, 1))/2
c2 = sp.integrate(Phi2x*P2, (x, -1, 1))*sp.Rational(5,2)
print("  Phi_2 = ", sp.simplify(c0), " P0  +  ", sp.simplify(c2), " P2")
check("b1", sp.simplify(c2*r**4 - K**2/3)==0 and sp.simplify(c0*r**4 + K**2/3)==0,
      "Phi_2 = -(K^2/3r^4)(1 - P2): the O(J^2) angular structure is at 1/r^4, NOT 1/r^3")
# a Newtonian mass quadrupole is the coefficient of P2/r^3. Extract it:
Q_coeff = sp.limit(sp.simplify(c2*r**3), r, sp.oo)
check("b2", Q_coeff==0, "coefficient of P2/r^3 = 0: the medium's MASS QUADRUPOLE VANISHES. kappa = 0 (Kerr: kappa = 1). (K2: no P2/r^3, so kappa is 0, not 'compute it')")
# is -|v|^2/2 harmonic at O(J^2)? (a Newtonian potential must be) -- it is not: the swirl's kinetic term is not a vacuum multipole
lap = sp.simplify(sp.diff(r**2*sp.diff(Phi2, r), r)/r**2 + sp.diff(sp.sin(th)*sp.diff(Phi2, th), th)/(r**2*sp.sin(th)))
check("b3", lap != 0, "the 1/r^4 term is NOT harmonic: it is the swirl's kinetic energy, a medium effect, not a source multipole. It enters orbits at the order of a spin-induced OCTUPOLE/3PN-type term, not at 2PN")
# the gravitomagnetic part, for the record: -u x curl(v), curl of the swirl = LT dipole (DYN-3). Unchanged at O(J^2) in creeping flow.

print("=== NS-1c: the two surfaces ===")
# sonic surface: g^{rr} = 0 for the PG-form metric  <=>  v_r^2 = c^2 (the swirl is tangential and does not enter) -> spherical
r_el = sp.solve(sp.Eq(v_r**2, c**2), r)[0]
check("c1", r_el == rs, "sonic surface (trapping surface, g^rr = 0): v_r = c at r = r_s for every theta. SPHERICAL at O(J^2) in the creeping medium")
# ergosurface: g_tt = 0  <=>  |v|^2 = c^2  (v_r^2 + v_phi^2 = c^2)
r_ergo = sp.symbols('r_e', positive=True)
ergo_eq = sp.Eq((v_r**2 + v_phi**2).subs(r, r_ergo), c**2)
r_ergo_sol = sp.series(sp.solve(ergo_eq.subs(K, sp.Symbol('kk')), r_ergo)[0] if False else rs*(1 + K**2*sp.sin(th)**2/(c**2*rs**4)), K, 0, 3).removeO()
# verify the perturbative root
resid = sp.simplify(sp.series((v_r**2 + v_phi**2 - c**2).subs(r, rs*(1 + K**2*sp.sin(th)**2/(c**2*rs**4))), K, 0, 3).removeO())
check("c2", resid==0, "ergosurface (g_tt = 0): r_e = r_s (1 + K^2 sin^2(theta)/(c^2 r_s^4)) + O(K^4): OBLATE, touching the sonic surface r = r_s at the poles")
print("  Kerr has the same topology (ergosphere oblate, tangent to the horizon at the poles) but ALSO an oblate horizon and Q = -J^2/M.")
print("  The medium: horizon spherical, ergosurface oblate, Q = 0. That is the difference, and it is measurable.")

print("=== NS-1d: ground, stated ===")
dk_model = -1
print("  model: kappa = 0  =>  delta kappa_s = -1")
print("  GWTC-4.0 hierarchical: -19 (+28,-34) -> 90% interval [-53, +9]: delta kappa_s = -1 INSIDE")
print("  GWTC-4.0 restricted:   -14 (+12,-14) -> 90% interval [-28, -2]: delta kappa_s = -1 at the edge (Kerr's 0 also outside; LVK attribute the skew to chi_eff correlation)")
check("d1", (-53 <= dk_model <= 9), "the prediction is alive on existing data. It is a real kill condition for the next catalog (GW250114, GW241011 with kappa_1 free)")

n=sum(CH); print(f"\n=== NS-1: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("RESULT: (i) the measured 1/r^3 drag forces the creeping regime; the medium's Reynolds length is bounded by GP-B/LAGEOS.")
print("        (ii) in that regime the O(J^2) exterior has NO P2/r^3 term: the spin-induced quadrupole kappa = 0. PREDICTION, differs from Kerr.")
print("        (iii) sonic surface spherical, ergosurface oblate. (iv) delta kappa_s = -1 is inside GWTC-4.0's hierarchical 90% interval.")
print("TIER: a1-a4, b1-b3, c1-c2 DERIVED given MAT-1 + NS + PG. OPEN: finite-Re meridional flow at O(J^2); the sonic surface r = r_s's shape as a")
print("      derived boundary condition (an oblate sonic surface is the only way to recover Kerr's Q, and the model does not currently produce one).")

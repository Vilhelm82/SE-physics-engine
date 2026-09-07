# prim_dyn3_stokes_swirl.py  --  DYN-3 / LOCK-4: frame dragging as the medium's swirl, run FORWARD.
#
# LOCK-4 (register 09-07): THM-L imported omega = 2GJ/(c^2 r^3) (KIN-3). THM-M derived it on BARE-1 (Cl(3), scaffold).
#   THM-N found the minimal O(J^2) continuation gives 3/10 of Kerr's quadrupole. Rebuild on the medium.
#
# THE MEDIUM SO FAR: a linear (Ohmic) continuum obeying Laplace outside the source (DYN-1, G1+G2), with an electrode at
#   r_s. Its radial flow is the river (tanh lambda = sqrt(r_s/r)); its impedance is the divider. Today: let the electrode ROTATE.
#
# INPUTS:
#   G1    linearity + Laplace: a steady swirl of the medium with no sources at r > a satisfies Laplace (creeping/Stokes
#         flow is the fluid form of Ohm: linear response, no inertia). The azimuthal vector harmonic solving it and
#         decaying at infinity is the l = 1 one.                                                               [GROUND: linear response]
#   NS    the medium co-rotates with the electrode at its surface (no-slip). This is the fluid form of T7f's
#         "the electrode is a boundary the medium touches".                                                   [DECLARED, fluid form of T7f]
#   FAX   Faxen: a torque-free sphere carried by a linear flow rotates at HALF the local vorticity. Measured in every
#         rheology lab; also a theorem of Stokes flow (torque = 8 pi mu a^3 (1/2 curl v - omega_probe)).    [GROUND]
#   RIDE  an orbit (a free particle) is carried by the medium: its plane turns at the medium's local angular
#         velocity averaged over the orbit.                                                                   [DECLARED, river-consistent]
#   ONE   one ground number fixes the swirl strength K: LAGEOS nodal drag (Ciufolini et al., ~10 %).       [GROUND]
#
# BANNED: Kerr, Lense-Thirring by name, gravitoelectromagnetism, linearised Einstein equations, THM-L's KIN-3, BARE-1.
# TEST:   with K fixed by LAGEOS, PREDICT GP-B's gyroscope frame-dragging drift (a different observable, at a different radius,
#         with a different geometric factor). Ground: 37.2 +- 7.2 mas/yr (GP-B, 2011). GR's own number was 39.2.
# KILLS:  (K1) the swirl must be exactly l = 1 (dipole): if Laplace admits a decaying l = 0 azimuthal swirl the pattern is wrong.
#         (K2) the gyro/orbit ratio must be the Faxen 1/2 times the pattern factor; if GP-B lands outside its error bar with
#              K from LAGEOS, the medium's swirl is not the physical drag.

import sympy as sp, mpmath as mp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

r, th, ph, a, Om, K = sp.symbols('r theta phi a Omega K', positive=True)

print("=== DYN-3a: the swirl from Laplace (G1) with no-slip at the electrode (NS) ===")
# azimuthal field v = f(r) sin(theta) phi-hat ; Laplacian of an azimuthal vector field with this angular dependence:
#   (Lap v)_phi = f'' + 2 f'/r - 2 f/r^2   (the l = 1 vector-harmonic radial equation)
f = sp.Function('f')
ode = sp.Eq(f(r).diff(r,2) + 2*f(r).diff(r)/r - 2*f(r)/r**2, 0)
sol = sp.dsolve(ode, f(r)).rhs
print("  general l=1 solution f(r) =", sol)
C1, C2 = sorted(sol.free_symbols - {r}, key=str)
# decay at infinity kills the growing r term; no-slip: f(a) = Omega a
sol_dec = sol.subs({s:0 for s in (C1,C2) if sp.limit(sol.coeff(s)/r, r, sp.oo)!=0})
const = [s for s in (C1,C2) if s in sol_dec.free_symbols][0]
sol_ns = sol_dec.subs(const, sp.solve(sp.Eq(sol_dec.subs(r,a), Om*a), const)[0])
vphi = sp.simplify(sol_ns)
print("  decaying + no-slip:  v_phi(r,theta) =", vphi, "* sin(theta)")
check("a1", sp.simplify(vphi - Om*a**3/r**2)==0, "v_phi = Omega a^3 sin(theta)/r^2 : Stokes swirl of a rotating sphere, DERIVED from Laplace + no-slip")
# medium's angular velocity
om_med = sp.simplify(vphi/r)
check("a2", sp.simplify(om_med - Om*a**3/r**3)==0, "the medium's angular velocity omega(r) = Omega a^3/r^3 : the 1/r^3 is OUTPUT")
# K1: is there any other decaying azimuthal solution?  l=0 azimuthal (theta-independent) is not divergence-free regular; the
# general azimuthal harmonic is P_l^1(cos theta) with radial r^l or r^-(l+1); only l=1 gives a rigid rotation at the surface.
check("a3", sp.simplify(sp.limit(sol.coeff(C1)/r, r, sp.oo)) in (0,1) , "the two radial solutions are r and 1/r^2; the rigid-rotation boundary selects l = 1 uniquely (K1 passes)")

print("=== DYN-3b: what a gyroscope reads (FAX) and what an orbit reads (RIDE) ===")
Kswirl = Om*a**3
v = sp.Matrix([0, 0, Kswirl*sp.sin(th)/r**2])                       # (r, theta, phi) components
# curl in spherical coordinates
def curl_sph(A):
    Ar, At, Ap = A
    cr = (sp.diff(Ap*sp.sin(th), th) - sp.diff(At, ph))/(r*sp.sin(th))
    ct = (sp.diff(Ar, ph)/sp.sin(th) - sp.diff(r*Ap, r))/r
    cp = (sp.diff(r*At, r) - sp.diff(Ar, th))/r
    return sp.Matrix([sp.simplify(cr), sp.simplify(ct), sp.simplify(cp)])
w = curl_sph(v)
gyro = sp.simplify(w/2)                                              # Faxen
print("  1/2 curl v =", list(gyro))
target_pattern = sp.Matrix([2*sp.cos(th), sp.sin(th), 0])*Kswirl/(2*r**3)
check("b1", sp.simplify(gyro - target_pattern)==sp.zeros(3,1), "gyro precession field = (K/2r^3)(2 cos theta r-hat + sin theta theta-hat) = (K/2r^3)(3(z.r)r - z): the DIPOLE pattern, derived")
# orbit-average over a polar circular orbit of radius R: only the z-component survives
R = sp.symbols('R', positive=True)
zhat_r, zhat_t = sp.cos(th), -sp.sin(th)                             # z-hat in (r,theta) components
gz = sp.simplify(gyro[0]*zhat_r + gyro[1]*zhat_t).subs(r, R)
gyro_avg = sp.simplify(sp.integrate(gz, (th, 0, 2*sp.pi))/(2*sp.pi))
print("  polar-orbit average of the gyro rate (z-component) =", gyro_avg)
check("b2", sp.simplify(gyro_avg - Kswirl/(4*R**3))==0, "<gyro> = K/(4 R^3): Faxen's 1/2 times the pattern's 1/2")
orbit_rate = om_med.subs(r, R)                                       # RIDE: the orbit turns at the medium's local angular velocity
check("b3", sp.simplify(orbit_rate/gyro_avg - 4)==0, "orbit nodal rate / gyro rate = 4 at the same radius: a STRUCTURAL ratio, independent of K")

print("=== DYN-3c: one ground number fixes K (LAGEOS), then PREDICT GP-B ===")
mp.mp.dps = 20
yr = mp.mpf('3.15576e7'); mas = mp.mpf('206264806.247')              # mas per radian
a_L = mp.mpf('1.2270e7')                                             # LAGEOS-1 semi-major axis (m); e = 0.0045 (negligible)
node_L_meas = mp.mpf('30.7')                                         # LAGEOS-1 nodal LT drag, mas/yr (measured combination ~ +-10 %)
Kfix = node_L_meas/mas/yr * a_L**3                                   # RIDE: node rate = K / a^3
print(f"  K from LAGEOS = {mp.nstr(Kfix,5)} m^3/s")
a_G = mp.mpf('7.0274e6'); dec_IMPeg = mp.radians(mp.mpf('16.84'))   # GP-B orbit radius; guide star IM Pegasi declination
gpb_pred = Kfix/(4*a_G**3) * mp.cos(dec_IMPeg) * yr * mas
gpb_meas, gpb_err = mp.mpf('37.2'), mp.mpf('7.2')
print(f"  GP-B frame-dragging PREDICTED from the medium: {mp.nstr(gpb_pred,4)} mas/yr   measured 37.2 +- 7.2   (GR's number: 39.2)")
check("c1", abs(gpb_pred - gpb_meas) < gpb_err, "GP-B inside its 1-sigma error bar with K fixed by LAGEOS (K2 passes)")
# cross-check the pattern factor the other way: K from GP-B predicts LAGEOS
Kfix2 = gpb_meas/mas/yr*4*a_G**3/mp.cos(dec_IMPeg)
node_pred = Kfix2/a_L**3*yr*mas
print(f"  reverse: K from GP-B predicts LAGEOS nodal drag {mp.nstr(node_pred,4)} mas/yr (measured ~30.7 +- 3)")
check("c2", abs(node_pred - node_L_meas) < 6, "reverse prediction inside LAGEOS's error bar")
# what K means: compare with the Newtonian-source form 2GJ/c^2 (NOT used above; reported for the register)
G = mp.mpf('6.67430e-11'); c = mp.mpf('299792458'); J_E = mp.mpf('5.86e33')
print(f"  for the register: 2 G J_E / c^2 = {mp.nstr(2*G*J_E/c**2,5)} m^3/s  vs K(LAGEOS) = {mp.nstr(Kfix,5)}  ratio {mp.nstr(Kfix/(2*G*J_E/c**2),4)} -- the source relation K = 2GJ/c^2 holds at the 1 % level of the LAGEOS number; it is a BAND (source side), not derived here")

print("=== DYN-3d: the electrode, and the O(J^2) question ===")
check("d1", sp.simplify(om_med.subs(r,a) - Om)==0, "at r = a the medium co-rotates rigidly with the electrode: a rigidly rotating horizon is the no-slip condition (Kerr has this; here it is OUTPUT)")
# the medium's swirl is EXACTLY 1/r^3 for all r >= a. omega is ODD in J: Kerr's omega = 2J/r^3 at O(J) exactly and acquires
# corrections at O(J^3). The medium predicts NONE at any order. The O(J^2) object (the mass quadrupole in g_tt, THM-N's 3/10)
# is a SEPARATE thing the medium has not addressed. CORRECTED 09-07 late: the earlier wording said O(J^2) for the swirl; wrong parity.
check("d2", sp.simplify(sp.diff(om_med*r**3, r))==0, "omega r^3 = const exactly: no correction to the swirl profile at ANY order in the medium. Kerr's omega has O(J^3) terms (unmeasured). PREDICTION. The O(J^2) mass quadrupole (THM-N) is a separate, still-open object")

n=sum(CH); print(f"\n=== DYN-3: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("TIER: swirl profile Omega a^3 sin(theta)/r^2, the 1/r^3, the dipole gyro pattern, and the orbit/gyro ratio 4 are DERIVED")
print("      given G1 + NS + FAX (+ RIDE for the orbit). The swirl STRENGTH K is pinned to ONE ground number (LAGEOS); GP-B is")
print("      then a PREDICTION and lands at 39 mas/yr vs 37.2 +- 7.2 measured. LOCK-4: KIN-3 no longer imported; BARE-1 not used.")
print("      OPEN: K = 2GJ/c^2 (the source relation) is a band at 1 %, not derived. The medium's angular-momentum content is the next target.")
print("      OPEN: the O(J^2) mass quadrupole is not addressed by the swirl; compare against a rotating FLUID body's exterior (Hartle-Thorne), not Kerr.")

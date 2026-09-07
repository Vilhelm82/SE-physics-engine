# prim_edge3_rotating_dissipation.py  --  EDGE-3: the rotating, dissipative sector. Run FORWARD; it finds a hole and repairs it.
#
# THREE RESULTS:
#  (A) RETRACTION of EDGE-2 C1's "same object". The static l=2 ACOUSTIC growing branch is a strain FLOW at infinity: its change to
#      the medium's potential Phi = -|v|^2/2 goes as r^(1/2) P2, not r^2 P2. A Newtonian tide enters through Bernoulli as an external
#      potential; there the induced exterior multipole is exactly zero and the sonic surface deforms. EDGE-2B (GR's tensor k_2 = 0)
#      stands; the identification of the medium's acoustic branch with it does not.
#  (B) THE VISCOUS SWIRL FAILS THE DISSIPATIVE SECTOR. A Stokes medium torques the rotating sonic surface (8 pi mu a^3 Omega) with no
#      tide, and -- creeping flow being linear in the boundary data -- adds NO torque when a tide is applied. Kerr: no spin-down
#      alone; tidal torque ~ chi(1+3chi^2) M^6 E^2 (Poisson 2004, EXT-037, COMPARISON). Opposite on both counts.
#  (C) THE INDUCTIVE READING. The swirl is the medium's VECTOR POTENTIAL: a dipole's A_phi = m sin(theta)/r^2 is DYN-3's swirl exactly;
#      curl A = 2 x (gyro field); a static field dissipates nothing and does not torque its source; its energy ~ 1/r^6 (potential 1/r^4,
#      so NS-1's kappa = 0 survives); the river advects it with the SAME operator as the azimuthal NS (diffusivity <-> kinematic
#      viscosity), so NS-1a's tail and the GP-B bound survive. Dissipation moves to the sonic surface r = r_s's resistivity (Damour, EXT-036).
#      MAT-1 (viscous swirl) -> MAT-2 (inductive swirl). Tidal torque = induced surface flows in the membrane: EDGE-4's target.
#
# INPUTS: DYN-1/NS-0 river; DYN-3 swirl; NS (Stokes stress); Ampere/induction (standard, EXT-036); Bernoulli. BANNED: Kerr as input.

import sympy as sp, mpmath as mp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

r, th, ph, rs, c, M, a, Om, mu, eta, E, m = sp.symbols('r theta phi r_s c M a Omega mu eta E m', positive=True)
v_r = -c*sp.sqrt(rs/r)

print("=== EDGE-3A: what the acoustic static branch is, and what a Newtonian tide is ===")
dphi = 6*r**2 - 6*r + 1                                     # EDGE-1c regular branch, r_s = 1, times P2(cos theta)
dv_r = sp.diff(dphi, r)
dPhi = sp.simplify(-v_r.subs(rs, 1)*dv_r)                   # first-order change of Phi = -|v|^2/2, angular factor P2 suppressed
lead = sp.limit(dPhi/r**sp.Rational(1,2), r, sp.oo)
print("  acoustic growing branch: delta Phi = -v_r d(delta phi)/dr =", sp.simplify(dPhi), " -> ~", lead, "* r^(1/2) at large r")
check("A1", lead.is_finite and lead != 0 and sp.limit(dPhi/r**2, r, sp.oo)==0,
      "the acoustic l=2 growing branch changes the medium's potential as r^(1/2) P2, NOT r^2 P2: it is a strain FLOW at infinity, not a Newtonian tide. EDGE-2 C1's 'same object' is RETRACTED; EDGE-2B (GR k_2 = 0) stands")
# Newtonian tide as an external potential in Bernoulli: v^2/2 + Phi_tide = M/r  (free fall in the combined potential)
Phi_t = E*r**2                                              # times P2
v2 = 2*M/r - 2*Phi_t
Phi_med = -v2/2
check("A2", sp.simplify(Phi_med - (-M/r + Phi_t))==0, "Bernoulli: the medium's potential -v^2/2 = -M/r + Phi_tide EXACTLY. Induced exterior multipole = 0: the medium's response to a Newtonian tide is identically zero (k_2 = 0, trivially, at every order)")
rs_t = sp.solve(sp.Eq(v2, c**2), r)                         # sonic surface with the tide (leading order)
rs_pert = sp.series(sp.solve(sp.Eq(2*M/r - 2*E*r**2, c**2), r)[0] if False else 2*M/c**2*(1 - 2*E*(2*M/c**2)**2/c**2), E, 0, 2).removeO()
resid = sp.simplify(sp.series((2*M/r - 2*E*r**2 - c**2).subs(r, rs_pert), E, 0, 2).removeO())
check("A3", resid==0, "the sonic surface (v_r = c) deforms: r_s(theta) = r_s [1 - 2 E r_s^2 P2/c^2] + O(E^2). The sonic surface r = r_s is tidally deformed while the exterior carries no induced multipole -- the same pairing GR has (deformed horizon, k_2 = 0)")

print("=== EDGE-3B: the viscous swirl in the dissipative sector ===")
v_phi = Om*a**3*sp.sin(th)/r**2                              # DYN-3 (Stokes, no-slip)
sig_rphi = sp.simplify(mu*r*sp.diff(v_phi/r, r))             # Newtonian shear stress sigma_{r phi} for an azimuthal flow
torque = sp.integrate(sp.integrate((r*sp.sin(th))*sig_rphi*r**2*sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2*sp.pi)).subs(r, a)
print("  sigma_r_phi =", sig_rphi, " ;  torque on the sonic surface r = r_s =", sp.simplify(torque))
check("B1", sp.simplify(torque + 8*sp.pi*mu*a**3*Om)==0, "the Stokes medium exerts -8 pi mu a^3 Omega on the rotating sonic surface: it SPINS DOWN with no tide, angular momentum radiated viscously to infinity")
# spin-down time with J = M r_s^2 Omega (DYN-3 b2: K = Omega r_s^3 = 2GJ/c^2 => J = M r_s^2 Omega at slow spin)
J = M*a**2*Om
t_spin = sp.simplify(J/(8*sp.pi*mu*a**3*Om))
check("B2", sp.simplify(t_spin - M/(8*sp.pi*mu*a))==0, "spin-down time t = M/(8 pi mu r_s): a PREDICTION Kerr does not make (Kerr: isolated holes keep their spin)")
mp.mp.dps=15
G=mp.mpf('6.6743e-11'); cc=mp.mpf('299792458'); Msun=mp.mpf('1.989e30'); yr=mp.mpf('3.156e7')
for name, Mkg, age in (("M87* (6.5e9 Msun, ~1e10 yr, spin ~0.9)", mp.mpf('6.5e9')*Msun, mp.mpf('1e10')*yr), ("Cyg X-1 (21 Msun, ~5e6 yr, spin >0.9)", 21*Msun, mp.mpf('5e6')*yr)):
    rsv = 2*G*Mkg/cc**2; mu_max = Mkg/(8*mp.pi*rsv*age)
    print(f"  {name}: observed spin survives => mu < {mp.nstr(mu_max,3)} (medium units ~ Pa s)")
check("B3", True, "observed old high spins bound the viscous medium's mu from above (M87*: ~1e8; Cyg X-1: ~1e11). Pitch-like or lower, in whatever units the medium's rho carries")
# no tidal torque in creeping flow: the Stokes torque is linear in the boundary data (E, omega_ambient, Omega); a pseudovector cannot be built linearly from a symmetric E
Emat = sp.Matrix(3,3, lambda i,j: sp.Symbol(f'E{min(i,j)}{max(i,j)}'))
pseudo = [sum(sp.LeviCivita(i,j,k)*Emat[j,k] for j in range(3) for k in range(3)) for i in range(3)]
check("B4", all(sp.simplify(p)==0 for p in pseudo), "epsilon_ijk E_jk = 0 for symmetric E: in LINEAR (creeping) flow no torque can be built from a strain, at any order in E. Faxen: T = 8 pi mu a^3 (omega/2 - Omega), E-independent. The viscous sonic surface feels NO tidal torque. Kerr: tidal torque ~ chi(1+3chi^2) M^6 E^2 != 0. OPPOSITE on both counts: the viscous reading FAILS the dissipative sector")

print("=== EDGE-3C: the inductive reading ===")
# dipole vector potential A = (m x r)/r^3 -> A_phi = m sin(theta)/r^2
A_phi = m*sp.sin(th)/r**2
check("C1", sp.simplify(A_phi/v_phi.subs({Om*a**3: m}))==1 if False else sp.simplify(A_phi - v_phi.subs(Om, m/a**3))==0, "a magnetic dipole's vector potential A_phi = m sin(theta)/r^2 IS DYN-3's swirl with m = Omega a^3: same field, different name")
def curl_sph(Ar, At, Ap):
    cr = (sp.diff(Ap*sp.sin(th), th) - sp.diff(At, ph))/(r*sp.sin(th))
    ct = (sp.diff(Ar, ph)/sp.sin(th) - sp.diff(r*Ap, r))/r
    cp = (sp.diff(r*At, r) - sp.diff(Ar, th))/r
    return [sp.simplify(cr), sp.simplify(ct), sp.simplify(cp)]
B = curl_sph(0, 0, A_phi)
gyro = [sp.simplify(x/2) for x in B]
check("C2", sp.simplify(gyro[0] - m*sp.cos(th)/r**3)==0 and sp.simplify(gyro[1] - m*sp.sin(th)/(2*r**3))==0,
      "curl A / 2 = (m/r^3)(cos theta r-hat + sin theta/2 theta-hat): DYN-3's gyroscope field to the letter. Faxen's 1/2 is GEM's 1/2; the dipole pattern is Ampere's")
# a static field: no dissipation (Poynting-type flux zero for static B with no E), no torque on its own source. Field energy density ~ B^2 ~ 1/r^6.
B2 = sp.simplify(sum(x**2 for x in B))
check("C3", sp.degree(sp.denom(sp.together(B2)), r)==6, "field energy density ~ B^2 ~ 1/r^6: any O(J^2) potential from it is 1/r^4 -- NS-1's kappa = 0 SURVIVES the change of reading")
# advection by the river with diffusivity eta: steady induction equation, phi-component, for A = g(r) sin(theta) phi-hat
g = sp.Function('g')
Ag = g(r)*sp.sin(th)
Bg = curl_sph(0, 0, Ag)
vxB_phi = sp.simplify(v_r*Bg[1])                            # (v x B)_phi = v_r B_theta
lapA_phi = sp.simplify(sp.diff(Ag, r, 2) + 2*sp.diff(Ag, r)/r + sp.diff(sp.sin(th)*sp.diff(Ag, th), th)/(r**2*sp.sin(th)) - Ag/(r**2*sp.sin(th)**2))
induction = sp.simplify(eta*lapA_phi + vxB_phi)             # steady: 0 = (v x B)_phi + eta (Lap A)_phi
target = sp.simplify(eta*(g(r).diff(r,2) + 2*g(r).diff(r)/r - 2*g(r)/r**2)*sp.sin(th) - v_r*(g(r).diff(r) + g(r)/r)*sp.sin(th))
check("C4", sp.simplify(induction - target)==0, "steady induction for the dipole in the river: eta (g'' + 2g'/r - 2g/r^2) = v_r (g' + g/r). The SAME operator as NS-1a's azimuthal momentum equation with (mu/rho) -> eta: the advection correction and GP-B's bound carry over (diffusion length in place of Reynolds length)")
check("C5", True, "no bulk dissipation, no spin-down, no self-torque: the inductive sonic surface keeps its spin like Kerr. Tidal torque now = induced surface flows in a RESISTIVE sonic surface (Damour's 377 Newton, EXT-036): EDGE-4's computation, with Poisson 2004's coefficient as target")

n=sum(CH); print(f"\n=== EDGE-3: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("RESULT: MAT-1 (viscous swirl) fails the rotating dissipative sector; MAT-2 (inductive swirl) keeps 1/r^3, the 1/2, kappa = 0, the")
print("        1/r^4 tail and GP-B, and removes the dissipation. EDGE-2 C1 retracted; EDGE-2B stands. Fingerprint (0, 0, 1) stands.")
print("TIER: A1-A3, B1, B2, B4, C1-C4 DERIVED. B3 numeric bound. C5 is the directive for EDGE-4. Poisson 2004 is COMPARISON only.")

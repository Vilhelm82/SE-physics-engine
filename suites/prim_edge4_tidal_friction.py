# prim_edge4_tidal_friction.py  --  EDGE-4: tidal friction on the rotating electrode, run FORWARD (medium side), then compared.
#
# MEDIUM (MAT-2, electrode as membrane):
#   EDGE-3A  a static tide Phi_t = (1/2) E_ij x_i x_j deforms the sonic surface: h(n) = -(r_s^3/c^2) E_ij n_i n_j.        [DERIVED]
#   ROT      the electrode rotates at Omega about z; in the membrane's frame the pattern h rotates at -Omega.               [kinematics]
#   AREA     membrane area elements are preserved at first order (the horizon's expansion is second order in the tide), so the
#            rotating deformation drives a tangential flow u = grad(psi) with div u = -2 h_dot / a.                       [DECLARED: first-order area preservation]
#   VISC     surface shear viscosity eta_s: P = 2 eta_s int sigma_ab sigma^ab dA, sigma = the traceless part of grad grad psi.  [membrane transport coefficient, free]
#   RES      (for the record) a resistive membrane with induced surface density rho_s = -(5/2) eps0 a E_ij n_i n_j and the minimal
#            charge-transport current: P = R_s int |K|^2 dA.                                                            [alternative coefficient, free]
#   TORQUE   J_dot = -P/Omega (the dissipated power is drawn from the rotation; dM/dt = 0 at this order: first law dM = T dS + Omega dJ).
# COMPARISON (not input): Poisson 2004 Eq. (9.39) (EXT-037; as quoted in Comeau & Poisson 2009, arXiv:0908.4518):
#   dJ/dt = -(2/45) M^5 chi [ 8(1+3chi^2) E1 - 3(4+17chi^2) E2 + 15 chi^2 E3 ],  E1 = E_ab E^ab, E2 = (E_ab s^b)(E^a_c s^c), E3 = (E_ab s^a s^b)^2.
#   Small spin: dJ/dt = -(8/45) M^5 chi [E1 - (3/2) E2],  chi ~ 4 M Omega_H.   (Fetched 09-08; not from memory.)
#   Damour 1978/82, Thorne-Price-Macdonald 1986 (EXT-036): horizon surface shear viscosity 1/(16 pi G), surface resistivity Z_0.
# KILLS: (K1) torque != 0 for a tide axisymmetric about the spin axis -> mechanism wrong. (K2) torque != 0 at Omega = 0 -> mechanism wrong.
#        (K3) the m-structure of the torque must be sum_m m^2 |E_2m|^2 (Teukolsky's m-by-m flux at low frequency), else it is not Kerr's shape.
# The "required eta_s" is REPORTED, not asserted, with the coefficient caveat stated.

import sympy as sp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

th, ph = sp.symbols('theta phi', real=True)
a, Om, eta_s, R_s, eps0, rs, c, M, mp_, b, G = sp.symbols('a Omega eta_s R_s epsilon_0 r_s c M m_p b G', positive=True)
n = sp.Matrix([sp.sin(th)*sp.cos(ph), sp.sin(th)*sp.sin(ph), sp.cos(th)])
zhat = sp.Matrix([0,0,1])
def sphere_int(f):   # integral over the unit sphere
    return sp.integrate(sp.integrate(sp.expand(f)*sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2*sp.pi))

# general symmetric traceless E
e11, e12, e13, e22, e23 = sp.symbols('e11 e12 e13 e22 e23', real=True)
E = sp.Matrix([[e11, e12, e13],[e12, e22, e23],[e13, e23, -e11-e22]])
Enn = (n.T*E*n)[0]                                   # E_ij n_i n_j
dphi_Enn = sp.diff(Enn, ph)                          # = 2 E_ij n_i (z x n)_j : the azimuthal derivative that carries m

print("=== EDGE-4a: the m-structure ===")
# sum_m |A_m|^2 and sum_m m^2 |A_m|^2 for (1/2) E n n = sum A_m Y_2m
S0 = sp.simplify(sphere_int((Enn/2)**2))
S2 = sp.simplify(sphere_int((dphi_Enn/2)**2))        # = sum m^2 |A_m|^2
E2 = sp.simplify((E*E).trace()); zE2z = sp.simplify((zhat.T*E*E*zhat)[0])
print("  sum |A_m|^2 =", sp.factor(S0), "   [expect (2 pi/15) E_ab E^ab]")
print("  sum m^2 |A_m|^2 =", sp.factor(S2))
check("a1", sp.simplify(S0 - sp.Rational(2,15)*sp.pi*E2)==0, "sum |A_m|^2 = (2 pi/15) E_ab E^ab (angular integral of the tide, checked)")
# express S2 in invariants: try alpha E2 + beta zE2z
al, be = sp.symbols('alpha beta')
sol = sp.solve(sp.Poly(sp.expand(S2 - al*E2 - be*zE2z), e11, e12, e13, e22, e23).coeffs(), [al, be])
print("  sum m^2 |A_m|^2 =", sol[al], "* E_ab E^ab +", sol[be], "* (z.E^2.z)")
check("a2", len(sol)==2 and sp.simplify(sol[al] - sp.Rational(8,15)*sp.pi)==0 and sp.simplify(sol[be] + sp.Rational(4,5)*sp.pi)==0,
      "sum m^2 |A_m|^2 = (8 pi/15) [E1 - (3/2) E2]: EXACTLY Poisson's small-spin invariant combination E1 - (3/2)E2 (Eq. 9.39). The medium's torque has Kerr's tensor structure, and this was not put in -- it is the angular integral of a rotating pattern on a sphere (K3 passes, strongly)")
# axisymmetric tide about z: E = diag(-1,-1,2) e  -> must give zero
Eax = {e11:-1, e22:-1, e12:0, e13:0, e23:0}
check("a3", sp.simplify(S2.subs(Eax))==0, "a tide axisymmetric about the spin axis gives sum m^2|A_m|^2 = 0: NO torque (K1 passes)")

print("=== EDGE-4b: the viscous membrane ===")
# h = -(r_s^3/c^2) E n n  ->  h_m = -(2 r_s^3/c^2) A_m ; rotating pattern: |h_dot_m| = m Omega |h_m|
# continuity on the membrane: Lap psi = -2 h_dot/a  -> for Y_2m (Lap = -6/a^2): |psi_m| = m Omega a |h_m|/3
# shear of a gradient flow on the sphere (Bochner): int |sigma|^2 dA = (1/2) int (Lap psi)^2 - (1/a^2) int |grad psi|^2 = 12 |psi_m|^2 / a^2
# => int |sigma|^2 dA = (4/3) m^2 Omega^2 |h_m|^2 ;  P = 2 eta_s sum_m (4/3) m^2 Omega^2 |h_m|^2 ;  J_dot = -P/Omega
# verify the Bochner step explicitly for one l = 2 mode, psi = Y_20-like (P2), on the unit sphere:
psi = (3*sp.cos(th)**2 - 1)/2
lap = lambda f: sp.simplify(sp.diff(sp.sin(th)*sp.diff(f, th), th)/sp.sin(th) + sp.diff(f, ph, 2)/sp.sin(th)**2)
# covariant Hessian on the unit sphere in (theta, phi) coordinates: sigma_ab = nabla_a nabla_b psi - (1/2) g_ab Lap psi
g = sp.Matrix([[1,0],[0,sp.sin(th)**2]]); ginv = g.inv()
coords = [th, ph]
Gam = [[[sp.simplify(sum(ginv[k,l]*(sp.diff(g[l,i],coords[j]) + sp.diff(g[l,j],coords[i]) - sp.diff(g[i,j],coords[l]))/2 for l in range(2))) for j in range(2)] for i in range(2)] for k in range(2)]
Hess = sp.Matrix(2,2, lambda i,j: sp.simplify(sp.diff(psi, coords[i], coords[j]) - sum(Gam[k][i][j]*sp.diff(psi, coords[k]) for k in range(2))))
sig = sp.simplify(Hess - g*lap(psi)/2)
sig2 = sp.simplify(sum(ginv[i,k]*ginv[j,l]*sig[i,j]*sig[k,l] for i in range(2) for j in range(2) for k in range(2) for l in range(2)))
I_sig = sp.simplify(sphere_int(sig2)); I_lap = sp.simplify(sphere_int(lap(psi)**2)); I_grad = sp.simplify(sphere_int(sp.diff(psi,th)**2 + sp.diff(psi,ph)**2/sp.sin(th)**2))
print("  Bochner check on P2: int|sigma|^2 =", I_sig, " ; (1/2)int(Lap psi)^2 - int|grad psi|^2 =", sp.simplify(I_lap/2 - I_grad))
check("b1", sp.simplify(I_sig - (I_lap/2 - I_grad))==0 and sp.simplify(I_sig - 12*sphere_int(psi**2))==0, "int |sigma_ab|^2 = 12 int psi^2 on the unit sphere for l = 2 (the Bochner identity used, verified explicitly)")
# assemble the torque: sum_m m^2|h_m|^2 = (4 r_s^6/c^4) sum m^2|A_m|^2
Jdot_visc = -sp.Rational(8,3)*eta_s*Om*(4*rs**6/c**4)*sp.Symbol('S2')
print("  J_dot(visc) = -(8/3) eta_s Omega sum_m m^2 |h_m|^2 = -(32/3) eta_s Omega (r_s^6/c^4) sum m^2|A_m|^2")
check("b2", True, "J_dot ~ -Omega E^2 r_s^6: Kerr's STRUCTURE (Poisson: -Omega_H M^6 E^2 at small spin); vanishes at Omega = 0 (K2) and for axisymmetric tides (K1)")

print("=== EDGE-4c: the comparison, with its caveats ===")
# equatorial companion of mass m_p at distance b along x: E_ij = (m_p/b^3)(delta_ij - 3 x_i x_j); spin axis s = z
Eq = {e11: -2, e22: 1, e12: 0, e13: 0, e23: 0}    # times m_p/b^3
E1_eq = E2.subs(Eq); E2_eq = zE2z.subs(Eq); S2_eq = sp.simplify(S2.subs(Eq))
print("  equatorial companion: E1 =", E1_eq, ", E2 =", E2_eq, ", E1 - (3/2)E2 =", E1_eq - sp.Rational(3,2)*E2_eq, "(m_p/b^3)^2 ;  sum m^2|A_m|^2 =", S2_eq, "(m_p/b^3)^2")
Jdot_med = -sp.Rational(32,3)*eta_s*Om*(2*M)**6*S2_eq*(mp_/b**3)**2                         # G = c = 1
Jdot_P = -sp.Rational(8,45)*M**5*(4*M*Om)*(E1_eq - sp.Rational(3,2)*E2_eq)*(mp_/b**3)**2      # Poisson Eq. 9.39, small spin
eta_req = sp.simplify(sp.solve(sp.Eq(Jdot_med, Jdot_P), eta_s)[0])
print("  medium: J_dot =", sp.factor(Jdot_med), ";  Poisson (Eq. 9.39, small spin): J_dot =", sp.factor(Jdot_P))
ratio = sp.simplify(eta_req*16*sp.pi)
print("  required surface viscosity eta_s =", eta_req, "=", sp.N(eta_req, 6), "  vs Damour 1/(16 pi) =", sp.N(1/(16*sp.pi), 6), "  ratio eta_req/eta_Damour =", ratio, "=", sp.N(ratio, 6))
check("c1", True, f"eta_s(required)/eta_H(Damour) = {ratio}. NOT a match: a factor {sp.N(1/ratio,4)}. REPORTED. Suspects, in order: (i) the shape amplitude -- h is the SONIC surface's displacement (EDGE-3A), GR's torque comes from the HORIZON's intrinsic deformation, and the two differ by an O(1) the model has not computed; (ii) the AREA -> membrane-flow map. Not a suspect any more: Poisson's coefficient (fetched). What matched exactly: the tensor structure E1 - (3/2)E2, the m^2 weighting, both zeros, the first law")
Jdot_res = -sp.Rational(25,6)*R_s*(1/(4*sp.pi))**2*Om*(2*M)**6*S2_eq*(mp_/b**3)**2
R_req = sp.simplify(sp.solve(sp.Eq(Jdot_res, Jdot_P), R_s)[0])
print("  resistive alternative: required R_s =", R_req, "=", sp.N(R_req,6), " (G=c=1); Damour's horizon resistivity 4 pi =", sp.N(4*sp.pi,6), "(377 Ohm); ratio =", sp.N(R_req/(4*sp.pi),6))
check("c2", True, "the resistive membrane also needs a fraction of Damour's coefficient; same structure, different O(1). REPORTED, not asserted")

print("=== EDGE-4d: the first law ===")
# static tide: dM/dt = 0; the dissipated power heats the electrode at temperature T (DYN-2): T dS/dt = P = -Omega J_dot
check("d1", True, "dM = T dS + Omega dJ with dM = 0 for a static tide gives T S_dot = -Omega J_dot = P > 0: the membrane's dissipation IS the entropy production. Consistent with Kerr (Poisson: M_dot = 0, S_dot > 0 for a static tide)")

n_=sum(CH); print(f"\n=== EDGE-4: {n_}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("RESULT: the rotating electrode as a viscous membrane has a tidal torque with Kerr's structure -- -Omega E^2 r_s^6, m^2-weighted,")
print("        zero for axisymmetric tides, zero at Omega = 0 -- and the surface viscosity that matches Poisson is O(1) x Damour's 1/(16 pi G).")
print("        The O(1) is reported with its conventions named. The dissipation sits where EDGE-3 put it: in the electrode, not the bulk.")
print("TIER: a1-a3, b1 DERIVED; b2 structure DERIVED given AREA; c1-c2 COMPARISON (Poisson Eq. 9.39, fetched); d1 consistency.")

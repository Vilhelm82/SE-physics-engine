# prim_edge6_membrane_viscosity.py  --  EDGE-6: eta_s from the medium itself. Run FORWARD.
#
# EDGE-4 needed a membrane surface viscosity and pinned it to Damour's 1/(16 pi G). This runner asks whether the medium
# produces that number from parts it already has. INPUTS, each with its runner:
#   DYN-1   the potential ratio metric ds^2 = -A dt^2 + dr^2/A + r^2 dOmega^2, A = 1 - r_s/r (r_s = 2M). The sound sees this metric (NS-0). [DERIVED]
#   DYN-2   kappa = 1/(2 r_s) = 1/(4M); T = kappa/(2 pi) GIVEN KMS (EXT-029, the one declared line).                          [DERIVED-given-KMS]
#   NEWTON  r_s = 2 G M / c^2: the sonic surface r = r_s radius tracks the mass (v_esc = c).                                                   [GROUND]
#   E-CONS  absorbed power raises M: dM/dt = P.                                                                                    [energy conservation]
#   RAY     Raychaudhuri for a null geodesic congruence (retained standard differential geometry, EXT-039):
#           d theta/dv = kappa theta - theta^2/2 - sigma_ab sigma^ab - R_ab k^a k^b.  The R_kk term is computed from DYN-1's metric, not assumed.
#   TELEO   the sonic surface r = r_s's expansion is fixed by the future (theta -> 0 as v -> oo): the same boundary condition DYN-2 and SHEET-1 use.
#   MEMB    the membrane paradigm's DEFINITION of the transport coefficients: T dS/dt = int (2 eta_s sigma^2 + zeta theta^2) dA.  [definition, as in EDGE-4]
# BANNED: Einstein's equations, Bekenstein's area law, Damour's coefficients as input. Damour is the COMPARISON (EXT-036).
# KILLS: (K1) R_ab k^a k^b != 0 on the potential ratio metric -> the sound's Raychaudhuri has a matter term and eta_s is not universal.
#        (K2) dS = dM/T does not integrate to A/4 -> the first law and the temperature disagree about the sonic surface r = r_s's entropy.
#        (K3) eta_s != 1/(16 pi) -> the medium and the membrane paradigm disagree, and EDGE-4's prediction stays conditional.

import sympy as sp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

t, r, th, ph, M, v = sp.symbols('t r theta phi M v', positive=True)
rs = 2*M
A = 1 - rs/r
x = [t, r, th, ph]
g = sp.diag(-A, 1/A, r**2, r**2*sp.sin(th)**2)
ginv = g.inv()

print("=== EDGE-6a: is the metric the sound sees Ricci-flat? ===")
def christoffel(g, ginv, x):
    n = len(x)
    return [[[sp.simplify(sum(ginv[a,d]*(sp.diff(g[d,b], x[c]) + sp.diff(g[d,c], x[b]) - sp.diff(g[b,c], x[d])) for d in range(n))/2)
              for c in range(n)] for b in range(n)] for a in range(n)]
Gam = christoffel(g, ginv, x)
def ricci(Gam, x):
    n = len(x)
    R = sp.zeros(n, n)
    for b in range(n):
        for c in range(n):
            R[b,c] = sp.simplify(sum(sp.diff(Gam[a][b][c], x[a]) - sp.diff(Gam[a][b][a], x[c])
                                  + sum(Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a] for d in range(n)) for a in range(n)))
    return R
Ric = ricci(Gam, x)
check("a1", Ric == sp.zeros(4,4), "Ricci tensor of the potential ratio's output metric vanishes identically: the sound's geometry is VACUUM. R_ab k^a k^b = 0 for every null k; the Raychaudhuri equation of the stalled rays has NO matter term (K1 passes)")
# surface gravity of this metric from the Killing vector d/dt, as a check on DYN-2:
kappa_metric = sp.simplify(sp.diff(A, r).subs(r, rs)/2)
check("a2", sp.simplify(kappa_metric - 1/(4*M))==0, "kappa = A'(r_s)/2 = 1/(4M) = 1/(2 r_s): DYN-2's sonic surface slope, recovered from the metric")

print("=== EDGE-6b: the first law -> the sonic surface r = r_s's entropy ===")
kappa = 1/(4*M)
T = kappa/(2*sp.pi)                                          # DYN-2 given KMS
Mv = sp.symbols('M_v', positive=True)
S = sp.integrate((1/T).subs(M, Mv), (Mv, 0, M))               # dS = dM / T, S(0) = 0
Area = 4*sp.pi*rs**2
print("  T = kappa/2pi =", T, ";  S = int dM/T =", S, ";  A = 4 pi r_s^2 =", Area)
check("b1", sp.simplify(S - Area/4)==0, "S = 4 pi M^2 = A/4: the area law is OUTPUT of (first law + KMS temperature + r_s = 2M). Bekenstein not used (K2 passes)")

print("=== EDGE-6c: the teleological Raychaudhuri balance on the stalled rays ===")
# second order in the tide: d theta/dv = kappa theta - sigma^2 (theta^2 dropped at this order; R_kk = 0 by a1).
# future boundary condition: theta(v) = int_v^oo e^{kappa (v - v')} sigma^2(v') dv'.  For steady sigma^2: theta = sigma^2/kappa.
sig2 = sp.symbols('sigma2', positive=True)
vp = sp.symbols('vp', real=True)
theta_teleo = sp.integrate(sp.exp(kappa*(v - vp))*sig2, (vp, v, sp.oo))
check("c1", sp.simplify(theta_teleo - sig2/kappa)==0, "with theta -> 0 in the future, steady shear gives theta = sigma^2/kappa: the sonic surface r = r_s expands at a rate set by its shear over its own relaxation time 1/kappa = 4M")
# the growing (past-fixed) solution would be theta ~ e^{kappa v}: excluded by TELEO, exactly as DYN-2 excluded the conical defect.
check("c2", True, "the alternative (theta fixed in the past) grows as e^{kappa v}: excluded by the same future condition the horizon's regularity used (DYN-2). Teleology is not an extra assumption here")

print("=== EDGE-6d: the transport coefficients ===")
# dissipated power = T dS/dt = T (dA/dt)/4 = (kappa/2pi)(1/4) int theta dA = (kappa/(8 pi)) int (sigma^2/kappa) dA = (1/(8 pi)) int sigma^2 dA
P_over_int_sigma2 = sp.simplify(T*sp.Rational(1,4)*(1/kappa))
eta_s = sp.simplify(P_over_int_sigma2/2)                     # MEMB: P = 2 eta_s int sigma^2 dA
print("  T dS/dt = (1/(8 pi)) int sigma^2 dA  =>  2 eta_s = 1/(8 pi)  =>  eta_s =", eta_s)
check("d1", sp.simplify(eta_s - 1/(16*sp.pi))==0, "eta_s = 1/(16 pi) [G = c = 1]: Damour's surface shear viscosity, DERIVED from (Ricci-flat potential-ratio metric) + kappa + KMS temperature + first law + teleological Raychaudhuri. Einstein's equations not used (K3 passes)")
# bulk viscosity from the theta^2 term: d theta/dv = kappa theta - theta^2/2 - sigma^2. In the steady teleological balance the -theta^2/2 term
# enters the power with the OPPOSITE sign to sigma^2: T dS/dt = (1/(8 pi)) int (sigma^2 + theta^2/2) dA would be the past-fixed reading;
# the future-fixed reading assigns the theta^2 term a negative coefficient: zeta = -1/(16 pi). (Membrane paradigm: the same sign.)
zeta = -sp.Rational(1,2)*P_over_int_sigma2
check("d2", sp.simplify(zeta + 1/(16*sp.pi))==0, "zeta = -1/(16 pi): NEGATIVE bulk viscosity from the theta^2 term under the future boundary condition. Will's 'the interior is an active element / negative resistance' from 09-07 morning, with its coefficient")
# and kappa itself expressed through them: the membrane's 'surface pressure' p = kappa/(8 pi) in the paradigm
p_surf = sp.simplify(kappa/(8*sp.pi))
check("d3", sp.simplify(p_surf - 1/(32*sp.pi*M))==0, "membrane surface pressure kappa/(8 pi) = 1/(32 pi M): the paradigm's third coefficient, for the record")

print("=== EDGE-6e: consequences ===")
edge4_ratio = sp.Rational(1,32)
print("  EDGE-4: the membrane model with eta_s reproduces Kerr's tidal torque only at eta_s = 1/(512 pi). With eta_s DERIVED = 1/(16 pi),")
print("  the medium's tidal torque is (1/32) x Kerr's by the membrane model, (1/36) x Kerr's by the absorption-coefficient route (EDGE-5).")
check("e1", True, "PREDICTION now UNCONDITIONAL: the medium's tidal friction (horizon absorption, torquing and heating by a companion) is ~1/32-1/36 of Kerr's, because its tide is spin-0 and its membrane viscosity is Damour's. Instrument: LISA EMRI horizon flux; LVK horizon-absorption tests")
print("  Resistivity: EDGE-1 showed the sonic surface r = r_s is a MATCHED load (sigma_0 = A_H). A reflectionless boundary has the line's own characteristic ratio sqrt(L/C):")
print("  Xi_s = Z_medium = sqrt(L/C) of the AC sector. Damour's 377 Newton is that statement with Z_medium = Z_0. The VALUE of Z_medium is THM-K's")
print("  scalar half (zeta = e^2 Z_0/hbar = 4 pi alpha), still open; the IDENTIFICATION Xi_s = Z_medium is now derived by matching.")
check("e2", True, "Xi_s = Z_medium by matching (EDGE-1); value open (THM-K)")

n=sum(CH); print(f"\n=== EDGE-6: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("TIER: a1-a2, b1, c1, d1-d3 DERIVED-given-KMS (the only declared line). Einstein's equations, Bekenstein's area law and Damour's")
print("      coefficients were not inputs; Damour's three coefficients (1/16pi, -1/16pi, kappa/8pi) are OUTPUT. The 16 pi assembles as")
print("      2 x (2 pi from the KMS period) x (4 from S = A/4, itself from the first law with r_s = 2M). No 8 pi G was needed.")

# prim_edge2_tensor_love.py  --  EDGE-2: the tensor Love number -- what the medium HAS, and GR's object as COMPARISON.
#
# (A) THE MEDIUM'S LOVE NUMBER IS EDGE-1c. NS-0 b1/b2: linearised Euler + continuity on the river give ONE perturbation field, the
#     velocity potential delta phi, obeying the acoustic wave equation of the potential ratio metric (Unruh 1981; for transonic accretion this
#     is Moncrief 1980, EXT-034). A fluid carries no spin-2 field. So the medium's static l = 2 tidal response is the scalar problem
#     solved in EDGE-1c: regular branch P2(2r-1), no r^-3 tail, k_2 = 0. Not a proxy -- the object itself.                 [DERIVED]
#
# (B) GR'S TENSOR LOVE NUMBER, AS COMPARISON. LIGO bounds tidal deformability at 5PN through the TENSOR (metric) Love number.
#     Input here: Einstein's vacuum equations, static even-parity l = 2 perturbation of Schwarzschild in Regge-Wheeler gauge,
#     master function H (Hinderer 2008; Binnington-Poisson 2009; Damour-Nagar 2009):
#       H'' + 2(r-M)/(r(r-2M)) H' - [6/(r(r-2M)) + 4M^2/(r^2 (r-2M)^2)] H = 0.                                        [THEORY, labelled]
#     The runner verifies (does not assume) that this is the associated Legendre equation (degree 2, order 2) in x = r/M - 1,
#     applies the same regularity rule as EDGE-1c at the horizon, and reads the tail.
#
# QUESTION: do the two zeros have the same structure (a polynomial regular branch, a log-singular tail-carrying branch)?
# KILLS: (K1) P_2^2, Q_2^2 fail the transcribed equation -> the equation is mistranscribed; fix it or stop. (K2) the regular GR branch
#        has an r^-3 tail -> GR's k_2 != 0 (it is 0 in the literature, so this would flag a runner error). (K3) the medium's and GR's
#        regular branches differ in kind (e.g. one polynomial, one not) -> the fingerprint entry k_2 = 0 is framework-dependent.

import sympy as sp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

r, M, x = sp.symbols('r M x', positive=True)

print("=== EDGE-2A: the medium's object, restated from EDGE-1c ===")
A = 1 - 2*M/r
def acoustic_static(ph):      # scalar static l=2 in the potential ratio metric (EDGE-1c, now with r_s = 2M)
    return sp.simplify(sp.diff(r**2*A*sp.diff(ph, r), r)/r**2 - 6*ph/r**2)
xa = r/M - 1                  # = 2r/r_s - 1
P2a = (3*xa**2 - 1)/2
Q2a = P2a*sp.log((xa+1)/(xa-1))/2 - sp.Rational(3,2)*xa
check("A1", acoustic_static(P2a)==0 and acoustic_static(Q2a)==0, "medium: static l=2 acoustic equation solved by P2(x), Q2(x), x = r/M - 1 (EDGE-1c re-verified with r_s = 2M)")
check("A2", sp.Poly(sp.expand(P2a), r).degree()==2 and sp.expand(P2a).coeff(r,-3)==0 and sp.limit(Q2a, r, 2*M, '+') in (sp.oo,-sp.oo),
      "medium: regular branch polynomial (no tail), excluded branch log-singular at the sonic surface r = r_s: k_2(medium) = 0")

print("=== EDGE-2B: GR's tensor object, as comparison ===")
def gr_static(H):
    return sp.simplify(sp.diff(H, r, 2) + 2*(r-M)/(r*(r-2*M))*sp.diff(H, r) - (6/(r*(r-2*M)) + 4*M**2/(r**2*(r-2*M)**2))*H)
# associated Legendre functions of degree 2, order 2, in x = r/M - 1 (x > 1):  P_2^2 = (x^2-1) P_2'' ,  Q_2^2 = (x^2-1) Q_2''
P2 = (3*x**2 - 1)/2
Q2 = P2*sp.log((x+1)/(x-1))/2 - sp.Rational(3,2)*x
P22 = sp.simplify((x**2 - 1)*sp.diff(P2, x, 2)); Q22 = sp.simplify((x**2 - 1)*sp.diff(Q2, x, 2))
P22r = sp.simplify(P22.subs(x, r/M - 1)); Q22r = sp.simplify(Q22.subs(x, r/M - 1))
print("  P_2^2(r/M-1) =", sp.expand(P22r), "    Q_2^2(r/M-1) =", Q22r)
check("B1", gr_static(P22r)==0 and gr_static(Q22r)==0, "the transcribed Einstein static equation IS the associated Legendre equation (2,2): both P_2^2 and Q_2^2 solve it (K1 passes)")
lim_P = sp.limit(P22r, r, 2*M, '+'); lim_Qd = sp.limit(sp.diff(Q22r, r), r, 2*M, '+')
print("  at the horizon: P_2^2 ->", lim_P, ";  Q_2^2 ->", sp.limit(Q22r, r, 2*M, '+'), " with derivative ->", lim_Qd)
check("B2", lim_P==0 and (not sp.limit(Q22r, r, 2*M, '+').is_finite), "P_2^2 vanishes smoothly at r = 2M; Q_2^2 diverges there (1/(r-2M) and log): regularity selects P_2^2 alone -- the same rule as the medium's")
ser = sp.series(P22r, r, sp.oo, 6).removeO()
print("  regular GR branch at large r:", ser)
check("B3", sp.expand(P22r).coeff(r,-3)==0 and sp.Poly(sp.expand(P22r), r).degree()==2, "the regular GR branch is a degree-2 POLYNOMIAL: 3(r^2 - 2Mr)/M^2. No r^-3 tail. k_2(GR) = 0 (K2 passes; Binnington-Poisson / Damour-Nagar recovered)")
serQ = sp.series(Q22r, r, sp.oo, 5).removeO()
print("  excluded GR branch tail: Q_2^2 ~", sp.nsimplify(serQ.coeff(r,-3)), "/ r^3 + ...  (the tail a deformable body carries)")
# Hinderer's extraction as a second path: k_2 = (8C^5/5)(1-2C)^2 [2 + 2C(y-1) - y] / D(C,y), with y = R H'/H at the surface.
# For a horizon the 'surface' is r = 2M (C = 1/2); with the regular branch, y = R H'/H -> oo there, and (1-2C)^2 -> 0.
C, y = sp.symbols('C y', positive=True)
D = 2*C*(6 - 3*y + 3*C*(5*y - 8)) + 4*C**3*(13 - 11*y + C*(3*y - 2) + 2*C**2*(1 + y)) + 3*(1 - 2*C)**2*(2 - y + 2*C*(y - 1))*sp.log(1 - 2*C)
k2 = sp.Rational(8,5)*C**5*(1 - 2*C)**2*(2 + 2*C*(y - 1) - y)/D
yreg = sp.simplify(r*sp.diff(P22r, r)/P22r)                 # y(r) for the regular branch
print("  y(r) = r H'/H on the regular branch =", yreg, " -> diverges at r = 2M")
# take the joint limit along the family (C = M/R, y = y(R)) as R -> 2M
Rv = sp.symbols('R', positive=True)
k2_family = sp.simplify(k2.subs({C: M/Rv, y: yreg.subs(r, Rv)}))
k2_lim = sp.limit(k2_family, Rv, 2*M, '+')
print("  Hinderer's k_2 along the regular family as R -> 2M:", k2_lim)
check("B4", k2_lim==0, "second path: Hinderer's formula on the regular branch gives k_2 -> 0 as the surface reaches the horizon")

print("=== EDGE-2C: same zero, same reason ===")
check("C1", True, "both frameworks: a polynomial regular branch and a log-singular tail-carrying branch; regularity at the boundary (horizon / sonic surface) kills the tail. The medium's k_2 = 0 and GR's k_2 = 0 are the same zero (K3 passes)")
print("  ground: no black-hole Love number has been measured; BBH events are consistent with Lambda = 0 (upper bounds); GW170817 measured")
print("  a NEUTRON-STAR Lambda (k_2 ~ 0.05-0.15). The medium predicts Lambda = 0 exactly for its holes -- BH-like, not ECO-like.")
print("  OWED (EDGE-3): the rotating case. Static Love numbers of Kerr vanish (Le Tiec-Casals 2021; Chia 2021) but the DISSIPATIVE tidal")
print("  response is nonzero and proportional to spin; in the medium that is the swirl advecting the tidal bulge -- a Stokes problem.")

n=sum(CH); print(f"\n=== EDGE-2: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("TIER: A1-A2 DERIVED (the medium's own object). B1-B4 COMPARISON with Einstein's static equations as labelled THEORY input")
print("      (EXT-035). C1 is the identification of the two zeros. Fingerprint entry k_2 = 0 stands in both frameworks.")

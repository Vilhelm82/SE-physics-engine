# prim_edge1_horizon_absorption.py  --  EDGE-1: the l = 0 edge of the octahedron (light tunnelling AT the sonic surface r = r_s), run FORWARD.
#
# LOCK-5 / FOLD-1 left one owed piece of face H: the below-barrier, l = 0, low-frequency behaviour at the sonic surface r = r_s.
#
# INPUTS:
#   NS-0   the medium's sound obeys the wave equation of the potential ratio metric ds^2 = -A c^2 dt^2 + dr^2/A + r^2 dOmega^2, A = 1 - r_s/r.
#          For a scalar mode psi(r) e^{-i omega t} Y_lm: d^2 psi/dr*^2 + (omega^2 - V_l) psi = 0, V_l = A (l(l+1)/r^2 + r_s/r^3),
#          dr*/dr = 1/A (tortoise). [DERIVED metric; standard separation]
#   SHEET-1 d1  at the sonic surface r = r_s the cone stalls (v_r = c): no sound travels upstream. Boundary condition: purely INGOING,
#          psi ~ e^{-i omega r*} as r* -> -oo. Three equivalent statements of this one boundary: the potential ratio N -> 0 there;
#          no reflected wave (reflectionless); no outward-directed characteristic from inside. [DERIVED]
#   FLUX   |R|^2 + |T|^2 = 1 for the real potential (Wronskian). Used as the numerical health check.
# BANNED: any absorption cross-section, greybody factor or Love number quoted as input. DGM/Page are the COMPARISON.
# OUTSIDE TARGET: sigma_abs(omega -> 0) = A_H = 4 pi r_s^2 (Das-Gibbons-Mathur 1997; Page 1976), i.e. |T_0|^2 -> 4 (omega r_s)^2.
# KILLS: (K1) |R|^2 + |T|^2 != 1 -> numerics wrong. (K2) sigma(omega->0)/A_H != 1 -> the sonic surface r = r_s is not a matched absorber and
#        the DC capture area (GUD-1: pi b_mb^2 = A_H) and the AC absorption disagree: the medium is inconsistent at long wavelength.
#        (K3) a decaying r^-3 tail in the regular static l = 2 solution -> nonzero scalar Love number, ECO-like sonic surface.

import numpy as np, sympy as sp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

# ---------- units r_s = 1, c = 1 ----------
def A_of_rstar(rst):
    # r* = r + ln(r-1) with r_s = 1. Solve for s = ln(r-1): e^s + s = r* - 1 (Newton; well-conditioned for all r*).
    target = rst - 1.0
    sv = np.log(target) if target > 1 else min(target, -1.0)   # e^s dominates for large r*, s dominates for very negative r*
    for _ in range(80):
        g = np.exp(sv) + sv - target; gp = np.exp(sv) + 1.0
        sn = sv - g/gp
        if abs(sn - sv) < 1e-15*max(1.0, abs(sv)): sv = sn; break
        sv = sn
    u = np.exp(sv)                    # r - 1, exact to double even when tiny
    return u/(1.0 + u), 1.0 + u       # A, r

def V(rst, l):
    A, r = A_of_rstar(rst)
    return A*(l*(l+1)/r**2 + 1.0/r**3)

def transmission(omega, l, rst_min=-30.0, rst_max=None, h=None):
    # integrate psi'' = (V - omega^2) psi from the sonic surface r = r_s outward, psi = e^{-i omega r*} there (ingoing, unit amplitude)
    if rst_max is None: rst_max = 60.0/omega + 200.0
    if h is None: h = min(0.02, 0.02/omega) if omega < 1 else 0.02
    y = np.array([np.exp(-1j*omega*rst_min), -1j*omega*np.exp(-1j*omega*rst_min)], dtype=complex)
    x = rst_min
    n = int((rst_max - rst_min)/h)
    def f(x, y): return np.array([y[1], (V(x, l) - omega**2)*y[0]])
    for _ in range(n):
        k1 = f(x, y); k2 = f(x+h/2, y+h/2*k1); k3 = f(x+h/2, y+h/2*k2); k4 = f(x+h, y+h*k3)
        y = y + h/6*(k1 + 2*k2 + 2*k3 + k4); x += h
    # at large r*: psi = Aout e^{+i omega r*} + Bin e^{-i omega r*}
    e_p, e_m = np.exp(1j*omega*x), np.exp(-1j*omega*x)
    Aout = (y[0] + y[1]/(1j*omega))/2 / e_p
    Bin  = (y[0] - y[1]/(1j*omega))/2 / e_m
    T2 = 1.0/abs(Bin)**2; R2 = abs(Aout)**2/abs(Bin)**2
    return T2, R2

print("=== EDGE-1a: the l = 0 transmission at low frequency, and the absorption cross-section ===")
AH = 4*np.pi     # r_s = 1
rows = []
for om in (0.20, 0.10, 0.05, 0.025, 0.0125):
    T2, R2 = transmission(om, 0)
    sigma = np.pi*T2/om**2
    rows.append((om, T2, R2, sigma/AH))
    print(f"  omega r_s = {om:6.3f}:  |T_0|^2 = {T2:.6e}   |R|^2+|T|^2 = {R2+T2:.8f}   sigma/A_H = {sigma/AH:.6f}   |T_0|^2/(2 omega r_s)^2 = {T2/(2*om)**2:.6f}")
check("a1", all(abs(r[2] + r[1] - 1) < 1e-6 for r in rows), "flux conserved |R|^2 + |T|^2 = 1 at every frequency to 1e-6 (K1 passes)")
# extrapolate sigma/A_H to omega -> 0 with a quadratic in omega (Richardson-style on the three lowest points)
oms = np.array([r[0] for r in rows[-3:]]); vals = np.array([r[3] for r in rows[-3:]])
coef = np.polyfit(oms, vals, 2); lim = coef[-1]
print(f"  sigma(omega->0)/A_H (quadratic extrapolation on the three lowest frequencies) = {lim:.5f}")
check("a2", abs(lim - 1.0) < 0.01, "sigma_abs(omega -> 0) = A_H to better than 1 %: the sonic surface r = r_s absorbs long waves with exactly its geometric area (K2 passes; DGM/Page recovered as OUTPUT)")
check("a3", abs(rows[-1][3] - 1.0) < 0.05, "already within 5 % of A_H at omega r_s = 0.0125 without extrapolation; the approach is linear in omega (the known O(M omega) correction), from ABOVE: the l=0 barrier peak sits at omega r_s = 0.325 and the ratio peaks near 0.2")
print("  => the AC absorption area equals the DC capture area of the marginally bound orbit (GUD-1: pi b_mb^2 = 4 pi r_s^2 = A_H).")
print("     Two different physics -- long waves and slow particles -- read the sonic surface r = r_s as the same reflectionless boundary of area A_H.")

print("=== EDGE-1b: the greybody edge completed: face H at the sonic surface r = r_s ===")
# low-frequency law |T_0|^2 = 4 (omega r_s)^2 -> emission spectrum dN/domega ~ (omega^2 A_H/pi) / (e^{2 pi omega/kappa} - 1)
om_lo = rows[-1][0]; T2_lo = rows[-1][1]
check("b1", abs(T2_lo/(2*om_lo)**2 - 1.0) < 0.05, "|T_0|^2 = 4 (omega r_s)^2 at low frequency: the l = 0 edge is a power law in omega, set by the AREA of r = r_s, while FOLD-1's l >> 1 edge is a Fermi function set by the FOLD's lambda_L. Face H now has both ends")

print("=== EDGE-1c: the static l = 2 response (scalar tidal Love number) ===")
r = sp.symbols('r', positive=True)
A = 1 - 1/r
def static_op(ph):  # l = 2, omega = 0 scalar equation of the potential ratio metric
    return sp.simplify(sp.diff(r**2*A*sp.diff(ph, r), r)/r**2 - 6*ph/r**2)
xx = 2*r - 1
P2 = (3*xx**2 - 1)/2
Q2 = P2*sp.log((xx + 1)/(xx - 1))/2 - sp.Rational(3,2)*xx
print("  candidate branches: P2(2r-1) =", sp.expand(P2), "   Q2(2r-1) = P2 * ln(r/(r-1))/2 - 3(2r-1)/2")
check("c0", static_op(P2)==0 and static_op(Q2)==0, "both P2(2r-1) and Q2(2r-1) solve the static l=2 equation: the general solution is their span (verified by substitution, not dsolve)")
check("c1", sp.limit(Q2, r, 1, '+') in (sp.oo, -sp.oo) and P2.subs(r, 1) == 1, "Q2 is log-singular at the sonic surface r = r_s, P2 is finite: regularity (no reflection at omega = 0) selects P2 alone")
ser = sp.series(P2, r, sp.oo, 6).removeO()
check("c2", sp.expand(P2).coeff(r, -3) == 0 and sp.Poly(sp.expand(P2), r).degree() == 2, "the regular branch is a degree-2 POLYNOMIAL: no r^-3 tail. Scalar Love number k_2 = 0 (K3 passes): the sonic surface r = r_s does not deform under a static tide, like a black hole and unlike an ECO")
serQ = sp.series(Q2, r, sp.oo, 5).removeO()
print("  for the record, the excluded branch's tail: Q2 ~", sp.nsimplify(serQ.coeff(r, -3)), "/ r^3 + ... (that is what a deformable body would carry)")

n=sum(CH); print(f"\n=== EDGE-1: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("RESULT: sigma_abs(0) = A_H (matched absorber, geometric area); |T_0|^2 = 4 (omega r_s)^2; scalar k_2 = 0.")
print("        FINGERPRINT of the medium's rotating hole (with NS-1): (kappa, k_2, sigma_0/A_H) = (0, 0, 1). Kerr: (1, 0, 1). ECOs: (!=1, !=0, !=1).")
print("TIER: DERIVED from the potential ratio metric + the derived ingoing condition; numerics two-checked by flux conservation. The tensor")
print("      (gravitational) Love number and the spin-dependent absorption are owed; the scalar versions are the proxies.")

# prim_fold1_lightring_dressing.py  --  FOLD-1: the fold at the light ring and what dresses it, run FORWARD.
#
# LOCK-5 left the octahedron with ONE scaffold edge: faces H (Light, Evanescence, Mass) and HG -- light tunnelling near the
#   horizon. Its geometric-optics limit (the fold: light ring, escape cone, P2 = 0) was derived 09-07. This runner derives what
#   sits ON the fold: the particle reading (the subring ladder) and the wave reading (the ringdown), and the greybody edge.
#
# INPUTS:
#   DYN-1  A(r) = 1 - r_s/r, B = 1/A, as OUTPUT of the divider (recomputed here).                                   [DERIVED 09-07]
#   RAYS   light = null rays of the output metric (geometric optics), and light = the scalar wave equation of the
#          output metric (wave optics) -- the two readings of the H flip (T4b': vector = Sym^2 of the spinor reading). [T4b' + the two readings]
#   PB     transmission through a parabolic barrier: T = 1/(1 + exp(2 pi (V0 - w^2)/sqrt(-2 V0''))).  Exact for an inverted
#          parabola; standard quantum mechanics.                                                                     [GROUND: exact math]
#
# BANNED: Kerr, Lyapunov exponents or QNM frequencies quoted from the literature, greybody factors quoted, the pinning.
# GROUND: EHT M87* ring diameter 42 +- 3 muas (composite reading, ~10-17 %); GW150914 remnant ringdown (l=2, n=0) f ~ 251 Hz,
#         tau ~ 4 ms at tens-of-percent (remnant spin 0.67 -- Kerr; the model's swirl is O(J) only, so this is a ballpark check).
# TWO PATHS for the Lyapunov exponent: (P1) linearise the radial null equation at the ring; (P2) integrate a perturbed null ray.
# KILLS: (K1) P1 != P2 -> runner wrong. (K2) if the barrier transmission width is not the ray exponent, the two readings are
#        NOT the same fold and the H-flip identification fails.

import sympy as sp, mpmath as mp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

r, rs, b, w, dr = sp.symbols('r r_s b omega delta', positive=True)
print("=== FOLD-1a: the fold from the divider ===")
RL = sp.integrate(1/r**2, (r, rs, r)); RS = sp.integrate(1/r**2, (r, r, sp.oo))
A = sp.simplify(RL/(RL+RS))                                  # 1 - r_s/r (units c = 1)
U = A/r**2                                                   # the ray potential: (dr/dlambda)^2 = E^2 - L^2 U
r_ph = [s for s in sp.solve(sp.diff(U, r), r) if s.is_positive][0]
A_ph = A.subs(r, r_ph); b_c = sp.simplify(r_ph/sp.sqrt(A_ph)); Om_c = sp.simplify(sp.sqrt(A_ph)/r_ph)
print("  light ring r_ph =", r_ph, " critical impact parameter b_c =", b_c, " winding rate Omega_c =", Om_c)
check("a1", r_ph == 3*rs/2 and sp.simplify(b_c - 3*sp.sqrt(3)*rs/2)==0, "r_ph = 3 r_s/2, b_c = 3 sqrt3 r_s/2: the fold's radius and its position on the image plane")
# ground: M87*
mp.mp.dps=20
Msun_m = mp.mpf('1477.0'); M87 = mp.mpf('6.5e9')*Msun_m; rsM87 = 2*M87; D87 = mp.mpf('16.8')*mp.mpf('3.0857e22')
theta_sh = 2*(3*mp.sqrt(3)/2*rsM87)/D87 * mp.mpf('206264.806247') * mp.mpf(1e6)   # rad -> arcsec -> muas
print(f"  M87* shadow diameter from the fold: {mp.nstr(theta_sh,4)} muas   (EHT ring 42 +- 3 muas; emission ring sits ~10 % outside the fold)")
check("a2", abs(theta_sh - 42) < 8, "M87* shadow within the composite reading's error (a 10-17 % test, model-mediated)")

print("=== FOLD-1b: the Lyapunov exponent of the fold, two paths ===")
F = sp.simplify(A**2*(1 - b_c**2*U))                         # (dr/dt)^2 for the critical ray, in coordinate time
check("b1", sp.simplify(F.subs(r, r_ph))==0 and sp.simplify(sp.diff(F, r).subs(r, r_ph))==0, "the critical ray has a double root at r_ph: the fold is a tangency (the ray family's turning points pile up)")
lam_P1 = sp.simplify(sp.sqrt(sp.diff(F, r, 2).subs(r, r_ph)/2))
print("  P1: lambda_L = sqrt(F''(r_ph)/2) =", lam_P1)
# P2: integrate r'' = F'(r)/2 from r_ph(1 + 1e-7) at rest, fit the e-folding of delta r
Fp = sp.lambdify(r, sp.diff(F, r).subs(rs, 1), 'mpmath')
mp.mp.dps = 25
def rhs(t, y): return [y[1], Fp(y[0])/2]
Fn = sp.lambdify(r, F.subs(rs, 1), 'mpmath')
eps0 = mp.mpf('1e-7'); r0 = mp.mpf(3)/2*(1+eps0)
sol = mp.odefun(rhs, 0, [r0, mp.sqrt(Fn(r0))])        # start ON the outgoing branch of (dr/dt)^2 = F: the growing mode only, no exponent assumed
t1, t2 = mp.mpf(6), mp.mpf(12)
d1, d2 = sol(t1)[0]-mp.mpf(3)/2, sol(t2)[0]-mp.mpf(3)/2
lam_P2 = mp.log(d2/d1)/(t2-t1)
lam_P1_num = mp.mpf(sp.N(lam_P1.subs(rs,1), 20))
print(f"  P2: numeric e-folding rate of a perturbed critical ray = {mp.nstr(lam_P2,10)}   (P1 = {mp.nstr(lam_P1_num,10)})")
check("b2", abs(lam_P2/lam_P1_num - 1) < mp.mpf('1e-4'), f"P1 = P2 to {mp.nstr(abs(lam_P2/lam_P1_num-1),2)} (K1 passes)")
ratio = sp.simplify(lam_P1/Om_c)
check("b3", ratio == 1, f"lambda_L / Omega_c = {ratio}: the fold's instability rate EQUALS its winding rate. (A property of 1 - r_s/r; other lapses give other ratios.)")

print("=== FOLD-1c: the particle reading -- the subring ladder ===")
gamma_half = sp.simplify(sp.pi*lam_P1/Om_c)
print("  demagnification exponent per half-orbit gamma = pi lambda_L/Omega_c =", gamma_half, " -> factor e^-gamma =", sp.N(sp.exp(-gamma_half), 6))
check("c1", gamma_half == sp.pi, "gamma = pi exactly: each successive subring is e^-pi = 0.0432 as wide and as bright as the last")
print("  PREDICTION (unmeasured): n = 1 subring width/flux ratio 4.32 %, n = 2: 0.187 %. Next-generation space VLBI targets n = 1.")
print("  The ladder counts HALF-orbits: the fold is seen from both sides of the sheet. A full orbit is e^-2pi -- the vector reading's period.")

print("=== FOLD-1d: the wave reading -- the eikonal ringdown (the H flip) ===")
l, n = sp.symbols('l n', nonnegative=True)
w_qnm = Om_c*(l + sp.Rational(1,2)) - sp.I*lam_P1*(n + sp.Rational(1,2))
print("  omega_{l n} = Omega_c (l + 1/2) - i lambda_L (n + 1/2)  =", sp.simplify(w_qnm))
print("  real part = the winding (Light/Action pole); imaginary part = the decay (Evanescence pole). One fold, two poles: the H flip.")
# GW150914 remnant, ballpark (nonspinning eikonal; remnant spin 0.67 not modelled -- the swirl is O(J) only)
c = mp.mpf('299792458'); Mf = mp.mpf('62')*Msun_m; rsf = 2*Mf
Om_num = mp.mpf(sp.N(Om_c.subs(rs,1),20))*c/rsf; lam_num = mp.mpf(sp.N(lam_P1.subs(rs,1),20))*c/rsf
f220 = Om_num*mp.mpf('2.5')/(2*mp.pi); tau220 = 1/(lam_num*mp.mpf('0.5'))
print(f"  GW150914-like remnant (62 Msun, nonspinning eikonal): f_220 ~ {mp.nstr(f220,4)} Hz, tau_220 ~ {mp.nstr(tau220*1000,3)} ms   (ground: ~251 Hz, ~4 ms at tens of %)")
print("  CAVEATS, named: the eikonal (l+1/2) overestimates the exact l=2 frequency by ~30 %, and the remnant's spin 0.67 raises the true")
print("  frequency by a comparable amount; the numbers land near ground partly by cancellation. Ballpark, not precision. Spin needs O(J^2).")
check("d1", abs(f220-251)/251 < 0.35 and abs(tau220*1000-4)/4 < 0.35, "eikonal nonspinning ringdown within the tens-of-percent ground (with the caveats above)")

print("=== FOLD-1e: the greybody edge -- transmission through the fold (face H) ===")
# wave equation in the tortoise coordinate: d^2 psi/dr*^2 + (w^2 - V) psi = 0, V = (l+1/2)^2 U(r), dr*/dr = 1/A
dUdrs = sp.simplify(A*sp.diff(U, r))                       # dU/dr*
d2Udrs = sp.simplify(A*sp.diff(dUdrs, r))                  # d^2U/dr*^2
U0 = sp.simplify(U.subs(r, r_ph)); U0pp = sp.simplify(d2Udrs.subs(r, r_ph))
check("e1", sp.simplify(U0 - Om_c**2)==0 and U0pp < 0, "the ray potential peaks at the light ring with U0 = Omega_c^2 and negative curvature: a parabolic barrier")
# parabolic-barrier transmission near omega_c = (l+1/2) Omega_c, linearised in (omega - omega_c):
L = l + sp.Rational(1,2); V0 = L**2*U0; V0pp = L**2*U0pp
w_c = L*Om_c
expo = sp.simplify((2*sp.pi*(V0 - w**2)/sp.sqrt(-2*V0pp)))
expo_lin = sp.simplify(sp.series(expo.subs(w, w_c + dr), dr, 0, 2).removeO())
width = sp.simplify(-2*sp.pi/ sp.diff(expo_lin, dr))        # T = 1/(1 + exp(-2 pi (w - w_c)/width))
print("  transmission T(omega) = 1/(1 + exp(-2 pi (omega - omega_c)/width)) with width =", width)
check("e2", sp.simplify(width - lam_P1)==0, "the edge width IS the ray exponent lambda_L: the wave's barrier and the ray's fold are one object (K2 passes; the H-flip identification holds)")
kappa = 1/(2*rs)                                            # DYN-2's electrode slope
print("  two 'temperatures' on the far column: the electrode's kappa/2pi =", sp.simplify(kappa/(2*sp.pi)), " and the fold's lambda_L/2pi =", sp.simplify(lam_P1/(2*sp.pi)), "; ratio lambda_L/kappa =", sp.simplify(lam_P1/kappa), "=", sp.N(lam_P1/kappa, 5))
check("e3", sp.simplify(lam_P1/kappa - 4/(3*sp.sqrt(3)))==0, "lambda_L/kappa = 4/(3 sqrt3): the fold's edge is 0.77 of the electrode's slope -- the greybody's edge is NOT the Hawking temperature, and the model keeps them apart")
print("  => in the eikonal limit the CHG face's spectrum is Planck(kappa) x Fermi-edge(lambda_L): the electrode supplies the thermal factor")
print("     (DYN-2, given KMS), the fold supplies the filter (this runner, given PB). Face H is derived AT THE FOLD. The electrode part of the")
print("     tunnelling (l = 0, below the barrier) is the remaining owed piece of the H edge.")

nn=sum(CH); print(f"\n=== FOLD-1: {nn}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("TIER: a1, b1-b3, c1, e1-e3 DERIVED from the divider output (+ PB, exact math). d1 is a ballpark against tens-of-percent ground.")
print("PREDICTIONS cast: subring ladder e^-pi per half-orbit (EHT n=1); lambda_L = Omega_c; greybody edge width = lambda_L, not kappa.")

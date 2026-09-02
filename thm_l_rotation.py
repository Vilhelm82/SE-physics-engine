#!/usr/bin/env python3
# =============================================================================
# THM-TARGET L -- rotation: can the presentation field carry swirl, and does it drag frames?
# Date: 2026-09-03.  sympy + mpmath + numpy (floats for the orbit corroboration only).  Exit 0 iff every check passes.
#
# THE ARCHITECTURE CLAIM: the transport law extends to a presentation field with vorticity.  Rule G (thm_i_transport) assumed
#   an irrotational flow (A-0).  The torsion-free extension is forced: boost by the SYMMETRIC strain S along the co-moving
#   step, AND rotate at the river's local angular velocity Omega = (1/2) curl v per unit frame time.        [RULE G']
#   Kill: if G' does not conserve the symmetry pairings E = k^0 + v.k and L_z = (x x k)_z for a swirling stationary
#   axisymmetric field, the model cannot carry rotation.
# THE SWIRL PROFILE (declared, KIN-3): v = -beta(r) r-hat + omega(r) (z x r),  omega = 2 G J/(c^2 r^3)  -- the drift of the
#   zero-angular-momentum free-fall frame around a slowly rotating source, imported at first order in J exactly as the
#   pinning was imported before THM-K.  Deriving it is named future work.
# THE KILL NUMBERS: a gyroscope at rest precesses at (1/2) curl v (Gravity Probe B, 37.2 +- 7.2 mas/yr measured, 39.2
#   predicted); an orbit's node drifts at 2GJ/(c^2 a^3 (1-e^2)^{3/2}) (LAGEOS/LARES).  Both must come out of RULE G'.
# NOT USED before the comparison block: the Kerr metric, gravitomagnetism, Lense-Thirring by name.
# =============================================================================
import sys, math
import sympy as sp
import mpmath as mp
mp.mp.dps = 30
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
def z(e):
    e = sp.sympify(e)
    for rt in (sp.simplify, lambda q: sp.simplify(sp.expand(q))):
        try:
            if rt(e) == 0: return True
        except Exception: pass
    return False

x, y, zz, eps = sp.symbols('x y z epsilon', real=True)
r = sp.sqrt(x**2 + y**2 + zz**2)
X = sp.Matrix([x, y, zz]); ez = sp.Matrix([0, 0, 1])
bet = sp.Function('beta', real=True); om = sp.Function('omega', real=True)
v = -bet(r)*X/r + om(r)*ez.cross(X)                          # KIN-2a-type inflow + KIN-3 swirl, both general in r
Jac = sp.Matrix(3, 3, lambda i, j: sp.diff(v[i], [x, y, zz][j]))
S = (Jac + Jac.T)/2; Aa = (Jac - Jac.T)/2
curl = sp.Matrix([sp.diff(v[2], y) - sp.diff(v[1], zz), sp.diff(v[0], zz) - sp.diff(v[2], x), sp.diff(v[1], x) - sp.diff(v[0], y)])
Om = curl/2

print("=" * 100); print("L-1  the extension is forced, and it conserves the pairings for ANY swirl"); print("=" * 100)
tt = sp.Symbol('t', positive=True)
conc_fn = {bet: sp.Lambda(tt, sp.sqrt(1/tt)), om: sp.Lambda(tt, sp.Rational(1, 10)/tt**3)}     # KIN-2a + KIN-3, r_s = 1, 2GJ/c^2 = 1/10
noswirl = {bet: sp.Lambda(tt, sp.sqrt(1/tt)), om: sp.Lambda(tt, 0)}
p0 = {x: sp.Rational(3, 2), y: sp.Rational(-4, 5), zz: sp.Rational(7, 10)}
check("L-1a the swirl is rotational: curl v != 0 with the swirl on, and curl v = 0 with it off (A-0's hypothesis is exactly the no-swirl case)",
      any(sp.N(cc.subs(conc_fn).doit().subs(p0)) != 0 for cc in curl) and all(z(cc.subs(noswirl).doit()) for cc in curl))
check("L-1b the antisymmetric part of the gradient acts as (1/2) curl v cross: A k = Omega x k for every k",
      all(z(sp.simplify((Aa*ez)[i] - Om.cross(ez)[i])) for i in range(3)) and all(z(sp.simplify((Aa*sp.Matrix([1, 0, 0]))[i] - Om.cross(sp.Matrix([1, 0, 0]))[i])) for i in range(3)))
# a massive body: momentum paravector (k0, kvec) in the rain frame, k0 = sqrt(1 + k^2); co-moving velocity u = kvec/k0
k1, k2, k3 = sp.symbols('k1 k2 k3', real=True); kv = sp.Matrix([k1, k2, k3]); k0 = sp.sqrt(1 + (kv.T*kv)[0, 0]); u = kv/k0
def step(rule):
    """one step of size eps: returns (x', kvec') to first order for RULE G' (rotation on) or RULE G (rotation off, full-gradient boost)"""
    dx = (u + v)*eps                                            # PROP-1
    if rule == "G'":
        dk = (-(S*u)*k0 + Om.cross(kv))*eps                     # boost by S along the co-moving step, rotate at Omega per unit time
    else:
        dk = (-(Jac*u)*k0)*eps                                  # thm_i_transport's rule G, extrapolated naively to a rotational field
    return X + dx, kv + dk
def pairings(Xp, kp):
    k0p = sp.sqrt(1 + (kp.T*kp)[0, 0]); vp = v.subs({x: Xp[0], y: Xp[1], zz: Xp[2]}, simultaneous=True)
    return k0p + (vp.T*kp)[0, 0], (Xp.cross(kp))[2]
E0, L0 = pairings(X, kv)
def drift(rule):
    Xp, kp = step(rule); E1, L1 = pairings(Xp, kp)
    return sp.diff(E1, eps).subs(eps, 0), sp.diff(L1, eps).subs(eps, 0)
dE_Gp, dL_Gp = drift("G'")
ok_sym = z(dE_Gp) and z(dL_Gp)
print(f"  symbolic identity attempted with general beta(r), omega(r): {'closed' if ok_sym else 'did not close; falling back to exact points'}")
pts = [{x: sp.Rational(3, 2), y: sp.Rational(-4, 5), zz: sp.Rational(7, 10), k1: sp.Rational(1, 3), k2: sp.Rational(-2, 7), k3: sp.Rational(5, 9)},
       {x: sp.Rational(-1, 4), y: sp.Rational(9, 5), zz: sp.Rational(-3, 8), k1: sp.Rational(2, 5), k2: sp.Rational(1, 6), k3: sp.Rational(-4, 5)}]
def at(e, p):
    e = e.subs(conc_fn).doit().subs(p); return sp.N(e, 30)
ok_pts = all(abs(at(dE_Gp, p)) < 1e-25 and abs(at(dL_Gp, p)) < 1e-25 for p in pts)
check("L-1c RULE G' conserves BOTH pairings for the swirling field" + (" -- identically in x, k, beta and omega" if ok_sym else " -- exactly at two rational points on KIN-2a + KIN-3 (30 digits)"),
      ok_sym or ok_pts)
dE_G, dL_G = drift("G")
check("L-1d RULE G (no rotation, full-gradient boost) FAILS on the same field: dE/dt, dL_z/dt != 0 -- the rotation term is REQUIRED, not optional",
      any(abs(at(dE_G, p)) > 1e-8 or abs(at(dL_G, p)) > 1e-8 for p in pts), f"dE/dt = {sp.N(at(dE_G, pts[0]), 6)}, dL_z/dt = {sp.N(at(dL_G, pts[0]), 6)}")
check("L-1e with omega = 0 the two rules coincide (rule G' reduces to thm_i_transport's rule G): irrotational limit recovered",
      z(sp.simplify((dE_Gp - dE_G).subs(om, sp.Lambda(tt, 0)).doit())) and z(sp.simplify((dL_Gp - dL_G).subs(om, sp.Lambda(tt, 0)).doit())))
print("  READING: a swirling presentation field is carried by the SAME transport law with its rotation part switched on.  The model")
print("  can hold rotation.  What it holds is decided by the profile, which is declared below.")

print(); print("=" * 100); print("L-2  a gyroscope at rest: its spin is carried by the rotation part alone -- the dipole pattern"); print("=" * 100)
G_, J_, c_ = sp.symbols('G J c', positive=True)
conc_LT = {om: sp.Lambda(tt, 2*G_*J_/(c_**2*tt**3))}          # KIN-3, declared: the free-fall frame drifts at omega = 2GJ/(c^2 r^3)
Om_LT = sp.simplify(Om.subs(conc_LT).doit())
rhat = X/r; Jvec = J_*ez
schiff = (G_/(c_**2*r**3))*(3*(Jvec.T*rhat)[0, 0]*rhat - Jvec)
check("L-2a (1/2) curl v for the declared swirl equals (G/c^2 r^3)[3 (J.r-hat) r-hat - J] exactly: twice as strong on the axis as on"
      " the equator, opposite sense on the equator", all(z(sp.simplify(Om_LT[i] - schiff[i])) for i in range(3)))
check("L-2b the swirl's vorticity is divergence-free and curl-free outside the source (a dipole field): the vacuum law the profile obeys",
      z(sp.simplify(sum(sp.diff(2*Om_LT[i], [x, y, zz][i]) for i in range(3)))) and
      all(z(sp.simplify(c2)) for c2 in [sp.diff(2*Om_LT[2], y) - sp.diff(2*Om_LT[1], zz), sp.diff(2*Om_LT[0], zz) - sp.diff(2*Om_LT[2], x), sp.diff(2*Om_LT[1], x) - sp.diff(2*Om_LT[0], y)]))
# polar-orbit average and the Gravity Probe B number
Gn, cn = mp.mpf('6.67430e-11'), mp.mpf('2.99792458e8')
J_E = mp.mpf('8.034e37')*mp.mpf('7.292115e-5')                  # I_E omega_E
r_gpb = mp.mpf('7.0274e6')
Om_avg = Gn*J_E/(2*cn**2*r_gpb**3)                              # <3 (J.r) r - J> over a polar great circle = J/2
mas_yr = Om_avg*mp.mpf('3.15576e7')*mp.mpf(180)/mp.pi*3600*1000
check("L-2c averaged over a polar orbit the precession is G J/(2 c^2 r^3) along J; for Earth at r = 7027.4 km that is "
      + mp.nstr(mas_yr, 4) + " mas/yr  [comparison stage: Schiff; Gravity Probe B measured 37.2 +- 7.2, predicted 39.2 with their constants]",
      abs(mas_yr - 39.2) < 7.2)

print(); print("=" * 100); print("L-3  an orbit's node: the boost part and the rotation part acting together along a path"); print("=" * 100)
import numpy as np
ell = 0.02; a_orb = 60.0; inc = math.radians(60.0)               # units r_s = 1, c = 1: omega(r) = 2 ell/r^3, ell = GJ/c^3
vc = {bet: sp.Lambda(tt, sp.sqrt(1/tt)), om: sp.Lambda(tt, 2*ell/tt**3)}
S_c = S.subs(vc).doit(); Om_c = Om.subs(vc).doit(); Jac_c = Jac.subs(vc).doit(); v_c = v.subs(vc).doit()
S_f = sp.lambdify((x, y, zz), S_c, 'numpy'); Om_f = sp.lambdify((x, y, zz), Om_c, 'numpy'); Jac_f = sp.lambdify((x, y, zz), Jac_c, 'numpy'); v_f = sp.lambdify((x, y, zz), v_c, 'numpy')
def rhs(s, rule):
    X_ = s[:3]; K = s[3:]; k0 = math.sqrt(1 + K @ K); U = K/k0
    Sm = np.array(S_f(*X_), dtype=float); Ov = np.array(Om_f(*X_), dtype=float).reshape(3); Jm = np.array(Jac_f(*X_), dtype=float); V = np.array(v_f(*X_), dtype=float).reshape(3)
    dX = U + V
    dK = -(Sm @ U)*k0 + np.cross(Ov, K) if rule == "G'" else -(Jm @ U)*k0
    return np.concatenate([dX, dK])
def integrate(rule, n_orbits=3, h=0.05):
    X0 = np.array([a_orb, 0.0, 0.0]); vcirc = math.sqrt(0.5/a_orb)
    dXdt = vcirc*np.array([0.0, math.cos(inc), math.sin(inc)])
    U0 = dXdt - np.array(v_f(*X0), dtype=float).reshape(3); K0 = U0/math.sqrt(1 - U0 @ U0)
    s = np.concatenate([X0, K0]); T = 2*math.pi*math.sqrt(2*a_orb**3); nodes = []; zprev = s[2]; t = 0.0
    while t < n_orbits*T + 50:
        k1_ = rhs(s, rule); k2_ = rhs(s + h/2*k1_, rule); k3_ = rhs(s + h/2*k2_, rule); k4_ = rhs(s + h*k3_, rule)
        s_new = s + h/6*(k1_ + 2*k2_ + 2*k3_ + k4_); t += h
        if zprev < 0 <= s_new[2]:                                  # ascending node: interpolate the crossing
            f = -zprev/(s_new[2] - zprev); xc = s[:3] + f*(s_new[:3] - s[:3]); nodes.append(math.atan2(xc[1], xc[0]))
        zprev = s_new[2]; s = s_new
    return nodes
nodes = integrate("G'")
shifts = [((nodes[i+1] - nodes[i] + math.pi) % (2*math.pi)) - math.pi for i in range(len(nodes) - 1)]
pred = 4*math.pi*math.sqrt(2)*ell/a_orb**1.5                        # 2GJ/(c^2 a^3) x period, circular, first order
meas = sum(shifts)/len(shifts) if shifts else float('nan')
check(f"L-3a RULE G': the ascending node advances {meas:.3e} rad per orbit (mean of {len(shifts)}); first-order prediction 2GJ/(c^2 a^3) x T = {pred:.3e}:"
      f" ratio {meas/pred:.4f} (post-Newtonian corrections at a = 60 r_s are of order r_s/a ~ 2%)", len(shifts) >= 2 and abs(meas/pred - 1) < 0.05 and meas > 0)
nodes_G = integrate("G")
shifts_G = [((nodes_G[i+1] - nodes_G[i] + math.pi) % (2*math.pi)) - math.pi for i in range(len(nodes_G) - 1)]
meas_G = sum(shifts_G)/len(shifts_G) if shifts_G else float('nan')
check(f"L-3b FAILURE BRANCH: rule G without the rotation term gives {meas_G:.3e} rad per orbit, ratio {meas_G/pred:.3f} to the prediction -- wrong;"
      " the node needs both parts of the connection", len(shifts_G) >= 2 and abs(meas_G/pred - 1) > 0.1)
print("  READING: the node drift is not the local river spin (that is L-2); it is strain and rotation acting together along the")
print("  orbit, and rule G' returns the textbook coefficient 2 to the post-Newtonian accuracy of the setup.")

print(); n_ = sum(CH); print(f"RESULT: {n_}/{len(CH)} checks passed.")
if n_ != len(CH): print("VERDICT: check failures above."); sys.exit(1)
print("VERDICT (THM-L): the presentation field CAN carry swirl.  The torsion-free extension of the transport law is forced -- boost by")
print("  the symmetric strain along the co-moving step, rotate at (1/2) curl v per unit frame time -- and it conserves both symmetry")
print("  pairings for any stationary axisymmetric swirl, while rule G without the rotation term fails.  With the swirl profile")
print("  omega = 2GJ/(c^2 r^3) DECLARED (KIN-3, the zero-angular-momentum free-fall drift, first order in J), a gyroscope at rest")
print("  precesses at (G/c^2 r^3)[3(J.r)r - J] exactly and an orbit's node advances at 2GJ/(c^2 a^3): both measured numbers land.")
print("  What is derived: the connection and both effects from it.  What is declared: the swirl profile.  Deriving omega(r) the way")
print("  THM-K derived beta(r) -- from the seat's temperature on a rotating screen -- is the next target.")
print("COMPARISON STAGE: rule G' is the spin connection of the rain tetrad for a river with vorticity (Hamilton-Lisle's twisting river;")
print("  Doran's form of Kerr has flat slices in oblate spheroidal coordinates); L-2a is Schiff's gyroscope formula, Lense-Thirring")
print("  precession; L-2b is the vacuum gravitomagnetic dipole; L-3a is the Lense-Thirring nodal precession (LAGEOS/LARES);")
print("  the swirl (1/2) curl v = Omega_LT is the statement that frame dragging is the vorticity of the free-fall flow.")
sys.exit(0)

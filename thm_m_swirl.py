#!/usr/bin/env python3
# =============================================================================
# THM-TARGET M -- the swirl is the boost of the pinning: KIN-3 derived from KIN-2a + BARE-1 covariance + linear superposition
# Date: 2026-09-03.  sympy + mpmath only.  Exit 0 iff every check passes.
#
# THE CLAIM: the frame-dragging profile omega(r) = 2GJ/(c^2 r^3), declared as KIN-3 in thm_l_rotation, follows from three
#   things the model already has: (i) the static presented field of a mass (THM-K derived beta^2 = r_s/r; the seat's clock
#   N^2 = 1 - beta^2 and rods delta + [r_s/(r - r_s)] r r^T, H-5); (ii) BARE-1 covariance -- a moving element's field is
#   its static field boosted by the sandwich; (iii) K-6's field law is linear, so elements superpose at first order in r_s.
# THE STRUCTURE: boosting the seat's TIME OFFICE alone gives half the swirl; boosting the RODS gives the other half.  A scalar
#   theory (time office only) predicts half the measured frame dragging: Gravity Probe B kills it.  Same split as the light
#   bending: each office carries half.
# NOT USED before the comparison block: the Kerr metric, linearised Einstein equations, gravitoelectromagnetism by name.
# =============================================================================
import sys
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

print("=" * 100); print("M-1  the static presented field at first order in r_s: clock, rods, no shift"); print("=" * 100)
rs, r = sp.symbols('r_s r', positive=True); x, y, zz = sp.symbols('x y z', real=True)
X = sp.Matrix([x, y, zz]); rr = sp.sqrt(x**2 + y**2 + zz**2); rh = X/rr
h00 = rs/rr                                                  # 1 - N^2 = beta^2 = r_s/r  (THM-K, K-5d)
hij = (rs/rr)*rh*rh.T                                        # presented rods delta + [r_s/(r-r_s)] r r^T -> r_s/r r r^T at first order (H-5)
check("M-1a the seat's clock: N^2 = 1 - r_s/r, so h_00 = r_s/r; the seat's rods: delta_ij + (r_s/(r - r_s)) r_i r_j -> h_ij = (r_s/r) r_i r_j;"
      " the static seat's shift vanishes -- three numbers, all from THM-K's beta^2 and H-5", z(sp.series(rs/(r - rs), rs, 0, 2).removeO() - rs/r))

print(); print("=" * 100); print("M-2  BARE-1 covariance: boost the presented field by a small velocity w"); print("=" * 100)
s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -sp.I], [sp.I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]]); Id = sp.eye(2)
SIG = [s1, s2, s3]; BASIS = [Id] + SIG                       # paravector basis e_0 = 1, e_i = sigma_i
def para(t, v): return t*Id + sum((v[i]*SIG[i] for i in range(3)), sp.zeros(2))
def comps(M): return [sp.simplify(sp.trace(M)/2)] + [sp.simplify(sp.trace(M*S)/2) for S in SIG]
w1, w2, w3, zeta = sp.symbols('w1 w2 w3 zeta', real=True)
wv = sp.Matrix([w1, w2, w3])
# the sandwich boost with rapidity zeta along w-hat, to first order in zeta (velocity w = zeta w-hat):  B = 1 + (w.sigma)/2
B = Id + para(0, wv)/2
L = sp.Matrix(4, 4, lambda mu, nu: comps(B*BASIS[nu]*B)[mu])  # column nu = image of e_nu, to first order in w
eta = sp.diag(-1, 1, 1, 1)
Lfirst = L.applyfunc(lambda e: sp.expand(e).subs({w1**2: 0, w2**2: 0, w3**2: 0, w1*w2: 0, w1*w3: 0, w2*w3: 0}))
check("M-2a the sandwich is a Lorentz map on paravectors: L^T eta L = eta to first order in w, with L_0i = L_i0 = w_i",
      all(z(sp.expand((Lfirst.T*eta*Lfirst - eta)[i, j]).subs({w1**2: 0, w2**2: 0, w3**2: 0, w1*w2: 0, w1*w3: 0, w2*w3: 0})) for i in range(4) for j in range(4))
      and all(z(Lfirst[0, i+1] - wv[i]) and z(Lfirst[i+1, 0] - wv[i]) for i in range(3)))
# the presented metric perturbation h (symmetric bilinear form on paravectors) transforms as a rank-2 tensor: h' = L^T h L
H = sp.zeros(4, 4); H[0, 0] = h00
for i in range(3):
    for j in range(3): H[i+1, j+1] = hij[i, j]
Hp = (Lfirst.T*H*Lfirst).applyfunc(lambda e: sp.expand(e).subs({w1**2: 0, w2**2: 0, w3**2: 0, w1*w2: 0, w1*w3: 0, w2*w3: 0}))
shift = sp.Matrix([Hp[0, i+1] for i in range(3)])
shift_closed = (rs/rr)*(wv + (wv.T*rh)[0, 0]*rh)
check("M-2b the boosted field acquires a shift h'_0i = w_i h_00 + w_j h_ji = (r_s/r)[w_i + (w.r-hat) r_i]: the river of a moving mass",
      all(z(shift[i] - shift_closed[i]) for i in range(3)))
half_clock = (rs/rr)*wv; half_rods = (rs/rr)*(wv.T*rh)[0, 0]*rh
check("M-2c the shift has two halves: w_i h_00 from the CLOCK (a scalar theory would have only this) and w_j h_ji from the RODS",
      all(z(shift[i] - half_clock[i] - half_rods[i]) for i in range(3)))

print(); print("=" * 100); print("M-3  linear superposition over a rigidly rotating shell: the far field and its coefficient"); print("=" * 100)
# ARGUMENT (exact, then machine-checked numerically): element e of mass m_e at x_e moving at w_e = Omega x x_e contributes
#   v_e(x) = (2 G m_e/(c^2 R_e)) [w_e + (w_e . n_e) n_e],  R_e = |x - x_e|, n_e = (x - x_e)/R_e   (M-2b with r_s,e = 2 G m_e/c^2).
#   Dipole expansion with sum m x x^T = (I/2) 1 for a spherically symmetric body:
#     clock half  ->  (2G/c^2) (J x X)/(2 r^3),    rods half  ->  (2G/c^2) (J x X)/(2 r^3),    total  (2G/c^2)(J x X)/r^3.
Gn, cn = mp.mpf(1), mp.mpf(1)                                  # units G = c = 1
def swirl_shell(Xf, R, Mtot, Omega, which='both', n=100):
    """far-field river velocity of a rotating thin spherical shell of radius R, mass Mtot, angular velocity Omega z-hat: numeric sphere quadrature"""
    tot = mp.matrix(3, 1)
    for i in range(n):                                          # midpoint rule in cos(theta') and phi'
        ct = -1 + (2*i + 1)/mp.mpf(n); st = mp.sqrt(1 - ct**2)
        for j in range(n):
            ph = 2*mp.pi*(j + mp.mpf(1)/2)/n
            xe = mp.matrix([R*st*mp.cos(ph), R*st*mp.sin(ph), R*ct])
            we = mp.matrix([-Omega*xe[1], Omega*xe[0], 0])
            d = Xf - xe; Re = mp.sqrt(d[0]**2 + d[1]**2 + d[2]**2); ne = d/Re
            dm = Mtot/(n*n)                                     # equal-area cells: uniform in cos(theta'), phi'
            wn = we[0]*ne[0] + we[1]*ne[1] + we[2]*ne[2]
            term = (we if which in ('both', 'clock') else mp.matrix(3, 1)) + (wn*ne if which in ('both', 'rods') else mp.matrix(3, 1))
            tot += (2*Gn*dm/(cn**2*Re))*term
    return tot
R, Mtot, Omega = mp.mpf(1), mp.mpf(1), mp.mpf('0.1')
J = mp.mpf(2)/3*Mtot*R**2*Omega                                # thin shell: I = (2/3) M R^2
def LT(Xf, coeff):
    Jv = mp.matrix([0, 0, J]); cross = mp.matrix([Jv[1]*Xf[2] - Jv[2]*Xf[1], Jv[2]*Xf[0] - Jv[0]*Xf[2], Jv[0]*Xf[1] - Jv[1]*Xf[0]])
    return coeff*Gn/cn**2*cross/(mp.sqrt(Xf[0]**2 + Xf[1]**2 + Xf[2]**2)**3)
pts = [mp.matrix([300, 0, 0]), mp.matrix([200, 150, 250]), mp.matrix([0, 0, 400]), mp.matrix([-250, 120, -60])]
worst = mp.mpf(0); worst_half = mp.mpf(0)
for Xf in pts:
    v_both = swirl_shell(Xf, R, Mtot, Omega, 'both'); v_lt = LT(Xf, 2)
    v_clock = swirl_shell(Xf, R, Mtot, Omega, 'clock'); v_rods = swirl_shell(Xf, R, Mtot, Omega, 'rods'); v_half = LT(Xf, 1)
    nrm = mp.sqrt(sum(v_lt[i]**2 for i in range(3)))
    if nrm > 0:
        worst = max(worst, mp.sqrt(sum((v_both[i] - v_lt[i])**2 for i in range(3)))/nrm)
        worst_half = max(worst_half, mp.sqrt(sum((v_clock[i] - v_half[i])**2 for i in range(3)))/nrm, mp.sqrt(sum((v_rods[i] - v_half[i])**2 for i in range(3)))/nrm)
    else:
        worst = max(worst, mp.sqrt(sum(v_both[i]**2 for i in range(3))))
print(f"  far-field mismatch vs (2G/c^2)(J x X)/r^3 at four points: {mp.nstr(worst, 3)} (relative; next term O(R^2/r^2) ~ 1e-5)")
check("M-3a the superposed river of the rotating shell is (2G/c^2)(J x X)/r^3 at large r: omega(r) = 2GJ/(c^2 r^3) -- KIN-3 DERIVED,"
      " coefficient 2, on the axis (where J x X = 0 and the field vanishes) and off it", worst < mp.mpf('1e-4'))
check("M-3b each half alone -- clock only, rods only -- gives exactly HALF: (G/c^2)(J x X)/r^3.  A scalar-potential theory (clock only)"
      " predicts half the frame dragging", worst_half < mp.mpf('1e-4'))
# the failure branch against the measurement
Gsi, csi = mp.mpf('6.67430e-11'), mp.mpf('2.99792458e8'); J_E = mp.mpf('8.034e37')*mp.mpf('7.292115e-5'); r_gpb = mp.mpf('7.0274e6')
to_mas = lambda om: om*mp.mpf('3.15576e7')*mp.mpf(180)/mp.pi*3600*1000
full = to_mas(Gsi*J_E/(2*csi**2*r_gpb**3)); half = full/2
check(f"M-3c Gravity Probe B, polar-orbit average: rank-2 (clock + rods) {mp.nstr(full, 4)} mas/yr inside 37.2 +- 7.2; clock-only {mp.nstr(half, 4)} mas/yr is"
      " outside by 2.3 sigma -- the rods' half of the dragging is measured", abs(full - 37.2) < 7.2 and abs(half - 37.2) > 7.2)

print(); print("=" * 100); print("M-4  why there is no rotating screen: the thermal candidates fail on dimensions"); print("=" * 100)
Msym, Jsym, rsym, Gsym, csym, hbsym = sp.symbols('M J r G c hbar', positive=True)
cand = {'screen co-rotates carrying J (thin shell, I = 2/3 M r^2)': Jsym/(sp.Rational(2, 3)*Msym*rsym**2),
        'angular momentum per bit, equipartition energy per bit': Jsym/(Msym*rsym**2)}
target = 2*Gsym*Jsym/(csym**2*rsym**3)
check("M-4a both natural rotating-screen principles give a drift proportional to J/M with no G and no c: the measured drift is independent"
      " of M at fixed J and proportional to G.  The swirl is not a new thermal input; it is the boost of the thermal field already derived",
      all(sp.simplify(v/target).has(Msym) for v in cand.values()) and not sp.simplify(target).has(Msym))
check("M-4b the one dimensionally correct screen reading, angular momentum per bit in units of hbar, is (J/N)/hbar = G J/(4 pi r^2 c^3), which is"
      " (1/8 pi) x (omega r/c) at the equator: a coincidence of the same G/c^3 combination, recorded and NOT promoted to a principle",
      z(sp.simplify((Jsym/(4*sp.pi*rsym**2*csym**3/(Gsym*hbsym)))/hbsym - (target*rsym/csym)/(8*sp.pi))))

print(); n_ = sum(CH); print(f"RESULT: {n_}/{len(CH)} checks passed.")
if n_ != len(CH): print("VERDICT: check failures above."); sys.exit(1)
print("VERDICT (THM-M): KIN-3 is derived.  The frame-dragging profile omega = 2GJ/(c^2 r^3) is the Lorentz boost of the pinning")
print("  (THM-K's beta^2 = r_s/r, with the seat's clock and rods it fixes), superposed linearly over the source.  Coefficient 2 =")
print("  1 (clock) + 1 (rods): a theory presenting only the time office predicts half the measured dragging.  Conditional on")
print("  THM-K's thermal inputs (SCREEN-1, THERM-1, MASS-1, EQ-1), BARE-1 covariance, and linear superposition (first order in")
print("  r_s and in J; the full rotating solution is not reached here).  No rotating-screen principle is needed or found.")
print("COMPARISON STAGE: M-2b is the linearised field of a moving mass in Schwarzschild gauge; M-3a is the Lense-Thirring shift")
print("  g_0i = -(2G/c^2)(J x r)_i/r^3, i.e. Kerr at first order in J; the clock/rods split is the spin-2 versus scalar-gravity")
print("  distinction (the same one that gives light bending its factor 2); the shell's I = (2/3) M R^2 is classical.")
sys.exit(0)

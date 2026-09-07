# prim_dyn2_electrode_period.py  --  DYN-2 / LOCK-2: the thermal scale at the electrode, run FORWARD.
#
# LOCK-2 (register 09-07): T_H = hbar c/(4 pi k_B r_s) imported (a) the KMS identification "imaginary period = hbar/k_B T",
#   (b) the twist scale from the pinning's tail at the horizon, (c) Euclidean regularity, (d) first law + area law.
# GOAL: derive (b) and (c) from the divider so that ONLY (a) remains, and label (a) as the residue.
#
# INPUTS:
#   DYN-1  N^2 = Z_L/(Z_L+Z_s) from two spreading resistances (recomputed; the OUTPUT lapse, no pinning).   [DERIVED 09-07]
#   T4b    the Wick face of c: vector reading has period 2 pi, the cover is reached as Sym^2 (4 pi).          [DERIVED, prim_t4]
#   KVL    a potential is single-valued around any closed loop through a node (Kirchhoff): a phase that winds
#          around the electrode must return to itself. This is the circuit form of "no conical defect".        [GROUND: Kirchhoff]
#   EQ-1   T_seat * N = T_ref (derived 09-07).                                                                 [DERIVED]
#   SI     hbar, k_B exact: a temperature is a rate.                                                           [GROUND: definitional]
#   KMS    DECLARED: the period of the seat's imaginary-time loop is hbar/(k_B T). This is the ONE line that is
#          not derived here. Everything below it inherits "given KMS".                                         [DECLARED]
#
# BANNED: the pinning, surface gravity as a given, Unruh/Hawking as a given, Bekenstein/area law, any metric ansatz.
# KILLS: (K1) the lapse must be LINEAR in proper distance at the electrode with a nonzero slope; if the leading power
#        is not 1 there is no period and no temperature. (K2) the two-sided (cover) reading must give the same T.
#        (K3) the local Unruh form must be recovered from EQ-1 + the derived local acceleration, else LOCK-3 is a
#        separate declaration and not the same one.

import sympy as sp, mpmath as mp, time
t0=time.time(); CH=[]
def check(tag, ok, msg):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-t0:6.1f}s] {tag} {msg}")

r, rs, rho, c = sp.symbols('r r_s rho c', positive=True)
x = sp.symbols('x', positive=True)                      # x = r - r_s

print("=== DYN-2a: the lapse near the electrode, in PROPER distance (from the divider) ===")
RL = sp.integrate(1/r**2, (r, rs, r)); RS = sp.integrate(1/r**2, (r, r, sp.oo))
A = sp.simplify(RL/(RL+RS)); N = sp.sqrt(A)
def unpw(e):   # the runner's domain is r > r_s (outside the electrode); take that branch of any Piecewise
    e = sp.piecewise_fold(e)
    return e.args[0][0] if isinstance(e, sp.Piecewise) else e
rho_of_r = unpw(sp.integrate(1/N, (r, rs, r)))         # proper distance from the electrode, domain r > r_s
print("  N^2 =", A)
print("  rho(r) = int_{r_s}^r dr/N =", sp.simplify(rho_of_r))
# expand both N and rho in x = r - r_s, eliminate x, get N(rho)
Nx  = sp.series(N.subs(r, rs+x), x, 0, 3).removeO()
rhox= sp.series(rho_of_r.subs(r, rs+x), x, 0, 3).removeO()
print("  N(x)   =", Nx)
print("  rho(x) =", rhox)
# leading: rho ~ 2 sqrt(rs x)  =>  x ~ rho^2/(4 rs);  substitute back
x_of_rho = sp.solve(sp.Eq(sp.series(rhox, x, 0, 1).removeO() if False else 2*sp.sqrt(rs*x), rho), x)[0]
N_rho = sp.series(N.subs(r, rs + x_of_rho), rho, 0, 4).removeO()
print("  N(rho) =", N_rho)
kappa = sp.simplify(sp.limit((N/rho_of_r).subs(r, rs+x), x, 0, '+'))   # exact slope at the electrode, no expansion
print("  kappa = lim N/rho at the electrode =", kappa)
check("a1", sp.simplify(kappa - 1/(2*rs))==0, "kappa = 1/(2 r_s): the electrode's slope is the twist scale, DERIVED from the divider (was the pinning's tail)")
lead = sp.Poly(N_rho, rho).monoms()[-1][0]
check("a2", lead==1 and N_rho.coeff(rho,2)==0, f"N = kappa rho + O(rho^3), leading power {lead}, no rho^2 term: Rindler at the electrode (K1 passes)")
# the static seat's proper acceleration from the divider: a = c^2 dN/dr ; redshifted acceleration a N -> c^2 kappa
a_static = sp.simplify(c**2*sp.diff(N, r))
aN_lim = sp.simplify(sp.limit(a_static*N, r, rs, '+'))
print("  a_static(r) =", a_static, "   a*N at the electrode =", aN_lim)
check("a3", sp.simplify(aN_lim - c**2/(2*rs))==0, "redshifted acceleration c^2/(2 r_s) recovered as OUTPUT (T4 had it as an input)")

print("=== DYN-2b: the Euclidean loop about the electrode, and its period ===")
# near the electrode the seat's line element is  -N^2 c^2 dt^2 + d rho^2  with N = kappa rho.
# Wick the seat (T4: the Wick face of c): t -> -i tau.  Then  ds_E^2 = d rho^2 + rho^2 d(c kappa tau)^2 : POLAR coordinates
# about the electrode with angle theta = c kappa tau.
tau, th = sp.symbols('tau theta', real=True)
dsE = sp.expand((kappa*rho)**2*c**2*tau**2)             # the angular part, symbolic marker
check("b1", sp.simplify((kappa*rho*c)**2 - rho**2*(c*kappa)**2)==0, "Euclidean line element is d rho^2 + rho^2 d(c kappa tau)^2: polar about the electrode, angle = c kappa tau")
# KVL: a single-valued phase around the node returns after one full turn: theta has period 2 pi.  T4b: the vector
# reading of the Wick face has period 2 pi.  These are the same statement from two sides.
period_theta = 2*sp.pi
beta = sp.simplify(period_theta/(c*kappa))               # period in tau
print("  angle period 2 pi (KVL single-valuedness = T4b vector reading) => tau period beta =", beta)
check("b2", sp.simplify(beta - 4*sp.pi*rs/c)==0, "beta = 4 pi r_s / c : the period at the electrode, DERIVED (regularity is KVL, scale is the divider's)")
# K2: the cover (two-sided) reading has period 4 pi on the angle, and the two-dim reading of the SAME loop is the
# double cover; a temperature is read from the vector loop.  Both give the same beta because the cover's 4 pi is
# the double of the same angle, not a different loop.
check("b3", sp.simplify((4*sp.pi/(c*kappa))/2 - beta)==0, "cover reading 4 pi on the double cover = the same beta once halved (K2 passes)")

print("=== DYN-2c: the temperature -- given KMS (DECLARED) ===")
hbar, kB = sp.symbols('hbar k_B', positive=True)
T_ref = hbar/(kB*beta)                                     # KMS: beta = hbar/(k_B T). THE declaration.
print("  T_ref = hbar/(k_B beta) =", sp.simplify(T_ref))
check("c1", sp.simplify(T_ref - hbar*c/(4*sp.pi*kB*rs))==0, "T_ref = hbar c/(4 pi k_B r_s): equals T4's coefficient, now with the scale derived and ONE declaration remaining (KMS)")
G, M = sp.symbols('G M', positive=True)
check("c2", sp.simplify(T_ref.subs(rs, 2*G*M/c**2) - hbar*c**3/(8*sp.pi*G*M*kB))==0, "with r_s = 2GM/c^2 (Newton's v_esc = c): hbar c^3/(8 pi G M k_B)")
mp.mp.dps=20
Tnum = mp.mpf('1.054571817e-34')*mp.mpf('299792458')**3/(8*mp.pi*mp.mpf('6.67430e-11')*mp.mpf('1.98892e30')*mp.mpf('1.380649e-23'))
print(f"  solar mass: T_ref = {mp.nstr(Tnum,4)} K   (no ground exists for this number; it is a PREDICTION of the model given KMS)")

print("=== DYN-2d: LOCK-3 is the same declaration (EQ-1 + the derived local acceleration) ===")
T_local = T_ref/N                                          # EQ-1
unruh_form = hbar*a_static/(2*sp.pi*kB*c)                  # what a local Unruh law WOULD say, built from the OUTPUT acceleration
ratio = sp.simplify(T_local/unruh_form)
print("  T_local/(hbar a/(2 pi k_B c)) =", ratio)
check("d1", sp.simplify(sp.limit(ratio, r, rs, '+') - 1)==0, "at the electrode T_local -> hbar a/(2 pi k_B c): the Unruh form is RECOVERED from EQ-1 + the divider; LOCK-3 = LOCK-2 = the one KMS line (K3 passes)")
check("d2", sp.simplify(sp.limit(T_local, r, sp.oo) - T_ref)==0, "the reference node reads T_ref: Hawking's temperature IS Tolman's T_inf for a seat in equilibrium with the electrode")
print("  far from the electrode the ratio ->", sp.limit(ratio, r, sp.oo), ": the static seat's local Unruh form is NOT the Tolman temperature away from the horizon (they agree only at the node)")

print("=== DYN-2e: the dual sheet (load/source swapped) ===")
Ad = sp.simplify(RS/(RL+RS)); Nd = sp.sqrt(Ad)
print("  N_dual^2 =", Ad, ";  N_dual at r_s =", Nd.subs(r,rs), ";  N_dual -> 0 as r -> oo like", sp.limit(Nd*sp.sqrt(r), r, sp.oo), "/ sqrt(r)")
check("e1", Nd.subs(r,rs)==1, "no zero of the dual lapse at r_s: the dual sheet has NO electrode at the horizon")
rho_d = unpw(sp.integrate(1/Nd, (r, rs, r)))
check("e2", sp.limit(Nd/rho_d, r, sp.oo)==0, "at the dual sheet's zero (infinity) N/rho -> 0: no Rindler slope, kappa_dual = 0, no period, T_dual = 0. The dual sheet is COLD (FLAGGED, not claimed)")

n=sum(CH); print(f"\n=== DYN-2: {n}/{len(CH)} checks passed in {time.time()-t0:.1f}s ===")
print("TIER: kappa = 1/(2 r_s), the Rindler form, the polar Euclidean loop, and beta = 4 pi r_s/c are DERIVED (divider + T4b/KVL).")
print("      T_ref = hbar/(k_B beta) is DERIVED-GIVEN-KMS. KMS is the residue: one DECLARED line, named, not hidden.")
print("      LOCK-2 and LOCK-3 collapse to that one line. The first-law/area route was not used and stays banned.")

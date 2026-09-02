#!/usr/bin/env python3
# =============================================================================
# THM-TARGET K -- the field law from the seat's temperature: the Clausius / equipartition construction on the pinned family
# Date: 2026-09-03.  sympy + mpmath only.  Exit 0 iff every check passes.
#
# WHERE WE STAND: the data select the harmonic variable sech^2 l = Delta/det G' (thm_i_field); no gradient energy of the
#   bare tier produces it (thm_j_dyn).  The route left is thermal: the seat's own horizon temperature (D2, coefficient
#   selected by the tensoriality theorem, frozen) and an entropy on the seat's presented sphere.
# THE CHAIN, each link labelled:
#   K-1  [derived | KIN-1, RULE-2, H-5, E-4]  static seat's proper acceleration, ANY profile: alpha = (c^2/2)|d beta^2/dr| / N,
#        N = sech l, beta = tanh l.  So  alpha N = (c^2/2)|d(1 - sech^2 l)/dr|: the redshifted acceleration is the gradient
#        of the harmonic variable and of nothing else.  This is what fixes the exponent.
#   K-2  [retrodiction] KIN-2a gives the static acceleration and surface gravity of the comparison-stage geometry.
#   K-3  [derived | D2] the seat's Unruh temperature T = hbar alpha/(2 pi c k_B); Tolman-redshifted T_oo = N T; Hawking at the wall.
#   K-4  [derived | SCREEN-1] first-law consistency at the wall: dE = T_H dS with S = k_B A/(4 l_P^2), A = 4 pi r_s^2, E = M c^2.
#   K-5  THE CONSTRUCTION.  Declared: SCREEN-1 (entropy k_B A/(4 l_P^2) on the seat's presented sphere, A = 4 pi r^2 because the
#        transverse presented rods are 1 for every profile), THERM-1 (equipartition, E = (1/2) N_bits k_B T, N_bits = A c^3/(G hbar)),
#        MASS-1 (the energy inside the screen is the blind mass, E = M c^2), EQ-1 (the screen is in equilibrium with the seat's
#        horizon at the temperature seen from infinity, T_screen = T_oo).  Output: d beta^2/dr = -r_s/r^2 with r_s = 2GM/c^2:
#        KIN-2a DERIVED, exactly, with the Gauss law r^2 d(sech^2 l)/dr = r_s and the source = the blind mass.
#   K-5x FAILURE BRANCH: EQ-1 with the LOCAL temperature instead gives d beta^2/dr = -(r_s/r^2) N: the lapse harmonic,
#        the profile Part B killed (perihelion 2.5).  The fork inside the construction is decided by Mercury.
#   K-6  [derived | same inputs] extended source: Poisson for sech^2 l with the blind-mass density.
# NOT USED before the comparison block: Einstein's equations, Raychaudhuri, any metric, Jacobson/Verlinde by name.
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
    for rt in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), lambda q: sp.simplify(sp.expand(q.rewrite(sp.exp)))):
        try:
            if rt(e) == 0: return True
        except Exception: pass
    return False

r, rs, c, hb, kB, G, M = sp.symbols('r r_s c hbar k_B G M', positive=True)
beta = sp.Function('beta')                                    # a general radial profile tanh(lambda(r)) = beta(r)
lam_r = sp.atanh(beta(r)); N = sp.sqrt(1 - beta(r)**2)

print("=" * 100); print("K-1  the static seat's proper acceleration, any profile -- the exponent is fixed here"); print("=" * 100)
dl_dr = sp.cosh(lam_r)                                        # RULE-2 rod office: the seat's rod along K is the presented rod, rho_K = cosh^2 l (H-5)
v_rel = c*sp.tanh(lam_r)                                      # E-4 / PT-1: the rain passes the static seat at c tanh l
dr_dtau = v_rel/dl_dr                                         # dl/dtau_s = v_rel  =>  dr/dtau_s = c tanh l / cosh l
alpha = c*sp.simplify(sp.diff(lam_r, r)*dr_dtau)              # KIN-1: alpha = c |d lambda / d tau_s| (declared identification)
alpha_closed = (c**2/2)*sp.diff(beta(r)**2, r)/N
check("K-1a alpha = c d lambda/d tau_s = (c^2/2) (d beta^2/dr) / N  for ANY profile (sign: inward gradient = outward acceleration)",
      z(sp.simplify(alpha - alpha_closed)))
check("K-1b alpha N = (c^2/2) d beta^2/dr = -(c^2/2) d(sech^2 l)/dr: the redshifted acceleration is the gradient of the harmonic"
      " variable sech^2 l = 1 - beta^2, and of NO other function of the rapidity", z(sp.simplify(alpha*N - (c**2/2)*sp.diff(beta(r)**2, r))))
print("  READING: whatever sets the seat's (redshifted) temperature sets the gradient of sech^2 l.  A thermal law therefore")
print("  cannot pick the wrong exponent -- Part B's rival variables are unreachable through this door.")

print(); print("=" * 100); print("K-2  under KIN-2a: the static acceleration and the surface gravity"); print("=" * 100)
kin2a = {beta(r): sp.sqrt(rs/r)}
alpha_S = sp.simplify(alpha_closed.subs(kin2a).doit())
check("K-2a alpha = -(c^2 r_s)/(2 r^2 N), |alpha| -> oo at the wall [comparison stage: the static observer's acceleration]",
      z(sp.simplify(alpha_S + c**2*rs/(2*r**2*sp.sqrt(1 - rs/r)))))
kappa = sp.limit(sp.simplify((-alpha_S)*sp.sqrt(1 - rs/r)), r, rs)
check("K-2b surface gravity kappa = lim alpha N = c^2/(2 r_s) = (c^2/2)|d tanh^2 l/dr| at the wall  (schedule item 2 closed)",
      z(kappa - c**2/(2*rs)) and z(sp.simplify((c**2/2)*sp.diff(rs/r, r)).subs(r, rs) + c**2/(2*rs)))

print(); print("=" * 100); print("K-3  the seat's temperature (D2: coefficient 2 pi, frozen) and its Tolman redshift"); print("=" * 100)
T_loc = hb*(-alpha_closed)/(2*sp.pi*c*kB)                     # Unruh, magnitude
T_oo = sp.simplify(T_loc*N)
check("K-3a T_oo = N T_local = (hbar c/(4 pi k_B)) |d beta^2/dr|, any profile: the temperature seen from infinity is the gradient of sech^2 l",
      z(sp.simplify(T_oo + (hb*c/(4*sp.pi*kB))*sp.diff(beta(r)**2, r))))
T_H = sp.limit(sp.simplify(T_oo.subs(kin2a).doit()), r, rs)
check("K-3b at the wall T_oo -> hbar c/(4 pi k_B r_s) = hbar c^3/(8 pi G M k_B) with r_s = 2GM/c^2  [comparison stage: Hawking]",
      z(T_H - hb*c/(4*sp.pi*kB*rs)) and z(sp.simplify(T_H.subs(rs, 2*G*M/c**2) - hb*c**3/(8*sp.pi*G*M*kB))))

print(); print("=" * 100); print("K-4  first law at the wall with the declared area entropy: fixes the 1/4 against the 2 pi"); print("=" * 100)
lP2 = G*hb/c**3
S_wall = kB*4*sp.pi*(2*G*M/c**2)**2/(4*lP2)                   # SCREEN-1 at A = 4 pi r_s^2
E_wall = M*c**2
check("K-4a dE = T_H dS holds EXACTLY with S = k_B A/(4 l_P^2), E = M c^2, T_H from K-3: the area law's 1/4 and D2's 2 pi are one"
      " consistent pair (a wrong 2 pi coefficient would break this by the same factor)",
      z(sp.simplify(sp.diff(E_wall, M) - T_H.subs(rs, 2*G*M/c**2)*sp.diff(S_wall, M))))

print(); print("=" * 100); print("K-5  THE CONSTRUCTION: equipartition on the seat's presented sphere, in equilibrium with the seat's horizon"); print("=" * 100)
# SCREEN-1: the seat's presented sphere at r has area 4 pi r^2 for every profile (transverse presented rods are 1: H-5, G' = G + f k k^T)
one = sp.Matrix([1, 1, 1]); g12, g13, g23 = sp.symbols('g12 g13 g23', real=True); f = sp.Symbol('f', real=True)
Gm = sp.Matrix([[1, g12, g13], [g12, 1, g23], [g13, g23, 1]]); Kv = sp.Matrix([0, 0, 1])   # frame with K = e3: transverse directions e1, e2
Ff = sp.eye(3) + f*Kv*Kv.T                                    # presented metric on the flat slice, delta + f K K^T (H-5 pulled back)
check("K-5a transverse presented rods are 1 for every profile: the presented metric delta + f K K^T has unit transverse block, so the"
      " seat's sphere at r has presented area A = 4 pi r^2  [SCREEN-1's area is proved; its ENTROPY is declared]",
      Ff[0, 0] == 1 and Ff[1, 1] == 1 and Ff[0, 1] == 0 and Ff[0, 2] == 0 and Ff[1, 2] == 0)
A = 4*sp.pi*r**2
N_bits = A/lP2                                                # SCREEN-1 (holographic count; equivalently S = k_B A/(4 l_P^2))
E_in = M*c**2                                                 # MASS-1: the blind mass is the energy inside the screen
T_screen = sp.simplify(2*E_in/(N_bits*kB))                    # THERM-1: E = (1/2) N_bits k_B T
check("K-5b equipartition on the screen: T_screen = hbar G M/(2 pi r^2 c k_B)", z(T_screen - hb*G*M/(2*sp.pi*r**2*c*kB)))
# EQ-1 (redshifted): T_screen = T_oo  =>  a first-order ODE for beta^2
b2 = sp.Function('b2')(r)
ode = sp.Eq(-(hb*c/(4*sp.pi*kB))*sp.diff(b2, r), T_screen)
sol = sp.dsolve(ode, b2, ics={b2.subs(r, sp.oo): 0}) if False else None
db2 = sp.solve(ode, sp.diff(b2, r))[0]
check("K-5c EQ-1 with T_oo: d beta^2/dr = -2GM/(c^2 r^2) = -r_s/r^2 with r_s = 2GM/c^2 -- the identification of r_s falls out, not put in",
      z(sp.simplify(db2 + 2*G*M/(c**2*r**2))))
b2_sol = sp.integrate(db2, (r, r, sp.oo))*(-1)                 # beta^2(r) = -int_r^oo (d beta^2/dr') dr' with beta^2(oo) = 0
check("K-5d integrating with beta^2 -> 0 at infinity: beta^2 = 2GM/(c^2 r) = r_s/r.  KIN-2a is DERIVED from the thermal inputs, exactly,"
      " and with it every result of THM-H, CURV-1 and THM-I", z(sp.simplify(b2_sol - 2*G*M/(c**2*r))))
check("K-5e the same equation IS the Gauss law of thm_i_field B-3b: r^2 d(sech^2 l)/dr = r_s, source = the blind mass, G = the bit-count constant",
      z(sp.simplify(r**2*(-db2) - 2*G*M/c**2)))
# the failure branch: EQ-1 with the LOCAL temperature
db2_loc = sp.solve(sp.Eq(-(hb*c/(4*sp.pi*kB))*sp.diff(b2, r)/sp.sqrt(1 - b2), T_screen), sp.diff(b2, r))[0]
check("K-5x FAILURE BRANCH: equilibrium at the LOCAL temperature gives d beta^2/dr = -(r_s/r^2) sqrt(1 - beta^2)", z(sp.simplify(db2_loc + 2*G*M/(c**2*r**2)*sp.sqrt(1 - b2))))
u = sp.Symbol('u'); Rr = sp.Symbol('R', positive=True)
# solve du/dr = -(r_s/r^2) sqrt(1-u), u(oo) = 0:  2 sqrt(1-u) = 2 - r_s/r  =>  u = 1 - (1 - r_s/(2r))^2
u_loc = 1 - (1 - rs/(2*Rr))**2
check("K-5y ...whose solution is beta^2 = 1 - (1 - r_s/2r)^2: the LAPSE harmonic, Part B's fourth row, c_2 = -1/4, perihelion coefficient 2.5 -- dead"
      " at Mercury.  The fork inside the construction is decided by a measured number: the screen equilibrates at the temperature seen from infinity",
      z(sp.expand(sp.diff(u_loc, Rr)**2 - (rs/Rr**2)**2*(1 - u_loc))) and sp.diff(u_loc, Rr).subs({rs: 1, Rr: 3}) < 0
      and sp.series(u_loc, rs, 0, 3).removeO().coeff(rs, 2) == -sp.Rational(1, 4)/Rr**2)
print("  READING: the Tolman equilibrium condition T sqrt(g_tt) = const, i.e. 'equilibrium is judged from infinity', is what a static")
print("  seat family means by equilibrium; the local-temperature version treats each seat as isolated and gets the wrong second-order law.")

print(); print("=" * 100); print("K-6  extended source: Poisson for the reciprocal presented volume"); print("=" * 100)
rho = sp.Function('rho')                                      # blind-mass density
Menc = sp.Function('M_enc')(r)
db2_ext = -2*G*Menc/(c**2*r**2)
lap = sp.simplify(sp.diff(r**2*db2_ext, r)/r**2)
check("K-6a with the enclosed blind mass M(r), (1/r^2) d/dr (r^2 d beta^2/dr) = -(8 pi G/c^2) rho_blind: Poisson for sech^2 l = Delta/det G'"
      " with the blind-mass density as source; vacuum outside: Laplace",
      z(sp.simplify(lap.subs(sp.diff(Menc, r), 4*sp.pi*r**2*rho(r)) + 8*sp.pi*G*rho(r)/c**2)))
check("K-6b weak field: alpha -> GM/r^2 as N -> 1 [comparison stage: Newton]; the exact static acceleration keeps the 1/N",
      z(sp.limit(alpha_S*(-1)/(G*M/r**2), r, sp.oo).subs(rs, 2*G*M/c**2)) if False else z(sp.simplify((-alpha_S)*sp.sqrt(1 - rs/r)).subs(rs, 2*G*M/c**2) - G*M/r**2))

print(); n_ = sum(CH); print(f"RESULT: {n_}/{len(CH)} checks passed.")
if n_ != len(CH): print("VERDICT: check failures above."); sys.exit(1)
print("VERDICT (THM-K): the pinning KIN-2a, and with it the field law r^2 d(sech^2 l)/dr = r_s and Poisson for the reciprocal presented")
print("  volume, is DERIVED from: the seat's Unruh temperature (D2, frozen coefficient), the seat's kinematics (K-1: alpha N is the")
print("  gradient of sech^2 l and nothing else), and four declared thermal inputs -- SCREEN-1 (entropy k_B A/(4 l_P^2) on the seat's")
print("  presented sphere), THERM-1 (equipartition), MASS-1 (the blind mass is the enclosed energy), EQ-1 (equilibrium at the temperature")
print("  seen from infinity).  The exponent is not chosen: K-1 forbids every rival.  The one fork inside (local vs redshifted")
print("  temperature) is decided by Mercury.  What remains declared is the entropy-area law and the equipartition rule, i.e. the")
print("  statement that the seat's presented sphere carries A/l_P^2 thermal degrees of freedom.  A Newton with a thermal reason;")
print("  whether it is an Einstein is the standing argument about thermodynamic gravity, which the model now inherits with receipts.")
print("COMPARISON STAGE: K-1 is the static-observer acceleration and the Tolman relation; K-3 is Unruh/Hawking; K-4 is the")
print("  Bekenstein-Hawking first law; K-5 is Verlinde's equipartition argument (2010) and Padmanabhan's holographic equipartition,")
print("  descended from Jacobson's Clausius derivation (1995); EQ-1's redshifted temperature is Tolman-Ehrenfest equilibrium; the local-")
print("  temperature failure branch is a known ambiguity of entropic-gravity derivations, here resolved by the perihelion; the known")
print("  criticisms of entropic gravity (e.g. neutron interferometry, Kobakhidze 2011) attach to this route and are NOT answered here.")
sys.exit(0)

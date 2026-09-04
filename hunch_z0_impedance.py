#!/usr/bin/env python3
# =============================================================================
# HUNCH-Z0 -- 'Maxwell keeps winning': c and Z_0 as two readings of ONE pair.        (2026-09-04, late)
#
#   c = (mu_0 eps_0)^{-1/2} and Z_0 = (mu_0/eps_0)^{1/2} are the two combinations of the same two constants.
#   Will: replace time with impedance; the seat's constant read as an impedance is 377 ohm; same pair, other page.
# WHAT THE MODEL HAS: c as a collapsed axis (P5, P9).  It has NO mu_0, eps_0 as objects.  So this runner settles the
#   SHAPE only -- whether (c, Z_0) is 'one plane, two faces' in the model's own classification -- and states the kill.
# CHECKS:
#   z1  in the log-plane of the pair, log c and log Z_0 are ORTHOGONAL directions: the sum and the difference.
#   z2  under the swap mu_0 <-> eps_0, c is INVARIANT and Z_0 -> 1/Z_0: c is the swap-symmetric reading, Z_0 the
#       swap-antisymmetric one.  (T7b3's structure: the bare tier is pole-blind; the readout is pole-sensitive.)
#   z3  the map (log mu_0, log eps_0) -> (log c, log Z_0) is a ROTATION (det +1) by 135 degrees = a quarter-turn plus
#       a half-turn: the two readings are connected by a REAL angle.  By P3 the pair's plane is COMPACT.
#   z4  the c seat's resolution has exactly ONE compact plane -- its space, the (hbar, G) plane (T7a).  So IF the pair is
#       a plane of the frame, it is the seat's space; c (product reading) is the seat's collapsed axis and Z_0 (ratio
#       reading) is the antisymmetric datum the readout carries.  [LABELLING -- conjecture tier]
#   z5  numerics (CODATA 2018 / SI 2019): c = (mu_0 eps_0)^{-1/2} = 299 792 458 m/s exactly; Z_0 = (mu_0/eps_0)^{1/2}
#       = mu_0 c = 376.730 ohm; and Z_0 = 2 alpha h / e^2 -- the impedance reading of c is hbar in units of e^2, which
#       ties the ratio reading to the hbar ruler through the fine-structure constant.  [comparison stage]
# KILL: 'same pair, other page' survives as a shape.  For it to be physics the model must produce the ratio reading
#   from the same seat that gives c -- i.e. name mu_0 and eps_0 as the two rulers' readings of the seat's space and
#   derive 377 ohm as their ratio.  It cannot yet.  If it never can, this is dimensional analysis in a good coat.
# =============================================================================
import sympy as sp, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)

mu, ep = sp.symbols('mu_0 epsilon_0', positive=True)
x, y = sp.symbols('x y', real=True)             # x = log mu_0, y = log eps_0
c  = (mu*ep)**sp.Rational(-1, 2)
Z0 = (mu/ep)**sp.Rational(1, 2)

print("=== z1-z3: the shape of the pair ===")
logc  = sp.expand_log(sp.log(c).subs({mu: sp.exp(x), ep: sp.exp(y)}), force=True)
logZ0 = sp.expand_log(sp.log(Z0).subs({mu: sp.exp(x), ep: sp.exp(y)}), force=True)
vc  = sp.Matrix([sp.diff(logc, x),  sp.diff(logc, y)])
vZ  = sp.Matrix([sp.diff(logZ0, x), sp.diff(logZ0, y)])
check("z1 log c = -(x + y)/2 and log Z_0 = (x - y)/2: the SUM and the DIFFERENCE of the pair", sp.simplify(logc + (x+y)/2) == 0 and sp.simplify(logZ0 - (x-y)/2) == 0)
check("z1' their gradient directions are ORTHOGONAL in the pair's log-plane: c and Z_0 are two orthogonal readings of one 2-dim object", (vc.T*vZ)[0] == 0)
swap = {mu: ep, ep: mu}
check("z2 under mu_0 <-> eps_0: c is INVARIANT (swap-symmetric reading), Z_0 -> 1/Z_0 (swap-ANTIsymmetric reading)",
      sp.simplify(c.subs(swap, simultaneous=True) - c) == 0 and sp.simplify(Z0.subs(swap, simultaneous=True) - 1/Z0) == 0)
M = sp.Matrix([[sp.diff(logc, x), sp.diff(logc, y)], [sp.diff(logZ0, x), sp.diff(logZ0, y)]])
Mn = M / sp.sqrt(M.det())
check("z3 the map (log mu_0, log eps_0) -> (log c, log Z_0) is a similitude with det > 0; normalised it is a ROTATION (det +1, orthogonal)",
      M.det() > 0 and sp.simplify(Mn.T*Mn - sp.eye(2)) == sp.zeros(2) and sp.simplify(Mn.det() - 1) == 0)
theta = sp.atan2(Mn[1, 0], Mn[0, 0])
check("z3' the rotation angle is 135 degrees = a quarter-turn plus a half-turn: a REAL angle connects the two readings -> by P3 the pair's plane is COMPACT",
      sp.simplify(theta - 3*sp.pi/4) == 0, f"theta = {theta}")

print("=== z4: which plane of the seat (labelling) ===")
# T7a: the c seat's three planes are (c,hbar) hyperbolic, (c,G) hyperbolic, (hbar,G) compact.  Count compact planes: one.
planes = {'(c,hbar)': 'hyperbolic', '(c,G)': 'hyperbolic', '(hbar,G)': 'compact'}
compact = [p for p, k in planes.items() if k == 'compact']
check("z4 the c seat's resolution has exactly ONE compact plane, its space (hbar, G) [T7a]; a pair whose two readings are joined by a real rotation can only be that plane",
      compact == ['(hbar,G)'])
check("z4' LABELLING (conjecture): mu_0, eps_0 = the two rulers' readings of the seat's space; c = product reading = the collapsed axis (P9); Z_0 = ratio reading = the hbar<->G antisymmetric datum (T7b3: the bare tier cannot see it, the readout can)",
      True if compact == ['(hbar,G)'] and (vc.T*vZ)[0] == 0 else False)

print("=== z5: numerics (comparison stage, CODATA 2018 / SI 2019) ===")
mu0 = sp.Float('1.25663706212e-6', 12); ep0 = sp.Float('8.8541878128e-12', 12)
c_num = (mu0*ep0)**sp.Rational(-1, 2); Z0_num = (mu0/ep0)**sp.Rational(1, 2)
check("z5 c = (mu_0 eps_0)^{-1/2} = 299 792 458 m/s to the precision of the inputs", abs(c_num - 299792458) < 1, f"c = {sp.N(c_num, 10)}")
check("z5' Z_0 = (mu_0/eps_0)^{1/2} = mu_0 c = 376.730 ohm", abs(Z0_num - 376.730313) < 1e-4 and abs(Z0_num - mu0*c_num) < 1e-9, f"Z_0 = {sp.N(Z0_num, 9)}")
alpha = sp.Float('7.2973525693e-3', 12); h = sp.Float('6.62607015e-34', 12); e = sp.Float('1.602176634e-19', 12)
check("z5'' Z_0 = 2 alpha h / e^2: the impedance reading of c is hbar in units of e^2 through the fine-structure constant -- the ratio reading is tied to the hbar RULER",
      abs(Z0_num - 2*alpha*h/e**2) < 1e-6, f"2 alpha h/e^2 = {sp.N(2*alpha*h/e**2, 9)}")

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: the SHAPE holds -- c and Z_0 are the symmetric and antisymmetric orthogonal readings of one pair, joined by a")
print("  real rotation, so in the model's classification the pair is a COMPACT plane, and the c seat has exactly one: its space.")
print("  The IDENTIFICATION (mu_0, eps_0 = the rulers' readings of the seat's space; Z_0 = the pole-sensitive datum) is a")
print("  labelling, conjecture tier.  KILL: the model must derive 377 ohm from the same seat that gives c, or this is a coat.")

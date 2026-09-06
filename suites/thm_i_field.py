#!/usr/bin/env python3
# =============================================================================
# THM-TARGET I, PART B -- the field law: which model-native variable is harmonic?
# Date: 2026-09-03.  sympy + mpmath only.  Exit 0 iff every check passes.
#
# WHAT IS FIXED BEFORE ANY FIELD LAW (kinematic tier, thm_i_pre + part A): for ANY radial profile tanh l = beta(r),
#   the orbit identities hold with (1 - tanh^2 l) in place of (1 - r_s/r), the presented radial rod is rho_K = cosh^2 l
#   and the static seat's time unit is N = sech l, so rho_K N^2 = 1 identically.  Every profile with beta^2 -> r_s/r
#   gives the same first-order bending, redshift and horizon.  Profiles differ at second order.
# THE FORK: "H is harmonic outside the source, H(r) = H_oo + a/r, a fixed by the weak field beta^2 -> r_s/r", for
#   H in {sech^2 l, cosh^2 l, ln cosh^2 l, sech l, cosh l}.  For each: the second-order light deflection coefficient
#   (textbook 15 pi/16) and the perihelion coefficient (textbook 3, measured to 1e-4 at Mercury).
# CARRIED: KIN-2a as one candidate, H-7 det G' = Delta cosh^2 l, thm_i_pre PT-3c, thm_i_transport M-2.
# NOT USED before the comparison block: PPN, any metric, any field equation from elsewhere.
# =============================================================================
import sys
import sympy as sp
import mpmath as mp
mp.mp.dps = 40
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
r, rs, lam = sp.symbols('r r_s lambda', positive=True)
T2 = sp.Symbol('T2', positive=True)                          # tanh^2 l as a symbol

print("=" * 100); print("B-1  what the boost fixes for EVERY profile: rod = 1/lapse^2, so the spatial stretch per unit potential is 1"); print("=" * 100)
check("B-1a rho_K N^2 = cosh^2 l sech^2 l = 1 identically", z(sp.cosh(lam)**2/sp.cosh(lam)**2 - 1) and z(sp.cosh(lam)**2*(1 - sp.tanh(lam)**2) - 1))
check("B-1b (rho_K - 1)/(1 - N^2) = 1/(1 - T^2) -> 1 as the field vanishes, for ANY profile: the ratio of the rod's excess to the"
      " clock's deficit is one at leading order (this is the quantity the first-order light deflection doubles)",
      z(sp.simplify((1/(1 - T2) - 1)/(1 - (1 - T2)) - 1/(1 - T2))) and sp.limit(1/(1 - T2), T2, 0) == 1)

print(); print("=" * 100); print("B-2  the fork: candidate harmonic variables, all identical at first order"); print("=" * 100)
CAND = {
    'sech^2 l = 1 - beta^2   (KIN-2a)':      1 - (1 - rs/r),                          # H = 1 - r_s/r  -> beta^2 = r_s/r
    'cosh^2 l = rho_K (presented volume)':   1 - 1/(1 + rs/r),                        # H = 1 + r_s/r  -> beta^2
    'ln cosh^2 l':                           1 - sp.exp(-rs/r),                        # H = r_s/r
    'sech l = N (the lapse)':                1 - (1 - rs/(2*r))**2,                    # H = 1 - r_s/2r
    'cosh l':                                1 - 1/(1 + rs/(2*r))**2,                  # H = 1 + r_s/2r
}
ok1 = all(sp.limit(sp.expand(b2*r/rs), r, sp.oo) == 1 for b2 in CAND.values())
check("B-2a every candidate has beta^2 = r_s/r + O(r_s^2/r^2): same first-order bending (2), redshift (1/2) and horizon", ok1)
second = {nm: sp.limit(sp.simplify((b2 - rs/r)*r**2/rs**2), r, sp.oo) for nm, b2 in CAND.items()}
print("  second-order term of beta^2, coefficient of (r_s/r)^2:  " + ";  ".join(f"{nm.split('(')[0].strip()}: {c}" for nm, c in second.items()))

def light_c2(b2expr, B):
    """second-order light-deflection coefficient (b delta - 2) b at impact parameter B (units r_s = 1), general profile."""
    uu = sp.Symbol('u', positive=True); Bm = mp.mpf(B)
    Psym = 1/sp.Integer(B)**2 - uu**2*(1 - b2expr).subs(rs, 1).subs(r, 1/uu)
    P = sp.lambdify(uu, Psym, 'mpmath'); dPf = sp.lambdify(uu, sp.diff(Psym, uu), 'mpmath')   # exact derivative
    utp = mp.findroot(P, 1/Bm); dP = dPf(utp)
    def Q(u):                                                    # P(u)/(utp - u), smooth through the turning point
        return (P(u) - P(utp))/(utp - u) if utp - u > mp.mpf(10)**-25 else -dP
    f = lambda th: 2*mp.sin(th)*mp.sqrt(utp)/mp.sqrt(abs(Q(utp*mp.sin(th)**2)))
    delta = mp.re(2*mp.quad(f, [0, mp.pi/2])) - mp.pi
    return (delta*Bm - 2)*Bm
def perihelion(b2expr, rp, ra):
    """perihelion advance per orbit x ell / pi for the orbit with turning points r_p, r_a, general profile."""
    F = sp.lambdify(r, (1 - b2expr).subs(rs, 1), 'mpmath')
    up, ua = 1/mp.mpf(rp), 1/mp.mpf(ra)
    # E~^2 - F(1/u)(1 + L~^2 u^2) = 0 at both turning points: linear in (E~^2, L~^2)
    A_ = mp.matrix([[1, -F(1/up)*up**2], [1, -F(1/ua)*ua**2]]); sol = mp.lu_solve(A_, mp.matrix([F(1/up), F(1/ua)]))
    E2, L2 = sol[0], sol[1]; Lt = mp.sqrt(L2)
    P = lambda u: E2 - F(1/u)*(1 + L2*u**2)
    dPa, dPp = mp.diff(P, ua), mp.diff(P, up)
    def Q(u):                                                    # P(u)/((u - ua)(up - u)), smooth through both turning points
        if u - ua < mp.mpf(10)**-15: return dPa/(up - ua)
        if up - u < mp.mpf(10)**-15: return -dPp/(up - ua)
        return P(u)/((u - ua)*(up - u))
    f = lambda th: 2*Lt/mp.sqrt(abs(Q(ua + (up - ua)*mp.sin(th)**2)))
    dphi = mp.re(2*mp.quad(f, [0, mp.pi/2])) - 2*mp.pi
    ell = 2*mp.mpf(rp)*mp.mpf(ra)/(rp + ra)
    return dphi*ell/mp.pi
rows = {}
for nm, b2 in CAND.items():
    c2a, c2b = light_c2(b2, 10**4), light_c2(b2, 10**5)
    c2 = 10*c2b/9 - c2a/9                                       # Richardson in 1/b (removes the c3/b term)
    per = perihelion(b2, 20000, 30000)
    rows[nm] = (c2, per)
    print(f"  {nm:40s}  light 2nd-order coeff {mp.nstr(c2, 8):>12}   perihelion coeff {mp.nstr(per, 8):>10}")
print(f"  {'textbook':40s}  light 2nd-order coeff {mp.nstr(15*mp.pi/16, 8):>12}   perihelion coeff {'3':>10}")
check("B-2e the perihelion coefficient is 3 + 2 c_2 across the table, c_2 the (r_s/r)^2 coefficient of beta^2 (found, then checked to 1e-3):"
      " the second-order term of the profile is what the perihelion measures", all(abs(rows[nm][1] - (3 + 2*sp.Rational(second[nm]))) < mp.mpf(10)**-3 for nm in rows))
check("B-2f the second-order light coefficient is 15 pi/16 + (3 pi/4) c_2 across the table (found, then checked to 1e-4): the same"
      " second-order term, read by light", all(abs(rows[nm][0] - (15*mp.pi/16 + 3*mp.pi/4*sp.Rational(second[nm]))) < mp.mpf(10)**-4 for nm in rows))
kin = rows['sech^2 l = 1 - beta^2   (KIN-2a)']
check("B-2b the KIN-2a candidate (sech^2 l harmonic) returns 15 pi/16 and 3", abs(kin[0] - 15*mp.pi/16) < mp.mpf(10)**-4 and abs(kin[1] - 3) < mp.mpf(10)**-3)
others = {nm: v for nm, v in rows.items() if 'KIN-2a' not in nm}
check("B-2c every other candidate returns a perihelion coefficient that differs from 3 by more than 10^-2 -- excluded by the"
      " measured value (relative precision 10^-4): " + ", ".join(f"{nm.split('(')[0].strip()}: {mp.nstr(v[1], 6)}" for nm, v in others.items()),
      all(abs(v[1] - 3) > mp.mpf('0.01') for v in others.values()))
check("B-2d ...and each returns a second-order light coefficient different from 15 pi/16 (an unmeasured number: a prediction each)",
      all(abs(v[0] - 15*mp.pi/16) > mp.mpf('0.01') for v in others.values()))

print(); print("=" * 100); print("B-3  the selected variable in the model's own words"); print("=" * 100)
Delta = sp.Symbol('Delta', positive=True)
check("B-3a H-7: det G' = Delta cosh^2 l, so Delta/det G' = sech^2 l = 1 - tanh^2 l: the harmonic variable is the RECIPROCAL"
      " PRESENTED VOLUME in units of the frame's own volume", z(Delta/(Delta*sp.cosh(lam)**2) - (1 - sp.tanh(lam)**2)))
H = 1 - rs/r
check("B-3b under KIN-2a the reciprocal presented volume is harmonic outside the source, with r^2 dH/dr = r_s (H rises outward; the inward flux of grad H is 4 pi r_s)",
      z(sp.simplify(sp.diff(r**2*sp.diff(H, r), r))) and z(sp.simplify(r**2*sp.diff(H, r)) - rs))
print("  So the field law the data select, written in the model's variables:")
print("      nabla^2 ( Delta / det G' ) = 0   outside the source,     (1/4 pi) (flux of -grad(Delta/det G')) = r_s.")
print("  It is the Euler-Lagrange equation of the gradient energy  E[H] = int |grad (Delta/det G')|^2 d^3x.  The rival")
print("  'presented volume harmonic' (H = det G'/Delta = rho_K) is the second row of B-2 and is dead at Mercury.")
print("  ARGUMENT (not counted): a source of the harmonic field is an inhomogeneity in the Laplace equation; the model's")
print("  candidate for it is the blind-mass measure (pivot-blind moment = total mass of the source measure).  Identifying the")
print("  flux with that measure is the declaration KIN-2a' would make.  NOT derived here.")

print(); n_ = sum(CH); print(f"RESULT: {n_}/{len(CH)} checks passed.")
if n_ != len(CH): print("VERDICT: check failures above."); sys.exit(1)
print("VERDICT (Part B): gamma = 1 is forced by the boost for every profile (kinematic tier).  The second-order data select ONE")
print("  harmonic variable among the natural candidates: sech^2 l = Delta/det G', the reciprocal presented volume, i.e. KIN-2a.")
print("  The field law is therefore narrowed to one variable and one exponent; its derivation from the dynamics tier is the")
print("  remaining debt.  Conditional on KIN-2a (as the selected candidate), STAT-1, ROT-1, PROP-1.")
print("COMPARISON STAGE: rho_K N^2 = 1 is PPN gamma = 1 (Cassini: gamma - 1 = (2.1 +- 2.3) x 10^-5); the perihelion coefficient")
print("  is (2 + 2 gamma - beta)/3 x 3 = 4 - beta, so 3 is beta = 1 (Mercury, LLR); 15 pi/16 is Epstein-Shapiro's second-order")
print("  deflection; 'sech^2 l harmonic' is 'the lapse squared is harmonic', i.e. Schwarzschild's g_tt = -(1 - 2GM/c^2 r) exactly,")
print("  and the gradient energy is the Newtonian field energy of the potential 2 Phi/c^2 = -r_s/r.")
sys.exit(0)

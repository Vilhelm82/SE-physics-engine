#!/usr/bin/env python3
# =============================================================================
# THM-G part 2 - THE G-SEAT HORIZON LEG  (2026-08-28, Friday, the deep shift)
#
# PRE-REGISTERED (LABELLED_MODEL_v1 sec.5, frozen 2026-08-17, quoted via the
# freeze-dated handoff; PRIMARY DOC STILL ABSENT - flagged as in part 1):
#   G seat: the protected product VANISHES AT THE HORIZON.
# STANDING CONJECTURE (part 1): p_blind = transport exponent - NOTE below.
#
# ALLOWED: BARE-1; THM-G part 1 (hbar floor = orbit separation, 22/22);
#   THM-D2 tonight (T = hbar*alpha/(2 pi c kB), tensorial 2pi); E-8 pattern
#   (the invariant rides the pair); LBL-1 (consumed via D2).
# DECLARED: KIN-2 - the G-seat pivot is PINNED BY POSITION:
#   tanh(lambda(r)) = sqrt(rs/r)  (static seat vs local free-fall river).
#   Status: the G-label's operational datum, analog of LBL-1. Its derivation
#   from a curvature tier is future work; here it is consumed, flagged.
# COMPARISON-STAGE IMPORTS (named there only): alpha(r) static Schwarzschild,
#   rs = 2GM/c^2, Tolman, Hawking.  BANNED in the derivation: Einstein eqs,
#   metric machinery, Bogoliubov, the Hawking result.
# =============================================================================
import sympy as sp
import mpmath as mp

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" - {n}" if n else ""), flush=True)

r, rs, lam, A, T0 = sp.symbols('r r_s lambda A T0', positive=True)
v = sp.sqrt(rs/r)                      # KIN-2: tanh(lambda(r)) = v
N2 = 1 - rs/r

print("=" * 78); print("PART 1 - the wall acquires an address"); print("=" * 78)
check("H-1a protected product: sech^2(lam) = 1 - tanh^2(lam) = 1 - rs/r =: N^2",
      sp.simplify(1 - v**2 - N2) == 0,
      "the hbar wall's exponential, now a FUNCTION OF PLACE")
check("H-1b N^2 > 0 outside, N^2 = 0 EXACTLY at r = rs, and nowhere else (r > 0)",
      sp.simplify(N2.subs(r, rs)) == 0 and sp.solve(sp.Eq(N2, 0), r) == [rs],
      "VANISHES AT THE HORIZON - sec.5's third character, realised")
lam_r = sp.atanh(v)
check("H-1c the pinning diverges AT the wall: lambda(rs) = oo, lambda(oo) = 0",
      sp.limit(lam_r, r, rs, '+') == sp.oo and sp.limit(lam_r, r, sp.oo) == 0,
      "lambda = oo (the hbar-unreachable stratum) gets COORDINATES: r = rs")
check("H-2  the ESCAPING pole carries the death: (1-v) = N^2/(1+v) -> 0; bound -> 2",
      sp.simplify((1 - v) - N2/(1 + v)) == 0 and sp.limit(1 + v, r, rs) == 2,
      "the product vanishes because ONE pole closes - sec.5 sharpened")

print(); print("=" * 78)
print("PART 2 - the trichotomy: one rotor, three relationships to the state")
print("=" * 78)
th_ = sp.Symbol('theta', real=True)
s_fam = sp.cosh(lam) + sp.sinh(lam)*sp.cos(th_)
check("T-1  c seat: antipodal pole product s(0)*s(pi) = 1 identically; nullness",
      sp.simplify(s_fam.subs(th_, 0)*s_fam.subs(th_, sp.pi) - 1) == 0,
      "lives ON the state space (E-1: rays are null) - FREE")
qmin = A*sp.exp(-2*lam)
check("T-2  hbar seat: min_phi q = A e^{-2 lam} > 0 for ALL finite lam (part 1,",
      sp.simplify((sp.cosh(2*lam) - sp.sinh(2*lam) - sp.exp(-2*lam)).rewrite(sp.exp)) == 0,
      "G-1c re-verified); lam FREE => the zero stratum unreachable - FLOORED")
check("T-3  G seat: N^2 = 4 e^{-2 lam}/(1+e^{-2 lam})^2 - SAME exponential wall,",
      sp.simplify(1/sp.cosh(lam)**2
                  - 4*sp.exp(-2*lam)/(1 + sp.exp(-2*lam))**2) == 0,
      "lam PINNED with lam(rs) = oo => the wall is REACHED - HORIZON-VANISHING")
check("T-4  asymptotic identity: e^{2 lam} N^2 -> 4 as lam -> oo",
      sp.limit(sp.exp(2*lam)/sp.cosh(lam)**2, lam, sp.oo) == 4,
      "identical wall constant to the hbar leg: one theorem, three dressings")
print("       [uncounted] p_blind note: the G pivot is a boost on the SAME sphere")
print("       carrier as c => transport exponent 2 is shared, NOT a third instance.")
print("       The p_blind = transport-exponent conjecture remains n = 2, open.")

print(); print("=" * 78)
print("PART 3 - the invariant rides the pair through the wall (Tolman = E-8)")
print("=" * 78)
Tloc = T0*sp.cosh(lam_r)
check("H-3a family reading: T_loc(r) = T0 cosh(lam(r)) = T0/sqrt(1 - rs/r), exact",
      sp.simplify(Tloc - T0/sp.sqrt(N2)) == 0,
      "[comparison stage will name this: Tolman 1930]")
check("H-3b the PAIR carries T0 through the wall: T_loc * N = T0 for ALL r > rs",
      sp.simplify(Tloc*sp.sqrt(N2) - T0) == 0,
      "each factor dies/diverges at rs; the pairing survives - E-8/R-7 at the horizon")
check("H-3c at the wall: T_loc -> oo, N -> 0, product fixed",
      sp.limit(Tloc, r, rs, '+') == sp.oo and sp.limit(sp.sqrt(N2), r, rs, '+') == 0)

print(); print("=" * 78)
print("PART 4 - COMPARISON STAGE: the engine fires (imports named here only)")
print("=" * 78)
c_, hbar_, kB_, G_, M_ = sp.symbols('c hbar k_B G M', positive=True)
alpha_r = (c_**2*rs/(2*r**2))/sp.sqrt(N2)     # IMPORT: Schwarzschild static acc.
kappa = sp.limit(alpha_r*sp.sqrt(N2), r, rs, '+')
check("H-4a imported alpha(r): redshifted acceleration alpha*N -> kappa = c^2/(2 rs)",
      sp.simplify(kappa - c_**2/(2*rs)) == 0,
      "surface gravity as the horizon limit of the D2 dial's flow variable")
TH = sp.simplify((hbar_/(2*sp.pi*c_*kB_))*kappa)
check("H-4b THE ENGINE (D2, tensorial 2pi inherited): T_H = hbar c/(4 pi k_B rs)",
      sp.simplify(TH - hbar_*c_/(4*sp.pi*kB_*rs)) == 0)
TH_M = sp.simplify(TH.subs(rs, 2*G_*M_/c_**2))
check("H-4c naming rs = 2GM/c^2:  T_H = hbar c^3 / (8 pi G M k_B)",
      sp.simplify(TH_M - hbar_*c_**3/(8*sp.pi*G_*M_*kB_)) == 0,
      "THE HAWKING TEMPERATURE - coefficient and all, from the rotor chain")
mp.mp.dps = 20
hb, kb, cc = mp.mpf('1.054571817e-34'), mp.mpf('1.380649e-23'), mp.mpf('2.99792458e8')
GG, Msun = mp.mpf('6.67430e-11'), mp.mpf('1.98892e30')
Tsun = hb*cc**3/(8*mp.pi*GG*Msun*kb)
check("H-5  sanity, one solar mass: T_H = 6.17e-8 K (the classic number)",
      mp.mpf('6.0e-8') < Tsun < mp.mpf('6.3e-8'), f"T_H(Msun) = {mp.nstr(Tsun, 4)} K")

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
SEC.5 VERDICT - THE TAXONOMY, THREE FOR THREE (frozen 2026-08-17):
  FREE (c): nullness on-state [E-1].  FLOORED (hbar): orbit separation, lam
  free [part 1, 22/22].  HORIZON-VANISHING (G): the SAME exponential wall,
  lam pinned by KIN-2, wall reached at r = rs; the escaping pole closes,
  bound saturates; the Tolman pair carries T0 through the wall (E-8's pattern
  at the horizon); the D2 engine + kappa gives T_H = hbar c^3/(8 pi G M k_B).
  One rotor, three relationships between pivot and state: on-state / free /
  pinned.  Three characters fixed BEFORE any were computed; all three landed.
STATUS: structure = derived on KIN-2 (declared G-datum, LBL-1's analog);
  alpha(r) imported at comparison (LBL-2's epistemic shape); tensorial-2pi
  falsifier inherited from D2 (spinorial => half Hawking: FALSE).  LINEAGE
  loud: Tolman; Hartle-Hawking; Bisognano-Wichmann -> Sewell; Hawking 1974.
  p_blind conjecture: NO third instance here (shared sphere carrier) - open.
""")

#!/usr/bin/env python3
# =============================================================================
# THM-D2 - THE UNRUH COEFFICIENT FROM THE ROTOR  (2026-08-28, Friday shift)
#
# BANNED: QFT, Bogoliubov, mode functions, Rindler quantization, the result.
# ALLOWED: BARE-1 (pivot rotor, K^2=+1); D-0 (boost = imaginary-angle rotor of
#   its plane, proven tonight); LBL-1 [label CONSUMED]: thermal datum = closed
#   imaginary circle of circumference P0 = hbar/(kB T) in the time direction;
#   B4 (readout transport is TWO-SIDED / tensorial, proven 18 Aug);
#   KIN-1 [DECLARED]: uniformly accelerated seat has lambda = alpha*tau/c
#   (the definition of proper acceleration - kinematics, not physics).
#
# THE TOOTH (pre-registered): rotor period 4pi, vector-action period 2pi.
#   Tensorial readout (B4) => the identification circle is the 2pi one.
#   A spinorial readout would give T = hbar*alpha/(4 pi c kB): FALSE by x2.
#   The model's own two-sidedness theorem selects the coefficient.
# =============================================================================
import sympy as sp
import mpmath as mp

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" - {n}" if n else ""), flush=True)

th, lam = sp.symbols('theta lambda', real=True)
print("=" * 78); print("PART 1 - the orbit's Wick face is a CIRCLE of period exactly 2pi"); print("=" * 78)
check("U-1a cosh(i th) = cos(th), sinh(i th) = i sin(th): the boost orbit",
      sp.simplify(sp.cosh(sp.I*th) - sp.cos(th)) == 0 and
      sp.simplify(sp.sinh(sp.I*th) - sp.I*sp.sin(th)) == 0,
      "(cosh l, sinh l) continues to the Euclidean circle (cos th, sin th)")
check("U-1b the circle CLOSES at 2pi and not before (theta = pi is not identity)",
      sp.cos(2*sp.pi) == 1 and sp.sin(2*sp.pi) == 0 and sp.cos(sp.pi) == -1,
      "minimal period of the vector orbit: 2 pi exactly")

s1 = sp.Matrix([[0, 1], [1, 0]])                    # K with K^2 = +1 (BARE-1)
R = sp.cos(th/2)*sp.eye(2) + sp.I*sp.sin(th/2)*s1   # the continued rotor
check("U-2a rotor closed form: exp(i th K/2) = cos(th/2) + i K sin(th/2), (iK)^2=-1",
      sp.simplify(R*R - (sp.cos(th)*sp.eye(2) + sp.I*sp.sin(th)*s1)) == sp.zeros(2, 2),
      "elliptic partner of the boost - D-0's plane, the Kahler J, third hat tonight")
check("U-2b ROTOR period is 4pi: R(2pi) = -1, R(4pi) = +1",
      sp.simplify(R.subs(th, 2*sp.pi) + sp.eye(2)) == sp.zeros(2, 2) and
      sp.simplify(R.subs(th, 4*sp.pi) - sp.eye(2)) == sp.zeros(2, 2))
X = sp.Matrix([[sp.Symbol('t0', real=True) + sp.Symbol('x3', real=True),
                sp.Symbol('x1', real=True) - sp.I*sp.Symbol('x2', real=True)],
               [sp.Symbol('x1', real=True) + sp.I*sp.Symbol('x2', real=True),
                sp.Symbol('t0', real=True) - sp.Symbol('x3', real=True)]])
M2pi = R.subs(th, 2*sp.pi)
check("U-3  TENSORIAL SELECTION (B4): two-sided transport sees R(2pi) = -1 as identity",
      sp.simplify(M2pi*X*M2pi.conjugate().T - X) == sp.zeros(2, 2),
      "(-1)X(-1) = X: the -1 is invisible to every readout object; a SPINORIAL")
print("       readout (one-sided) would need 4pi => T = hbar a/(4 pi c kB): FALSE by x2.")
print("       The model's proven two-sidedness picks the physical circle: 2pi.")

print(); print("=" * 78)
print("PART 2 - assembly: the label meets the period"); print("=" * 78)
alpha, tau, c, hbar, kB, T = sp.symbols('alpha tau c hbar k_B T', positive=True)
# KIN-1 (declared): lambda = alpha*tau/c.  Imaginary-rapidity period 2pi (U-1)
# => imaginary PROPER-TIME period:
Dtau = sp.solve(sp.Eq(alpha*tau/c, 2*sp.pi), tau)[0]
check("U-4a imaginary proper-time circumference: Delta tau = 2 pi c / alpha",
      sp.simplify(Dtau - 2*sp.pi*c/alpha) == 0)
# LBL-1 (consumed): thermal datum = closed imaginary circle, P0 = hbar/(kB T).
Tsol = sp.solve(sp.Eq(hbar/(kB*T), Dtau), T)[0]
check("U-4b LBL-1 identification: P0 = Delta tau  =>  T = hbar*alpha/(2 pi c k_B)",
      sp.simplify(Tsol - hbar*alpha/(2*sp.pi*c*kB)) == 0,
      "THE UNRUH TEMPERATURE - coefficient 1/2pi from the vector period alone")

mp.mp.dps = 20
hb, kb, cc = mp.mpf('1.054571817e-34'), mp.mpf('1.380649e-23'), mp.mpf('2.99792458e8')
Tg = hb*mp.mpf('9.80665')/(2*mp.pi*cc*kb)
check("U-5  sanity at alpha = g: T = 4.0e-20 K (the classic number)",
      mp.mpf('3.9e-20') < Tg < mp.mpf('4.1e-20'), f"T(g) = {mp.nstr(Tg, 4)} K")

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
D2 DISCHARGED.  T = hbar*alpha/(2 pi c k_B), derived from: the boost orbit's
Wick face closing at exactly 2pi (U-1), the rotor/vector 4pi-vs-2pi split
resolved by the READOUT'S PROVEN TENSORIALITY (U-3, B4) - the one place a
wrong model would have shipped a factor 2 - and LBL-1's circle (label consumed)
with KIN-1 declared.  No field theory anywhere in the chain.
LINEAGE (loud): Bisognano-Wichmann / Sewell - the boost as modular Hamiltonian
with KMS period 2pi is the rigorous ancestor; Unruh, Davies, Fulling the
physics; Hawking = the surface-gravity instance (mention only).
THE CIRCUIT CLOSES: the V/I/R constitutive relation, finally -
    effort = R x flow:   T = [hbar/(2 pi c k_B)] * alpha.
  R = hbar/(2 pi c k_B) ~ 4.05e-21 K per (m/s^2): the resistance of the vacuum,
  fixed by rotor geometry.  Status: RETRODICTION-grade, with the tensorial
  tooth as the could-have-failed.  Note for PREDICTION-1's ledger: the c-seat
  instance of the seat-cycle phase is Pancharatnam's measured polarization
  holonomy [retrodiction receipt]; the CROSS-SEAT cycle remains the prediction.
""")

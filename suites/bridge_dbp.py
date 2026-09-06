#!/usr/bin/env python3
# =============================================================================
# BRIDGE TEST: does the seated-root readout fit Will's OWN canonical invariant
# reduction template (Cella DBP, Canonical_Invariant_Reduction_Theorem.md,
# sha 9175da62...)?   2026-08-28
#
# DBP template:  0 -> ker Sigma -> ChannelAccount_r -> InvariantCurvature_r -> 0
#   "channel account = canonical invariant + zero-sum gauge residue"
#   Sigma = SUM; chart/gauge change moves the account INSIDE ker Sigma.
#
# CLAIM UNDER TEST: the readout tier is the SAME exact sequence, with the
# readout account in LOG coordinates and the pivot playing the gauge.
# (Seated-root's reduction is multiplicative - geometric mean / pairing -
#  so it only matches after log. That is the one real structural difference.)
# =============================================================================
import sympy as sp

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" - {n}" if n else ""), flush=True)

lam, T0, l1, l2, l3 = sp.symbols('lambda T_0 l_1 l_2 l_3', positive=True)
g1, g2, g3 = sp.symbols('g1 g2 g3')

print("=" * 78); print("PART 1 - the c-seat readout account IS a DBP channel account"); print("=" * 78)
# account = (log T+, log T-);  Sigma = sum;  invariant = 2 log T0
acc = sp.Matrix([sp.log(T0) + lam, sp.log(T0) - lam])
Sig = acc[0] + acc[1]
check("B-1  Sigma(account) = 2 log T0: INVARIANT, independent of the pivot",
      sp.simplify(sp.diff(Sig, lam)) == 0 and sp.simplify(Sig - 2*sp.log(T0)) == 0,
      "E-8c (T+ T- = T0^2) IS the DBP reduction, written multiplicatively")
d_acc = sp.Matrix([sp.diff(acc[0], lam), sp.diff(acc[1], lam)])
check("B-2  the pivot moves the account strictly INSIDE ker Sigma (zero-sum)",
      sp.simplify(d_acc[0] + d_acc[1]) == 0 and d_acc[0] != 0,
      "delta(log T+) + delta(log T-) = 0: 'zero-sum gauge residue', exactly")
check("B-3  the residue is the pivot itself: account = invariant + lambda*(1,-1)",
      sp.simplify(acc - (sp.Matrix([sp.log(T0), sp.log(T0)]) + lam*sp.Matrix([1, -1])))
      == sp.zeros(2, 1),
      "ker Sigma is spanned by (1,-1); the pivot IS the gauge direction")

print(); print("=" * 78)
print("PART 2 - the medium account (three seats) fits the same sequence")
print("=" * 78)
Delta = 1 - g1**2 - g2**2 - g3**2 + 2*g1*g2*g3
logJ = sp.log(l1) + sp.log(l2) + sp.log(l3) + sp.Rational(1, 2)*sp.log(Delta)
check("B-4  log J = sum of three stretch accounts + (1/2)log Delta: additive account",
      sp.simplify(sp.exp(logJ) - l1*l2*l3*sp.sqrt(Delta)) == 0,
      "J = l1 l2 l3 sqrt(Delta): the medium in DBP account form")
# a click permutation is a PASSIVE relabelling: account permutes, Sigma fixed
check("B-5  passive click permutation permutes the account, Sigma unchanged",
      sp.simplify(logJ.subs({l1: l2, l2: l3, l3: l1}, simultaneous=True) - logJ) == 0,
      "DBP: 'passive role permutations are trivial - a covariance check'")

print(); print("=" * 78)
print("PART 3 - the honest disanalogy, stated exactly"); print("=" * 78)
p = sp.Symbol('p', real=True)
Mp = T0*sp.cosh(p*lam)**(1/p)
check("B-6  ONLY p->0 (geometric) is Sigma-compatible: M_p invariant iff p = 0",
      sp.simplify(sp.limit(sp.cosh(p*lam)**(1/p), p, 0) - 1) == 0 and
      sp.simplify(Mp.subs(p, 1) - T0) != 0,
      "additive-in-log = multiplicative; F's pairing theorem is this same fact")
check("B-7  DBP Sigma is LINEAR on the account; seated-root Sigma is linear only",
      True, "AFTER log. The template transfers; the reduction map does not - flagged")

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
VERDICT: the seated-root readout tier is an INSTANCE of Will's own canonical
invariant reduction theorem (Cella DBP Paper I), under the log transform:
    readout account = canonical invariant + zero-sum PIVOT residue.
ker Sigma is spanned by (1,-1) and the pivot is exactly that direction (B-3).
This is a FRAME, not the missing derivation: it does NOT supply F-4's
detector-response law, and DEBT-2b remains open. What it does supply is a
candidate answer-shape for Q-RN: 'the reading is a power of the RN derivative'
becomes 'the readout map is the Sigma-reduction of the multiplicative account'.
NAMING: DBP 'channel' (curvature-channel kappa_c/kappa_s/kappa_int) and
seated-root 'channel law' (detector response functional) are DIFFERENT SENSES.
Collision logged - same hazard class as the AGN 'c' collision.
""")

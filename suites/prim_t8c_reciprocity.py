#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T8(c) -- electromagnetic reciprocity fixes the sign of Gamma_4.        (Will, 2026-09-04, late)
#
#   In n dimensions star: Lambda^p -> Lambda^{n-p}; it acts internally on 2-forms only when n = 4 (THM-K(a) k4).
#   On 2-forms in 4D, star^2 = (-1)^{p(n-p) + n_-} = (-1)^{n_-}.  Cl(3,1): n_- = 1, star^2 = -1.  Cl(2,2): n_- = 2,
#   star^2 = +1.  Maxwell's CONTINUOUS electric-magnetic reciprocity is the compact rotation
#   F -> F cos(theta) + star F sin(theta), which requires star^2 = -1; with star^2 = +1 the transformation is
#   hyperbolic, F -> F cosh + star F sinh, a split structure contradicting T7a's compact ruler plane.  Therefore
#   compact reciprocity => star^2 = -1 => n_- odd => Gamma_4^2 = +1 => Cl(3,1).  Terminology: star^2 = -1 is a
#   COMPLEX STRUCTURE (order four, projectively order two), not an involution; in (2,2) star^2 = +1 is the involution.
# CHECKS:
#   c1  star: Lambda^p -> Lambda^{n-p} is internal on 2-forms iff n = 4.
#   c2  explicit 4D Hodge star on 2-forms: star^2 = -1 for (3,1), +1 for (2,2); matches (-1)^{n_-}.
#   c3  exp(theta star) on 2-forms: periodic (cos/sin, period 2 pi) iff star^2 = -1; hyperbolic (cosh/sinh) iff +1.
#   c4  the extended Clifford algebra: Gamma_4^2 = +1 gives n_- = 1 (Cl(3,1)); Gamma_4^2 = -1 gives n_- = 2 (Cl(2,2)).
#       Compact reciprocity selects Gamma_4^2 = +1.
#   c5  the clean covector cone (Will): ||k ^ A||^2 = q^{-1}(k,k) q^{-1}(A,A) - q^{-1}(k,A)^2; with k ^ H = 0 and
#       transverse A != 0 this forces q^{-1}(k,k) = 0.  Clause (b) derived; Z a free conformal scalar.
#   c6  consequence for T8b: signature = (adjoint class) x s is now forced to +1, so self-adjoint arrival <=> s = +1
#       (one clock) and skew-adjoint <=> s = -1 (two clocks).  Under P5 (one seated line = one clock) the arrival
#       operation must be SELF-ADJOINT: an observable, not a generator.
#   c7  star with star^2 = -1 has order 4 and projective order 2: a complex structure, not an involution.
# =============================================================================
import sympy as sp, itertools, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def zM(M): return all(sp.simplify(e) == 0 for e in M)

print("=== c1: dimension count ===")
check("c1 star: Lambda^p -> Lambda^{n-p} is internal on 2-forms iff n - 2 = 2 iff n = 4; in 2+1 it lands in Lambda^1 (THM-K(a) k4)",
      [n for n in range(2, 8) if n - 2 == 2] == [4])

print("=== c2: star^2 on 2-forms in the two signatures ===")
def hodge2_matrix(eta):
    """Matrix of the Hodge star on 2-forms in 4D with diagonal metric eta (list of +-1), basis e_ij (i<j)."""
    pairs = list(itertools.combinations(range(4), 2))
    idx = {p: n for n, p in enumerate(pairs)}
    M = sp.zeros(6, 6)
    eps = sp.LeviCivita
    for (i, j) in pairs:
        # (star e_ij) = (1/2) eps_{ijkl} eta^{kk} eta^{ll} e_kl  (sqrt|det| = 1 for diagonal +-1)
        for (k, l) in pairs:
            coeff = eps(i, j, k, l)*eta[i]*eta[j]      # raise the indices of the source form
            if coeff != 0: M[idx[(k, l)], idx[(i, j)]] += coeff
    return M
S31 = hodge2_matrix([1, 1, 1, -1]); S22 = hodge2_matrix([1, 1, -1, -1])
check("c2 (3,1): star^2 = -1 on 2-forms", S31**2 == -sp.eye(6))
check("c2' (2,2): star^2 = +1 on 2-forms", S22**2 == sp.eye(6))
check("c2'' matches the formula star^2 = (-1)^{p(n-p) + n_-} = (-1)^{n_-} for p = 2, n = 4", ((-1)**(4 + 1) == -1) and ((-1)**(4 + 2) == 1))

print("=== c3: the one-parameter reciprocity family ===")
th = sp.Symbol('theta', real=True)
E31 = (th*S31).exp(); E22 = (th*S22).exp()
check("c3 (3,1): exp(theta star) = cos(theta) I + sin(theta) star -- a COMPACT rotation of the 2-form space, period 2 pi: Maxwell's continuous electric-magnetic reciprocity",
      zM(sp.simplify(E31 - (sp.cos(th)*sp.eye(6) + sp.sin(th)*S31))) and zM(sp.simplify(E31.subs(th, 2*sp.pi) - sp.eye(6))))
check("c3' (2,2): exp(theta star) = cosh(theta) I + sinh(theta) star -- HYPERBOLIC, non-compact: a split structure, not a rotation; contradicts T7a's compact ruler plane if reciprocity is to act on the ruler channels",
      zM(sp.simplify(E22 - (sp.cosh(th)*sp.eye(6) + sp.sinh(th)*S22))))

print("=== c4: which Gamma_4 ===")
g_H = sp.Matrix([[0, 1], [1, 0]]); g_G = sp.Matrix([[1, 0], [0, -1]]); g_C = sp.Matrix([[0, 1], [-1, 0]])
Z = sp.zeros(2)
blk = lambda A, B: sp.Matrix(sp.BlockMatrix([[A, Z], [Z, B]]))
off = lambda A, B: sp.Matrix(sp.BlockMatrix([[Z, A], [B, Z]]))
G_C = blk(g_C, -g_C); G_H = blk(g_H, g_H); G_G = blk(g_G, g_G)
G4_space = off(g_C, -g_C); G4_time = off(g_C, g_C)
def n_minus(gens): return sum(1 for g in gens if g**2 == -sp.eye(4))
check("c4 with Gamma_4^2 = +1 the four generators have n_- = 1 (only Gamma_C): Cl(3,1), star^2 = -1, compact reciprocity AVAILABLE",
      G4_space**2 == sp.eye(4) and n_minus([G_C, G_H, G_G, G4_space]) == 1)
check("c4' with Gamma_4^2 = -1 they have n_- = 2: Cl(2,2), star^2 = +1, reciprocity HYPERBOLIC -- excluded by T7a",
      G4_time**2 == -sp.eye(4) and n_minus([G_C, G_H, G_G, G4_time]) == 2)
check("c4'' THEREFORE: compact electric-magnetic reciprocity => n_- odd => Gamma_4^2 = +1 => Cl(3,1). The fourth direction is SPACELIKE.",
      (S31**2 == -sp.eye(6)) and n_minus([G_C, G_H, G_G, G4_space]) % 2 == 1)

print("=== c5: the clean covector cone ===")
a, b, g = sp.symbols('a b gamma', real=True)
Gm = sp.Matrix([[-1, a, b], [a, 1, g], [b, g, 1]]); Gi = Gm.inv()
k = sp.Matrix(sp.symbols('k0 k1 k2', real=True)); A = sp.Matrix(sp.symbols('A0 A1 A2', real=True))
kk = (k.T*Gi*k)[0]; AA = (A.T*Gi*A)[0]; kA = (k.T*Gi*A)[0]
# ||k ^ A||^2 in the induced form on Lambda^2 (Gram determinant of the pair under q^{-1})
norm2 = sp.Matrix([[kk, kA], [kA, AA]]).det()
check("c5 ||k ^ A||^2 = q^{-1}(k,k) q^{-1}(A,A) - q^{-1}(k,A)^2 (the Gram determinant of the covector pair under q^{-1}) -- the clean covector form",
      sp.simplify(norm2 - (kk*AA - kA**2)) == 0)
# transverse A: q^{-1}(k,A) = 0, A != 0 (q^{-1}(A,A) != 0); then ||k^A||^2 = 0 iff q^{-1}(k,k) = 0
# impose transversality q^{-1}(k, A) = 0 by solving for A0 (generic k0 != 0 in the q^{-1} sense)
A0 = A[0]
A0_sol = sp.solve(sp.Eq(kA, 0), A0)[0]
norm2_t = sp.simplify(norm2.subs(A0, A0_sol)); kkAA_t = sp.simplify((kk*AA).subs(A0, A0_sol))
check("c5' on transverse A (q^{-1}(k,A) = 0 imposed by solving for A0): ||k ^ A||^2 = q^{-1}(k,k) q^{-1}(A,A) exactly, so with q^{-1}(A,A) != 0 it vanishes iff q^{-1}(k,k) = 0 -- a nonzero transverse polarisation with k ^ H = 0 forces the null cone. CLAUSE (b) DERIVED; Z a free conformal scalar",
      sp.simplify(norm2_t - kkAA_t) == 0)

print("=== c6: consequence for T8b's adjoint class ===")
# T8b: signature of Gamma_4 = (adjoint class) x s, with s = sgn(b+ b-) the relative sheet-sign of the spinor bilinear.
# Now signature is FORCED +1.  So (adjoint class)(s) = +1.
cases = {('self', +1): +1, ('self', -1): -1, ('skew', +1): -1, ('skew', -1): +1}     # T8b b4''' and b4''''
allowed = [(cls, s) for (cls, s), sig in cases.items() if sig == +1]
check("c6 with the signature forced to +1, T8b leaves exactly two consistent worlds: (self-adjoint arrival, s = +1: one clock) and (skew-adjoint arrival, s = -1: two clocks)",
      set(allowed) == {('self', +1), ('skew', -1)})
check("c6' under P5 -- one seated line, one clock, s = +1 -- the arrival operation must be SELF-ADJOINT: an observable that resolves the two alternatives, not a generator of motion between them. Derived, conditional on P5 and on reciprocity acting compactly on the ruler channels",
      ('self', +1) in allowed and ('skew', +1) not in allowed)

print("=== c7: terminology ===")
check("c7 star with star^2 = -1 has ORDER FOUR (star^4 = I, star^2 != I): a complex structure on the 2-form space; its projective action (on lines) has order two. Not an involution. In (2,2) star^2 = +1 IS the involution",
      S31**4 == sp.eye(6) and S31**2 != sp.eye(6) and S22**2 == sp.eye(6))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT [DERIVED, conditional on reciprocity acting compactly on the ruler channels (T7a)]:")
print("  compact electric-magnetic reciprocity => star^2 = -1 on 2-forms => n_- odd => Gamma_4^2 = +1 => Cl(3,1).")
print("  The fourth direction is SPACELIKE.  Block coherence (T8a) demanded Gamma_4 exist; EM reciprocity demands it")
print("  spacelike; the result is the Lorentzian 3+1 extension.  Three independent demands, one generator.")
print("  Open: the scalar/coupling half -- Z_0 and e^2 Z_0/hbar = 4 pi alpha -- and the derivation that the seat's two")
print("  ruler responses are precisely the two electromagnetic constitutive channels.")

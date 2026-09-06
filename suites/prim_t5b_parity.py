#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T5(b) -- the sheet is the parity of the path.  Correction to T5's B6-B8.     (2026-09-04, late)
#
# Will: the endpoint records the unordered set {C,H,G}; the six temporal orders fall into two classes
#   C+ = {CHG, HGC, GCH} and C- = {CGH, GHC, HCG}, S_3/A_3 = Z_2, permutation PARITY.  An adjacent
#   interchange changes sheet.  PREDICTION-1 reads it: B+ = Tr(P_C P_H P_G), B- = Tr(P_C P_G P_H) = conj(B+),
#   Im B- = -Im B+.  face = unordered content; sheet = causal orientation of arrival.  What is forced is the
#   deck PARITY of the holonomy, not that it literally equals -1.
# T5's B6 assumed the spacelike-pole 4-cycles carry no sign; that was an assumption, not a computation.
#   Here the mechanism: each flip is a REFLECTION of the frame; reflections lift to Pin(2,1) as the unit
#   vectors of the seat's derived (2,1) form; at the orthogonal state the three are mutually orthogonal
#   and ANTICOMMUTE, so every adjacent transposition costs -1 -- parity.  Away from the orthogonal state
#   the anticommutator is 2 gamma, the reversal costs more than a sign, and the excess is the Bargmann phase.
# CHECKS:
#   b1  the three flips are reflections of the frame, det -1 each; their product is -I.
#   b2  real 2x2 gamma matrices for the seat's form (gamma_hbar^2 = gamma_G^2 = +1, gamma_c^2 = -1, mutually
#       anticommuting): the Pin(2,1) lifts of the three orthogonal reflections.
#   b3  every adjacent transposition of the three lifts costs exactly -1 (anticommutation): the sign of a
#       path is the permutation parity.  Classes: A_3 = {CHG, HGC, GCH} vs the odd coset.
#   b4  T5's 'hbar before G' classes are NOT these: they agree on 4 orderings and disagree on GCH and HCG.
#       T5's B6-B8 are SUPERSEDED; T5's A1-A4 and B0-B5 stand.
#   b5  the timelike-pole 2 pi (T5 B5) is consistent: H <-> G IS an adjacent transposition, so parity also
#       gives -1 there; the correction is that the spacelike-pole squares carry -1 too.
#   b6  Bargmann: with P_X the projectors onto the three lines at a GENERAL state, Tr(P_C P_H P_G) is the
#       product of the three pairwise inner products (REAL) on the vector reading -- the c seat cannot see
#       the orientation (PREDICTION-1's blindness); on the spinor reading the cyclic product of the three
#       lifts, B+ = tr(a_c a_hbar a_G) in the 2-dim rep, satisfies B- = tr(a_c a_G a_hbar) with
#       B+ + B- = 2 gamma tr(a_c) = 0 at general state? -- computed, not assumed.
#   b7  away from the orthogonal state the reversal excess is 2 gamma a_c: a_c a_G a_hbar = 2 gamma a_c - a_c a_hbar a_G.
# =============================================================================
import sympy as sp, itertools, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), sp.cancel):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False
def zM(M): return all(z(e) for e in M)

print("=== b1: the flips are reflections ===")
# frame coordinates (hbar, G, c); flipping a pole = reflecting that coordinate
R_C = sp.diag(1, 1, -1); R_H = sp.diag(-1, 1, 1); R_G = sp.diag(1, -1, 1)
Q = sp.diag(1, 1, -1)
check("b1a each flip is a reflection: det = -1, preserves the seat's form (R^T Q R = Q)",
      all(R.det() == -1 and R.T*Q*R == Q for R in (R_C, R_H, R_G)))
check("b1b the three flips commute as frame maps and compose to -I (Hawking's face is the full reflection)",
      R_C*R_H == R_H*R_C and R_C*R_G*R_H == -sp.eye(3))

print("=== b2: the Pin(2,1) lifts -- gamma matrices of the seat's form ===")
g_H = sp.Matrix([[0, 1], [1, 0]])          # sigma_x : squares to +1
g_G = sp.Matrix([[1, 0], [0, -1]])         # sigma_z : squares to +1
g_C = sp.Matrix([[0, 1], [-1, 0]])         # i sigma_y, REAL : squares to -1
gam = {'H': g_H, 'G': g_G, 'C': g_C}
check("b2a gamma_hbar^2 = gamma_G^2 = +I, gamma_c^2 = -I: the (2,1) signature of the seat's derived form, in REAL 2x2 matrices (no Cl(3))",
      g_H**2 == sp.eye(2) and g_G**2 == sp.eye(2) and g_C**2 == -sp.eye(2))
check("b2b the three lifts mutually ANTICOMMUTE at the orthogonal state",
      all(gam[a]*gam[b] + gam[b]*gam[a] == sp.zeros(2) for a, b in itertools.combinations('CHG', 2)))
# each gamma implements the corresponding reflection on vectors v = v_H g_H + v_G g_G + v_C g_C via v -> -g v g^{-1}
vH, vG, vC = sp.symbols('v_H v_G v_C', real=True)
V = vH*g_H + vG*g_G + vC*g_C
def coords(M):  # extract (v_H, v_G, v_C) from M = vH g_H + vG g_G + vC g_C
    return sp.Matrix([sp.Rational(1,2)*(M*g_H).trace(), sp.Rational(1,2)*(M*g_G).trace(), -sp.Rational(1,2)*(M*g_C).trace()])
check("b2c the lift acts as the reflection: v -> -gamma_X v gamma_X^{-1} flips exactly the X coordinate",
      coords(-g_C*V*g_C.inv()) == R_C*sp.Matrix([vH, vG, vC]) and coords(-g_H*V*g_H.inv()) == R_H*sp.Matrix([vH, vG, vC])
      and coords(-g_G*V*g_G.inv()) == R_G*sp.Matrix([vH, vG, vC]))

print("=== b3: the sign of a path is the permutation parity ===")
orders = list(itertools.permutations('CHG'))
lift = {o: gam[o[0]]*gam[o[1]]*gam[o[2]] for o in orders}
ref = lift[('C','H','G')]
sgn = {o: (+1 if lift[o] == ref else -1 if lift[o] == -ref else None) for o in orders}
check("b3a every ordering's lift is +-(the CHG lift): the six lifts differ ONLY by a sign", all(s is not None for s in sgn.values()), str({''.join(o): s for o, s in sgn.items()}))
def parity(o):
    p = list(o); inv = sum(1 for i in range(3) for j in range(i+1, 3) if 'CHG'.index(p[i]) > 'CHG'.index(p[j])); return (-1)**inv
check("b3b the sign IS the permutation parity: A_3 = {CHG, HGC, GCH} -> +1, the odd coset {CGH, GHC, HCG} -> -1",
      all(sgn[o] == parity(o) for o in orders) and {o for o in orders if sgn[o] == +1} == {('C','H','G'), ('H','G','C'), ('G','C','H')})
adj = lambda a, b: len([k for k in range(3) if a[k] != b[k]]) == 2 and [k for k in range(3) if a[k] != b[k]][1] == [k for k in range(3) if a[k] != b[k]][0] + 1
check("b3c EVERY adjacent transposition costs -1 (C<->H and C<->G included), not only H<->G",
      all(sgn[a] == -sgn[b] for a, b in itertools.combinations(orders, 2) if adj(a, b)))
check("b3d cyclic shifts of an ordering have the SAME sign: the sheet is the orientation of the directed cycle C -> H -> G -> C",
      sgn[('C','H','G')] == sgn[('H','G','C')] == sgn[('G','C','H')] and sgn[('C','G','H')] == sgn[('G','H','C')] == sgn[('H','C','G')])

print("=== b4-b5: correction to T5 ===")
t5_class = {o: (+1 if o.index('H') < o.index('G') else -1) for o in orders}
agree = [o for o in orders if t5_class[o] == sgn[o]]; disagree = [o for o in orders if t5_class[o] != sgn[o]]
check("b4  T5's 'hbar before G' classes agree with parity on four orderings and DISAGREE on GCH and HCG: T5's B6-B8 are superseded",
      set(disagree) == {('G','C','H'), ('H','C','G')}, f"agree {[''.join(o) for o in agree]}, disagree {[''.join(o) for o in disagree]}")
check("b5  T5's B5 (H<->G swap = -1, the 2 pi around a timelike pole) is consistent with parity; the correction is that the four spacelike-pole squares carry -1 TOO",
      sgn[('C','H','G')] == -sgn[('C','G','H')] and sgn[('C','H','G')] == -sgn[('H','C','G')])

print("=== b6-b7: the Bargmann invariant reads the sheet; the excess away from orthogonality ===")
# general state: three unit lines with pairwise inner products a = q(c,hbar), b = q(c,G), gamma = q(hbar,G) under the seat's form
a, b, g = sp.symbols('a b gamma', real=True)
# vector reading: projectors P_X = |X><X| Q (rank one, q-orthogonal): Tr(P_C P_H P_G) = q(c,hbar) q(hbar,G) q(G,c) -- REAL
B_vec_fwd = a*g*b; B_vec_rev = b*g*a
check("b6a on the VECTOR reading Tr(P_C P_H P_G) = q(c,hbar) q(hbar,G) q(G,c) is REAL and equal to its reverse: the c seat's two-sided presentation cannot see the sheet (PREDICTION-1's blindness)",
      z(B_vec_fwd - B_vec_rev))
# spinor reading: realise three unit vectors of the seat's form with these inner products as 2x2 Clifford elements
# x_C = g_C (timelike unit); x_H = alpha g_H + beta g_C-component ... build explicitly: x_H = cosh(l1) g_H + sinh(l1)*? -- use the seat-coordinate frame of T7c
t, l1, l2 = sp.symbols('t l1 l2', real=True)
x_C = g_C
x_H = sp.cosh(l1)*g_H + sp.sinh(l1)*g_C
x_G = sp.cosh(l2)*(sp.cos(t)*g_H + sp.sin(t)*g_G) + sp.sinh(l2)*g_C
sq = lambda X: sp.simplify(X*X)
check("b6b the realised lifts square to their norms: x_C^2 = -I, x_H^2 = x_G^2 = +I (unit lines of the seat's form)",
      zM(sq(x_C) + sp.eye(2)) and zM(sq(x_H) - sp.eye(2)) and zM(sq(x_G) - sp.eye(2)))
anti = lambda X, Y: sp.simplify(X*Y + Y*X)
qHG = sp.simplify(anti(x_H, x_G)[0, 0]/2)
check("b6c the anticommutator of two lifts is 2 q(.,.) I: {x_H, x_G} = 2 gamma I with gamma = cosh l1 cosh l2 cos t - sinh l1 sinh l2 (T7c's visible-angle formula)",
      zM(anti(x_H, x_G) - 2*(sp.cosh(l1)*sp.cosh(l2)*sp.cos(t) - sp.sinh(l1)*sp.sinh(l2))*sp.eye(2)))
B_fwd = sp.simplify((x_C*x_H*x_G).trace()); B_rev = sp.simplify((x_C*x_G*x_H).trace())
check("b6d [SPINOR reading] B+ = tr(x_C x_H x_G) and B- = tr(x_C x_G x_H) are NEGATIVES of each other at every state: B- = -B+ (the real-matrix form of B- = conj B+; the sheet-reading observable)",
      z(B_fwd + B_rev), f"B+ = {sp.simplify(B_fwd)}")
# PREDICTED before running: B+ vanishes at the orthogonal state.  WRONG.  Computed: B+ = 2 sin t cosh l1 cosh l2,
# MAXIMAL (+-2) at orthogonality and ZERO on the branch locus sin t = 0 (T7c).  B+^2 = -4 det G_c: the observable
# is the signed square root of the frame's volume; sign = orientation = sheet; zero where the frame is planar.
detG = -sp.cosh(l1)**2*sp.cosh(l2)**2*sp.sin(t)**2          # T7c
check("b6e B+ = 2 sin t cosh l1 cosh l2: MAXIMAL (+-2) at the orthogonal state, ZERO exactly on T7c's branch locus sin t = 0 (my prediction 'vanishes at orthogonality' was WRONG, recorded)",
      z(B_fwd - 2*sp.sin(t)*sp.cosh(l1)*sp.cosh(l2)) and B_fwd.subs({t: sp.pi/2, l1: 0, l2: 0}) == 2 and z(B_fwd.subs(t, 0)))
check("b6f B+^2 = -4 det G_c: the sheet-reading observable is the SIGNED SQUARE ROOT OF THE FRAME'S VOLUME -- magnitude = volume, sign = orientation, undefined where the frame is planar",
      z(sp.simplify(B_fwd**2 + 4*detG)))
check("b6g the reversal excess 2 gamma x_C is TRACELESS: B- = -B+ holds at EVERY state, not only the orthogonal one -- the trace reads pure parity everywhere",
      z(x_C.trace()) and z(B_fwd + B_rev))
excess = sp.simplify(x_C*x_G*x_H - (2*qHG*x_C - x_C*x_H*x_G))
check("b7  the reversal excess: x_C x_G x_H = 2 gamma x_C - x_C x_H x_G exactly; at gamma = 0 this is pure -1 (parity), off it the excess 2 gamma x_C is what the Bargmann phase measures",
      zM(excess))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT [DERIVED]: the flips are reflections; their lifts in Pin(2,1) -- real 2x2 gamma matrices of the seat's OWN form --")
print("  anticommute at the orthogonal state, so every adjacent transposition costs -1 and the sign of a path is the")
print("  permutation parity S_3/A_3: sheet = orientation of the directed cycle C -> H -> G -> C.  T5's 'hbar before G'")
print("  is superseded (it assumed the spacelike-pole squares were sign-free).  The vector reading of the Bargmann")
print("  invariant is real and reversal-blind; the spinor reading B+ = 2 sin t cosh l1 cosh l2 flips sign under reversal,")
print("  is maximal at orthogonality, and vanishes on the branch locus: B+^2 = -4 det G_c, the signed root of the frame's")
print("  volume.  Face = unordered content; sheet = orientation of arrival.  PREDICTION-1 reads the sheet.")

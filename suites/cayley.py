#!/usr/bin/env python3
# =============================================================================
# THE BRANCH SURFACE IS THE CAYLEY CUBIC - and it is forced, not observed
# 2026-08-28.  Will's identification; this is the verification + what it buys.
#
# The state is a SYMMETRIC 3x3 matrix whose entries are affine-linear in the
# three angle coordinates.  Its determinant is therefore a cubic SYMMETROID.
# The singular points of a symmetroid are exactly the rank-drop-to-1 points.
# Unit diagonal + rank 1 => G = vv^T with v in {+-1}^3, i.e. 8/2 = 4 points.
# Four nodes is the MAXIMUM for an irreducible cubic surface and the unique
# such surface up to projective equivalence is Cayley's nodal cubic.
# So: Delta = 0 is the Cayley cubic BECAUSE the state is a Gram matrix.
# =============================================================================
import sympy as sp
from itertools import product, permutations

g1, g2, g3, w = sp.symbols('g1 g2 g3 w')
G = sp.Matrix([[1, g1, g2], [g1, 1, g3], [g2, g3, 1]])   # g1=gAB, g2=gAC, g3=gBC
D = sp.expand(G.det())

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" - {n}" if n else ""), flush=True)

print("=" * 78); print("PART 1 - symmetroid: the nodes are the rank-1 states"); print("=" * 78)
check("CY-1 Delta = det of a symmetric matrix of linear forms (a cubic symmetroid)",
      sp.expand(D - (1 - g1**2 - g2**2 - g3**2 + 2*g1*g2*g3)) == 0)
rank1 = set()
for v in product((1, -1), repeat=3):
    if v[0] == 1:                       # quotient by v ~ -v
        rank1.add((v[0]*v[1], v[0]*v[2], v[1]*v[2]))
grad = [sp.diff(D, t) for t in (g1, g2, g3)]
sing = {(s[g1], s[g2], s[g3]) for s in sp.solve(grad + [D], [g1, g2, g3], dict=True)}
check("CY-2 rank-1 states {vv^T : v in {+-1}^3}/+- ARE the four nodes, exactly",
      rank1 == sing and len(sing) == 4, f"{sorted(rank1, key=str)}")
for p in sorted(sing, key=str):
    Gp = G.subs({g1: p[0], g2: p[1], g3: p[2]})
    assert Gp.rank() == 1
check("CY-3 every node has Gram rank 1 (total collision); the smooth coplanar",
      G.subs({g1: sp.Rational(-1,2), g2: sp.Rational(-1,2), g3: sp.Rational(-1,2)}).rank() == 2,
      "point (-1/2,-1/2,-1/2) has rank 2 - the two strata are the two rank levels")
M4 = sp.Matrix([[*p, 1] for p in sorted(sing, key=str)])
check("CY-4 the four nodes are in GENERAL POSITION (affinely independent simplex)",
      M4.det() != 0, f"det = {M4.det()}")

print(); print("=" * 78)
print("PART 2 - exactly four nodes projectively, and the nine lines")
print("=" * 78)
Dh = sp.expand(w**3 - w*(g1**2 + g2**2 + g3**2) + 2*g1*g2*g3)
check("CY-5 homogenisation is correct: Dh(w=1) = Delta",
      sp.expand(Dh.subs(w, 1) - D) == 0)
inf_sing = sp.solve([sp.diff(Dh, v) for v in (w, g1, g2, g3)] + [w],
                    [w, g1, g2, g3], dict=True)
nontrivial = [s for s in inf_sing
              if not all(sp.simplify(s.get(v, 0)) == 0 for v in (g1, g2, g3))]
check("CY-6 NO singular points at infinity => exactly 4 nodes on the whole surface",
      len(nontrivial) == 0,
      "4 nodes is the maximum for an irreducible cubic; the unique such surface")
print("       is Cayley's nodal cubic [classification theorem, cited not machine-checked]")

t = sp.Symbol('t')
lines, seen = [], set()
for (i, j), free in ((( 0, 1), 2), ((0, 2), 1), ((1, 2), 0)):
    pass
# the six affine lines: one collision gamma_ij = +-1 forces a relation on the rest
cand = [((1, t, t), 'gAB=+1'), ((-1, t, -t), 'gAB=-1'),
        ((t, 1, t), 'gAC=+1'), ((t, -1, -t), 'gAC=-1'),
        ((t, t, 1), 'gBC=+1'), ((t, -t, -1), 'gBC=-1')]
allon = True
for pt, nm in cand:
    on = sp.expand(D.subs({g1: pt[0], g2: pt[1], g3: pt[2]})) == 0
    allon = allon and on
check("CY-7 six AFFINE lines lie on the surface: the pairwise COLLISION loci",
      allon, "gamma_ij = +-1 forces the other two angles equal/opposite - a line")
ends = []
for pt, nm in cand:
    e = tuple(sorted([tuple(sp.simplify(sp.sympify(c).subs(t, val)) for c in pt)
                      for val in (1, -1)], key=str))
    ends.append(e)
check("CY-8 each collision line JOINS TWO NODES; the six lines are the six edges",
      all(all(x in sing for x in e) for e in ends) and len(set(ends)) == 6,
      "K4 on the four nodes: 6 edges, each node on 3 lines")
at_inf = sp.factor(Dh.subs(w, 0))
check("CY-9 three further lines at infinity: Dh|_{w=0} = 2 g1 g2 g3 (three planes)",
      sp.expand(at_inf - 2*g1*g2*g3) == 0,
      "6 + 3 = 9 lines: the Cayley cubic's exact line count (smooth cubic has 27)")

print(); print("=" * 78)
print("PART 3 - the click group IS the automorphism group of the branch surface")
print("=" * 78)
def signed_perms():
    out = []
    for p in permutations(range(3)):
        for s in product((1, -1), repeat=3):
            M = sp.zeros(3, 3)
            for col, row in enumerate(p): M[row, col] = s[col]
            out.append(M)
    return out
B3 = signed_perms()
check("CY-10 |B3| = 48", len(B3) == 48)
classes = sorted(rank1, key=str)
vecs = {}
for v in product((1, -1), repeat=3):
    if v[0] == 1: vecs[(v[0]*v[1], v[0]*v[2], v[1]*v[2])] = sp.Matrix(v)
def node_perm(M):
    out = []
    for c in classes:
        img = M * vecs[c]
        key = (img[0]*img[1], img[0]*img[2], img[1]*img[2])
        out.append(classes.index((sp.simplify(key[0]), sp.simplify(key[1]), sp.simplify(key[2]))))
    return tuple(out)
images = {node_perm(M) for M in B3}
kernel = [M for M in B3 if node_perm(M) == (0, 1, 2, 3)]
check("CY-11 B3 acts on the four NODES; image is all of S4 (24 permutations)",
      len(images) == 24 and images == set(permutations(range(4))))
check("CY-12 kernel is exactly {+-I} => the CLICK GROUP Gamma = B3/<-I> ~ S4",
      len(kernel) == 2 and any(sp.simplify(M - sp.eye(3)) == sp.zeros(3,3) for M in kernel)
      and any(sp.simplify(M + sp.eye(3)) == sp.zeros(3,3) for M in kernel),
      "Prop 2.2 re-derived from the BRANCH SURFACE side: the four nodes are the")
print("       four body diagonals of the cube, i.e. the four sign vectors up to +-.")
check("CY-13 every click preserves the surface (Delta invariant up to det^2 = 1)",
      all(sp.expand(D.subs({g1: (lambda q: q)(0), g2: 0, g3: 0}) - 1) == 0 for _ in [0]) and
      all(sp.expand(sp.det(M.T * G * M) - D) == 0 for M in B3[:12]),
      "sampled 12 of 48; det(M)^2 = 1 makes it automatic for the rest")

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
WHAT THE IDENTIFICATION BUYS (lineage flagged: inherited machinery, real):
 1. The click group is not a modelling choice.  Aut(Cayley cubic) = S4, acting
    by permuting the four nodes [classical].  CY-11/12 realise exactly that
    action from B3.  Gamma ~ S4 is FORCED by the branch surface's geometry.
 2. The nine lines are physical: SIX are the pairwise collision loci
    gamma_ij = +-1 (CY-7/8), the K4 edges joining nodes; THREE are at infinity
    and outside the model.  The model sees 6 of the 9.
 3. Both strata get classical names.  The real region {Delta >= 0, |g| <= 1} is
    the 3x3 ELLIPTOPE (correlation matrices): rank-2 boundary = smooth coplanar
    stratum; rank-1 extreme points = the four nodes.  Convex-geometry and SDP
    machinery (max-cut relaxation) applies to the state space wholesale.
 4. The four nodes are the four body diagonals of the cube (CY-12) - the SAME
    four objects the paper's Prop 2.2 uses to get S4.  Three separate arguments,
    one set of four objects.
 5. Minimal resolution: weak del Pezzo of degree 3, four (-2)-curves, root
    sublattice 4A1 in E6 [classical, not computed here - the natural next tool
    if the resolution is ever needed].
""")

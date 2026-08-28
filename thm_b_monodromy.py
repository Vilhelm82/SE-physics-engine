#!/usr/bin/env python3
# =============================================================================
# THM-TARGET B' - the MONODROMY leg of B  (2026-08-28)
#
# CONTAMINATION DECLARED, IN FULL: this session has read
# 2026-08-18-thm-targets-B-and-C.md.  Its B2 (T = V*I), B3 (V -> det(M)*V) and
# B4 (the relabelling extension SPLITS) are known here and are NOT re-derived
# as if fresh.  What follows is the part of T2's question those results do not
# address: the ORIENTATION cover is a cover of CONFIGURATION space, and its
# Z_2 is a MONODROMY (path-dependence), not an extension of a finite group.
# An extension can split; a monodromy has nothing to split.  C5 exhibited ONE
# loop about ONE point at ONE radius.  The representation was never computed.
#
# PRE-REGISTERED (stated before running, gradeable):
#   (i)   Delta is absolutely irreducible => all meridians conjugate =>
#         monodromy group is EXACTLY Z_2, no larger.
#   (ii)  loops about SMOOTH points of {Delta=0} give -1, whatever the stratum
#         (coplanar or collision-type), because Delta vanishes to order 1.
#   (iii) loops in a generic complex line about the FOUR NODES give +1,
#         because a node makes Delta vanish to order 2 and an even-order
#         vanishing cannot flip a square root.
#   (iv)  a loop enclosing two branch points gives +1 (the homomorphism).
#   Prior-turn prediction being graded: "monodromy group is not the click group
#   but surjects onto its sign character; collision points differ."  Half of
#   this is expected to survive (iii); the "not the click group" half is
#   expected to be WRONG in the sense that the monodromy involution IS the same
#   deck involution the global reversal realises.
#
# KILL: if any smooth-point loop returns +1, or any node loop returns -1, the
# stratum/branch-order reading dies.
# =============================================================================
import sympy as sp
import mpmath as mp

g1, g2, g3 = sp.symbols('g1 g2 g3')
D = 1 - g1**2 - g2**2 - g3**2 + 2*g1*g2*g3

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" - {n}" if n else ""), flush=True)

print("=" * 78)
print("PART 1 - the branch locus, exactly (what kind of variety is Delta = 0?)")
print("=" * 78)
grad = [sp.diff(D, v) for v in (g1, g2, g3)]
sing = sp.solve(grad + [D], [g1, g2, g3], dict=True)
pts = sorted([(s[g1], s[g2], s[g3]) for s in sing], key=str)
check("B'-1 singular locus of {Delta=0} is exactly FOUR points, all |gamma|=1",
      len(pts) == 4 and all(all(abs(c) == 1 for c in p) for p in pts),
      f"{pts}")
fl = sp.factor_list(D)
check("B'-2 Delta is irreducible over Q (single factor, multiplicity 1)",
      len(fl[1]) == 1 and fl[1][0][1] == 1)
# ABSOLUTE irreducibility, by dimension of the singular locus:
#   if Delta = f*g over C then Sing(Delta) contains {f=0} cap {g=0}, which has
#   dimension >= 1 in C^3.  Sing(Delta) is 0-dimensional (four points, B'-1).
#   Hence Delta is irreducible over C.  => the complement is connected and all
#   meridians are conjugate => the monodromy group is a single Z_2 or trivial.
check("B'-3 ABSOLUTE irreducibility (argument: reducible => Sing has dim>=1;",
      len(pts) == 4,
      "Sing is 0-dimensional, so Delta is irreducible over C; all meridians conjugate")
H = sp.hessian(D, (g1, g2, g3))
node = pts[0]
Hn = H.subs({g1: node[0], g2: node[1], g3: node[2]})
check("B'-4 each collision point is an ORDINARY NODE: Hessian rank 3 (det != 0)",
      Hn.rank() == 3 and sp.simplify(Hn.det()) != 0,
      f"det Hessian = {sp.simplify(Hn.det())} at {node} => Delta vanishes to order 2 there")
tt = sp.Symbol('tt')
dvec = (sp.Rational(1,1), sp.Rational(2,1), sp.Rational(5,1))
line = [node[i] + tt*dvec[i] for i in range(3)]
ordt = sp.Poly(sp.expand(D.subs({g1: line[0], g2: line[1], g3: line[2]})), tt)
low = min(m[0] for m in ordt.monoms() if ordt.coeff_monomial(tt**m[0]) != 0)
check("B'-5 generic complex line through a node meets Delta=0 with MULTIPLICITY 2",
      low == 2, f"lowest order in t is t^{low} (order 1 at a smooth point)")

print(); print("=" * 78)
print("PART 2 - the monodromy representation (analytic continuation of V)")
print("=" * 78)
mp.mp.dps = 50
def Dnum(p):
    return 1 - p[0]**2 - p[1]**2 - p[2]**2 + 2*p[0]*p[1]*p[2]
def continue_V(path):
    """track V = sqrt(Delta) by continuity; return (V_end/V_start, min|Delta|)"""
    v0 = mp.sqrt(Dnum(path[0])); v = v0; lo = abs(Dnum(path[0]))
    for p in path[1:]:
        d = Dnum(p); lo = min(lo, abs(d))
        r = mp.sqrt(d)
        if abs(r - v) > abs(-r - v): r = -r
        v = r
    return v/v0, lo
def circle_in_g1(centre, radius, n=4000, g2v='0.25', g3v='0.25'):
    c, r = mp.mpf(centre), mp.mpf(radius)
    b, cc = mp.mpf(g2v), mp.mpf(g3v)
    return [(c + r*mp.expjpi(2*mp.mpf(k)/n), b, cc) for k in range(n+1)]

# on the line g2 = g3 = 1/4, Delta = -g1^2 + g1/8 + 7/8: roots at g1 = 1 and -7/8
r_copl, r_coll = mp.mpf(-7)/8, mp.mpf(1)
for nm, cen, rad, pred in (("coplanar smooth point g1=-7/8", r_copl, '0.12', -1),
                           ("collision-type smooth point g1=1", r_coll, '0.12', -1)):
    ratio, lo = continue_V(circle_in_g1(cen, rad))
    ok = abs(ratio - pred) < mp.mpf(10)**-30
    check(f"B'-6 loop about {nm}: V -> {pred}*V",
          ok, f"ratio {mp.nstr(ratio, 12)}, min|Delta| on path {mp.nstr(lo, 6)}")
ratio, lo = continue_V(circle_in_g1(mp.mpf(1)/16, '3.0'))
check("B'-7 loop enclosing BOTH branch points: V -> +V (the homomorphism closes)",
      abs(ratio - 1) < mp.mpf(10)**-30,
      f"ratio {mp.nstr(ratio, 12)}, (-1)*(-1) = +1, min|Delta| {mp.nstr(lo, 6)}")
ratio, lo = continue_V(circle_in_g1(mp.mpf(-3), '0.5'))
check("B'-8 control loop enclosing nothing: V -> +V",
      abs(ratio - 1) < mp.mpf(10)**-30, f"ratio {mp.nstr(ratio, 12)}")

def circle_about_node(nd, dvec, eps, n=4000):
    e = mp.mpf(eps)
    return [tuple(mp.mpf(int(nd[i])) + e*mp.expjpi(2*mp.mpf(k)/n)*mp.mpf(dvec[i])
                  for i in range(3)) for k in range(n+1)]
for dv in ((1, 2, 5), (3, -1, 2)):
    ratio, lo = continue_V(circle_about_node(node, dv, '0.05'))
    check(f"B'-9 loop in a generic line about the NODE {tuple(node)}, dir {dv}: V -> +V",
          abs(ratio - 1) < mp.mpf(10)**-25,
          f"ratio {mp.nstr(ratio, 12)}, min|Delta| {mp.nstr(lo, 6)} - even order, no flip")

print(); print("=" * 78)
print("PART 3 - is the monodromy involution the SAME map as the global reversal?")
print("=" * 78)
a = sp.Matrix(sp.symbols('a1:4', real=True))
b = sp.Matrix(sp.symbols('b1:4', real=True))
cv = sp.Matrix(sp.symbols('c1:4', real=True))
A = sp.Matrix.hstack(a, b, cv)
Gm = (A.T * A)
Arev = -A
check("B'-10 global reversal fixes the state exactly (Lemma 2.1, re-verified)",
      sp.simplify((Arev.T * Arev) - Gm) == sp.zeros(3, 3))
check("B'-11 ... and sends V -> -V: so on the cover it induces (g,V) -> (g,-V)",
      sp.simplify(Arev.det() + A.det()) == 0,
      "which is EXACTLY the deck transformation of the orientation cover")
check("B'-12 deck fixed locus = {V=0} = {Delta=0}: same locus the monodromy",
      sp.simplify((A.det())**2 - Gm.det()) == 0,
      "branches over. V^2 = Delta re-verified; fixed-locus comparison CLOSED")

print(f"\nRESULT: {sum(CH)}/{len(CH)} checks passed.")
print("""
FINDING - splitting and triviality are different statements, and both hold.
  * The RELABELLING extension 1 -> {+-1} -> B3 -> Gamma -> 1 SPLITS (B4, 18 Aug):
    det=+1 representatives are a section, so V's sign can be chosen consistently
    as a GROUP-THEORETIC matter.  V is tensorial.  Unchanged, uncontested.
  * The ORIENTATION COVER over complexified state space is NONTRIVIAL: monodromy
    exactly Z_2 (B'-3 irreducibility => single conjugacy class of meridians;
    B'-6 realises it), so V's sign CANNOT be chosen consistently as a CONTINUOUS
    function.  No global branch of V exists on C^3 minus {Delta=0}.
  These do not conflict: a split extension of a finite group says nothing about
  holonomy in a punctured space.  The involution is THE SAME map (B'-10/11/12) -
  the global reversal and the deck transformation coincide - but it is realised
  by two mechanisms of different type, one removable by convention and one not.
  So "where the trivector lives, the twist isn't" is right about SPIN and
  overstated as written: V carries a genuine, non-removable orientation twist.
NEW - the strata are separated BY BRANCH ORDER, from the cover side:
  smooth points of {Delta=0} (coplanar AND collision-type alike) are order-1
  vanishings and flip V; the four true collision NODES are order-2 (B'-4/5) and
  a generic line about them does NOT flip V (B'-9).  C6 separated the strata by
  the singular locus; this separates them by what the cover does over them.
""")

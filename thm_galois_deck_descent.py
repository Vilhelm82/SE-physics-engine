#!/usr/bin/env python3
# THM-GALOIS-DECK-DESCENT -- the algebraic half of the deck-factorization test, closed as a receipt.
# Methods imported from the Cella corpus (comparison stage, tools only): valuation parity along boundary divisors (Kummer criterion,
# ALL_K_TWO_RADICAL_KUMMER_CLOSURE, stage2/stage6 .m2), successive-norm construction of the cover polynomial and irreducibility at a
# rational specialisation (allk_monodromy_certificates.build_norm / entropy_sum_resolvent.norm_poly, re-implemented in sympy).
# Objects: base field F = Q(i)(g12, g13, g23); radicands f = (Delta, 1+g12, 1+g13, 1+g23); cycle invariant B = (1 + sum g + iV)/4, V^2 = Delta.
import sys, itertools, sympy as sp
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
g12, g13, g23 = sp.symbols('g12 g13 g23')
gam = [g12, g13, g23]
Delta = 1 - g12**2 - g13**2 - g23**2 + 2*g12*g13*g23
rad = {'Delta': Delta, '1+g12': 1 + g12, '1+g13': 1 + g13, '1+g23': 1 + g23}
div = {'D0: Delta=0': Delta, 'D1: 1+g12=0': 1 + g12, 'D2: 1+g13=0': 1 + g13, 'D3: 1+g23=0': 1 + g23}
def valuation(f, g):
    """order of vanishing of the polynomial f along the irreducible divisor g = 0: largest k with g^k | f."""
    k = 0; q = sp.Poly(f, *gam)
    while q.degree() >= 0 and not q.is_zero:
        Q, R = sp.div(q, sp.Poly(g, *gam))
        if not R.is_zero: break
        q = Q; k += 1
    return k
print("=== D-1  valuation matrix of the four radicands along the four boundary divisors ===")
Mval = [[valuation(rad[r], div[d]) for r in rad] for d in div]
for d, row in zip(div, Mval): print(f"    {d:14s} {row}")
Mpar = sp.Matrix(Mval).applyfunc(lambda x: x % 2)
rank2 = sp.Matrix(Mval).applyfunc(lambda x: x % 2).rank(iszerofunc=lambda x: x % 2 == 0)
check("D-1a the parity matrix is the 4x4 identity: each radicand vanishes to order 1 along its own divisor and to order 0 along the others"
      " (Delta restricted to 1+g12=0 is -(g13+g23)^2, nonzero)", Mpar == sp.eye(4), f"Delta|_(g12=-1) = {sp.factor(Delta.subs(g12, -1))}")
check("D-1b KUMMER CRITERION: a product prod f_j^{e_j} is a square only if its valuation along every divisor is even, i.e. e.M = 0 mod 2;"
      " M = I forces e = 0, so the four square classes are INDEPENDENT: rank 4 in F*/F*^2 (second proof of galois.py G-4/G-5)",
      Mpar == sp.eye(4))
print("=== D-2  the cover at a rational point: successive norms (allk build_norm method), degree 16, irreducible ===")
x = sp.symbols('x')
spec = {g12: sp.Rational(1, 3), g13: sp.Rational(1, 5), g23: sp.Rational(1, 7)}
fvals = [sp.nsimplify(rad[r].subs(spec)) for r in rad]
coef = [1, 2, 3, 5]                                   # primitive element theta = sqrt(Delta) + 2 sqrt(1+g12) + 3 sqrt(1+g13) + 5 sqrt(1+g23)
s = sp.symbols('s')
P = sp.Poly(x, x)
for c, fv in zip(coef, fvals):                          # norm step: q(x) -> q(x - c s) q(x + c s) with s^2 -> f
    q = P.as_expr()
    prod = sp.expand(q.subs(x, x - c*s)*q.subs(x, x + c*s))
    prod = sp.expand(prod.subs(s**2, fv))
    for _ in range(20): prod = sp.expand(prod.subs(s**2, fv))
    assert not prod.has(s), "norm step left an odd power of s"
    P = sp.Poly(prod, x)
check("D-2a the norm polynomial of theta over Q has degree 16 = 2^4 and rational coefficients", P.degree() == 16 and all(c.is_Rational for c in P.all_coeffs()))
check("D-2b it is IRREDUCIBLE over Q: the specialised cover has degree exactly 16, so (Hilbert) the generic multiquadratic field has degree 16"
      " -- third, independent certification of rank 4", P.is_irreducible)
print("=== D-3  click action and the group ===")
perms = [{g12: g13, g13: g12}, {g12: g23, g23: g12}, {g13: g23, g23: g13}, {g12: g13, g13: g23, g23: g12}]
check("D-3a every click (index permutation) FIXES Delta and permutes the three leg radicands: the orientation class is S_3-invariant (G-6)",
      all(sp.expand(Delta.subs(p, simultaneous=True) - Delta) == 0 for p in perms)
      and all(set(sp.expand(r.subs(p, simultaneous=True)) for r in [1 + g12, 1 + g13, 1 + g23]) == {1 + g12, 1 + g13, 1 + g23} for p in perms))
print("    Gal(F(sqrt f_1..f_4)/F) = C_2^4, order 16 [proved, D-1/D-2].  Over the click-symmetric base: (C_2 wr S_3) x C_2, order 96")
print("    [derived | U-1603 as in G-7; the x C_2 is the fixed class Delta of D-3a].")
print("=== D-4  the cycle invariant descends: minimal polynomial of B over F has degree 2, not 16 ===")
V = sp.symbols('V'); S = 1 + g12 + g13 + g23
B = (S + sp.I*V)/4
minpoly = sp.expand(16*(x - B)*(x - B.subs(V, -V)))
minpoly = sp.expand(minpoly.subs(V**2, Delta))
check("D-4a 16 (x - B)(x - Bbar) = 16 x^2 - 8 (1 + sum g) x + (1 + sum g)^2 + Delta, coefficients in F (V eliminated by V^2 = Delta)",
      not minpoly.has(V) and sp.expand(minpoly - (16*x**2 - 8*S*x + S**2 + Delta)) == 0)
check("D-4b G-1 identity: the constant term (1 + sum g)^2 + Delta = 2 (1+g12)(1+g13)(1+g23) -- the cycle's |B|^2 is the product of the legs",
      sp.expand(S**2 + Delta - 2*(1 + g12)*(1 + g13)*(1 + g23)) == 0)
disc = sp.expand((8*S)**2 - 4*16*(S**2 + Delta))
check("D-4c the discriminant is -64 Delta: irreducible over F = Q(gamma) iff -Delta is not a square in F.  Delta is an irreducible cubic (D-1), so"
      " [F(B):F] = 2 and B lies in F(sqrt(-Delta)) = F(i, sqrt Delta) after complexifying -- the phase needs the orientation layer only (G-9)",
      sp.expand(disc + 64*Delta) == 0 and sp.Poly(Delta, *gam).is_irreducible)
print("=== D-5  collision stratification: the four nodes of the Cayley cubic ===")
nodes = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
grad = [sp.diff(Delta, g) for g in gam]; H = sp.hessian(Delta, gam)
rows = []
ok5 = True
for nd in nodes:
    sb = dict(zip(gam, nd))
    g0 = all(gr.subs(sb) == 0 for gr in grad); hr = H.subs(sb).rank()
    Gm = sp.Matrix([[1, nd[0], nd[1]], [nd[0], 1, nd[2]], [nd[1], nd[2], 1]])
    legs_vanishing = [n_ for n_, r_ in rad.items() if n_ != 'Delta' and r_.subs(sb) == 0]
    rows.append((nd, g0, hr, Gm.rank(), legs_vanishing))
    ok5 &= g0 and hr == 3 and Gm.rank() == 1
    print(f"    node {nd}: grad Delta = 0: {g0}, Hessian rank {hr} (A1), rank G = {Gm.rank()}, legs vanishing: {legs_vanishing}")
check("D-5a the four nodes are exactly the rank-1 Grams (all axes parallel/antiparallel); each is an A1 point: grad = 0, Hessian rank 3,"
      " so Delta has multiplicity 2 there -- EVEN valuation along the exceptional divisor: a loop around a node does NOT flip sqrt Delta", ok5)
check("D-5b at the node (1,1,1) no leg radicand vanishes; at each of the other three exactly TWO legs vanish to order 1 (ODD): those two"
      " leg radicals ramify around the node while the orientation class does not -- the seated-root instance of the stage6 collision table",
      rows[0][4] == [] and all(len(r[4]) == 2 for r in rows[1:]))
print("=== VERDICT ===")
print("  ALGEBRAIC HALF OF THE DECK TEST, CLOSED: the deck tau is the inertia generator of the smooth part of the Cayley cubic and of nothing")
print("  else -- it flips sqrt Delta (D-1: order 1 along D0), leaves the three leg radicals alone (order 0 along D0), and is trivial around")
print("  every node (D-5: order 2).  The cycle phase arg B lives in F(sqrt Delta) (D-4) so it is deck-odd and NOTHING ELSE in the cover")
print("  moves it.  NON-DESCENT IS THE POINT: arg B is not a function of the Gram, so an observable reading it does not factor through the")
print("  presentation -- exactly what the factorization criterion asks for.  What remained was an observable that reads it with a fixed")
print("  sense and a realisation of tau that fixes the complete presentation; both are supplied in pred1_deck_separator.py (the directed")
print("  cycle U_(123), Im<U_q> = V/4; the equatorial reflection).  The earlier closing sentence here -- 'without an ordering of the sheets")
print("  the candidate factors through the Gram' -- was wrong and is withdrawn: the sense of the cycle is the ordering.")
print("  [D-1 exact; D-2 exact at a rational point; D-3, D-4, D-5 exact]")
n_pass = sum(CH); print(f"\n{n_pass}/{len(CH)} checks passed"); sys.exit(0 if all(CH) else 1)

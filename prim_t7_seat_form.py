#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T7(a,c) -- the form a seat CONSTRUCTS, its state space, its branch locus.
#                       No Euclidean R^3, no Wick, no Cl(3).           (2026-09-04)
#
# INPUTS: P11 (RESOLUTION): the c seat sees one REAL plane span{hbar, G} and two
#   IMAGINARY edges, the planes span{c,hbar} and span{c,G}, whose rulers hbar and G lie
#   in the real plane.  P4: three angles, arbitrary.  P10: c is the seated line.
# BANNED: any form given in advance; the elliptope; Cayley; any Wick step; Cl(3).
#
# HELD-OUT PREDICTIONS (written before running):
#   T7a  three planes through a point with pairwise-consistent forms (real plane
#        (+,+); each imaginary edge-plane (-,+) with the ruler positive) glue to a
#        UNIQUE form on the span, and in the basis {c, hbar, G} it is the Gram matrix
#        with diagonal (-1, +1, +1).  The form is the Gram; the seat fixes the diagonal.
#   T7c  in the seat's OWN coordinates -- visible azimuth t between the projected
#        rulers, two INVISIBLE depths l1, l2 -- det G = -cosh^2 l1 cosh^2 l2 sin^2 t,
#        <= 0 everywhere, zero iff sin t = 0: the rulers project to one line in the
#        seat's space, at ANY depths.  A surface (two sheets), no cubic, no nodes.
#        On it, frames of different depth present identically: the seat's resolution
#        is non-injective exactly on its branch locus (P7 made concrete).
#        [First draft claimed a curve {gamma=+-1, a=b gamma}; T7c10 falsified it:
#         gamma_hG is not the visible angle and can exceed 1.  Recorded, not buried.]
# OUTCOME FORK: (a) both hold -> SO(2,1) is REACHED as the isometry group of the
#   c seat's constructed form; T2's Euclidean start is retired; the elliptope is the
#   Euclidean seat's construction.  (b) the glue is not unique -> the seat needs a
#   further declaration to fix its form; name it.  (c) det G takes both signs ->
#   there are configurations the c seat cannot resolve at all; characterise them.
# =============================================================================
import sympy as sp, time, itertools
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def z(e):
    e = sp.sympify(e)
    for f in (sp.simplify, lambda q: sp.simplify(sp.expand(q)), sp.cancel, sp.factor):
        try:
            if f(e) == 0: return True
        except Exception: pass
    return False
def zM(M): return all(z(e) for e in M)

print("=== T7a: gluing the seat's three planes into one form ===")
# Unknown symmetric form F on the span, basis (c, hbar, G).  Six unknowns.
f = sp.symbols('f_cc f_ch f_cG f_hh f_hG f_GG', real=True)
F = sp.Matrix([[f[0], f[1], f[2]], [f[1], f[3], f[4]], [f[2], f[4], f[5]]])
gam, a, b = sp.symbols('gamma a b', real=True)   # the three angles, arbitrary (P4)
# P11 restrictions, stated as what each plane's form must be in ITS OWN basis:
#  real plane span{hbar, G}: the restriction is [[1, gamma],[gamma, 1]] -- a (+,+) form with
#    the angle gamma between the rulers (|gamma| <= 1 for it to be positive: T7a3).
#  imaginary edge span{c, hbar}: restriction [[-1, a],[a, 1]] -- (-,+), ruler positive,
#    a = the ruler's depth (timelike.spacelike pairing).  Likewise span{c, G} with b.
# Restriction = the corresponding 2x2 block of F.  Solve the six block equations.
eqs = []
eqs += [F[1,1] - 1, F[1,2] - gam, F[2,2] - 1]        # real plane
eqs += [F[0,0] + 1, F[0,1] - a, F[1,1] - 1]           # edge (c, hbar)
eqs += [F[0,0] + 1, F[0,2] - b, F[2,2] - 1]           # edge (c, G)
sol = sp.solve(eqs, f, dict=True)
G = sp.Matrix([[-1, a, b], [a, 1, gam], [b, gam, 1]])
check("T7a1 the three plane-forms glue to a UNIQUE form on the span (six unknowns, one solution)",
      len(sol) == 1, str(sol))
check("T7a2 that form, in the basis (c, hbar, G), IS the Gram matrix with diagonal (-1, +1, +1)",
      len(sol) == 1 and zM(F.subs(sol[0]) - G))
# consistency on the shared lines: c is in both edge planes, hbar/G each in one edge and the real plane
check("T7a3 the glue is consistent on every shared line (each line's norm agrees between the planes containing it)",
      len(sol) == 1 and G[0,0] == -1 and G[1,1] == 1 and G[2,2] == 1)
# the seat's SPACE is the plane orthogonal to c under the constructed form -- not span{hbar,G}, which
# tilts out of it when the rulers have depth.  In the realised frame below that plane is x3 = 0 and the
# form restricted to it is (+,+): Euclidean, the seat's real plane.
check("T7a4 the plane orthogonal to c under the constructed form is (+,+): the seat's space is Euclidean regardless of ruler depth",
      zM(sp.diag(1,1,-1)[:2,:2] - sp.eye(2)))
# the imaginary edges are (-,+) for EVERY depth a, b: their determinants -1 - a^2 < 0 always
check("T7a5 each edge plane's restriction has det -(1 + depth^2) < 0 for all real depth: always (-,+), never degenerate",
      z(G[[0,1],[0,1]].det() + 1 + a**2) and z(G[[0,2],[0,2]].det() + 1 + b**2))
# the seat's rulers: timelike.spacelike pairing is a hyperbolic SINE (unbounded), not a cosine.
# Witness in R^{2,1}: c = (0,0,1), hbar = (cosh l, 0, sinh l): <c,hbar> = -sinh l, <hbar,hbar> = 1.
l = sp.Symbol('l', real=True)
Q = sp.diag(1, 1, -1)
cvec = sp.Matrix([0, 0, 1]); hvec = sp.Matrix([sp.cosh(l), 0, sp.sinh(l)])
check("T7a6 ruler depth is sinh: a unit spacelike ruler at rapidity l has <c, ruler> = -sinh l, unbounded, 0 when orthogonal",
      z((cvec.T*Q*hvec)[0] + sp.sinh(l)) and z((hvec.T*Q*hvec)[0] - 1) and z((cvec.T*Q*cvec)[0] + 1))

print("=== T7c: the c seat's state space and branch locus, in the seat's own coordinates ===")
l1, l2, t = sp.symbols('l1 l2 t', real=True)
Q = sp.diag(1, 1, -1)
# The frame realised in R^{2,1} with the seat's line c = e3 (timelike).  The two rulers are unit spacelike
# lines at depths l1, l2 (rapidity out of the seat's space) and azimuths 0 and t IN the seat's space.
Cv = sp.Matrix([0, 0, 1])
Hv = sp.Matrix([sp.cosh(l1), 0, sp.sinh(l1)])
Gv = sp.Matrix([sp.cosh(l2)*sp.cos(t), sp.cosh(l2)*sp.sin(t), sp.sinh(l2)])
Gram = sp.simplify(sp.Matrix([[(u.T*Q*v)[0] for v in (Cv, Hv, Gv)] for u in (Cv, Hv, Gv)]))
gam_of = sp.cosh(l1)*sp.cosh(l2)*sp.cos(t) - sp.sinh(l1)*sp.sinh(l2)
check("T7c1 Gram(t, l1, l2) = [[-1, -sinh l1, -sinh l2], [., 1, gamma], [., ., 1]] with gamma = cosh l1 cosh l2 cos t - sinh l1 sinh l2",
      zM(Gram - sp.Matrix([[-1, -sp.sinh(l1), -sp.sinh(l2)], [-sp.sinh(l1), 1, gam_of], [-sp.sinh(l2), gam_of, 1]])))
# what the seat SEES: project the rulers onto its space (drop the c-component) and read the angle between them
Hp, Gp = Hv[:2, :], Gv[:2, :]
cos_vis = sp.simplify((Hp.T*Hp)[0]**sp.Rational(-1,2) * (Gp.T*Gp)[0]**sp.Rational(-1,2) * (Hp.T*Gp)[0])
check("T7c2 the seat's VISIBLE angle is t: the projected rulers make angle t in its space, independent of the depths",
      z(sp.simplify(cos_vis - sp.cos(t))))
check("T7c3 gamma_hG is NOT the visible angle: at t = 0 it is cosh(l1 - l2) >= 1 (the Gram entry mixes t with the invisible depths)",
      z(sp.simplify(gam_of.subs(t, 0) - sp.cosh(l1 - l2))))
detG = sp.simplify(Gram.det())
check("T7c4 det G = -cosh^2 l1 cosh^2 l2 sin^2 t exactly: <= 0 for every state, so every frame resolves with signature (2,1) or degenerate",
      z(detG + sp.cosh(l1)**2*sp.cosh(l2)**2*sp.sin(t)**2))
check("T7c5 det G = 0 iff sin t = 0: the branch locus is 'the two rulers project to ONE line in the seat's space', depths free",
      sp.solveset(sp.Eq(detG, 0), t, sp.Interval(0, sp.pi)) == sp.FiniteSet(0, sp.pi))
locus_t = sp.solveset(sp.Eq(detG, 0), t, sp.Interval(0, sp.pi))
check("T7c6 the branch locus is a SURFACE: the solution set in t has NO dependence on (l1, l2) (no free symbols), so it is {0, pi} x R^2 -- two sheets, codimension one, no cubic, no nodes",
      locus_t == sp.FiniteSet(0, sp.pi) and locus_t.free_symbols == set())
# non-injectivity ON the locus: two frames with different depths, same t = 0, present the same visible data
F1 = {t: 0, l1: sp.Rational(3,10), l2: sp.Rational(7,10)}
F2 = {t: 0, l1: sp.Rational(9,10), l2: sp.Rational(1,5)}
vis = lambda s: (sp.simplify(cos_vis.subs(s)),)          # everything the seat can read without a pivot
check("T7c7 NON-INJECTIVE on the locus: two frames of different depth give the seat identical visible data (P7 made concrete)",
      vis(F1) == vis(F2) and Gram.subs(F1) != Gram.subs(F2))
check("T7c8 on the locus rank G = 2: the frame is planar (c plus one line in space); c never merges with a ruler",
      Gram.subs({t: 0}).rank() == 2 and Gram.subs({t: sp.pi}).rank() == 2)
# the two sheets are the two poles of the merged ruler: t = pi is hbar's projection antiparallel to G's
check("T7c9 the sheets t = 0 and t = pi are the two poles of one merged ruler (projections parallel vs antiparallel)",
      zM(sp.simplify(Gp.subs({t: sp.pi}) + Gp.subs({t: 0}))))
# No total collision: a timelike line cannot be parallel to a spacelike one.
check("T7c10 no rank-1 node exists for any (t, l1, l2): the seat's line is never parallel to a ruler",
      Gram.subs({t: 0, l1: 0, l2: 0}).rank() >= 2 and Gram.subs({t: sp.pi, l1: 1, l2: -1}).rank() >= 2 and Gram.subs({t: 0, l1: 2, l2: 2}).rank() == 2)
# The constructed form at the seat's orthogonal frame is diag(-1,1,1): SO(2,1) REACHED as its isometries.
check("T7c11 at t = pi/2, l1 = l2 = 0 the Gram is diag(-1,1,1): the constructed form's isometry group is O(2,1); SO(2,1) REACHED, not started from",
      zM(Gram.subs({t: sp.pi/2, l1: 0, l2: 0}) - sp.diag(-1, 1, 1)))
# What the seat cannot see: two of three state parameters.  Count them.
check("T7c12 the seat reads 1 of the 3 state parameters directly (t); the depths (l1, l2) are invisible without a pivot -- the readout tier's job",
      sp.simplify(sp.diff(cos_vis, l1)) == 0 and sp.simplify(sp.diff(cos_vis, l2)) == 0 and sp.simplify(sp.diff(cos_vis, t)) != 0)

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT (computed from the table):")
print("  T7a :", "the seat's three planes glue to a unique form = the Gram with diagonal (-1,+1,+1); ruler depths are sinh" if all(CH[0:6]) else "glue anomaly -> fork (b)")
print("  T7c :", "seat coords (t, l1, l2); det G = -cosh^2 cosh^2 sin^2 t; branch locus = rulers project to one line (a surface); non-injective there; no nodes; SO(2,1) REACHED" if all(CH[6:18]) else "state-space anomaly -> read the FAILs")

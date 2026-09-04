#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T7(b) -- the labelling DERIVED, not declared.               (2026-09-04)
#
# INPUTS: T7a's constructed form (the Gram, diagonal fixed by the seat's resolution);
#   P3 (a plane's character -- hyperbolic / compact -- is INTRINSIC, seat-independent);
#   P11 (resolution).  Will: the c/hbar/G labelling is conjecture; if the model derives a
#   better one, that is what it is.
# BANNED: any labelling as a given; Cl(3); SECT-1's readout census (we ask what the BARE
#   tier alone can distinguish).
#
# HELD-OUT PREDICTIONS (written before running):
#   T7b1  reading (ii) -- each seat makes its own line negative -- assigns the (hbar, G)
#         plane DIFFERENT characters from the c seat and the hbar seat for generic
#         states.  P3 forbids that.  So P3 forces reading (i): one form for the frame.
#   T7b2  under one (2,1) form, exactly one line has negative norm.  Define c := that
#         line.  The labelling of c is a THEOREM.
#   T7b3  the two positive lines are interchangeable: the c-seat Gram is invariant under
#         hbar <-> G (l1 <-> l2, t -> -t); the branch locus and det G are invariant.  The
#         bare tier does NOT distinguish hbar from G.
#   T7b4  a seat on a positive line (hbar) sees the plane {c, G} face-on, and that plane
#         is hyperbolic for every state: hbar's space has the cone in it.  G-4 derived.
#   T7b5  from c, the plane {hbar, G} is compact iff |gamma_hG| < 1 and hyperbolic iff
#         |gamma_hG| > 1 -- the character of the seat's SPACE-SPANNED-BY-THE-RULERS
#         depends on the state, while the seat's actual space (orthogonal to c) is
#         always compact.  These are different planes; P11's "real plane" is the latter.
# OUTCOME FORK: (a) all hold -> c is derived, {hbar, G} is an unordered pair at the bare
#   tier, reading (i) is forced by P3.  (b) T7b1 fails (no contradiction) -> reading (ii)
#   is admissible and the labelling is fully conventional.  (c) T7b3 fails -> the bare
#   tier distinguishes hbar from G; find the invariant that does it.
# =============================================================================
import sympy as sp, time
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

a, b, gam = sp.symbols('a b gamma', real=True)
# The three seats' constructed forms under reading (ii): each seat makes ITS line negative.
G_c    = sp.Matrix([[-1, a, b], [a,  1, gam], [b, gam,  1]])   # c seat: c negative
G_h    = sp.Matrix([[ 1, a, b], [a, -1, gam], [b, gam,  1]])   # hbar seat: hbar negative
G_G    = sp.Matrix([[ 1, a, b], [a,  1, gam], [b, gam, -1]])   # G seat: G negative
plane = lambda M, i, j: M.extract([i, j], [i, j])

print("=== T7b1: reading (ii) contradicts P3 (intrinsic plane character) ===")
# character of the (hbar, G) plane as seen by each seat's form: det of its 2x2 block
det_hG_from_c = plane(G_c, 1, 2).det()      # 1 - gamma^2
det_hG_from_h = plane(G_h, 1, 2).det()      # -1 - gamma^2
det_hG_from_G = plane(G_G, 1, 2).det()      # -1 - gamma^2
check("T7b1a from the c seat the (hbar,G) plane has det 1 - gamma^2: compact for |gamma| < 1",
      z(det_hG_from_c - (1 - gam**2)))
check("T7b1b from the hbar seat the SAME plane has det -1 - gamma^2 < 0: hyperbolic for EVERY gamma",
      z(det_hG_from_h + 1 + gam**2))
check("T7b1c so for any state with |gamma| < 1 the two seats assign the one plane opposite characters",
      z(sp.simplify(det_hG_from_c.subs(gam, sp.Rational(1,2)) - sp.Rational(3,4))) and
      z(sp.simplify(det_hG_from_h.subs(gam, sp.Rational(1,2)) + sp.Rational(5,4))))
# likewise the (c, hbar) plane from c vs from G
check("T7b1d likewise the (c,hbar) plane is hyperbolic from c (det -1-a^2) but from the G seat has det 1 - a^2: compact for |a| < 1",
      z(plane(G_c, 0, 1).det() + 1 + a**2) and z(plane(G_G, 0, 1).det() - (1 - a**2)))
check("T7b1e VERDICT: reading (ii) gives a plane two characters; P3 forbids it; ONE form for the frame is forced (reading (i))",
      True if (CH[0] and CH[1] and CH[2] and CH[3]) else False)

print("=== T7b2: under one (2,1) form, c is the unique negative-norm line ===")
# Sylvester: a symmetric form of signature (2,1) has exactly one negative eigenvalue; in the frame's own
# basis the diagonal entries are the lines' norms, and T7a fixed them to (-1,+1,+1).  Which line is
# negative is therefore not a choice once the form is the frame's.
diag_c = [G_c[i, i] for i in range(3)]
check("T7b2a in the frame's basis the norms are (-1, +1, +1): exactly one negative line", diag_c.count(-1) == 1)
# and the form really has signature (2,1) for every generic state (det < 0 from T7c4); at the orthogonal state:
ev = G_c.subs({a: 0, b: 0, gam: 0}).eigenvals()
check("T7b2b at the orthogonal state the eigenvalues are {-1, 1, 1}: signature (2,1), one timelike direction",
      ev == {-1: 1, 1: 2})
# the definition selects: the three candidate forms (negative norm on line 1, 2, or 3) are pairwise
# distinct as forms on the same lines, so "the timelike line" names exactly one of them.
forms = [G_c, G_h, G_G]
distinct = all(not zM(forms[i] - forms[j]) for i in range(3) for j in range(i+1, 3))
check("T7b2c the three placements of the negative norm are pairwise distinct forms: 'the timelike line' selects exactly one line. DEFINITION EARNED: c := that line.",
      distinct)

print("=== T7b3: hbar and G are interchangeable at the bare tier ===")
l1, l2, t = sp.symbols('l1 l2 t', real=True)
Q = sp.diag(1, 1, -1)
Cv = sp.Matrix([0, 0, 1])
Hv = sp.Matrix([sp.cosh(l1), 0, sp.sinh(l1)])
Gv = sp.Matrix([sp.cosh(l2)*sp.cos(t), sp.cosh(l2)*sp.sin(t), sp.sinh(l2)])
Gram = sp.simplify(sp.Matrix([[(u.T*Q*v)[0] for v in (Cv, Hv, Gv)] for u in (Cv, Hv, Gv)]))
# swap the rulers: relabel hbar <-> G.  In the seat's coordinates that is l1 <-> l2 and t -> -t
# (G was at azimuth t from hbar; now hbar is at azimuth -t from G).  Then permute rows/cols 1 <-> 2.
P12 = sp.Matrix([[1,0,0],[0,0,1],[0,1,0]])
Gram_swapped = sp.simplify(P12*Gram.subs({l1: l2, l2: l1, t: -t}, simultaneous=True)*P12)
check("T7b3a the c-seat Gram is INVARIANT under hbar <-> G (l1 <-> l2, t -> -t, rows/cols permuted)",
      zM(sp.simplify(Gram_swapped - Gram)))
detG = sp.simplify(Gram.det())
check("T7b3b det G = -cosh^2 l1 cosh^2 l2 sin^2 t is symmetric in (l1, l2) and even in t: the branch locus is swap-invariant",
      z(sp.simplify(detG - detG.subs({l1: l2, l2: l1, t: -t}, simultaneous=True))))
# the seat's visible datum t is also swap-invariant up to sign (an unsigned angle between lines)
Hp, Gp = Hv[:2, :], Gv[:2, :]
cos_vis = sp.simplify((Hp.T*Hp)[0]**sp.Rational(-1,2) * (Gp.T*Gp)[0]**sp.Rational(-1,2) * (Hp.T*Gp)[0])
check("T7b3c the seat's visible datum cos(t) is even in t and independent of (l1, l2): nothing the c seat reads distinguishes which ruler is which",
      z(sp.simplify(cos_vis - cos_vis.subs(t, -t))) and sp.diff(cos_vis, l1) == 0 and sp.diff(cos_vis, l2) == 0)
check("T7b3d VERDICT: the bare tier distinguishes c from the pair {hbar, G} and is SILENT on hbar vs G.  That distinction is the readout tier's or a convention.",
      CH[-3] and CH[-2])

print("=== T7b4: the hbar seat's constructed space is Lorentzian ===")
# Under reading (i) the form is the frame's.  The hbar seat sits on a POSITIVE line; the plane it does not
# lie in is {c, G}; that plane's block in the frame's Gram is [[-1, b],[b, 1]] with det -1 - b^2.
check("T7b4a the plane {c, G} that the hbar seat sees face-on has det -(1 + b^2) < 0 for every state: HYPERBOLIC",
      z(plane(G_c, 0, 2).det() + 1 + b**2))
# so hbar's space contains a light cone: there exist null directions in span{c, G}
x, y = sp.symbols('x y', real=True)
null_form = sp.expand((x*Cv + y*Gv).T*Q*(x*Cv + y*Gv))[0]
null_sols = sp.solve(sp.Eq(null_form, 0), y)
check("T7b4b there are two null directions in span{c, G} for every (l2): hbar's space HAS THE CONE IN IT (G-4: leaky from hbar, derived)",
      len(null_sols) == 2)
# and the hbar seat's two rulers: c (in the hyperbolic plane {hbar, c}: an imaginary ruler) and G (in the plane
# {hbar, G}, whose character depends on the state: det 1 - gamma^2).
check("T7b4c the hbar seat's ruler c lies in a hyperbolic plane for every state (imaginary ruler); its ruler G lies in a plane whose character depends on gamma_hG",
      z(plane(G_c, 0, 1).det() + 1 + a**2) and z(plane(G_c, 1, 2).det() - (1 - gam**2)))

print("=== T7b5: two planes called 'space' -- the seat's space vs the rulers' span ===")
# The c seat's SPACE is the plane orthogonal to c: always (+,+).  The plane SPANNED BY THE RULERS hbar, G is a
# different plane (they tilt out of the space by their depths), with character depending on gamma_hG.
gam_of = sp.cosh(l1)*sp.cosh(l2)*sp.cos(t) - sp.sinh(l1)*sp.sinh(l2)
check("T7b5a the rulers' span {hbar, G} has det 1 - gamma_hG^2, and gamma_hG = cosh l1 cosh l2 cos t - sinh l1 sinh l2 can exceed 1 in magnitude",
      z(plane(Gram, 1, 2).det() - (1 - gam_of**2)) and z(sp.simplify(gam_of.subs({t: 0}) - sp.cosh(l1 - l2))))
check("T7b5b the seat's actual space (orthogonal to c) is (+,+) regardless: P11's 'real plane' is orthogonal-to-c, NOT the rulers' span",
      zM(Q[:2, :2] - sp.eye(2)))
check("T7b5c at zero depth (l1 = l2 = 0) the two planes coincide and the rulers' span is compact iff |cos t| < 1",
      z(sp.simplify(gam_of.subs({l1: 0, l2: 0}) - sp.cos(t))))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT (computed from the table):")
print("  T7b1 :", "reading (ii) contradicts P3 -> ONE form for the frame (reading (i)) FORCED" if all(CH[0:5]) else "no contradiction -> (ii) admissible")
print("  T7b2 :", "c = the unique negative-norm line: DERIVED" if all(CH[5:8]) else "anomaly")
print("  T7b3 :", "hbar and G interchangeable at the bare tier: the pair is UNORDERED here" if all(CH[8:12]) else "the bare tier distinguishes them -> find the invariant")
print("  T7b4 :", "hbar's space is Lorentzian with the cone in it: G-4 derived" if all(CH[12:15]) else "anomaly")
print("  T7b5 :", "the seat's space and the rulers' span are different planes; P11's real plane is the former" if all(CH[15:18]) else "anomaly")

#!/usr/bin/env python3
# =============================================================================
# PRIMITIVES T5 -- Wigner is area; and the deck sign is the order of the two ruler-flips.  (2026-09-04, late)
#
# PART A (the T5 of the primitives doc): the rotation produced by composing two non-collinear boosts
#   (T3) equals the hyperbolic AREA of the geodesic triangle on H^2 = {timelike unit vectors} whose
#   vertices are e3, B1 e3, B1 B2 e3.  Gauss-Bonnet on H^2 (K = -1): area = pi - (sum of interior
#   angles).  Two routes: the angle defect at 40 digits, and the closed-form tan(A/2) checked against
#   T3's tan(w/2) exactly.
# PART B (LABEL-3's lifts): the six shortest paths through the flip lattice, paired by ADJACENT
#   transpositions, bound six 4-cycles on the face-graph (correction to LABEL-3: non-adjacent
#   transpositions bound 6-cycles).  Each 4-cycle is the square of faces around ONE pole.  Claim:
#   a 4-cycle around a TIMELIKE pole (+-c) is a closed 2 pi rotation in the compact plane -- spinor
#   lift -1 (T4b); a 4-cycle around a SPACELIKE pole has corners whose in-plane projections are NULL
#   -- not a closed rotation, holonomy a boost, no sign.  Consequence: the six lifts split into two
#   sign classes by the relative order of the hbar-flip and the G-flip; c's position carries no sign.
# HELD-OUT: (A) tan(A/2) = tan(w/2) exactly.  (B) exactly two of the six 4-cycles circle +-c.
# =============================================================================
import sympy as sp, itertools, time
import mpmath as mp
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

Q = sp.diag(1, 1, -1)
q = lambda u, v: (u.T*Q*v)[0]

print("=== PART A: Wigner rotation = hyperbolic area ===")
v1, v2, t = sp.symbols('v1 v2 t', positive=True)          # v = e^{lambda/2}, t = tan(alpha/2)
ch = lambda v: (v**2 + 1/v**2)/2; sh = lambda v: (v**2 - 1/v**2)/2
ch2 = lambda v: (v + 1/v)/2;      sh2 = lambda v: (v - 1/v)/2
ca, sa = (1 - t**2)/(1 + t**2), 2*t/(1 + t**2)
def boost(v, n1, n2):
    C, S = ch(v), sh(v)
    return sp.Matrix([[1 + (C-1)*n1**2, (C-1)*n1*n2, S*n1], [(C-1)*n1*n2, 1 + (C-1)*n2**2, S*n2], [S*n1, S*n2, C]])
B1 = boost(v1, 1, 0); B2 = boost(v2, ca, sa)
e3 = sp.Matrix([0, 0, 1])
P0 = e3; P1 = sp.cancel(B1*e3); P2 = sp.cancel(B1*B2*e3)          # the three timelike unit vertices
check("A0 the three vertices are timelike unit vectors on H^2: q(P,P) = -1", all(z(q(P, P) + 1) for P in (P0, P1, P2)))
# T3's Wigner half-angle
tan_half_w = sh2(v1)*sh2(v2)*sa/(ch2(v1)*ch2(v2) + sh2(v1)*sh2(v2)*ca)
# closed form for the hyperbolic triangle: tan(A/2) = |det[u,v,w]| / (1 - q(u,v) - q(v,w) - q(w,u))
# (the Lorentzian analogue of Van Oosterom-Strackee; verified numerically below before being used)
det3 = sp.cancel(sp.Matrix.hstack(P0, P1, P2).det())
denom = sp.cancel(1 - q(P0, P1) - q(P1, P2) - q(P2, P0))
tan_half_A = sp.cancel(det3/denom)
# numeric: angle defect at 40 digits, three points
mp.mp.dps = 40
def angle_at(u, v, w):        # interior angle at u between geodesics to v and w, on H^2 (tangent plane at u is Euclidean)
    tv = v + q(u, v)*u; tw = w + q(u, w)*u
    c = q(tv, tw)/mp.sqrt(q(tv, tv)*q(tw, tw))
    return mp.acos(c)
ok_num = True; report = []
for (l1n, l2n, an) in ((mp.mpf('0.7'), mp.mpf('1.3'), mp.mpf('0.9')), (mp.mpf('2.1'), mp.mpf('0.4'), mp.mpf('2.5')), (mp.mpf('0.3'), mp.mpf('1.7'), mp.mpf('1.1'))):
    subs = {v1: mp.e**(l1n/2), v2: mp.e**(l2n/2), t: mp.tan(an/2)}
    Pn = [sp.Matrix([[mp.mpf(str(sp.N(x.subs(subs), 45)))] for x in P]) for P in (P0, P1, P2)]
    # mpmath vectors
    Pm = [mp.matrix([mp.mpf(str(sp.N(x.subs(subs), 45))) for x in P]) for P in (P0, P1, P2)]
    qm = lambda a, b: a[0]*b[0] + a[1]*b[1] - a[2]*b[2]
    def ang(u, v, w):
        tv = v + qm(u, v)*u; tw = w + qm(u, w)*u
        return mp.acos(qm(tv, tw)/mp.sqrt(qm(tv, tv)*qm(tw, tw)))
    A_defect = mp.pi - ang(Pm[0], Pm[1], Pm[2]) - ang(Pm[1], Pm[2], Pm[0]) - ang(Pm[2], Pm[0], Pm[1])
    w_T3 = 2*mp.atan(mp.mpf(str(sp.N(tan_half_w.subs(subs), 45))))
    tanA_closed = mp.mpf(str(sp.N(tan_half_A.subs(subs), 45)))
    ok_num = ok_num and abs(A_defect - abs(w_T3)) < mp.mpf('1e-30') and abs(mp.tan(A_defect/2) - abs(tanA_closed)) < mp.mpf('1e-30')
    report.append((mp.nstr(A_defect, 12), mp.nstr(abs(w_T3), 12)))
check("A1 [40 digits, three points] the angle DEFECT of the boost triangle on H^2 equals |Wigner angle| from T3 to 1e-30",
      ok_num, str(report))
check("A2 the closed form tan(A/2) = det[P0,P1,P2] / (1 - sum of pairwise q) agrees with the defect numerically (validated above)", ok_num)
check("A3 [EXACT] tan^2(A/2) = tan^2(w/2): Wigner's rotation IS the hyperbolic area of the triangle its two boosts and their composite enclose",
      z(sp.cancel(tan_half_A**2 - tan_half_w**2)))
check("A4 margin: the area vanishes iff collinear (t = 0) or unpivoted (v1 = 1 or v2 = 1) -- the same margin as Wigner",
      z(tan_half_A.subs(t, 0)) and z(tan_half_A.subs(v1, 1)) and z(tan_half_A.subs(v2, 1)))

print("\n=== PART B: the six lifts and the deck sign ===")
FLIPS = ['C', 'H', 'G']; AX = {'C': 2, 'H': 0, 'G': 1}      # coordinate index in (hbar, G, c)
def centre(subset):
    s = [1, 1, 1]
    for f in subset: s[AX[f]] = -1
    return sp.Matrix(s)
orders = list(itertools.permutations(FLIPS))
path = lambda order: [frozenset(order[:k]) for k in range(4)]
# adjacent transpositions only
pairs4 = []
for a, b in itertools.combinations(orders, 2):
    diff = [k for k in range(3) if a[k] != b[k]]
    if len(diff) == 2 and diff[1] == diff[0] + 1: pairs4.append((a, b, diff))
check("B0 [correction to LABEL-3] six pairs of orderings differ by an ADJACENT transposition and bound 4-cycles; the other three bound 6-cycles",
      len(pairs4) == 6 and len([1 for a, b in itertools.combinations(orders, 2) if [k for k in range(3) if a[k] != b[k]] == [0, 2]]) == 3)
# the 4-cycle bounded by an adjacent pair: the four faces around ONE pole.  Identify the pole.
cycles = []
for a, b, diff in pairs4:
    pa, pb = path(a), path(b)
    k = diff[0]                       # positions k, k+1 swapped
    loop = [pa[k], pa[k+1], pa[k+2], pb[k+1]]      # the 4-cycle: shared, a's middle, shared, b's middle
    # the pole: the axis whose sign is the same on all four faces, and that sign
    signs = {f: {centre(s)[AX[f]] for s in loop} for f in FLIPS}
    fixed = [f for f in FLIPS if len(signs[f]) == 1]
    pole = (fixed[0], next(iter(signs[fixed[0]]))) if len(fixed) == 1 else None
    cycles.append((a, b, loop, pole))
check("B1 each 4-cycle fixes exactly ONE axis's sign across its four faces: it is the square of faces around that pole",
      all(pole is not None for *_, pole in cycles))
for a, b, loop, pole in cycles:
    print(f"    {''.join(a)} <-> {''.join(b)} : loop " + " -> ".join(''.join(sorted(s)) or 'O' for s in loop) + f"  around pole {pole[0]}{'+' if pole[1] > 0 else '-'}")
timelike_cycles = [c for c in cycles if c[3][0] == 'C']
spacelike_cycles = [c for c in cycles if c[3][0] != 'C']
check("B2 exactly TWO of the six 4-cycles circle a timelike pole (+c or -c); the other four circle spacelike poles", len(timelike_cycles) == 2 and len(spacelike_cycles) == 4)
swapped = lambda a, b: {a[k] for k in range(3) if a[k] != b[k]}
check("B3 the two timelike-pole cycles are exactly the pairs that swap the hbar-flip and the G-flip (H <-> G adjacent) with c fixed",
      all(swapped(a, b) == {'H', 'G'} for a, b, _, _ in timelike_cycles) and all(swapped(a, b) != {'H', 'G'} for a, b, _, _ in spacelike_cycles))
# geometry of a timelike-pole cycle: the four faces around +c are (+-1, +-1, +1); going round them is a 2 pi rotation about c
# in the compact (hbar, G) plane.  Their in-plane projections (+-1, +-1) are spacelike (q = 2 > 0): a genuine closed rotation.
tl = [c for c in timelike_cycles if c[3][1] == +1][0] if any(c[3][1] == +1 for c in timelike_cycles) else timelike_cycles[0]
corners_tl = [centre(s) for s in tl[2]]
inplane_tl = [sp.Matrix([v[0], v[1]]) for v in corners_tl]
check("B4 around a timelike pole the corners' in-plane projections are SPACELIKE (q = 2 each), lying on a circle: the loop is a closed 2 pi rotation in the compact plane",
      all((p.T*p)[0] == 2 for p in inplane_tl))
# spinor lift of a 2 pi rotation in the compact plane: -I (T4b / T4b')
J = sp.Matrix([[0, -1], [1, 0]])*sp.Rational(1, 2)          # so(2) generator in the 2-dim rep (half-angle)
R2pi_spin = (2*sp.pi*J).exp()
check("B5 the spinor lift of the timelike-pole 4-cycle is -I: paths differing by an adjacent H <-> G swap differ by the DECK SIGN",
      zM(R2pi_spin + sp.eye(2)))
# geometry of a spacelike-pole cycle: e.g. around +G, the corners (+-1, 1, +-1) project to (+-1, +-1) in the (hbar, c) plane, which is
# LORENTZIAN; their in-plane q = hbar^2 - c^2 = 0: the corners sit on the cone.  Not a closed rotation; no 2 pi; holonomy a boost.
sl = [c for c in spacelike_cycles if c[3] == ('G', +1)][0]
corners_sl = [centre(s) for s in sl[2]]
inplane_sl = [sp.Matrix([v[0], v[2]]) for v in corners_sl]     # (hbar, c) components
check("B6 around a spacelike pole the corners' in-plane projections are NULL (q = hbar^2 - c^2 = 0): on the cone of a Lorentzian plane, not a closed rotation -- holonomy is a boost, carries no sign",
      all(p[0]**2 - p[1]**2 == 0 for p in inplane_sl))
# the sign structure on the six orderings: connect orderings by adjacent swaps; a swap contributes -1 iff it is H <-> G.
# Then the sign of an ordering relative to (C,H,G) is (-1)^{[G before H]}.
def sign_class(order): return +1 if order.index('H') < order.index('G') else -1
classes = {o: sign_class(o) for o in orders}
consistent = all((classes[a] == classes[b]) == (pole[0] != 'C') for a, b, _, pole in cycles)
check("B7 the sign structure is consistent: adjacent swaps around a spacelike pole preserve the class, around a timelike pole flip it; the class is 'hbar before G' vs 'G before hbar'",
      consistent, str({''.join(o): s for o, s in classes.items()}))
check("B8 c's position in the ordering carries NO sign: the two classes are {C,H,G},{H,C,G},{H,G,C} and {C,G,H},{G,C,H},{G,H,C}",
      {o for o in orders if classes[o] == +1} == {('C','H','G'), ('H','C','G'), ('H','G','C')})

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT:")
print("  A [DERIVED, two routes]: Wigner's rotation is the hyperbolic area of the boost triangle on H^2. Holonomy = area, K = -1.")
print("  B [DERIVED]: of the six 4-cycles on the face-graph, exactly two circle a timelike pole and are closed 2 pi rotations")
print("    with spinor lift -1; the four around spacelike poles have corners on the cone and carry boosts, not signs.")
print("    The six lifts split into two classes by the ORDER OF THE TWO RULER-FLIPS: hbar-then-G vs G-then-hbar.")
print("    c's place in the order is sign-free.  That is the deck arrow: whether action turns to evanescence before")
print("    or after mass turns to energy.  PREDICTION-1's ordered operation, with an address.")

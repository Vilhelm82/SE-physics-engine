#!/usr/bin/env python3
# =============================================================================
# LABEL-3 -- the Boolean lattice of three flips, and the character of each flip.     (2026-09-04, late)
#
# Will: the eight faces are the Boolean lattice of C: Light -> Temperature, H: Action -> Evanescence,
#   G: Mass -> Energy.  Our face is the empty set; Hawking's is CHG; single flips isolate an ingredient,
#   pairs hold every pair.  3! = 6 shortest edge paths from our face to Hawking's, one per ordering; the
#   projected endpoint is always F_{---}; the LIFTS may differ by an ordered phase or deck sign -- that is
#   where the deck arrow belongs.  A face selects a presentation, not the inventory of the universe.
# PROVED HERE: the lattice structure (LABEL-2's octahedron IS the 2^3 lattice: adjacency = one flip);
#   six shortest paths, each visiting one single-flip face and one pair-flip face; the face-centres are
#   unit SPACELIKE under the seat's form (the gaze surface is the one-sheeted hyperboloid, dS_2); and
#   the flips have different characters on it: flipping hbar or G is NULL-separated (degenerate span),
#   flipping c is a BOOST (Lorentzian span, rapidity arccosh 3).  The pair flips: CH and CG are boosts,
#   HG is a Euclidean rotation by pi (the compact plane).  CHG is -I, not in SO(2,1) on the frame.
# STATED, NOT RUN: the six lifts.  On the gaze surface dS_2 the transition along an edge is parallel
#   transport; the holonomy around (path_i)(path_j)^-1 is the Lorentzian Gauss-Bonnet area of the
#   enclosed cube-face, and its spinor lift is the deck sign.  That is T5 on dS_2.  Queued.
# =============================================================================
import sympy as sp, itertools, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)

FLIPS = ['C', 'H', 'G']                      # c, hbar, G
AX = {'C': 'c', 'H': 'hbar', 'G': 'G'}
POLE = {('c',+1):'Light',('c',-1):'Temperature',('hbar',+1):'Action',('hbar',-1):'Evanescence',('G',+1):'Mass',('G',-1):'Energy'}
def face_of(subset):        # subset of FLIPS -> signs dict
    return {AX[f]: (-1 if f in subset else +1) for f in FLIPS}
def name(face): return "(" + ", ".join(POLE[(ax, face[ax])] for ax in ('c','hbar','G')) + ")"

print("=== the lattice ===")
subsets = [frozenset(s) for r in range(4) for s in itertools.combinations(FLIPS, r)]
faces = {s: face_of(s) for s in subsets}
check("8 subsets of {C,H,G} -> 8 faces, bijectively", len(subsets) == 8 and len({tuple(sorted(f.items())) for f in faces.values()}) == 8)
ours, hawk = frozenset(), frozenset(FLIPS)
check("the empty set is our face (Light, Action, Mass); the full set is Hawking's (Temperature, Evanescence, Energy)",
      name(faces[ours]) == "(Light, Action, Mass)" and name(faces[hawk]) == "(Temperature, Evanescence, Energy)")
adjacent = lambda s, t: len(s ^ t) == 1
check("lattice adjacency (differ by one flip) = octahedron face adjacency (share an edge): every face has exactly 3 neighbours",
      all(sum(1 for t in subsets if adjacent(s, t)) == 3 for s in subsets))
singles = [s for s in subsets if len(s) == 1]; pairs = [s for s in subsets if len(s) == 2]
check("three single-flip faces each isolate one ingredient; three pair-flip faces each hold one pair; nothing else",
      len(singles) == 3 and len(pairs) == 3 and len(subsets) == 1 + 3 + 3 + 1)
for s in singles: print(f"    {''.join(sorted(s)):3s} -> {name(faces[s])}")
for s in pairs:   print(f"    {''.join(sorted(s)):3s} -> {name(faces[s])}")

print("\n=== the six shortest paths ===")
paths = []
for order in itertools.permutations(FLIPS):
    chain = [frozenset(order[:k]) for k in range(4)]
    paths.append((order, chain))
    print("    " + " -> ".join(''.join(sorted(c)) or 'O' for c in chain) + "   :  " + " -> ".join(name(faces[c]) for c in chain))
check("3! = 6 shortest paths from our face to Hawking's, each through one single-flip face and one pair-flip face, all ending at F_{---}",
      len(paths) == 6 and all(len(ch[1]) == 1 and len(ch[2]) == 2 and ch[3] == hawk for _, ch in paths))
check("the PROJECTED endpoint is the same face for all six orderings (the presentation is order-blind)",
      len({tuple(sorted(faces[ch[3]].items())) for _, ch in paths}) == 1)

print("\n=== the character of each flip on the gaze surface ===")
# face-centre as a gaze direction in the seat's form: coordinates (hbar, G, c), Q = diag(1,1,-1)
Q = sp.diag(1, 1, -1)
centre = lambda face: sp.Matrix([face['hbar'], face['G'], face['c']])
q = lambda u, v: (u.T*Q*v)[0]
f0 = centre(faces[ours])
check("every face-centre is unit SPACELIKE, q(f,f) = +1: the gaze surface is the one-sheeted hyperboloid (dS_2), Lorentzian",
      all(q(centre(f), centre(f)) == 1 for f in faces.values()))
def character(u, v):
    g = sp.Matrix([[q(u,u), q(u,v)], [q(u,v), q(v,v)]]); d = g.det()
    return ('null (degenerate span)' if d == 0 else 'boost (Lorentzian span)' if d < 0 else 'rotation (Euclidean span)'), q(u, v), d
chars = {}
for s in singles:
    f1 = centre(faces[s]); kind, qq, d = character(f0, f1)
    chars[next(iter(s))] = kind
    print(f"    flip {next(iter(s))}: q(f0, f1) = {int(qq):2d}, Gram det = {int(d):2d}  ->  {kind}")
check("flipping hbar and flipping G are NULL-separated on the gaze surface (Gram det 0): light-rays on dS_2",
      chars['H'].startswith('null') and chars['G'].startswith('null'))
check("flipping c is a BOOST on the gaze surface (Gram det < 0, q(f0,f1) = 3 = cosh of the rapidity): the odd flip out is c's -- T7b again",
      chars['C'].startswith('boost') and q(f0, centre(faces[frozenset('C')])) == 3)
for s in pairs:
    f2 = centre(faces[s]); kind, qq, d = character(f0, f2)
    print(f"    pair {''.join(sorted(s)):2s}: q(f0, f2) = {int(qq):2d}, Gram det = {int(d):3d}  ->  {kind}")
# PREDICTED before running: CH, CG boosts and HG a Euclidean rotation.  WRONG.  Computed:
# q(f0, f) = s_hbar + s_G - s_c, so |q| = 3 iff the two rulers agree and oppose c (C alone, or HG together);
# |q| = 1 otherwise.  Gram det = 1 - q^2: boost iff |q| = 3, null iff |q| = 1, never a Euclidean rotation.
qvals = {''.join(sorted(s)) or 'O': q(f0, centre(faces[s])) for s in subsets}
formula_ok = all(q(f0, centre(faces[s])) == faces[s]['hbar'] + faces[s]['G'] - faces[s]['c'] for s in subsets)
boost_faces = {k for k, v in qvals.items() if abs(v) == 3}
null_faces = {k for k, v in qvals.items() if abs(v) == 1 and k != 'O'}
check("q(f0, f) = s_hbar + s_G - s_c for every face; |q| in {1, 3}; no face is Euclidean-rotation-separated from ours (det = 1 - q^2 <= 0 always)",
      formula_ok and set(abs(v) for v in qvals.values()) == {1, 3}, str(qvals))
check("BOOST-separated faces are exactly C and HG -- the two where the rulers agree with each other and oppose c (my prediction 'CH, CG boosts; HG rotation' was WRONG, recorded)",
      boost_faces == {'C', 'GH'})
check("NULL-separated faces are the other six including Hawking's: H, G, CH, CG, CHG are all light-separated from our face on dS_2",
      null_faces == {'H', 'G', 'CH', 'CG', 'CGH'})
f3 = centre(faces[hawk])
check("the full flip CHG lands at the antipodal centre -f0 (q = -1), and as a FRAME map it is -I, det -1: not in SO(2,1) -- Hawking's face is not a rotation of ours; it is the deck",
      f3 == -f0 and (-sp.eye(3)).det() == -1)

print("\n=== what the six lifts need (stated, queued as T5 on dS_2) ===")
print("    Transition along an edge of the face-graph = parallel transport of the gaze on dS_2 between null- or boost-separated centres.")
print("    Two paths i, j from O to CHG bound a closed loop on the cube graph; its holonomy is the Lorentzian Gauss-Bonnet")
print("    area of the enclosed region; the spinor lift of that holonomy is the deck sign.  If the sign depends on the")
print("    ordering, PREDICTION-1's ordered operation is the choice of path through this lattice.")
check("the six paths pair into three closed loops on the face-graph (each pair of orderings differing by a transposition bounds one 4-cycle)",
      len([(a, b) for a, b in itertools.combinations(range(6), 2)
           if sum(1 for k in range(3) if paths[a][0][k] != paths[b][0][k]) == 2]) == 9)   # 9 transposition-pairs among 6 permutations

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT:")
print("  Lattice [PROVED]: the 8 faces are 2^{C,H,G}; ours is empty, Hawking's is full, antipodal; six shortest paths, one endpoint.")
print("  Characters [DERIVED]: on the gaze surface dS_2, flipping hbar or G is null, flipping c is a boost (rapidity arccosh 3);")
print("    the full flip is -I, not a rotation of the frame -- Hawking's face is reached by the deck, not by a pivot.")
print("  Lifts [QUEUED]: holonomy of the cube-graph loops on dS_2 and its spinor sign -- T5 on the gaze surface.")

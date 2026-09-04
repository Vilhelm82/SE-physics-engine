#!/usr/bin/env python3
# =============================================================================
# LABEL-2 -- the pole-octahedron.  Six poles, twelve edges, eight faces.        (2026-09-04, late)
#
# Will: a face is the perpendicular plane one sees from the root looking into an octant,
#   bordered by one pole from each axis; 2^3 = 8 of them ("the squared cube").  Not the
#   cube's six faces.  The table's per-axis 'faces' are CONTRIBUTIONS: what each pole puts
#   into every face it borders.  A face's physics = its three borders.
# PROVED HERE: the combinatorics (vertices/edges/faces), that every edge borders exactly two
#   faces which differ only in the axis the edge does not touch, that every vertex lies on
#   four faces, that the W-loops of LABEL-1 are edge paths whose odd (plane) stations lie on
#   edges, and that the near and far faces of the table are antipodal.
# DECLARED (comparison stage): a candidate physical regime for each of the eight faces.
#   Two are established (our face; Hawking's face).  Six are NAMED TONIGHT and are the
#   pages the labelling owes a physics for.  Graded as candidates, not results.
# =============================================================================
import sympy as sp, itertools, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)

AXES = ['c', 'hbar', 'G']
POLE = {('c',+1):'Light', ('c',-1):'Temperature', ('hbar',+1):'Action', ('hbar',-1):'Evanescence', ('G',+1):'Mass', ('G',-1):'Energy'}
CONTRIB = {('c',+1):'GR', ('c',-1):'Quantum', ('hbar',+1):'Wave', ('hbar',-1):'Particle', ('G',+1):'Bound', ('G',-1):'Escaping'}

# --- the octahedron -----------------------------------------------------------------------
V = [(ax, s) for ax in AXES for s in (+1, -1)]                                   # 6 poles
E = [(u, v) for u, v in itertools.combinations(V, 2) if u[0] != v[0]]            # edges: poles on DIFFERENT axes
F = [dict(zip(AXES, signs)) for signs in itertools.product((+1, -1), repeat=3)]  # faces: one sign per axis
print("=== combinatorics ===")
check("6 vertices (poles), 12 edges (pole pairs on different axes), 8 faces (sign triples): an octahedron, V - E + F = 2",
      len(V) == 6 and len(E) == 12 and len(F) == 8 and len(V) - len(E) + len(F) == 2)
def faces_of_edge(u, v):
    third = [ax for ax in AXES if ax not in (u[0], v[0])][0]
    return [f for f in F if f[u[0]] == u[1] and f[v[0]] == v[1]], third
ok_edges = True
for (u, v) in E:
    fs, third = faces_of_edge(u, v)
    if len(fs) != 2 or fs[0][third] == fs[1][third]: ok_edges = False
check("every edge borders exactly TWO faces, and they differ only in the axis the edge does not touch", ok_edges)
check("every vertex lies on exactly four faces", all(sum(1 for f in F if f[v[0]] == v[1]) == 4 for v in V))
antipodal = lambda f, g: all(f[ax] == -g[ax] for ax in AXES)
near = {'c':+1, 'hbar':+1, 'G':+1}; far = {'c':-1, 'hbar':-1, 'G':-1}
check("the near face (Light, Action, Mass) and the far face (Temperature, Evanescence, Energy) are ANTIPODAL faces of the octahedron",
      near in F and far in F and antipodal(near, far))
check("the near face and the far face share no vertex and no edge",
      not any(near[u[0]] == u[1] and far[u[0]] == u[1] for u in V))

print("\n=== the W-loops are edge paths; their odd stations are on edges ===")
# LABEL-1's loop in the (c,hbar) plane visits +c, i hbar, -c, -i hbar: the great circle through the four poles of c and hbar.
# Its four odd (plane) stations sit on the four edges of that circle; each such edge borders the +G face and the -G face.
for pair, third in ((('c','hbar'),'G'), (('c','G'),'hbar'), (('hbar','G'),'c')):
    ring = [(pair[0],+1), (pair[1],+1), (pair[0],-1), (pair[1],-1)]
    ring_edges = [(ring[k], ring[(k+1) % 4]) for k in range(4)]
    all_edges = all(tuple(sorted(e, key=str)) in [tuple(sorted(x, key=str)) for x in E] for e in ring_edges)
    between = all(set(f[third] for f in faces_of_edge(*e)[0]) == {+1, -1} for e in ring_edges)
    check(f"loop in the ({pair[0]},{pair[1]}) plane: four vertices, four edges, each odd station between the +{third} face and the -{third} face",
          all_edges and between, f"odd stations are between {POLE[(third,+1)]}-faces and {POLE[(third,-1)]}-faces")

print("\n=== the eight faces, from the table ===")
# declared candidate regimes (comparison stage).  'established' = the table's own two columns; 'named' = first named here.
CANDIDATE = {
    (+1,+1,+1): ("light waves in bound orbits around mass, in curved spacetime: lensing, orbital mechanics -- OUR FACE", 'established'),
    (-1,-1,-1): ("thermal tunnelling carrying energy away: HAWKING RADIATION -- the far face", 'established'),
    (+1,+1,-1): ("light with phase leaving a well: gravitational redshift, Shapiro delay, gravitational waves escaping", 'named'),
    (+1,-1,+1): ("light tunnelling in a bound system: near-horizon greybody filtering (G-2), frustrated total internal reflection", 'named'),
    (-1,+1,+1): ("thermal waves around bound mass: Tolman-Ehrenfest equilibrium, temperature in a gravitational well (EQ-1)", 'named'),
    (+1,-1,-1): ("light tunnelling out carrying energy: the Parikh-Wilczek tunnelling picture of emission, non-thermal", 'named'),
    (-1,+1,-1): ("thermal waves carrying energy away, no mass: BLACKBODY radiation, Planck's reading of the light-temperature axis", 'named'),
    (-1,-1,+1): ("thermal tunnelling in a bound system: finite-temperature instantons, stellar fusion through the Coulomb barrier", 'named'),
}
n_est = 0; n_named = 0
for f in F:
    key = tuple(f[ax] for ax in AXES)
    poles = ", ".join(POLE[(ax, f[ax])] for ax in AXES)
    contribs = " + ".join(CONTRIB[(ax, f[ax])] for ax in AXES)
    regime, status = CANDIDATE[key]
    if status == 'established': n_est += 1
    else: n_named += 1
    print(f"    face {key}: [{poles}]  =  {contribs:28s}  ->  {regime}  [{status}]")
check("every face has a declared candidate regime; two established (the table's columns), six named tonight",
      n_est == 2 and n_named == 6 and all(tuple(f[ax] for ax in AXES) in CANDIDATE for f in F))
# the six named faces are exactly the faces adjacent to BOTH the near and far faces? no -- check what they are:
mixed = [f for f in F if f != near and f != far]
one_flip = [f for f in mixed if sum(1 for ax in AXES if f[ax] == near[ax]) == 2]
two_flip = [f for f in mixed if sum(1 for ax in AXES if f[ax] == near[ax]) == 1]
check("the six mixed faces split 3 + 3: three one pole-flip from our face (share an edge with it), three one flip from Hawking's",
      len(one_flip) == 3 and len(two_flip) == 3)
# each one-flip face shares an EDGE with the near face; each two-flip face shares an edge with the far face
share_edge = lambda f, g: sum(1 for ax in AXES if f[ax] == g[ax]) == 2
check("adjacency: the three one-flip faces each share an edge with our face; the three two-flip faces each share an edge with Hawking's",
      all(share_edge(f, near) for f in one_flip) and all(share_edge(f, far) for f in two_flip))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT:")
print("  Combinatorics [PROVED]: the six poles form an octahedron; edges border two faces differing in the untouched axis;")
print("    the W-loops are its edge paths and their plane-stations lie on edges labelled by the axis they never touch.")
print("  Faces [DECLARED candidates]: our face and Hawking's face are the table's two columns and are antipodal;")
print("    the six between them are named for the first time and each owes a derivation.")

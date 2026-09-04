#!/usr/bin/env python3
# =============================================================================
# LABEL-1 -- the Wick loop run on every direction from every pole.          (2026-09-04, late)
#
# INPUTS: the seat's (2,1) form from T7a with c = e3 negative, rulers e1 (hbar), e2 (G);
#   the Wick face exp(i theta K) from T4 for the two hyperbolic planes; the real rotation
#   exp(theta J) for the compact plane.  Will's labelling table (CONJECTURE tier, 2026-09-04):
#       axis   faces               poles
#       c      GR <-> Quantum      Light <-> Temperature
#       hbar   Wave <-> Particle   Action <-> Evanescence
#       G      Bound <-> Escaping  Mass <-> Energy
# COMPARISON-STAGE TABLE (declared, not derived): known physical relations between quantities.
#   The test asks whether every step of every loop lands on a known relation, and names the
#   shared-station consistency conditions the labelling implies.
# WHAT IS PROVED HERE: the geometry of the loops (stations land on +-axes, quarter-turn
#   structure, the two Wick loops share c's poles, each pair of loops shares one axis).
# WHAT IS CHECKED AGAINST A DECLARED TABLE: the physics at each step.  That part is a
#   labelling test, not a derivation, and is graded as such.
# =============================================================================
import sympy as sp, itertools, time
T0 = time.time(); CH = []
def check(t, ok, n=""):
    CH.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}][{time.time()-T0:5.1f}s] {t}" + (f" -- {n}" if n else ""), flush=True)
def zM(M): return all(sp.simplify(e) == 0 for e in M)

I = sp.I; th = sp.Symbol('theta', real=True)
e = {'c': sp.Matrix([0,0,1]), 'hbar': sp.Matrix([1,0,0]), 'G': sp.Matrix([0,1,0])}
Q = sp.diag(1, 1, -1)

# --- the labelling (Will's table) ----------------------------------------------------------
POLE = {('c', +1): 'Light',  ('c', -1): 'Temperature',
        ('hbar', +1): 'Action', ('hbar', -1): 'Evanescence',
        ('G', +1): 'Mass',   ('G', -1): 'Energy'}
FACE = {('c', +1): 'GR',    ('c', -1): 'Quantum',
        ('hbar', +1): 'Wave', ('hbar', -1): 'Particle',
        ('G', +1): 'Bound', ('G', -1): 'Escaping'}

# --- the declared physics table (comparison stage) -----------------------------------------
# key: frozenset of two pole labels -> (relation, kind).  kind: 'wick' = a 90-degree imaginary
# step (real <-> imaginary counterpart); 'anti' = 180 degrees (the two ends of one line);
# 'real' = a 90-degree real rotation between rulers (dimensional conversion in the seat's space).
KNOWN = {
    frozenset({'Light','Action'}):        ('photon phase = S/hbar; e^{iS/hbar}', 'wick'),
    frozenset({'Action','Temperature'}):  ('Euclidean action = beta E; e^{iS} -> e^{-S_E}: QM <-> stat mech', 'wick'),
    frozenset({'Temperature','Evanescence'}): ('thermal decay e^{-beta E}; finite-T instantons', 'wick'),
    frozenset({'Evanescence','Light'}):   ('tunnelling out as radiation (Parikh-Wilczek Hawking)', 'wick'),
    frozenset({'Light','Mass'}):          ('E = m c^2 with E the photon energy; lensing', 'wick'),
    frozenset({'Mass','Temperature'}):    ('T_H = hbar c^3/(8 pi G M k_B)', 'wick'),
    frozenset({'Temperature','Energy'}):  ('E = k_B T', 'wick'),
    frozenset({'Energy','Light'}):        ('E = h nu', 'wick'),
    frozenset({'Action','Mass'}):         ('S = -m c^2 integral d tau', 'real'),
    frozenset({'Mass','Evanescence'}):    ('Compton / Yukawa decay length hbar/(m c)', 'real'),
    frozenset({'Evanescence','Energy'}):  ('tunnelling rate ~ e^{-2 kappa d}, kappa = sqrt(2m(V-E))/hbar', 'real'),
    frozenset({'Energy','Action'}):       ('S = integral E dt', 'real'),
    # antipodes
    frozenset({'Light','Temperature'}):   ('Unruh: the light-seat under acceleration reads a temperature', 'anti'),
    frozenset({'Action','Evanescence'}):  ('S -> i S_E: oscillatory <-> decaying', 'anti'),
    frozenset({'Mass','Energy'}):         ('dE = -c^2 dM: what the bound pole loses the escaping pole carries', 'anti'),
}

# --- the three loop operators ------------------------------------------------------------------
def boost_gen(u, v):      # symmetric generator mixing the negative line u with a positive line v
    return (v*u.T + u*v.T)*sp.Matrix([[1]])[0,0]
K_h = sp.Matrix([[0,0,1],[0,0,0],[1,0,0]])   # (c, hbar) plane: hyperbolic
K_G = sp.Matrix([[0,0,0],[0,0,1],[0,1,0]])   # (c, G) plane: hyperbolic
J   = sp.Matrix([[0,-1,0],[1,0,0],[0,0,0]])  # (hbar, G) plane: compact
for name, M, kind in (('K_h', K_h, 'boost'), ('K_G', K_G, 'boost'), ('J', J, 'rotation')):
    check(f"{name} is an so(2,1) generator ({kind})", zM(M.T*Q + Q*M))
E_h = sp.eye(3) + (sp.cos(th) - 1)*K_h**2 + I*sp.sin(th)*K_h      # exp(i theta K_h), T4a
E_G = sp.eye(3) + (sp.cos(th) - 1)*K_G**2 + I*sp.sin(th)*K_G      # exp(i theta K_G)
R_J = sp.eye(3) + sp.sin(th)*J + (1 - sp.cos(th))*J**2            # exp(theta J), real rotation
check("E_h = exp(i theta K_h): E' = i K_h E, E(0) = I", zM(E_h.diff(th) - I*K_h*E_h) and zM(E_h.subs(th,0) - sp.eye(3)))
check("R_J = exp(theta J): R' = J R, R(0) = I", zM(R_J.diff(th) - J*R_J) and zM(R_J.subs(th,0) - sp.eye(3)))

LOOPS = {
    'A (c,hbar) WICK': (E_h, ['c', 'hbar'], 'wick'),
    'B (c,G)    WICK': (E_G, ['c', 'G'],    'wick'),
    'C (hbar,G) REAL': (R_J, ['hbar', 'G'], 'real'),
}

def station(vec):
    """Identify a station vector as (axis, sign, factor) where vec = factor * sign * e[axis], factor in {1, i}."""
    for ax, basis in e.items():
        for sgn in (+1, -1):
            for fac in (1, I):
                if zM(vec - fac*sgn*basis): return ax, sgn, fac
    return None

print("\n=== geometry: every loop, every start, both directions ===")
all_ok = True
transcripts = []
for lname, (Op, axes, kind) in LOOPS.items():
    starts = [(ax, sgn) for ax in axes for sgn in (+1, -1)]
    for (ax0, s0) in starts:
        for direction in (+1, -1):
            v0 = s0*e[ax0]
            seq = [(ax0, s0)]
            for k in (1, 2, 3, 4):
                vk = sp.simplify(Op.subs(th, direction*k*sp.pi/2)*v0)
                st = station(vk)
                if st is None: all_ok = False; seq.append(None); continue
                seq.append((st[0], st[1]))
            closed = seq[-1] == seq[0]
            axes_visited = {s[0] for s in seq[:4] if s}
            two_each = all(sum(1 for s in seq[:4] if s and s[0] == a) == 2 for a in axes)
            ok = closed and axes_visited == set(axes) and two_each
            all_ok = all_ok and ok
            labels = [POLE[s] for s in seq[:4]]
            transcripts.append((lname, ax0, s0, direction, kind, seq[:4], labels))
check("every loop from every pole in both directions: four quarter-turn stations, two per axis, alternating axes, closing at 2 pi",
      all_ok, f"{len(transcripts)} walks")
# antipode structure: station 2 is minus the start; the two intermediate stations are +-(i or 1) * the other axis
anti_ok = all(seq[2] == (seq[0][0], -seq[0][1]) for (_,_,_,_,_,seq,_) in transcripts)
check("station 2 of every walk is the ANTIPOLE of the start (i^2 = -1): the loop passes through the far pole of the starting axis", anti_ok)
# the intermediate stations in the WICK loops carry the factor i; in the REAL loop they do not
def factors(lname, Op, ax0, s0):
    v1 = sp.simplify(Op.subs(th, sp.pi/2)*(s0*e[ax0])); st = station(v1); return st[2]
wick_i = all(factors(l, LOOPS[l][0], a, s) == I for l in ('A (c,hbar) WICK','B (c,G)    WICK') for a in LOOPS[l][1] for s in (1,-1))
real_1 = all(factors('C (hbar,G) REAL', R_J, a, s) == 1 for a in ('hbar','G') for s in (1,-1))
check("in the two WICK loops the intermediate stations are on the complex plane (factor i); in the REAL loop they are real", wick_i and real_1)

print("\n=== the physics at every step (declared comparison table) ===")
missing = []; kinds_ok = True
seen = set()
for (lname, ax0, s0, direction, kind, seq, labels) in transcripts:
    for k in range(4):
        L1, L2 = labels[k], labels[(k+1) % 4]
        key = frozenset({L1, L2})
        if key in seen: continue
        seen.add(key)
        if key not in KNOWN:
            missing.append((lname, L1, L2)); continue
        rel, rkind = KNOWN[key]
        # adjacent steps in a Wick loop should be 'wick', in the real loop 'real'
        expected = 'wick' if kind == 'wick' else 'real'
        if rkind != expected: kinds_ok = False
        print(f"    {lname}: {L1:12s} -> {L2:12s}  [{rkind}]  {rel}")
check("every adjacent step of every loop has a declared known relation", len(missing) == 0, str(missing))
check("every adjacent relation's kind matches its loop (Wick loops: imaginary counterparts; real loop: dimensional conversions)", kinds_ok)
# antipodes
anti_missing = []
for ax in ('c','hbar','G'):
    key = frozenset({POLE[(ax,+1)], POLE[(ax,-1)]})
    if key not in KNOWN or KNOWN[key][1] != 'anti': anti_missing.append(ax)
    else: print(f"    antipode on {ax:5s}: {POLE[(ax,+1)]:12s} <-> {POLE[(ax,-1)]:12s}  {KNOWN[key][0]}")
check("every axis's two poles have a declared 'antipode' relation", len(anti_missing) == 0)

print("\n=== shared stations: the consistency conditions the labelling IMPLIES ===")
# Loops A and B share the c line: {Light, Temperature}.  Loops A and C share hbar.  Loops B and C share G.
shared = {('A','B'): ('c',   'Temperature reached via hbar (Unruh: action -> Euclidean action -> Boltzmann) must equal Temperature reached via G (Hawking: mass -> horizon).  DERIVED: T4 (same 2 pi, same kappa).'),
          ('A','C'): ('hbar','Action/Evanescence reached from c (Wick) must equal Action/Evanescence reached from G (real rotation).  OWED.'),
          ('B','C'): ('G',   'Mass/Energy reached from c (Wick: E = mc^2 at the light pole) must equal Mass/Energy reached from hbar (real rotation: S = -mc^2 tau).  OWED.')}
for (l1, l2), (ax, cond) in shared.items():
    print(f"    loops {l1} & {l2} share axis {ax}: {cond}")
pair_axes = {}
names = list(LOOPS.keys())
for i in range(3):
    for j in range(i+1, 3):
        pair_axes[(names[i][0], names[j][0])] = set(LOOPS[names[i]][1]) & set(LOOPS[names[j]][1])
check("the three loops pairwise share exactly one axis each (computed from the loop definitions): three consistency conditions, one per axis",
      all(len(v) == 1 for v in pair_axes.values()) and set(map(lambda s: next(iter(s)), pair_axes.values())) == {'c','hbar','G'},
      str({k: next(iter(v)) for k, v in pair_axes.items()}))
# geometric content of (A,B): from Light, both Wick loops reach the SAME vector at station 2 (-c).  The numeric
# equality of the temperature reached each way is T4 (prim_t4_hawking_period.py, 30/30), cited not re-run.
stA = sp.simplify(E_h.subs(th, sp.pi)*e['c']); stB = sp.simplify(E_G.subs(th, sp.pi)*e['c'])
check("condition (A,B), geometric half: from Light both Wick loops land on the same station vector -c at 2 quarter-turns; the numeric half (Unruh T = Hawking T) is T4",
      zM(stA - stB) and zM(stA + e['c']))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT:")
print("  Geometry [PROVED]: three loops, four stations each, two per axis, antipole at i^2, Wick loops on the complex plane, real loop real.")
print("  Labelling [CONJECTURE, checked against a DECLARED table]: every step and every antipode names a known relation, with the right kind.")
print("  Consistency [one DERIVED, two OWED]: the shared-station condition on c is T4; on hbar and on G it has not been run.")

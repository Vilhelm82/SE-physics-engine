#!/usr/bin/env python3
# CENSUS-C: is the hbar transition phase arg B realised by the forced pivot circle's spin lift (cap angle 2 pi (1-h))?
# Runnable form: the only hbar-flavoured angle the record defines ON THE GRAM CELL is arg B = arg(1 + sum gamma + iV) (PRED-1,
# bargmann.py).  Step 1 identifies it geometrically (no name used until the comparison line).  Step 2 compares with the cap.
# Step 3 states the relation that DOES hold.  Kill (Will, VIEW-1 doc s.3): identification dies if the two disagree on the cell.
import sys, sympy as sp, numpy as np
CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else "")); sys.stdout.flush()
def unit(v): return v/np.linalg.norm(v)
def solid_angle_gb(a, b, c):
    """solid angle of the geodesic triangle by Gauss-Bonnet: sum of interior angles minus pi (independent of any formula)."""
    def ang(p, q, r):                      # interior angle at p between geodesics to q and r
        tq = unit(q - np.dot(q, p)*p); tr = unit(r - np.dot(r, p)*p)
        return np.arccos(np.clip(np.dot(tq, tr), -1, 1))
    return ang(a, b, c) + ang(b, c, a) + ang(c, a, b) - np.pi
def argB(a, b, c):
    V = np.dot(a, np.cross(b, c)); S = 1 + np.dot(a, b) + np.dot(b, c) + np.dot(c, a)
    return np.angle(S + 1j*V)
rng = np.random.default_rng(11)
frames = []
while len(frames) < 300:
    f = [unit(v) for v in rng.normal(size=(3, 3))]
    if np.dot(f[0], np.cross(f[1], f[2])) < 0: f[1], f[2] = f[2], f[1]          # orient V > 0 (the deck is sgn V; fixed here)
    # keep the pole inside the triangle (acute spherical triangle) for the unsigned decomposition; the signed one is a corollary
    A = (np.cross(f[0], f[1]) + np.cross(f[1], f[2]) + np.cross(f[2], f[0]))/2; n = unit(A)
    if solid_angle_gb(f[0], f[1], f[2]) < 2*np.pi and all(np.dot(n, np.cross(f[i], f[(i+1) % 3])) > 0 for i in range(3)):
        frames.append(f)
print("=== C-1  what arg B is, geometrically ===")
ok1 = all(abs(2*argB(*f) - solid_angle_gb(*f)) < 1e-9 for f in frames)
check("C-1a on 300 random frames 2 arg B = the solid angle of the GEODESIC TRIANGLE (a_c, a_hbar, a_G) computed by Gauss-Bonnet"
      " from its interior angles: arg B is the half solid angle of the three great-circle legs", ok1)
th = sp.symbols('theta', positive=True)                              # exact: the trine
e = [sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0]), sp.Matrix([0, 0, 1])]
Vt = e[0].dot(e[1].cross(e[2])); St = 1 + sum(e[i].dot(e[j]) for i in range(3) for j in range(i+1, 3))
check("C-1b exact at the trine: arg B = atan2(1, 1) = pi/4, twice it is pi/2 = the octant", sp.simplify(sp.atan2(Vt, St) - sp.pi/4) == 0)
print("=== C-2  compare with the forced circle's cap 2 pi (1 - h), h = V/2A ===")
def cap(f):
    A = (np.cross(f[0], f[1]) + np.cross(f[1], f[2]) + np.cross(f[2], f[0]))/2
    h = np.dot(f[0], np.cross(f[1], f[2]))/(2*np.linalg.norm(A)); return 2*np.pi*(1 - h), h
ratios = np.array([solid_angle_gb(*f)/cap(f)[0] for f in frames])
h_trine = 1/sp.sqrt(3); cap_trine = 2*sp.pi*(1 - h_trine)
check("C-2a exact at the trine: cap = 2 pi (1 - 1/sqrt 3) = %.6f, triangle = pi/2 = %.6f: NOT EQUAL" % (float(cap_trine), float(sp.pi/2)),
      sp.simplify(cap_trine - sp.pi/2) != 0)
check("C-2b across the cell the ratio triangle/cap is not constant (min %.3f, max %.3f): the transition phase is NOT the cap holonomy,"
      " not even up to a fixed factor.  CENSUS-C as an identification of angles is KILLED" % (ratios.min(), ratios.max()),
      ratios.max() - ratios.min() > 0.05 and ratios.max() < 1)
print("=== C-3  the relation that does hold: cap = triangle + three leg segments ===")
def decomposition(f):
    A = (np.cross(f[0], f[1]) + np.cross(f[1], f[2]) + np.cross(f[2], f[0]))/2; n = unit(A); h = np.dot(n, f[0])
    segs = []; dphis = []
    for i in range(3):
        p, q = f[i], f[(i+1) % 3]
        tp = unit(p - h*n); tq = unit(q - h*n)                              # azimuthal directions about the pole n
        dphi = np.arctan2(np.dot(n, np.cross(tp, tq)), np.dot(tp, tq)) % (2*np.pi)
        seg = dphi*(1 - h) - solid_angle_gb(n, p, q)                         # cap sector minus pole triangle
        segs.append(seg); dphis.append(dphi)
        gam = np.dot(p, q)
        assert abs(np.cos(dphi) - (gam - h**2)/(1 - h**2)) < 1e-9          # each segment is a function of (gamma_ij, h) only
    return sum(segs), sum(dphis), segs
ok3 = True; ok3b = True
for f in frames:
    ssum, dsum, segs = decomposition(f)
    ok3 &= abs(cap(f)[0] - solid_angle_gb(*f) - ssum) < 1e-9 and abs(dsum - 2*np.pi) < 1e-9
    ok3b &= all(s > -1e-12 for s in segs)
check("C-3a cap(h) - triangle = sum of three circular segments, each a function of ONE leg gamma_ij and h (cos dphi_ij = (gamma_ij - h^2)/(1 - h^2));"
      " azimuths sum to 2 pi.  NOTE (galois.py G-8/G-9): every angle here -- cap, triangle, segments -- lives in the sqrt(Delta) class;"
      " the conjugate triple sqrt(1+gamma_ij) enters only the legs' MODULI |<a_i|a_j>|, never an angle.  The split is moduli vs phase", ok3)
check("C-3b every segment is non-negative: the cap strictly contains the triangle; the hbar phase is LESS than the circle's holonomy by the legs", ok3b)
print("=== C-4  parity under the deck (the horizon paper's principle: invariants descend, odd data does not) ===")
# tau realised as a_i -> -a_i (V -> -V, Gram fixed).  Classify each cell function by parity.
def parity(fn):
    res = set()
    for f in frames:
        g = [-v for v in f]; x, y = fn(f), fn(g)
        res.add('odd' if abs(x + y) < 1e-9 else ('even' if abs(x - y) < 1e-9 else 'mixed'))
    return res
tab = {'arg B': parity(lambda f: argB(*f)),
       'h = V/2A': parity(lambda f: cap(f)[1]),
       'cap 2pi(1-h)': parity(lambda f: cap(f)[0]),
       '|B|': parity(lambda f: abs(np.exp(1j*argB(*f)))*np.sqrt(np.prod([1 + np.dot(f[i], f[j]) for i in range(3) for j in range(i+1, 3)])/8)),
       'segments (sum)': parity(lambda f: decomposition(f)[0])}
check("C-4a arg B and h are PURE ODD under the deck; |B| is EVEN; the cap 2 pi (1-h) is MIXED (2 pi even part, 2 pi h odd part) -- exactly the"
      " Theta / Theta^2 / T = Theta-over-root pattern of the horizon-pair paper.  The Gram quotient forgets sgn V; a model observable"
      " built from arg B DESCENDS only if the model supplies an ordering of the two sheets (the horizon analogue is S_+ > S_-)",
      tab['arg B'] == {'odd'} and tab['h = V/2A'] == {'odd'} and tab['|B|'] == {'even'} and tab['cap 2pi(1-h)'] == {'mixed'}, f"{tab}")
print("=== VERDICT ===")
print("  CENSUS-C as stated -- 'the hbar character is realised by C's spin lift with angle 2 pi (1-h)' -- is KILLED on the cell: the")
print("  record's hbar transition phase arg B is half the solid angle of the geodesic triangle through the three alignments, and the")
print("  cap of the forced circle exceeds it by three leg segments that never vanish off the branch locus.  What survives: (i) the")
print("  doubled-angle structure IS realised on the cell -- arg B is a half-angle of a solid angle (spinor of a vector), no circle needed;")
print("  (ii) cap, triangle and segments all live in the sqrt(Delta) class (G-9); the conjugate triple lives in the legs' moduli (G-8);")
print("  the circle carries the class datum h alone, the triangle carries the same class through Sum gamma as well; (iii) both arg B and")
print("  h are deck-odd, so VIEW-1's circle and the transition phase are the SAME sector's objects -- the orientation sector -- and the")
print("  census's 2:1 is the square class itself, not a locus.  The deck test now reads: find the model's ordering of the two sheets,")
print("  or arg B's sign does not descend (galois_horizon_cover: the sign of Theta is the datum the Galois quotient forgets; extremal")
print("  chirality = death of the odd sector at the branch point = R -> -I at Delta = 0 here).  [C-1a, C-3, C-4 numeric on 300 frames]")
n_pass = sum(CH); print(f"\n{n_pass}/{len(CH)} checks passed"); sys.exit(0 if all(CH) else 1)

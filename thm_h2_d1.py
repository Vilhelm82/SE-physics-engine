#!/usr/bin/env python3
# =============================================================================
# THM-TARGET H, part 2 -- D1 two-seat interior on the seat-layer pinning
# Date: 2026-09-02
#
# QUESTION: does ANY rod-free click-invariant diverge at the centre mu = 0?
# Named falsifier of "no singularity exists" (Will, 2026-09-02).
#
# ALLOWED INPUTS (each carried with its label):
#   BARE-1   pivot = hyperbolic rotor A = exp(lambda K/2), two-sided X -> A X A [proved, D1]
#   KIN-2a   pinning tanh(lambda) = sqrt(r_s/r)                                [DECLARED]
#   CONT-1   interior sheet lambda = mu + i pi/2, holomorphic continuation      [DECLARED]
#   RULING-2 (Will, 2026-09-02): K_C ~ a1+a2+a3 is the WORLD layer (the hole);
#            K_A = root is the SEAT layer (the infalling seat), projected out
#            of C. Part 2 runs A; C runs alongside as control, since every
#            identity below is K-independent in form.                         [DECLARED]
#   H-5/H-7  presented Gram G' = G + sinh^2(l) k k^T, det G' = Delta cosh^2(l) [proved, part 1]
#
# DEFINITION (rod-free): built from the rotor family, the frame's invariant
#   Gram, and angular separations, with NO division by a presented length,
#   by r, by dr, or by r*eps.
#
# BANNED (comparison stage only): tidal tensors, geodesic deviation, curvature
#   invariants, Kruskal/Penrose language, Kerr, regular-BH metrics.
#
# TWO PATHS for the relative rotor: Path 1 Cl(3) in the Pauli rep; Path 2 the
#   4x4 Lorentz matrices, linked by tr(L) = |tr M|^2. Share nothing but sympy.
#
# OUTCOME FORK, declared in advance:
#   (a) some rod-free click-invariant diverges at mu = 0: the singularity is
#       IN the model; "no singularity" dies at home.
#   (b) every rod-free invariant is regular at mu = 0 and every divergence
#       there factors through a rod: "no singularity" is a theorem of the
#       model, conditional on the labels above.
#   (c) undecidable without a new rule: named debt.
# =============================================================================
import sys
import sympy as sp
import mpmath as mp

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" -- {n}" if n else ""))
def z(e):
    e = sp.sympify(e)
    for r in (sp.simplify, lambda q: sp.simplify(sp.expand(q)),
              lambda q: sp.simplify(sp.expand(q.rewrite(sp.exp)))):
        try:
            if r(e) == 0: return True
        except Exception: pass
    return False
def zM(M): return all(z(e) for e in M)
def zc(e):
    e = sp.sympify(e)
    for r_ in (lambda q: sp.simplify(sp.expand_complex(sp.expand(q.rewrite(sp.exp)))),
               lambda q: sp.simplify(sp.expand_complex(sp.expand_trig(q))),
               lambda q: sp.simplify(q.rewrite(sp.exp))):
        try:
            if r_(e) == 0: return True
        except Exception: pass
    return False
def zcM(M): return all(zc(e) for e in M)

lam = sp.Symbol('lambda', positive=True)      # exterior rapidity
mu  = sp.Symbol('mu', positive=True)          # interior real rapidity
r, rs = sp.symbols('r r_s', positive=True)
g12, g13, g23 = sp.symbols('g12 g13 g23', real=True)
G = sp.Matrix([[1, g12, g13], [g12, 1, g23], [g13, g23, 1]])
Delta = sp.expand(G.det())
one = sp.Matrix([1, 1, 1])
K_OPTS = {
    "A (K = root a3)":      sp.Matrix([g13, g23, 1]),
    "C (K ~ a1+a2+a3)":     G * one / sp.sqrt((one.T * G * one)[0, 0]),
}
f = rs / (r - rs)                              # the ONE coefficient, both sheets

print("=" * 78); print("PART A -- one formula through the horizon"); print("=" * 78)
check("D-1  exterior: sinh^2(l) = r_s/(r - r_s) under tanh^2(l) = r_s/r",
      z(sp.simplify(sp.sinh(lam)**2 - (rs/(r-rs)).subs(r, rs/sp.tanh(lam)**2))))
check("D-2  interior: -cosh^2(mu) = r_s/(r - r_s) under r = r_s tanh^2(mu)  (CONT-1)",
      z(sp.simplify(-sp.cosh(mu)**2 - (rs/(r-rs)).subs(r, rs*sp.tanh(mu)**2))),
      "so G'(r) = G + [r_s/(r - r_s)] k k^T on BOTH sheets: one rational formula")
for name, k in K_OPTS.items():
    check(f"D-3  |K|^2 = k^T G^-1 k = 1   [{name}]", z(sp.simplify((k.T * G.inv() * k)[0, 0] - 1)))
    Gp = G + f * k * k.T
    check(f"D-4  det G'(r) = Delta r/(r - r_s), both sheets   [{name}]",
          z(sp.simplify(Gp.det() - Delta * r / (r - rs))),
          "simple zero at the centre, simple pole at the horizon, nothing else")
    e1 = sp.simplify(Gp.trace()); e2 = sp.simplify((Gp.trace()**2 - (Gp*Gp).trace()) / 2); e3 = sp.simplify(Gp.det())
    def only_pole_at_rs(e):
        den = sp.Poly(sp.denom(sp.together(e)), r)
        return set(sp.roots(den).keys()) <= {rs}
    poles_ok = all(only_pole_at_rs(e) for e in (e1, e2, e3))
    check(f"D-5  click-invariants tr, e2, det of G'(r) are rational with NO pole at r = 0   [{name}]", poles_ok,
          f"values at the centre: tr={sp.simplify(e1.subs(r,0))}, e2={sp.factor(e2.subs(r,0))}, det={e3.subs(r,0)}")
    d1 = sp.diff(Gp, r).subs(r, 0); d2 = sp.diff(Gp, r, 2).subs(r, 0)
    check(f"D-6  dG'/dr(0) = -k k^T/r_s and d2G'/dr2(0) = -2 k k^T/r_s^2: per-rod, yet FINITE (G' is even in mu)   [{name}]",
          zM(d1 + k * k.T / rs) and zM(d2 + 2 * k * k.T / rs**2))

print(); print("=" * 78); print("PART B -- the fold at the centre, and its deck involution"); print("=" * 78)
mu_of_r = sp.atanh(sp.sqrt(r / rs))
dmu = sp.simplify(sp.diff(mu_of_r, r))
check("D-7  fold: mu(r) = artanh(sqrt(r/r_s)) ~ sqrt(r/r_s); dmu/dr ~ r^(-1/2) at the centre",
      z(sp.simplify(sp.limit(sp.sqrt(r) * dmu, r, 0) - 1 / (2 * sp.sqrt(rs)))),
      "sqrt(r) * dmu/dr -> 1/(2 sqrt r_s): the divergence IS the square-root branch of r = r_s tanh^2(mu)")
s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -sp.I], [sp.I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]]); Id = sp.eye(2)
SIG = [s1, s2, s3]
def vec(v): return sum((v[i] * SIG[i] for i in range(3)), sp.zeros(2))
def rotor(L, Kv): return sp.cosh(L / 2) * Id + sp.sinh(L / 2) * vec(Kv)
e3 = sp.Matrix([0, 0, 1])
Ain = rotor(mu + sp.I * sp.pi / 2, e3)
deck = rotor(-mu + sp.I * sp.pi / 2, e3)
conj_flip = rotor(mu + sp.I * sp.pi / 2, -e3).applyfunc(lambda e: sp.conjugate(e))
check("D-8  deck involution of the fold, mu -> -mu, acts on the rotor as (I -> -I) composed with (K -> -K)",
      zM((deck - conj_flip).applyfunc(lambda e: sp.expand(sp.expand_trig(e).rewrite(sp.exp)))),
      "the centre's cover is the orientation Z2 and the pole-reversal Z2 acting TOGETHER")
check("D-9  the presented Gram is even in mu (deck-invariant): it lives on the base, the rotor on the cover",
      z(sp.simplify((G - sp.cosh(mu)**2 * K_OPTS["A (K = root a3)"] * K_OPTS["A (K = root a3)"].T)
                    .subs(mu, -mu).det() + Delta * sp.sinh(mu)**2)))

print(); print("=" * 78); print("PART C -- two seats: naive product vs frame-transported relative rotor"); print("=" * 78)
eps = sp.Symbol('epsilon', real=True)
lam1, lam2 = sp.symbols('lambda1 lambda2', real=True)
R = sp.cos(eps / 2) * Id - sp.sin(eps / 2) * s3 * s1            # rotor in the plane e3^e1
Rt = sp.cos(eps / 2) * Id + sp.sin(eps / 2) * s3 * s1           # its reverse
K1 = e3
K2m = sp.simplify(R * vec(K1) * Rt)                              # K2 = R K1 R~
K2 = sp.Matrix([sp.simplify(sp.trace(K2m * S) / 2) for S in SIG])
check("C-0  R rotates K1 = e3 to K2 = (sin e, 0, cos e), |K2| = 1, R R~ = 1",
      zM(K2 - sp.Matrix([sp.sin(eps), 0, sp.cos(eps)])) and zM(sp.simplify(R * Rt - Id)))
A1, A2 = rotor(lam, K1), rotor(lam, K2)
Nv = sp.simplify(A2 * A1.inv())                                  # NAIVE: frames pretended coincident
scalarN = sp.simplify(sp.trace(Nv) / 2)
closedN = 1 + 2 * sp.sinh(lam / 2)**2 * sp.sin(eps / 2)**2
check("D-10 naive product A2 A1^-1: scalar part = 1 + 2 sinh^2(l/2) sin^2(e/2)   [path 1]",
      z(sp.simplify(scalarN - closedN)))
def L4(L, Kv):
    ch, sh = sp.cosh(L), sp.sinh(L); M = sp.zeros(4, 4); M[0, 0] = ch
    for i in range(3):
        M[0, i+1] = sh * Kv[i]; M[i+1, 0] = sh * Kv[i]
        for j in range(3): M[i+1, j+1] = (1 if i == j else 0) + (ch - 1) * Kv[i] * Kv[j]
    return M
trL = sp.simplify((L4(lam, K2) * L4(lam, K1).inv()).trace())
check("D-11 tr(L2 L1^-1) = 4 (1 + 2 sinh^2(l/2) sin^2(e/2))^2   [path 2, 4x4; tr L = |tr M|^2]",
      z(sp.simplify(trL - 4 * closedN**2)))
mp.mp.dps = 40; worst = mp.mpf(0)
for lv, ev in ((mp.mpf('0.7'), mp.mpf('0.9')), (mp.mpf('2.3'), mp.mpf('0.31')), (mp.mpf('4.1'), mp.mpf('1.7'))):
    c, s = mp.cosh(lv), mp.sinh(lv); k2 = [mp.sin(ev), 0, mp.cos(ev)]; k1 = [0, 0, 1]
    def L4n(k):
        M = mp.matrix(4, 4); M[0, 0] = c
        for i in range(3):
            M[0, i+1] = s * k[i]; M[i+1, 0] = s * k[i]
            for j in range(3): M[i+1, j+1] = (1 if i == j else 0) + (c - 1) * k[i] * k[j]
        return M
    t = sum((L4n(k2) * L4n(k1)**-1)[i, i] for i in range(4))
    worst = max(worst, abs(t - 4 * (1 + 2 * mp.sinh(lv / 2)**2 * mp.sin(ev / 2)**2)**2))
check("D-12 40-digit numerics of D-11 at three (l, e) points", worst < mp.mpf(10)**-25, f"worst {mp.nstr(worst, 3)}")

# --- the frame-transported relative rotor (the D1 fence, dissolved along the family) ---
A2K2 = rotor(lam2, K2); A1K1 = rotor(lam1, K1)
check("D-13 frame transport: A(l, K2) = R A(l, K1) R~  (a seat whose K is rotated is the rotated seat)",
      zM(sp.simplify(rotor(lam, K2) - R * rotor(lam, K1) * Rt)))
Rel = sp.simplify(A2K2 * R * A1K1.inv())
check("D-14 relative rotor of two seats = R * A(l2 - l1, K1): rotation by their angle x COLLINEAR boost by their rapidity gap",
      zM((Rel - R * rotor(lam2 - lam1, K1)).applyfunc(lambda e: sp.simplify(sp.expand(e.rewrite(sp.exp))))),
      "no non-collinear composition ever arises along the pinned family: the D1 fence does not bite here")
mu1, mu2 = sp.symbols('mu1 mu2', positive=True)
RelIn = sp.simplify(rotor(mu2 + sp.I*sp.pi/2, K2) * R * rotor(mu1 + sp.I*sp.pi/2, K1).inv())
check("D-15 INTERIOR: two interior seats are related by R * A(mu2 - mu1, K1): REAL rotation, REAL boost; the i pi/2 cancels",
      zM((RelIn - R * rotor(mu2 - mu1, K1)).applyfunc(lambda e: sp.simplify(sp.expand(e.rewrite(sp.exp))))))
check("D-16 at the CENTRE (mu1, mu2 -> 0) the relative rotor -> R: two seats near r = 0 are a pure rotation apart",
      zM(sp.simplify(RelIn.subs({mu1: 0, mu2: 0}) - R)))
c_q = (Id + sp.I * vec(K1)) / sp.sqrt(2)
x = sp.Symbol('x', real=True)
check("D-17 ACROSS the horizon (l1 real, l2 = mu + i pi/2): relative rotor = R * A(mu - l1, K1) * c, c = (1 + I K1)/sqrt2",
      zcM(rotor(x + sp.I*sp.pi/2, K1) - c_q * rotor(x, K1)),
      "one quarter turn, the same for every depth: the crossing cost is fixed, not a function of where you end up")
scalarN_in = closedN.subs(lam, mu + sp.I*sp.pi/2)
ghost_closed = 1 + (sp.I*sp.sinh(mu) - 1) * sp.sin(eps/2)**2
check("D-18 OFFICE GHOST #4: the naive product's scalar on the interior = 1 + (i sinh mu - 1) sin^2(e/2); finite at the centre (= cos^2(e/2)), divergent at the horizon",
      zc(scalarN_in - ghost_closed) and z(ghost_closed.subs(mu, 0) - sp.cos(eps/2)**2)
      and sp.limit(sp.im(ghost_closed).subs(eps, 1), mu, sp.oo) == sp.oo,
      "composing rotors of two seats WITHOUT the frame transport is 'two offices, one product'; its divergence is the ghost")

print(); print("=" * 78); print("PART D -- what actually happens at the centre on the seat layer (K = root)"); print("=" * 78)
kA = K_OPTS["A (K = root a3)"]
GpA = G + f * kA * kA.T
check("D-19 K = root: the root's presented vector is cosh(l) * a3 (pure rescaling); |v_root|^2 = r/(r - r_s)",
      z(sp.simplify(GpA[2, 2] - r / (r - rs))) and zM(sp.simplify(GpA.row(2).T - (r/(r-rs)) * sp.Matrix([g13, g23, 1]))),
      "zero at the centre, negative (timelike) for 0 < r < r_s, pole at the horizon")
DpA = sp.simplify(GpA.det() / (GpA[0, 0] * GpA[1, 1] * GpA[2, 2]))
DpA_closed = Delta / ((1 + f * g13**2) * (1 + f * g23**2))
gp12 = (g12 - g13 * g23) / sp.sqrt((1 - g13**2) * (1 - g23**2))
check("D-20 K = root: normalised Delta_pres = Delta / [(1 + f g13^2)(1 + f g23^2)] -- REGULAR and POSITIVE at the centre",
      z(sp.simplify(DpA - DpA_closed)) and z(sp.simplify(DpA_closed.subs(r, 0) - (1 - gp12**2))),
      "Delta_pres(centre) = 1 - gamma'_12^2 > 0: the presented DIRECTIONS never degenerate at r = 0; only the root's LENGTH passes through zero")
check("D-21 CORRECTION to part 1: the coplanar landing at the centre (H-19/20) is a generic-K feature; for K = root the seat's presented frame is non-degenerate there",
      not z(DpA_closed.subs(r, 0)),
      "part 1's centre verdict stands for option C (control: v_i = a_i - k_i K coplanar) and is REPLACED for option A by D-19/20")
kC = K_OPTS["C (K ~ a1+a2+a3)"]
GpC = G + f * kC * kC.T
check("D-22 control, K = click-symmetric: det G' -> 0 at the centre with all three presented lengths finite: the coplanar landing is real for C",
      z(sp.simplify(GpC.det().subs(r, 0))) and all(sp.simplify(GpC[i, i].subs(r, 0)) != 0 for i in range(3)))

print(); print("=" * 78); print("PART E -- every divergence at the centre is a rod in the denominator"); print("=" * 78)
check("E-1  dmu/dr diverges as r^(-1/2); sqrt(r) * dmu/dr is finite   [the fold]",
      z(sp.simplify(sp.limit(sp.sqrt(r) * dmu, r, 0) - 1/(2*sp.sqrt(rs)))))
check("E-2  angle per transverse rod separation, e/(r e) = 1/r, diverges; r * (1/r) = 1   [the rod itself]", True, "trivial, recorded so the list is complete")
check("E-3  the root's presented length |v_root| = sqrt(r/(r - r_s)) -> 0; anything normalised by it is ~ r^(-1/2) times finite",
      z(sp.simplify(sp.limit(sp.sqrt(r) / sp.sqrt(-GpA[2, 2]), r, 0) - sp.sqrt(rs))))
print("  ARGUMENT (prose, not counted): on the interior sheet the rotor family A(mu + i pi/2, K) is ENTIRE in mu,")
print("  the frame's invariant Gram is constant, R(eps) is entire in eps. Every polynomial click-invariant of relative")
print("  rotors and presented Grams is therefore entire in mu near 0, hence finite at the centre. A divergence at the")
print("  centre requires a denominator vanishing there; the only such denominators along the family are rod lengths:")
print("  r itself (the fold, E-1/E-2) and presented axis lengths (E-3). The horizon is different: it is Re(lambda) -> oo,")
print("  infinite rapidity FOLDED to finite r by the pinning -- a rod-free quantity can diverge there because the")
print("  static seat itself ceases to exist (its rotor becomes projectively null), not because any frame quantity blows up.")

print(); n_ = sum(CH); print(f"RESULT: {n_}/{len(CH)} checks passed.")
print("=" * 78)
if n_ != len(CH):
    print("VERDICT: check failures above -- do NOT read the verdict below as earned."); sys.exit(1)
print("VERDICT: OUTCOME (b). No rod-free click-invariant diverges at the centre, on the seat layer OR the world layer.")
print("  One rational formula G'(r) = G + [r_s/(r - r_s)] k k^T runs through the horizon on both sheets (D-1..D-6).")
print("  Two seats along the family are related by (rotation by their angle) x (collinear boost by their rapidity gap);")
print("  inside, both are REAL and the relative rotor -> a pure rotation at the centre (D-14..D-16). Crossing the")
print("  horizon costs one fixed quarter turn (D-17). On the seat layer the presented directions are non-degenerate at")
print("  r = 0 (D-20): only the root's rod reading passes through zero (D-19). Every divergence at the centre is a rod")
print("  in a denominator (E-1..E-3). Conditional on KIN-2a, CONT-1, RULING-2 (all declared).")
print("=" * 78)
print("COMPARISON STAGE (banned names spoken here only):")
print("  - 'tidal divergence at r = 0' (geodesic deviation ~ 1/r^3): NOT reproduced and NOT computed -- the model has no")
print("    acceleration; its per-rod divergences at the centre are r^(-1/2) (fold) and 1/r (angular), both Jacobians.")
print("  - the seat-layer centre is a place where one rod reading is zero and everything else is finite: the shape of")
print("    a coordinate singularity, not a curvature one. Whether that is (a) artifact or (b) the model being too coarse")
print("    is now a single question: does the model own a rod-free curvature reading? It does not yet. Named debt: CURV-1.")
sys.exit(0)

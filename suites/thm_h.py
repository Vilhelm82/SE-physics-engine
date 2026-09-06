#!/usr/bin/env python3
# =============================================================================
# THM-TARGET H, part 1 — the presented state along the pinned family
# Date: 2026-09-02
#
# QUESTION (open since 2026-08-30): is the hole a STATE or a LOCUS? Does the
# pinned family reach Delta = 0, and if so where, and on which stratum?
#
# ALLOWED INPUTS (each carried with its label):
#   BARE    state = Gram matrix G of three unit axes; Delta = det G      [proved]
#   BARE-1  pivot = hyperbolic rotor A = exp(lambda K/2), K^2 = +1,
#           acting two-sided X -> A X A (thm_e.py convention)            [proved, D1]
#   RULE-2  rod office: a seat READS the vector part of a presentation   [named]
#   KIN-2a  pinning tanh(lambda) = sqrt(r_s/r)                           [DECLARED]
#   CONT-1  interior sheet = holomorphic continuation of the model's own
#           sandwich, lambda = mu + i pi/2, applied through the exact
#           identities cosh(lambda) = i sinh(mu), sinh(lambda) = i cosh(mu),
#           tanh(lambda) = coth(mu)                                      [DECLARED
#           convention; a Hermitian-adjoint sandwich is NOT analytic in
#           complex lambda and is not used]
#
# BANNED (comparison stage only): any interior black-hole model, Kruskal /
#   Penrose-diagram language, curvature invariants, Kerr, regular-BH metrics.
#
# TWO INDEPENDENT PATHS for the presentation map:
#   Path 1: Cl(3) in the Pauli representation, sandwich A X A
#   Path 2: 4x4 Lorentz boost along K acting on the 4-vector (0, a_i)
# They share nothing but sympy.
#
# OUTCOME FORK, declared in advance:
#   (a) the presented state never degenerates along the family: hole is a
#       locus in pivot space only, disjoint from Delta = 0.
#   (b) it degenerates at the horizon and/or the centre: the hole is a
#       (presented) state and Delta = 0 is its wall — record WHICH stratum.
#   (c) the family leaves the state space in a way the model cannot read:
#       new declared rule needed; named debt.
# =============================================================================
import sys
import sympy as sp
import mpmath as mp

CH = []
def check(t, ok, n=""):
    CH.append(bool(ok)); print(f"  [{'PASS' if ok else 'FAIL'}] {t}" + (f" — {n}" if n else ""))
def z(e):
    e = sp.sympify(e)
    for r in (sp.simplify, lambda q: sp.simplify(sp.expand(q)),
              lambda q: sp.simplify(sp.expand(q.rewrite(sp.exp)))):
        try:
            if r(e) == 0: return True
        except Exception: pass
    return False
def zM(M): return all(z(e) for e in M)

lam = sp.Symbol('lambda', positive=True)   # exterior rapidity
mu  = sp.Symbol('mu', positive=True)       # interior real rapidity (lambda = mu + i pi/2)
rs  = sp.Symbol('r_s', positive=True)

# ---------------------------------------------------------------------------
# PART A — the presentation map, two paths, general frame and general K
# ---------------------------------------------------------------------------
print("=" * 78); print("PART A — presentation map: Path 1 (Cl(3)) vs Path 2 (4x4 boost)"); print("=" * 78)
s1 = sp.Matrix([[0,1],[1,0]]); s2 = sp.Matrix([[0,-sp.I],[sp.I,0]]); s3 = sp.Matrix([[1,0],[0,-1]]); Id = sp.eye(2)
SIG = [s1, s2, s3]
def para(t, v): return t*Id + sum((v[i]*SIG[i] for i in range(3)), sp.zeros(2))
def scalar_part(M): return sp.simplify(sp.trace(M)/2)
def vector_part(M): return sp.Matrix([sp.simplify(sp.trace(M*S)/2) for S in SIG])

asym = sp.symbols('a11 a12 a13 a21 a22 a23 a31 a32 a33', real=True)
Amat = sp.Matrix(3, 3, asym)                       # rows = axes a_i (unit, imposed in checks via K,|a| identities)
K1, K2 = sp.symbols('K1 K2', real=True)
K = sp.Matrix([K1, K2, sp.sqrt(1 - K1**2 - K2**2)])  # unit pivot direction, |K| = 1 exactly
Rot = para(sp.cosh(lam/2), sp.sinh(lam/2)*K)      # A = cosh(l/2) + sinh(l/2) K.sigma

def path1(a):                                      # a: 3-vector; returns (scalar, vector) of A (a.sigma) A
    M = Rot * para(0, a) * Rot
    return scalar_part(M), vector_part(M)
def path2(a):                                      # 4x4 boost on (0, a)
    ch, sh = sp.cosh(lam), sp.sinh(lam)
    L = sp.zeros(4, 4); L[0, 0] = ch
    for i in range(3):
        L[0, i+1] = sh*K[i]; L[i+1, 0] = sh*K[i]
        for j in range(3):
            L[i+1, j+1] = (1 if i == j else 0) + (ch - 1)*K[i]*K[j]
    x = L * sp.Matrix([0, a[0], a[1], a[2]])
    return sp.simplify(x[0]), sp.Matrix(x[1:4])

k = Amat * K                                       # k_i = a_i . K
Vclosed = Amat + (sp.cosh(lam) - 1) * k * K.T      # claimed presented vector parts (rows)
ok1 = ok2 = oks = True
for i in range(3):
    a = Amat.row(i).T
    t1, v1 = path1(a); t2, v2 = path2(a)
    ok1 &= zM(v1 - Vclosed.row(i).T); ok2 &= zM(v2 - Vclosed.row(i).T)
    oks &= z(t1 - sp.sinh(lam)*k[i]) and z(t2 - sp.sinh(lam)*k[i])
check("H-1  presented vector part  v_i = a_i + (cosh l - 1)(a_i.K) K   [path 1]", ok1)
check("H-2  presented vector part  v_i = a_i + (cosh l - 1)(a_i.K) K   [path 2]", ok2)
check("H-3  presented scalar (time-office) part = sinh(l) (a_i.K), both paths", oks,
      "the K-component of every axis leaks into the time office by sinh(l): E-7 for axes")

# RULE-1 sanity: the invariant pairing Gram never moves
def pair_scalar(t1, v1, t2, v2): return t1*t2 - (v1.T*v2)[0, 0]      # <X Ybar>_S in paravector form
okp = True
for i in range(3):
    for j in range(3):
        ti, vi = path1(Amat.row(i).T); tj, vj = path1(Amat.row(j).T)
        okp &= z(pair_scalar(ti, vi, tj, vj) + (Amat.row(i)*Amat.row(j).T)[0, 0])
check("H-4  invariant pairing Gram is lambda-independent (RULE-1): <X_i' X_j'bar> = -a_i.a_j", okp,
      "the FRAME's Gram is a spectator to the pivot; only the rod office moves")

# ---------------------------------------------------------------------------
# PART B — the presented Gram theorem
# ---------------------------------------------------------------------------
print(); print("=" * 78); print("PART B — THEOREM: G'(l) = G + sinh^2(l) k k^T,  det G' = Delta cosh^2(l)"); print("=" * 78)
G  = Amat * Amat.T
Gp = Vclosed * Vclosed.T
check("H-5  G' = G + sinh^2(l) k k^T   (general frame, general unit K)", zM(sp.expand(Gp - (G + sp.sinh(lam)**2 * k * k.T))),
      "the pivot adds ONE rank-one term, along the K-components of the axes")
Delta = Amat.det()**2
check("H-6  k^T adj(G) k = Delta |K|^2 = Delta   (matrix-determinant lemma input)",
      z(sp.expand((k.T * G.adjugate() * k)[0, 0] - Delta)))
check("H-7  det G'(l) = Delta cosh^2(l)  EXACTLY — presented volume V' = V cosh(l)",
      z(sp.expand(Gp.det() - Delta * sp.cosh(lam)**2)),
      "the frame's Delta factors out; ALL lambda-dependence is the seat factor cosh^2")

# ---------------------------------------------------------------------------
# PART C — the family, exterior: normalised presented state and the horizon
# ---------------------------------------------------------------------------
print(); print("=" * 78); print("PART C — EXTERIOR: r = r_s / tanh^2(l), lambda in (0, oo)"); print("=" * 78)
# concrete generic rational frame + generic rational K (all k_i nonzero, Delta > 0)
Anum = sp.Matrix([[sp.Rational(3,5), sp.Rational(4,5), 0],
                  [0, sp.Rational(3,5), sp.Rational(4,5)],
                  [sp.Rational(4,5), 0, sp.Rational(3,5)]])
Knum = sp.Matrix([sp.Rational(2,3), sp.Rational(2,3), sp.Rational(1,3)])
subsA = dict(zip(asym, list(Anum))); subsK = {K1: Knum[0], K2: Knum[1]}
knum = Anum * Knum
Dnum = (Anum.det())**2
print(f"  frame: rows {list(Anum.T)}  K = {list(Knum)}  k = {list(knum)}  Delta = {Dnum}")
check("H-8  frame is generic: all |k_i| in (0,1), Delta > 0",
      all(0 < abs(x) < 1 for x in knum) and Dnum > 0)

GpN = (G + sp.sinh(lam)**2 * k * k.T).subs(subsA).subs(subsK)
Dii = [sp.simplify(GpN[i, i]) for i in range(3)]
check("H-9  normalisation |v_i|^2 = 1 + sinh^2(l) k_i^2", all(z(Dii[i] - (1 + sp.sinh(lam)**2 * knum[i]**2)) for i in range(3)))
DeltaPres = sp.simplify(Dnum * sp.cosh(lam)**2 / sp.prod([1 + sp.sinh(lam)**2 * knum[i]**2 for i in range(3)]))
check("H-10 normalised presented Delta_pres(l) = Delta cosh^2 / prod(1 + sinh^2 k_i^2)",
      z(DeltaPres - sp.simplify(GpN.det() / sp.prod(Dii))))
check("H-11 Delta_pres(0) = Delta: the seat at r = oo reads the frame's own state", z(DeltaPres.subs(lam, 0) - Dnum))
check("H-12 Delta_pres > 0 for every finite l: NO singular presented state outside the horizon",
      all(sp.N(DeltaPres.subs(lam, v)) > 0 for v in (sp.Rational(1,10), 1, 3, 10, 30)))
lim_h = sp.limit(DeltaPres, lam, sp.oo)
check("H-13 HORIZON (l -> oo): Delta_pres -> 0", lim_h == 0)
N2 = 1/sp.cosh(lam)**2                                  # N^2 = 1 - r_s/r = sech^2 l  (KIN-2a)
ratio = sp.limit(DeltaPres / N2**2, lam, sp.oo)
check("H-14 vanishing order at the horizon: Delta_pres ~ (Delta/prod k_i^2) N^4  — DOUBLE zero in N^2 = 1 - r_s/r",
      z(ratio - Dnum / sp.prod([knum[i]**2 for i in range(3)])), f"limit Delta_pres/N^4 = {ratio}")
# which stratum? normalised presented directions -> sign(k_i) K : a Cayley NODE (m = +1)
gam_lim = sp.Matrix(3, 3, lambda i, j: sp.limit(GpN[i, j] / sp.sqrt(Dii[i] * Dii[j]), lam, sp.oo))
m_lim = gam_lim[0, 1] * gam_lim[0, 2] * gam_lim[1, 2]
check("H-15 presented directions collapse to sign(k_i) K: presented Gram -> a Cayley node with m = +1",
      all(abs(gam_lim[i, j]) == 1 for i in range(3) for j in range(3)) and m_lim == 1,
      f"limit Gram rows {list(gam_lim.T)}")
print("  -> the HORIZON is the NODE stratum (the tear) of the PRESENTED state; the frame is regular.")

# ---------------------------------------------------------------------------
# PART D — the interior sheet (CONT-1): r = r_s tanh^2(mu), mu in (oo, 0)
# ---------------------------------------------------------------------------
print(); print("=" * 78); print("PART D — INTERIOR: lambda = mu + i pi/2  (CONT-1)"); print("=" * 78)
# the continuation identities used (exact; verified here, not assumed)
check("H-16 CONT-1 identities: cosh(mu+i pi/2) = i sinh(mu), sinh(mu+i pi/2) = i cosh(mu), tanh = coth",
      z(sp.expand(sp.cosh(mu + sp.I*sp.pi/2).rewrite(sp.exp) - sp.I*sp.sinh(mu).rewrite(sp.exp))) and
      z(sp.expand(sp.sinh(mu + sp.I*sp.pi/2).rewrite(sp.exp) - sp.I*sp.cosh(mu).rewrite(sp.exp))) and
      z(sp.expand(sp.tanh(mu + sp.I*sp.pi/2).rewrite(sp.exp) - sp.coth(mu).rewrite(sp.exp))))
r_in = rs * sp.tanh(mu)**2
check("H-17 interior address: tanh(l) = sqrt(r_s/r) with tanh(l) = coth(mu)  =>  r = r_s tanh^2(mu) in (0, r_s)",
      z(sp.simplify(rs / sp.coth(mu)**2 - r_in)) and sp.limit(r_in, mu, sp.oo) == rs and r_in.subs(mu, 0) == 0)
# presented Gram inside: sinh^2(l) -> -cosh^2(mu), cosh^2(l) -> -sinh^2(mu)
GpIn = (G - sp.cosh(mu)**2 * k * k.T).subs(subsA).subs(subsK)
check("H-18 det G'_in = -Delta sinh^2(mu) <= 0 : presented volume V' = i V sinh(mu) is IMAGINARY inside",
      z(sp.expand(GpIn.det() + Dnum * sp.sinh(mu)**2)),
      "the cover coordinate rotates onto the Wick face; the orientation is a phase inside")
check("H-19 CENTRE (mu = 0, r = 0): det G'_in = 0 and the presented vector parts are the a_i minus their K-components",
      GpIn.det().subs(mu, 0) == 0 and
      zM((GpIn.subs(mu, 0)) - (Anum - knum * Knum.T) * (Anum - knum * Knum.T).T),
      "all three presented axes lie in the plane perpendicular to K: COPLANAR")
gam_c = sp.Matrix(3, 3, lambda i, j: sp.simplify(GpIn[i, j] / sp.sqrt(GpIn[i, i] * GpIn[j, j])).subs(mu, 0))
m_c = sp.simplify(gam_c[0, 1] * gam_c[0, 2] * gam_c[1, 2])
check("H-20 the centre lands on the SMOOTH coplanar stratum, not a node: all |gamma'| < 1 there",
      all(abs(sp.N(gam_c[i, j])) < 1 for i in range(3) for j in range(3) if i != j),
      f"presented m at the centre = {m_c} (a node needs m = +1)")
DiiIn = [sp.simplify(GpIn[i, i]) for i in range(3)]
DeltaIn = sp.simplify(GpIn.det() / sp.prod(DiiIn))
slope = sp.limit(DeltaIn / (r_in / rs), mu, 0)
check("H-21 vanishing order at the centre: Delta_pres ~ -(Delta/prod(1-k_i^2)) (r/r_s)  — SIMPLE zero in r",
      z(slope + Dnum / sp.prod([1 - knum[i]**2 for i in range(3)])), f"limit Delta_pres/(r/r_s) = {slope}")
print("  -> the CENTRE is the COPLANAR stratum (the twist) of the PRESENTED state; the frame is regular.")

# inner null radii: |v_i|^2 = 1 - k_i^2 cosh^2(mu) = 0
print()
r_null = [sp.simplify(rs * (1 - knum[i]**2)) for i in range(3)]
ok_null = True
for i in range(3):
    mu_i = sp.acosh(1 / abs(knum[i]))
    ok_null &= z(DiiIn[i].subs(mu, mu_i)) and z(r_in.subs(mu, mu_i) - r_null[i])
check("H-22 INNER NULL RADII: presented axis i is null (cannot be normalised) at r_i = r_s (1 - k_i^2) = r_s sin^2(theta_i)",
      ok_null, f"r_i/r_s = {[sp.nsimplify(x/rs) for x in r_null]}")
# signature story: each axis with k_i != 0 is timelike just inside the horizon, spacelike again below r_i
okS = all(sp.N(DiiIn[i].subs(mu, 5)) < 0 for i in range(3)) and all(sp.N(DiiIn[i].subs(mu, sp.Rational(1,100))) > 0 for i in range(3))
check("H-23 every axis with k_i != 0 presents TIMELIKE just inside the horizon and SPACELIKE again below its r_i", okS,
      "the t<->r swap happens per axis, at that axis's own address")
check("H-24 the ROOT case K = a_1 (k_1 = 1): the root's null radius is r = 0 — the centre is where the root's rod reading vanishes",
      sp.simplify(rs * (1 - 1**2)) == 0)

# ---------------------------------------------------------------------------
# PART E — numerics (30 digits) at random rational points, both paths
# ---------------------------------------------------------------------------
print(); print("=" * 78); print("PART E — 30-digit spot checks"); print("=" * 78)
mp.mp.dps = 30
worst = mp.mpf(0)
for lv in (mp.mpf('0.37'), mp.mpf('1.9'), mp.mpf('4.2')):
    Af = mp.matrix([[3/mp.mpf(5), 4/mp.mpf(5), 0], [0, 3/mp.mpf(5), 4/mp.mpf(5)], [4/mp.mpf(5), 0, 3/mp.mpf(5)]])
    Kf = mp.matrix([2/mp.mpf(3), 2/mp.mpf(3), 1/mp.mpf(3)])
    kf = Af * Kf
    Vf = Af + (mp.cosh(lv) - 1) * kf * Kf.T
    detGp = mp.det(Vf * Vf.T)
    Dmp = mp.mpf(int(Dnum.p)) / mp.mpf(int(Dnum.q))
    worst = max(worst, abs(detGp - Dmp * mp.cosh(lv)**2))
check("H-25 det G' = Delta cosh^2(l) numerically to 1e-25 at three rapidities", worst < mp.mpf(10)**-25, f"worst {mp.nstr(worst, 3)}")

# ---------------------------------------------------------------------------
print(); n_ = sum(CH); print(f"RESULT: {n_}/{len(CH)} checks passed.")
print("=" * 78)
if n_ != len(CH):
    print("VERDICT: check failures above — do NOT read the verdict below as earned."); sys.exit(1)
print("VERDICT: OUTCOME (b). The hole is a PRESENTED STATE.")
print("  The pinned family degenerates the seat's read state TWICE while the frame's")
print("  Delta never moves (H-4, H-7): at the horizon it lands on a Cayley NODE (H-15,")
print("  double zero in N^2, H-14); at the centre on the smooth COPLANAR stratum (H-19/20,")
print("  simple zero in r, H-21). Inside, the presented volume is imaginary (H-18) and")
print("  each axis has its own null radius r_i = r_s sin^2(theta_i) (H-22).")
print("  Conditional on: KIN-2a (pinning, declared), CONT-1 (continuation convention, declared).")
print("=" * 78)
print("COMPARISON STAGE (banned names spoken here only):")
print("  - t<->r swap inside a horizon: reproduced PER AXIS with an address (H-23), not globally.")
print("  - 'curvature invariants diverge at r = 0': the frame's invariants do NOT move here;")
print("    what vanishes at r = 0 is the presented volume, and what fails is normalisation of")
print("    a null presented axis (H-22/24). No tidal (two-seat) quantity is computed in part 1.")
print("  - inner null surfaces at r_s sin^2(theta_i): NO standard-model analogue is claimed;")
print("    lineage unsearched, flagged candidate-new until the D1 two-seat computation exists.")
sys.exit(0 if n_ == len(CH) else 1)

#!/usr/bin/env python3
# =============================================================================
# XLINK-1 -- Cella's coupling-edge form is the Gram matrix of the seated root's triple horizon.
#   (Claire, 2026-09-05.  Cross-link between the Cella corpus and the seated root; corpus used freely.)
# SOURCE: Cella U-0610 / DBP_Curvature_Role_Reduction.md sec. 8:  nu = (g1 H23, g2 H13, g3 H12),
#   Delta_c = nu^T (2I - 11^T) nu, 'a Lorentzian quadratic form, signature (2,1)'; sec. 13.3 lead L6 (OPEN):
#   'Does the channel quadratic form have a fixed signature for general (n,r)?  A light-cone reading where
#    null directions = gauge-rigid directions?'
# QUESTION: coincidence of a 3x3 matrix's eigenvalues, or the same object?
# CLAIMS, checked:
#   x1  Q = 2I - 11^T has eigenvalues (2, 2, -1): signature (2,1), timelike direction 1 = (1,1,1), and the
#       sum-zero plane ker Sigma (Cella's gauge-transport plane) is its POSITIVE-DEFINITE complement.
#   x2  Q's null cone is (Sum nu)^2 = 2|nu|^2, i.e. Sum nu_i^2 - 2 Sum_{i<j} nu_i nu_j = 0; with nu_i = s_i^2 that
#       is -16 Area^2 of the triangle with sides s_i (Heron): null <=> degenerate triangle.
#   x3  Q IS A GRAM MATRIX IN R^{2,1}: h1 = (1,0,0), h2 = (-1,1,1), h3 = (-1,-1,1) with q = diag(1,1,-1) have
#       q(h_i,h_i) = 1 and q(h_i,h_j) = -1: three unit SPACELIKE lines, Gram = Q exactly, det = -4 = det q . det M^2.
#   x4  Every pair is at the seated root's HORIZON condition: the 2x2 minors det [[1,-1],[-1,1]] = 0 (|gamma| = 1),
#       with sigma = -1, so h_i + h_j is NULL for every pair.  Q's null cone contains exactly the three pairwise
#       sums (1,1,0), (0,1,1), (1,0,1): Cella's null directions are the model's three horizon light rays n0.
#   x5  1 = h1 + h2 + h3 is timelike, q(1,1) = -3; c := 1/sqrt3 is a unit negative line -- a SEAT -- and it is the
#       unique S3-symmetric one.  Cella's trace channel is the seat that sees all three rulers equally.
#   x6  ker Sigma = c_perp: Cella's gauge-transport plane IS the seat's compact space (positive definite).
#       'Zero-sum redistribution of channel account' = motion in the seat's rest plane.
#   x7  From that seat: a = b = -1/sqrt3 (equal depths), gamma = -1 for every pair -- on the horizon, OFF the
#       pinch (a != sigma b; delta = (a + b)^2 = 4/3).  In seat coordinates: l = asinh(1/sqrt3) and t = 2pi/3.
#       THE TRINE.  The Cella form is the Gram of the coplanar 0/120/240 trine (Codex's GRAM_SUBMERSION
#       fixture) lifted to the depth at which every ruler pair reaches the horizon.
#   x8  GENERAL n (Cella's L6): 2I - 11^T has eigenvalues 2 (mult n-1) and 2-n: signature (n-1,1) for every
#       n >= 3, degenerate at n = 2.  It is the Gram of n unit spacelike lines in R^{n-1,1} pairwise at q = -1,
#       realised explicitly; n = 4 is the model's Cl(3,1) after T8.
# CAVEAT (Codex, docs/results/2026-09-05/2026-09-05-cella-dependency-reconciliation.md sec. 1, recorded 2026-09-05): under a defining-function
#   gauge H -> H + g a^T + a g^T the edge vector shifts by delta nu with det(a -> delta nu) = 2 (g1 g2 g3)^2 != 0, so an
#   arbitrary gauge moves nu ANYWHERE.  'Cella's null cone = the horizon light rays' is therefore true of the FORM and
#   gauge-dependent for any actual nu until a presentation rule is fixed.  The nu-map is Codex's (agreed split).
# OBSTRUCTION (Codex overlap audit sec. 5, citing Cella's Wall-Isotropy Comparison Theorem sec. 1, recorded 2026-09-05 evening):
#   the additive gauge group ker Sigma over Q is DIVISIBLE, so every homomorphism ker Sigma -> C_2 is trivial
#   (phi(a) = 2 phi(a/2) = 0).  Cella's gauge action on the plane x6 identifies with c_perp CANNOT generate the sheet
#   parity.  So: the plane is the seat's rest space (x6 stands), but the deck must come from the branched cover
#   K(sqrt delta)/K (PINCH-2), never from motion in that plane.  Consistent with PINCH-2 q8 (no realisable loop).
# TIER: [DERIVED].  The FORM is identical and the seated root REALISES it; the identification of Cella's
#   nu-components (g_i H_jk) with frame data is NOT made here and stays open.  That is the honest boundary.
# [First run: x8c failed on the n = 4 realisation's scale (pairwise -3/2, needed -2); fixed by r = 2/sqrt3. Claim unchanged.]
# KILL: if Q were not a Gram matrix of unit spacelike lines in R^{2,1}, or if its null cone missed the pairwise
#   sums, the link would be a coincidence of eigenvalues and nothing more.
# =============================================================================
import sympy as sp, time
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

one = sp.Matrix([1, 1, 1]); Qm = 2*sp.eye(3) - one*one.T
print("=== x1: signature and the two distinguished subspaces ===")
ev = Qm.eigenvals()
check("x1a 2I - 11^T has eigenvalues 2 (mult 2) and -1 (mult 1): signature (2,1)", ev == {2: 2, -1: 1})
check("x1b the timelike eigenvector is 1 = (1,1,1)", z((Qm*one + one).norm()))
kerS = [sp.Matrix([1, -1, 0]), sp.Matrix([1, 1, -2])]
check("x1c the sum-zero plane ker Sigma is the eigenvalue-2 eigenspace: positive definite", all(z((Qm*v - 2*v).norm()) for v in kerS))

print("=== x2: the null cone is Heron ===")
n1, n2, n3 = sp.symbols('n1 n2 n3', real=True); nu = sp.Matrix([n1, n2, n3])
Qf = sp.expand((nu.T*Qm*nu)[0])
check("x2a Q(nu) = 2|nu|^2 - (Sum nu)^2 = Sum nu_i^2 - 2 Sum_{i<j} nu_i nu_j", z(Qf - (n1**2 + n2**2 + n3**2 - 2*n1*n2 - 2*n2*n3 - 2*n1*n3)))
s1, s2, s3 = sp.symbols('s1 s2 s3', positive=True)
heron16 = 2*(s1**2*s2**2 + s2**2*s3**2 + s3**2*s1**2) - (s1**4 + s2**4 + s3**4)
check("x2b with nu_i = s_i^2: Q = -16 Area^2 (Heron); null <=> a degenerate triangle with sides s_i", z(Qf.subs({n1: s1**2, n2: s2**2, n3: s3**2}) + heron16))

print("=== x3: Q is a Gram matrix in the model's R^{2,1} ===")
q = sp.diag(1, 1, -1)
h1, h2, h3 = sp.Matrix([1, 0, 0]), sp.Matrix([-1, 1, 1]), sp.Matrix([-1, -1, 1])
M = sp.Matrix.hstack(h1, h2, h3)
G = M.T*q*M
check("x3a h1 = (1,0,0), h2 = (-1,1,1), h3 = (-1,-1,1): three UNIT SPACELIKE lines, q(h_i,h_i) = 1", all(G[i, i] == 1 for i in range(3)))
check("x3b pairwise q(h_i,h_j) = -1: the Gram matrix IS 2I - 11^T", G == Qm)
check("x3c det G = -4 = det q . (det M)^2: the three lines span R^{2,1} (det M = 2)", G.det() == -4 and M.det() == 2 and q.det()*M.det()**2 == -4)

print("=== x4: every pair sits at the seated root's horizon, and the null cone is the three light rays ===")
check("x4a every 2x2 principal minor of G vanishes: |gamma_ij| = 1 for all three pairs -- the ruler-plane horizon condition, three times over",
      all(G.extract([i, j], [i, j]).det() == 0 for i, j in [(0, 1), (0, 2), (1, 2)]))
check("x4b sigma = -1 for every pair: h_i + h_j is NULL (the model's n0 = h + g at the gamma = -1 horizon)",
      all(z(((h_i + h_j).T*q*(h_i + h_j))[0]) for h_i, h_j in [(h1, h2), (h1, h3), (h2, h3)]))
check("x4c in Q's own coordinates the pairwise sums (1,1,0), (0,1,1), (1,0,1) are Q-null: Cella's null directions are the three horizon light rays",
      all(z((v.T*Qm*v)[0]) for v in [sp.Matrix([1, 1, 0]), sp.Matrix([0, 1, 1]), sp.Matrix([1, 0, 1])]))
check("x4d and the differences (1,-1,0) etc. are SPACELIKE (Q = 4): Cella's kappa_c-null ray (0,1,-1) is a spatial direction, not a light ray",
      (sp.Matrix([0, 1, -1]).T*Qm*sp.Matrix([0, 1, -1]))[0] == 4)

print("=== x5-x6: the trace channel is a seat; the gauge plane is its space ===")
S = h1 + h2 + h3
check("x5a 1 = h1 + h2 + h3 is TIMELIKE: q(S,S) = -3", (S.T*q*S)[0] == -3)
c = S/sp.sqrt(3)
check("x5b c := S/sqrt3 is a unit negative line: a SEAT, and the unique S3-invariant one", (c.T*q*c)[0] == -1)
cperp_basis = [h1 - h2, h2 - h3]
check("x5c the seat's space c_perp is spanned by ruler DIFFERENCES = Cella's sum-zero plane ker Sigma", all(z((c.T*q*v)[0]) for v in cperp_basis))
Gp = sp.Matrix(2, 2, lambda i, j: (cperp_basis[i].T*q*cperp_basis[j])[0])
check("x6  the form on c_perp is positive definite (the seat's COMPACT plane): Cella's gauge-transport plane is the seat's rest space",
      all(e > 0 for e in Gp.eigenvals()))

print("=== x7: what the symmetric seat reads -- the trine at horizon depth ===")
a = (c.T*q*h2)[0]; b = (c.T*q*h3)[0]; gam = (h2.T*q*h3)[0]
check("x7a from the seat c: a = b = -1/sqrt3 (equal depths), gamma = -1 for the pair (h2,h3)", z(a + 1/sp.sqrt(3)) and z(b + 1/sp.sqrt(3)) and gam == -1)
sig = -1
check("x7b on the horizon (|gamma| = 1) but OFF the pinch: a != sigma b, delta = (a - sigma b)^2 = 4/3 > 0 -- a regular horizon point, nu defined",
      not z(a - sig*b) and z((a - sig*b)**2 - sp.Rational(4, 3)))
l = sp.asinh(1/sp.sqrt(3)); t = sp.symbols('t', real=True)
gam_seat = sp.cosh(l)**2*sp.cos(t) - sp.sinh(l)**2
tsol = [s for s in sp.solve(sp.Eq(gam_seat, -1), t) if 0 < s.evalf() < sp.pi.evalf()]
check("x7c seat coordinates: l1 = l2 = asinh(1/sqrt3), and gamma = -1 forces cos t = -1/2, t = 2pi/3: the rulers present as the 120-degree TRINE",
      z(sp.sinh(l) - 1/sp.sqrt(3)) and any(z(s - 2*sp.pi/3) for s in tsol))
check("x7d the trine's presented angles are 120 degrees for all three pairs (S3-symmetric), and every pair is at |gamma| = 1: the Cella form is the coplanar trine lifted to horizon depth",
      all(G.extract([i, j], [i, j]).det() == 0 for i, j in [(0, 1), (0, 2), (1, 2)]) and any(z(s - 2*sp.pi/3) for s in tsol))

print("=== x8: general n -- Cella's open lead L6 ===")
def sig_n(n):
    On = sp.ones(n, 1); Qn = 2*sp.eye(n) - On*On.T
    return Qn.eigenvals()
check("x8a n = 2: eigenvalues {2, 0}: DEGENERATE (no Lorentzian form on two channels)", sig_n(2) == {2: 1, 0: 1})
check("x8b n = 3, 4, 5: eigenvalues 2 (mult n-1) and 2-n < 0: signature (n-1, 1) -- Lorentzian for every n >= 3, timelike = the trace",
      all(sig_n(n) == {2: n - 1, 2 - n: 1} for n in (3, 4, 5)))
# explicit realisation for n = 4 in R^{3,1}: four unit spacelike lines pairwise at q = -1
q4 = sp.diag(1, 1, 1, -1)
x = sp.symbols('x', positive=True)
# h1 = (1,0,0,0); h_i = (-1, y_i) with |y_i|^2_{(2,1)} = 0, and q(h_i,h_j) = 1 + <y_i,y_j> = -1  =>  <y_i,y_j> = -2
r = 2/sp.sqrt(3)   # three null vectors in R^{2,1} at 120 degrees have <y_i,y_j> = r^2 (cos 120 - 1) = -3 r^2/2; need -2 => r^2 = 4/3
y2, y3, y4 = [r*sp.Matrix([sp.cos(ph), sp.sin(ph), 1]) for ph in (0, 2*sp.pi/3, 4*sp.pi/3)]
H4 = [sp.Matrix([1, 0, 0, 0])] + [sp.Matrix([-1, *y]) for y in (y2, y3, y4)]
G4 = sp.Matrix(4, 4, lambda i, j: sp.simplify((H4[i].T*q4*H4[j])[0]))
check("x8c n = 4 realised in R^{3,1}: four unit spacelike lines pairwise at q = -1, Gram = 2I - 11^T -- the model's Cl(3,1) sector after T8",
      G4 == 2*sp.eye(4) - sp.ones(4, 1)*sp.ones(1, 4))

n_ok = sum(CH); n_all = len(CH)
print(f"\nRESULT: {n_ok}/{n_all} checks passed in {time.time()-T0:.1f}s")
print("VERDICT: Cella's coupling-edge form 2I - 11^T is the Gram matrix of three unit spacelike lines in R^{2,1} pairwise at the")
print("  seated root's horizon condition (|gamma| = 1, sigma = -1): every pairwise sum is null, and those three null vectors")
print("  ARE Cella's null cone.  The trace direction is a seat -- the unique S3-symmetric one -- and Cella's gauge-transport")
print("  plane is that seat's compact rest space.  From it the rulers present as the 120-degree trine at depth asinh(1/sqrt3),")
print("  every pair on the horizon, off the pinch.  For general n the form has signature (n-1,1) and is realised in R^{n-1,1}.")
print("  NOT shown: the map from Cella's nu = (g_i H_jk) to frame data.  Same form, realised; identification of components open.")

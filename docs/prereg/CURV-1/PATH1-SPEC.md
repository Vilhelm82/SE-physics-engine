# CURV-1, Path 1 — the presented seat triangle and its excess along the pinned family
Written 2026-09-02 against HEAD `b4ac768`. For an independent runner. Do not read
`PREDICTIONS.md` in this folder before your suite is green and committed.

## 0. What is computed
A seat pinned at radius r reads its three axes as three presented vectors v_1, v_2, v_3.
Normalising them puts three points on the unit sphere; those points span a spherical
triangle, and the area of that triangle — its spherical excess Omega — is a reading
that uses angles only: rescaling any v_i leaves it unchanged. This spec computes that
excess, for every choice of poles of the presented axes, as an exact function of r
on both sheets, and takes its limits at the centre and at the horizon. The output
is the profile Omega(r), the endpoint limits, and their classification.

## 1. Inputs (frozen; cite, do not re-derive)
Symbols: r, r_s > 0; g12, g13, g23 real, G = [[1,g12,g13],[g12,1,g23],[g13,g23,1]]
positive definite; Delta = det G.

    f(r)    = r_s/(r - r_s)              exterior r > r_s: f > 0;  interior 0 < r < r_s: f < -1
    G'(r)   = G + f(r) k k^T             det G'(r) = Delta r/(r - r_s)             [D-1..D-4]
    rho_i(r) = G'_ii = 1 + f k_i^2 = (r - r_s(1 - k_i^2))/(r - r_s)   presented rods
    r_i     = r_s (1 - k_i^2)            null radii, rho_i(r_i) = 0, inside the horizon  [H-22]

Pivot options, k_i = a_i . K:

    A  (seat layer,  K = root = a_3)    k = (g13, g23, 1);   rho_3 = r/(r - r_s), r_3 = 0   [D-19]
    C  (world layer, K ~ a1+a2+a3)      k = G 1 / sqrt(1^T G 1),   1 = (1,1,1)^T

Where the rapidity is needed (route V, section 3): exterior l(r) = artanh sqrt(r_s/r);
interior mu(r) = artanh sqrt(r/r_s), l = mu + i pi/2 [CONT-1]. Presented vector part
v_i = a_i + (cosh l - 1) k_i K [H-1/H-2], with cosh l = i sinh mu inside. Dot products
of complex vectors are bilinear (sum of products, no conjugation): v.v = rho.

## 2. The excess orbit (the object)
Three lines through a point cut the unit sphere into eight triangles, in four
antipodal pairs; the pairs are indexed by the choice of pole on each axis. For a
sign vector eps = (e1, e2, e3), e_i = +-1, define on the presented Gram

    Gam_ij  = G'_ij / (sqrt(rho_i) sqrt(rho_j))                         presented cosines
    V'      = sqrt(det G') / (sqrt(rho_1) sqrt(rho_2) sqrt(rho_3))       presented volume
    S_eps   = 1 + e1 e2 Gam_12 + e1 e3 Gam_13 + e2 e3 Gam_23
    V_eps   = (e1 e2 e3) V'
    B_eps   = S_eps + i V_eps
    W_eps   = (S_eps + i V_eps) / (S_eps - i V_eps)                     = e^{i Omega_eps}
    Omega_eps = -i log W_eps                                            branch: section 4

S_{-eps} = S_eps and V_{-eps} = -V_eps, so Omega_{-eps} = -Omega_eps. The four members
with e1 e2 e3 = +1 are the unsigned excesses of the four triangles. Outside the
horizon everything is real and the four sum to 2 pi.

Square roots are the principal branch, taken for each rho_i separately
(sqrt(rho_1) sqrt(rho_2), never sqrt(rho_1 rho_2)). Inside the horizon some rho_i are
negative and the corresponding sqrt(rho_i) are imaginary. Choosing the other branch
of sqrt(rho_i) is the same as e_i -> -e_i: it permutes the orbit. The multiset
{B_eps : all eight eps} is therefore independent of every branch choice, and every
classification in this spec is made on that multiset, not on a single member.

Exact identity (branch-free, both sheets): with the product over the four members
having e1 e2 e3 = +1,

    prod B_eps = -4 (1 - Gam_12^2)(1 - Gam_13^2)(1 - Gam_23^2)      hence   prod W_eps = 1.
    (sign corrected 2026-09-02 evening: the original draft had +4; verified exact in curv1_path1.py P1-11a/b)

## 3. Two implementations, sharing only sympy
Route G (Gram). Build G'(r) as a rational matrix in r from section 1 and evaluate
section 2 directly. Everything before the square roots is rational in r.

Route V (vectors). Build the rotor A(l, K) = cosh(l/2) + sinh(l/2) K.sigma in the
Pauli representation (or the 4x4 boost), sandwich each concrete unit axis a_i to get
v_i(l), substitute the pinned rapidity of section 1, and compute the excess from the
vectors with the three-vector form

    tan(Omega/2) = [v1, v2, v3] / ( |v1||v2||v3| + (v1.v2)|v3| + (v1.v3)|v2| + (v2.v3)|v1| )

with |v| = sqrt(v.v) principal branch and [v1,v2,v3] = det. The eps signs are applied
to the vectors themselves, v_i -> e_i v_i. Route V never forms G'.

Agreement (must be able to fail): B_eps identical between routes at three rational
exterior r and three rational interior r — exactly where sympy closes, otherwise to
30 digits — for both K options and every eps.

## 4. The profile and the branch
Outside (r > r_s) everything is real and Omega_eps = 2 atan2(V_eps, S_eps).
Moving inward along the real r axis, each null radius r_i is a point where sqrt(rho_i)
passes through zero and turns imaginary. There the labelling of the orbit forks:
continuing sqrt(rho_i) through zero one way or the other differs by e_i -> -e_i.
Record the fork; do not resolve it. Track the eight B_eps(r) as continuous curves in
the complex plane (small steps, nearest-neighbour matching between consecutive r,
principal branch at each point) and report the label permutation at each r_i.

Report, for each (frame, K option):
  (a) the exact orbit at r -> oo (it must be the frame's own orbit);
  (b) the exact limits r -> r_s+ and r -> r_s- of every B_eps and W_eps. When S_eps
      and V_eps both vanish there, take the limit of W_eps along the family (sympy
      limit in r); if it does not exist, say so and print both expressions;
  (c) the exact limit r -> 0+ of every B_eps and W_eps, same 0/0 rule;
  (d) a numeric profile on r/r_s in {3, 2, 1.5, 1.2, 1.05, 1.01, 0.99, 0.95, 0.5,
      0.1, 0.01, 0.001} plus eight points inside each band between consecutive null
      radii, listing Re Omega_eps and Im Omega_eps for the four e1 e2 e3 = +1 members;
  (e) on each interior band, for each member, whether Omega_eps is real, purely
      imaginary, of the form pi + i(real), or general complex — decided exactly:
      |W_eps| = 1 means Omega real; W_eps real positive means Omega purely imaginary;
      W_eps real negative means Omega = pi + i(real).

## 5. Checks that can fail (minimum set; add, never remove)
P1-1  G'(r) reproduces D-4 (det), D-19 (rho_3 for option A) and H-22 (null radii).
P1-2  Omega_eps is unchanged by v_i -> c_i v_i for symbolic c_i > 0 (rod-free), and the
      orbit multiset is unchanged by any relabelling of the three axes.
P1-3  Exterior: all B_eps real; the four e1e2e3 = +1 excesses lie in (0, 2 pi) and sum
      to 2 pi (30 digits at three rational r); the identity prod B_eps = -4 prod (1 - Gam^2)
      holds as a rational identity in r.
P1-4  The same product identity holds on the interior sheet, exactly in r.
P1-5  Routes G and V agree (section 3).
P1-6  r -> oo returns the frame's own orbit.
P1-7  Option A: at r = 0 the root's rod is zero, so Gam_13 and Gam_23 are 0/0 there.
      Evaluate them as limits in r, never by substitution, and print the limit values.
P1-8  Option C: a1+a2+a3 is not a function of the three lines, so K_C moves with the
      pole labelling of the frame. Run the four pivot directions K ~ e1 a1 + e2 a2 + e3 a3
      (e1e2e3 = +1) on the original frame as four separate option-C runs, and report
      whether the endpoint classifications agree across them. Option A needs no such
      run: reversing the root leaves k k^T unchanged, and reversing a visible axis is
      a member of the orbit.
P1-9  Every endpoint classification is printed with its value, and the suite contains
      at least one check whose PASS/FAIL depends on which class came out — written so
      that the other classes make it FAIL and print the member responsible.

## 6. Frames
Run the general symbolic frame wherever sympy returns in reasonable time (the
endpoint limits usually do), and always these concrete rational frames:

    Frame P (all cosines positive, unequal)   a1 = (1,0,0)  a2 = (3/5,4/5,0)   a3 = (2/7,3/7,6/7)
                                              g12 = 3/5, g13 = 2/7, g23 = 18/35, Delta = (24/35)^2
    Frame Q (one obtuse angle)                a1 = (1,0,0)  a2 = (-3/5,4/5,0)  a3 = (2/7,3/7,6/7)
                                              g12 = -3/5, g13 = 2/7, g23 = 6/35
    Frame H (thm_h.py's frame, equal angles)  rows (3/5,4/5,0), (0,3/5,4/5), (4/5,0,3/5)
                                              g12 = g13 = g23 = 12/25   (control for coincidences)

For option A the root is a3 in every frame. For option C see P1-8.

## 7. The deliverable
One row per (frame, K option or pivot direction, endpoint), columns:
class | Omega_eps for the four e1e2e3 = +1 members | spin sign of each real B_eps
(sign of B_eps when it is real) | reality class on the adjacent band | forks crossed.
Endpoints: r -> oo, r -> r_s+, r -> r_s-, r -> 0+.
Write the table to `docs/prereg/CURV-1/PATH1-RESULTS.md`, the raw values to
`curv1_path1_results.json`, and commit both with the suite `curv1_path1.py` and its log.

## 8. Not used here
No black-hole metric, curvature invariant, tidal tensor, or interior model of any
kind enters this computation or its checks. Those belong to the comparison stage,
a separate document written after this suite is green.

## 9. Outcome fork for Path 1 (declared before running)
  (a) some member is unbounded at r -> 0+: the model owns a diverging rod-free
      curvature reading, and the centre is a singularity of the model itself.
  (b) every member is finite at r -> 0+: record identity / finite-nontrivial and the
      spin signs. The centre is then a finite curvature event and the recorded values
      are the model's own reading of it, one per layer (A and C).
  (c) indeterminate limits: the excess is not the right reading; named debt.

A (b) with identity on both layers and no -1 spin sign reads "no curvature event at
the centre". A (b) with a nontrivial value or a -1 sign is a finite event with an
address. Only (a) is a divergence. The three are distinguishable by the output of
sections 4(b)-(c); nothing in this spec presumes which one arrives.

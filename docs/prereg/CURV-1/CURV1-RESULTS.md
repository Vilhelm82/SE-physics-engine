# CURV-1 — RESULTS.  Does the model own a rod-free curvature reading, and what does it do at the centre?
Written 2026-09-02, evening.  Repo HEAD before this commit: `37424f7`.  All three suites green on the 7800X3D:
`curv1_path1.py` **40/40** (`run_curv1_p1.log`, ~5 min), `curv1_path23.py` **15/15** (`run_curv1_p23.log`, ~3 min).
Raw values: `curv1_path1_results.json`.  Re-running is the verification.

## Verdict in one paragraph
**Fork (b) at the centre, on every path, both layers.**  No rod-free curvature reading of the model is unbounded at r = 0.
The centre carries two finite readings, one per layer: the infalling seat (K = root) reads the frame's own dihedral angle
at the root, φ₃ = arccos[(γ₁₂ − γ₁₃γ₂₃)/(sinθ₁₃ sinθ₂₃)]; the hole (K ~ a₁+a₂+a₃) reads a full turn, 2π, whose spin lift
is −1 — the deck ℤ₂ — and every loop of seats around the centre returns exactly −1 in SU(2) at every radius.  The
divergence CURV-1 was looking for exists, but not at the centre: the presented excess has logarithmic singularities at
interior radii where two presented directions become antipodal, on the seat layer at exactly r* = r_s·Δ/(1 − γ₁₂²),
below both null radii.  The horizon is invisible to the excess (finite on both sides, sum rule intact); its cost is the
fixed rotor quarter-turn, which only Path 2 sees.

## What was run
Inputs, all re-derived in-suite: CARRY-1 G′(r) = G + f kkᵀ with f = r_s/(r − r_s), det G′ = Δ·r/(r − r_s); KIN-2a pinning;
CONT-1 interior sheet; RULING-2 layers A (root) and C (a₁+a₂+a₃); the area-tangent identity; CARRY-7 null radii;
CARRY-8 rod-free (frozen).  Grams N1 (¼, −⅓, ⅕), N2 near-coplanar (9/10, 9/10, 5/8), N3 near-orthogonal, N4 negative-sum
(−9/20 ×3); control (9/10, 9/10, 3/5) rejected, Δ = −1/125.  Symbolic theorems in the general Gram where they close.

Two functionals on the presented Gram, four pole-orbit members each (the four triangles cut out by three presented lines):
- Ω_A raw: S = 1 + Σ e_ie_j G′_ij, V = √det G′ — obeys CARRY-8's letter (no division by a rod), not the excess of any triangle.
- Ω_B normalised: S = 1 + Σ e_ie_j Γ′_ij, Γ′ = G′_ij/(√ρ_i√ρ_j), V̂ = √det G′/∏√ρ_i — the excess of the presented directions,
  degree 0 in every presented length; divides by rods, so not rod-free under CARRY-8's letter.
Ω = −i log W, W = (S+iV)/(S−iV), principal branches; the multiset over the four members is branch-free (a branch flip of
√ρ_i is e_i → −e_i).  Classification: IDENTITY (= baseline at r → ∞) / FINITE-NONTRIVIAL / UNBOUNDED / NO-LIMIT, plus
the half-angle sign where B is real.

## Path 1 — the seat-triangle excess (40/40)

| point | layer A, Ω_B (seat) | layer A, Ω_A (raw) | layer C, Ω_B (hole) | layer C, Ω_A (raw) |
|---|---|---|---|---|
| r → ∞ | frame's own orbit | same | same | same |
| r → r_s± | finite, real, sum 2π; all-agreeing member → 0, mixed members finite non-zero (lune limits, 0/0 resolved along the family), two sides equal | 0 | same shape | 0 |
| r → 0⁺ | **{φ₃, π−φ₃, π−φ₃, φ₃}** real, exact at N1–N4, theorem T1 for every Gram | **0**, B a positive rational, every Gram (T1f) | **{2π, 0, 0, 0}**: (+++) member's B real negative, half-angle **−1**; three members 0 (T2, hemisphere) | 0; sign of B is a frame sign bit (T2c), vanishing at the orthogonal frame |

Nothing is UNBOUNDED or NO-LIMIT at r → 0⁺ on any (layer, functional, member, Gram).  P1-7e is the check that would
have printed the diverging member; it printed none.

Theorems (general symbolic Gram, in-suite):
- T1: on layer A, Γ′₁₃, Γ′₂₃ → 0 like √r, Γ′₁₂ → cos φ₃, V̂ → sin φ₃; (γ₁₂ − γ₁₃γ₂₃)² + Δ = (1−γ₁₃²)(1−γ₂₃²) puts the centre
  value on the unit circle; 2 arg(1 + cos φ + i sin φ) = φ.  So Ω_B(+++)(0) = φ₃ exactly.  Picture: the root's presented vector
  vanishes, the two visible axes project to the plane perpendicular to the root; the presented triangle is the polar one
  whose excess is its apex angle.
- T1f: on layer A, S_A(0) = 1 + γ₁₂ − γ₁₃γ₂₃ = 1 + cos φ₃ sinθ₁₃ sinθ₂₃ > 0 for every admissible Gram, so the raw functional
  reads 0 at the centre with no frame content.
- T2: on layer C, G′(0)·(1,1,1)ᵀ = 0 — the presented vectors sum to zero (projections of the axes onto the plane ⊥ K, with
  K ∝ their sum).  Three plane vectors summing to zero surround the origin, so S < 0 and V = 0: the (+++) triangle is the
  hemisphere, Ω = 2π, B negative (argument + exact at N1–N4).
- T3: on layer A the two visible presented directions collide (G′₁₂² = ρ₁ρ₂) at r* = r_s Δ/(1 − γ₁₂²), with r*/r₁ =
  sin²(dihedral at a₁) < 1, so below both null radii, where both visible rods are positive: a real collision inside.

Interior structure (exact radii, both layers):
- Null radii r_i = r_s(1 − k_i²) (H-22): the orbit multiset is continuous through each one; individual labels swap (P1-8b).
- Presented collisions: layer A one, at r* = r_sΔ/(1−γ₁₂²) (N1: 2711/3375); layer C three, one per pair.  At each, the members
  for which the pair is antipodal have S + iV → 0: Im Ω → ∞, a logarithmic singularity of the excess with no rod in any
  denominator (P1-8c).  The raw functional has its own exact rational interior zeros of S_A² + det G′ (P1-13a), where
  Im Ω_A → ∞ — divergences of a functional with no geometric content.
- Reality by band: exterior real; just inside the horizon real (all rods negative, all cosines real); between null radii
  complex; on layer C below all three null radii W is real, so Ω is purely imaginary or π + imaginary; on layer A the
  root's rod is negative on the whole interior, so Ω is complex there and only becomes real at the centre.

Cross-checks: route V (rotor sandwich in the Pauli representation at the pinned rapidity, +i sheet, never forming G′) =
route G to 5.6e-40 on both sheets; l'Huilier agrees with Ω_B in the exterior to 3.7e-40 (mod 2π); the exact orbit identity
∏_ε(S_ε + iV) = **−4**(1−Γ₁₂²)(1−Γ₁₃²)(1−Γ₂₃²) holds identically in r on both sheets (the spec had the sign wrong; fixed);
the four exterior excesses sum to 2π at 40 digits; Ω_B is invariant under v_i → c_iv_i and Ω_A is not; the classifier
proves it can return UNBOUNDED (det G′ at the horizon, exponent −1), IDENTITY, and a rejection.

Conditionality: KIN-2a (pinning), CONT-1 (continuation), RULING-2 (layers) — all declared; the branch convention (principal
square roots, +i sheet) — declared, with the mirror sheet a sign map; CARRY-8's letter for the raw functional.

## Path 2 — loop holonomy from the relative rotors (P2-1..P2-9, all pass)
- The transport built from part 2's relative rotor A_b R A_a⁻¹ is **flat**: every contractible loop of seats has holonomy
  exactly 1 on both sheets, at every size and base rapidity (exact: A(−L,K)A(L,K) = 1 and R(−e)R(e) = 1; rectangles and a
  5-seat polygon at three rational points, 40 digits).  The curvature density (H − 1)/(d·e) is 0 identically, including at
  the centre: **IDENTITY**.  There is no local curvature reading in the model's own variables along the family.
- Any closed polygon telescopes to A₀(∏R)A₀⁻¹: only the product of the geodesic rotations survives.
- A ring of seats **around the centre**, at any interior rapidity (0.5, 0.01, 2.25; N = 3, 4, 6) and on the exterior sheet, has
  holonomy **−1 exactly** in SU(2): identity in SO(3), spin sign −1.  Exact: ∏R = R(2π/N)^N = R(2π) = −1.  Shrinking the ring to
  the centre leaves it −1: **FINITE-NONTRIVIAL, the deck ℤ₂**, r-independent — the centre is where the pivot field's full
  turn can no longer be contracted away.  The frame returns to itself (conjugation by −1 is trivial); only a spinorial datum
  flips sign — the click-monodromy paper's Q₈ statement in dynamical form.
- Failure branch: with a twist of π/5 about the pivot at each transport step the ring holonomy is not −1 (scalar part
  −0.309).  The −1 is conditional on TRANSPORT-1: geodesic transport, D-13's R.
- The fold: continuing r once around the centre returns μ → −μ (P2-8); the deck acts on the rotor as a real boost by −2μ
  along K, → 1 at the centre, and equals (I → −I)(K → −K) (D-8 re-verified).
- The horizon: interior rotor = exterior rotor × c, c = (1 + iK)/√2 (D-17 re-verified); a loop in and back out along the
  same K sees c·c⁻¹ = 1.  The crossing is a cost, not a holonomy — W-3 in Path 2's language.

## Path 3 — every frozen candidate reading at the centre (S2, S3, S4 pass)
| candidate (freeze) | layer A (seat) | layer C (hole) |
|---|---|---|
| S1 excess, Path 1 (08-30 identity) | finite, φ₃ | finite, 2π, spin −1 |
| S2 Gram-submersion ‖F‖ ~ Δ_pres^(−m/2) (08-30) | m = 0, finite: Δ_pres(0) = Δ/((1−γ₁₃²)(1−γ₂₃²)) | m = 1, **unbounded ~ r^(−1/2)**: Δ_pres has a simple zero |
| S3 inertia node detector det𝕀 (08-30) | 8 − 2cos²φ₃ ≠ 0, not a node | ≠ 0, coplanar not rank-one; → 0 at the horizon (H-15 re-found) |
| S4 Kummer radicands 1 + Γ′ (08-28) | {1 + cos φ₃, 1, 1} finite, positive | finite; at r* the antipodal member's radicand → 0 exactly |
| Path 2 loop holonomy | 1 contractible, −1 encircling | same |

Agreement: every candidate that is rod-free under CARRY-8's letter is finite at the centre on both layers.  The one
unbounded candidate is S2 on the world layer.  Two honest readings of it, both recorded: (i) under CARRY-8, Δ_pres divides
by presented lengths, so the divergence "has a rod in the denominator" (the same √r as the fold, E-1); (ii) geometrically,
Δ_pres is the normalised presented state reaching the coplanar boundary of the elliptope, and ‖F‖ ~ Δ^(−1/2) is that
boundary's own curvature blowing up — a genuine event of the presented directions, not bookkeeping.  Which reading
governs is the CARRY-8 ruling (below), not a computation.

## Predictions scorecard (filed in PREDICTIONS.md before the runs)
| | prediction | outcome |
|---|---|---|
| W-1 | Ω imaginary inside because det G′ < 0 | **Partly.** True on the world layer below all three null radii (W real: Ω ∈ iℝ or π + iℝ) and for the raw functional everywhere inside; false on the seat layer, where the root's rod is negative on the whole interior and Ω is complex until the centre. |
| W-2 | centre holonomy finite, non-trivial, the deck ℤ₂ | **Yes.** Path 2: −1 in SU(2) for every loop around the centre, r-independent. Path 1, world layer: the hemisphere, spin −1. |
| W-3 | horizon crossing a fixed quarter turn | **Yes** (Path 2, D-17 re-verified). Path 1 cannot see it: the excess is finite and two-sided at r_s. |
| C-1 | reality splits by layer | **Yes**, as stated. |
| C-2 | seat-layer centre = dihedral angle at the root, {φ₃, φ₃, π−φ₃, π−φ₃} | **Yes**, exact, every Gram (T1). |
| C-3 | world-layer centre = hemisphere {2π, 0, 0, 0}, spin −1 | **Yes** (T2). |
| C-4 | horizon: no value predicted | mixed members finite non-zero, sum 2π; the excess does not see the crossing. |
| C-5 | visible pair collides at r_sΔ/sin²θ₁₂ | **Yes** (T3), and it is where the excess diverges. |

## Two things this run decides, and two it hands back
Decided:
1. CURV-1 fork **(b)**: the model's rod-free curvature readings are all finite at the centre.  The centre is a finite event
   with two addresses — an angle on the seat layer, a spin sign on the world layer — not a divergence.
2. The excess *does* diverge, at interior radii r* where presented directions become antipodal (exact formula on the seat
   layer).  These are not rods in denominators.  They are a new interior structure: a second node-type degeneration of the
   presented state between the centre and the null radii.  No standard analogue is claimed; lineage unsearched.

Handed back (rulings, not computations):
- **CARRY-8.**  Read literally, the frozen definition selects Ω_A, which reads 0 at the centre for every admissible frame on
  the seat layer and a sign bit on the world layer; read as "angular separations are allowed", it selects Ω_B, which reads
  the dihedral angle and the hemisphere.  The same ruling decides whether S2's world-layer divergence counts.
- **Identity semantics.**  The tables use the baseline-relative label (twin spec) with the half-angle sign as a column.
  Under 0-mod-2π semantics the world-layer centre is "identity with spin −1" instead of "finite-nontrivial, 2π"; same content.

## Comparison stage (names spoken here only)
- Ω_B is the spherical excess = the integrated curvature enclosed by the presented triangle (Gauss–Bonnet); Path 2's ring is
  the parallel-transport holonomy of a 2π turn of the pivot direction, i.e. a hemisphere's solid angle: 2π rotation, −1 spin.
  The two paths agree because they are the same theorem in two languages.
- A flat connection with ℤ₂ monodromy at a point is a topological defect (a vortex of the pivot field), not a curvature
  singularity; the Berry/Pancharatnam sign of the 2π turn is the deck.
- The seat-layer centre reading φ₃ is the apex angle of the polar triangle: the classical "excess of a triangle with two
  right angles equals its apex angle".
- The world-layer S2 divergence is the elliptope's boundary curvature (O'Neill), reached because the presented state is
  coplanar at the centre; standard tidal divergence (~1/r³) is NOT reproduced and was not computed — the model has no
  acceleration tier.
- The interior antipodal collisions and the ordering r* < r_i < r_s have no analogue named here.

## Files
`curv1_path1.py`, `run_curv1_p1.log`, `curv1_path1_results.json`; `curv1_path23.py`, `run_curv1_p23.log`;
`docs/prereg/CURV-1/{README, PATH1-SPEC, PREDICTIONS, CURV1-RESULTS}.md`.  The recovered twin spec's protocol (runner
isolation, harness self-tests, rod audit, both functionals, N1–N4) and this session's object (the branch-free orbit,
W = e^{iΩ}, route V, collisions) were merged in the suite itself; the separate-runner discipline was dropped by Will's
instruction on 2026-09-02 evening in favour of exact symbolic verification.

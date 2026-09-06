# MATH FOUNDATION AUDIT — every body of mathematics the seated root uses, and what it is doing there
Written 2026-09-05 against HEAD `45c7443` (paper v0.5, PRIMITIVES v0 draft, LABELLED-MODEL, CONJECTURE-COSMOLOGY).
**Status: AUDIT.** Nothing here changes a tier in the repo. It says which tiers are carrying the wrong label, which body of
mathematics is doing which job, where a stronger and more general frame already exists, and where the model would have to
build its own. Every runner cited was read in full or in its header block; `suites/prim_t7_seat_form.py` (18/18) and
`suites/prim_t8c_reciprocity.py` (14/14) were re-run during the audit and are green in this environment. Will's untracked
curvature notes (`GRAM_SUBMERSION_CURVATURE.md`, `INERTIA_NODE_DETECTOR.md`) were read from the main checkout. No other
corpus was used as an input.

**RULING (Will, 2026-09-05), applied throughout this audit:** nothing is banned. A derivation lists its hypotheses; anything it
uses that is not listed must be derived within it or declared with its cargo. Every "banned" in the repo's docs and runner
headers was shorthand for "not among the inputs" and has been corrected to that wording in the same session as this audit.

The question put to each body: (1) what job is it doing; (2) is it the best fit for that job; (3) is it conditioned to its
answer, by the repo's own test — *could the physics have failed to come out?*; (4) what is its first-principles home, and does
that home generalise the construction; (5) what does it look like with every physics label removed.

---

## 0. Verdict in one paragraph

The model's mathematics is sound where it computes and mislabelled where it narrates. Almost every "DERIVED" result in the
primitives layer is a correct theorem of standard mathematics — real forms of a complex Lie algebra, the Gram matrix of a
Lorentzian frame, Gauss–Bonnet on H², Cl(2,1) = M₂(ℝ) ⊕ M₂(ℝ), the Hodge star on 2-forms — and none of them could have failed,
because the declaration that makes them true is P3/P6/P11 (which planes are hyperbolic), not Cl(3). Removing Cl(3) from the inputs moved the
Lorentzian hypothesis from an algebra into three sentences; it did not remove it. That is fine, and the primitives are a *better*
place for it than Cl(3) was, but "reached, not declared" is the wrong label for T1, T2, T7a and T8c, and the ledger should say
so. The single most consequential finding is structural, not tiering: **the repository carries two state spaces that are not
the same object** — the paper's (S²)³ Euclidean Gram cell with its Cayley cubic, Kummer module, Bargmann phase, VIEW-1 and
PREDICTION-1, and the primitives' Lorentzian frame space with signature (2,1), no real nodes, and a different branch locus.
Everything the paper calls its own (§2, §9, PRED-1, CURV-1) lives on the first; everything the rebuild calls derived lives on
the second; nothing in the repo yet connects them, and they are related as *two real forms of one complex object*, which is
also what the "imaginary angle" and the Wick loops are. On the constructive side: the model can be stripped of every physics
label and stated as a small general framework (§5) — a frame of *n* lines in a space with a form of signature (p,q), a seat
as a line, a presentation as the projection off it, the pivot group as O(p,q), the sheet as the orientation ℤ₂, the horizon
and branch locus as rank-drop strata of principal minors of the Gram matrix — and three of its bodies have a first-principles
home that is strictly stronger than the version in the repo (real forms, §2.1; the Poisson kernel and boundary
representations for the whole readout tier, §2.9; determinantal strata and the vanishing-order calculus for every
singularity, §2.11). Original mathematics is warranted in exactly one place (§6): a uniform valuation calculus for readouts
at determinantal strata across real forms — the thing the trichotomy, the blind-mass theorem, the inertia exponent law and the
Kummer descent are each one corner of. It is not needed anywhere else, and the dynamics tier in particular should not be
invented: Will's impedance is a Hermitian structure on the seat's state space, and the Kähler dial already is one.

---

## 1. Inventory — what mathematics is in the repo, where, and what it is doing

| # | body of mathematics | where it lives | job in the model | standard name for what is actually computed |
|---|---|---|---|---|
| 1 | complexified rotations, real forms | `prim_t1_t3` | the pivot group; "imaginary angle" | real forms of SO(2,ℂ) and so(3,ℂ): SO(1,1), SO(2,1) vs SO(3) |
| 2 | Gram matrices of an indefinite form; gluing; tilt | `prim_t7*` | the seat's form, its state space, the lapse, degeneracy sets | linear algebra of a frame in ℝ^{2,1}; ADM lapse in linear form |
| 3 | hyperbolic geometry, Gauss–Bonnet | `prim_t5` | Wigner = area | holonomy of the Levi-Civita connection on H² |
| 4 | complexified one-parameter groups, Sym² | `prim_t4`, `thm_d2_unruh` | the 2π / 4π periods, "the cover reached" | SO(2,1)'s vector rep = Sym² of SL(2,ℝ)'s fundamental |
| 5 | real Clifford algebra, Pin(2,1), central idempotents | `prim_t5b`, `prim_t5c`, `prim_t8a`, `prim_t8b` | the sheet, the fourth direction, the spinor bilinear | Cl(2,1) ≅ M₂(ℝ) ⊕ M₂(ℝ); Cl(2,1) ⊂ Cl(3,1) or Cl(2,2); the charge-conjugation form |
| 6 | exterior algebra, Hodge star | `prim_t8c`, `thm_k_response_map` | signature of Γ₄; constitutive map | ⋆² = (−1)^{n₋} on 2-forms in n = 4; premetric electrodynamics (Hehl–Obukhov) |
| 7 | finite combinatorics | `label2`, `label3`, `label3b`, `prim_t5b` part B | the octahedron, the Boolean lattice, parity | the 3-cross-polytope; 2^{C,H,G}; S₃/A₃; hyperoctahedral group B₃ |
| 8 | real algebraic geometry of symmetroids, convex geometry | `suites/cayley.py`, paper §2 | the branch surface, the click group | Cayley's nodal cubic; the elliptope; Aut = S₄ |
| 9 | Kummer theory, Galois covers, monodromy | `suites/galois.py`, `thm_b_monodromy`, `thm_galois_deck_descent`, `correspondence_july` | the deck ℤ₂, the descent of arg B | multiquadratic extensions; odd-valuation rank criterion; inertia at nodes |
| 10 | U(1) graph cohomology, Bargmann invariants | `suites/bargmann.py`, `suites/graph_cocycle.py`, `pred1_*`, `pivot_map*`, `census_c`, `view1` | the phase, PREDICTION-1 | Bargmann/Pancharatnam; Li–Wagner–Zhang Thm V.1; Van Oosterom–Strackee |
| 11 | Möbius action on the boundary, Radon–Nikodym cocycles, power means, Mellin | `thm_g`, `thm_rn`, `bridge_dbp`, `graph_cocycle` part 1, spec I | the readout tier: presentation families, blind mass, mirror, trichotomy | Poisson kernel of H^n; conformal weights; Hellinger midpoint |
| 12 | Kähler geometry of (S²)³; Hamiltonian vs gradient flow | `suites/kahler.py`, `suites/doors.py` | the dynamics candidates; the dial | h = g + iω on (ℂP¹)³; metriplectic flow |
| 13 | Riemannian submersion, O'Neill, mechanical connection | Will's `GRAM_SUBMERSION_CURVATURE`, `INERTIA_NODE_DETECTOR`, THM-J | the curvature tier KIN-2b; the dead dynamics | shape-space geometry (Littlejohn–Reinsch); ‖F‖ ∼ Δ^{−m/2} |
| 14 | Cl(3) paravector boosts; PG river; Gauss–Codazzi on a flat slice; geodesics | `thm_h*`, `thm_i*`, `thm_l`, `thm_m`, `thm_n`, `thm_o`, `curv1_*`, `sect1` | the gravitational sector | Lorentz boosts of a 4-vector; ADM constraints in Painlevé–Gullstrand gauge; Schwarzschild/Lense–Thirring |
| 15 | thermodynamic declarations | `thm_k_clausius` | the pinning from temperature | Jacobson/Verlinde equipartition with the Bekenstein–Hawking entropy declared |
| 16 | independent-oscillator elimination, FDT | `debt2`, `debt2b` | the seat's dissipator | Ford–Kac–Mazur / Caldeira–Leggett |
| 17 | impedance / circuit reading | PRIMITIVES "IMPEDANCE", `hunch_z0` | the shape of the dynamics tier | linear response with a complex (dissipative + reactive) coefficient |

Rows 1–7 are the primitives layer (post-09-04). Rows 8–13 are the paper's bare and readout tiers (08-28 to 09-04), all built on
the Euclidean state space. Rows 14–16 are the gravitational and thermal sectors, all built on Cl(3). Row 17 is a hunch.

---

## 2. Body-by-body

### 2.1 Real forms — the pivot group (row 1; T1–T3)

**What it computes.** T1: conjugating a rotation R(iλ) by diag(1, i) gives a real boost; the only invariant symmetric form is
diag(1, −1). T2: making one coordinate of ℝ³ imaginary turns two of the three so(3) generators symmetric and leaves one
antisymmetric; the realified algebra preserves diag(1,1,−1); over all eight realifications the number of boost planes is 0 or 2.
T3: the Wigner closed form in SO(2,1). All exact, all correct.

**What it is.** T1 and T2 are the classification of real forms of SO(2,ℂ) and so(3,ℂ). A complex simple Lie algebra has finitely
many real forms; for so(3,ℂ) they are so(3) (compact) and so(2,1) (split), and a real form realised on ℝ³ carries a form of
signature (3,0) or (2,1) — nothing else exists, which is T2f's {0, 2}. "Rotation by an imaginary angle" is the Weyl unitary
trick run backwards: multiplying a generator by *i* moves between the compact and the split real form inside the complex group.
The abstract algebra {L, iL} ≅ so(1,3) of T2g is so(3,ℂ) as a real Lie algebra (correction #1 in the handoff already says this).

**Could it have failed?** No. Once P3 says a plane is hyperbolic and P6 says a pivot is an imaginary-angle rotation, the
signature and the group are fixed by the classification. The runner's T1c ("the real pivot does NOT preserve Q_E: signature
change is forced") is not a discovery about the model; it is the statement that SO(1,1) ≠ SO(2). By the repo's own test, T1,
T2 and T7a are constructions scoped to their answer.

**Best fit?** As a *computation*, yes — cheap, exact, and it produces the one non-trivial fact (T2f: no (1,3) on three lines).
As a *foundation*, the honest statement is stronger and shorter than the runner's: *the model declares a real form.* P3 + P6
is exactly the declaration "the frame's form has signature (2,1)". That is a better declaration than Cl(3), because it
carries no representation, no spinors and no fourth direction with it — the rebuild was worth doing — but it is a declaration,
and the ledger entry should read **signature (2,1): declared by P3/P11**, with T1–T2 as its consequences, not its derivation.

**First-principles home, and the generalisation.** Real forms and Cartan decompositions. For *n* lines through a root the same
primitives give so(n,ℂ) with real forms so(p,q), p+q = n, and "one imaginary coordinate" gives so(n−1,1). Nothing in T1–T3 is
specific to three. The label-free statement of P6 is: *the pivot group is the identity component of O(q) for the frame's form q,
and the complexified group O(q_ℂ) contains every real form; a pivot by an imaginary angle is a path in the complex group between
real forms.* This is also what the Wick loops of LABEL-1 are (§2.7), so one sentence covers both.

### 2.2 The seat's form — gluing, seat coordinates, tilt, degeneracy sets (row 2; T7a–g)

**T7a is a definition, not a theorem.** The six block equations assign the six entries of a symmetric 3×3 matrix directly
(each diagonal entry twice, consistently); "unique glue" is the statement that a symmetric 3×3 matrix is determined by its
entries. What T7a actually says is: *a frame consisting of one negative line and two positive lines has Gram diagonal (−1, +1,
+1).* That is P11 written as a matrix. T7c11 ("SO(2,1) REACHED, not started from") is the same fact read at the orthogonal
state. The tier should be **[declared: P11] → Gram diagonal (−1,+1,+1)**.

**T7c–T7g are real theorems, and they are the model's best work in the primitives layer.** The seat coordinates (t, l₁, l₂),
det G = −cosh²l₁ cosh²l₂ sin²t, the branch locus {sin t = 0} as a surface with no nodes (T7c10 — see §4 for what this does
to the Cayley cubic), the tilt η = vᵀS⁻¹v with sech λ the lapse, the two degeneracy sets D_plane ⊋ B_form, the bounded lapse,
the horizon off the branch locus, the region table, the two horizons ħ ∓ G. These are the linear algebra of a Lorentzian frame,
correctly done, and they are where the model says something *about seats* rather than about signature.

**Their first-principles home: determinantal strata of the Gram matrix.** Every object in T7c–g is a rank condition on a
principal minor of G:

| model name | condition | which minor |
|---|---|---|
| branch locus B_form | det G = 0 | the full 3×3 |
| D_plane (rulers' span degenerate; the horizon η → ∞ lives here) | det S = 1 − γ² = 0 | the (ħ,G) 2×2 |
| the forbidden zone | signature (1,2) | sign pattern of the minors |
| the interior sector | det S < 0, det G < 0 | same minors, other signs |
| CURV-1's interior collisions | Γ′₁₂² = ρ₁ρ₂ | a presented 2×2 |

So the horizon, the branch locus and the interior antipodal collisions are three strata of one stratification — the space of
symmetric matrices with a fixed diagonal, stratified by the ranks and signs of its principal minors. This is general (any *n*,
any signature), it is label-free, and it is the setting in which "the singularity is a locus with computable structure"
(paper §1) is literally true: the singular loci are algebraic subvarieties with a vanishing order, and §2.11 says what the
vanishing order buys. Recommended as the label-free definition of *horizon* in §5.

**One consequence already available.** In T7f the interior is a *real* sector (|γ| > 1, rulers' span Lorentzian) and
N² = (1 − γ²)/(−det G) passes through zero continuously at |γ| = 1 and is negative inside. That is CONT-1 (λ = μ + iπ/2,
lapse imaginary) as a real path through the state space rather than a declared analytic continuation. If the pinning's fibre
(two free functions of r, T7f5) can be chosen so that η(r) reaches D_plane and the path continues into |γ| > 1, CONT-1 moves
from **[declared]** to **[derived | T7f]**. One runner. Not done here.

### 2.3 Hyperbolic geometry — Wigner is area (row 3; T5 part A)

Correct, and classical: the rotation produced by composing two boosts is the holonomy of the Levi-Civita connection on
H² = SO(2,1)/SO(2), and by Gauss–Bonnet holonomy equals enclosed area at curvature −1. The runner's Lorentzian
Van Oosterom–Strackee formula is the right closed form. **[retrodiction]** as mathematics; nothing conditioned. It generalises
verbatim to H^n (the Thomas–Wigner rotation is the area of the geodesic triangle for any n). Its role in the model is to identify
the seat-cycle holonomy of §9 with the boost-triangle area — which is a *Lorentzian* statement on the (2,1) frame space, and
therefore the natural replacement for the Euclidean solid-angle statement of `suites/bargmann.py` when §9 is migrated (§4).

### 2.4 Complexified one-parameter groups — the periods and the cover (row 4; T4, THM-D2)

**What holds.** exp(iθK) on the vector rep has period 2π and equals diag(−1,1,−1) at π; on the 2-dim rep it has period 4π;
the vector rep is Sym² of the 2-dim rep, so the −1 squares away. All standard: SO(2,1)'s vector representation is the adjoint,
the symmetric square of SL(2,ℝ)'s fundamental. "The cover is REACHED as Sym²" is exact. The tensoriality theorem is
correctly retired.

**What is conditioned.** T4c–T4g take the pinning tanh λ = √(r_s/r) as input. The pinning is derived in THM-K from SCREEN-1,
which declares S = k_B A/(4ℓ_P²). The check T4f3 — "the first law with S = A/4 returns T_from_vector and NOT T_from_cover" —
is then a consistency check between two *declared* conventions (Bekenstein's ¼ and the KMS circle), not an independent
selection of 2π: Bekenstein's ¼ was fixed by Hawking's 2π. The genuine content of the selection is elsewhere and is structural:
**P4 makes the state a set of angles between lines, so a seat reads the vector (projective) representation, and the vector
representation's period is 2π.** That sentence is model-internal and label-free, and it is the one that should carry the
result. The first-law check should be re-labelled **[consistency of two declared conventions]**.

**Label-free form.** The Wick face is the compact real form's orbit inside the complex group; its period on the representation
the seat reads is the period of that representation. The thermal reading of the period (LBL-1: circumference ħ/k_BT) is a
label and the paper says so. Generalises: for any real form and any representation the period is a representation-theoretic
invariant, and "which representation a seat reads" is a primitive-level question (P4/P5).

### 2.5 Real Clifford algebra — the sheet and the fourth direction (row 5; T5b, T5c, T8a, T8b)

**Correct throughout, and honest about where it stops.** Cl(2,1) ≅ M₂(ℝ) ⊕ M₂(ℝ) with central pseudoscalar I, I² = +1; the
two blocks are the two central idempotents (1 ± I)/2; the raw spinor trace tr(Γ_CΓ_HΓ_G) = 2D is the pseudoscalar part of a
product of three vectors, i.e. the determinant; the flips are reflections, their Pin(2,1) lifts anticommute at the orthogonal
state, so the sheet is the parity in S₃/A₃; the anticommutant of the three generators in M₄(ℝ) is a 2-parameter family;
adjoining Γ₄ makes the algebra simple; the per-block spinor bilinear is forced symplectic; the relative sheet sign *s* is not
fixed by Cl(2,1). T8b's a6 ("my prediction was WRONG: the blocks are inequivalent irreps") is exactly right and is recorded.

**What it is, label-free.** The sheet is the orientation ℤ₂ of the frame, read through the spinor representation because
the vector representation is orientation-blind (T5c 2f: the projector Bargmann invariant is even in D). The two blocks are the
two connected components of the space of oriented frames; Γ₄ is a generator of the larger algebra whose conjugation is the grade
involution — the *outer* automorphism of Cl(2,1) made inner. This generalises: for any (p,q) with p − q ≡ 3 mod 4, Cl(p,q) is
a direct sum of two simple algebras (I central, I² = +1) and the same story holds; for the other residues the sheet is not a
central idempotent and the model's sheet mechanism would need restating. That is a real constraint on generalisation and
should be recorded: **the sheet as a block label exists only in signatures with p − q ≡ 3 (mod 4)**, of which (2,1) and
(3,0) are the smallest.

**The dependency worth naming.** Everything here is standard Clifford theory *given* the (2,1) form, i.e. given P3/P11. The
readout's role is the one honest input: which block a one-sided reading selects, and the sign *s* of T8b. Neither is
conditioned; both are declared, and labelled as such.

### 2.6 Hodge star and reciprocity — the signature of Γ₄ (row 6; T8c, THM-K(a))

**The computation is right; the tier omits a hypothesis the argument uses.** ⋆² = (−1)^{n₋} on 2-forms in
four dimensions; exp(θ⋆) is compact iff ⋆² = −1; Cl(3,1) has n₋ = 1. All exact. But the premise — "Maxwell's electric-magnetic
reciprocity is a compact rotation, therefore n₋ is odd" — is a hypothesis about a named physical theory, and it is not
listed among the runner's inputs. It is also not an inference: "duality on 2-forms is compact" and "n₋ is odd" are the
*same* statement, so the argument selects Lorentzian signature by asserting it in another vocabulary. The link that would make
it model-internal — that 2-form duality acts on the seat's compact ruler plane as a real rotation — is admitted in P13 to be
"argued, not computed". Tier should be **[derived | Maxwell duality compact] → Γ₄² = +1**: conditional on a declared input, not
unconditionally [derived]. Nothing forbids the hypothesis; the ledger has to carry it.

**A model-internal route exists and is already in the runners.** T8b proves signature(Γ₄) = (adjoint class) × s. T8c's c6′
then uses P5 ("one seated line = one clock") to set s = +1 *after* fixing the signature by Maxwell. Run the inference the other
way: P5 ⇒ s = +1 (Will's "time comes as a pair; one clock"); declare that the arrival operation is an *observable* (self-adjoint
under the seat's own bilinear — a P12-type principle, not a physical theory); then Γ₄² = +1 follows. That replaces a hypothesis about a named
theory with a model principle and moves Maxwell to the comparison stage. Recommended (§7).

**THM-K(a)** is honest and general: the Hodge star of q is the unique q-built constitutive map up to a scalar, which is
premetric electrodynamics; the cone clause follows; the scalar is free. Label-free, generalises to any signature and dimension.

### 2.7 Finite combinatorics — the octahedron and the lattice (row 7; LABEL-2, LABEL-3, T5b part B)

**Mathematically trivial, correctly labelled as such, and it generalises immediately.** Six poles, twelve edges, eight faces
is the 3-cross-polytope; the faces are the sign vectors {±1}³; adjacency = one flip; the six shortest paths are S₃; the sheet
is S₃/A₃. For *n* lines: 2n poles, 2n(n−1) edges, 2ⁿ faces (orthants), the Boolean lattice 2ⁿ, the sheet S_n/A_n — the model's
labelling scheme is the n = 3 case of the cross-polytope, and nothing about it is three-specific. The paper's relabelling group
B₃ (order 48) is the hyperoctahedral group, the full symmetry group of this octahedron; the click group S₄ = B₃/±I is its
rotation group. So §2 of the paper and LABEL-2 are the same finite object seen from the state space and from the pole set.

**A structural coincidence, recorded as a hunch and no more.** The 6 poles ±eᵢ are the short roots of the B₃ root system;
the 12 edges, read as the sums ±eᵢ ± eⱼ of their two poles, are the 12 long roots; the 8 faces (±1,±1,±1) are, up to the
factor ½, the weights of the 8-dimensional spin representation of Spin(7). The counts are forced by the combinatorics and are
therefore not evidence of anything; the repo's own rule on triality applies ("dim 8 is a coincidence until a construction
exists"). Worth exactly one runner asking whether the sheet parity of T5b is the Spin(7) chirality; not worth more.

**The Wick loops (LABEL-1)** are paths in the complex group SO(3,ℂ) through its real forms — quarter-turns by imaginary angle
in a hyperbolic plane land on the antipode. The physics table they are checked against is declared, and the runner says so.
Label-free: the loop structure is the action of the complex group on the pole set.

### 2.8 The Euclidean state space — Cayley, elliptope, Kummer, monodromy, Bargmann, VIEW-1, PRED-1 (rows 8–10)

**All correct as mathematics; all built on an object the rebuild no longer owns.** The paper's state is three unit vectors in
S² with a positive-definite Gram, Δ = det G ≥ 0, the elliptope, the Cayley cubic with four real nodes, the click group S₄ as
its automorphisms, the multiquadratic field with four independent radicands, the descent B ∈ F(i, √Δ), the deck as the inertia
of the smooth cubic, the geodesic triangle's solid angle as 2 arg B, the forced pivot circle, the pivoted observation map, the
ordered-register operation reading sgn V. Under the primitives the frame carries one (2,1) form (T7b1), has one negative line,
and its Gram symmetroid has *no real nodes* (T7c10). See §4 for what survives migration and what does not. Here: the
Kummer/Galois layer (row 9) is over the field ℚ(γ), i.e. it is a statement about the *complex* surface and survives the change
of real form intact; the elliptope and the real nodes (row 8) do not; the Bargmann/Pancharatnam layer (row 10) must be redone
with the Lorentzian triangle (T5 part A already contains the n = 3 case). Row 10 is a retrodiction (Li–Wagner–Zhang) and the
repo says so; nothing in rows 8–10 is conditioned to a physical answer, but PRED-1 is conditioned to a state space that is not
the model's.

### 2.9 The readout tier — cocycles, blind mass, power means, the trichotomy (row 11)

**This is the one body where a stronger first-principles frame exists, is already verified in the repo's own checks, and
has not been named.** The facts the repo has proved: the c-seat reading family is s(u) = cosh λ + sinh λ·u with the forced
measure du/2 on the sphere of directions; the pivot acts on directions by aberration and du′/du = s^{−2} (RN-1, Q = 2); the
ħ-seat has q(φ) = cosh 2λ + sinh 2λ cos 2φ with dφ′/dφ = q^{−1} (RN-2, Q = 1); the blind moment is the total mass of the
source measure; the mirror I(p,m) = I(−p, 1−m); the moments are Legendre polynomials P_{n−1}(cosh 2λ); the bath is a flat
band; DEBT-2b forces c(ω)²g₀(ω) ∝ ω².

**What these are.** In the ball model of hyperbolic n-space, the Poisson kernel at hyperbolic distance λ from the origin,
evaluated at a boundary point at angle θ, is P = (cosh λ − sinh λ cos θ)^{−1}, and the boundary Jacobian of the isometry that
moves the origin by λ is P^{\,n−1} in the round measure. So:

- the c-seat's family s(θ) is the reciprocal Poisson kernel of **H³** (boundary S², the view sphere), and Q = 2 = n − 1;
- the ħ-seat's family q(φ) is the reciprocal Poisson kernel of **H²** (boundary S¹, or ℝP¹ in the doubled angle), Q = 1;
- "the measure is the plane" (Archimedes: uniform in u = cos θ) is the round measure on ∂H³;
- the blind-mass theorem ∫ s^{Q} dμ′ = 1 is the statement that the Poisson kernel integrates to one (harmonic measure);
- the Legendre polynomials are the zonal spherical functions of SO(2,1);
- the mirror (p, m) ↦ (−p, 1−m) is the weight duality w ↔ Q − w of the boundary (principal-series) representations, the
  α ↔ 1−α symmetry `suites/graph_cocycle.py` already downgraded it to;
- Q-RN ("why is a reading a power of a Radon–Nikodym derivative?") has an answer: because a seat's readings are sections of
  homogeneous line bundles over the boundary of its symmetric space, and the weight is the representation label. Q is the
  carrier dimension — the 08-28 handoff's "structural bet", confirmed.

**Two consequences the repo has not drawn.** (i) The paper's Q_c = 2 presupposes a *sphere* of directions, i.e. 3+1. On the
primitives' own 2+1 frame the c-seat's directions form a circle and Q_c would be 1, like ħ's; the transport census (2,1,2) is
therefore a consequence of the fourth direction (T8), not of the bare tier — a consistency the repo can bank. (ii) The horizon
η → ∞ is λ → ∞, i.e. the state going to the boundary of the seat's H²; the seat's readout carrier *is* that boundary. The sky
(debt2's bath) and the horizon are the same boundary read two ways. This is the clean, label-free version of "the horizon is a
ruler turning into light" (P9 sharpened), and it is where Hawking's period (§2.4) lives without any thermal label: the
complexified boundary action's period.

**Best fit?** The current cocycle/Mellin formulation is correct and general (spec I of 08-28 is a good abstract statement of
it), but it is the representation theory of the pivot group on its boundary written without saying so. Naming it costs
nothing and buys the trichotomy (free / floored / pinned = the three ways a positive cocycle's zero set can sit relative to a
group orbit, a conjugation class, and a base map — the 08-28 "trichotomy classification" item) and the aggregation-exponent
question (a distinguished weight is a distinguished representation). **Nothing here is conditioned**: the reading family is
forced by isotropy, and the checks RN-1/RN-2 are exactly the Poisson-kernel Jacobian identities. Kill: if s^{−Q} is not the
boundary Jacobian of the pivot on ∂H^{Q+1} in the round measure, the identification is wrong — RN-1 and RN-2 already say it is.

### 2.10 The Kähler dial and the doors — dynamics candidates (row 12)

Correct on (S²)³ = (ℂP¹)³: J_a v = a × v, ω(u,v) = g(Ju,v), Hamiltonian flow = gradient flow rotated by J, the dial e^{Jθ}.
Standard Kähler geometry; nothing conditioned; the doors experiment was honest (RULING-1 is a ruling, and says so). Two things
follow. First, it is on the Euclidean state space and must be restated on the (2,1) frame space — where the natural Kähler
object is the seat's H² (a Hermitian symmetric space whose invariant Kähler structure is unique up to scale), not (S²)³.
Second, and more useful: **the dial is Will's impedance.** A complex coefficient Z = R + iX in V = IZ is a linear response
with an in-phase (dissipative) part and a quadrature (reactive) part, and that is exactly a Hermitian form h = g + iω on the
state space: g is the resistance, ω the reactance, the dial angle θ the impedance phase. So "the load's impedance as a function
of state" — named in PRIMITIVES as the one unknown behind the dynamics tier — is *a Hermitian structure on the state space*,
which the repo already has in `suites/kahler.py`. The open question is not what mathematics to invent but which Hermitian structure on
which quotient of the (2,1) frame space is forced; on H² there is one candidate up to a scalar. This is a hunch with a runner
behind it, and it is the reason §6 recommends against inventing a dynamics calculus.

### 2.11 The curvature tier — the submersion and the inertia tensor (row 13; Will's notes, THM-J)

**Correct, and it contains the model's most original general statement.** F = 2𝕀^{−1}(Σ Xᵢ × Yᵢ); ‖F‖ ∼ C(b)·Δ^{−m/2} with m
the local vanishing order of Δ; det 𝕀 = 8 − 2Σγ² − 2γ₁₂γ₁₃γ₂₃ = 0 iff rank-one; C(b) non-universal with a closed form. THM-J
then correctly showed the harmonic-map energy of this metric is not the field law (perihelion 1.5–1.87). All on the Euclidean
state space (round metric on (S²)³, elliptope base).

What is general and label-free in it is not the metric but the **exponent law**: the parity of the vanishing order flips the
sheet, the value of the vanishing order is the curvature exponent, and the same integer m is readable off a pointwise algebraic
invariant (det 𝕀) instead of an approach analysis. That is a statement about determinantal strata of a Gram symmetroid, and it
is the same mechanism as the monodromy dichotomy (smooth points flip √Δ, nodes do not), the Kummer valuation criterion (odd
valuation along a divisor), and the trichotomy (where a cocycle's zero sits). Four results, one calculus, not yet written as
one. This is the candidate for original mathematics in §6.

### 2.12 The gravitational sector — Cl(3), the river, Gauss–Codazzi, geodesics (row 14)

**Correctly labelled by the paper as bands; the audit adds only what the rebuild changes.** THM-H through THM-O use BARE-1
(a Cl(3) rotor acting on a paravector = a Lorentz boost of a 4-vector), the pinning, CONT-1, STAT-1/ROT-1 (Killing
conservation), PROP-1 (the river), NORMAL-1 (the Hamiltonian and momentum constraints), one boost per seat (flat rain rods),
K-6 linearity. Every physical number is Schwarzschild, PG, Lense–Thirring, and the paper's band/edge/position/surplus register
says so. By the repo's own test the sector could not have failed at first order; it *did* fail at second order in the spin
(THM-N, THM-O fork (b)), which is the sector's only evidential content and is honestly recorded.

**What the rebuild changes.** The paper's three axes a_c, a_ħ, a_G are three *spacelike* unit vectors (Bloch vectors) with
time supplied by the paravector scalar; the primitives' three lines are one timelike (c) and two spacelike. These are not the
same object in a different basis (§4). So the gravitational sector's "presented Gram" G′ = G + sinh²λ kkᵀ is a Euclidean Gram
of three spatial axes under a boost, and its interior continuation, its null radii and CURV-1's centre readings are all
statements about that object. They are consistent physics (they are PG's geometry) but they are not yet statements about the
primitives' frame. Migration route: the presented object in the (2,1) frame is the frame itself under a pivot — T7's
(t, l₁, l₂) with η(r) pinned — and the "hole" is a path in the state space crossing D_plane (§2.2); the two free tangential
functions of T7f5 are the sector's genuine debt, correctly named in the handoff.

### 2.13 The thermal declarations (row 15; THM-K)

SCREEN-1 declares S = k_B A/(4ℓ_P²), i.e. the Bekenstein–Hawking entropy with its ¼; THERM-1 equipartition; MASS-1 blind mass
as enclosed energy; EQ-1 Tolman equilibrium. Given Unruh's 2π these return r_s = 2GM/c² — this is Jacobson/Verlinde and the
paper cites them. The model's own content in THM-K is K-1: αN = (c²/2)|d(sech²λ)/dr| for any profile, which fixes the *exponent*
(and which T4e re-derives from λ alone). That identity is label-free (it is the lapse gradient of a static PG metric) and
general. The declarations are physics, correctly labelled, and cannot be stripped: they are the labels. The only tiering
correction is the one in §2.4 (the first-law check is a consistency of two declared conventions).

### 2.14 The bath (row 16) and the impedance reading (row 17)

Row 16: Ford–Kac–Mazur elimination and classical FDT, lineage loud, coupling derived in DEBT-2b from the Jacobian — honest.
Label-free it is: the seat's environment is the boundary of its symmetric space with the round measure pushed through the
pivot (§2.9), which is why the band is flat. Row 17: N² = 1/(1+η) as a divider ratio is algebra (any number in (0,1] is a
divider ratio); its content is the identification η = Z_source/Z_load, which is a claim about a response map not yet built.
Its mathematical form is §2.10's Hermitian structure. `suites/hunch_z0_impedance.py` is honest that it is dimensional analysis until the
charged load is constructed.

---

## 3. Register of conditioned mathematics — ranked by how much rests on it

By the repo's test ("could the physics have failed to come out?"). "Conditioned" is not "wrong"; several entries below are
the *right* place to put a declaration. The problem in each case is the label; in one case (row 3) the label omits a hypothesis.

| rank | item | what it presupposes | current label | audit label | what rests on it |
|---|---|---|---|---|---|
| 1 | P3 + P6 + P11 → (2,1) | the signature; a split real form | primitives (words) | **declared: signature** | T1, T2, T7a, T7c11, all of Cl(2,1) |
| 2 | the Euclidean state space (S²)³ | Cl(3)'s vector part; three spacelike axes | paper §2 "bare tier" | **[Cl(3) residue]** — not any seat's construction (T7b1) | Cayley, elliptope, S₄, VIEW-1, CENSUS-C, pivot_map, PRED-1's protocol, CURV-1, THM-H's G′ |
| 3 | T8c: Maxwell duality compact ⇒ Γ₄² = +1 | a named physical theory as arbiter; and it restates n₋ odd | [derived] | **[derived \| Maxwell duality compact]**; model-internal route via P5 + "arrival is an observable" available | P13, Cl(3,1), the 3+1 reading, THM-K(a) clause (a) |
| 4 | SCREEN-1 (S = A/4ℓ_P²) with LBL-1 (KMS circle) | Hawking's ¼ and ħ | [declared] | [declared] — correct; but T4f3 is a consistency of the two, not a selection of 2π | KIN-2a, THM-K–O, T4c–g |
| 5 | KIN-2a as T4's input | Schwarzschild's escape velocity, via row 4 | "the model's own profile" | [derived \| SCREEN-1, THERM-1, MASS-1, EQ-1] — as the paper says; T4's header should say so too | T4c–g, Hawking's coefficient |
| 6 | CONT-1 | an analytic continuation chosen for one rational formula | [declared] | dischargeable by T7f's real interior sector (§2.2) | THM-H interior, CURV-1, the null radii |
| 7 | NORMAL-1, STAT-1, ROT-1, PROP-1, one boost per seat, K-6 | Einstein's constraints, Killing symmetry, the PG river, flat rain rods, linearity | [declared], named loudly | [declared] — correct; nothing to change | the whole gravitational band |
| 8 | VIEW-1 | one pivot orbit presents all three alignments | [declared] | [declared] — and Euclidean (row 2) | the circle, CENSUS-C |
| 9 | the LABEL-1 physics table | the known relations | [declared] | [declared] — correct | the labelling's kill conditions |
| 10 | SEAT-hb (the ħ state is a quadratic form) | the one-sided (spinor) reading | [declared] | [declared]; label-free it is "reads the fundamental rep of the cover" (§2.4) | E-8-X, thm_g, the census |

Rows 1–3 are the ones the ledger currently mislabels. Row 2 is structural and is §4. Row 3 is the one result whose tier omits a
hypothesis it uses.

---

## 4. The two state spaces — the structural finding

**The paper and the primitives do not describe the same object.**

| | paper (v0.5, §2/§9, THM-H..O, CURV-1, PRED-1) | primitives (T7 onward) |
|---|---|---|
| the state | three unit vectors in S²; Gram G positive semidefinite | three lines in ℝ^{2,1}; Gram diagonal (−1,+1,+1) |
| symmetry of the state space | SO(3) (and the click group S₄ on labels) | SO(2,1) (and the flips B₃) |
| Δ = det G | ≥ 0; the elliptope; boundary Δ = 0 | ≤ 0 (T7c4); det G = −D²; branch locus D = 0 |
| singular points of the branch surface | four real nodes (Cayley's cubic, Aut = S₄) | **no real nodes** (T7c10) |
| the deck | V ↦ −V, V² = Δ | D ↦ −D, D² = −det G (T5c) |
| the phase / sheet meter | arg B, B = (1 + Σγ + iV)/4 on S² | V_spin = 2D in a Cl(2,1) block |
| the seat | occupies one of three spacelike axes | occupies the one timelike line |
| where time comes from | the paravector scalar (Cl(3)) | the negative line c (T7b2) |
| the hole | a boosted Euclidean Gram G′ = G + sinh²λ kkᵀ, continued through CONT-1 | a path in (t,l₁,l₂) crossing D_plane (§2.2) |

T7b1 proves one form for the frame; the frame therefore has exactly one negative line, and no seat's resolution is the
positive-definite Gram. So the Euclidean cell is not "the Euclidean seat's construction" (as the PRIMITIVES corrections say) —
there is no such seat. It is the vector part of the Cl(3) paravector: the import, surviving in every file that predates T7.
That includes the banked prediction. PRED-1's algebraic half (`suites/thm_galois_deck_descent.py`) and its operational protocol
(`suites/pred1_operational_chain.py`, `suites/pred1_physical_protocol.py`) are computed on S².

**What survives migration unchanged, what changes, what dies.** The two state spaces are the two real forms of one complex
object: the complex symmetroid {det G = 0} in the space of complex symmetric 3×3 matrices with the diagonal fixed up to signs,
acted on by SO(3,ℂ). Sign changes of the diagonal are a change of real form of the *same* complex surface; the four nodes
are still there over ℂ, at v = (i, ±1, ±1) — they are just not real on the (2,1) slice. Hence:

- **survives verbatim**: everything stated over ℚ(γ) or ℂ — the irreducibility of Δ, the rank-4 square classes, the conjugate
  Kummer module and its Galois group, the descent B ∈ F(i, √Δ), the ℤ₂ monodromy and its node dichotomy (loops about complex
  nodes do not flip). Migration is a re-parametrisation: the radicand becomes −det G = D² in seat coordinates.
- **changes**: the branch locus is a surface with no real nodes, so the click group's identification with Aut(Cayley) has no
  real content on the (2,1) slice; the geodesic-triangle solid angle becomes the hyperbolic area of T5A; the spinor lift is
  Pin(2,1), not SU(2); the Bargmann phase becomes the block-resolved spinor trace, which T5c already computed for n = 3.
- **dies or must be re-run**: VIEW-1's circle (a Euclidean circumcircle on S²), CENSUS-C, pivot_map's readings, PRED-1's three
  physical protocols (Hadamard test on qubits, spin chirality, the tritter — all realise the SU(2) lift), CURV-1's centre
  readings (dihedral angle, hemisphere), THM-H's null radii. Each has a Lorentzian analogue; none has been computed.

**This is the audit's first recommendation (§7): decide which state space the model owns and migrate the other.** The
PRIMITIVES doc has already decided (one form, (2,1)); the paper has not been told.

---

## 5. Can the physics labels be stripped? Yes. The label-free core

Everything the model computes in the primitives layer, and everything in the readout and curvature tiers, is an instance of
the following. Nothing in it names c, ħ, G, time, light, temperature, mass, or a horizon.

**Definitions.**

1. **Frame.** A real vector space V of dimension n with a nondegenerate symmetric bilinear form q of signature (p, q), and an
   ordered set of n lines ℓ₁,…,ℓₙ through the origin spanning V. The **state** is the Gram matrix G of unit representatives of
   the lines (P4); it is determined up to the sign flips {±1}ⁿ (the poles) and lives in the space of symmetric matrices with
   diagonal (±1,…,±1) fixed by q. [The seated root: n = 3, (p,q) = (2,1).]
2. **Pivot group.** O(q), acting on frames; its complexification O(q_ℂ) contains every real form O(p′,q′) with p′+q′ = n.
   A pivot by an imaginary angle is a path in O(q_ℂ) between real forms (P6). [T1–T3, LABEL-1's Wick loops.]
3. **Seat.** A choice of one line ℓ (P5). The seat's **resolution** is the q-orthogonal decomposition V = ℓ ⊕ ℓ^⊥; its
   **space** is ℓ^⊥ with the restricted form; its **rulers** are the projections of the other lines into ℓ^⊥; its
   **constant** is ℓ itself (P9, P11). The seat's stabiliser in O(q) is O(q|ℓ) × O(q|ℓ^⊥); the seat's **symmetric space**
   is O(q)/stabiliser. [T7a, T7c: the seat reads the projected angle t and not the depths.]
4. **Presentation.** The map Π_ℓ: state ↦ (Gram of the projected rulers in ℓ^⊥) — non-injective by construction (P7); its
   fibres are the depths. The **tilt** of a seat is the q-angle between ℓ and the normal to the rulers' span, cosh λ =
   √(1+η), η = vᵀS⁻¹v (T7d); it is a point of the seat's symmetric space H = O(q)/K.
5. **Readout.** A function of the state that factors through a presentation, transported by the pivot group. The natural
   class: sections of homogeneous line bundles over the boundary ∂H of the seat's symmetric space, of weight w; the pivot
   acts through the boundary Jacobian, which is a power of the Poisson kernel. [§2.9: the presentation family, the blind mass,
   the mirror, Q = dim ∂H.]
6. **Sheet.** The orientation class of the frame, sgn det[ℓ₁…ℓₙ]; read by the spinor (one-sided) representation of the pivot
   group and invisible to the vector (two-sided) one; when Cl(p,q) has a central pseudoscalar with I² = +1 (p − q ≡ 3 mod 4),
   it is the block label. [T5b, T5c, T8a.]
7. **Degeneracy strata.** The subvarieties of the state space where a principal minor of G drops rank, stratified by which
   minor and by the local vanishing order m of the corresponding determinant. Named strata for n = 3: the full minor (the
   **branch locus**, where the frame is planar); the rulers' minor (the **rulers' degeneracy**, where the seat's normal cannot
   be constructed — the seat's tilt goes to ∂H); their intersection; and, for presented sub-frames, the **presented collisions**.
   [T7c, T7e, T7f, T7g; CURV-1's r*.]
8. **Vanishing-order invariants.** On a stratum of order m: the parity of m decides whether a loop about the stratum flips
   the sheet; the value of m fixes the growth exponent of the connection curvature of the state-space submersion; m is readable
   pointwise from the locked inertia determinant. [Will's July theorem; THM-B′; GRAM_SUBMERSION; INERTIA_NODE_DETECTOR.]
9. **Extension.** A generator anticommuting with the frame's Clifford generators, unique up to two scalars, exchanging the
   sheets; its square is the one datum the readout supplies (T8a, T8b).

**The dictionary back to the labels** is then a separate, short document: c := the unique negative line under (2,1) (T7b2);
{ħ, G} := the unordered positive pair (T7b3); *time* := an ordered path in the state space read through the tilt (P10);
*horizon* := the rulers' degeneracy stratum reached by a path along which the tilt diverges (T7e); *the interior* := the
sector where the rulers' minor is negative (T7f); *the singularity* := whichever stratum a chosen path meets next, with its
order m; *temperature* := the period of the complexified boundary action in the representation the seat reads (T4); *the
fourth direction* := item 9 with square +1; *Hawking's face* := the antipodal orthant; *the sheet at the horizon* := the path
holonomy of item 6. Every physical output of the repo is one of these words applied to n = 3, (p,q) = (2,1).

**What generalises without change:** n arbitrary; (p,q) arbitrary for items 1–5 and 7–8; item 6 needs p − q ≡ 3 (mod 4) for
the block reading and otherwise needs restating; the octahedron becomes the n-cross-polytope with the Boolean lattice 2ⁿ and
sheet S_n/A_n (§2.7). **What is specific to the seated root and should be recorded as its content, not as mathematics:** the
choice n = 3, (2,1); the identification of the three lines with three constants; the labelling table; the thermal declarations;
the pinning.

---

## 6. Where original mathematics is warranted — and where it is not

**Not warranted (existing frameworks are the best fit).** The pivot group (real forms and symmetric spaces); the seat's form
and its degeneracies (linear algebra and determinantal varieties); the sheet and the extension (real Clifford algebra); the
readout tier (boundary representations of the pivot group, the Poisson kernel, principal-series weights); the octahedron (the
cross-polytope and the hyperoctahedral group); the Galois layer (Kummer theory, already the author's own); the dynamics tier
(a Hermitian structure on the seat's symmetric space — §2.10 — not a new calculus). In each case the repo either already uses
the framework without naming it, or uses a weaker hand-rolled version of it. Naming is the upgrade. Inventing here would be
the failure mode the 09-04 handoff warns against in the other direction: reaching for original mathematics where the standard
treatment is not only available but strictly stronger.

**Warranted in one place: a uniform valuation calculus for readouts at determinantal strata, across real forms.** The model
has four results that are each one corner of an unwritten theory:

| corner | statement | file |
|---|---|---|
| monodromy | odd vanishing order of Δ flips the sheet; even does not | `thm_b_monodromy` |
| curvature | ‖F‖ ∼ C(b)·Δ^{−m/2}, m the vanishing order; C(b) from the inertia tensor | GRAM_SUBMERSION |
| Kummer | square classes are independent iff the valuations along the boundary divisors are | `galois`, `thm_galois_deck_descent` |
| trichotomy | free / floored / pinned = where a positive cocycle's zero set sits relative to a group orbit, a conjugation class, a base map | `sect1`, 08-28 handoff item 2 |

What is missing is the theory that contains all four: for a Gram symmetroid over a field, stratified by rank and by vanishing
order of its principal minors, and for a readout of weight w on the boundary of the pivot group's symmetric space, *what does
the readout do at each stratum, as a function of (w, m, which minor, which real form)?* The seated root needs this in exactly
one place — the black hole: what a seat reads at the horizon (rulers' degeneracy, m = 1), at the branch locus (full
degeneracy), at a presented collision (CURV-1's logarithmic singularity), and whether the answer depends on the real form
(the interior sector). That is the model's founding question ("the singularity as a locus with computable structure") stated
as a theorem target rather than a wish, and no existing body of mathematics answers it as posed, because it mixes
representation-theoretic weights (§2.9) with algebraic strata (§2.11) at their meeting point. It is the best case for
original work because (i) all four corners are already proved in the repo, (ii) two of them are the author's own prior theorems,
(iii) the general statement is label-free and dimension-free, and (iv) its physics payload is precisely singularities and
horizons. **Kill condition for the proposal:** if the readout's behaviour at a stratum is determined by m alone (no dependence
on w or on the real form), the "calculus" is the existing vanishing-order theorem and there is nothing new to write.

---

## 7. Recommendations, in order

1. **Rule on the state space** (§4). The PRIMITIVES doc has already ruled (one form, (2,1)); write the ruling into the paper
   as v0.6's scaffold entry: the Euclidean cell is BARE-1's residue, not the bare tier. Then migrate in this order, each a
   runner: (a) the Kummer/descent layer in seat coordinates (radicand −det G = D²; nothing should change); (b) the Bargmann
   phase → the block-resolved spinor trace and the Lorentzian triangle area (T5A + T5c already contain it); (c) PRED-1's
   operational reading restated in Pin(2,1) — decide whether the three qubit protocols have Lorentzian analogues or whether
   PRED-1 is a 3+1 statement that needs Γ₄ first; (d) THM-H/CURV-1's hole as a path in (t, l₁, l₂) crossing D_plane, with
   CONT-1 replaced by T7f's real interior (§2.2).
2. **Relabel three ledger entries** (§3 rows 1–3): signature (2,1) as declared by P3/P11 with T1–T2/T7a as consequences;
   T8c as derived | Maxwell duality compact, with the P5 route named as the removal of that hypothesis; T4f3 as a consistency of two declared conventions.
   None of these change a number.
3. **Run the P5 route for Γ₄²** (§2.6): s = +1 from P5, "arrival is an observable" as a declared model principle, signature
   from T8b's identity. If it lands, T8c's tier drops the Maxwell hypothesis.
4. **Name the readout tier** (§2.9): one short note, no new computation, stating the c- and ħ-seat families as the reciprocal
   Poisson kernels of H³ and H², Q = dim ∂H, the blind mass as harmonic measure, the mirror as the weight duality. Then bank
   the consistency that Q_c = 2 presupposes Γ₄, and re-derive DEBT-2b's flat band from the boundary measure.
5. **Write the label-free core** (§5) as its own document, with the dictionary as a separate section, so that every future
   runner's header can cite a definition instead of a physical word. This is cheap and it is what makes the framework usable
   for any n and (p,q).
6. **Do not invent a dynamics calculus.** Take the Kähler dial to the (2,1) frame space as a Hermitian structure on the seat's
   H², identify it with the impedance (§2.10), and ask which structure is forced. One runner decides whether the impedance is
   a free function of state or a scalar.
7. **Open the valuation-calculus target** (§6) as a THM-TARGET with the four corners as its regression tests and the horizon
   readout as its first theorem. This is the one place the audit endorses original mathematics.
8. **Housekeeping.** Record in §2.5 that the block-label sheet needs p − q ≡ 3 (mod 4); record the B₃/Spin(7) count in
   LABEL-2 as a coincidence with a one-runner kill; give T4's header the pinning's actual provenance.

---

## 8. Kill conditions for this audit's own claims

- §2.1/§2.2: if a runner can produce signature (2,1) from primitives that do *not* assign a character to any plane, the
  "declared signature" finding is wrong and T1/T7a are derivations. (P3 as written assigns characters; the finding stands until
  P3 is weakened.)
- §2.6: if 2-form duality is shown to act on the seat's compact ruler plane as a real rotation *from the model's own
  objects*, T8c's premise becomes model-internal and the Maxwell hypothesis comes off its tier.
- §2.9: if s^{−Q} ≠ the round-measure boundary Jacobian of the pivot on ∂H^{Q+1} for either seat, the Poisson-kernel
  identification dies. (RN-1 and RN-2 are that check; they pass.)
- §4: if a seat's resolution can produce a positive-definite Gram of the three lines under P3/P11, the Euclidean cell is a
  seat's construction and the two-state-space finding softens to "two seats". T7b1 says no seat does.
- §6: if the readout's behaviour at every stratum is a function of m alone, there is no calculus to write.

# THM-TARGET J — does the model's own geometry produce the field law?  (schedule item 0b)
Written 2026-09-03.  `thm_j_dyn.py` **7/7**, `run_thm_j_dyn.log`.  Read after `docs/2026-09-03-thm-target-I.md`.

## Verdict: NOT DERIVED — and the reason is structural, not numerical
1. **Theorem (D-1).**  Any gradient-quadratic static energy Σ M_ab(g)∇g_a·∇g_b of click-invariant fields, restricted to the
   spherically symmetric pinned family, is F(λ)λ′²; its Euler–Lagrange equation is (r² ds/dr)′ = 0 with s = ∫√F dλ.
   A gradient energy makes its own ARC LENGTH along the family harmonic.  So "produces KIN-2a" means "the energy's arc
   length is affine in sech²λ"; to fourth order, s₄/s₂ = −2/3.
2. **Obstruction (D-2).**  An affine function of sech²λ has a pole where the presented volume vanishes (the centre) and is
   bounded at the horizon.  Every polynomial invariant of the presented cell (ρ_K, sinh²λ, ρ_i, √ρ_i, det G′) does the
   opposite.  So no polynomial gradient energy can produce KIN-2a: the energy must be built from the RECIPROCAL of the
   volume invariant, Δ/det G′.
3. **The model's own candidate (D-3).**  The only frozen, choice-free geometry the model owns is the state-space metric of
   `docs/GRAM_SUBMERSION_CURVATURE` (round metric on (S²)³, mechanical connection, horizontal lift, quotient metric on
   the elliptope).  Its harmonic-map energy along the pinned family:

   | frame / layer | s₄/s₂ (need −0.6667) | c₂ | perihelion coefficient (Mercury: 3) |
   |---|---|---|---|
   | N1 / A (seat) | +0.0457 | −0.712 | 1.575 |
   | N1 / C (hole) | −0.0992 | −0.567 | 1.865 |
   | N3 / A | +0.0829 | −0.750 | 1.501 |
   | N3 / C | −0.0876 | −0.579 | 1.842 |

   Fit converged to 1e-8; including the vertical part of the presented lift changes nothing; the arc length is not
   proportional to tanh²λ.  The profile the frozen metric would produce is definite, frame- and layer-dependent, and
   excluded at Mercury by a factor of two.

## What this settles
- The bare tier's geometry is a metric on **directions**.  The field law lives in the **reciprocal volume**.  A metric on
  directions cannot see the volume's reciprocal, and the numbers say it does not.
- The energy that produces KIN-2a is therefore exactly the gradient energy of Δ/det G′ — one variable, one exponent — and
  it is not in the bare tier.  It has to be declared: **KIN-2a″**: E_static = ∫|∇(Δ/det G′)|² d³x, sourced by the
  blind-mass measure.  That is a Newton's law in the model's own variable.  It is not an Einstein: nothing in the model
  yet says why the reciprocal presented volume, rather than any other function of the cell, is the thing whose gradient
  costs energy.

## What would change the verdict
A dynamics-tier energy (nine generators, stretch–angle couplings, the dissipator's channel law) whose static reduction is
a function with a pole at det G′ = 0.  The obstruction D-2 says where to look: not among polynomials.  The one frozen
object with the right shape is the capacity reading of Δ (scale-blind, a volume) — inverted.  A test target: does the
seat dissipator's pivot-invariant bath (c(ω)²g₀(ω) = Cω², DEBT-2b) have a static limit whose energy is the reciprocal
presented volume?  That is the next place the model could earn the law rather than declare it.

## Comparison stage
The Dirichlet energy of N² = 1 − 2Φ/c² is the Newtonian field energy; the elliptope quotient metric is the shape-space
metric of three unit vectors (mechanical connection, Littlejohn–Reinsch); harmonic maps from ℝ³ into a round quotient
target are sigma-model configurations, and sigma models do not give 1/r potentials — which is what the table shows.

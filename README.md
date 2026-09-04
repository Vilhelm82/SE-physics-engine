# SE-physics-engine — the seated root

A derivational-economy study. One root, three lines through it, three planes between
them, and an observer *seated* on one line. Everything else — the metric signature, the
pivot group, the light cone, the horizon, the double cover — is derived from that, or it
is written down as a declared input with its removal route.

The point is not to predict something general relativity does not. It is to find out how
little has to be assumed before the usual structures become forced, and to keep an honest
ledger of what was assumed anyway.

## The rule that shapes the repo

**No representation that already contains the target.**

The paper's first draft used Cl(3) as its algebra. Cl(3) is the Pauli algebra: its
paravectors carry the Minkowski norm and its rotors *are* the Lorentz group. So every
"the c seat reproduces relativity" result had been written in relativity's own algebra
and could not have failed. That is not a derivation; it is a construction steered by a
known answer.

The test applied since: *could the physics have failed to come out?* If not, the
mathematics has been scoped to the answer. Cl(3) is banned as an input. Everything the
model now uses — the (2,1) form, the double cover, the four-generator extension — is
reached from the primitives or declared with its full cargo.

## The primitives

`docs/2026-09-04-PRIMITIVES-v0.md` is canonical and is still a draft.

| | |
|---|---|
| **P1** | A point: the root. |
| **P2** | Three lines through it. Each line has two ends (poles). |
| **P3** | Three planes, each spanned by two lines, each with fixed intrinsic character (hyperbolic or compact) that never deforms. |
| **P4** | The state is the three angles between the lines. Nothing else. |
| **P5** | A seat is one line taken as reference. The seat cannot see its own line; both poles collapse to the root. |
| **P6** | A pivot is a rotation of the frame by an *imaginary* angle. Not a spatial motion. |
| **P7** | A presentation is what a seat reads after a pivot. Presentations are non-injective: distinct states can present identically. |
| **P8** | Axes cannot be moved in isolation. |
| **P9** | A constant is what a view collapses to a point. Constancy is a property of the view. |
| **P10** | *c* is the negative line of the frame's form — derived, not declared. "Timelike" is the form's word; "time" needs a path through ordered states. |
| **P11** | A seat resolves the frame into one point (its constant), two imaginary lines (its rulers), one real plane (its space). |
| **P12** | Laws are seat-derived. Another seat derives another physics from the same frame. |
| **P13** | The fourth direction is the generator that exchanges the two Cl(2,1) sheets. Derived. |

## Results

Each is a Python file with a matching `.log`. Re-running the file *is* the verification;
every suite finishes in seconds.

**The pivot group.** Three lines and "rotation by an imaginary angle" give SO(2,1). The
count of boost planes over all realifications is `{0, 2}` — never three. On its own
primitives the model is **2+1**; Cl(3) had been supplying a fourth direction silently.
→ `prim_t1_t3_pivot_group.py`

**Hawking's coefficient.** The Wick face has period 2π on the vector reading and 4π on
the cover, and the vector reading is Sym² of the spinor one — so the cover is *reached*,
not imported. The twist is logarithmic with scale 4·r_s; κ from the twist and κ from the
seat's redshifted acceleration agree; T_H follows by three independent routes.
→ `prim_t4_hawking_period.py`

**The seat's form.** The seat's three planes glue to a unique Gram matrix with diagonal
(−1, +1, +1). In seat coordinates det G = −cosh²l₁ cosh²l₂ sin²t; the branch locus is a
surface, and on it frames of different depth present identically (P7, concrete). The seat
reads one of its three parameters directly.
→ `prim_t7_seat_form.py`

**Time dilation as tilt.** The seat decomposes *c* against its rulers' span; η = vᵀS⁻¹v
and sech λ is the lapse. The paper's Schwarzschild pinning is a value of the seat's own
tilt invariant. The lapse is bounded; the horizon is η → ∞, and it is *not* the branch
locus.
→ `prim_t7d_tilt.py`, `prim_t7e_two_degeneracies.py`, `prim_t7f_regions.py`

**Wigner is area.** The rotation from composing two non-collinear boosts equals the
hyperbolic area of the geodesic triangle they enclose on H². Two routes: angle defect at
40 digits, and an exact identity in rational coordinates.
→ `prim_t5_area_and_lifts.py`

**The sheet.** Cl(2,1) = M₂(ℝ) ⊕ M₂(ℝ); the two blocks are two sheets, superselected. The
observable that reads them is the raw spinor trace 𝒱 = tr(Γ_c Γ_ħ Γ_G) = 2D, with
D² = −det G — the *signed* volume of the frame. Reversal sends D → −D. The block label is
the sheet, and the full trace is sheet-blind.
→ `prim_t5b_parity.py`, `prim_t5c_corrections.py`

**The fourth direction.** A generator anticommuting with the seat's three exists, is
unique up to two scalars, and exchanges the sheets — lifting the superselection so they
can interfere. Its signature is fixed by electromagnetic reciprocity: ⋆ acts on 2-forms
only in four dimensions, ⋆² = (−1)^n₋, and Maxwell's *compact* duality rotation needs
⋆² = −1. Hence Γ₄² = +1: spacelike, Cl(3,1).
→ `prim_t8a_fourth_direction.py`, `prim_t8b_adjoint.py`, `prim_t8c_reciprocity.py`

**The labelling.** The six poles form an octahedron whose eight faces are the Boolean
lattice of three flips. Our face is the empty set; Hawking's is the full set, its exact
antipode. The sheet above a face is the parity of the path that arrived there.
→ `label2_octahedron.py`, `label3_lattice.py`, `label3b_corrections.py`

## Tiers

Every claim carries one, and they are never blurred:

- **proved** — an exact symbolic identity, verified in a runner
- **derived | X** — follows from the primitives given declared input X
- **retrodiction** — reproduces a known result
- **conjecture** — stated, with a kill condition
- **declared** — assumed, with its cargo and its removal route named

A result without a stated kill condition is recorded as having none. That is information
too.

Failures are kept. `THM-J` is the first model-internal candidate the model itself killed;
it is in the repo because a framework that cannot lose is not saying anything. Corrections
are recorded inside the runners that made them, next to the check that caught them.

## Running a suite

```bash
python3 prim_t7_seat_form.py        # prints PASS/FAIL per check and a verdict
python3 prim_t8c_reciprocity.py
```

Requires `sympy` and `mpmath`. Each file is self-contained and states its inputs, its
checks, and what it does *not* establish in a header block. The `.log` beside it is the
committed output.

`build_html.py` renders `PAPER-seated-root-v0.5.md` to self-contained HTML.

## State

The paper is at v0.5 and is known to carry an undeclared import (Cl(3), above); v0.6 owes
that scaffold entry. The primitives document is a draft. The mathematics below it is
green.

**Open:**

- **The scalar and the coupling.** A response map must be built from the seat and a
  charged load *before* either positive ruler is named, and must yield the vacuum
  impedance — unit-independently, ζ = e²Z₀/ħ = 4πα. Done so far: the cone. Open: the
  scalar, the coupling, and the theorem that the seat's two ruler responses *are* the two
  electromagnetic constitutive channels. Prove it and ħ and G are separated by a
  measurement rather than a label.
- **The dynamics tier.** The frame is static; the seat constructs a rate but nothing
  propagates. The candidate shape is V = IZ — the tilt as potential, the seat's load as
  impedance, the state-rate as current — with the load's impedance as a function of state
  the single unknown. The tangential sector, the horizon sheet and the propagator are all
  that one unknown.
- **The horizon sheet.** Whether trapped-surface formation selects one orientation of the
  frame path. The mathematics has supplied the variable; the collapse dynamics has to say
  which sign is the horizon.

## Provenance

Written by Will Lloyd. The newest handoff file is the state of the project; every older
one is provenance. Commit messages are the detailed record — they are written as receipts,
and they are more complete than this file.

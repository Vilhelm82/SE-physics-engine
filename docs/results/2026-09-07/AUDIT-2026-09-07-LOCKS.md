# LOCK REGISTER — 2026-09-07

**Directive (Will, 09-07):** any part of the model structurally dependent on theory rather than
measured observation is a LOCK. Locks are to be removed and the quantity re-derived on the model's
own steam, or demoted to DECLARED / CONJECTURE until it is.

**Ground rule used here:** GROUND = a number somebody measured, with its precision.
THEORY = a formula that fits ground and extends past it. Extension past the measured order is
scaffold regardless of how good the mathematics is.

Audit basis: `HANDOFF-2026-09-04-evening.md`, `docs/LABELLED-MODEL.md`, `docs/CONJECTURE-COSMOLOGY.md`,
runner headers (`prim_t4_hawking_period.py`, `prim_t7f_regions.py`, `thm_n_kerr_quadrupole.py`).
Remote HEAD 2026-09-05 12:16. Nothing local was visible.

---

## Ground available to the model (measured, with precision)

| Quantity | Measured how | Precision | Order of the metric it pins |
|---|---|---|---|
| γ_PPN (space part, 1st order) | Cassini Shapiro delay | γ−1 = (2.1±2.3)×10⁻⁵ | g_ij at O(U) |
| Redshift (lapse, 1st order) | GP-A, Galileo eccentric, optical clocks | few ×10⁻⁵ | g_00 at O(U) |
| β_PPN (lapse, 2nd order) | Mercury perihelion + MESSENGER | ~10⁻⁴ | g_00 at O(U²) |
| Frame dragging | GP-B, LAGEOS | 19 % / few % | g_0i at O(U), Lense–Thirring form |
| Orbital decay (radiation) | double pulsar | ~10⁻⁴ | PN dynamics, low η |
| S2 precession + redshift | GRAVITY | ~15 % | g_00 at η≈7×10⁻⁴ |
| Light-ring scale | EHT M87*, Sgr A* | 10–17 %, mass-degenerate | η ≈ 2 region, model-mediated |
| Merger/ringdown | LIGO/Virgo | tens of % | near light ring, **template-mediated (Kerr NR)** |
| α | atomic/QED | 10⁻¹⁰ | — |
| c, ħ | SI exact by definition | exact | — (conventions, not measurements) |
| G | CODATA | 2×10⁻⁵ | — |
| Blackbody spectrum | Planck's law, lab | very high | station only (Planck's face) |
| Thomas precession | atomic fine structure | high | Wigner rotation, low velocity |
| CMB flatness / isotropy | Planck satellite | Ω_k = 0.000±0.005 | cosmology targets only |

**No ground exists for:** the lapse beyond O(U²); the exact light-ring radius; the horizon (inferred by
absence of surface emission); anything at η > ~3; the interior; r = 0; the white-hole / past-horizon
region; Hawking or Unruh temperature; Tolman–Ehrenfest equilibrium; Parikh–Wilczek tunnelling; greybody
factors; ISCO (Kerr-assumed X-ray fits); marginally bound orbit.

---

## The locks

### LOCK-1 — The pinning η = r_s/(r − r_s)  [T7f]
- **Imports:** Schwarzschild g_00 = 1 − r_s/r, the whole series, via the presented Gram coefficient.
- **Ground it actually has:** lapse to O(U²) (β), space part to O(U) (γ). Everything past that is scaffold.
- **Depends on it:** N² = 1/(1+η); horizon ↔ η→∞; interior ↔ η<−1; centre ↔ η=−1; the r-labels on T7g's
  region table; photon sphere at η=2; tilt = free-fall rapidity (the "river"); ISCO at η=½; the marginally
  bound point η=1; redshift √3; T4's twist scale 4r_s (LOCK-2 inherits); EQ-1's r-dependence; the whole
  r-dictionary in CONJECTURE-COSMOLOGY.md.
- **What survives without it:** the sector geometry itself (f0–f4, forbidden zone, D_plane, B_form, the
  seat switch, sin t = 0 as the ground with two components, three-seat factorisation η_seat = q(seat,seat))
  — these are facts about the Gram and carry no r.
- **Removal route:** derive η(r) from the dynamics tier (frame path; V = IZ; the flow tanh λ = f(r) as
  output). Check against β, γ at measured order. The O(U³) lapse term and the exact light-ring radius are
  then *predictions*, testable by next-generation EHT/GW.
- **Tier after removal:** the pinning becomes **SCAFFOLD (measured to O(U²) lapse, O(U) space; unmeasured
  beyond)**. All downstream r-statements inherit that tier until the dynamics tier reproduces them.
- **Kill:** if the model's own flow is not 1 − r_s/r to O(U²), the model is wrong at measured order and dies
  there — the honest way.

### LOCK-2 — Hawking coefficient T_H = ħc/(4π k_B r_s)  [T4]
- **Imports:** (a) the identification "Wick period ↔ 1/k_BT" (KMS condition; equilibrium stat mech —
  theory); (b) the twist scale 4r_s (surface gravity κ = 1/(2r_s), i.e. LOCK-1's tail at the horizon);
  (c) the Euclidean-regularity route (theory); (d) the first-law route uses dM = T dS with the area law
  (Bekenstein — theory).
- **Ground:** none. No gravitational Hawking quantum has been observed; none will be at astrophysical mass.
- **What survives:** the *structure* — the Wick face exists, W is an eighth-turn, vector period 2π /
  2-dim 4π, the cover REACHED as Sym² (T4b'), three routes agreeing with each other. That is a derived
  fact about the octahedron and is kept.
- **Removal route:** demote the temperature identification to DECLARED (a hypothesis in the list, per
  "nothing is banned"). Drop the 4r_s scale until LOCK-1 is removed. T_H's *number* becomes CONJECTURE.
  If the dynamics tier produces its own scale at the open circuit, T_H is re-derived; otherwise it stays
  a named conjecture with three internally consistent routes and no ground.
- **Tier after removal:** structure DERIVED; identification DECLARED; number CONJECTURE.

### LOCK-3 — Unruh temperature (ħ seat under acceleration)
- **Imports:** same KMS identification as LOCK-2. **Ground:** none (never observed).
- **Removal route:** same as LOCK-2. Keep only the shared-station consistency (Unruh route = Hawking
  route) as a DERIVED relation between two DECLARED identifications.

### LOCK-4 — Frame dragging "reproduced"  [thm_l / thm_m / thm_n, Kerr quadrupole]
- **Imports:** Kerr / Lense–Thirring as the comparison, and (per runner names) as input to the swirl.
- **Ground:** GP-B 19 %, LAGEOS few %, leading order only.
- **Status:** a BAND on Lense–Thirring at measured order; anything beyond leading order is scaffold.
  These are Cl(3)-era THM files — check whether the swirl was *reached* or *read in* before crediting.
- **Removal route:** produce g_0i from the frame path; compare to the measured leading coefficient.

### LOCK-5 — Octahedron far column (Temperature / Evanescence / Energy)  [LABEL-2/3]
- **Imports by station:** Tolman–Ehrenfest (theory, unmeasured); Parikh–Wilczek (theory); greybody (theory);
  finite-T instantons (theory, though fusion rates are ground); Hawking CHG (theory).
- **Ground:** exactly one station — Planck's face (CG), blackbody, measured.
- **What survives:** the lattice as *geometry* (8 faces = 2^{C,H,G}, null/boost classes, antipode) is
  DERIVED and carries no physics import. The *names* are provisional labels, which is fine.
- **Removal route:** relabel: each station carries GROUND / THEORY. "Hawking radiation is the unique
  conjunction" is a labelling fact, not physics, until LOCK-2 is removed.

### LOCK-6 — EQ-1 Tolman equilibrium
- **Open ruling from earlier sessions: "Will's principle or Mercury's fact?"** Resolved by this directive:
  **neither — it is a target.** T√g_00 = const is theory (GR + thermodynamics), unmeasured.
- **Removal route:** derive equilibrium from the divider (V = IZ). The condition "no work through the load
  ⇒ no time" at the open circuit is the model's candidate. Until derived: CONJECTURE.

### LOCK-7 — The white-hole sheet  [CONJECTURE-COSMOLOGY.md]
- **Imports:** the maximal analytic extension (Kruskal) as the thing the sheet exchange Γ_C → −Γ_C
  corresponds to; all r-language via LOCK-1.
- **Ground:** none for the past horizon; GR itself does not claim it is realised.
- **What survives:** the sgn D computation on the I = −1 block is a Gram fact and stands on its own.
  The *targets* (flatness Ω_k, isotropy, arrow) have ground; the *route* is scaffold on scaffold.
- **Tier:** stays CONJECTURE / DIRECTIVE. No change to the computation order.

### LOCK-8 (soft) — EM reciprocity ⇒ Cl(3,1)  [T8b]
- **Imports:** ⋆² = −1 on 2-forms (Hodge structure of Maxwell). Lorentz reciprocity is measured (every
  antenna lab); the *Hodge formulation* is theory.
- **Ruling:** keep as GROUND-ADJACENT; state explicitly that the input is "reciprocity holds", not
  "Maxwell's equations". Low priority.

### Not locks (checked)
- T1–T3 (Cl(2,1) reached from primitives) — clean.
- T5/T5b/T5c (Wigner = hyperbolic area; V_spin = 2D) — clean; Thomas precession is ground at low v.
- T7a–e (form constructed by the seat; branch; sectors as geometry) — clean.
- PPN γ = 1 forced — the model predicting ground. Exemplar of what the directive wants.
- HUNCH-Z0 ζ = 4πα — exact algebra on measured α.
- Bargmann / Galois descent — pure mathematics.
- c, ħ, G as axes — DECLARED (a choice), values are SI/CODATA ground.
- Today's structural surplus (three-seat factorisation, null edges, lit face, P₂ fold, magic plane = CHG
  face-centre tangent) — Gram facts, no import. The r-values attached to them are LOCK-1.

---

## Order of removal

1. **LOCK-1** (everything else waits on it). Deliverable: η(r) as output of the dynamics tier.
   Pass condition: 1 − r_s/r to O(U²) in the lapse, γ = 1 in the space part. Prediction: O(U³).
2. LOCK-6 alongside (the divider *is* the dynamics candidate; EQ-1 is its first theorem or its first casualty).
3. LOCK-2/3: demote now (one-line tier edits), re-derive after 1.
4. LOCK-5: relabel now.
5. LOCK-4: read the three THM runners for reached-vs-read-in; band them at measured order.
6. LOCK-7: unchanged; sgn D runs on schedule.

## Edits owed to the repo (when Will is at the workstation)
- LABELLED-MODEL.md: add a GROUND/THEORY column to the station table; retier T_H and T_U.
- HANDOFF: replace "pinning" with "scaffold pinning (measured to O(U²))" wherever it is load-bearing.
- CONJECTURE-COSMOLOGY.md: header note that all r-language is LOCK-1.
- PRIMITIVES-v0: still DRAFT. Freeze first; none of the above is the model's until P1–P13 are.

---

# LOCK REGISTER / STATION STATUS — UPDATE 2026-09-08

Appended to `AUDIT-2026-09-07-LOCKS.md` and `LOCK5-STATIONS-2026-09-07.md`. See
`2026-09-08-EDGE-chain.md` for the full results and the four retractions (R-1 … R-4).

## Register changes

- **LOCK-4 (frame dragging).** Re-tiered: DERIVED given **MAT-2** (inductive swirl) + Ampère/induction (EXT-036),
  strength K still pinned to LAGEOS at 1 %. The earlier "DERIVED given G1/NS/FAX" wording is superseded: Faxén's ½ is
  retained as GEM's ½, and the Stokes constitutive reading is retired (R-2). NS-1's finite-Re correction carries over
  under the induction operator with diffusivity for kinematic viscosity.
- **LOCK-5 (station labels).** The scaffold edge is **closed**. Face H was the last: its geometric-optics fold (09-07)
  now has the wave dressing at both ends — the electrode's power law |T₀|² = 4(ωr_s)² (EDGE-1) and the fold's Fermi edge
  of width λ_L (FOLD-1). Face HG (Parikh–Wilczek) remains THEORY: the model has not derived a tunnelling *emission*
  channel, only the transmission edges. Retier HG explicitly rather than letting H's closure cover it.
- **New pin, named:** the membrane surface viscosity η_s (EDGE-4) is a **free transport coefficient**. The tidal-friction
  prediction is conditional on Damour's 1/(16πG) until the medium derives it. Same status as K's pinning to LAGEOS:
  honest because it is pinned to something outside, dishonest if quoted as derived.

## Station status changes (`docs/LABELLED-MODEL.md` §4)

| flips | was (09-07) | now (09-08) | receipt |
|---|---|---|---|
| H | THEORY, with a DERIVED geometric-optics limit | **DERIVED** (both ends) | FOLD-1 (fold, Fermi edge λ_L) + EDGE-1 (electrode, 4(ωr_s)²); DYN-2 supplies the thermal factor between them |
| HG | THEORY | **THEORY** (unchanged, and now stated separately) | no emission channel derived; H's closure does not cover it |
| CHG | DERIVED-given-KMS · PREDICTION | unchanged, **plus** the greybody factor is now derived at both ends | the spectrum is Planck(κ) × [Fermi-edge(λ_L) at high l, 4(ωr_s)² at l = 0] |

## Residue and pins, consolidated

One declared line (KMS, EXT-029) and three pinned-to-ground quantities: d = 3 (T8c, reciprocity ground-adjacent),
K = 2GJ/c² (LAGEOS, 1 %), η_s (Damour, EDGE-4). Everything else in the DYN/EDGE/GUD/NS chain is derived from Ohm,
Laplace, Euler, induction, and the primitives.

## Predictions register (new section; four live, all with instruments)

1. κ = 0 → δκ_s = −1 (GWTC-4.0 cannot distinguish from Kerr's 0; needs σ ≲ 0.3, factor 30–100).
2. Tidal friction ~1/36 of Kerr's, conditional on η_s = Damour (LISA EMRI horizon flux; LVK horizon absorption).
3. Subring ladder e^{−π} per half-orbit; λ_L = Ω_c; greybody edge width λ_L ≠ κ (next-gen space VLBI).
4. No swirl correction at any order (Kerr has O(J³)); isotropic g₀₀ at O(U³) = −3/2 (unmeasured).

# HANDOFF — 2026-09-07 (evening)  — LIVE STATE

Continues `handoffs/HANDOFF-2026-09-05.md` (the single current session handoff): this file supersedes ONLY its
"original black-hole programme" continuation item (native material selection); the quantum targets there stand unchanged.
Written while Will was away from the workstation, against the GitHub remote (HEAD 09-05 12:16, flat tree, 61 commits behind
local). **I worked all day from HANDOFF-2026-09-04-evening, two versions stale.** Reconciliation with the 09-05/06 local state
is in the section below and was done BEFORE anything was placed in the tree.

## Reconciliation with 09-05/06 (read first)

- **Material law.** 09-05 left "a solved passive crossing example exists; a unique material law selected by the primitives does not."
  Today's DYN-1/DYN-3/NS-0 supply a DECLARED material law — an Ohmic/Euler continuum (Laplace outside the source, Euler river,
  Stokes swirl, no-slip electrode) — of the same tier as Codex's declared passive load, but tested against ground: beta = gamma = 1,
  Mercury 42.9807", deflection 1.7512", GP-B 39.1 from LAGEOS. It is NOT selected by the primitives. File it as MAT-1 (declared).
- **The pinch.** What I called "B_form, the corner of the t = 0 sheet where l1 = l2" this afternoon IS the PINCH of PINCH-1/2/3
  (P_sigma = {gamma = sigma, a = sigma b}, codim 2, 81 checks, Kummer coordinate). Cite those; do not re-narrate.
- **SHEET-1 qualifier.** The river path never touches the pinch, so SHEET-1's "no fluid path changes sheet" is "no GENERIC path."
  Codex's load transport (09-05) crosses the pinch and changes sheet: a codim-2 route. Both hold. The white sheet is a reading
  along generic paths and a destination only through the pinch.
- **HG-1 is not GUD-1.** HG-1's swap is hbar <-> G (D -> -D). GUD-1/2's swap is load <-> source (eta -> 1/eta, D unchanged).
  Distinct objects; both exist.
- **P10/P11 correction (to Will, tonight).** Both 09-05 foundation audits: signature (2,1) is DECLARED by P3/P11 with T1/T2/T7a as
  consequences. I said T7a/T7b derive P10/P11; that overclaimed. T7b2a shows one negative line GIVEN the constructed form.
- **Ledger.** AGENTS.md (09-05) requires every adopted external result to be recorded with source, use, assumptions, and whether it
  is retained standard mathematics. Today's imports are EXT-027..EXT-033 (appended to docs/EXTERNAL-MATHEMATICS-DEBT.md). The lock
  register (AUDIT-2026-09-07-LOCKS.md) is the same discipline from the physics side; it does not replace the ledger.
- **"Do not mark foundational debt closed merely because a check passes."** Every "LOCK-n -> DERIVED-given-X" below is a
  conditional derivation with X named. None closes a foundational debt.
- **Placement.** Runners -> suites/ (with .log). Docs -> docs/results/2026-09-07/. This file -> handoffs/.

## What happened (one paragraph)

Will directed: any part of the model structurally dependent on theory rather than measurement is a LOCK; remove it and let the
model derive the quantity itself. A lock register was drawn (8 locks). The removal tool turned out to be Will's own
lapse-as-voltage-divider: with the horizon as an electrode and Ohm's law in a 3-D medium, the exact Schwarzschild lapse is the
divider of two spreading resistances. That single object then derived Tolman equilibrium, the Hawking scale (leaving one
declared line), frame dragging as the medium's Stokes swirl, the sheet holonomy of a collapsing fluid, and the dressing of the
light-ring fold. Seven runners, 93/93 checks, sandbox only.

## Runners produced today (all in /mnt/user-data/outputs; NOT on the workstation yet)

| runner | checks | result | tier |
|---|---|---|---|
| prim_dyn1_ohmic_divider.py | 16/16 | A(r) = 1 − r_s/r as OUTPUT of Ohm + T8c(d=3) + THM-I(B=1/A); β = γ = 1 read off; Mercury 42.9807″ (two paths, 2e-7); deflection 1.7512″; d=2 control fails; U³ coeff −3/2 | DERIVED given Ohm + two T7f declarations |
| prim_eq1_tolman_from_divider.py | 12/12 | T·N = T_∞ by mode balance and Nyquist balance; no static equilibrium inside (theorem); horizon T→∞; dual sheet swaps the hot end (flagged) | DERIVED given DYN-1 + SI + definition of equilibrium |
| prim_dyn2_electrode_period.py | 12/12 | κ = 1/(2r_s), Rindler form N = κρ − ρ³/(16r_s³), polar Euclidean loop, β = 4πr_s/c derived; T_H recovered GIVEN KMS; LOCK-3 = LOCK-2; dual sheet cold (flagged) | DERIVED-given-KMS |
| prim_dyn3_stokes_swirl.py | 10/10 | swirl Ωa³sinθ/r² from Laplace + no-slip; dipole gyro pattern via Faxén; node/gyro = 4 structural; K pinned to LAGEOS predicts GP-B 39.1 vs 37.2±7.2; ωr³ = const at all orders (Kerr has O(J³)); **wording corrected late: swirl parity is odd in J; the O(J²) mass quadrupole (THM-N) is separate and open** | DERIVED given G1+NS+FAX+RIDE; strength pinned to ground |
| prim_sheet1_block_and_holonomy.py | 17/17 | sheet swap Γ_C→−Γ_C unique (reflection alternative breaks anticommutation); I=−1 reverses the interior clock → **kill 2 does not fire**; fluid path from the river has D = (r₀/r)^{3/2}, never 0, → ∞ at the centre: sheet inherited, white sheet is a READING not a destination; θ₊=0 at river=c, θ₋=0 on I=−1 | DERIVED; fluid→state map given RIDE-2 + P11 ratio convention |
| prim_fold1_lightring_dressing.py | 10/10 | λ_L = Ω_c (two paths, 2.6e-6); subring ladder e^{−π} per half-orbit (PREDICTION); eikonal QNM = the H flip; greybody edge = parabolic-barrier Fermi function, width λ_L ≠ κ, λ_L/κ = 4/(3√3); M87* shadow 39.7 μas vs 42±3 | DERIVED from the divider output + one exact formula; ringdown ballpark only |
| (docs) AUDIT-2026-09-07-LOCKS.md, LOCK5-STATIONS-2026-09-07.md | — | the register; the station relabel | — |

Sandbox re-runs today of existing runners: label3_lattice 14/14, label3b_corrections 10/10, label2_octahedron 11/11 (geometry intact).

## Lock register — status tonight

- LOCK-1 pinning → **DERIVED-given-Ohm** (DYN-1). All r-language downstream inherits this.
- LOCK-2 T_H, LOCK-3 Unruh → **merged; DERIVED-given-KMS** (DYN-2). KMS = the ONE declared line on the far column.
- LOCK-4 frame dragging → **DERIVED given G1/NS/FAX/RIDE; K pinned to LAGEOS** (DYN-3). KIN-3 retired; THM-M (on BARE-1) archivable.
- LOCK-5 station labels → **relabelled**; two scaffold faces (H, HG) share one edge; H derived at the fold (FOLD-1); l=0 electrode tunnelling still owed.
- LOCK-6 EQ-1 → **DERIVED**; the open ruling ("Will's or Mercury's") is closed: neither, it is a theorem of the divider.
- LOCK-7 white sheet → **CONJECTURE (unchanged tier); kill 2 survived; kill 1 (the number) is the only kill with teeth; kill 3 (signalling) runner owed.**
- LOCK-8 EM reciprocity → soft; state the input as "reciprocity holds" not "Maxwell". Clerical.

Pins remaining (ground-pinned, not derived): d = 3 (T8c, reciprocity ground-adjacent); K = 2GJ/c² (LAGEOS, 1 %); KMS (declared).

## Declarations spent today (name them or they become imports)
- T7f-load: the load is the shell from the open circuit to the seat. (Will, pre-existing)
- T7f-divider: N² = Z_L/(Z_L+Z_s). (Will, pre-existing; now tested once and passed)
- NS: the medium co-rotates with the electrode (no-slip).
- RIDE / RIDE-2: an orbit rides the medium; the seat's rulers are material lines of the medium.
- P11-ratio: cosh l = radial/transverse ruler length ratio (backed by THM-I's ρ_K = cosh²l).
- KMS: the period of the imaginary loop is ħ/(k_B T). The residue.

## Witness disclosure
The DYN-1 shell assignment was found BACKWARD in chat (both tried, match kept) before the runner fixed it by declaration. The
runner header records this. Everything after DYN-1 was run forward.

## Order of what is next (critical first)

0. **Will:** freeze PRIMITIVES-v0 (P1–P13). Line zero.
1. **Will:** re-run the seven on the 7800X3D; commit each with its .log and the draft receipts (in chat log / below).
2. **Docs:** LABELLED-MODEL.md edits (five items in LOCK5-STATIONS); CONJECTURE-COSMOLOGY header (all r-language is LOCK-1;
   kill 2 survived; sheet = reading); HANDOFF: this file. Paper v0.6: scaffold entry BARE-1/Cl(3) + the pinning re-tiered.
3. **Kill 3 (signalling)** — a cheap runner with a "dies immediately" outcome; must run before any further cosmology.
4. **The AC sector = the 377 Ω chase, reframed.** The medium so far is DC: resistive only. A purely resistive continuum carries no
   waves; the model must carry gravitational waves at c (GW170817: |c_gw − c|/c < 10⁻¹⁵) and must have a light cone. That needs
   inductance and compliance — the L and C of V = IZ — and the wave impedance of such a line is √(L/C). Damour's horizon
   resistivity 377 Ω = Z₀ is the matched termination of that line at the electrode. So "chase 377 Ω" and "give the medium
   waves" are ONE task, and it is THM-K's scalar half (ζ = e²Z₀/ħ = 4πα). External target: reproduce Damour 1978 / Znajek 1978.
   Header must be written blind; this is the most seductive item on the list and therefore the one most likely to be steered.
5. Remaining edge of the octahedron: l = 0 below-barrier tunnelling at the electrode.
6. Unpin K: the medium's angular-momentum content. Compare the O(J²) exterior against Hartle–Thorne (fluid body), not Kerr.
7. Structural runners from the afternoon, unrun: ħ- and G-seat readings of the ground (six readings, one surface); lit face on the
   octahedron vertices; the η → 1/η involution as an object on the Gram (is it a rotor?).

## Draft commit receipts
- `DYN-1 ohmic divider 16/16: lapse derived from Ohm+T8c+THM-I; beta=gamma=1 read off output; Mercury 42.9807 by two paths; d=2 control fails; U^3 = -3/2 PREDICTION`
- `EQ-1 tolman-from-divider 12/12: T*N = T_inf by mode balance and Nyquist balance; no static equilibrium inside (theorem); horizon T->oo; dual sheet swaps the hot end (flagged)`
- `DYN-2 electrode period 12/12: kappa=1/(2r_s), Rindler form, beta=4 pi r_s/c derived from divider + T4b/KVL; T_H recovered given KMS (the one declared line); LOCK-3 = LOCK-2; dual sheet cold (flagged)`
- `DYN-3 stokes swirl 10/10: frame dragging derived as medium swirl (Laplace + no-slip + Faxen); dipole pattern and node/gyro = 4 structural; K pinned to LAGEOS predicts GP-B 39.1 vs 37.2±7.2; swirl profile uncorrected at any order (Kerr O(J^3)); O(J^2) quadrupole open`
- `SHEET-1 block+holonomy 17/17: sheet swap unique; I=-1 reverses interior clock, kill 2 does not fire; fluid path D=(r0/r)^{3/2} never zero, sheet inherited; theta_+=0 at river=c, theta_-=0 on I=-1; white sheet is a reading not a destination`
- `FOLD-1 light-ring dressing 10/10: lambda_L = Omega_c (two paths); subring ladder e^-pi (PREDICTION); eikonal QNM = H flip; greybody edge width lambda_L != kappa; face H derived at the fold`

## Working with Will — additions from today
- The circuit was right five times in one session (divider, negative-resistance interior = Damour's negative bulk viscosity,
  KVL as Euclidean regularity, Faxén for gyroscopes, matched load at η = 1). Treat a circuit statement as a hypothesis with a
  strong prior; verify; expect to lose.
- "The singularity is a seat artefact" was corroborated by SHEET-1: the fluid meets D → ∞, the outside seat draws D = 0.
- Will's caution ("no contradiction found" ≠ "no test that could find one") was the right call; the register came from it.

---

# AMENDMENT — 2026-09-08 (afternoon/evening, Will at the workstation)

This section is appended to `HANDOFF-2026-09-07-evening.md` and makes it the live state through 09-08.
All runners below were executed on the workstation and committed with their logs; nothing is sandbox-only.

## What happened

The 09-07 register left three things owed on the physics side: the O(J²) discriminator, the l = 0 edge of the
octahedron, and the tensor Love number. All three ran. The chain produced one prediction with an instrument, closed
the last edge of the octahedron, and — in the course of the rotating case — **retracted one identification and replaced
one constitutive reading**. The retractions are the most important part of this amendment; they are at the front of
`docs/results/2026-09-08/2026-09-08-EDGE-chain.md`.

## Runners added (all on the workstation, all committed with logs)

| runner | checks | commit | result |
|---|---|---|---|
| `prim_ns1_swirl_inertia.py` | 10/10 | `300eb99` | finite-Re swirl = K′/r³ + εK/(4r⁴); GP-B/LAGEOS bound ε < ~1.3×10⁴ km ⇒ the measured 1/r³ drag forces the creeping regime; in it the O(J²) exterior has **no P₂/r³ term: κ = 0** (Kerr 1); electrode spherical, ergosurface oblate; δκ_s = −1 inside GWTC-4.0's hierarchical 90 % interval |
| `prim_edge1_electrode_absorption.py` | 7/7 | `fcf183b` | σ_abs(ω→0)/A_H = 0.9993 (DGM/Page as OUTPUT); \|T₀\|² = 4(ωr_s)²; equals the DC capture area πb_mb² (GUD-1); static l=2 regular branch polynomial ⇒ scalar k₂ = 0. **Face H closed at both ends** |
| `prim_edge2_tensor_love.py` | 7/7 | `f8e054a` | GR's tensor Love number computed as labelled COMPARISON: static even-parity l=2 verified to be associated Legendre (2,2), regular branch polynomial, Hinderer → 0. Same zero, same mechanism |
| `prim_edge3_rotating_dissipation.py` | 12/12 | `796ed88` | **R-1** retraction of EDGE-2 C1; **R-2** the viscous swirl fails the dissipative sector; **R-3** MAT-1 → MAT-2 (swirl as vector potential), every earned result preserved |
| `prim_edge4_tidal_friction.py` | 8/8 | `eac22cc` | membrane tidal friction: Σ_m m²\|A_m\|² = (8π/15)[E₁ − (3/2)E₂] = **exactly Poisson 2004 Eq. 9.39's small-spin structure**; both zeros; first law closes; coefficient off by 32 |
| `prim_edge5_spin_ratio.py` | 7/7 | `cc60d99` | spin-0 and spin-2 l=2 absorption on the same background: 1/2025 and 4/225 as OUTPUT; **ratio 36.006** = Page's factorial. **R-4**: EDGE-4's 1/(512π) is Damour's 1/(16π) seen through a spin-0 coupling |

Ledger: EXT-034 (Moncrief), EXT-035 (Hinderer/Binnington–Poisson, comparison), EXT-036 (induction/MAT-2),
EXT-037 + amendment (Poisson Eq. 9.39, **fetched not remembered**), EXT-038 (Regge–Wheeler + Page, comparison).

## State changes

- **MAT-1 → MAT-2.** The medium's constitutive reading changed: the swirl is inductive, not viscous. The bulk is
  inviscid (river: potential flow, Euler with p = 0); dissipation lives at the electrode. Every DYN-3/NS-1 result
  survives because the *equation* was right and only the constitutive word was wrong. Record this as the second
  cautionary case after Cl(3), with the opposite moral: an equation reached from the primitives survived a change of
  interpretation that a result reached *for* would not have.
- **The octahedron is closed.** Face H has both ends: the electrode's power law \|T₀\|² = 4(ωr_s)² (EDGE-1) and the
  fold's Fermi edge of width λ_L (FOLD-1), with DYN-2's electrode temperature as the thermal factor between them.
- **The fingerprint** (κ, k₂, σ₀/A_H, tidal friction) = (0, 0, 1, ~1/36 of Kerr). Two BH-like, two not, none ECO-like.

## Method notes worth keeping

- **Fetch, don't remember.** EDGE-4's first draft quoted Poisson's coefficient from memory. Fetching it changed the
  comparison from "an O(1)" to "a factor 32", which is what made EDGE-5 findable. The ledger records which numbers were
  fetched.
- **Structure matching with a coefficient mismatch is a signature, not a failure.** Exact tensor structure + exact zeros
  + wrong amplitude ⇒ correct kinematics, wrong coupling. That reading found the spin conflation in one runner.
- **Same integrator, two fields.** EDGE-5's design — compute both coefficients with one code so normalisation cancels —
  is the pattern to reuse whenever a ratio is the claim.

## Next, in order

1. **η_s from the medium itself** (makes prediction 2 unconditional). Two routes, both named in the results page; the
   resistivity route *is* the 377 Ω chase / THM-K's scalar half.
2. **Does the AC sector's frame twist (T8) give the medium a spin-2 tidal channel?** Decides whether prediction 2 stands.
3. Kill 3 of CONJECTURE-COSMOLOGY (signalling) — still owed, still cheap, still ahead of any further cosmology.
4. NS-1's open items: finite-Re meridional flow at O(J²); the electrode's shape as a derived boundary condition.
5. The 09-07 structural runners, still unrun: ħ- and G-seat readings of the ground; lit face on the octahedron
   vertices; the η → 1/η involution as an object on the Gram.

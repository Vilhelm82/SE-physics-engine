# THM-TARGET O — the strain law: Gauss and Codazzi on the seat's flat rods, at THM-N's gate
Written 2026-09-04, early.  `thm_o_strain_law.py` **33/33** (`714bf5c`, `run_o.log`, 114.7 s on the workstation).  Read after `docs/2026-09-03-thm-target-N.md`.
Re-running is the verification; each check flushes as it lands, so `tail -f run_o.log` shows progress.

## In one paragraph
The handoff asked for the closure THM-N demands — e₂(strain) = 0 — to be derived as the unique admissible invariant from a model
principle before use, or it is ADM in a hat.  It is derived, and the derivation is one line: one boost per seat means the rain seat's
rods are flat, and Gauss's theorema egregium says the normal-normal Einstein component of any 4-geometry containing a flat slice with
second fundamental form K is exactly e₂(K).  Under a stationary shift v with unit lapse, K is the symmetric strain of v.  Codazzi gives the
rest the rods can see: G_ni = −½ curl curl v.  "No source through the seat's normal" is therefore e₂(sym ∇v) = 0 and curl curl v = 0 —
not one of nine generators but the only invariant flat rods present.  The identity is proved; the vanishing in vacuum is declared
(NORMAL-1).  The pinning is derived a second time (exponent only; the coefficient needs MASS-1), c₂ = 0 is forced, Lense–Thirring's
1/r³ exponent follows from Codazzi (the amplitude is still THM-M's), and the first-order drag on flat rods is exact vacuum in all ten
components.  At second order in the spin the closure is a determined linear problem: Codazzi makes the correction a gradient, Gauss is
the anisotropic operator 4∂²ᵣ + (2/r)∂ᵣ + Δ_S/r² sourced by the drag's strain-squared, and its exterior exponents are (ℓ+1)/2 and −ℓ/2 —
half of Laplace's.  Only ℓ = 0 (mass) and ℓ = 1 (a boost) coincide with Newton.  **Gate 1 passes by structure: no P₂·r^{1/2} mode
exists, so no P₂/r term can arise.  Gate 2 is not reachable: the orbital invariant ½(Ω₊²+Ω₋²)(R_c) has no R⁻⁵ slot, where Kerr carries
3J²/(2M); the model's ℓ = 2 freedom is a half-integer tail R^{−9/2}.**  The tangential Einstein components fail at O(ε²) for every value
of the two free constants.  Fork (b): flat rods fail at second order in the spin; the multipoles of any source live in the bending of
the rods, and the bending is rank two.

## The chain, with tiers
| step | statement | tier |
|---|---|---|
| O-1a | the presentation −dt² + \|dx − v dt\|² has unit lapse, flat rods, n·n = −1, det g = −det h | proved |
| O-1b | **Gauss:** G_nn = e₂(K), K = sym ∇v, identically for three arbitrary functions v(r,θ) (theorema egregium, R(h) = 0) | proved |
| O-1c/d | **Codazzi:** G_ni = −(D_jK^j_i − D_iK) = −½[curl curl v]_i identically (Cartesian route for the second equality, any v(x,y,z)) | proved |
| O-2a–e | of e₁, e₂, e₃ of the strain cell only e₂ = 0 returns rβ² = C (e₁ = 0 → β ∝ r⁻², e₃ = 0 → β′ = 0, both dead); c₂ = 0 forced | derived \| NORMAL-1 + flat rods (second route to the pinning, disjoint from THM-K) |
| O-3a–d | Codazzi ⇒ rw″ + 4w′ = 0 ⇒ w = A/r³ (two routes); the O(ε) cross term in e₂ vanishes for any β, w; e₂(S₁) = −(rw′)² sin²θ/4 is the second-order source | proved (exponent); amplitude A = 2GJ/c² remains THM-M's |
| O-4a | poloidal axisymmetric v₂ with curl curl v₂ = 0 and a regular axis is a gradient (ω_φ = C/(r sinθ) forces C = 0) | proved |
| O-4b/c | polarised Gauss: −(β′+β/r)∇²ψ + (β′−β/r)ψ_rr for any β, ψ; with the pinning, −(β/2r)[4f″ + 2f′/r − ℓ(ℓ+1)f/r²]P_ℓ (ℓ = 0..3, spherical route) | proved |
| O-4d/e/f | **exponent theorem:** p ∈ {(ℓ+1)/2, −ℓ/2}; the Newtonian clock multipole needs ψ ∝ r^{1/2−ℓ}, in the spectrum for ℓ = 0, 1 only; no P₂r^{1/2} mode | proved |
| O-5a/b | exterior: ψ₂ = a₀√r − (A²/10√r_s) r^{−5/2} + [b₂/r + (A²/8√r_s) r^{−5/2}]P₂ satisfies Gauss and Codazzi through O(ε²) for every a₀, b₂; the particular constants are forced (explicitly perturbed field leaves an O(ε²) residual) | proved |
| O-5c | clock at O(ε²): −√r_s a₀/r, 2√r_s b₂P₂ r^{−5/2}, (A²/6 − A²P₂/24)/r⁴; **no P₂/r, no P₂/r³** | proved |
| O-5d | orbital invariant ½(Ω₊²+Ω₋²)(R_c): r_s/(2R³) − √r_s a₀/(2R³) − (5/4)√r_s b₂ R^{−9/2} + **0·R⁻⁵** + (7/8)A² R⁻⁶ | derived \| transport = geodesic (THM-I/L) |
| O-6a/b | the presentation is vacuum in all ten components at O(ε⁰) and — new — at O(ε¹): flat rods + curl-free drag is an exact first-order field | proved |
| O-6c | at O(ε²) the components the rods fix vanish (corollary of O-1b/c + O-5a) | proved |
| O-6d | at O(ε²) G_rr, G_θθ, G_φφ ≠ 0; a₀ drops out exactly (a mass shift is still vacuum); the b₂ that kills G_θθ at one point leaves the others nonzero | proved (corroboration under Einstein's tangential law; the rods-bend conclusion does NOT rest on it — see 4 below) |
| O-7a–c | independent numeric path (mpmath, nested central differences, 40 digits): e₂ = O(ε⁴) with the forced constants (ratios 16.0), O(ε²) with them doubled (3.99); curl curl v = 0 to 10⁻¹⁵ | numeric corroboration |
| C-1/2 | **Kerr:** same invariant, same variable: M/R³ + (3J²/2M)R⁻⁵ + 6J²R⁻⁶.  Gate 1 passed; Gate 2 not reachable; fork (b) | comparison |

Declared: NORMAL-1 (no source through the seat's normal in vacuum).  Conditional throughout on one boost per seat (flat rain rods,
unit lapse) and BARE-1/M-2 (the presentation metric is the η-sandwich of the boosted frame).

## Failure or disagreement — sorted by the model's own standard
The criterion: a *failure* contradicts a model theorem, a number measured in a model-owned observable, or the model's own axioms; a
*disagreement* is a definite, integer-structured prediction differing from GR where GR is unmeasured.  GR is never the arbiter.
1. **Gate 1 (P₂/r): failure of THM-N's construction, internal.**  K-6 puts ℓ = 2 at 1/r³; MASS-1 makes the 1/r term the scalar blind
   mass.  The closure restores consistency; no import.
2. **The exponent theorem: failure of flat rods, internal + ground.**  The half-integer tail contradicts K-6's integer multipoles
   (K-6 is a scaffold — linearity declared) and, through the same operator on a static oblate source, J₂'s measured a^{−7/2} scaling
   (the ground: on flat rods the ℓ = 2 clock would fall as P₂r^{−5/2} and the nodal regression as a⁻³).  The model's static rods for a
   non-spherical source are already a sum of rank-ones (THM-M), not flat — so this fails the *formulation*, not the model.
3. **The R⁻⁶ coefficient (7/8 A² vs Kerr's 6J²): neither yet.**  Superseded when rods bend; not comparable.
4. **The rods bend at O(J²) — rests on O-4e + K-6 + J₂, not on O-6d.**  O-6d adds only: Einstein's tangential law also rejects flat rods.
5. **The tangential sector: where the model is allowed to disagree.**  Gauss and Codazzi are the complete *normal* projections of the presentation (G_nn, G_ni) (O-6c); the tangential G_ij are computable from (h, v) in the stationary case (O-6d) but no model law owns their value.  The presentation restriction h = δ is observer-side, not an ontological claim (see `docs/2026-09-04-axis-planes-rods-clarification.md`).  GR has a law for the six components; the model is silent, not wrong.  Silence, not error.
6. **Kerr's Q = −J²/M: not yet a legal comparison.**  Flat rods produce no quadrupole at all.  A tangential law will produce some κ_model
   (Q = −κJ²/M); integer-structured and κ ≠ 1 is PREDICTION-2 with kill numbers (LIGO spin-induced-quadrupole bounds, order unity to
   order ten; ringdown no-hair tests ~10–20%).  A hole's Q also needs regularity at the wall (THM-H at O(J²)); a shell's needs source
   matching and is compared with GR's shell, not Kerr (de la Cruz & Israel 1968; Pfister & Braun 1985).
7. **Rank two (conjecture, labelled).**  The O(J²) bending and the TT wave (handoff item 4) are both rank-two objects one boost cannot
   present.  A second boost gives PSD deformations; the ħ squeeze gives trace-free indefinite ones.  Whichever the tangential law
   demands is which pivot the model has.  Kill: if the demanded deformation is PSD, the squeeze is dead before GW-1.

## Band, named loudly; edge, in the model's words
A reproduction is a band of the spectrum; it is worth recording only with its edge.  Three fields per reproduced result: the band, the
edge (model-internal), the position (which seat principle puts it there).  A special-case identification has a fourth: the surplus.
**Band:** O-1 on flat rods with unit lapse IS the ADM Hamiltonian–momentum constraint pair in Painlevé–Gullstrand gauge; O-2 is
Schwarzschild in PG form; O-3/O-6b are the Doran–Hamilton–Lisle river at first order in J.  **Edge (no GR name):** flat rods carry a
mass and a boost and no other multipole, for any source (O-4e).  Garat & Price 2000 (PRD 61, 124011) and Valiente Kroon 2004 are the
shadow of that edge, narrower and about Kerr only.  **Position:** one boost per seat.  **Surplus:** none yet — the general case needs
the tangential law.  Two-thirds of a special-case identification.

**Scaffold register** (ground / declared / scaffold, with removal route):
- NORMAL-1 — declared; Einstein's constraint renamed.  Removal: sourced Gauss from the blind-mass theorem (MASS-1): if the pivot-blind
  moment is what passes through the seat's normal, G_nn = 8πρ_blind falls out and the beam becomes structure.
- K-6 linearity — declared; the integer-multipole kill of the tail rests on it.  Ground beneath it: J₂'s a^{−7/2}.  Removal: the
  tangential law, whatever it is.
- Drag amplitude 2 = 1 + 1 — THM-M superposition.  Removal: Codazzi sourced by the shell's momentum density (one section, not run).
- EQ-1 Tolman — Mercury is the ground, Tolman the coat; ruling still owed.
- Ground, never removed: Mercury 3, GP-B 37.2, J₂ scaling.
- The tangential sector has no scaffold because it has nothing: the empty lot.

## Gate 2, re-declared
THM-N's Gate 2 as written — Q = −J²/(Mc²) *from a rotating shell* — is Kerr's no-hair number.  Two gates replace it: **(i) hole** — Q
fixed by regularity at the wall, THM-H at O(J²), compared with Kerr; **(ii) shell** — Q fixed by the source's O(ε²) energy and stress,
compared with GR's shell.  Neither arises on flat rods.

## Process note
The suite was built in the assistant's container and run there first; that container is a 4 GB sandbox whose processes are reaped
between calls, and the run took forty minutes and three kills.  On the workstation it took 114.7 s.  Rule for the next instance:
the container is for syntax checks; every real run is on the workstation, tee'd to a log the operator can tail.

## Next, in order
1. **Sourced Gauss from MASS-1** — NORMAL-1's removal route.  Cheap; a scaffold comes out or the beam is named load-bearing.
2. **Sourced Codazzi for the drag amplitude** — does the shell's momentum density return 2 = 1 + 1?  Removes K-6 from the drag.
3. **Tangential-law candidates on a STATIC OBLATE source at O(ε⁰)** before any spin: THM-I's reciprocal presented volume Δ/det G′ as
   harmonic variable; THM-K's Tolman equilibrium on a rotating screen; integrability of rule G′.  Any candidate that does not return
   Poisson's integer multipoles there is dead for a dollar.  Survivors go to O(J²): (a) Einstein's tangential equations — reproduction,
   named loudly; (b) something else, integer-structured — PREDICTION-2 with κ_model; (c) something else failing J₂ or Mercury — failure,
   naming which principle to drop.
4. **BEND-1**: the presented Gram at O(J²) — which deformation G′ → G′ + ε²D the tangential law demands; second boost (PSD) vs squeeze
   (trace-free indefinite).  Shares the pivot with GW-1.
5. ħ-seat reading as queued (handoff item 5); the floor A₀² as a derived number rather than a unit is the target Will named.

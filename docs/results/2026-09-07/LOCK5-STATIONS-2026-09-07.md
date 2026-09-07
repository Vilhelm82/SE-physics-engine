# LOCK-5 — STATION STATUS RELABEL — 2026-09-07

Directive: every station of the octahedron carries a GROUND / THEORY / DERIVED tag. The lattice geometry
(8 faces = 2^{C,H,G}, null/boost classes, antipode = deck) is import-free and was re-run today before any label
was touched: `label3_lattice.py` 14/14, `label3b_corrections.py` 10/10, `label2_octahedron.py` 11/11 (sandbox).
Nothing below changes a face. It changes what each face is *allowed to claim*.

Tags:
- **GROUND** — a measured number exists at the station, with a precision.
- **DERIVED** — reached from the primitives / the divider today or earlier, runner cited.
- **DERIVED-given-X** — derived except for one named declaration X.
- **THEORY** — a formula with a name and no measurement; scaffold.
- **DECLARED** — a hypothesis in the list, by choice.

---

## Poles

| pole | axis | tag | receipt |
|---|---|---|---|
| Light | c+ | GROUND | c exact by SI definition; light propagation measured everywhere |
| Temperature | c− | structure DERIVED · identification DECLARED | Wick face, eighth-turn, cover reached (T4b). "Period = temperature" is the KMS line (DYN-2) |
| Action | ħ+ | GROUND | ħ exact by SI; wave phase e^{iS/ħ} measured (interference) |
| Evanescence | ħ− | GROUND | tunnelling measured (α-decay lifetimes, STM, Josephson) |
| Mass | G+ | GROUND | G to 2×10⁻⁵ (CODATA); binding measured |
| Energy | G− | GROUND | E = mc² via mass defect, ~10⁻⁷ |

## Antipodal relations

| axis | relation | tag | receipt |
|---|---|---|---|
| c | Unruh | DERIVED-given-KMS | DYN-2 d1: EQ-1 + derived local acceleration recover ħa/(2πk_Bc) at the electrode. No ground (unobserved). LOCK-3 merged into LOCK-2 |
| ħ | S → iS_E | DERIVED (structure) · GROUND-adjacent (rates) | T4 eighth-turn; WKB tunnelling rates measured |
| G | dE = −c² dM | GROUND | mass defect |

## Faces

| flips | physics | tag before 09-07 | tag after 09-07 | receipt |
|---|---|---|---|---|
| ∅ | our face: light waves in bound orbits, curved spacetime | established | **GROUND + DERIVED** | PPN β = γ = 1 read off the divider output; Mercury 42.9807″ two paths; deflection 1.7512″ (DYN-1 16/16) |
| G | redshift, Shapiro, escaping GW | named | **GROUND + DERIVED** | first-order lapse −2U from the divider (DYN-1 d3); redshift measured to few×10⁻⁵; Cassini γ; GW detected |
| H | near-horizon greybody filter (G-2) | named | **THEORY**, with a DERIVED geometric-optics limit | full greybody unmeasured. The geometric-optics capture edge — light ring at 3r_s/2, escape cone sin ψ = (3√3/2)ηN³, its fold at P₂(cos θ)=0 — is DERIVED from the divider (DYN-1 f1 + 09-07 chat). EHT shadow consistent at 10–17 % |
| C | Tolman–Ehrenfest (EQ-1) | named | **DERIVED** | T·N = T_ref by mode balance and Nyquist balance (EQ-1 12/12). Ground-adjacent: the clock gradient −g/c² is measured; the thermal step is definitional |
| HG | Parikh–Wilczek tunnelling emission | named | **THEORY** | unmeasured; no model derivation yet |
| CG | Planck's face: blackbody | named | **GROUND** | blackbody / CMB spectrum to ~10⁻⁴. The one far-column station the world has measured |
| CH | finite-T instantons; stellar fusion | named | **GROUND (rates) · THEORY (formalism)** | fusion rates measured (Gamow peak, solar neutrinos); instanton machinery is scaffold |
| CHG | Hawking radiation | established | **DERIVED-given-KMS · PREDICTION** | DYN-2 12/12: κ = 1/(2r_s), β = 4πr_s/c derived; T_ref = ħc/(4πk_Br_s) given the one KMS line. No ground; 6.17×10⁻⁸ K for the Sun is a prediction |

---

## What the relabel says

**Earned today.** Three far-side faces moved: C (Tolman) is DERIVED; CHG (Hawking) is DERIVED-given-KMS; the
geometric-optics edge of H is DERIVED. All three came off the same tool — the divider with the horizon as electrode.

**Still scaffold: exactly two faces, and they share an edge.** H and HG are the two faces containing *Light + Evanescence*
— light tunnelling near the horizon. The greybody and Parikh–Wilczek are the same physics seen with and without the Energy
pole. The model owes *one* derivation, not two: how light tunnels at the electrode. Everything the model has said about
that edge so far is geometric optics (the fold). The wave part — Airy fringes on the fold, in the rainbow's language — is
the open item, and it is the natural next use of the eighth-turn: T4b's Sym² relation between the two-sided (wave) and
one-sided (particle) readings is precisely the H flip.

**The "unique conjunction" sentence.** "Hawking radiation is the unique conjunction thermality + tunnelling + escape" is
PROVED as lattice geometry (CHG = −f₀, the deck) and DECLARED as physics. After DYN-2 the physics half is one KMS line
away from derived. The sentence keeps both tags.

**The residue, in one place.** Across LOCK-2, LOCK-3, LOCK-5 there is now a single undischarged declaration: *the period of
the seat's imaginary-time loop is ħ/(k_BT)*. Every "temperature" on the far column inherits it. It is named on the
LABELLED-MODEL page, once, and nowhere else.

---

## Edits owed to `docs/LABELLED-MODEL.md`
1. §2 poles: add the tag column above.
2. §4 faces: replace the `status` column with the "after 09-07" column and receipts.
3. §4 footnote: replace "Three of the named faces are things the paper already touches…" with the paragraph
   "Still scaffold: exactly two faces…" above.
4. Add a **Residue** line under §4: the KMS declaration, cited to DYN-2.
5. Runner list: append `prim_dyn1_ohmic_divider.py` (16/16), `prim_eq1_tolman_from_divider.py` (12/12),
   `prim_dyn2_electrode_period.py` (12/12) — after workstation re-run.

Note kept from `label3b_corrections.py` for the next runner: D ≠ 0 at both horizons (1.1752, 4.2880), so the sheet at a
horizon is the path/block holonomy χ[Γ], not a sign flip of D. The sheet-sign computation is a holonomy computation.

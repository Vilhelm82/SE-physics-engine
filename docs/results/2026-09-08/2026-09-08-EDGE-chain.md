# 2026-09-08 — The EDGE chain: the electrode's response, and the fingerprint

Runners: `suites/prim_edge1_electrode_absorption.py` (7/7), `prim_edge2_tensor_love.py` (7/7),
`prim_edge3_rotating_dissipation.py` (12/12), `prim_edge4_tidal_friction.py` (8/8), `prim_edge5_spin_ratio.py` (7/7).
All run on Will's workstation; logs committed beside each runner. Ledger: EXT-034 … EXT-038.

---

## 1. Corrections and retractions (front of the page, per the tiering rule)

**R-1 (EDGE-3A retracts EDGE-2 C1).** EDGE-2 identified the medium's static l = 2 *acoustic* growing branch with GR's
tensor tidal field. Wrong. That branch changes the medium's potential Φ = −|v|²/2 as r^{1/2}P₂, not r²P₂: it is a strain
**flow** at infinity, not a Newtonian tide. EDGE-2B (GR's tensor k₂ = 0, computed as comparison) stands untouched.
The medium's actual tidal response is EDGE-3A's: a tide enters through Bernoulli, the induced exterior multipole is
identically zero at every order, and the sonic surface deforms. k₂ = 0 becomes *trivially exact* rather than computed.

**R-2 (EDGE-3B retires MAT-1's constitutive word).** DYN-3 derived the swirl *viscously* (Stokes + no-slip). A viscous
medium torques the rotating electrode at −8πμa³Ω — it spins holes down with no tide, on t = M/(8πμr_s) — and, creeping
flow being linear in its boundary data, admits **no tidal torque at all** (ε_ijk E_jk = 0 for symmetric E). Kerr is the
opposite on both counts. The viscous reading fails the rotating dissipative sector.

**R-3 (MAT-1 → MAT-2, the repair).** The swirl is the medium's **vector potential**: a dipole's A_φ = m sinθ/r² is
DYN-3's field exactly, with m = Ωa³. ½∇×A is the gyroscope field (Faxén's ½ = GEM's ½; the dipole pattern is Ampère's);
the field energy ∝ 1/r⁶ so NS-1's κ = 0 survives; and the river advects it by a steady induction equation whose
φ-component is *the same operator* as NS-1a's azimuthal momentum equation with diffusivity for kinematic viscosity —
so the 1/r⁴ tail and the GP-B/LAGEOS bound survive. Nothing earned is lost. Dissipation moves out of the bulk and into
the electrode's resistivity, which is where the membrane paradigm always had it.

**R-4 (EDGE-5 explains EDGE-4's coefficient).** EDGE-4 needed a membrane viscosity of 1/(512π) = Damour's 1/(16π) ÷ 32.
That was a spin conflation: the medium's tide is spin-0, Poisson's formula is for a spin-2 tide, and on the same horizon
a spin-2 l = 2 field is absorbed **36×** more strongly. Corrected, the membrane needs 1.125 × Damour's coefficient —
inside EDGE-4's two declared O(1)s. **1/(512π) is not a model number; Page 1976's factorial owns it.**

---

## 2. The electrode's response, derived

**EDGE-1 — absorption.** The medium's sound in the divider metric (NS-0), with the ingoing condition derived from the
cone stalling at the sonic point (SHEET-1 d1). Circuit readings of that one boundary: **short** in the DC divider
(Z_L → 0), **matched** in AC, **open** from inside. Results:
- σ_abs(ω → 0)/A_H = **0.9993** (extrapolated); |T₀|² = 4(ωr_s)². Das–Gibbons–Mathur / Page as *output*.
- Flux |R|² + |T|² = 1 to 10⁻⁸ at every frequency.
- The approach is linear in ω and from **above**: 1.48 → 1.20 → 1.09 → 1.04 as ωr_s halves from 0.1.
- **The AC absorption area equals the DC capture area** πb_mb² = 4πr_s² (GUD-1). Long waves and slow particles read the
  electrode as the same matched load of the same area.

**EDGE-1c / EDGE-2 / EDGE-3A — the Love number.** k₂ = 0, three ways: the acoustic regular branch is a polynomial
(P₂(2r−1), no r⁻³ tail); GR's tensor branch is a polynomial too (3(r²−2Mr)/M², Hinderer's formula → 0 at R → 2M),
verified to be the associated Legendre equation (2,2) rather than assumed; and after R-1, the medium's response to a
Newtonian tide is identically zero by Bernoulli. Same zero, same mechanism: a polynomial regular branch, a log-singular
tail-carrying branch, regularity at the boundary killing the tail.

**EDGE-4 — tidal friction.** The rotating electrode as a membrane: the tidal deformation of the sonic surface rotates
under it, area preservation at first order forces a tangential flow, the flow shears, the shear dissipates, J̇ = −P/Ω.
- Structure: J̇ ∝ −Ω E² r_s⁶, weighted by Σ_m m²|A_m|² = (8π/15)[E₁ − (3/2)E₂].
- **That invariant combination is exactly Poisson 2004 Eq. (9.39)'s small-spin structure** — from the angular integral of
  a rotating quadrupole on a sphere, not put in.
- Zero for a tide axisymmetric about the spin; zero at Ω = 0; first law closes (dM = 0, T Ṡ = −Ω J̇ = P > 0).
- Bochner identity used in the shear integral verified explicitly, not trusted.

**EDGE-5 — the spin ratio.** Same background, same integrator, same boundary condition, two fields:
spin-0 l = 2 → 4.959×10⁻⁴ (Page 1/2025 = 4.938×10⁻⁴); spin-2 l = 2 (Regge–Wheeler, comparison) → 1.785×10⁻²
(Page 4/225 = 1.778×10⁻², the classic 256/225 (Mω)⁶). **Ratio 36.006**, against Page's 36.

---

## 3. The fingerprint

| quantity | medium | Kerr | ECOs | status |
|---|---|---|---|---|
| spin-induced quadrupole κ | **0** | 1 | ≠1 (NS 2–14; boson stars 10–150; gravastars <0) | NS-1, DERIVED |
| tidal Love number k₂ | **0** | 0 | ≠0 | EDGE-3A trivially exact; EDGE-2B for GR |
| σ_abs(0)/A_H | **1** | 1 | ≠1 | EDGE-1, DERIVED |
| tidal friction | **~1/36 of Kerr** | 1 | varies | EDGE-4+5, given Damour's η_s |

Two BH-like entries, two not, none ECO-like. The GWTC-4.0 paper's expectation — that κ ≠ 1 implies nonzero tidal
deformability — is exactly where the medium breaks from the ECO profile: **κ = 0 with black-hole tidal properties is a
combination nothing in that literature produces.**

---

## 4. Predictions, with instruments

1. **κ = 0 → δκ_s = −1.** GWTC-4.0: hierarchical 90 % [−53, +9], restricted [−28, −2]. The model and Kerr differ by
   **one unit on a ruler whose tick is 26–62**; current data cannot distinguish them. Needs σ(δκ_s) ≲ 0.3 — a factor 30–100.
   Named movers: GW241011 (κ₁ free), GW250114 (loudest), 3G/LISA.
2. **Tidal friction ~1/36 of Kerr's**, if the medium's tide stays spin-0 and the membrane viscosity is Damour's.
   Instrument: LISA EMRI horizon-flux phasing (tens of radians in GR → ~1); LVK horizon-absorption tests, currently
   unconstraining. *Escape route, open:* the medium may acquire a spin-2 tide through the AC sector's frame twist (T8).
3. **Λ = 0 exactly** for the medium's holes (BH-like, not ECO-like). No BH Love number has ever been measured.
4. Carried from 09-07: subring ladder e^{−π} per half-orbit; λ_L = Ω_c; greybody edge width λ_L ≠ κ; no swirl correction
   at any order (Kerr has O(J³)); isotropic g₀₀ at O(U³) = −3/2.

---

## 5. Open, named

- **η_s from the medium itself.** EDGE-4 treats the surface viscosity as a free transport coefficient. Until the medium
  derives it, the tidal-friction prediction is conditional on Damour's value. Two routes: the electrode's own resistivity
  (the 377 Ω chase, THM-K), or the horizon's intrinsic deformation vs the sonic surface's (Poisson 2005; Poisson–Vlasov 2010).
- **Spin-2 in the medium.** Does the AC sector's frame twist give the medium a spin-2 tidal channel? If yes, prediction 2
  moves; if no, it stands.
- **Finite-Re meridional flow at O(J²)** (NS-1), and the **electrode's shape as a derived boundary condition** — an oblate
  electrode is the only way the medium recovers Kerr's Q, and the model does not currently produce one.
- **AREA** (first-order area preservation) and the shape→flow map are DECLARED in EDGE-4; they are the residual O(1).

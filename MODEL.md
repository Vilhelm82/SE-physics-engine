# The model

Two parts. **Part I** (§0–6) is the labelling of the octahedron — the geometry is proved, the names on it are a
labelling and are graded as one. **Part II** (§7–13) is the dynamics tier: what the model derives on that geometry,
with every claim tiered and every runner cited. Part II supersedes any status note elsewhere in the repo.

---

# Part I — The labelled octahedron

**Part I (§0–6) — the labelling of the geometry. Status: CONJECTURE tier throughout, except where marked.** The geometry it sits on is proved
(LABEL-2, LABEL-3, T5b, T5c). The *names* on the geometry are the labelling, and the labelling is
graded as a labelling: it has been checked against known physics (LABEL-1) and it has kill
conditions, but it is not derived. Written 2026-09-05 from the runner headers and the evening
handoff, so that it exists in one place.

The 08-17 labelled model (item codes c-1…, ħ-1…, G-1…, SECT-1) is superseded by this one where
they differ and is not in the repo.

---

## 0. What the octahedron is

Three lines through a root (P2), each with two ends. Six ends = six poles. Take the six poles as
vertices: they form an **octahedron**. Its 12 edges join poles on *different* axes. Its 8 faces are
the octants — each bordered by exactly one pole from each axis. Every edge borders exactly two
faces, and they differ only in the axis the edge does not touch. [PROVED, `suites/label2_octahedron.py`]

A face is what the seat sees looking into an octant: the plane perpendicular to that line of
sight, bordered by three poles. "Faces" here are these eight — not the six faces of a cube. [Will]

---

## 1. The three axes

| axis | near pole | far pole | near face | far face |
|---|---|---|---|---|
| **c** | Light | Temperature | GR | Quantum |
| **ħ** | Action | Evanescence | Wave | Particle |
| **G** | Mass | Energy | Bound | Escaping |

- A **pole** is the quantity at the end of the line. A **face** (per axis) is the physics you see
  looking at that end — what that pole *contributes* to every octant it borders. [Will]
- **Near** = toward the c seat, reachable by real operations. **Far** = behind the seat, reachable
  only through the imaginary angle. From c we see light propagate and matter in bound orbits; the
  far faces leak — tunnelling, Hawking. (G-4: sealed from c, leaky from ħ.)
- **c's near face is GR.** SR sits inside it as the zero-tilt case (T1–T3 derive SR; T7d derives
  the tilt; GR is the near face with $\eta \neq 0$).
- **ħ's faces are Wave/Particle**, not Classical/Quantum: action *is* wave phase ($e^{iS/\hbar}$);
  evanescence is where the tunnelling quantum is detected as one thing. Wave = the two-sided
  reading, Particle = the one-sided reading, and T4b′ says the vector reading is Sym² of the
  spinor reading. [labelling, with T4b′ as structural support]
- **G's poles are Mass/Energy.** From c they are identified by $c^2$, which is why the c seat
  cannot see G's pole structure — the seat collapses it, exactly as it collapses its own.
- **The bare tier cannot tell ħ from G** (T7b3, PROVED). Which of the two positive rulers is which
  is not in the geometry. See §6.

Read down the near column: GR, Wave, Bound, Light, Action, Mass — our physics. Read down the far
column: Quantum, Particle, Escaping, Temperature, Evanescence, Energy — **Hawking radiation in six
words.** [Will]

---

## 2. The six poles (vertices)

| pole | axis | end | what it is |
|---|---|---|---|
| Light | c | + | the null direction; what moves at the seat's constant |
| Temperature | c | − | the Wick face of c; the thermal circle (T4) |
| Action | ħ | + | real $S$; oscillatory phase |
| Evanescence | ħ | − | imaginary $S$; exponential decay, tunnelling |
| Mass | G | + | the source; what binds; what stays |
| Energy | G | − | what leaves; what radiation carries; what a hole loses |

Antipodes, and the known relation across each (LABEL-1, declared table):

| axis | antipodal relation |
|---|---|
| c | Unruh: the light-seat under acceleration reads a temperature |
| ħ | $S \to iS_E$: oscillatory ↔ decaying |
| G | $dE = -c^2\,dM$: what the bound pole loses the escaping pole carries |

**Poles are not symmetric.** Left and right halves of an axis are different identities; the c root
is one side of its antipole, not the same line with a sign. The pivot group is the identity
component $SO^+(2,1)$; the deck $V \to -V$ is *not* a symmetry of the model, it is the map between
identities. [Will; T5b makes it exact]

---

## 3. The twelve edges

Every edge joins two poles on different axes. Every one names a known relation, and its *kind*
matches its loop — imaginary counterparts in the two Wick loops, dimensional conversions in the
real loop. [LABEL-1, 13/13, checked against a DECLARED comparison table; graded as a lookup]

**Loop A — the (c, ħ) plane, hyperbolic, a Wick loop**

| edge | relation |
|---|---|
| Light → Action | photon phase $= S/\hbar$ |
| Action → Temperature | Euclidean action $= \beta E$; $e^{iS} \to e^{-S_E}$: QM ↔ stat mech |
| Temperature → Evanescence | thermal decay $e^{-\beta E}$; finite-$T$ instantons |
| Evanescence → Light | tunnelling out as radiation (Parikh–Wilczek) |

**Loop B — the (c, G) plane, hyperbolic, a Wick loop**

| edge | relation |
|---|---|
| Light → Mass | $E = mc^2$ with $E$ the photon energy; lensing |
| Mass → Temperature | $T_H = \hbar c^3 / (8\pi G M k_B)$ |
| Temperature → Energy | $E = k_B T$ |
| Energy → Light | $E = h\nu$ |

**Loop C — the (ħ, G) plane, compact, a real loop (the seat's space)**

| edge | relation |
|---|---|
| Action → Mass | $S = -mc^2 \int d\tau$ |
| Mass → Evanescence | Compton / Yukawa decay length $\hbar / mc$ |
| Evanescence → Energy | tunnelling rate $\sim e^{-2\kappa d}$, $\kappa = \sqrt{2m(V-E)}/\hbar$ |
| Energy → Action | $S = \int E\,dt$ |

**The loops.** A Wick step is an eighth-turn: $W^8 = 1$ on the vector reading, $W^{16} = 1$ on the
cover; even powers of $W$ are axes (poles), odd powers are planes. [PROVED, inline 09-04] Each loop
has eight stations; the four odd ones sit on the octahedron's edges, and an edge is labelled by
the axis the loop never touches — the seat's space showing up as what you are not looking at.
[PROVED, `suites/label2_octahedron.py`]

**Shared-station consistency.** Loops A and B both pass through Light and Temperature: the
temperature reached via ħ (Unruh) must equal the one reached via G (Hawking). **DERIVED — T4.**
Loops A and C share Action/Evanescence; B and C share Mass/Energy: those two conditions are
**OWED.** [LABEL-1]

---

## 4. The eight faces

Each face is a triple — one pole per axis — and its physics is the sum of its three borders.
The faces form the Boolean lattice $2^{\{C,H,G\}}$ of three flips: $C$: Light → Temperature,
$H$: Action → Evanescence, $G$: Mass → Energy. Ours is $\varnothing$; Hawking's is $CHG$, its exact
antipode $-f_0$. [PROVED, `suites/label3_lattice.py`, `suites/label3b_corrections.py`]

| flips | poles | contributions | physics (DECLARED candidate) | status |
|---|---|---|---|---|
| $\varnothing$ | Light, Action, Mass | GR + Wave + Bound | light waves in bound orbits around mass, in curved spacetime — **our face** | established |
| $G$ | Light, Action, Energy | GR + Wave + Escaping | light with phase leaving a well: gravitational redshift, Shapiro delay, escaping GW | named 09-04 |
| $H$ | Light, Evanescence, Mass | GR + Particle + Bound | light tunnelling in a bound system: near-horizon greybody filter (G-2) | **DERIVED 09-07 late, both ends**: fold's Fermi edge width $\lambda_L$ (FOLD-1) and sonic surface's $|T_0|^2 = 4(\omega r_s)^2$ (EDGE-1) |
| $C$ | Temperature, Action, Mass | Quantum + Wave + Bound | thermal waves around bound mass: Tolman–Ehrenfest (EQ-1) | named 09-04 |
| $HG$ | Light, Evanescence, Energy | GR + Particle + Escaping | light tunnelling out carrying energy: Parikh–Wilczek emission, pre-thermal | **THEORY** (09-07 late): no emission channel derived; H's closure does not cover it |
| $CG$ | Temperature, Action, Energy | Quantum + Wave + Escaping | thermal waves carrying energy, no mass: **blackbody, Planck's face** | named 09-04 |
| $CH$ | Temperature, Evanescence, Mass | Quantum + Particle + Bound | thermal tunnelling in a bound system: finite-$T$ instantons, stellar fusion | named 09-04 |
| $CHG$ | Temperature, Evanescence, Energy | Quantum + Particle + Escaping | thermal tunnelling carrying energy away — **Hawking radiation** | established |

The three single-flip faces isolate one ingredient each; the three pair-flip faces hold every
pair; the antipodal face holds all three. Hawking radiation is the unique conjunction *thermality
+ tunnelling + escape*, opposite *light propagation + action phase + bound mass*. [Will]

Three of the "named" faces are things the paper already touches without having known which face
they were on: the greybody (G-2), Tolman (EQ-1), Parikh–Wilczek. Planck's face sits one edge from
Hawking's — the blackbody measurement was a reading taken from the face next door to the horizon.

**On the gaze surface** (face-centres are unit spacelike; the surface is $dS_2$):
$q(f_0, f) = s_\hbar + s_G - s_c \in \{\pm1, \pm3\}$. Boost-related to our centre: $C$ ($q=3$).
Directly null-related ($f - f_0$ null): $CH$, $CG$. Null-related to *Hawking's* centre: $H$, $G$.
Boost-related to Hawking's: $GH$. The antipode: $CHG = -f_0$ exactly. No face is
rotation-separated from ours. The lattice is symmetric under $f \to -f$. [PROVED, LABEL-3/3b]

**Witnesses, not inventories.** A face selects a *presentation*. Gravitational redshift is the
operational signature of the $G$-flipped face; it is not a claim that Mass ceased to exist.
[Will]

---

## 5. The sheet above a face

The endpoint records the unordered set of flips. The six temporal orders of three flips fall into
two classes by permutation parity, $S_3/A_3 = \mathbb Z_2$:

$$\mathcal C_+ = \{CHG, HGC, GCH\}, \qquad \mathcal C_- = \{CGH, GHC, HCG\}.$$

Face = unordered physical content. **Sheet = causal orientation of arrival.** [Will; PROVED as
the Pin(2,1) lift, `suites/prim_t5b_parity.py`: the flips are orthogonal reflections, their lifts
anticommute, every adjacent transposition costs $-1$]

The observable that reads the sheet is the raw spinor trace

$$\mathcal V_{\rm spin} = \operatorname{tr}_2(\Gamma_C \Gamma_H \Gamma_G) = 2D, \qquad D^2 = -\det G_c, \qquad \tau: D \mapsto -D,$$

with $D = \det[\hbar, G, c]$ the frame's *oriented* volume. Magnitude = volume, sign = orientation,
defined and zero on the branch locus where only $\operatorname{sgn} D$ is undefined. This is **not**
the projector Bargmann invariant, which is real and even in $D$. [PROVED, `suites/prim_t5c_corrections.py`]

$\mathrm{Cl}(2,1) = M_2(\mathbb R) \oplus M_2(\mathbb R)$; the two blocks give $\mathcal V_{\rm spin} =
\pm 2D$; **the block label is the sheet**; the full $4\times4$ trace is sheet-blind. The
"one-sided reading" is the choice of block. [PROVED, T5c]

**The stellar reading** [CONJECTURE, dynamical]: a neutron star and a black hole are the *same far
face* — thermal, degenerate, radiating. They differ by the sheet: whether quantum support
established first ($H \prec G$) or gravitational collapse did ($G \prec H$). Which sheet is the
trapped surface must be derived from collapse dynamics; $\operatorname{sgn} D$ is constant along
any regular trajectory and $D \neq 0$ at the horizon, so the sheet there is the path/block
holonomy $\chi[\Gamma]$, to be compared with $\theta_+ = 0 \iff 2Gm/c^2R = 1$. [LABEL-3b]

---

## 6. What the labelling does not fix, and the route that would

**ħ versus G.** The geometry gives two positive lines and no way to say which converts what
(T7b3). The octahedron carries it as a contextual ordering. The labelling above *assigns* ħ to
Action/Evanescence and G to Mass/Energy; that assignment is the conjecture.

**The route** [Will; target THM-K]: put a charged load on the seat's compact plane. Its even
response is $c$; its odd response is $Z_0$. The odd response couples in units of $h/e^2$
($e^2 Z_0/\hbar = 4\pi\alpha$, unit-independent), which identifies the ħ-ruler; the remaining
positive ruler is G. The symmetry broken by a **measurement**, not a label. Requires the response
map to be built *before* either ruler is named. Done: the cone (THM-K(a)). Open: the scalar, the
coupling, and the theorem that the seat's two ruler responses are the two electromagnetic
constitutive channels.

**Kill conditions** carried by the labelling:
- Any edge with no known relation of the right kind → that edge's labels are wrong. (LABEL-1: none.)
- The two owed shared-station conditions (ħ, G) fail → the loop labelling is inconsistent.
- The frame path makes the escaping face carry something other than the energy the bound face
  loses → Mass/Energy is wrong.
- $\theta_+ = 0$ and $\chi = \chi_{\rm BH}$ fail to coincide under the collapse map → the stellar
  reading dies.
- The model cannot derive $4\pi\alpha$ from a charged load → HUNCH-Z0 is dimensional analysis in
  a coat, and the ħ/G assignment stays a label.

---

## Sources

`suites/label1_wick_loops.py` (13/13) · `suites/label2_octahedron.py` (11/11) · `suites/label3_lattice.py` (14/14) ·
`suites/label3b_corrections.py` (10/10) · `suites/prim_t5_area_and_lifts.py` (14/14) · `suites/prim_t5b_parity.py`
(19/19) · `suites/prim_t5c_corrections.py` (19/19) · `suites/prim_t7b_labelling.py` (18/18) ·
`suites/hunch_z0_impedance.py` (8/8) · `suites/thm_k_response_map.py` (7/7). The table in §1 is Will's,
2026-09-04.

---

# Part II — The dynamics tier

Built 2026-09-07 under Will's directive: *any part of the model structurally dependent on theory rather than measured
observation is a LOCK, to be removed and re-derived on the model's own steam.* Part I above is the labelling of the
geometry; Part II is the physics the model derives on it. Every runner listed here was executed on the workstation and
committed with its log.

---

## 7. The lapse

**N² = 1 − φ(r)/φ(r_s)**, with φ = −GM/r the Newtonian potential and the normalisation taken where escape velocity = c
(so φ(r_s) = −c²/2). Equivalently A(r) = 1 − r_s/r.

Inputs: the Newtonian potential (ground); conservation of current in d spatial dimensions, giving the shell integral of
$1/r^{d-1}$ (ground); **d = 3 from T8c** (derived, reciprocity ground-adjacent); B = 1/A from THM-I; and Will's
declaration that N² is the ratio of the inner shell to the whole. Einstein's equations are not used.

**The claim is the tail.** First order in φ is automatic for anything matching Newton. What the model asserts is that
the exact lapse is linear in the Newtonian potential **to all orders**. Read off the output and never fed in:

| read off | value | ground |
|---|---|---|
| γ_PPN | 1 | Cassini, γ−1 = (2.1 ± 2.3)×10⁻⁵ |
| β_PPN | 1 | Mercury/MESSENGER, ~10⁻⁴ |
| solar-limb deflection | 1.7512″ | 1.7516″ |
| Mercury precession | 42.9807″/century, by two independent paths agreeing to 2×10⁻⁷ | 42.98 ± 0.04 |
| isotropic g₀₀ at O(U³) | −3/2 | **unmeasured — prediction** |

Controls: in d = 2 the outer shell integral diverges and N² = 0 everywhere, so the dimension is load-bearing and is
T8c's, not the medium's. [DERIVED, `suites/prim_dyn1_lapse_from_potential.py`, 16/16]

**The tilt.** η = r_s/(r − r_s) = sinh²λ, where λ is the rapidity of free fall from rest at infinity relative to a
seat holding station. So η is a presentation quantity: the tilt of the rulers against a **fixed** seat, not a motion of
the seat. N = sech λ is simultaneously the gravitational lapse and the boosted clock rate — the equivalence principle as
an identity in the tilt.

**The σ involution.** Swapping which shell is taken as the numerator sends η → 1/η and N² ↔ v²_esc, exactly. This is an
involution of the radial line with fixed point η = 1 at r = 2r_s. [DERIVED, DYN-1 b1]

---

## 8. Equilibrium

**T·N = T_∞.** Temperature is a rate (ħ and k_B are exact by SI), and a rate divides by the lapse. Derived by two
independent routes — mode-by-mode occupation balance, and Johnson–Nyquist power balance in a band — which agree and are
independent of both the mode and the bandwidth, so equilibrium is well defined. Tolman–Ehrenfest is the output, not an
input.

Consequences: T → ∞ at N → 0 (a station-keeping thermometer has no time to count a rate against); T = √3 T_∞ at the
light ring; and **no static equilibrium exists inside**, since N² < 0 there — a theorem of the potential ratio, not
imported. [DERIVED, `suites/prim_eq1_tolman_from_lapse.py`, 12/12]

---

## 9. The medium

The model's continuum, MAT-2. Three components, each with its own runner.

**The river** — radial flow at the escape velocity, v = −c√(r_s/r). It is a **potential flow** (irrotational,
φ = −2c√(r r_s)) whose convective acceleration (v·∇)v is exactly −GM/r²: Euler with p = 0 and no body force. **Gravity
is the flow's own inertia.** Its Bernoulli integral is ½v² = GM/r, the energy equation of free fall.
[DERIVED, `suites/prim_ns0_river_euler.py`, 6/6]

**The swirl** — the l = 1 harmonic of a rotating boundary. Read as the medium's **vector potential**, not as a viscous
velocity: a dipole's A_φ = m sinθ/r² is the field exactly, ½∇×A is the gyroscope precession field, its energy goes as
1/r⁶, and the river advects it by an induction equation whose φ-component is the same operator as the azimuthal
momentum equation with diffusivity in place of kinematic viscosity. Strength K pinned to one ground number (LAGEOS);
GP-B is then a prediction and lands at 39.1 mas/yr against 37.2 ± 7.2 measured. The node-to-gyroscope ratio is 4,
structural and independent of K. [DERIVED, `suites/prim_dyn3_stokes_swirl.py` 10/10, `suites/prim_edge3_rotating_dissipation.py` 12/12]

*Why vector potential and not viscosity:* a viscous medium torques a rotating boundary with no tide present (spinning
it down on t = M/8πμr_s) and, being linear in its boundary data, admits no tidal torque at all (ε_ijk E_jk = 0 for
symmetric E). Both are contrary to observation. A static field dissipates nothing and does not torque its source, and
every derived result survives the change. Dissipation lives at the boundary, §11.

**Sound** — the medium's only propagating field. Linearised Euler plus continuity on the river give an acoustic metric:
sound moves at ±c_s **relative to the flow**, so the cone is carried by the medium. With c_s = c, a shift of the time
coordinate takes that acoustic metric to **exactly the lapse metric of §7**: g_TT = −c²(1 − r_s/r), g_rr = 1/(1 − r_s/r).
[DERIVED, NS-0 b1, b2]

**This is the AC sector's construction principle, and it has a consequence the model must own: the medium carries a
spin-0 field and no spin-2 field.** Everything read at the boundary in §11 is sound. The characteristic ratio √(L/C) of
this sector is ρc_s; its *value* is open (THM-K's scalar half, §6).

**O(J²).** The measured 1/r³ frame dragging bounds the medium's advection length (GP-B against LAGEOS: ≲ 1.3×10⁴ km),
which puts the transverse dynamics in the creeping regime. There the flow is linear, river and swirl superpose, and the
O(J²) exterior potential −|v|²/2 has angular structure at 1/r⁴, **not** 1/r³. The coefficient of P₂/r³ is exactly zero:
**the spin-induced quadrupole κ = 0** (Kerr: 1). The trapping surface stays spherical; the ergosurface is oblate.
[DERIVED, `suites/prim_ns1_swirl_inertia.py`, 10/10]

---

## 10. Orbits and the fold

All three circular-orbit radii come out of the §7 metric, none put in: light ring 3r_s/2, ISCO 3r_s, marginally bound
2r_s.

**The σ involution maps the orbit families.** As an identity in r,

$$(L/E)_{\rm circ}(r) \;=\; b_{\rm turn}(\sigma r)/\sqrt2 ,$$

so the ISCO condition is the image of the light-ring condition under σ, by the chain rule — not a coincidence of two
roots. σ is **not** a symmetry of the geodesic flow (the E² term picks up (ρ−r_s)⁴/r_s⁴), so this is an orbit-family
symmetry of a spacetime that has no such metric symmetry. Lineage: weaker and different from the Couch–Torrence
inversion of extremal Reissner–Nordström. Novelty unverified against the literature.
[DERIVED, `suites/prim_gud1_complement.py` 8/8, `suites/prim_gud2_inversion_map.py` 11/11]

In the tilt angle θ with η = tan²θ, the two families are one function with sine and cosine exchanged:
b_turn = r_s/(sin²θ cos θ), b_circ = r_s/(√2 sin θ cos²θ). The light ring sits at P₂(cos θ) = 0 (θ = 54.74°), the ISCO
at P₂(sin θ) = 0 (35.26°), the marginally bound orbit at θ = 45° — the eighth-turn W of §3, at a radius. And
N² + v² = 1: the lapse and the river are the two legs of one right triangle whose angle is the tilt.

**The fold.** The light ring is a fold of the ray family, with λ_L = Ω_c — the instability rate equals the winding rate
(two independent paths, agreeing to 2.6×10⁻⁶). Consequences: the photon-ring subring ladder is **e^{−π} per half-orbit**
(4.32 %, unmeasured — next-generation space VLBI); eikonal quasinormal modes ω = Ω_c(l+½) − iλ_L(n+½), whose real part
is the winding (Light/Action) and imaginary part the decay (Evanescence) — face H's flip; and the greybody edge is a
Fermi function in frequency of width λ_L, which is **not** κ (ratio 4/(3√3)). M87*'s shadow from the fold: 39.7 μas
against a 42 ± 3 μas ring. [DERIVED, `suites/prim_fold1_lightring_dressing.py`, 10/10]

---

## 11. The boundary at r = r_s

Where the river runs at c the sound cone stalls: no outward-directed characteristic. Three equivalent statements of one
boundary — N → 0 there, no reflected wave, no outward characteristic.

**Absorption.** σ_abs(ω→0)/A_H = 0.9993, with |T₀|² = 4(ωr_s)². The long-wavelength absorption area equals the
capture area of the marginally bound orbit, πb²_mb = 4πr_s². Flux conserved to 10⁻⁸ at every frequency.
[DERIVED, `suites/prim_edge1_horizon_absorption.py`, 7/7]

**Tidal response.** A tide enters the medium through Bernoulli as an external potential: −v²/2 = −M/r + Φ_tide exactly,
so the induced exterior multipole is identically zero at every order while the sonic surface deforms. **k₂ = 0.** GR's
tensor Love number, computed as a labelled comparison, is zero by the same mechanism — a polynomial regular branch and a
log-singular tail-carrying branch, with regularity at the boundary killing the tail.
[DERIVED, `suites/prim_edge2_tensor_love.py` 7/7, EDGE-3A]

**Transport coefficients, derived.** The §7 metric is **Ricci-flat** (computed, not assumed), so the stalled sound rays
obey the vacuum Raychaudhuri equation with no matter term. With κ = 1/(2r_s) from the metric, T = κ/2π (the one declared
line, KMS), the first law with r_s = 2M giving **S = A/4 as output** (Bekenstein not used), and the teleological steady
state θ = σ²/κ, the dissipated power is (1/8π)∫σ²dA. Hence

$$\eta_s = \tfrac{1}{16\pi}, \qquad \zeta = -\tfrac{1}{16\pi}, \qquad p = \kappa/8\pi .$$

Damour's three membrane coefficients as output, with no Einstein equations. The 16π assembles as 2 × (2π from the
thermal period) × (4 from S = A/4). ζ < 0 — the boundary is an active element. [DERIVED-given-KMS,
`suites/prim_edge6_membrane_viscosity.py`, 10/10]

**Tidal friction.** The rotating boundary as a membrane: the tidal deformation of the sonic surface rotates under it,
area preservation forces a tangential flow, the flow shears, the shear dissipates, J̇ = −P/Ω. The torque's tensor
structure is Σ_m m²|A_m|² = (8π/15)[E₁ − (3/2)E₂] — **exactly** the invariant combination of Poisson's small-spin
result, from the angular integral alone. Zero for a tide axisymmetric about the spin; zero at Ω = 0; the first law
closes. [DERIVED, `suites/prim_edge4_tidal_friction.py`, 8/8]

**Magnitude.** Reproducing GR's torque with this membrane would need η_s/32; the derived η_s is 32× larger. The
discrepancy is the field, not the mechanism: a spin-2 l = 2 field is absorbed **36×** more strongly than a spin-0 one by
the same boundary — Page's factorial, recovered here by computing both coefficients with one integrator (spin-0: 1/2025;
spin-2 via Regge–Wheeler as comparison: 4/225, the classic 256/225 (Mω)⁶; ratio 36.006). The model's tide is spin-0
(§9), so its tidal friction is ~1/32–1/36 of Kerr's. [DERIVED, `suites/prim_edge5_spin_ratio.py`, 7/7]

**The thermal scale.** κ = 1/(2r_s) is the metric's slope at the boundary; near it N = κρ − ρ³/16r_s³, the Rindler form.
Wick the seat and the line element is polar about the boundary; single-valuedness of a potential around the loop fixes
the period at β = 4πr_s/c, and one circuit rotates the frame by 2π — +1 on vectors, **−1 on spinors**, forced by
smoothness, which is §3's W⁸/W¹⁶. Given KMS, T_ref = ħc/(4πk_B r_s). The Unruh form is the same declaration read at the
seat rather than at the boundary, so it is not a second assumption. [DERIVED-given-KMS,
`suites/prim_dyn2_horizon_period.py`, 12/12]

**The Euclidean section.** A < 0 inside, so the section covers r ≥ r_s only: **r = 0 is absent from it**. The tip is a
regular interior point where the thermal circle shrinks. The section is a cigar (disc × S², simply connected, radial
direction non-compact) — not a torus. The Lorentzian interior and r = 0 still exist; "no centre" is true of the
Euclidean section only.

---

## 12. What the model says that Kerr does not

| quantity | model | Kerr | ECOs | instrument |
|---|---|---|---|---|
| spin-induced quadrupole κ | **0** (δκ_s = −1) | 1 | ≠1 | GWTC-4.0 gives [−53, +9] hierarchical; one unit apart on a 26–62 unit ruler — **cannot yet distinguish**. Needs σ(δκ_s) ≲ 0.3 |
| tidal Love number k₂ | 0 | 0 | ≠0 | no BH Love number ever measured |
| σ_abs(0)/A_H | 1 | 1 | ≠1 | — |
| tidal friction | **~1/32–1/36 of Kerr** | 1 | varies | LISA EMRI horizon-flux phasing; LVK horizon absorption |
| swirl profile | ωr³ = const at all orders | O(J³) corrections | — | ringdown no-hair |
| photon-ring subrings | e^{−π} per half-orbit | same | — | next-gen space VLBI |
| isotropic g₀₀ at O(U³) | −3/2 | same | — | unmeasured |

Two entries BH-like, two not, none ECO-like. The combination κ = 0 with black-hole tidal properties is not produced by
anything in the ECO literature.

**Conditional on:** the medium's tide remaining spin-0 (§9). If the AC sector's construction gives the medium a spin-2
channel, the tidal-friction entry moves.

---

## 13. Residue

**One declared line:** KMS — the period of the seat's imaginary-time loop is ħ/k_BT. Everything thermal inherits it.

**Two quantities pinned to ground:** d = 3 (T8c, reciprocity ground-adjacent) and K = 2GJ/c² (LAGEOS, 1 %).

Everything else in Part II is derived from the Newtonian potential, Laplace, Euler, induction and the primitives.

**Open constructions**, in order:

1. **The AC sector — L and C.** One construction carrying: the value of the characteristic ratio (THM-K's scalar half,
   ζ = e²Z₀/ħ = 4πα); whether the medium has a spin-2 tidal channel (§12's condition); and kill 1 of
   `docs/CONJECTURE-COSMOLOGY.md`.
2. ~~SHEET-1's c-block.~~ **Rebuilt** as `suites/prim_sheet1c_two_charts.py` (11/11): the infalling observer's
   reference frame is the **untilted** presentation, η = 0 at every r, so its D is the frame's D and is constant — the
   old D → ∞ came from RIDE-2's strain and is gone with it. The station-keeping presentation is the tilted one,
   η = sinh²λ = r_s/(r−r_s), and its family **terminates at r_s** (holding position below needs v > c relative to the
   river), so its interior readings are a formula continued past its last occupant: η < −1 and N² < 0 is an imaginary
   clock rate for a seat that is not there. det G′ vanishes at r = r_s(1 − kᵀG⁻¹k) — a locus whose position is set by
   the seat's own direction, which is what an artefact does and a place does not. The block label references no η, so
   the sheet reads the same in both charts. No new declaration; RIDE-2 not used.
   *Still open:* what the infalling chart reads **at** r = 0 is a question about the frame, not about either
   presentation, and the runner does not answer it — only that neither chart diverges there.
3. **Unpin K:** the medium's own angular-momentum content.
4. Finite-Reynolds meridional flow at O(J²); the trapping surface's shape as a derived boundary condition (an oblate
   one is the only way the medium recovers Kerr's Q, and the model does not produce one).
5. Structural, unrun: the ħ- and G-seat readings of the ground; the lit face on the octahedron's vertices; σ as an
   object on the Gram rather than on the pinning.

---

## Sources — Part II

`suites/prim_dyn1_lapse_from_potential.py` (16/16) · `prim_eq1_tolman_from_lapse.py` (12/12) ·
`prim_dyn2_horizon_period.py` (12/12) · `prim_dyn3_stokes_swirl.py` (10/10) · `prim_sheet1_block_and_holonomy.py`
(17/17, c-block withdrawn) · `prim_fold1_lightring_dressing.py` (10/10) · `prim_gud1_complement.py` (8/8) ·
`prim_gud2_inversion_map.py` (11/11) · `prim_ns0_river_euler.py` (6/6) · `prim_ns1_swirl_inertia.py` (10/10) ·
`prim_edge1_horizon_absorption.py` (7/7) · `prim_edge2_tensor_love.py` (7/7) · `prim_edge3_rotating_dissipation.py`
(12/12) · `prim_edge4_tidal_friction.py` (8/8) · `prim_edge5_spin_ratio.py` (7/7) ·
`prim_edge6_membrane_viscosity.py` (10/10). External imports: `docs/EXTERNAL-MATHEMATICS-DEBT.md`, EXT-027…039.

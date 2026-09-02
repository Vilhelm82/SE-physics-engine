# THM-TARGET I — light, matter and the field law in the pinned family
Written 2026-09-03.  Three suites, all green on the 7800X3D, all committed with logs:
`thm_i_pre.py` **19/19** (`baf7401`), `thm_i_transport.py` **11/11** (`67f573f`), `thm_i_field.py` **10/10** (`3fedde7`).
Read after `docs/prereg/CURV-1/CURV1-RESULTS.md`.  Re-running is the verification.

## The question and what it became
Opened 2026-09-02 evening as "light bending from the presented geometry alone: factor 1 (Newton) or 2 (Einstein)?".
The pre-tests replaced the planned seat-to-seat integration with an exact identity, the transport probe found the
one local rule that reproduces it, and the field probe turned "where does r_s come from" into "which variable is
harmonic", with a measured number selecting the answer.  Every result below is a **structural retrodiction**
conditional on the declared inputs named in each section.  Nothing here derives r_s from matter; the last section
says exactly what is left.

## Inputs
Carried and re-derived in-suite: BARE-1 rotor and sandwich; RULE-1 pairing; E-4 presented scalar; KIN-2a pinning
tanh λ = √(r_s/r); CONT-1; the c-seat null locus; H-7 det G′ = Δ cosh²λ.
New, named this target:
- **STAT-1** — the pinned rotor field does not depend on time, and a ray's frequency with respect to the frame's
  common time office (its pairing with the static seat's presented time unit A·1·A, divided by that unit's scalar
  part) is constant along the ray.  [declared principle]
- **ROT-1** — the pinned family is invariant under rotations about the centre, and the pairing with the rotation
  generator (0, ẑ × r⃗) is constant along the ray.  [declared principle]
- **PROP-1** — position bookkeeping: the ray advances at c along its direction in the local frame while the frame
  drifts at the pinned velocity toward the centre relative to the static seats.  [construction]
- **RULE G** — the local transport law found in part A (below).  [derived from STAT-1/ROT-1: the unique candidate]

## 1. Pre-tests (`thm_i_pre.py`)
- **Redshift.**  The static seat reads a ray by pairing with A·1·A; the frozen E-4 convention gives e^{−λ} for an
  outgoing ray (redshift; sign not chosen).  With STAT-1, ω_static ∝ cosh λ, so between static seats
  ω₂/ω₁ = √((1−r_s/r₁)/(1−r_s/r₂)) = 1 − (r_s/2)(1/r₁ − 1/r₂) + …: coefficient ½.
- **Two local rules die in 1D.**  Boosting the ray between neighbouring rain frames by their velocity difference is
  wrong by the factor (1+v) — office ghost #5.  D-14's static relative rotor is first order in √r_s where the redshift
  is second order: it relates seats, it does not propagate rays.
- **Horizon.**  c(1 − tanh λ) = 0 exactly where λ → ∞, where ρ_K and det G′ have their pole, where e^{−λ} → 0.
- **The orbit identity.**  STAT-1 + ROT-1 + the null locus give, for ANY radial profile,
  (dr/dφ)² = r⁴/b² − r²(1 − tanh²λ), b = L/ω_t.  KIN-2a turns 1 − tanh²λ into 1 − r_s/r.
- **The deflection.**  b·δ → 2.000003 at b/r_s = 10⁶; second-order coefficient 2.945248 against 15π/16 = 2.945243.
  The factor 2 has an address: the rod office supplies 1 − r_s/r, the time office supplies the drift in ω_t.

## 2. The local transport rule (`thm_i_transport.py`, part A)
Three candidate local maps X(x+dx) = Λ X Λ†, symbolic first-order drift of (ω_t, L), then integrated as ODEs:

| rule | what it is | d ω_t/dt, dL/dt at a generic point | swept angle vs exact (b = 100 r_s, r₀ = 1000 r_s) |
|---|---|---|---|
| G | boost by the gradient of the pinned velocity field along the ray's own co-moving step, dλ⃗ = (n̂·∇)v⃗ dt, passive sign | **0, 0 identically**, every profile, every direction (A-1) | agrees to 7.8×10⁻⁷, pairings conserved to 10⁻¹⁵ |
| N | naive product of seat rotors A(x+dx)A(x)⁻¹ (the D1 composition, with Wigner rotation) | 0.0417, 0.392 | off by 0.0345 — more than the whole deflection (0.02) |
| F | D-14 frame-transported relative rotor R(dθ)A(dλ,K) | 0.0759, −0.0344 | off by 0.75; pairings drift to 10³⁶ |

The only property of the field rule G uses is that ∇v⃗ is symmetric (the pinned flow is irrotational, A-0).
Propagation is **local**; the local law is the derivative of the pivot field, not a comparison of two seats against
a background.  That is the stopwatch objection with a positive half: the rotors relate seats to the frame;
frame-to-frame transport needs the field's gradient, which carries no background.
(Log cosmetic: the A-6 "deflection" numbers are swept angle minus π at finite r₀, not asymptotic deflections.)

**Massive bodies.**  With a timelike paravector the pairings are E = m(cosh ζ − sinh ζ tanh λ cos α) and
L = m r sinh ζ sin α, and the orbit closes as the identity
(dr/dφ)² = (r⁴/L̃²)[Ẽ² − (1 − tanh²λ)(1 + L̃²/r²)] for any profile (M-2); rule G conserves E and L with the body's
own step e⃗ = u⃗ dt (M-3); the perihelion advance for r_p = 2×10⁴, r_a = 3×10⁴ r_s gives Δφ·ℓ/π = **3.0003** (M-4).

## 3. The field law (`thm_i_field.py`, part B)
**Forced for every profile (kinematic tier):** ρ_K·N² = cosh²λ·sech²λ = 1.  The presented radial rod and the static
clock are one Lorentz factor, so the spatial stretch per unit potential is 1 — PPN γ = 1 with no field law at all.

**The fork.**  "H is harmonic outside the source, H = H_∞ + a/r, with a fixed by the weak field β² → r_s/r", for the
natural candidates.  All agree at first order (bending 2, redshift ½, horizon).  They differ at second order:

| harmonic variable | c₂ = (r_s/r)² term of tanh²λ | 2nd-order light coeff | perihelion coeff |
|---|---|---|---|
| **sech²λ = Δ/det G′  (KIN-2a)** | 0 | **2.9452431** (15π/16) | **3.0003** |
| cosh²λ = ρ_K (presented volume) | −1 | 0.5890 | 1.0000 |
| ln cosh²λ | −½ | 1.7671 | 2.0001 |
| sech λ (the lapse) | −¼ | 2.3562 | 2.5002 |
| cosh λ | −¾ | 1.1781 | 1.5000 |
| textbook | — | 2.9452431 | 3 |

Two laws found in the table and then checked: perihelion coefficient = 3 + 2c₂ (B-2e), second-order light
coefficient = 15π/16 + (3π/4)c₂ (B-2f).  Mercury's perihelion (relative precision 10⁻⁴) selects c₂ = 0.  So the
harmonic variable is sech²λ = Δ/det G′: **the reciprocal presented volume, in units of the frame's own volume**.
In the model's words the field law the data demand is

    ∇²(Δ/det G′) = 0 outside the source,   r² d(Δ/det G′)/dr = r_s,

the Euler–Lagrange equation of the gradient energy ∫|∇(Δ/det G′)|² d³x.  The rival "presented volume harmonic"
(ρ_K = 1 + r_s/r) is the second row and is dead at Mercury.  The other three second-order light coefficients are
predictions of theories that no longer exist.

## Predictions scorecard for THM-I (from the 2026-09-02 evening exchange)
| | | |
|---|---|---|
| factor 2 not 1 | **yes** | analytically (PT-3c/d), numerically (PT-3e) |
| the presentation carries Einstein's extra half | **yes** | rod office 1 − r_s/r, time office the drift; each half |
| seat-to-seat transport imports a convention | **yes, and worse** | both available conventions give wrong numbers (PT-1i/j, A-2/3/6) |
| a local rule exists at all | **yes** | the gradient rule; unique among the candidates (A-1, A-5) |
| γ = 1 structural | **yes** | ρ_K N² = 1 (B-1) |

## What is now decided, and what is not
Decided (conditional on KIN-2a, STAT-1, ROT-1, PROP-1):
1. The ray and orbit kinematics of the pinned family are the full Schwarzschild geodesic structure, null and timelike,
   second order included — the model has no freedom left once the rapidity field is given, because a seat's rods and
   clocks are one boost.
2. The local transport law is the gradient of the pivot field along the transported object's own co-moving step.
3. The field law is narrowed to one variable: the reciprocal presented volume Δ/det G′ is harmonic outside the source.

Not decided:
- **Where r_s comes from.**  The flux r² dH/dr = r_s is the source strength; identifying it with the blind-mass measure
  (pivot-blind moment = total mass of the source measure) is the declaration KIN-2a′ would make.  Not derived.
- **Whether the dynamics tier produces the law.**  The nine generators, the stretch–angle couplings and the dissipator
  have not been asked whether ∫|∇(Δ/det G′)|² is their static energy.  That is the next target, and it is now a
  question with one variable and one exponent in it rather than a wish.

## Comparison stage (names spoken here only)
STAT-1 and ROT-1 are the Killing conservation laws of a stationary, spherically symmetric field; PROP-1 is the
Painlevé–Gullstrand / river picture; rule G is the spin connection of the rain tetrad in PG form (boost part = velocity
gradient on the co-moving displacement, rotation part zero for irrotational flow); rule N's error is boosting
neighbouring river elements against a background; rule F is the static-observer relation.  PT-3d and M-2 are
Schwarzschild's null and timelike geodesic equations; 2r_s/b is Einstein 1915; 15π/16 (r_s/b)² is the second-order
deflection (Epstein–Shapiro 1980); the third-order residual seen in the b = 10⁴ run, c₃ ≈ 5.4, is consistent with the
textbook 16/3.  ρ_K N² = 1 is PPN γ = 1 (Cassini); the perihelion coefficient is (2 + 2γ − β)/3 × 3 = 4 − β, so 3 is
β = 1 (Mercury, LLR); "sech²λ harmonic" is "the lapse squared is harmonic", Schwarzschild's g_tt exactly, and the
gradient energy is the Newtonian field energy of 2Φ/c².  The weak-field ½ is Pound–Rebka.

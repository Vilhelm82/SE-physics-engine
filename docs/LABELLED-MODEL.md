# The labelled model — the whole octahedron

**Status: CONJECTURE tier throughout, except where marked.** The geometry it sits on is proved
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
faces, and they differ only in the axis the edge does not touch. [PROVED, `label2_octahedron.py`]

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
[PROVED, `label2_octahedron.py`]

**Shared-station consistency.** Loops A and B both pass through Light and Temperature: the
temperature reached via ħ (Unruh) must equal the one reached via G (Hawking). **DERIVED — T4.**
Loops A and C share Action/Evanescence; B and C share Mass/Energy: those two conditions are
**OWED.** [LABEL-1]

---

## 4. The eight faces

Each face is a triple — one pole per axis — and its physics is the sum of its three borders.
The faces form the Boolean lattice $2^{\{C,H,G\}}$ of three flips: $C$: Light → Temperature,
$H$: Action → Evanescence, $G$: Mass → Energy. Ours is $\varnothing$; Hawking's is $CHG$, its exact
antipode $-f_0$. [PROVED, `label3_lattice.py`, `label3b_corrections.py`]

| flips | poles | contributions | physics (DECLARED candidate) | status |
|---|---|---|---|---|
| $\varnothing$ | Light, Action, Mass | GR + Wave + Bound | light waves in bound orbits around mass, in curved spacetime — **our face** | established |
| $G$ | Light, Action, Energy | GR + Wave + Escaping | light with phase leaving a well: gravitational redshift, Shapiro delay, escaping GW | named 09-04 |
| $H$ | Light, Evanescence, Mass | GR + Particle + Bound | light tunnelling in a bound system: near-horizon greybody filter (G-2) | named 09-04 |
| $C$ | Temperature, Action, Mass | Quantum + Wave + Bound | thermal waves around bound mass: Tolman–Ehrenfest (EQ-1) | named 09-04 |
| $HG$ | Light, Evanescence, Energy | GR + Particle + Escaping | light tunnelling out carrying energy: Parikh–Wilczek emission, pre-thermal | named 09-04 |
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
the Pin(2,1) lift, `prim_t5b_parity.py`: the flips are orthogonal reflections, their lifts
anticommute, every adjacent transposition costs $-1$]

The observable that reads the sheet is the raw spinor trace

$$\mathcal V_{\rm spin} = \operatorname{tr}_2(\Gamma_C \Gamma_H \Gamma_G) = 2D, \qquad D^2 = -\det G_c, \qquad \tau: D \mapsto -D,$$

with $D = \det[\hbar, G, c]$ the frame's *oriented* volume. Magnitude = volume, sign = orientation,
defined and zero on the branch locus where only $\operatorname{sgn} D$ is undefined. This is **not**
the projector Bargmann invariant, which is real and even in $D$. [PROVED, `prim_t5c_corrections.py`]

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

`label1_wick_loops.py` (13/13) · `label2_octahedron.py` (11/11) · `label3_lattice.py` (14/14) ·
`label3b_corrections.py` (10/10) · `prim_t5_area_and_lifts.py` (14/14) · `prim_t5b_parity.py`
(19/19) · `prim_t5c_corrections.py` (19/19) · `prim_t7b_labelling.py` (18/18) ·
`hunch_z0_impedance.py` (8/8) · `thm_k_response_map.py` (7/7). The table in §1 is Will's,
2026-09-04.

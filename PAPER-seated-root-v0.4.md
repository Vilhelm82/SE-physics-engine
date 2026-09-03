# The Seated Root
## How much physics follows from three axes through a point: a derivational-economy study with machine-verified receipts

**W. Lloyd** — draft v0.4, 2026-09-04
*v0.4 (2026-09-04, later): PRED-1's Presentation paragraph carries the transition closed form T_ij = R_i R_j~, the pole-triangle dressing with SIGNED holonomy, the corrected parity bookkeeping (the pole triangle does not flip under the deck; the dressing is mixed through h alone), and the exact statement that V is the only primitive deck-odd generator — the seat supplies a deck-even, sign-varying gain and no independent deck bit. v0.3 retained in git history.*
*v0.3 (2026-09-04): PREDICTION-1 rewritten around presentation-versus-operation — the pivoted observation map derived cold (bilinear pairings: deck-even scalar office, mixed transition readings with closed forms, the pole-triangle dressing), the ordered-register operation reading sgn V with three exact physical readouts, the algebraic half certified three ways, and the cycle protocols recorded as an existence witness. Owed in this version and not yet folded: THM-O (the strain law as Gauss and Codazzi on flat rods; fork (b) at O(J^2)) into §6, VIEW-1 and CENSUS-C into §2/§9, and the scaffold register into §10 — see HANDOFF-2026-09-04.md. v0.2 (2026-09-03) is retained in git history.*
*v0.2 (2026-09-03): the pinned family worked through — the hole as a presented state (THM-H), what sits at the
centre (CURV-1), light and matter in the family (THM-I), and the field law narrowed and shown not to follow from
the bare tier (THM-J). §6 gains four blocks; §10's ledger is rescoped and records the programme's first falsified
model-internal candidate. Amended the same evening: THM-K derives the pinning from the seat's temperature, THM-L carries rotation, THM-M derives
the swirl as the boost of the pinning; §6's last block flipped from "not derived" to "derived from thermal inputs".
v0.1 (2026-08-28) is retained in git history.*
*All results in this paper are backed by exact symbolic computation committed to the project repository; commit hashes appear in §11.*

---

### Abstract

We study a deliberately minimal object: three polarized axes through a common
root, with state given by the Gram matrix $G$ of the axes, observers modelled
as *seats* (an observer occupies one axis, which becomes its unmeasurable
reference), and changes of perspective modelled as hyperbolic rotors
("pivots") in the Clifford algebra $\mathrm{Cl}(3)$. Under a strict
banned-input discipline — each derivation declares what it may not use, and
predictions are frozen with hashes before computation — this structure
reproduces a surprising amount of standard physics: the classical relativistic
temperature-transformation family and its resolution (with explicit ancestry
in Biró–Ván), the uncertainty floor of conjugate quadratures, the
horizon-vanishing structure of static gravitational observers, the Unruh and
Hawking temperatures *with their coefficients selected by an internal theorem
where a wrong coefficient was available*, an emergent-time architecture
consonant with Wheeler–DeWitt and Page–Wootters, and a derived (not posited)
seat-level dissipator obeying the fluctuation–dissipation theorem. Pinning the
rotor's rapidity by position (KIN-2a, declared) presents the frame's rods as the
static geometry of a black hole: the horizon and centre are degenerations of
the *presented* state while the frame's own state never moves; the centre is a
finite event on every rod-free curvature reading (an angle for the infalling
seat, a spin sign for the hole); the family's ray and orbit kinematics is
Schwarzschild's exactly, second order included, because a seat's rods and
clocks are one boost; and the field law behind the pinning is narrowed to one
variable — the reciprocal presented volume must be harmonic — shown *not* to
follow from the bare tier's own geometry (the first model-internal candidate
this programme proposed and had killed by a measured number), and then derived
from its *temperature*: the seat's presented sphere, given an entropy
proportional to its area and equilibrated with the seat's own Unruh
temperature as seen from infinity, returns the pinning exactly, with the
Schwarzschild radius falling out. Boosting that field gives frame dragging with
the measured coefficient, half from the seat's clock and half from its rods,
so the model's gravitational sector now rests on one declared idea — the
thermal content of a seat's sphere — and no imported profile. The paper's
claim is **not** new physics: nearly every physical output is a retrodiction
with named lineage. The claim is *economy*: the input set is small, the
banned lists are enforced, several freeze-dated internal theorems later
selected correct coefficients unprompted, and the framework yields a small
number of new mathematical identities plus one prediction-shaped statement —
the medium of the model is readable as a Bargmann holonomy around the seat
cycle. Every theorem is machine-verified in exact arithmetic; every numeric
claim is corroborated to stated precision; every conjecture is labelled.

---

### 1. Introduction, and what kind of paper this is

The project's founding motivation is the representation of singular behaviour:
to build a state space in which "the singularity" is a locus with computable
structure rather than a place where one writes $\infty$. The resulting
programme has an unusual epistemology which we state up front, because the
paper is unreadable without it.

Every result below carries a tier label: **[proved]** (exact symbolic theorem,
machine-checked), **[derived | X]** (exact, conditional on a declared input X),
**[retrodiction]** (recovers known physics; lineage cited), **[conjecture]**,
or **[declared]** (an input, frozen and named, never silently assumed). Each
computation preregisters its predictions and kill conditions before running;
several documents in the corpus carry SHA-256 hashes frozen prior to
execution. Where possible, results are verified by two independent
computational routes sharing nothing but the computer-algebra system. The
outcome of any derivation is booked into one of three classes: (a) it
reproduces known physics (a retrodiction receipt), (b) it *differs* from known
physics in a testable way (a prediction), or (c) it cannot close without a new
declared rule (a named debt). The honest headline of this paper is that
class (a) dominates overwhelmingly; the interesting exceptions are catalogued
in §10.

### 2. The bare tier: state space, symmetry, and the branch surface

**The object.** Three unit axes $a_1, a_2, a_3 \in S^2$ through a common
origin (the *root*). The intrinsic state is the Gram matrix
$G_{ij} = a_i \cdot a_j$, with unit diagonal and off-diagonal angle
coordinates $\gamma_{12}, \gamma_{13}, \gamma_{23}$. The capacity is
$\Delta := \det G = 1 - \gamma_{12}^2 - \gamma_{13}^2 - \gamma_{23}^2
+ 2\gamma_{12}\gamma_{13}\gamma_{23} \geq 0$, vanishing exactly at rank loss.
A *seat* is the occupation of one axis by an observer; the seated axis is the
observer's reference and does not appear in its own measurements (proved at
the two-axis tier as ROOT-1: the seat projection annihilates the seated axis).
Relabelings form the signed permutation group $B_3$ (order 48); the *click
group* is $\Gamma = B_3/\{\pm I\} \cong S_4$ **[proved]**.

**The branch surface is the Cayley cubic, and this is forced.** Because the
state is a symmetric matrix of affine-linear forms, $\Delta$ is a cubic
*symmetroid*; the singular points of a symmetroid are its rank-one points;
unit diagonal plus rank one forces $G = vv^{T}$ with $v \in \{\pm 1\}^3$,
i.e. exactly four nodes — the maximum for an irreducible cubic surface, whose
unique representative is Cayley's nodal cubic **[proved]**. The four nodes are
the four total-collision states and simultaneously the four body diagonals of
the cube; the classical fact $\mathrm{Aut}(\text{Cayley}) = S_4$ then shows
the click group is not a modelling choice but the automorphism group of the
boundary geometry **[proved + classical]**. Of the surface's nine lines, six
are the pairwise collision loci $\gamma_{ij} = \pm 1$ (the $K_4$ edges on the
nodes); three lie at infinity, outside the model. The physical region
$\{\Delta \ge 0, |\gamma_{ij}| \le 1\}$ is the $3\times 3$ **elliptope** of
correlation matrices: its smooth rank-two boundary is the coplanar stratum,
its four rank-one vertices are the collisions. Convex-geometric and
semidefinite-programming machinery therefore applies to the state space
wholesale [lineage: classical algebraic geometry; SDP literature].

![Figure 1](fig1_elliptope.png)

**Orientation, cover, and monodromy.** The signed volume $V = \det[a_1 a_2
a_3]$ satisfies $V^2 = \Delta$: oriented frames form a double cover of the
Gram world, branched over $\Delta = 0$. Because $\Delta$ is absolutely
irreducible (its singular locus is zero-dimensional), the monodromy group of
the cover over the complexified state space is exactly $\mathbb{Z}_2$
**[proved]**: loops about smooth branch points (coplanar *and*
collision-type alike) send $V \to -V$; loops about the four nodes, where
$\Delta$ vanishes to order two, do **not** flip $V$ **[proved + verified to
30 digits]**. This dichotomy is an instance of a general mechanism the author had already
proved, six weeks earlier and for a different family: the *inertia
localization theorem* for axial $k$-ellipse norms (2026-07-10) establishes
that all monodromy is generated at the critical locus while sheet-collision
components are nodes with trivial inertia. The shared mechanism is the parity
of the local vanishing order — odd order flips the sheet, even order does not
— and the quantitative echo is exact at $k=3$: that note's genus law gives
$g=1$, and a generic line section of the seated-root cover is $w^2 = $ (a
squarefree cubic), an elliptic curve **[verified]**. The two families are
different objects (a curve cover over an $M$-line; a surface cover over the
Gram cell), so the mechanism is shared but no identification is claimed.

Two distinct facts coexist here and must not be conflated: the
relabeling extension $1 \to \{\pm 1\} \to B_3 \to \Gamma \to 1$ **splits**
(so $V$'s sign is removable as a group-theoretic convention — $V$ is
tensorial), while the orientation cover is topologically **nontrivial** (no
continuous global branch of $V$ exists). A split extension says nothing about
holonomy in a punctured space. The same involution is realised by two
mechanisms of different type, one removable and one not **[proved]**.

**The covers, algebraically: a conjugate Kummer module.** The topological
statement above has an exact field-theoretic counterpart. Over the state field
$F = \mathbb{Q}(\gamma_{12},\gamma_{13},\gamma_{23})$ the constructions of this
paper introduce four radicands: the click-invariant $\Delta$ (the orientation
cover) and the three $1+\gamma_{ij}$ required by the half-angle spinor lift of
§9, which the click group permutes. All four are irreducible and pairwise
coprime, so by the odd-valuation criterion their square classes are independent
in $F^{*}/F^{*2}$ — rank four, multiquadratic group $(C_2)^4$ **[proved]**.
Because the click group permutes the three conjugate classes while fixing
$[\Delta]$, the covers assemble into a *conjugate Kummer module* with normal
closure
$$\mathrm{Gal}(L/F) \;\cong\; (C_2 \wr S_3) \times C_2, \qquad |{\cdot}| = 96,$$
**[derived]**. One polynomial relation ties the invariant radicand to the
conjugate triple:
$$(1+\gamma_{12}+\gamma_{13}+\gamma_{23})^2 + \Delta \;=\;
2(1+\gamma_{12})(1+\gamma_{13})(1+\gamma_{23}) \quad\textbf{[proved]},$$
and it is the mechanism behind the descent theorem of §9. Note that this
$\mathbb{Z}_2$ — a Kummer square class — is *not* the $\mathbb{Z}_2$ of the
non-split spin extension: same order, different category, and the two must not
be conflated.

**The two-axis ground floor.** A parallel, fully preregistered campaign
(AGNOSTIC-1, tasks T1–T6; SHA-frozen, byte-stable, independently re-run) had
already established at $(k,d)=(2,2)$: the invariant-ring lattice, the
orientation cover $\omega^2 = \Delta$, the conserved circulation momentum for
*any* rotation-invariant potential, and the selection theorems: $(2,2)$ is
the unique planar spanning rank with a single angle and *protected* chirality;
the three-axis model of this paper is the $(3,3)$ rung of the same ladder,
where protection is restored with three angles **[proved]**.

### 3. Pivots, composition, and the Kähler structure

**Pivots.** A change of perspective is a hyperbolic rotor
$B(\lambda) = \exp(\lambda K/2)$ with $K^2 = +1$ in $\mathrm{Cl}(3)$
**[declared: BARE-1; collinear leg discharged separately]**. Complexifying
the rotor angle identifies boosts as imaginary-angle rotations of their
planes: $-I\,(e_2 e_3) = e_1$ exactly, where $I = e_1e_2e_3$ is the central
pseudoscalar **[proved]**. Non-collinear composition was the model's own kill
condition ("the non-commutativity must produce the Wigner rotation, or the
labelling fails"): two boosts of rapidities $\lambda_1, \lambda_2$ along axes
at angle $\alpha$ compose into a boost times a rotation of angle
$$\tan(\omega/2) \;=\; \frac{s_1 s_2 \sin\alpha}{c_1 c_2 + s_1 s_2\cos\alpha},
\qquad c_i = \cosh(\lambda_i/2),\; s_i = \sinh(\lambda_i/2),$$
verified exactly in the rotor algebra and cross-checked entrywise against
$4\times 4$ Lorentz matrices **[proved; retrodiction: Wigner]**. The margin
$s_1 s_2 \sin\alpha$ — vanishing iff collinear or unpivoted — is the first of
five appearances of a single motif: every forced deviation in this model
carries a $\sinh$-product margin with an exact vanishing locus (the others:
the thermal identification's off-office leak, the medium-current's flux, the
uncertainty excess of misaligned quadratures, and the dissipation gate of §8).

**The Kähler dial.** (Figure 3.) The frame space $(S^2)^3 = (\mathbb{CP}^1)^3$ is Kähler;
the complex structure at axis $a$ acts on tangent vectors as
$J_a v = a \times v$, whose generator is the bivector $P_a = I a$ — the
axis's own plane. Since $\omega(u,v) = g(J u, v)$ exactly, Hamiltonian flow is
the gradient flow rotated by one application of $J$: the reversible
"circulation" dynamics and the dissipative "settling" dynamics driven by the
*same* potential differ by a quarter-turn of the generator,
$\dot{x} = -e^{J\theta}\,\mathrm{grad}\,V$, with descent rate
$\cos\theta$ **[proved]**. At $\theta = \pi/2$ the flow conserves the
potential and the total axis vector exactly; at intermediate $\theta$ the
modes carry complex frequency — oscillation with decay, i.e. *evanescence*,
which is the model's own frozen name for the $\hbar$-seat's second pole
[lineage: Landau–Lifshitz–Gilbert; metriplectic/GENERIC; imaginary-time
propagation].

![Figure 3](fig3_kahler_dial.png)

### 4. The readout layer at the $c$ seat

**Presentation families.** A pivoted seat reading an isotropic reference does
not read a number; it reads a *family* over directions,
$s(\theta) = \cosh\lambda + \sinh\lambda\cos\theta$, with the direction
measure forced to be the plane measure $du/2$, $u = \cos\theta$, by the
Jacobian of the pivot ("the measure is the plane") **[proved]**. Power-mean
aggregates $M_p = T_0[\langle s^p\rangle]^{1/p}$ over antipodal poles
reproduce the three historical relativistic-temperature laws at
$p = -1, 0, +1$ (Planck–Einstein, Landsberg, Ott) **[retrodiction]**; the
underlying unification — one family, one selector — has explicit ancestry in
Biró and Ván (2009–2010), who select by an internal heat-current assumption
where this model selects by a (measure, exponent) pair. The model's own
additions at this seat are exact identities: the mirror
$\langle s^{p}\rangle_{\text{piv}} = \langle s^{-p}\rangle_{\text{unpiv}}$
and its continuous generalisation $I(p, m) = I(-p, 1-m)$, whose fixed point
is the *half seat* $m = \tfrac12$ — the frame in which bath and observer move
equally and oppositely — where the whole power-mean spectrum is
self-reciprocal and the geometric mean equals the invariant by symmetry alone
**[proved; novelty status: unlocated in the literature, elementary]**.

**The invariant rides the pair.** No single pivoted reading returns the
invariant; the pairing of the presentation with the thermal identification
does — and the identification vector is, up to $\hbar/k_B$, the
inverse-temperature four-vector $\beta^\mu = u^\mu/T$ of van Kampen and
Israel **[retrodiction]**. The frame volume itself presents as the time
component of a four-current $n^\mu$: under a pivot the density reads
$\gamma J$ with a forced flux leak $\sinh\lambda\, J$ along the pivot axis,
and the seat-independent object is the current *norm*,
$(n^0)^2 - (n^1)^2 = J^2$ exactly **[proved]**. The medium of the model is
therefore a current, not a scalar; "nothing is created or destroyed" is a
continuity statement $\partial_\mu n^\mu = 0$, a conclusion forced
independently three times (by the transformation law, by the failure of any
scalar conservation in both minimal dynamics of §8, and by the additivity
structure of $\log J$, whose shared term is exactly the total correlation
$-\tfrac12\log\det G$ of the three axes) **[proved]**. Temperature-side
(effort, $\beta$) and medium-side (flow, $n$) form a conjugate pair with the
same $\sinh\lambda$ margin; their contraction is an action whose
dimensionless part is a Boltzmann exponent.

**The readout tier is an invariant-reduction exact sequence.** In logarithmic
coordinates the readout account $(\log T_+, \log T_-)$ carries a linear
reduction $\Sigma$ (the sum) whose value $2\log T_0$ is pivot-independent,
while the pivot moves the account strictly inside $\ker\Sigma$, the zero-sum
line spanned by $(1,-1)$ — so *the pivot is exactly the gauge direction*
**[proved]**. The antipodal identity $T_+T_- = T_0^2$ is this reduction written
multiplicatively, and the medium's $\log J$ is an additive account of the same
type. The readout layer thus instantiates the schema
$$0 \to \ker\Sigma \to \mathrm{Account} \xrightarrow{\ \Sigma\ }
\mathrm{Invariant} \to 0, \qquad
\text{account} = \text{invariant} + \text{zero-sum gauge residue},$$
with one structural difference from the additive case: here $\Sigma$ is linear
only after taking logarithms, because the reduction is multiplicative. That
difference is not incidental — it is the pairing theorem again, since the
geometric mean is the unique power mean that is multiplicative and hence the
unique one that is a $\Sigma$-reduction in log coordinates.

### 5. The $\hbar$ seat: the floor

The $\hbar$-seat pivot is a squeeze — the same hyperbolic rotor acting on a
quadrature plane. The reading family is
$q(\phi) = A(\cosh 2\lambda + \sinh 2\lambda\cos 2\phi)$, i.e. the $c$-seat
family pulled through the $2{:}1$ cover with doubled rapidity
($\Theta = 2\phi$, $\Lambda = 2\lambda$) **[proved]** — the orientation
double cover surfacing in the readout tier. The protected product
$\det C = A^2$ is preserved *identically* by the pivot group
($\det S = 1$, not an on-shell accident), and the zero stratum is a separate
orbit reachable only at infinite pivot: **"cannot vanish" is orbit
separation**; the geometry protects the product, the label supplies the
floor's value **[proved]**. The conjugate-product bound
$q(\phi)q(\phi + \pi/2) = A^2\!\left(1 + \sinh^2 2\lambda\,\sin^2
2\phi\right)$ recovers the Robertson–Schrödinger structure with saturation
exactly on principal axes **[proved; retrodiction]**; the pivoted moments are
Legendre polynomials, $\langle q^n\rangle = P_{n-1}(\cosh 2\lambda)$, and the
quadratic channel reads the squeezed-vacuum occupation $\sinh^2\lambda$
**[proved; retrodiction]**. Two structural discharges exceed retrodiction:
(i) the **blind-mass theorem** (W.L., 2026-08-28): if the reading family
satisfies $d\mu'/d\mu = f^{-n}$ on the seat's carrier, then
$\langle f^{n}\rangle_{\mu'} = \int d\mu = 1$ — the pivot-blind moment is
the *total mass of the source measure*, its exponent is the Jacobian power by
construction, and its value is $1$ for exactly that reason. The premise is
verified exactly per seat ($du'/du\,s^2 = 1$; $d\phi'/d\phi\,q = 1$), the
blind exponent is unique because only the degree-zero object of
Radon–Nikodym calculus is measure-free, and the Legendre index shift
$\langle q^{n}\rangle = P_{n-1}$ of this section is the Jacobian eating one
power of $q$ **[proved]**. This retires the former transport-exponent
conjecture; the surviving open question is sharper (Q-RN, §10): *why* the
reading is a power of the Radon–Nikodym derivative at all — the readout as
geometric rather than instrumental; (ii) the cover dictionary itself.

### 6. The $G$ seat: the horizon as the pinned wall

The $G$ seat runs the identical rotor with its rapidity **pinned by
position**: $\tanh\lambda(r) = \sqrt{r_s/r}$, the static seat's rapidity
against the local free-fall frame **[declared: KIN-2, the $G$-label's
operational datum; its derivation from a curvature tier is named future
work]**. Then $\mathrm{sech}^2\lambda(r) = 1 - r_s/r =: N^2$ exactly: the
$\hbar$ seat's exponential wall, previously unreachable because $\lambda$ was
free, acquires an *address* because the pinning diverges at finite $r = r_s$.
The reciprocal pair is (bound, escaping) $= (1+v,\, 1-v)$ with product
$N^2$; at the horizon the escaping pole alone closes,
$(1-v) = N^2/(1+v) \to 0$, while the bound pole saturates at $2$: the product
vanishes because one pole dies **[derived | KIN-2]**. The pair invariance
$T_{\mathrm{loc}}\cdot N = T_0$ holds for all $r > r_s$ while each factor
separately diverges or dies — the readout layer's pairing pattern standing at
the horizon [comparison-stage name: Tolman].

![Figure 2](fig2_trichotomy.png)

**The trichotomy theorem.** One rotor, three relationships between the pivot
parameter and the state, preregistered as a three-way character taxonomy
(frozen 2026-08-17, before any of the three was computed): the vanishing
locus lies **on-state** at the $c$ seat (null rays exist: *free*); at
$\lambda = \infty$ with $\lambda$ **free** at the $\hbar$ seat (*floored*);
at $\lambda = \infty$ with $\lambda$ **pinned** to a finite address at the
$G$ seat (*horizon-vanishing*). All three characters landed
**[proved / derived | KIN-2]**. The asymptotic wall constant is identical
across the $\hbar$ and $G$ legs ($e^{2\lambda} \times \text{product} \to 4$),
exhibiting the horizon as the floor's wall made reachable rather than as a
new mechanism.


**The pinned family as a state (THM-H, 2026-09-02).** Applying the pinned
rotor to the frame itself rather than to a probe presents the seat's rods:
$G' = G + \sinh^2\!\lambda\, k k^{\mathsf T}$ with $k_i = a_i\cdot K$ and
$\det G' = \Delta\cosh^2\lambda$ **[proved]**, the frame's own Gram a
spectator. Under the pinning, continued inside through
$\lambda = \mu + i\pi/2$ **[declared: CONT-1]**, both sheets collapse to one
rational formula, $G'(r) = G + [r_s/(r-r_s)]\,k k^{\mathsf T}$,
$\det G' = \Delta\, r/(r-r_s)$: a simple pole at the horizon, a simple zero at
the centre, nothing else **[derived | KIN-2a, CONT-1]**. Pulled back to the
frame, $G'$ is exactly $\delta_{ij} + [r_s/(r-r_s)]\,\hat x_i\hat x_j$, the
static observer's spatial rods, while the frame — the free-fall observer —
keeps flat rods [comparison-stage name: Painlevé–Gullstrand]. The hole is
therefore a *presented* state, not a locus in the frame: the horizon is where
the presented directions collapse onto the pivot (a node of the presented
state), the centre is where the presented state is coplanar (world layer) or
where the root's rod alone passes through zero (seat layer). RULING-2 assigns
$K\propto a_1+a_2+a_3$ to the hole and $K=$ root to the infalling seat
**[declared]**. KIN-2 is now split: KIN-2a, the pinning, remains declared;
KIN-2b, the curvature tier, was closed on 2026-08-30 by the Gram-submersion
curvature.

**What sits at the centre (CURV-1, 2026-09-02).** The founding question — is
"no singularity" a theorem or a blindness — was put to three readings of
curvature that use angles and rotors only. (i) The seat-triangle excess of the
presented directions, computed as the branch-free orbit over the four pole
choices, with the exact identity
$\prod_\epsilon(S_\epsilon + iV) = -4\prod_{i<j}(1-\Gamma_{ij}^2)$ on both
sheets **[proved]**: finite at the centre on both layers. The infalling seat
reads $\{\varphi_3,\ \pi-\varphi_3,\ \pi-\varphi_3,\ \varphi_3\}$ with
$\varphi_3$ the frame's dihedral angle at the root **[proved, every admissible
Gram]** — the root's presented vector vanishes, the visible axes project to its
perpendicular plane, and the presented triangle is the polar one whose excess
is its apex angle. The hole reads $\{2\pi,0,0,0\}$: the presented vectors sum
to zero ($G'(0)\,\mathbf 1 = 0$ **[proved]**), surround the origin, and the
triangle is the hemisphere, whose spin lift is $-1$: the deck $\mathbb Z_2$.
(ii) The loop holonomy of the two-seat relative rotors is exactly flat — every
contractible loop returns $1$ at every radius — and every ring around the
centre returns exactly $-1$ in $SU(2)$, $r$-independent, shrinking limit
included, conditional on geodesic frame transport (a twisted transport breaks
it) **[proved | TRANSPORT-1]**. (iii) The frozen candidate readings — the
Gram-submersion vanishing order, the inertia node detector, the Kummer
radicands — agree: finite at the centre on both layers; the one divergent
candidate, the Gram-submersion curvature on the world layer, diverges as
$\Delta_{\mathrm{pres}}^{-1/2}\sim r^{-1/2}$ because the presented state
reaches the coplanar boundary of the elliptope there. Fork (b): no rod-free
curvature reading diverges at the centre, which is a finite event with two
addresses, an angle on the seat layer and a spin sign on the world layer. The
divergence exists, but inside: at exact interior radii where two presented
directions become antipodal — on the seat layer $r^* = r_s\Delta/\sin^2\theta_{12}$,
below both null radii $r_i = r_s\sin^2\theta_i$ — the excess has logarithmic
singularities with no rod in any denominator **[derived | KIN-2a, CONT-1]**;
no standard analogue is claimed. Two rulings are owed: whether the frozen
rod-free definition ("no division by a presented length") admits normalised
cosines — the raw functional obeying its letter reads $0$ at the centre for
every admissible frame — and whether "identity" is baseline-relative or
mod $2\pi$.

**Light and matter in the pinned family (THM-I, 2026-09-03).** Two symmetry
pairings — a ray's frequency with respect to the frame's common time office,
$\omega_t = \omega(1-\tanh\lambda\cos\alpha)$, and its pairing with the
rotation generator, $L = \omega r\sin\alpha$ — are constant along a ray
**[declared: STAT-1, ROT-1]**; with the null locus and the drift bookkeeping
**[construction: PROP-1]** they close the orbit as an identity valid for *any*
radial profile,
$(dr/d\varphi)^2 = r^4/b^2 - r^2(1-\tanh^2\lambda)$, and for a massive body
$(dr/d\varphi)^2 = (r^4/\tilde L^2)\big[\tilde E^2 - (1-\tanh^2\lambda)(1+\tilde L^2/r^2)\big]$
**[proved]**. KIN-2a turns $1-\tanh^2\lambda$ into $1-r_s/r$ [comparison-stage:
Schwarzschild's null and timelike geodesic equations]. Receipts: gravitational
redshift $\cosh\lambda_2/\cosh\lambda_1$ with the sign forced by the frozen
Doppler convention and weak-field coefficient $\tfrac12$ [Pound–Rebka]; light
deflection $2r_s/b + (15\pi/16)(r_s/b)^2 + \dots$, both coefficients to $10^{-6}$
[Einstein 1915; Epstein–Shapiro]; perihelion coefficient $3.0003$ [Mercury].
The factor 2 has an address: the rod office supplies $1-r_s/r$ and the time
office the drift; each carries half, and the model cannot get it wrong because
a seat's rods and clocks are one boost ($\rho_K N^2 = 1$ for every profile
[PPN $\gamma=1$; Cassini]). The local transport law is the momentum paravector
boosted by the *gradient of the pivot field along its own co-moving step*,
$d\vec\lambda = (\hat n\cdot\nabla)\vec v\,dt$: it conserves both pairings
identically for every profile and integrates to the exact orbit; the two
seat-to-seat rules available in the corpus — the naive product of seat rotors
(the D1 non-collinear composition) and the frame-transported relative rotor —
both fail with order-one errors, the first by treating free-fall frames at
separated points as boosted copies against a background **[proved]**
[comparison-stage: the spin connection of the rain tetrad; Hamilton–Lisle's
river; Killing conservation]. The stopwatch objection of §9 is thereby a
theorem with a positive half: seats relate to the frame; frame-to-frame
transport needs the field's derivative, which carries no background.

**What the pinning is, and what it is not (THM-J, 2026-09-03).** Profiles that
agree with KIN-2a at first order differ at second, along two lines found and
then checked: perihelion coefficient $3 + 2c_2$ and second-order deflection
$15\pi/16 + (3\pi/4)c_2$, with $c_2$ the $(r_s/r)^2$ term of $\tanh^2\lambda$
[the PPN $\beta$-dependence, re-found]. Mercury selects $c_2 = 0$: the harmonic
variable is $\mathrm{sech}^2\lambda = \Delta/\det G'$, the reciprocal presented
volume in frame units, and the field law the data demand reads
$\nabla^2(\Delta/\det G') = 0$ outside the source with
$r^2\,d(\Delta/\det G')/dr = r_s$, the Euler–Lagrange equation of
$\int|\nabla(\Delta/\det G')|^2$. Whether the bare tier produces it: any
gradient-quadratic static energy makes its own arc length along the family
harmonic **[proved]**; an affine function of $\mathrm{sech}^2\lambda$ has a pole
at the presented volume's zero and is bounded at the horizon, which no
polynomial invariant of the presented cell is **[proved]**; and the
harmonic-map energy of the frozen elliptope metric of KIN-2b gives perihelion
coefficients $1.50$–$1.87$, frame- and layer-dependent **[derived, and
falsified]**. The field law is therefore not derivable from the bare tier's
*geometry*. It is derivable from the tier's *temperature*, which is the next
block.

**The pinning from the seat's temperature (THM-K, 2026-09-03).** For any
radial profile the static seat's proper acceleration is
$\alpha = (c^2/2)\,|d\beta^2/dr|/N$ **[derived | KIN-1, RULE-2, H-5, E-4]**, so
$\alpha N$, the acceleration seen from infinity, is the gradient of
$\mathrm{sech}^2\lambda = 1-\beta^2$ and of no other function of the rapidity:
whatever sets the seat's redshifted temperature sets the gradient of the
harmonic variable, and Part B's rival exponents are unreachable through this
door. The seat's Unruh temperature (D2, the $2\pi$ frozen), Tolman-redshifted,
is $(\hbar c/4\pi k_B)\,|d\beta^2/dr|$, Hawking's at the wall; the first law
$dE = T_H\,dS$ holds exactly with $S = k_B A/(4\ell_P^2)$, so the area law's
$\tfrac14$ and D2's $2\pi$ are one consistent pair. The construction: the
seat's presented sphere at $r$ has area $4\pi r^2$ for every profile (the
presentation stretches only along $K$) **[proved]**; give it the entropy
$k_B A/(4\ell_P^2)$ **[declared: SCREEN-1]**, equipartition **[declared:
THERM-1]**, the blind mass as the enclosed energy **[declared: MASS-1]**, and
equilibrium with the seat's horizon at the temperature seen from infinity
**[declared: EQ-1]**. Then $d\beta^2/dr = -2GM/(c^2 r^2)$, and with
$\beta^2 \to 0$ at infinity, $\beta^2 = r_s/r$ with $r_s = 2GM/c^2$ falling out
rather than put in: **KIN-2a derived**, together with the Gauss law
$r^2\,d(\mathrm{sech}^2\lambda)/dr = r_s$ and Poisson for the reciprocal
presented volume with the blind-mass density as source. The one fork inside —
equilibrate at the local temperature instead — gives $d\beta^2/dr =
-(r_s/r^2)\sqrt{1-\beta^2}$, the lapse harmonic, Part B's dead row (perihelion
$2.5$): Mercury decides it, in favour of Tolman's condition. Every result of
THM-H, CURV-1 and THM-I is thereby conditional on four thermal declarations
that say one thing — a seat's presented sphere carries $A/\ell_P^2$ thermal
degrees of freedom — instead of on an imported profile [comparison-stage:
Verlinde's equipartition argument and Padmanabhan's holographic equipartition,
descended from Jacobson's Clausius derivation; Tolman–Ehrenfest equilibrium;
the standing criticisms of entropic gravity, e.g. neutron interferometry, attach
to this route and are not answered here]. This is the third time a
freeze-dated theorem selected the correct answer to a question it was not
written for: a wrong Unruh coefficient would have broken $dE = T_H dS$ and
returned $r_s$ wrong by the same factor.

**Rotation (THM-L, 2026-09-03).** The transport law extends to a presentation
field with vorticity, and the extension is forced: boost by the *symmetric*
strain along the co-moving step, rotate at the river's local angular velocity
$\tfrac12\nabla\times\vec v$ per unit frame time (rule G$'$). It conserves both
symmetry pairings identically for any stationary axisymmetric field
$\vec v = -\beta(r)\hat r + \omega(r)\,\hat z\times\vec r$ **[proved]**; the
irrotational rule fails on the same field with explicit drift; with the swirl
off it reduces to THM-I's rule. With the swirl profile
$\omega = 2GJ/(c^2 r^3)$ (the drift of the zero-angular-momentum free-fall
frame, first order in $J$), a gyroscope at rest is carried by the rotation part
alone and precesses at $(G/c^2 r^3)[3(\vec J\cdot\hat r)\hat r - \vec J]$
exactly, the vorticity being a vacuum dipole field; averaged over a polar orbit,
$GJ/(2c^2 r^3)$, $40.8$ mas/yr for Earth at $7027$ km [comparison-stage:
Schiff, Lense–Thirring; Gravity Probe B measured $37.2\pm7.2$]. An orbit's node
advances at $2GJ/(c^2 a^3)$ — strain and rotation acting together along a
path, ratio $1.0000$ to the first-order value, zero without the rotation term
[LAGEOS/LARES]. Frame dragging is the vorticity of the free-fall flow.

**The swirl is the boost of the pinning (THM-M, 2026-09-03).** At first order
in $r_s$ the static presented field is the seat's clock, $h_{00} = r_s/r$, and
rods, $h_{ij} = (r_s/r)\hat r_i\hat r_j$, with no shift. The sandwich of BARE-1
is a Lorentz map ($L^{\mathsf T}\eta L = \eta$) and the presented metric a
rank-2 tensor under it, so a mass element moving at $\vec w$ acquires the
shift $h'_{0i} = w_i h_{00} + w_j h_{ji} = (r_s/r)[\vec w + (\vec w\cdot\hat
r)\hat r]$: half from the clock, half from the rods **[proved]**. Superposed
linearly over a rigidly rotating shell (K-6's law is linear), the far field is
$(2G/c^2)(\vec J\times\vec r)/r^3$, so $\omega = 2GJ/(c^2 r^3)$: **KIN-3
derived**, coefficient $2 = 1 + 1$ **[derived | THM-K's inputs, BARE-1, linear
superposition]**. A theory presenting only the time office predicts half the
dragging — $20.4$ mas/yr for Gravity Probe B, outside by $2.3\sigma$ — so the
rods' half is measured, the same split that gave light bending its factor two.
The natural rotating-screen thermal principles fail on dimensions (a drift
$\propto J/M$, no $G$): the swirl is not a new thermal input. Scope: first
order in $r_s$ and in $J$; the full rotating solution is not reached.

### 7. Thermal coefficients from rotor periods

**Unruh.** With no field theory anywhere in the chain — no Bogoliubov
transformations, no mode functions, no Rindler quantisation — the boost
orbit's Wick face is a circle of period exactly $2\pi$
($\cosh i\theta = \cos\theta$), while the *rotor* closes only at $4\pi$
($R(2\pi) = -1$). The fork is live: a spinorial readout would select $4\pi$
and predict $T = \hbar\alpha/4\pi c k_B$, which is false. What selects
$2\pi$ is an internal theorem proved ten days earlier for an unrelated
purpose: the readout layer transports two-sidedly (tensorially), so
$(-1)X(-1) = X$ and the $-1$ is invisible to every readout object
**[proved]**. With the frozen thermal label (the thermal datum is a closed
imaginary-time circle of circumference $\hbar/k_B T$) and the kinematic
definition $\lambda = \alpha\tau/c$, the assembly gives
$$T \;=\; \frac{\hbar\,\alpha}{2\pi c\, k_B}\,,$$
the Unruh temperature with its coefficient **[derived | LBL-1, KIN-1;
retrodiction; rigorous ancestry: Bisognano–Wichmann, Sewell]**. Read as a
constitutive law — effort proportional to flow — the "resistance of the
vacuum" is $R = \hbar/2\pi c k_B \approx 4.05\times 10^{-21}\,
\mathrm{K}/(\mathrm{m\,s^{-2}})$.

**Hawking.** Importing only the static-observer acceleration $\alpha(r)$ at
the comparison stage, the redshifted acceleration limits to the surface
gravity, $\alpha N \to \kappa = c^2/2r_s$, and the same engine yields
$$T_H \;=\; \frac{\hbar c^3}{8\pi G M k_B}\,,$$
numerically $6.17\times 10^{-8}\,$K at one solar mass **[derived | KIN-2 +
imported $\alpha(r)$; retrodiction: Hawking 1974]**. The factor-two
falsifier of the Unruh chain is inherited unchanged.

### 8. Closure, emergent time, and the derived dissipator

**The ruling.** Five independent computations (the two minimal dynamics run
side by side, the failure of scalar conservation in both, the sign-definite
additivity deficit, the bath argument, and the Kähler phase dial) converged on
one question — is the structure closed? — which was answered by a principle
frozen twelve days *before* any of those computations existed: temporal
predicates attach to proper subsystems only; no sentence, affirmative or
negative, places the totality in a temporal slot. The binding form of the
ruling is therefore not "the totality sits at $\theta = \pi/2$" but "the
totality is not on the dial": $\theta$ is a seat predicate **[declared:
RULING-1, on freeze-dated ground; consonant with Wheeler–DeWitt,
Page–Wootters (with the Giovannetti–Lloyd–Maccone repair and photonic
demonstrations), the decoherence programme, and Connes–Rovelli thermal
time]**. The model-native sharpening — new relative to Page–Wootters — is
that the partition is not free: there are exactly three seats, labelled by
$c$, $\hbar$, $G$, permuted by the click group, which §2 showed is the
automorphism group of the branch surface **[conjecture, the sharpest new one
on the board]**.

**The seat dissipator, derived.** A closed world at $\theta = \pi/2$ owes an
account of why seats age. The bath cannot be imported, so it is built from
owned parts: the seat's environment is *the sky* — every other cell's
unit-frequency angle mode, Doppler-presented through the readout, whose
frequencies under the forced plane measure land in an exactly **flat band**
$[\omega_0 e^{-\lambda}, \omega_0 e^{\lambda}]$ **[derived]**. Independent-
oscillator elimination [lineage: Ford–Kac–Mazur; Zwanzig; Caldeira–Leggett]
then gives the seat a generalized Langevin equation whose damping kernel and
noise carry the *same* coefficients — the fluctuation–dissipation identity,
verified symbolically: a seat that drinks must hiss **[proved]**. The gate:
at $\lambda = 0$ the band is degenerate, the kernel is an undamped cosine,
and no arrow exists — dissipation *requires* the pivot, margin
$\sinh\lambda$ (the motif's fifth appearance) **[proved]**. In full unitary
simulation the world's Hamiltonian is conserved to $10^{-7}$ while the seat's
energy relaxes to a floor equal to $k_B T$ (the hiss), a cold bath relaxes
toward zero, a degenerate band never relaxes, and a small bath *revives* —
the arrow is asymptotic in the crowd size, exactly as the ruling requires
**[numeric evidence]**. Residue: the per-mode coupling $c(\omega)$ is here
declared Ohmic-in-band, not derived; its derivation from the model's channel
law is the named open debt DEBT-2b.

### 9. The seat-cycle holonomy and the banked prediction

Lifting each axis to its minimal spinor (the half-angle map — the same
$2{:}1$ cover as §5's dictionary), the cycle
$B = \langle a_1|a_2\rangle\langle a_2|a_3\rangle\langle a_3|a_1\rangle$
admits the exact closed form
$$B \;=\; \frac{1 + \gamma_{12} + \gamma_{13} + \gamma_{23} + i\,V}{4}
\quad\textbf{[proved; verified to } 10^{-31}\textbf{]},$$
so that the Gram trace is its real part and the trivector its imaginary part:
the two invariants of the state are the two components of a single complex
number. Its phase is half the *signed solid angle* of the frame triangle
[lineage: Bargmann; Pancharatnam; Mukunda–Simon], recovering the classical
form $\tan(\Omega/2) = V/(1+\sum\gamma_{ij})$ [Van Oosterom–Strackee] as the
ratio of the two parts.
The phase generator is $V$ — **the medium is the phase** — with four verified
structural consequences: the sign of the phase is the orientation-cover
coordinate (an operational meter for the sheet, which §2 proved no continuous
function can supply); the meter nulls on the branch locus; seat circulation
preserves $B$ while reversals conjugate it (the click group's
untwisted/twisted split realised on a measurable number); and the branch
locus is $\mathbb{Z}_2$-graded — coplanar frames read $0$ on one side of the
symmetric ("trine") configuration and exactly $\pi$ on the other, with the
trine as the wall between phase classes **[proved / verified to 30 digits]**.

**Descent theorem** *(lineage located 2026-08-28, same night: this is
Li--Wagner--Zhang Thm V.1 — every single-qubit Bargmann invariant satisfies a
monic quadratic over the overlap field, determined up to one global
conjugation (arXiv:2507.14878); the cycle-space/recovery side is
Chien--Waldron (arXiv:1312.5393, SIAM JDM 2016); the $n=3$ closed form is
Zhang--Xie--Li, PRA 111, 042417 (2025). Their proof's property (A) is the
Binet--Cauchy engine used below. Everything in this subsection is therefore
**[retrodiction]** as mathematics; what remains ours is the instantiation —
$\arg B$ as the orientation-sheet meter of THIS cover, tied to the Cayley
$\Delta$ and the click-equivariant Kummer module.)* Each individual overlap $\langle a_i|a_j\rangle$ requires
its own radical $\sqrt{1+\gamma_{ij}}$, so a single leg of the cycle lives in
the degree-eight multiquadratic field of §2. The *closed cycle*, by the closed
form above, is rational in the Gram data together with $\sqrt{\Delta}$ alone:
the three conjugate radicals cancel around the loop and only the click-invariant
one survives, so the seat-cycle holonomy descends to the degree-two orientation
subfield $F(\sqrt{\Delta})$ **[proved]**. This upgrades the sheet-meter
property from an observed sign flip to a field-theoretic fact: $\arg B$ carries
exactly the orientation square class and nothing else, which is why no
continuous global branch of $V$ can supply what the phase supplies.

**PREDICTION-1 (banked; protocol assigned 2026-09-04).** $J = (1 + \gamma_{12} + \gamma_{13} +
\gamma_{23})\tan(\Omega/2)$ with $\Omega = 2\Phi/(2s)$, minimal lift
$s = \tfrac12$. The claim is stated as a separation between a *presentation*
and an *operation*, because that separation is what the model is for.

*Presentation.* The unrestricted pivoted observation map of one root —
VIEW-1's circumaxis as the pivot axis, the co-pivot rotor two-sided on the
axes and one-sided on the spinor roots, the roots geodesic from the
observer's own axis, and every bilinear bar-pairing among the presented
objects — was derived with its specification frozen before any output was
read (`pivot_map.py`). Its scalar office is deck-even for every pairing of
axes; two states never yield a pseudoscalar (the state lemma recovered); and
no reading of definite odd parity exists anywhere in it. What it does contain
are *mixed* readings on the transition pairings, with closed forms
(`pivot_map_closed_forms.py`): the axis–transition pseudoscalar
$\langle a_k\,\bar T_{ij}\rangle_I = -(E_{kij} + \epsilon_{kij}V)/(2(1+h))$,
$h = V/2|A|$, $E_{kij} = \hat n\cdot((a_j - a_i)\times a_k)$, whose deck-odd part
$(E_{kij}h - \epsilon_{kij}V)/(2(1-h^2))$ is $V$ times a Gram-rational function.
The transition itself has the closed form $T_{ij} = R_i\widetilde R_j$ with
$R_i = (1 + a_i\hat n)/\sqrt{2(1+h)}$ the shortest rotor from the observer's
axis to the root — undo $R_j$, apply $R_i$ — so that
$T_{ij} = [(1+2h+\gamma_{ij}) + i\,m_{ij}\cdot\sigma]/(2(1+h))$,
$m_{ij} = \hat n\times(a_j - a_i) + a_i\times a_j$. This is *not* the direct
geodesic rotor from $a_j$ to $a_i$: the route through $\hat n$ is a closed loop
$a_j \to \hat n \to a_i \to a_j$, and the two differ by a rotation about $a_j$
whose angle is the solid angle of the pole triangle $(\hat n, a_i, a_j)$ — the
loop's holonomy. The pair (scalar part of $T_{ij}$, $\hat n$-component of
$m_{ij}$) $= (1+2h+\gamma_{ij},\ \hat n\cdot(a_i\times a_j))$ is exactly that
triangle's Van Oosterom–Strackee pair **[proved on 150 frames, exact
identities]**. "Dressed by the observer's own axis" therefore has a closed
form: every transition an observer at one root presents is the direct geodesic
twisted by the solid angle its axis makes with the two roots, and the mixed
readings inherit that twist. The parity bookkeeping is exact and was first
stated wrongly here: the deck sends both roots to their antipodes, so the pole
triangle becomes $(\hat n, -a_i, -a_j)$ — its numerator $\hat n\cdot(a_i\times a_j)$
is deck-*even* and its denominator changes from $1+2h+\gamma_{ij}$ to
$1-2h+\gamma_{ij}$ — so the pole triangle does *not* flip sign; the dressing is
mixed through $h$ alone. The deck-odd content of every reading is the third
axis's component of $m_{ij}$, $a_k\cdot(a_i\times a_j) = \pm V$, the plain
triple product **[signed checks, 150 frames]**. The exact statement: $V$ is the
only primitive deck-odd generator — as the descent theorem requires, the odd
square class is $F(\sqrt\Delta)$ and has one generator. The observer-axis
geometry contributes no orientation datum of its own; after
antisymmetrisation it supplies a deck-*even*, seat-dependent coefficient
multiplying $V$, $[E_{kij}/2|A| - \epsilon_{kij}]/(2(1-h^2))$, which takes both
signs across the cell **[300 frames]**. So the seat changes the gain, and can
even flip the measured sign relative to $V$ through that even coefficient, but
it supplies no independent deck bit: it colours the reading and never signs it. The direct geodesic is what a pairwise measurement sees; $T_{ij}$
is what an observer at one root sees; their difference is the observer's
position on the view sphere, written as a solid angle. That is the seat, as
geometry.
So the model's *presentation* already separates a state from its deck image
— $\Pi(X) = \Pi(\tau X)$ on the Gram, $\mathcal O(X) \neq \mathcal O(\tau X)$
on the transitions — but reads the orientation only dressed by the observer's
own circumaxis, never cleanly. The orientation reference the presentation
carries is the sign of that circumaxis, i.e. the cyclic order in which the
three alignments are named. **[derived; the only declaration is the naming
order]**

*Operation.* A definite-parity reading of $\operatorname{sgn} V$ exists only as
an *ordered* operation on the presentation: at three ordered pivots the
observer records the presentation of the aligned axis into an ordered
register, $\varrho_X(i) = P(\Pi_{s,\theta_i}(X)) = (1 + a_i)/2$, and the
directed permutation $U_B$ of the three registers reads
$\operatorname{Im}\operatorname{Tr}[U_B\bigotimes_i \varrho_X(i)] = V/4$
(`pred1_operational_chain.py`). Two requirements are checked rather than
assumed: the pivot selects the presented axis without rotating the frame the
registers share (storing each presentation in its own pivoted frame breaks the
reading; identical registers read exactly $0$), and the source is unchanged by
observation, $X \mapsto X$ absolutely, so three ordered presentations of it can
be taken. The direction of the test is the temporal order of the
observations: reversing the order reads $-V/4$, cyclic re-orderings read the
same. Three physical readouts are exact (`pred1_physical_protocol.py`): the
controlled-$U_B$ Hadamard test ($\langle Y\rangle = V/4$), the spin-chirality
operator $\hat\chi = \sigma_1\!\cdot(\sigma_2\times\sigma_3)/2\sqrt3$ with
$U_B = e^{+i(2\pi/3)\hat\chi}$ exactly and $\langle\hat\chi\rangle = (2/\sqrt3)\operatorname{Im}B$,
and the three-photon tritter, $\operatorname{Im}B = \pm\tfrac{3\sqrt3}{2}\,(P_{(210)} - P_{(120)})$
with the sign fixed by the tritter's chirality. Pairwise SWAP/HOM tests return
$\Gamma_{ij} = 2\operatorname{Tr}(P_iP_j) - 1 = \gamma_{ij}$ and
$|\operatorname{Im}B| = \sqrt{\det\Gamma}/4$ only. In every realisation the
sign needs a directed device — a permutation order, a control line, a handed
splitter. **PROTOCOL-1**: observations are ordered; presentations are not.
That sentence is not an added axiom: it is the computed difference between the
bilinear map above and the ordered operation here.

*Algebraic half, certified three ways* (`thm_galois_deck_descent.py`):
the valuation matrix of $(\Delta, 1+\gamma_{ij})$ along the four boundary
divisors is the identity (rank-4 square classes); the degree-16 norm
polynomial at a rational point is irreducible; the minimal polynomial of $B$
over $F = \mathbb Q(\gamma)$ is $16x^2 - 8(1+\sum\gamma)x + 2\prod(1+\gamma_{ij})$,
discriminant $-64\Delta$, so $B \in F(\sqrt{-\Delta}) = F(i,\sqrt\Delta)$; the
four Cayley nodes are the rank-one Grams, $A_1$ points at which $\Delta$ has
multiplicity two, so a loop around a node does not flip $\sqrt\Delta$ while two
leg radicals do. The deck is the inertia of the smooth Cayley cubic and of
nothing else.

*Status, honestly.* The two "cycle" protocols first written for this
prediction (`pred1_protocol.py`: ordered idempotents; geodesic co-pivots
carrying the spinor root around the triangle) were built with the answer known
and are recorded as an **existence witness**, not a derivation; the cold map
does not produce them, because ordered products of idempotents are not among
the record's admitted operations. The cross-seat *content* of the prediction
is unchanged: the $c$-seat account (the even presentation — GR's account) cannot
read the orientation of three seat axes; the ordered $\hbar$-seat operation
can; the difference is the model's content. The lift charge $s$ and its tie to
the $\hbar$ quantum remain conjecture. Remainder: apparatus wiring — repeated
runs or stored temporal modes feeding three ordered registers.

### 10. Honest accounting

**Classical, with named lineage (the long list).** The Cayley cubic and its
automorphisms; the elliptope; double-cover monodromy of nodal hypersurfaces;
Wigner rotation; the binary-octahedral spin lift; Doppler, aberration,
Landsberg–Matsas; the temperature-family unification (Biró–Ván); the
inverse-temperature four-vector (van Kampen, Israel); Tolman; Robertson–
Schrödinger and squeezed states; Kähler geometry and Landau–Lifshitz–Gilbert;
GENERIC/metriplectic; Ford–Kac–Mazur/Zwanzig/Caldeira–Leggett and classical
FDT; Wheeler–DeWitt; Page–Wootters (+ Giovannetti–Lloyd–Maccone); decoherence
(Zeh, Zurek); Connes–Rovelli; Barbour; Bisognano–Wichmann, Sewell; Unruh,
Davies, Fulling; Hawking; Bargmann, Pancharatnam, Mukunda–Simon, Van
Oosterom–Strackee; Bronstein's $cG\hbar$ cube; Painlevé–Gullstrand and the Hamilton–Lisle river; Killing conservation; the spin connection of a tetrad; Einstein's deflection, Epstein–Shapiro's second-order term, Pound–Rebka, Cassini ($\gamma$), Mercury ($\beta$); Gauss–Bonnet and the polar triangle; the mechanical connection (Littlejohn–Reinsch); harmonic maps and sigma models; Jacobson's Clausius derivation, Verlinde's equipartition, Padmanabhan's holographic equipartition, Tolman–Ehrenfest; the Bekenstein–Hawking first law; Schiff and Lense–Thirring, Gravity Probe B, LAGEOS/LARES; Doran's form of Kerr and the Hamilton–Lisle twisting river; linearised Kerr and the spin-2/scalar distinction.

**Inherited from the author's own separate corpus (a distinct category).**
Three structural frames used above come not from the literature but from an
unrelated body of work by the same author (the Cella DBP programme on
role/channel/orbit geometry and its Galois–Kummer companion papers), predating
this project and frozen with content hashes: the *canonical invariant reduction
theorem* ("account = invariant + zero-sum gauge residue"), which §4 shows the
readout tier instantiates; the *conjugate Kummer module* and odd-valuation rank
criterion, which supply §2's cover group; and the wreath-closure family of which
§2's $(C_2\wr S_3)\times C_2$ is an instance. This inheritance is genuinely
useful — it replaced hand-rolled arguments with general theorems — but it must
be graded differently from external lineage: shared authorship means shared
assumptions and shared blind spots, so an apparent convergence between the two
corpora is weaker evidence than an independent confirmation would be. Readers
consulting both should also note a false friend: "channel" in the DBP corpus
means a curvature-channel component ($\kappa_c,\kappa_s,\kappa_{int}$), whereas
"channel law" here means a detector response functional. Same word, unrelated
referents.

**New or unlocated, honestly sized (the short list).** The half-seat
reflection $I(p,m) = I(-p,1-m)$ and the all-$p$ circle mirror (elementary,
unlocated); the trichotomy theorem (one rotor, on-state/free/pinned) as the
unifier of the frozen three-character taxonomy; the escaping-pole anatomy of
horizon-vanishing; the tensorial selection of the thermal $2\pi$ (a proven
internal theorem retiring a live factor-two error); the flat-band derivation
of the seat bath from the forced measure, and the $\sinh\lambda$ dissipation
gate; the $\mathbb{Z}_2$ grading of the branch locus by the seat-cycle phase;
the canonical three-seat partition menu as a sharpening of Page–Wootters
**[conjecture]**; the blind-mass theorem with its uniqueness clause (retiring
the former transport-exponent conjecture, no third carrier required); its
successor question Q-RN — derive the reading-as-RN-power premise from the
channel law, noting $n(c) = 2$ coincides with the quadratic channel exponent;
PREDICTION-1; the hole as a presented state and its one-formula interior; the centre's two readings (dihedral angle at the root; the hemisphere and its deck sign); the interior antipodal collisions with exact addresses $r_s\Delta/\sin^2\theta_{ij}$ **[candidate-new, lineage unsearched]**; the branch-free excess orbit and its product identity (elementary, unlocated); the arc-length theorem and the pole obstruction of THM-J; the finding that no seat-to-seat rotor rule propagates a ray, as a theorem; the K-1 identity that the redshifted acceleration is the gradient of the reciprocal presented volume for every profile (elementary, the static-observer relation in the model's variables); the resolution of the local-versus-redshifted fork in the equipartition argument by the perihelion; the clock/rods half-and-half split of frame dragging paired with the same split in light bending (each half known, the pairing unlocated); the dimensional failure of rotating-screen principles.

**Declared inputs and named debts.** LBL-1 (thermal circle), KIN-1 (proper
acceleration); KIN-2a and KIN-3 are now **derived** (THM-K, THM-M) from four thermal declarations that state one idea — a seat's presented sphere carries $A/\ell_P^2$ thermal degrees of freedom: SCREEN-1 (entropy $k_B A/4\ell_P^2$), THERM-1 (equipartition), MASS-1 (the blind mass is the enclosed energy), EQ-1 (equilibrium at the temperature seen from infinity); KIN-2b, the curvature tier, closed 2026-08-30; CONT-1 (interior continuation), RULING-2 (layer assignment), TRANSPORT-1 (geodesic frame transport), STAT-1 and ROT-1 (rays respect the pinned field's two symmetries), PROP-1 (drift bookkeeping); linear superposition (first order in $r_s$ and $J$); three rulings owed — the letter of the rod-free definition, the semantics of "identity", and whether EQ-1's Tolman choice is a principle the author owns or a fact Mercury owns; Ohmic-in-band coupling (DEBT-2b: derive $c(\omega)$ from the channel law); the Stage-3 role dictionary with its stabilizer-compatibility admissibility law; adjudication of the irreducibly-signed-observable hypothesis (H-B). Two
methodological facts deserve record. First, three times in this programme a theorem
frozen for one purpose later selected the correct answer to a question it was
not written for (the tensoriality theorem selecting $2\pi$; the closure rule
resolving the relax/circulate fork; the same $2\pi$ returning the correct $r_s$ in THM-K, where a wrong coefficient would have broken the first law) — the freeze-date discipline functioning
as the intended antidote to retrofit. This should be distinguished from the
cross-corpus inheritance catalogued above, which is a different event type:
preregistration working is evidence about *this* programme's discipline,
whereas one of an author's corpora supplying machinery to another is evidence
only that the two share an architecture — possibly because the architecture is
real, possibly because the author is. Conversely, the outcome ledger contains
as yet **no** class-(b) result against experiment: every physical output matches known physics. That statement now has a precise scope: given the four thermal declarations, the model's gravitational sector *is* Schwarzschild's geodesic structure plus slow-rotation Kerr (§6), so it can differ from the standard theory only beyond first order in $J$ or where the thermal declarations fail. The ledger contains, since 2026-09-03, its first falsified model-internal candidate — the bare tier's own state metric proposed as the dynamics, perihelion coefficients $1.50$–$1.87$, dead at Mercury — and two measured numbers that killed rival readings of the model's own construction (the local-temperature fork; the clock-only half of frame dragging). One candidate-new structure (the interior antipodal collisions) stands with no standard analogue claimed. A reformulation indistinguishable from the standard theory on all computed observables is either correct or unfalsifiable, and cannot tell which from inside; PREDICTION-1's cross-seat protocol and the thermal declarations' own tests (the entropic-gravity objections; rotation beyond first order) are the designated exits.

### 11. Methods and artifacts

All symbolic results are exact (SymPy), with adversarially chosen independent
second routes where feasible (hand-rolled Clifford algebra vs matrix
representations; spinor vs vector paths); numerics are mpmath at 20–50
digits or symplectic integration with conservation monitors reported inline.
Preregistrations and results carry SHA-256 hashes frozen before computation
in the AGNOSTIC-1 campaign; the present session's suites are committed
sequentially to the project repository (thm_g, thm_b_monodromy, cayley,
thm_d, doors, kahler, debt2, thm_d2_unruh, thm_g2, bargmann, thm_rn,
bridge_dbp, galois, correspondence_july; figures.py; RULING-1; PREDICTION-1; and, from 2026-09-02/03, thm_h, thm_h2_d1, curv1_path1, curv1_path23, thm_i_pre, thm_i_transport, thm_i_field, thm_j_dyn, thm_k_clausius, thm_l_rotation, thm_m_swirl), each with its full check output. The September suites use the numpy float path only for the ODE corroboration of the transport law; every claim there is also carried by an exact identity.
Head of the sequence at the time of writing: `68f3a9f` (THM-M), preceded by `d863efa` (THM-L), `5c58e8e` (THM-K), `2b03478` (THM-J), `a7a5350` (THM-I) and `90c0281` (CURV-1); at the previous writing `7e4e4e1` (Galois/Kummer
structure and the descent theorem), preceded by `685a89e` (invariant-reduction
bridge) and `2ba206f` (blind-mass theorem). Every FAIL encountered during
development is preserved in the transcript with its diagnosis; three were
substantive (a dead ladder conjecture; a non-resonance correction to the
seesaw mechanism; the branch-grading discovery), the remainder were
simplifier or harness defects, reported, not buried.
From 2026-09-04: thm_o_strain_law (33/33), view1_incidence (12/12), census_c (7/7), thm_galois_deck_descent (10/10), pred1_deck_separator (10/10), pred1_protocol (11/11; existence witness), pivot_map (6/6; spec frozen before output), pivot_map_closed_forms (10/10), pred1_physical_protocol (10/10), pred1_operational_chain (7/7); docs/2026-09-04-thm-target-O.md, the axis/planes/rods clarification and its addendum, HANDOFF-2026-09-04.md.



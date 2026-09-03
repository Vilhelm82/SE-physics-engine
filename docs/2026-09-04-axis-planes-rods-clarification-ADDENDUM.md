# ADDENDUM — ROOT SEATS, PIVOTED VIEWS, HIDDEN AXIAL MOTION, AND THE STATUS OF THE SUPERPOSITION CONJECTURE

**Purpose:** This document supersedes incorrect or overextended interpretations in `axis_planes_rods_clarification.md`. It records the model architecture clarified after that document, separates proved results from declarations and conjectures, and prevents conjectured measurement behaviour from entering the model as an input.

---

## 1. Status key

The following distinctions are mandatory:

- **[ontology]** — what the model says exists.
- **[proved identity]** — exact symbolic identity.
- **[declared law]** — rule inserted explicitly.
- **[derived | X]** — consequence conditional on declared input \(X\).
- **[conjecture]** — proposed interpretation or future consequence.
- **[prediction candidate]** — may be booked only if it falls out without being preconditioned.
- **[not run]** — proposed computation without a receipt.

---

# PART I — CORRECTIONS TO THE PREVIOUS CLARIFICATION

## 2. The fundamental planes never deform

This remains the model’s ontological statement:

\[
\boxed{
\delta \bar g_{P_i}=0
\qquad\text{for every fundamental plane }P_i.
}
\qquad\textbf{[ontology]}
\]

Any metric variation inferred by an observer is a variation of the **presented account**, not deformation of a fundamental plane:

\[
R\!\left[h^{(c)}\right]\neq0
\quad\not\Longrightarrow\quad
\delta\bar g_{P_i}\neq0.
\]

The phrase “the rods bend” is therefore admissible only as shorthand for:

\[
\boxed{
\text{the observer’s presented rod metric is nonflat.}
}
\]

It must never be used to assert ontic deformation of the underlying planes.

---

## 3. Non-invertibility does not itself imply curvature

The previous document correctly identified the seat map as non-invertible, but non-invertibility alone does not prove:

\[
R[h^{(c)}]\neq0.
\]

The actual open map is

\[
q_c
\longmapsto
h^{(c)}_{ij}
=
\mathcal F_{ij}[q_c,\partial q_c,\ldots],
\]

where \(q_c\) denotes the hidden position/orientation data of the \(c\)-axis representative within its owning plane.

The model must derive \(\mathcal F\). Curvature is then computed from its output.

Likewise,

\[
dE^a\neq0
\]

does not by itself imply curvature. A nonholonomic coframe can still describe flat geometry. For a presented coframe \(E^a\), one must compute

\[
dE^a+\omega^a{}_b\wedge E^b=0
\]

and then

\[
\mathcal R^a{}_b
=
d\omega^a{}_b
+
\omega^a{}_c\wedge\omega^c{}_b.
\]

Presented curvature requires

\[
\mathcal R^a{}_b\neq0.
\]

---

## 4. The tangential Einstein components are computable

The earlier wording could be read as saying that the observer cannot evaluate \(G_{ij}\). That is incorrect.

In the stationary unit-lapse presentation, once

\[
(h_{ij},v^i)
\]

is supplied, every component

\[
G_{ij}[h,v]
\]

is fully computable. O-6d explicitly computed them.

What is absent is not the tensor’s value as a mathematical functional. What is absent is a **model-owned law prescribing what that value must be**.

The current status is:

\[
G_{nn}=e_2(K)
\qquad\textbf{[proved identity on flat rods]},
\]

\[
G_{ni}
=
-\frac12(\nabla\times\nabla\times\vec v)_i
\qquad\textbf{[proved identity on flat rods, fixed convention]},
\]

followed by

\[
G_{nn}=0,
\qquad
G_{ni}=0
\qquad\textbf{[declared: NORMAL-1]}.
\]

No corresponding model principle has yet declared or derived:

\[
G_{ij}=\text{specified tensor}.
\]

GR supplies one possible law,

\[
G_{ij}=0
\]

in vacuum, but GR is not the model’s arbiter. O-6d establishes that GR’s tangential vacuum law cannot be satisfied by the complete flat-rod \(O(J^2)\) family. It does not establish that the model has failed.

The exact debt is:

\[
\boxed{
\text{derive the law governing the stationary tangential tensor }
G_{ij}[h,v].
}
\]

---

# PART II — AXES, PLANES, ROOTS, SEATS, AND VIEWS

## 5. An axis is an edge metric

Every axis is a plane seen edge-on from an adjacent plane.

The axis is not the plane itself and is not merely an intrinsic line drawn in the plane from which it is observed. It is the axial metric trace of its owning plane under the adjacent view.

For the \(c\)-plane:

\[
P_c
\longrightarrow
a_c
\longrightarrow
c,
\]

where:

- \(P_c\) is the hidden owning plane;
- \(a_c\) is its edge-on axial metric;
- \(c\) is the constant direct root reading when the observer is aligned with that axis.

Seating on the axial root makes the owning plane unavailable to direct observation. Motion within that plane therefore cannot be measured as a change in the root value.

---

## 6. One local root seat, many pivoted axial views

Treating the \(c\), \(\hbar\), and \(G\) views as three physical seats occupied in turn was a category error.

Let

\[
s=a_c\cap a_{\hbar}\cap a_G
\]

denote the observer’s local common root. The observer remains at \(s\) while the observation map pivots:

\[
\theta\in S^1.
\]

The complete local presentation is a family

\[
\Pi_{s,\theta}(X).
\]

The three constants correspond to distinguished axial alignments:

\[
\theta_c,\qquad
\theta_{\hbar},\qquad
\theta_G.
\]

Thus:

\[
\Pi_c=\Pi_{s,\theta_c},
\qquad
\Pi_{\hbar}=\Pi_{s,\theta_{\hbar}},
\qquad
\Pi_G=\Pi_{s,\theta_G}.
\]

These are three views from one local root, not three separate observer locations.

The legal operation is

\[
(s,\theta)\longmapsto(s,\theta+\varphi).
\]

This pivot is not ordinary spatial rotation in an assumed three-dimensional background. It changes which plane is presented edge-on through the fixed root.

Hence:

\[
\boxed{
\text{one root seat}
+
\text{one }360^\circ\text{ pivot orbit}
+
\text{three distinguished axial alignments}.
}
\]

---

## 7. Measuring a constant

When aligned with axis \(i\), the observer obtains its constant root value:

\[
\Pi^\parallel_{s,\theta_i}(P_i)=k_i,
\qquad
k_i\in\{c,\hbar,G\}.
\]

Hidden movement of that axis in its owning plane is annihilated in its direct channel:

\[
d\Pi^\parallel_{s,\theta_i}(\delta q_i)=0,
\]

although it can affect perpendicular presentations:

\[
d\Pi^\perp_{s,\theta_j}(\delta q_i)\neq0.
\]

For \(c\):

\[
\delta q_c\neq0
\]

is legal while

\[
\delta c=0
\]

is compulsory in the direct \(c\)-aligned view.

Therefore:

\[
\boxed{
\text{hidden motion of }c\text{ cannot be observed as varying }c;
\text{ it can only be inferred from relative perpendicular effects.}
}
\]

The same structural rule applies to each distinguished constant.

---

## 8. Other root positions still exist

The correction above does not mean that the model contains only one possible root position.

Let

\[
\mathcal S_c
\]

be the family of admissible root positions sharing the \(c\)-axial. An observer pivots views at their own local root \(s\). A complete model must nevertheless extrapolate the observable atlas at every other admissible root position:

\[
\mathfrak P_s(X)
=
\{\Pi_{s,\theta}(X):\theta\in S^1\},
\qquad
s\in\mathcal S_c.
\]

Thus two operations remain distinct:

1. **Pivoting at one local root**
   \[
   (s,\theta)\to(s,\theta').
   \]

2. **Comparing or extrapolating different root positions**
   \[
   s\to s'.
   \]

A complete model must determine both the local pivot family and the transition between root positions.

---

# PART III — MIDSAGITTAL MOVEMENT

## 9. Fixed root does not require fixed global intersection representative

The current model restricts all axes to one fixed global intersection point. That is stronger than root constancy requires.

Let \(Q\) denote a representative intersection point in the hidden owning plane and let \(\mathcal S_c\) denote the seat reduction. Root constancy requires only

\[
\mathcal S_c(Q)=c_0.
\]

It does not require

\[
Q=Q_0.
\]

The admissible hidden representatives form a fiber:

\[
\mathcal F_c
=
\mathcal S_c^{-1}(c_0).
\]

Midsagittal motion moves the hidden representative within this fiber:

\[
Q(\tau)\in\mathcal F_c,
\qquad
\mathcal S_c(Q(\tau))=c_0.
\]

Therefore the root remains constant even though its hidden representative position or relative angle changes.

---

## 10. Three descriptions of midsagittal movement

There are three candidate allocations of the same relative change.

### Type 1 — fixed observer–intersection distance

The observer traverses the axial direction while the representative intersection moves with the observer so that their relative distance remains fixed.

The observer continues to see the root as stationary, while hidden tilt or relative plane orientation changes the perpendicular presentation.

### Type 2 — fixed representative intersection

The representative intersection remains fixed while the observer moves midsagittally. The observer infers compression or expansion of the relative distance to the intersection point.

### Type 3 — fixed observer

The observer remains stationary while the representative intersection supplies the free motion. This is the reciprocal allocation of Type 2.

The conjecture is that all three become observationally equivalent when they induce the same complete perpendicular presentation:

\[
\Pi_{s,\theta}[O,Q,E]_1
=
\Pi_{s,\theta}[O,Q,E]_2
=
\Pi_{s,\theta}[O,Q,E]_3.
\]

This equivalence has not yet been proved.

The model must derive which combination of observer position \(O\), representative intersection \(Q\), and page orientation \(E\) survives the seat projection.

---

# PART IV — THM-O UNDER THE CLARIFIED ONTOLOGY

## 11. Exact proved and declared structure

For

\[
ds^2=-dt^2+\lvert d\vec x-\vec v\,dt\rvert^2
\]

with unit lapse and flat presented rods:

\[
G_{nn}=e_2(K)
\qquad\textbf{[proved identity]},
\]

\[
G_{ni}
=
-\frac12(\nabla\times\nabla\times\vec v)_i
\qquad\textbf{[proved identity]}.
\]

NORMAL-1 declares:

\[
G_{nn}=0,
\qquad
G_{ni}=0.
\]

Therefore:

\[
e_2(K)=0,
\qquad
\nabla\times\nabla\times\vec v=0
\qquad\textbf{[derived | NORMAL-1 + flat rods]}.
\]

The identity is proved. Its vacuum vanishing is declared through NORMAL-1.

---

## 12. Pinning: exponent versus coefficient

For radial inflow:

\[
e_2(K)
=
\frac1{r^2}\frac{d}{dr}(r\beta^2).
\]

NORMAL-1 gives:

\[
r\beta^2=C,
\qquad
\beta^2=\frac Cr
\qquad\textbf{[derived | NORMAL-1 + flat rods]}.
\]

This fixes the exponent and forces:

\[
c_2=0.
\]

It does not determine the source normalization:

\[
C=r_s=\frac{2GM}{c^2}.
\]

That coefficient comes from MASS-1/K-6 source normalization.

THM-K and THM-O are therefore independent derivations of the same exponent from different declared inputs:

\[
\begin{aligned}
\text{THM-K:}&\quad r^{-1}
&&\text{from the thermal route},\\
\text{THM-O:}&\quad r^{-1}
&&\text{from NORMAL-1 and flat rods}.
\end{aligned}
\]

THM-O is correctly classified as:

\[
\boxed{
\textbf{derived | NORMAL-1 + flat rods},
}
\]

not merely as a retrodiction.

---

## 13. Swirl: exponent versus amplitude

Codazzi gives:

\[
rw''+4w'=0,
\]

hence

\[
w(r)=\Omega_\infty+\frac{A}{r^3}.
\]

Removing rigid rotation at infinity leaves:

\[
w=\frac{A}{r^3}
\qquad\textbf{[derived | NORMAL-1 + flat rods + asymptotic office]}.
\]

This fixes the exponent only.

It does not derive:

\[
A=\frac{2GJ}{c^2}.
\]

That amplitude currently comes from THM-M’s boosted-source superposition, including the clock/rods split

\[
2=1+1.
\]

A sourced-Codazzi shell calculation could independently derive \(A\), but it has not been run:

\[
\boxed{
A=\frac{2GJ}{c^2}
\text{ from sourced Codazzi}
\qquad\textbf{[not run]}.
}
\]

---

## 14. Meaning of the flat-rod ceiling

The flat-rod ansatz means:

\[
h^{(c)}_{ij}=\delta_{ij}
\]

in the observer’s perpendicular account, with all gravitational information assigned to the shift \(\vec v\).

It does not mean that the fundamental planes have been physically tested for deformation.

THM-O proves that this presentation carries the Newton-compatible sectors

\[
\ell=0,1
\]

but no Newton-compatible exterior multipoles for

\[
\ell\ge2.
\]

At \(O(J^2)\), the flat presentation has:

\[
[\mathcal I]_{R^{-5}}=0
\]

and permits a half-integer tail

\[
[\mathcal I]_{R^{-9/2}}\neq0
\]

through the free \(\ell=2\) mode.

The half-integer mode is killed internally by K-6’s integer multipole structure and observationally by the static-oblate scaling test.

Therefore:

\[
\boxed{
\text{the flat observer presentation fails for higher multipoles;
the underlying planes do not deform or fail.}
}
\]

The required extension is an apparent rank-two presentation residue:

\[
h^{(c)}_{ij}
=
\delta_{ij}
+
q^{(c)}_{ij}.
\]

The open law is:

\[
q_c
\longmapsto
q^{(c)}_{ij}.
\]

This residue is conjectured to encode hidden midsagittal movement of the \(c\)-axis representative. That interpretation has not yet been derived.

---

## 15. The proposed explanation of the \(O(J^2)\) boundary

The fixed global intersection restriction is exact for the concurrent radial sector. The conjecture is that rotation requires a hidden family of nonconcurrent or skew representatives while preserving a constant seated root.

The distinction between Kerr-related congruences must remain explicit:

- the Doran/rain timelike flow;
- the spatial curl of the shift;
- the optical twist of Kerr’s principal null congruences;

are not the same object.

The principal null congruence and its ring caustic may suggest the loss of point concurrency, but they cannot be directly substituted for the rain congruence.

The proposed perturbative mechanism is:

\[
\text{hidden midsagittal displacement/tilt}=O(J),
\]

while the induced metric residue is even under spin reversal:

\[
q^{(c)}_{ij}=O(J^2).
\]

This matches the observed boundary:

\[
O(J^0),O(J^1):
\quad
\text{flat account closes},
\]

\[
O(J^2):
\quad
\text{rank-two presented residue required}.
\]

This is a structural conjecture, not yet a theorem.

---

# PART V — GR AS A SEAT-DERIVED THEORY

## 16. The conjecture

Let \(\mathcal X\) be the unseated state space and let

\[
\Pi_{s,\theta_c}:\mathcal X\to\mathcal P_c
\]

be the \(c\)-aligned observer presentation.

The conjecture is:

\[
\boxed{
\text{GR describes the derived }c\text{-aligned observer account,
not the complete unseated model.}
}
\]

The presented spacetime metric and curvature are lawful physical observations:

\[
g^{(c)}_{\mu\nu},
\qquad
R^{(c)}{}_{\mu\nu\rho\sigma},
\]

but they need not be ontic deformations of the fundamental planes.

THM-O supports the interpretation by showing that the normal Einstein projections arise as identities of the presented geometry. It does not by itself prove that GR is incomplete.

---

## 17. Factorization criterion

GR is physically incomplete relative to the model only if some model-owned observable does not factor through the complete \(c\)-seat presentation.

The required structure is:

\[
X_1\neq X_2,
\]

\[
\Pi_c(X_1)=\Pi_c(X_2),
\]

but

\[
\mathcal O(X_1)\neq\mathcal O(X_2).
\]

Equivalently:

\[
\nexists\,\widetilde{\mathcal O}
\quad\text{such that}\quad
\mathcal O=\widetilde{\mathcal O}\circ\Pi_c.
\]

If every observable factors through \(\Pi_c\), the hidden model is an ontological completion only.

---

## 18. Existing candidate pair: the deck involution

The record already contains a candidate:

\[
\tau:V\mapsto -V.
\]

The two oriented states satisfy:

\[
X_2=\tau X_1,
\]

while their Gram data agree:

\[
\pi_{\mathrm{Gram}}(X_1)
=
\pi_{\mathrm{Gram}}(X_2).
\]

The Bargmann invariant is:

\[
B
=
\frac{
1+\gamma_{12}+\gamma_{13}+\gamma_{23}+iV
}{4},
\]

so:

\[
B(\tau X)=\overline{B(X)},
\]

and generically:

\[
\arg B(\tau X)=-\arg B(X).
\]

Galois descent locates the phase in:

\[
F(\sqrt{\Delta}),
\]

while the Gram presentation lies in:

\[
F.
\]

E-8-X additionally places the phase on transitions rather than on a single presented state.

The candidate factorization test is therefore:

\[
X_1=X,
\qquad
X_2=\tau X,
\qquad
\mathcal O=\arg B.
\]

However, Gram equality alone is not yet enough. The computation must establish:

\[
\Pi_c(X)=\Pi_c(\tau X)
\]

for the complete relevant \(c\)-seat presentation, not only for its Gram coordinates.

It must also establish that \(\arg B\):

- is operationally accessible through a fixed-root pivot cycle;
- is not removed by an orientation gauge;
- is not already determined by the full \(c\)-seat presentation.

The candidate dies if:

\[
\arg B
\]

is gauge-equivalent between the deck states, is not measurable, or factors through the complete presentation despite not lying in the Gram field.

This is the first direct computation for the GR-incompleteness conjecture.

---

# PART VI — SUPERPOSITION AND MEASUREMENT: STRICT PREDICTION STATUS

## 19. Superposition is not a model condition

The following idea is a conjectured consequence only:

> Apparent superposition may arise from one fixed particle and one parasagittal shadow becoming indistinguishable under a special relative alignment of the \(c\)-view and the \(\hbar\)-plane.

It must never be inserted into the model as:

- a particle–shadow axiom;
- a collapse map;
- a two-branch state rule;
- a visibility selector;
- a measurement-record rule;
- a Born probability;
- a hand-chosen alignment gate;
- a requirement that one apparent support vanish;
- a rule guaranteeing interference.

The dependency must remain:

\[
\boxed{
\text{existing model}
\longrightarrow
\text{derived pivoted observation map}
\longrightarrow
\text{unexpected source/shadow structure, if present}
\longrightarrow
\text{quantum interpretation}.
}
\]

It must never run in reverse.

---

## 20. Exact prediction candidate

The conjecture is not that the particle occupies two absolute positions.

It is:

\[
\boxed{
\text{one fixed particle may possess a two-valued observer presentation:
the source and one parasagittal shadow.}
}
\]

The source and shadow would cast the same perpendicular projection under specific alignments. Before direct observation, an observer could infer two apparent spatial positions but could not determine which one carries the source.

Direct observation would pivot the view so that the shadow lies behind the particle or otherwise becomes indistinguishable from it. Only one apparent support would then remain.

Absolutely:

\[
X\longrightarrow X.
\]

Nothing moves, disappears, branches, or is destroyed.

The apparent collapse would be:

\[
\boxed{
\text{resolution of a non-injective observational projection,
not a change in the particle.}
}
\]

This entire mechanism remains a prediction candidate.

---

## 21. Alignment-gated, not universal

If the behaviour exists, it must occur only for specific relative configurations involving:

- the \(c\)-axis alignment;
- the provisional quantum/\(\hbar\)-plane alignment;
- the particle’s plane position;
- additional intrinsic variables derived by the model.

Observable superposition would therefore be a predicate of the complete configuration:

\[
\mathsf S=\mathsf S(X,s,\theta,\text{plane data}),
\]

not an intrinsic property of every particle:

\[
\mathsf S\neq\mathsf S(X).
\]

This could explain why superposition phenomena occur only under particular physical configurations.

But the model may not be instructed to produce this selectivity. The selector and its allowed locus must emerge from the already-owned rotor, seat, cover and transition structures.

Because the project intends to derive the \(\hbar\)-root from \(c\) and \(G\), no upstream superposition calculation may use a preconditioned \(\hbar\)-dependent angle, scale, commutator or probability rule.

---

## 22. Meaning of “global”

“Global” means global over admissible root positions sharing the \(c\)-axial:

\[
s\in\mathcal S_c.
\]

It does not mean a measurement record, state destruction or an imported instantaneous signal.

If the conjectured source–shadow presentation exists, every observer position sharing the \(c\)-axial must transport the same one-particle/one-shadow structure:

\[
T_{s's}(x_{\mathrm p}^{(s)})
=
x_{\mathrm p}^{(s')},
\]

\[
T_{s's}(x_{\mathrm{sh}}^{(s)})
=
x_{\mathrm{sh}}^{(s')}.
\]

Observers may occupy different positions in the perpendicular presented plane and assign different local box coordinates. Nevertheless:

- there is one fixed particle;
- there is one shadow;
- which apparent box contains the source is fixed absolutely;
- observation order does not move the particle;
- testing the shadow location must not detect a second source.

The model must derive this covariance if the structure appears.

---

## 23. The shadow cannot be a second source

If the conjecture emerges, the shadow must not independently carry:

\[
\text{mass},\qquad
\text{charge},\qquad
\text{energy},\qquad
\text{a second direct detection}.
\]

There must remain one conserved source.

At the same time, ordinary positional ambiguity is insufficient to reproduce quantum superposition. The model would also need to derive a transition quantity producing interference.

The existing candidate is the orientation-sensitive transition phase:

\[
\arg B,
\]

because:

\[
G(X)=G(\tau X)
\]

while

\[
B(\tau X)=\overline{B(X)}.
\]

E-8-X’s state lemma is compatible with this possibility because it puts phase on transitions rather than on individual state presentations.

Compatibility is not derivation.

---

## 24. Preregistered outcome fork

A future computation must begin only with already established model structures and derive the complete pivot-dependent observation map.

The allowed outcomes are:

### Outcome A — full structure emerges

The model independently produces:

- one source;
- one shadow presentation;
- alignment-specific spatial degeneracy;
- a transition phase;
- apparent collapse under direct alignment;
- no absolute state change.

Only then may the result be interpreted as a collapse mechanism.

### Outcome B — two-valued presentation without transition phase

The model produces classical projection ambiguity only, not quantum superposition.

### Outcome C — transition phase without spatial source/shadow degeneracy

The Bargmann phase survives as a transition observable, but the particle-shadow conjecture is dead.

### Outcome D — neither structure appears

The superposition conjecture is dead. No corrective projection rule may be inserted.

### Outcome E — the structure appears but contradicts observation

If the derived behaviour conflicts with model-owned measurements of interference, sequential measurement, erasure, Bell correlations or no-signalling, the model fails under the project’s criterion.

---

# PART VII — CURRENT COMPUTATIONAL PRIORITIES

## 25. Cheap calculations already identified

### A. Sourced Codazzi

Derive the amplitude \(A\) directly from a rotating shell’s momentum source.

Target:

\[
A=\frac{2GJ}{c^2}.
\]

Status:

\[
\textbf{[not run]}.
\]

Success would independently derive the drag amplitude and remove THM-M’s linear-superposition ancestry from that coefficient.

### B. Deck-factorization test

Determine whether:

\[
\Pi_c(X)=\Pi_c(\tau X)
\]

for the complete \(c\)-seat presentation while:

\[
\arg B(X)\neq\arg B(\tau X).
\]

This is the first direct test of whether the model contains an observable not factorizable through the GR-like \(c\)-seat account.

### C. Midsagittal projection law

Derive:

\[
(O,Q,E)
\longmapsto
(h^{(c)},v^{(c)})
\]

and test whether the three proposed allocations of midsagittal movement produce the same observer presentation.

The map must preserve:

\[
\delta\bar g_{P_i}=0.
\]

### D. Tangential law

Run candidate tangential principles first against a static nonspherical source. Any candidate that fails to produce integer exterior multipoles is dead before the \(O(J^2)\) rotating calculation.

A surviving law may then yield:

\[
Q=-\kappa_{\mathrm{model}}\frac{J^2}{M}.
\]

Only after the model derives \(\kappa_{\mathrm{model}}\) does comparison with Kerr become legal.

### E. Superposition search

Do not build a superposition suite by encoding source/shadow assumptions. First derive the unrestricted pivoted observation map. Search its output afterward for:

- non-injective spatial lifts;
- deck-even states with deck-odd transition phase;
- special alignment loci;
- one-source conservation;
- observer-position covariance.

---

# FINAL STATE

\[
\boxed{
\begin{aligned}
&\text{The fundamental planes never deform.}\\
&\text{An axis is the edge metric of its hidden owning plane.}\\
&\text{An observer occupies one local common root and pivots among axial views.}\\
&\text{Pivoting to }c,\hbar,\text{ or }G\text{ is not changing seats.}\\
&\text{A complete model must nevertheless extrapolate every admissible root position.}\\
&\text{Root constancy does not require one fixed global intersection representative.}\\
&\text{Midsagittal movement is legal inside the root’s invisible fiber.}\\
&\text{Flat rods are a restricted observer presentation, not ontic geometry.}\\
&\text{At }O(J^2)\text{ the flat presentation lacks rank-two multipole capacity.}\\
&\text{The tangential }G_{ij}[h,v]\text{ are computable; their model-owned law is missing.}\\
&\text{GR-as-a-seat-derived-theory remains a conjecture with a concrete deck-phase test.}\\
&\text{The particle-shadow account of superposition is prediction-only.}\\
&\text{It must fall out for free or be discarded.}
\end{aligned}
}
\]

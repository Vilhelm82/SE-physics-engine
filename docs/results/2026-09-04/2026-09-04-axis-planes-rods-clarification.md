## Required corrections

### 1. Separate ontology, presentation, and proved results

The document must distinguish:

\[
\bar g_{P_i}\ \text{fixed}
\qquad\textbf{[model ontology/axiom]},
\]

from

\[
h^{(c)}_{ij}[X]\ \text{variable}
\qquad\textbf{[observer presentation]},
\]

and from

\[
q_c\longmapsto h^{(c)}_{ij}
\qquad\textbf{[conjectured map; not yet derived]}.
\]

THM-O proves a limitation of the flat presentation. It does not yet prove that hidden motion of the \(c\)-page generates the required \(q^{(c)}_{ij}\).

### 2. Non-invertibility alone does not generate curvature

This implication is not established:

\[
\mathcal S_c\text{ non-invertible}
\quad\Longrightarrow\quad
R[h^{(c)}]\neq0.
\]

A non-invertible projection permits information loss, but the map producing the presented metric still has to be derived:

\[
h^{(c)}_{ij}
=
\mathcal F_{ij}[q_c,\partial q_c,\ldots].
\]

The open theorem is to find \(\mathcal F\) and prove that its curvature reproduces the required observer readings.

### 3. Correct the coframe criterion

This statement is false:

\[
dE^a\neq0
\quad\Longrightarrow\quad
R[h]\neq0.
\]

A nonholonomic coframe can describe flat geometry. For example, polar orthonormal coframes have \(dE^a\neq0\) while the metric remains flat.

The correct calculation is:

\[
dE^a+\omega^a{}_b\wedge E^b=0
\]

for the torsion-free connection, followed by

\[
\mathcal R^a{}_b
=
d\omega^a{}_b
+
\omega^a{}_c\wedge\omega^c{}_b.
\]

Apparent curvature requires

\[
\mathcal R^a{}_b\neq0,
\]

not merely \(dE^a\neq0\).

The sufficient flatness statement remains valid:

\[
E^a=dy^a
\quad\Longrightarrow\quad
R[h]=0
\]

locally, but its converse is false.

### 4. Tighten the Gauss–Codazzi scope

Do not say Gauss and Codazzi are everything the rods can say about the ambient geometry. What THM-O proves is that they supply the complete **normal constraint projections**:

\[
G_{nn}=e_2(K)
\]

on flat rods, and

\[
G_{ni}
=
-\frac12(\nabla\times\nabla\times\vec v)_i
\]

with your fixed conventions.

Their vanishing is conditional:

\[
G_{nn}=G_{ni}=0
\qquad\textbf{[declared: NORMAL-1]}.
\]

The tangential symmetric sector remains unowned by the model.

### 5. State the GR conjecture as a testable factorization claim

It is not enough that GR does not mention the proposed hidden state. To establish physical incompleteness, the model needs states \(X_1,X_2\) such that

\[
X_1\neq X_2,
\qquad
\Pi_c(X_1)=\Pi_c(X_2),
\]

and a model-owned observable satisfying

\[
\mathcal O(X_1)\neq\mathcal O(X_2).
\]

Equivalently, some observable must fail to factor through the \(c\)-seat presentation:

\[
\nexists\,\widetilde{\mathcal O}
\quad\text{such that}\quad
\mathcal O=\widetilde{\mathcal O}\circ\Pi_c.
\]

Without that, the hidden state is an ontological completion but not an empirically distinguishable completion.

## Clean replacement to give Claude

---

### Axis, plane, seat, and flat-rod clarification

#### Model ontology

Each fundamental plane \(P_i\) has a fixed intrinsic metric:

\[
\delta\bar g_{P_i}=0.
\]

The planes never deform.

An axis is not an intrinsic line contained in the plane from which it is observed. It is the edge-on metric trace of an adjacent plane. In particular, the \(c\)-axis is the edge presentation of \(P_c\).

Seating on \(c\) locks this edge as the observer’s root. Direct variation along its owning plane is annihilated by the seat:

\[
d\mathcal S_c(\delta q_c)=0,
\]

although the hidden configuration may vary:

\[
\delta q_c\neq0.
\]

If \(q_c\) includes orientation as well as position, its proper domain is a configuration or frame bundle over \(P_c\), not merely \(P_c\).

#### Observer presentation

The \(c\)-seat has access only to a reduced presentation:

\[
\Pi_c:X\longmapsto
\left(h^{(c)}_{ij},v^{(c)i},\ldots\right).
\]

Variation invisible in the \(c\)-channel is inferred from relative changes in the perpendicular presentation.

Consequently,

\[
R[h^{(c)}]\neq0
\]

means that the \(c\)-seat’s spatial account is curved. It does not mean that any fundamental plane has deformed:

\[
R[h^{(c)}]\neq0
\quad\not\Longrightarrow\quad
\delta\bar g_{P_i}\neq0.
\]

#### Meaning of flat rods in THM-O

The flat-rod ansatz is

\[
h^{(c)}_{ij}=\delta_{ij},
\]

with all observable gravitational content assigned to the shift:

\[
ds^2=-dt^2+\lvert d\vec x-\vec v\,dt\rvert^2.
\]

This is an observer-side presentation restriction, not an ontological claim about deformation.

THM-O proves that this presentation closes for the monopole and linear current dipole:

\[
O(J^0),\qquad O(J^1),
\]

but has no Newton-compatible multipole channel for

\[
\ell\ge2.
\]

At \(O(J^2)\), therefore, the \(c\)-seat requires a richer presented metric:

\[
h^{(c)}_{ij}
=
\delta_{ij}+q^{(c)}_{ij}.
\]

The tensor \(q^{(c)}_{ij}\) is an observer-side presentation residue. It is not physical deformation of an underlying plane.

#### Exact open debt

The model has not yet derived the projection law

\[
q_c
\longmapsto
q^{(c)}_{ij}.
\]

The required theorem is to construct

\[
q^{(c)}_{ij}
=
\mathcal F_{ij}[q_c,\partial q_c,\ldots]
\]

from the hidden movement of the seated root and then determine its presented curvature.

A coframe representation may be used:

\[
h^{(c)}_{ij}
=
\delta_{ab}E^a{}_iE^b{}_j,
\]

but nonholonomicity alone is not curvature. One must compute

\[
dE^a+\omega^a{}_b\wedge E^b=0,
\qquad
\mathcal R^a{}_b
=
d\omega^a{}_b+\omega^a{}_c\wedge\omega^c{}_b.
\]

#### GR conjecture

The conjecture is:

\[
\boxed{
\text{GR governs the }c\text{-seat’s derived presentation,
not the complete unseated state.}
}
\]

Formally, GR acts on the image of

\[
\Pi_c:X\rightarrow\mathcal G_c,
\]

where \(\mathcal G_c\) is the space of \(c\)-seat metric presentations.

THM-O supports this interpretation by showing that the normal Einstein projections arise as identities of the flat observer presentation. It does not yet prove that GR is incomplete.

Physical incompleteness requires a model-owned observable that does not factor through \(\Pi_c\):

\[
\Pi_c(X_1)=\Pi_c(X_2),
\qquad
\mathcal O(X_1)\neq\mathcal O(X_2).
\]

The candidate places to establish this are:

1. cross-seat observables;
2. an independently derived tangential closure;
3. a non-GR integer quadrupole coefficient;
4. derivation of the \(\hbar\)-seat from the hidden \(c\)- and \(G\)-root structure.

---


# The locked inertia tensor as a rank-one detector

Companion to `GRAM_SUBMERSION_CURVATURE.md`. Same object, stated as a boundary-classification result.

---

## 1. Objects

Three unit vectors $a_1,a_2,a_3\in S^2$, $A=[a_1\,a_2\,a_3]$.

$$G=A^{\mathsf T}A \quad(\text{intrinsic, } G_{ij}=\gamma_{ij}),\qquad M=AA^{\mathsf T}=\sum_i a_ia_i^{\mathsf T}\quad(\text{extrinsic})$$

$$\mathbb I \;=\; 3\,\mathrm{Id}-M \qquad\text{(locked inertia tensor of the simultaneous }SO(3)\text{ action)}$$

$\mathbb I$ arises as the metric on the $SO(3)$ orbit: $\sum_i\lVert\omega\times a_i\rVert^2=\omega^{\mathsf T}\mathbb I\,\omega$.

$G$ and $M$ are Gram/dual-Gram: same nonzero spectrum, $\operatorname{rank}G=\operatorname{rank}M$, $\det M=\det G=\Delta$, and $\operatorname{tr}M=3$ (unit axes).

$$\Delta=\det G=1-\gamma_{12}^2-\gamma_{13}^2-\gamma_{23}^2+2\gamma_{12}\gamma_{13}\gamma_{23}$$

---

## 2. Theorem (rank-one detection) [proved]

$$\boxed{\;\det\mathbb I\;\ge\;0\ \text{ on the whole physical region},\qquad \det\mathbb I=0\iff \operatorname{rank}G=1\;}$$

Hence $\{\det\mathbb I=0\}$ is exactly the four rank-one nodes $\gamma=(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)$, a proper subset of $\{\Delta=0\}$:

$$\{\det\mathbb I=0\}\;\subsetneq\;\{\Delta=0\}$$

**Proof.** $M$ is PSD with $\operatorname{tr}M=3$, so its eigenvalues satisfy $\mu_i\ge0$, $\sum\mu_i=3$, giving $\mu_i\le3$ with $\det\mathbb I=\prod_i(3-\mu_i)\ge0$. Equality holds iff some $\mu_i=3$, which forces the other two to vanish, i.e. $\operatorname{rank}M=1$, i.e. $\operatorname{rank}G=1$: all three axes parallel. $\square$

| stratum | $\operatorname{rank}G$ | $\Delta$ | $\operatorname{eig}(\mathbb I)$ | $\det\mathbb I$ |
|---|---|---|---|---|
| generic | 3 | $>0$ | $3-\mu_k>0$ | $>0$ |
| coplanar | 2 | $0$ | $3-\mu_1,\ 3-\mu_2,\ 3$ | $>0$ |
| node | 1 | $0$ | $0,3,3$ | $0$ |

$\mathbb I$ is blind to the coplanar stratum and degenerates only at genuine rank collapse.

---

## 3. $\det\mathbb I$ in Gram coordinates [proved]

$\mathbb I$ looks extrinsic but its determinant is a polynomial in the intrinsic data:

$$\boxed{\;\det\mathbb I \;=\; 3\,e_2(G)-\Delta \;=\; 8-2\big(\gamma_{12}^2+\gamma_{13}^2+\gamma_{23}^2\big)-2\,\gamma_{12}\gamma_{13}\gamma_{23}\;}$$

equivalently, with $e_3=\gamma_{12}\gamma_{13}\gamma_{23}$,

$$\det\mathbb I \;=\; 6+2\Delta-6e_3,\qquad\text{so on }\{\Delta=0\}:\quad \det\mathbb I=6\,(1-e_3).$$

From $\det\mathbb I=\prod(3-\mu_i)=27-9\!\sum\!\mu_i+3e_2-\Delta$ with $\sum\mu_i=3$, $e_2(M)=e_2(G)=3-\sum\gamma_{ij}^2$.

**No frame reconstruction is required.** The test is evaluated directly on measured overlaps.

**Values.** identity $\gamma=(0,0,0)$: $\det\mathbb I=8$. Trine $\gamma_{ij}=-\tfrac12$: $\det\mathbb I=\tfrac{27}{4}$. Nodes: $0$. Near-collinear coplanar $(0°,170°,190°)$: $0.5319$ — small but strictly positive.

**Verification.** Against explicit $3\times3$ matrices at three exact rational frames; all four nodes give exactly $0$; five exactly-coplanar rational points give $\Delta=0$ with $\det\mathbb I\in\{\tfrac{114054}{21025},\tfrac{1733664}{351125},\tfrac{34438176}{12880921},\tfrac{48534}{7225},\tfrac{6603438}{2628125}\}$, all $>0$. Two coincident axes with the third free ($\gamma=(1,\tfrac13,\tfrac13)$) gives $\tfrac{16}{3}>0$: **collinearity of a pair is not detected**, only total collapse.

---

## 4. Relation to $\Delta$

$\Delta$ does distinguish the strata, but through the **local vanishing order** $m$ ($m=1$ coplanar, $m=2$ node), which requires an approach analysis. $\det\mathbb I$ resolves the same distinction by a **pointwise algebraic evaluation**:

$$m=1\iff\det\mathbb I\neq0\ \text{ on }\{\Delta=0\},\qquad m=2\iff\det\mathbb I=0.$$

---

## 5. Curvature consequence [proved]

For the Riemannian submersion $\pi:(S^2)^3\to$ elliptope, $A\mapsto G$, the O'Neill curvature is

$$F(X,Y)=2\,\mathbb I^{-1}\Big(\sum_i X_i\times Y_i\Big),\qquad \lVert F\rVert\sim C(b)\,\Delta^{-m/2}=C(b)\,|V|^{-m},\quad V=\det A,\ V^2=\Delta.$$

| stratum | $m$ | exponent (exact-rational, target) | scaling |
|---|---|---|---|
| coplanar | 1 | $-1.00031\to-1$ in $\lVert F\rVert^2$ vs $\Delta$ | $\lVert F\rVert\sim|V|^{-1}$ |
| node | 2 | $-2.00135\to-2$ | $\lVert F\rVert\sim|V|^{-2}$ |

### Why the exponent doubles — kernel geometry

Let $M_0$ be the horizontal-lift system at the degenerate frame; every lift's divergent part lies in $\ker M_0$.

- **Coplanar ($m=1$).** $\dim\ker M_0=1$, purely out-of-plane: $N=(z_1e_z,z_2e_z,z_3e_z)$. Singular parts of any two lifts are $\alpha N,\beta N$, so per site the vectors are **parallel** and $\sum_iX_{-1,i}\times Y_{-1,i}=\alpha\beta\sum_iz_i^2(e_z\times e_z)=0$ identically. The $\varepsilon^{-2}$ term **cancels**; $\mathbb I$ is regular; $\lVert F\rVert=O(\Delta^{-1/2})$.

- **Node ($m=2$).** $\dim\ker M_0=4$, lying entirely **perpendicular** to the dying eigendirection $e_3$ of $\mathbb I$ (all kernel–$e_3$ inner products exactly zero). Perpendicularity means the cross-product pairing does **not** cancel, and its image is **parallel** to $e_3$ — precisely $\mathbb I$'s degenerating direction, where $\mathbb I_{33}=\tfrac{12t^2}{(1+t^2)^2}=3\sin^2\varepsilon$. Numerator $\varepsilon^{-2}$ along $e_3$ meets $\mathbb I^{-1}\sim(3\varepsilon^2)^{-1}$ on $e_3$: $\lVert F\rVert=O(\varepsilon^{-4})=O(\Delta^{-1})$ since $\Delta\sim\varepsilon^4$.

> **Precision point.** It is the *image of the cross-product pairing*, not the kernel, that aligns with the dying eigenvector. The kernel is perpendicular to it. Stating the kernel as "aligned" describes the $m=1$ case, which cancels and does not double.

### $\mathbb I$ also fixes the constant

On the smooth stratum, with $\mathbb I_\parallel$ the in-plane block,

$$C^2(b)=\frac{P(e_1,e_2)}{\big(\det\mathbb I_\parallel\big)^4},\qquad \det\mathbb I_\parallel=2(1-e_3),\qquad \det\mathbb I=3\det\mathbb I_\parallel\ \text{ on }\{\Delta=0\}$$

(verified exactly at three rational coplanar points). $\mathbb I$ therefore controls the exponent (via $\ker$), the doubling (via degeneracy), and the constant (via $\det$).

---

## 6. Scope and lineage

**Not new: the operator.** $\mathbb I$ is the locked inertia tensor of geometric mechanics (Marsden; Littlejohn–Reinsch; molecular gauge theory), standard for the $SO(3)$ action on shape space. The mechanical connection $\mathcal A=\mathbb I^{-1}J$ and O'Neill's formula are likewise classical.

**Candidate-new:** the characterization $\det\mathbb I=0\iff$ rank-one node with the Gram-coordinate closed form; the curvature exponent law $\lVert F\rVert\sim C(b)|V|^{-m}$ with the kernel-geometry mechanism; $C(b)$'s closed form. Novelty **unverified against the literature**.

**Terminology.** $G$ is the Gram matrix of three unit Bloch vectors — a $3\times3$ correlation matrix on the elliptope — not a density matrix. For pure qubits $\operatorname{rank}\rho=1$ identically; the relevant rank here is the Gram rank of the frame. Under the spinor lift, $V=\operatorname{Im}(4B)$ for the Bargmann invariant $B=\tfrac{1+\sum\gamma_{ij}+iV}{4}$, so $\{V=0\}=\{$coplanar$\}$ is the zero-imaginarity locus.

**Unverified claim.** Whether standard tomography and resource theories detect boundary states via $\operatorname{rank}\rho$ or $\Delta$ is a literature question not checked here; it is the load-bearing novelty premise.

---

## 7. Functions

```python
import sympy as sp

def det_inertia_from_gram(g12, g13, g23):
    """det I directly from the Gram entries. Zero iff rank-one node."""
    return 8 - 2*(g12**2 + g13**2 + g23**2) - 2*g12*g13*g23

def delta_from_gram(g12, g13, g23):
    return 1 - g12**2 - g13**2 - g23**2 + 2*g12*g13*g23

def inertia(A):
    """A is 3x3 with unit columns."""
    return 3*sp.eye(3) - A*A.T

def classify(g12, g13, g23, tol=0):
    D  = delta_from_gram(g12, g13, g23)
    dI = det_inertia_from_gram(g12, g13, g23)
    if D  > tol:  return "interior (rank 3)"
    if dI > tol:  return "coplanar boundary, m=1, ||F|| ~ |V|^-1"
    return "rank-one node, m=2, ||F|| ~ |V|^-2"
```

Checks: `verify_abstract.py` (8/8 with exact rationals).

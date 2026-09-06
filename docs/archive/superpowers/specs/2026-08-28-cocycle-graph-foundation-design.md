# Abstract Cocycle and Graph-Descent Foundation

**Status:** approved in chat on 2026-08-28; implementation awaits review of this written design.

## Purpose

Rebuild the mathematical spine of *The Seated Root* in the order

1. abstract seat-change cocycle;
2. graph and cycle-descent theorem;
3. the \(c\), \(\hbar\), and \(G\) seats as representations.

The existing three-seat results become corollaries and regression targets, not premises of the abstraction. The derivation must preserve the repository's separation between bare geometry, representation-specific structure, declared physical identifications, and comparison-stage retrodictions.

## Source of truth and working method

- **/home/williaml/seated-root/PAPER-seated-root-v0.1.md** remains the canonical paper source; the HTML remains generated output.
- Existing theorem scripts are evidence and regression tests. They do not substitute for an arbitrary-case proof.
- The new foundation is first written and checked as a standalone mathematical note. It is integrated into the paper only after its assumptions, theorem boundaries, and specialisations pass review and exact checks.
- No freeze, hash, or preregistration ceremony is added.

## Scope boundaries

This work will not:

- derive the numerical values or physical meanings of \(c\), \(\hbar\), or \(G\) from the abstract mathematics;
- discharge LBL-1, KIN-1, KIN-2, DEBT-2b, or the unassigned cross-seat measurement protocol;
- claim that an arbitrary groupoid action supplies a canonical measure without carrier geometry;
- identify an arbitrary graph vertex with a fundamental seat without a separate representation map;
- call a result novel merely because it has not yet been located in a targeted search.

## I. Abstract seat-change calculus

### Primitive data

A **carrier representation** of the seat-change groupoid consists of:

1. a groupoid \(\mathcal S\) of seats and reversible presentation changes;
2. for each object \(s\), a carrier \(X_s\) with geometry-induced finite probability measure \(\mu_s\);
3. for each arrow \(g:s\to t\), an invertible carrier map \(T_g:X_s\to X_t\);
4. the pullback Jacobian cocycle \(J_g:X_s\to\mathbb R_{>0}\), fixed by
   \[
   T_g^*\mu_t=J_g\,\mu_s;
   \]
5. when it exists, a positive presentation cocycle \(f_g\) and representation weight \(Q>0\) satisfying
   \[
   J_g=f_g^{-Q}.
   \]

For composable \(g:s\to t\) and \(h:t\to u\), this convention gives
\[
J_{h\circ g}(x)=J_g(x)J_h(T_gx),\qquad
f_{h\circ g}(x)=f_g(x)f_h(T_gx).
\]
Inversion gives
\[
J_{g^{-1}}(T_gx)=J_g(x)^{-1},\qquad
f_{g^{-1}}(T_gx)=f_g(x)^{-1}.
\]
These laws are derived from the chain rule, not separately assumed.

The phrase **the measure is not a free input** will mean precisely this: after the carrier geometry and its action have been specified, its induced measure and Jacobian are fixed. Carrier geometry itself remains representation data.

### Theorem ladder

**C0 -- Jacobian cocycle.** Prove the identity, composition, and inverse laws for \(J_g\), and hence for \(f_g\) wherever the positive \(Q\)-th root exists.

**C1 -- RN-root theorem at fixed weight.** Once the carrier representation fixes its density weight \(Q\), the unique positive presentation compatible with \(J_g=f_g^{-Q}\) is \(f_g=J_g^{-1/Q}\). This fixes the geometric candidate and removes a free measure choice. It does not derive \(Q\) from an arbitrary groupoid, nor does it prove that a particular detector couples to \(f_g\); those remain representation and labelled-coupling statements respectively.

**C2 -- Seat bridge and reflection.** For two presentations of the same carrier with
\[
d\mu_1=f^{-Q}d\mu_0,
\]
define the unnormalised geometric bridge
\[
d\nu_m=f^{-Qm}d\mu_0,
\qquad 0\le m\le1,
\]
and Mellin readout
\[
I_{01}(p,m)=\int f^p\,d\nu_m.
\]
Exchanging the seats sends \((f,m,p)\mapsto(f^{-1},1-m,-p)\) and leaves the bridge unchanged:
\[
\nu^{01}_m=\nu^{10}_{1-m},\qquad
I_{01}(p,m)=I_{10}(-p,1-m).
\]
The repository's \(I(p,m)=I(-p,1-m)\) follows when the two oriented descriptions are identified by seat reversal. At \(m=1/2\), \(\nu_{1/2}=\sqrt{d\mu_0d\mu_1}\) is the Hellinger/Bhattacharyya midpoint.

**C3 -- Blind-mass theorem.** Derive
\[
\int f^Q\,d\mu_1=\mu_0(X)=1.
\]
Prove the universal uniqueness statement by variation over positive \(f\) subject to \(\int f^{-Q}d\mu_0=1\): the only moments fixed for every admissible \(f\) are the normalisation moment \(p=0\) and the nontrivial blind moment \(p=Q\). Any additional blind exponent in a restricted family must be labelled accidental or representation-specific.

**C4 -- Finite invariant-measure support.** For a hyperbolic carrier map whose recurrent set equals its fixed-point set, use recurrence to prove that every finite invariant Borel measure is supported on those fixed points. Infinite or merely sigma-finite measures are explicitly excluded.

**C5 -- Pole-weight theorem.** Pivot invariance alone fixes the two-pole support but leaves arbitrary probability weights \(a\delta_++(1-a)\delta_-\). Invariance under the additional pole-exchange involution forces \(a=1/2\). Only then is the equal two-pole geometric mean canonical. This restores the qualifier already present in **thm_g.py** but compressed out of the paper prose.

## II. Graph and cycle-descent theorem

### Gauge-cohomology layer

Let \(\Gamma=(V,E)\) be a finite connected graph with one reference orientation chosen for each edge, on the open locus where every used overlap is nonzero. Associate an edge amplitude
\[
z_{ij}=\langle\psi_i\mid\psi_j\rangle\in\mathbb C^\times.
\]
Vertex rephasings act by
\[
z_{ij}\longmapsto h_i^{-1}z_{ij}h_j.
\]
Separate it canonically into magnitude and phase,
\[
z_{ij}=\rho_{ij}u_{ij},\qquad
\rho_{ij}=|z_{ij}|>0,\qquad u_{ij}\in U(1).
\]
The reverse edge satisfies \(\rho_{ji}=\rho_{ij}\) and \(u_{ji}=u_{ij}^{-1}\). Thus the normalized phases \(u\), not the full overlaps \(z\), form a multiplicative \(U(1)\) 1-cochain; vertex gauges are its coboundaries. An oriented closed walk \(C\) has phase holonomy
\[
\operatorname{Hol}_u(C)=\prod_{e\in C}u_e^{\epsilon_e},
\]
while its full Bargmann invariant is
\[
B_C=\left(\prod_{e\in C}\rho_e\right)\operatorname{Hol}_u(C).
\]

Using a spanning tree, prove the exact sequence
\[
1\longrightarrow C^0(\Gamma,U(1))/U(1)
\xrightarrow{\delta}C^1(\Gamma,U(1))
\xrightarrow{\mathrm{Hol}}\operatorname{Hom}(H_1(\Gamma;\mathbb Z),U(1))
\longrightarrow1.
\]
Thus a connected graph has \(\beta_1=|E|-|V|+1\) independent phase holonomies, alongside \(|E|\) gauge-invariant positive edge magnitudes subject to the state-space Gram constraints. Zero overlaps are recovered only as boundary strata after the nonzero theory is proved.

### Rank-three projector layer

For unit vectors \(a_i\in S^2\), let
\[
P_i=\frac{1+a_i\cdot\sigma}{2}=|\psi_i\rangle\langle\psi_i|.
\]
For every oriented cycle \(C=(i_1\ldots i_m)\), prove
\[
B_C=\operatorname{tr}(P_{i_1}\cdots P_{i_m}).
\]

Choose one noncoplanar reference triple with Gram determinant \(\Delta\ne0\) and orientation coordinate \(V^2=\Delta\). Let \(F\) be the function field generated by the full rank-three Gram data. Coordinates relative to the reference triple show that every scalar triple product is \(F\)-proportional to \(V\). Pauli reduction or orthogonal invariant theory then gives the **graph-wide descent form**
\[
B_C=A_C+iV D_C,
\qquad A_C,D_C\in F,
\]
for every closed cycle, independent of its length.

The proof will include the exact Pauli recurrence
\[
\alpha_{k+1}=\alpha_k+b_k\cdot a_{k+1},
\qquad
b_{k+1}=b_k+\alpha_k a_{k+1}+i\,b_k\times a_{k+1},
\]
for
\[
\prod_{j=1}^k(1+a_j\cdot\sigma)=\alpha_k1+b_k\cdot\sigma,
\]
with \(B_C=2^{1-m}\alpha_m\).

Taking the modulus yields the general cycle identity
\[
A_C^2+\Delta D_C^2
=2^{-m}\prod_{e\in C}(1+\gamma_e).
\]
For \(m=3\), this must reduce exactly to
\[
B_{123}=\frac{1+\gamma_{12}+\gamma_{23}+\gamma_{31}+iV}{4}
\]
and to the repository's cubic radicand identity.

### Kummer field statement

For edge radicands \(r_e=(1+\gamma_e)/2\), let \(R_\Gamma\subset F^*/F^{*2}\) be their generated square-class subgroup and let \(r_\Gamma=\dim_{\mathbb F_2}R_\Gamma\). The leg-magnitude field has degree
\[
[F(\sqrt{r_e}:e\in E):F]=2^{r_\Gamma}.
\]
No blanket independence claim is made without an odd-valuation proof for the chosen graph.

The full local-lift field may contain \(i\), \(V\), and all edge radicals. By contrast, the field generated by all graph-cycle holonomies is contained in \(F(i,V)\), so over \(F(i)\) its degree is at most two. It equals the orientation field whenever at least one cycle has \(D_C\ne0\); otherwise it descends to \(F(i)\). Graph automorphisms may be added as a semidirect action on the square-class module only after their action and splitting are proved. The existing \((C_2\wr S_3)\times C_2\) result is recovered specifically for the three-edge seat triangle.

## III. The three representations

Only after Sections I and II are complete will the physical labels appear.

### \(c\)-representation

Use
\[
T_\lambda(u)=\frac{u+\tanh\lambda}{1+\tanh\lambda\,u},
\qquad
s_\lambda(u)=\cosh\lambda+\sinh\lambda\,u.
\]
Prove the cocycle law and
\[
T_\lambda'(u)=s_\lambda(u)^{-2},
\]
so \(Q=2\). Recover the unpivoted, pivoted, and equal-two-pole power-mean rows, including seat reciprocity. The fixed-point theorem supplies support at \(u=\pm1\); pole exchange supplies equal weights.

### \(\hbar\)-representation

Use the symplectic squeeze \(S_\lambda=\operatorname{diag}(e^\lambda,e^{-\lambda})\) on a projective line. Derive
\[
q_\lambda(\phi)=\cosh2\lambda+\sinh2\lambda\cos2\phi,
\qquad
\frac{d\phi'}{d\phi}=q_\lambda(\phi)^{-1},
\]
so \(Q=1\). Recover the \(c\)-representation through \((\Theta,\Lambda)=(2\phi,2\lambda)\), then recover the determinant floor, orbit separation, mirror, blind moment, and Legendre index shift. The value assigned to the protected area remains comparison-stage physics.

### \(G\)-representation

Pull back the same rapidity cocycle along the declared pinning map
\[
\tanh\lambda(r)=\sqrt{r_s/r}.
\]
With \(N=\operatorname{sech}\lambda\), decompose the pole account as
\[
(\log(1+v),\log(1-v))
=\log N\,(1,1)+\lambda(1,-1).
\]
This separates the common lapse scale from the reciprocal cocycle residue. Pinning sends \(\lambda=\infty\) to the finite boundary \(r=r_s\), yielding the bound/escaping-pole anatomy and the free/floored/pinned trichotomy. KIN-2 remains declared; Tolman, Unruh, and Hawking names remain comparison-stage recoveries.

## Artifacts and integration order

1. **FOUNDATION-cocycle-graph.md** -- complete definitions and human proofs.
2. **thm_cocycle_graph.py** -- exact symbolic receipts and adversarial finite cases.
3. **run_cocycle_graph.log** -- captured verification output.
4. **PAPER-seated-root-v0.1.md** -- reordered theorem-first exposition with lineage and status labels.
5. **PAPER-seated-root.html** -- regenerated from the canonical Markdown.

The standalone note will distinguish established Radon--Nikodym, Hellinger, graph-cohomology, and Bargmann machinery from the graph-wide rank-three Kummer descent statement whose novelty remains to be audited.

## Verification

- Prove arbitrary-length statements in prose; do not infer them from testing finitely many cycles.
- Check composition, inversion, bridge reflection, blind mass, and both representation Jacobians symbolically.
- Check the Pauli recurrence and cycle identity exactly for symbolic triangles and quadrilaterals, then against independent high-precision spinor calculations for longer randomly generated cycles.
- Recover every current triangle identity and all nine cells of the seat-correct temperature table.
- Run every pre-existing theorem script unchanged and require all existing checks to remain green.
- Rebuild the self-contained HTML and inspect its structure for missing mathematics, figures, or references.
- Run repository diff and whitespace checks before any completion claim.

## Acceptance conditions

The work is complete only when:

1. no \(c\), \(\hbar\), or \(G\) assumption occurs in the abstract or graph proofs;
2. every representation-specific input is named at the point of specialisation;
3. finite invariant support and equal pole weighting are separate theorems;
4. the graph theorem states square-class rank conditionally rather than assuming generic independence;
5. the triangle result is obtained as a corollary of the arbitrary-cycle theorem;
6. Q-RN is split honestly into a mathematical RN-root theorem and an unresolved physical detector-coupling question;
7. all new and existing checks pass and the canonical HTML rebuilds successfully.

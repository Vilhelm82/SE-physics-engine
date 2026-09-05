# Cella sources for the horizon and load continuation

5 September 2026.

Target: derive the seat/load response and cut-area transport using the Cella constructions below.

## 1. Coupling form and native signature

**Existing Cella result.** The [Gauge Channel Transport Law](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/dbp_role_channel_and_orbit_geometry/Gauge_Channel_Transport_Law.md:191>) and [Canonical Invariant Reduction Theorem](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/dbp_role_channel_and_orbit_geometry/Canonical_Invariant_Reduction_Theorem.md:214>) give

\[
\nu=(g_1H_{23},g_2H_{13},g_3H_{12}),\qquad
\Delta_c=\nu^T B\nu,\qquad B=2I-\mathbf1\mathbf1^T.
\tag{1}
\]

Here \(g\) is the gradient of an implicit relation and \(H_c\) its off-diagonal Hessian in a declared role basis. Expanding the bordered determinant gives

\[
\det\begin{pmatrix}0&g^T\\g&H_c\end{pmatrix}
=\nu_1^2+\nu_2^2+\nu_3^2
-2\nu_1\nu_2-2\nu_1\nu_3-2\nu_2\nu_3,
\]

which proves (1).

Let \(e_0=(1,1,1)/\sqrt3\), and choose orthonormal \(e_1,e_2\) perpendicular to \(e_0\). Then \(Be_0=-e_0\) and \(Be_i=2e_i\). Consequently

\[
T=(e_0,e_1/\sqrt2,e_2/\sqrt2),\qquad
T^TBT=\operatorname{diag}(-1,1,1).
\tag{2}
\]

For a defining-function gauge with value one at the surface point and logarithmic gradient \(a\), Cella gives

\[
\widetilde H=H+ga^T+ag^T,\qquad
\delta\nu=
\begin{pmatrix}
g_1(g_2a_3+g_3a_2)\\
g_2(g_1a_3+g_3a_1)\\
g_3(g_1a_2+g_2a_1)
\end{pmatrix}.
\tag{3}
\]

For \(g_1g_2g_3\ne0\), the gauge map has determinant \(2(g_1g_2g_3)^2\), so it can move the edge vector arbitrarily. **Target:** construct native seat representatives and a presentation rule or covariant readout preserving their causal classification.

Path transport requires a rule relating successive relation jets: EXT-001/005.

## 2. Signed diagonal curvature

**Source:** [Local Curvature Calculus](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/local_curvature_and_black_hole_metrics/LOCAL_CURVATURE_CALCULUS_COMPLETE_v1.0.md:121>), directional channels and divisor valuations for positive diagonal metrics.

**Source:** [Weighted-jet theorem](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/local_curvature_and_black_hole_metrics/LEAD7_VARIABLE_TRANSVERSE_WEIGHTED_JET_THEOREM_v1.0.md:1>), arbitrary transverse coefficients with nonzero leading terms of either sign.

**Extension derived here.** On a patch where every diagonal component is nonzero, write

\[
g=\sum_i\varepsilon_i e^{2u_i}(dt^i)^2,\qquad
\varepsilon_i\in\{-1,1\}\ \text{constant}.
\]

With exactly Cella's channel polynomial

\[
\mathcal Q_j=
\sum_{i\ne j}\left[
\partial_j^2u_i+(\partial_ju_i)^2
-(\partial_ju_j)(\partial_ju_i)\right]
+\sum_{\substack{r<s\\r,s\ne j}}
(\partial_ju_r)(\partial_ju_s),
\]

the scalar curvature is

\[
\boxed{R[g]=-2\sum_j\varepsilon_j e^{-2u_j}\mathcal Q_j
=-2\sum_j\frac{\mathcal Q_j}{g_{jj}}.}
\tag{4}
\]

**Proof.** Put \(H_i=e^{u_i}\) and \(\beta_{ij}=\partial_iH_j/H_i\). The diagonal Christoffel formulas are

\[
\Gamma^i_{ii}=\partial_i u_i,\quad
\Gamma^i_{ij}=\partial_j u_i,\quad
\Gamma^i_{jj}=-\varepsilon_i\varepsilon_j
e^{2(u_j-u_i)}\partial_i u_j\quad(i\ne j).
\]

Substitution into the curvature definition gives the sectional curvature of the nondegenerate coordinate plane:

\[
K_{ij}=-\frac1{H_iH_j}
\left(\varepsilon_i\partial_i\beta_{ij}
+\varepsilon_j\partial_j\beta_{ji}
+\sum_{\ell\ne i,j}\varepsilon_\ell\beta_{\ell i}\beta_{\ell j}\right).
\]

Contracting \(R=2\sum_{i<j}K_{ij}\) and grouping by differentiation direction gives (4). This keeps the convention that the unit round two-sphere has scalar curvature \(+2\).

For \(g_{ii}=\varepsilon_i h_i z^{P_i}\), \(h_i>0\), the divisor decomposition in Cella's Theorem 3.1 becomes

\[
R=\sum_{a\le k}\varepsilon_a z^{-P_a-2e_a}F_a
+\sum_{\mu>k}\varepsilon_\mu z^{-P_\mu}F_\mu,
\tag{5}
\]

The coefficient germs \(F_i\), candidate support and jet dependence are unchanged. At coincident exponents, sum coefficients with the factors \(\varepsilon_i\); positivity-dependent corollaries require separate hypotheses.

**Application already calculated.** For the completed load witness,

\[
g=-dT^2+dr^2+\mathcal R^2d\Omega^2,\qquad
\mathcal R=r-T-\tfrac12 kT^2,\qquad \xi=1+kT,
\]

equation (4) gives

\[
R[g]=\frac{2\xi^2}{\mathcal R^2}-\frac{4k}{\mathcal R},
\qquad
R[g]\big|_{T=0}=\frac2{r^2}-\frac{4k}{r}.
\tag{6}
\]

Equation (6) agrees with direct Ricci contraction on \(\mathcal R>0\), using angular charts away from the poles.

**Target:** extend the channel and weighted-jet calculation to the off-diagonal metric or moving coframe. Native connection: EXT-005.

## 3. Area response and normal channels

**Sources:** [Mean Curvature Decomposition](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/geometric_fault_localization_and_decomposition/Mean Curvature Decomposition.md:1>) and [Three-Channel extension, section 2](</home/williaml/Cella Framework/Papers_Library/05_expository_companions_and_research_maps/dbp_role_channel_and_orbit_geometry/Three_Channel_KG_New_Math_Extension.md:125>). Their Euclidean hypersurface trace split extends as follows.

**Extension derived here.** Let a spacelike cut have positive metric \(q_{AB}\) and let \(X\) be a normal deformation. Define

\[
B_X=\tfrac12(\mathcal L_Xg)|_{TS},\qquad
\theta_X=\operatorname{tr}_q B_X.
\tag{7}
\]

Differentiating the determinant, by its multilinearity, gives

\[
\mathcal L_X\log dA
=\tfrac12 q^{AB}(\mathcal L_Xg)_{AB}
=\theta_X.
\]

Consequently any specified decomposition \(B_X=B_X^{(c)}+B_X^{(s)}\) gives

\[
\theta_X=\theta_X^{(c)}+\theta_X^{(s)}.
\tag{8}
\]

In the current normal basis,

\[
g(K,K)=-F,\quad g(K,J)=1,\quad g(J,J)=0,\qquad
D=\theta_K,\quad\rho=\theta_J>0,
\]

the area covector has metric dual and annihilating normal

\[
H_\theta=\rho K+(D+\rho F)J,\qquad
Z=\rho K-DJ.
\]

Direct multiplication gives

\[
\theta(Z)=0,\quad g(H_\theta,Z)=0,\quad
g(Z,Z)=-\rho^2F-2\rho D.
\tag{9}
\]

For each channel put \(D_a=\theta_K^{(a)}\), \(\rho_a=\theta_J^{(a)}\). Then

\[
Z=Z_c+Z_s,\qquad Z_a=\rho_aK-D_aJ.
\tag{10}
\]

Since \(\rho>0\), the area covector has a one-dimensional normal-plane kernel, proving uniqueness of the area-preserving direction. Norms of the channel sum include bilinear cross terms.

**Target:** derive \(B_K^{(c)},B_K^{(s)}\) from native cut transport and calculate \(D\). Regular coincidence requires \(D=F\Xi\) and the existing orientation condition. Native area and transport: EXT-002/003.

## 4. Constitutive dynamics and reduction

**Source:** [DIS chapter 10](</home/williaml/Cella Framework/research/derivations/DIS_chapters/10_bipolar_axes_and_constitutive_resistance.md:337>) proposes material/environment-dependent resistance and the candidate equation

\[
I_{\rm eff}(C,X)\ddot\alpha+
\gamma_{\rm eff}(C,X)\dot\alpha+
\partial_\alpha U(\alpha;C,X)
=F_{\rm ext}(C,X),
\]

with evolving material state and possible history dependence. [DIS chapter 12](</home/williaml/Cella Framework/research/derivations/DIS_chapters/12_internal_capacity_energy_and_final_hypotheses.md:180>) proposes stored energy and work/response pairings.

Status: constitutive and energy ansätze; coefficients and native identification remain to be derived.

**State extension.** Retain the internal variables required by the constitutive law:

\[
y=(a,x,d,z),\qquad
\dot y=(A,P,Q,Z)(y,X),
\tag{11}
\]

Here \(z\) contains internal variables and \(X\) environmental input. Use a finite state only when it represents the selected memory law; otherwise retain the history dependence.

At any twice differentiable transverse crossing, \(x=d=0,\ Q_0\ne0\), the completed pinch theorem still gives

\[
\xi_0=\frac{P_0}{Q_0},\qquad
\dot\xi_0=\frac{(\mathcal DP)_0Q_0-P_0(\mathcal DQ)_0}{2Q_0^2},
\tag{12}
\]

where, for a differentiable finite-state law,

\[
\mathcal D=
A\partial_a+P\partial_x+Q\partial_d+
Z\cdot\partial_z+\dot X\cdot\partial_X
\]

with an explicit time derivative added if present. Equation (12) follows by dividing the second-order Taylor expansions.

**Target:** derive the response, calculate the crossing rates (12) and obtain area transport (8). Apply the trapping criterion to the resulting trajectory.

### Reduction theorem

Let \(E\) be the full constitutive state manifold, \(N\) a candidate native state manifold, and \(\pi:E\to N\) a smooth surjective submersion. Let \(V(s,u)\) be a smooth constitutive vector field at fixed admissible input \(u\), whose input meaning has also been identified on the native side. A smooth native vector field \(v\) with

\[
D\pi_s V(s,u)=v(\pi(s),u)
\tag{13}
\]

exists **if and only if** the left-hand side is constant on each fiber \(\pi^{-1}(n)\).

**Proof.** Necessity follows from dependence on the fiber label alone. For sufficiency, define the reduced field by the common fiber value. Local smooth sections of the submersion give smoothness; the chain rule gives the projected evolution.

Example: \(\dot q=z,\dot z=0\) does not close under \(\pi(q,z)=q\), because the projected rate varies with \(z\). Retain \(z\) or derive a history-dependent reduction.

Work covectors transfer by the dual map:

\[
e_E=(D\pi_s)^T e_N
\quad\Longrightarrow\quad
\langle e_E,V\rangle=\langle e_N,D\pi_sV\rangle.
\tag{14}
\]

If \(H_E=H_N\circ\pi\), the inequality \(\dot H_E\le\langle e_E,V\rangle\) also descends. Otherwise retain the hidden energy/exchange terms.

For the displayed one-angle candidate, write \(m=I_{\rm eff}(C,X)\), \(\gamma=\gamma_{\rm eff}(C,X)\), \(v_\alpha=\dot\alpha\), and

\[
H=\tfrac12mv_\alpha^2+U(\alpha,C,X).
\]

Multiplying its equation of motion by \(v_\alpha\) and differentiating \(H\) gives exactly

\[
\dot H=F_{\rm ext}v_\alpha-\gamma v_\alpha^2
+\tfrac12\dot m\,v_\alpha^2
+U_C\cdot\dot C+U_X\cdot\dot X.
\tag{15}
\]

**Target:** derive the configuration response, field constitutive operator and cut transport from the same interaction. Equations (13)–(15) constrain reductions and energy exchange.

## Verification and dependencies

[Exact runner](/home/williaml/seated-root/cella_horizon_reuse.py) and [results](/home/williaml/seated-root/docs/cella-horizon-reuse-checks.json): 17 checks covering the coupling form, signed curvature and area/normal identities. Cella's weighted-jet replay: six assertions.

[Dependency ledger](/home/williaml/seated-root/docs/EXTERNAL-MATHEMATICS-DEBT.md): imported assumptions, native derivation targets and source provenance.

# Cella sources for the horizon and load continuation

5 September 2026.

**The earlier corpus check was insufficient. Cella supplies relevant proved algebra, a curvature calculus and an existing constitutive research direction. These change the work to do next.** This note records the sources actually inspected, proves two bounded extensions needed by the current horizon calculation, and connects them to the external-dependency ledger. No new online sources were used.

The first priority remains deriving the native seat/load response. The cited DIS chapters are legitimate physics research in the Cella corpus, as William has confirmed after checking their paths. Their material/environment response, storage and memory programme is an internal source for that derivation. Section 4 develops its connection to the native variables, including a reduction criterion where additional internal variables are eliminated. The Cella coupling-form, curvature and area-channel results remain as described below.

## 1. An internal algebraic source of the indefinite form

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

which proves (1). No spacetime metric is needed for this algebra.

Let \(e_0=(1,1,1)/\sqrt3\), and choose orthonormal \(e_1,e_2\) perpendicular to \(e_0\). Then \(Be_0=-e_0\) and \(Be_i=2e_i\). Consequently

\[
T=(e_0,e_1/\sqrt2,e_2/\sqrt2),\qquad
T^TBT=\operatorname{diag}(-1,1,1).
\tag{2}
\]

This is an explicit algebraic bridge to the native three-dimensional bilinear form. Cella already derives a form with the required inertia from its determinant structure. The isometry identifies the two bilinear spaces; selecting native clocks, rulers and their actual Cella representatives remains a substantive construction.

**Upgrade target.** Construct the native seat triplet from Cella relation/role data and carry its gauge transport along with it. Cella supplies the transformation law needed for this:

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

For \(g_1g_2g_3\ne0\), the determinant of the map \(a\mapsto\delta\nu\) is \(2(g_1g_2g_3)^2\). Thus an arbitrary defining-function gauge can move the edge vector arbitrarily. The concrete target is to derive the physical presentation rule, or a covariant enlarged state and readout, under which the native causal classification is well defined. Cella's exact pinning laws provide admissible transformations to investigate. Equation (2) alone does not choose that rule.

The source's transport is transport between defining-function presentations at a surface point. Turning this into transport along a physical path requires a rule relating successive relation jets. That is the specific connection target in EXT-001 and EXT-005.

## 2. Cella's curvature calculus extends to the current metric

**Existing Cella results.** The [completed local curvature calculus](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/local_curvature_and_black_hole_metrics/LOCAL_CURVATURE_CALCULUS_COMPLETE_v1.0.md:121>) proves an exact decomposition into differentiation-direction channels, followed by divisor valuations and cancellation laws. Its general manuscript states positive diagonal metrics.

The separate [variable-transverse weighted-jet theorem](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/local_curvature_and_black_hole_metrics/LEAD7_VARIABLE_TRANSVERSE_WEIGHTED_JET_THEOREM_v1.0.md:1>) already assumes nonzero leading coefficients without requiring their positivity. Its proof uses rational diagonal-metric formulas, and its replay retains arbitrary transverse coefficient functions. That existing result therefore already covers its stated fixed-signature germs.

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

with the same coefficient germs \(F_i\), finite candidate support and jet dependence as in the positive case. Signs enter the sums of coefficients at coincident exponents, so cancellation tests must use the signed sums. Any separate positivity-dependent conclusions need their own hypotheses.

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

This agrees with direct Christoffel/Ricci contraction. The horizon is within the regular patch \(\mathcal R>0\). The formulas are applied in angular charts away from coordinate poles. No change of spacetime signature is involved when the native rest normal becomes null.

**Use next.** Adopt these Cella channels and weighted jets for curvature orders, competing leading terms and cancellations in subsequent load-selected continuations. The current general null-coordinate family also allows off-diagonal terms; extending the channel calculus to its moving frame or full metric is a further target. Equation (4) covers the diagonal witness now. The transport connection's native derivation remains EXT-005.

## 3. Reuse the trace decomposition to derive the area response

**Existing Cella result.** [Mean Curvature Decomposition](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/geometric_fault_localization_and_decomposition/Mean Curvature Decomposition.md:1>) proves that the two Hessian contributions add linearly in mean curvature. [Three-Channel Mathematical Extension Notes, section 2](</home/williaml/Cella Framework/Papers_Library/05_expository_companions_and_research_maps/dbp_role_channel_and_orbit_geometry/Three_Channel_KG_New_Math_Extension.md:125>) express the result at the shape-operator trace level. Their particular implicit-hypersurface formulas use a Euclidean ambient metric and a declared shape-operator convention.

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

This is the trace mechanism from Cella, now stated directly on the positive cut metric in a Lorentzian ambient space. Its proof needs no Euclidean implicit-hypersurface formula. Selecting the two operators from native role/load data is the channel-construction target; the identity does not prescribe their values.

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

The total expansion and area-preserving normal add linearly; their squared norms include the ordinary bilinear cross terms. Since \(\rho>0\), the area covector is nonzero and its kernel in the normal plane is exactly one-dimensional. This proves the area-preserving direction without importing its interpretation from Anco.

**Use next.** Calculate \(B_K^{(c)},B_K^{(s)}\) from the transported native cut and derive \(D\) through (8). This replaces a freely specified area-response value with a channel calculation once that transport is supplied. The geometric coincidence condition at a regular boundary remains \(D=F\Xi\), with the orientation condition from the existing theorem.

The original Cella mean-curvature source explicitly declares the standard mean-curvature formula as a dependency. Its provenance stays recorded. The determinant/normal-plane derivation above is self-contained conditional on the declared cut metric and its transport; deriving those from native primitives remains EXT-002/003.

## 4. Continuing the DIS constitutive physics research

**Source identity, confirmed 5 September 2026.** The cited chapters 10 and 12 under Cella Framework/research/derivations/DIS_chapters are legitimate DIS physics research. The separate DIS-Engine project is a different body of work. This section continues the cited physics research and distinguishes its proposed constitutive laws from the results already derived.

In [DIS chapter 10](</home/williaml/Cella Framework/research/derivations/DIS_chapters/10_bipolar_axes_and_constitutive_resistance.md:337>), William explicitly makes axis resistance depend on external conditions, then [adds the material contained by the mass](</home/williaml/Cella Framework/research/derivations/DIS_chapters/10_bipolar_axes_and_constitutive_resistance.md:474>). The surrounding discussion develops the candidate

\[
I_{\rm eff}(C,X)\ddot\alpha+
\gamma_{\rm eff}(C,X)\dot\alpha+
\partial_\alpha U(\alpha;C,X)
=F_{\rm ext}(C,X),
\]

with evolving material state, and also considers dependence on prior exposure and deformation. [Chapter 12](</home/williaml/Cella Framework/research/derivations/DIS_chapters/12_internal_capacity_energy_and_final_hypotheses.md:180>) develops a stored-energy candidate and generalized work/response pairings.

These are research transcripts: the environmental/material dependence is William's stated direction; the displayed constitutive and energy ansätze are assistant formulations within that discussion. They supply a concrete programme to continue, rather than a proved unique dynamical law. Their proposed conservation and response conditions should be derived and assessed within that programme.

**Candidate shared state-space structure.** If native derivation requires additional internal variables, let the state include

\[
y=(a,x,d,z),\qquad
\dot y=(A,P,Q,Z)(y,X),
\tag{11}
\]

where \(z\) retains the material, velocity and internal variables required by the chosen constitutive law, and \(X\) is its environmental input. A finite-dimensional internal state is appropriate only when the selected memory law admits it; otherwise retain the history dependence explicitly.

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

and an explicit time derivative is added if present. Thus material evolution and environmental change contribute to the next pinch-direction coefficient through the total derivatives. Equation (12) follows by dividing the two second-order Taylor expansions; it does not require a memoryless linear response.

To continue this research into a native response calculation:

1. Construct the state, inputs, stored quantity and constitutive response from the DIS/native assumptions. If reducing a larger internal state to the seat variables, use the projection and power conditions below to determine whether that reduction preserves evolution and work.
2. From the resulting native law calculate \(P,Q,\mathcal DP,\mathcal DQ\) and the transported cut metric. Obtain \(D\) from its area variation, using the channel split where justified.
3. Insert those results into (12) and the existing trapping comparison. Determine which branch outcomes the derived response selects.

The passive matrix from the earlier continuation remains a solved admissible example. The DIS physics programme supplies an existing constitutive direction to develop beyond those freely chosen coefficients. Its response law and connection to the native variables remain the calculation to complete.

### Exact reduction target: dynamics and work must respect the state map

Let \(E\) be the full constitutive state manifold, \(N\) a candidate native state manifold, and \(\pi:E\to N\) a smooth surjective submersion. Let \(V(s,u)\) be a smooth constitutive vector field at fixed admissible input \(u\), whose input meaning has also been identified on the native side. A smooth native vector field \(v\) with

\[
D\pi_s V(s,u)=v(\pi(s),u)
\tag{13}
\]

exists **if and only if** the left-hand side is constant on each fiber \(\pi^{-1}(n)\).

**Proof.** Necessity follows because the right-hand side depends only on \(n\). For sufficiency, define \(v(n,u)\) by the common fiber value. Local smooth sections of a submersion show that this field is smooth, and the chain rule gives the projected evolution. The theorem is general and applies whenever the stated state map and vector field are supplied.

For example, \(s=(q,z)\), \(\dot q=z,\dot z=0\), and \(\pi(q,z)=q\) do not give a closed evolution on \(q\): different \(z\) values above the same \(q\) produce different rates. Retaining \(z\), or deriving a suitable history-dependent reduced law, is then the mathematical continuation target. This identifies which internal variables the constitutive theory must retain for a closed native description.

Work covectors transfer by the dual map:

\[
e_E=(D\pi_s)^T e_N
\quad\Longrightarrow\quad
\langle e_E,V\rangle=\langle e_N,D\pi_sV\rangle.
\tag{14}
\]

If the stored energy also descends, \(H_E=H_N\circ\pi\), then a full-state inequality \(\dot H_E\le\langle e_E,V\rangle\) transfers to the native system. If hidden variables store or exchange energy, retain that contribution before claiming the same inequality after reduction.

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

This identifies the exchange terms that an abstract response framework must carry when constitution and environment vary. It is a transferable energy-accounting identity. Which native quantities realize those terms remains a derivation target.

Equations (13)–(15) support development of the DIS constitutive theory and any justified reduction to the native seat variables. The proposed inertia response acts on configuration velocities, while THM-K's constitutive operator acts on differential forms. Connecting these response domains is a concrete derivation target. The scalar constitutive dynamics, field propagation and cut-area transport must be connected by the model's equations to determine a trapping prediction.

## 5. What was checked, and how the debt changes

The live Cella DAG was queried, relevant backing files were read, and text searches covered canonical papers, theorems, proof/audit material and the DIS derivation transcripts. The topic library uses symlinks; searching only its ordinary file listing is insufficient. Canonical backing paths were used here.

The DAG reports declared provenance and status. The conclusions above depend on the displayed source proofs and calculations, not on the status labels alone. Source SHA-256 values returned by the DAG matched the registered values for:

- Gauge Channel Transport Law: d7e31cb5b9041214a739f3491bb73b4f7ce92ca39cabb2abfc9d9a00f8b4fb0b.
- Mean Curvature Decomposition: ace98fbfd85a92cca4da4f348c410c81fe38a0409d57840915d98fa3359f9033.
- Completed local curvature calculus: 13aa19de2cc57f1bcd5cd87afadddf2ac66e9473d02eb65ddb24f75f3b368afb.
- Variable-transverse weighted-jet theorem: 35d29dc1cf1dfea9a074bd3d7e7d6f17b2504f23a059afe7b90cbc856fc15a18.

The [GR critical-exponent paper](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/general_relativity_and_gravity/whitepaper_gr_critical_exponents_v1.md>) and [GR applications note](</home/williaml/Cella Framework/Papers_Library/06_hypotheses_conjectures_and_research_programs/general_relativity_and_gravity/Gr_Applications.md>) provide diagnostics on specified GR relations. The [coupling field theory paper](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/cella_residue_and_coupling_theory/coupling_field_theory_paper.md>) names geometric flow, action and conservation as development targets. The DIS physics programme supplies the constitutive research direction developed in section 4.

The inspected sources did not supply a completed native spacetime soldering theorem, native physical null-propagation law, or premetric characteristic-reconstruction theorem. This is a bounded search result, not a claim of absence from every corpus artifact. Their next derivations now have identified internal inputs rather than an unexamined appeal to external frameworks.

Each of EXT-001 through EXT-007 now has a corpus reconciliation in the [dependency ledger](/home/williaml/seated-root/docs/EXTERNAL-MATHEMATICS-DEBT.md). Cella reuse retains both its own declared dependencies and the original external search receipts.

Verification: **17 exact checks pass** in [suites/cella_horizon_reuse.py](/home/williaml/seated-root/cella_horizon_reuse.py), with results in [cella-horizon-reuse-checks.json](/home/williaml/seated-root/docs/cella-horizon-reuse-checks.json). They check the coupling determinant and isometry, the gauge-map determinant, the signed curvature identity against direct Christoffel contraction for arbitrary diagonal functions in dimensions 2, 3 and 4, the existing load witness and signed reflection coefficient, and the area/normal identities. Cella's own weighted-jet replay was rerun and all six assertions passed.

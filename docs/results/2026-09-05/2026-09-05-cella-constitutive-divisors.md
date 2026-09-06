# Constitutive response at intersecting degeneracy surfaces

5 September 2026.

**Sources:** Cella's [local divisor/jet calculus](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/local_curvature_and_black_hole_metrics/LOCAL_CURVATURE_CALCULUS_COMPLETE_v1.0.md:248>), the DIS [material response](</home/williaml/Cella Framework/research/derivations/DIS_chapters/10_bipolar_axes_and_constitutive_resistance.md:474>) and [stored-energy programme](</home/williaml/Cella Framework/research/derivations/DIS_chapters/12_internal_capacity_energy_and_final_hypotheses.md:180>), and the [single-area response theorem](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-pinch-load-transport.md:197).

The extension classifies smooth passive response matrices at any finite collection of transverse degeneracy surfaces, including shared output constraints. It then gives the nonlinear integrability conditions, observable power ports, drift compensation and variable-constitution energy balance. Standard passivity, cotangent pullback and Hamiltonian energy accounting are retained mathematical tools.

## 1. Complete linear classification

Let \(z=(f_1,\ldots,f_k,\zeta)\) be local real coordinates, including both signs of each \(f_a\). The degeneracy divisor is \(\bigcup_a\{f_a=0\}\). Let

\[
j=M(z)e,\qquad M=R+J,\qquad
R=R^T\succeq0,\quad J=-J^T
\tag{1}
\]

be a smooth real response on \(n\) effort/flow channels. For each channel \(i\), prescribe a subset \(I_i\subseteq\{1,\ldots,k\}\) of surfaces on which that output must vanish for every effort:

\[
M_{i*}|_{f_a=0}=0\qquad(a\in I_i).
\tag{2}
\]

Set

\[
q_i=\prod_{a\in I_i}f_a,\quad Q=\operatorname{diag}(q_i),
\qquad u_{ij}=\prod_{a\in I_i\cup I_j}f_a.
\tag{3}
\]

**Theorem 1.** Conditions (1)–(2) hold exactly for matrices of the form

\[
\boxed{
R=Q\bar RQ,\qquad J_{ij}=u_{ij}\bar J_{ij},
\qquad \bar R=\bar R^T\succeq0,\quad
\bar J=-\bar J^T,
}
\tag{4}
\]

where the barred fields are smooth. The factors are unique by continuation from the complement of the divisor.

**Proof.** On \(f_a=0\), a constrained row gives \(R_{ii}=M_{ii}=0\). Positivity implies \(R_{ij}=0\) for every \(j\), so \(J_{ij}=0\) there as well. A smooth nonnegative scalar that vanishes on a two-sided coordinate hypersurface has at least a second-order zero. Applying smooth division successively gives

\[
R_{ii}=q_i^2r_i,\qquad r_i\ge0.
\]

The quotient remains nonnegative at each division. The positive-semidefinite Cauchy–Schwarz inequality now gives, on any sufficiently small compact neighbourhood,

\[
|R_{ij}|^2\le R_{ii}R_{jj}
=q_i^2q_j^2r_ir_j,
\qquad |R_{ij}|\le C|q_iq_j|.
\]

Fixing the other divisor coordinates away from zero shows that each required normal derivative of \(R_{ij}\) vanishes on the corresponding face. Continuity extends these vanishing derivatives across the intersections. Successive smooth division therefore gives \(R_{ij}=q_iq_j\bar R_{ij}\). Off the divisor, \(\bar R=Q^{-1}RQ^{-1}\succeq0\); continuity gives positivity on it.

Skew symmetry transfers a row condition to the corresponding column. Thus \(J_{ij}\) vanishes on every face in \(I_i\cup I_j\), and division gives exactly the factor \(u_{ij}\). Multiplication proves the converse. ∎

The dissipative and energy-neutral parts have different divisor laws. A face shared by channels \(i,j\) contributes **two** factors to \(R_{ij}\) and **one** to \(J_{ij}\).

For example,

\[
M(f)=f^2I+
f\begin{pmatrix}0&1\\-1&0\end{pmatrix}
\tag{5}
\]

is passive and both outputs vanish at \(f=0\). It cannot be written \(fI\,(\bar R+\bar J)\,fI\) with smooth \(\bar J\). This shared-face case is absent from the earlier single-area theorem.

For a weight \(w_a>0\), (4) gives the Cella valuation bounds

\[
\operatorname{ord}_wR_{ij}\ge
\sum_a(\mathbf1_{a\in I_i}+\mathbf1_{a\in I_j})w_a,
\qquad
\operatorname{ord}_wJ_{ij}\ge
\sum_{a\in I_i\cup I_j}w_a.
\tag{6}
\]

The first nonzero jets of the barred fields determine higher orders and cancellations. The power is exactly

\[
e^TMe=(Qe)^T\bar R(Qe).
\tag{7}
\]

Under a smooth invertible power-preserving port change \(j'=Tj,\ e'=T^{-T}e\), the response transforms as \(M'=TMT^T\). The vanishing output covectors transform by \(T^{-1}\); (4) is their expression in an adapted channel frame.

## 2. Nonlinear response and integrability

Let \(e\) range over an open convex domain containing zero, and let \(j(z,e)\) be smooth. Write

\[
b(z)=j(z,0),\qquad M(z,e)=\partial_ej.
\]

**Theorem 2.** A response is incrementally passive,

\[
(e-e')^T\bigl(j(z,e)-j(z,e')\bigr)\ge0,
\tag{8}
\]

and obeys \(j_i|_{f_a=0}=0\) for all efforts and \(a\in I_i\), exactly when:

1. \(b_i|_{f_a=0}=0\);
2. \(M\) has the factorization (4), with barred fields smooth in \(z,e\);
3. each row of \(M\,de\) is closed:

\[
\partial_{e_\ell}M_{ij}=\partial_{e_j}M_{i\ell}.
\tag{9}
\]

For such data the unique response with offset \(b\) is

\[
\boxed{j(z,e)=b(z)+\int_0^1M(z,te)e\,dt.}
\tag{10}
\]

**Proof.** Incremental passivity is equivalent to \(\operatorname{sym}\partial_ej\succeq0\), by differentiating (8) and integrating along effort segments. Boundary vanishing gives both the offset and Jacobian row conditions. Theorem 1 applies with \(e\) as extra smooth parameters. Conversely, (9) makes the row forms exact on the convex domain, so (10) has Jacobian \(M\). Its positive symmetric part proves (8), and its boundary rows vanish. ∎

This nonlinear theorem classifies **incremental** passivity. The linear theorem requires only the usual power condition \(e^TMe\ge0\).

A nonlinear potential example is

\[
\Psi(f,e)=\tfrac12(e_0+fe_1)^2+\tfrac14f^2e_1^4,
\qquad j=\nabla_e\Psi.
\]

Its Hessian factors with \(Q=\operatorname{diag}(1,f)\) and

\[
\bar R=
\begin{pmatrix}1&1\\1&1+3e_1^2\end{pmatrix}\succeq0.
\]

The second output vanishes for every effort at \(f=0\), while its nonlinear dependence is retained. Arbitrarily assigning barred matrices without (9) would not define a constitutive law.

## 3. Derived observable and area power ports

Suppose the full constitutive flow is \(j_s\), and measured rates have the affine form

\[
y=C(s)j_s+b(s).
\tag{11}
\]

For an observable effort \(e_y\), the state effort is \(e_s=C^Te_y\). Their power identity is

\[
\boxed{e_s^Tj_s=e_y^T(y-b).}
\tag{12}
\]

Thus \(b\) is an exchange/drift term with power \(e_y^Tb\), not part of the effort-driven mobility. For a linear state response,

\[
y-b=M_y e_y,\qquad M_y=CM_sC^T.
\tag{13}
\]

Positivity and skew symmetry are preserved. Theorem 1 applies to the observable mobility \(M_y\); full boundary matching also requires the total offset to vanish. For nonlinear state response use \(y=Cj_s(C^Te_y)+b\), whose Jacobian is \(C(\partial_{e_s}j_s)C^T\).

Use a positive **local cut area density** \(\mathcal A(s,r)\) in labels transported by the normal flow, and let \(K=C_{\rm clock}+wR_n\). Then

\[
D=\theta_K
=d_s\log\mathcal A\cdot\dot s+b_A,\qquad
b_A=wR_n\log\mathcal A+
(C_{\rm clock}\log\mathcal A)_{\text{explicit}}.
\tag{14}
\]

If a cut parametrization moves with velocity \(K+v_\parallel\), subtract \(\operatorname{div}_h v_\parallel\) in the drift to recover \(\theta_K\). The derivative of an integrated cut area gives its area-weighted average expansion; pointwise trapping uses the local density or the [implicit-cut trace formula](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-cella-normal-geometry.md), not that average.

The induced area effort pulls back as \(e_s=e_A\,d_s\log\mathcal A\). For spherical cuts \(\mathcal A\propto\mathcal R^2\),

\[
D=\frac{2\dot{\mathcal R}}{\mathcal R}
+\frac{2wR_n\mathcal R}{\mathcal R},\qquad
e_{\mathcal R}=\frac{2e_A}{\mathcal R}.
\tag{15}
\]

This derives the area work channel from the chosen geometric area functional and records its spatial drift. A moving cut or material background can add an offset even when applied effort is zero.

In the horizon comparison, \(\theta_L=D+\rho F/2\), with \(\rho>0\). For one area channel, Theorem 1 recovers the earlier factorization. With offset \(b_A\), universal boundary coincidence requires both the area row of \(M_y\) and the total offset to vanish at \(F=0\). Write

\[
D=F(\beta+c\cdot e_0+Fr e_A).
\]

At fixed geometry the orientation factor is

\[
1+\frac{2(\beta_H+c_H\cdot e_0)}{\rho_H}.
\tag{16}
\]

For all unbounded efforts it is positive exactly when \(c_H=0\) and \(1+2\beta_H/\rho_H>0\). A drift need not vanish to second order. For bounded efforts, test (16) on the admissible effort set. These conditions distinguish universal matching from a crossing selected by a particular applied load.

## 4. Variable constitution and momentum exchange

Use the DIS stored-energy candidate in momentum variables:

\[
H(q,p,z)=\tfrac12p^TM(q,z)^{-1}p+U(q,z),\qquad
v=\partial_pH.
\tag{17}
\]

Here \(M\) is positive definite and \(z\) includes material state and environmental parameters. The general forced canonical equations

\[
\dot q=v,\quad
\dot p=-\partial_qH-\Gamma v+Bu+\Pi,\quad
\dot z=W,\qquad\Gamma=\Gamma^T\succeq0
\]

have energy balance

\[
\boxed{\dot H=u^TB^Tv+v^T\Pi-v^T\Gamma v+
(\partial_zH)^TW.}
\tag{18}
\]

\(\Pi\) is momentum carried by constitutive exchange; the last term accounts for changes in material/environmental storage. These ports can be coupled to additional passive subsystems using the same power pairing.

For the source's one-angle equation

\[
m(C,X)\ddot\alpha+\gamma\dot\alpha+U_\alpha=f,
\]

\(p=m\dot\alpha\) gives \(\Pi=\dot m\,\dot\alpha\) relative to the canonical momentum equation. Its energy balance is therefore

\[
\dot H=f\dot\alpha-\gamma\dot\alpha^2
+\tfrac12\dot m\,\dot\alpha^2+U_C\dot C+U_X\dot X.
\tag{19}
\]

The alternative \(\Pi=0\) includes the term \(\dot m\,\dot\alpha\) in the angle equation. Specifying constitutive momentum exchange selects between these laws; changing inertia alone does not select it.

### Kinetic cost and constrained response

[KIN-1](/home/williaml/seated-root/kin1_cost_fork.py) supplies the cost-side alternatives. A kinetic form \(h:T\mathcal S\to T^*\mathcal S\) and a mobility \(M:T^*\mathcal S\to T\mathcal S\) have opposite types. Even the inverse of a smooth positive cost cannot by itself annihilate a nonzero constrained-output covector \(\ell\), since

\[
\ell h^{-1}\ell^T>0.
\]

Thus finite positive cost at the horizon is compatible with Theorem 1: the load response must supply its required divisor factors through the response or observable map. For any chosen smooth positive cost, \(Qh^{-1}Q\) is one admissible reciprocal response family, with the remaining rate scale and material selection explicit. This construction does not select among KIN-1's costs. The divisor theorem determines the boundary orders independently of that choice. Comparing the specific KIN-1 fixtures is a separate open-circuit test.

## 5. Smooth response across a weighted resolution

Cella's weighted approach calculus is already applied to the native pinch in [PINCH-3](/home/williaml/seated-root/pinch3_cella_valuation.py). Its response extension follows the power pairing through the resolved coordinates.

For a weighted chart \(z_i=r^{w_i}u_i\), \(w_i\) positive integers, use the original rate coframe

\[
dz_i-w_i r^{w_i-1}u_i\,dr=r^{w_i}du_i,\qquad dr=dr.
\]

Thus the resolved-to-original rate map is

\[
Q_w=\operatorname{diag}(r^{w_1},\ldots,r^{w_m},1).
\]

**Theorem 3.** Let \(M(r,u)\) be the original-coframe response pulled back to the punctured resolved chart. It has a smooth resolved response \(\widehat M\) exactly when

\[
\boxed{M=Q_w\widehat M Q_w}
\tag{20}
\]

with \(\widehat M\) smooth. Equivalently, each \(M_{ij}\) is smoothly divisible by \(r^{w_i+w_j}\), assigning weight zero to the radial channel. The lift is unique. Passivity on the punctured chart transfers to the lift by continuity.

**Proof.** Rate pushforward is \(j=Q_w\widehat j\); power preservation gives \(\widehat e=Q_w e\). Substitution gives (20). Off \(r=0\), its inverse is \(\widehat M=Q_w^{-1}MQ_w^{-1}\). The divisibility conditions are exactly its smooth extension conditions; congruence and continuity preserve positivity. ∎

Unlike boundary passivity alone, a smooth geometric lift constrains the energy-neutral entries to the same weight sums. The shared-face example (5) has no such weight-one lift for both channels.

For the native chart \(x=d\xi\), the contrast/radial coframe is

\[
dx-\xi\,dd=d\,d\xi,\qquad dd=dd.
\]

A constant resolved material response

\[
\begin{pmatrix}\dot\xi\\\dot d\end{pmatrix}
=
\widehat M
\begin{pmatrix}\widehat e_\xi\\\widehat e_d\end{pmatrix},
\qquad
\widehat M=
\begin{pmatrix}\mu&\beta\\\beta&\nu\end{pmatrix},
\quad \mu>0,\ \nu>0,\ \mu\nu>\beta^2
\tag{21}
\]

induces, in the original coordinate rates \((\dot x,\dot d)\),

\[
\boxed{
M_{\rm orig}=
\begin{pmatrix}
\mu d^2+2\beta d\xi+\nu\xi^2&\beta d+\nu\xi\\
\beta d+\nu\xi&\nu
\end{pmatrix},\qquad
\det M_{\rm orig}=d^2(\mu\nu-\beta^2).
}
\tag{22}
\]

The original effort \((e_x,e_d)=(0,1)\) pulls back to \((0,1)\). Starting at \(d=0,\xi=\xi_0\), the solution is

\[
d=\nu T,\qquad \xi=\xi_0+\beta T,\qquad
x=\nu\xi_0T+\nu\beta T^2.
\tag{23}
\]

For \(a=0,\xi_0=1,\nu=1,\beta=k\), this is the existing finite-time native crossing \(d=T,\ x=T+kT^2\). Its resolved material coefficients are constant; their original-coordinate dependence follows from the weighted chart.

The earlier full-rank matrix example and (22) generate the same chosen trajectory but differ for other efforts. The earlier matrix has no smooth full response lift at \(d=0\): its determinant remains one. Equation (22) supplies a constitutive continuation defined for every effort on the resolved state. These are original coordinate components **over the resolved state**; they need not descend to a single-valued tensor at the collapsed point. Distinct values of \(\xi\) remain distinct states.

**Corollary 3.1 — Descent obstruction for transverse response.** Let \(\widehat M(d,\xi)\) be any smooth response, without requiring symmetry or passivity. Suppose its coordinate pushforward

\[
M_{\rm orig}=J_\pi\widehat M J_\pi^T,\qquad
J_\pi=\begin{pmatrix}d&\xi\\0&1\end{pmatrix}
\]

descends continuously to a single matrix \(M_0\) at \((x,d)=(0,0)\), for every \(\xi\) in an open interval. Then

\[
\boxed{M_0=0,\qquad \widehat M_{dd}(0,\xi)=0,\qquad
\dot d|_{d=0}=0\quad\text{for every finite original effort}.}
\tag{24}
\]

**Proof.** At the exceptional fibre, set \(n(\xi)=\widehat M_{dd}(0,\xi)\). The pushforward is

\[
M_0=n(\xi)\begin{pmatrix}\xi^2&\xi\\\xi&1\end{pmatrix}.
\]

Its bottom-right entry makes \(n\) constant; its off-diagonal entry then forces that constant to vanish on an open interval. The pulled-back effort is \((d e_x,\xi e_x+e_d)^T\), hence \(\dot d=\widehat M_{d\xi}d e_x+\widehat M_{dd}(\xi e_x+e_d)=0\) at \(d=0\). ∎

For a smooth effort field the exceptional fibre is therefore invariant under the unique homogeneous resolved flow. If an affine offset must also descend, its boundary pushforward is \(\widehat b_d(\xi)(\xi,1)^T\); the same argument forces \(\widehat b_d=0\). A finite transverse passage such as (23), with \(\nu>0\), necessarily retains directional state or changes another stated hypothesis. The obstruction does not apply to a single prescribed ray, a singular original effort, or an offset whose descent is not required.

Area transport must still be calculated by section 3. Reproducing the native crossing does not impose \(D=0\).

## 6. Native use and closure

The divisor theorem gives the complete linear passive class for arbitrary channel/face incidence. The nonlinear theorem adds the exact integrability requirement. The weighted-lift theorem distinguishes these passive laws from those admitting a smooth resolved continuation and constructs the latter at the native pinch. The area port is a pullback of geometric work, with its drift carried explicitly. The momentum formulation completes the variable-constitution bookkeeping.

To evaluate a native prediction, supply the DIS energy/response coefficients, material exchange and geometric cut evolution. Extract \(P=\dot x,\ Q=\dot d\), then use

\[
\xi_0=P_0/Q_0,\qquad
\dot\xi_0=\frac{\dot P_0Q_0-P_0\dot Q_0}{2Q_0^2}
\]

at a twice differentiable transverse pinch. Equations (14)–(16) decide the corresponding area/trapping behavior. The theorems classify and connect these data; they do not select a unique material law.

**Dependencies:** EXT-004 retains standard passive/Hamiltonian structure; EXT-007 retains smooth division and local differential calculus. No new online source was used. The normal-crossing hypothesis is substantive: tangent or singular divisor intersections require their resolved charts before applying (4).

**Verification:** [runner](/home/williaml/seated-root/cella_constitutive_divisors.py), [results](/home/williaml/seated-root/docs/receipts/cella-constitutive-divisors-checks.json).

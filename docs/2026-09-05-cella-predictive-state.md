# Minimum state for predictive constitutive evolution

5 September 2026

**For the existing constant resolved response, both independent original effort channels make `(xi,d)` the exact minimum smooth state for predicting even the lapse `F`, with `a` and the material coefficients fixed.** This is two continuous coordinates, including at the pinch. The direction replaces the collapsed coordinate `x`; it is not a third independent coordinate. Restricting the available inputs can reduce this minimum. If all three constant material coefficients are unknown and may differ between preparations, five coordinates are necessary and sufficient for the same family's signed-volume predictions. Material memory has its own, calculable state count.

These are results for specified laws, inputs and outputs. The general construction below determines the minimum for another constitutive law; it does not choose that law. All trajectory statements are local, on a common admissible domain with unique full-state solutions.

## 1. Predictive equivalence and an effective closure test

Let `z` be a smooth finite-dimensional state, `h(z)` the outputs to predict, and

\[
\dot z=X_0(z)+\sum_{i=1}^m u_i X_i(z).
\tag{1}
\]

The same piecewise constant input is applied to compared states. The admissible input set has nonempty interior, so its affine span distinguishes all the fields in (1). If inputs transform with a presentation change, that transformation must instead be part of the comparison. Known environmental signals are inputs; an unknown evolving environment is additional state. Explicit time dependence can be included by adjoining `t`, with `dot t=1`.

Two states are **predictively equivalent** when every common admissible future input produces the same output trajectory. This relation is preserved by each common input evolution: otherwise concatenating the distinguishing continuation would distinguish the initial states. Every deterministic state reduction preserving the outputs and their evolution must separate inequivalent states. Conversely, the equivalence classes carry the induced input evolution and output. This supplies the exact minimum as a set with evolution; it need not be a smooth manifold.

The useful calculation is stronger than this characterization. Write `L_i f=df(X_i)` and form the observation family

\[
\mathcal H=\{h_j,L_{i_1}h_j,L_{i_2}L_{i_1}h_j,\ldots\}.
\tag{2}
\]

Equality of all future outputs implies equality of every member of (2). Indeed, compose short constant-input flows and differentiate with respect to their successive durations at zero. Mixed derivatives give the ordered Lie derivatives; affine spanning isolates `X_0,...,X_m`.

**Finite closure certificate.** Suppose a finite collection `q=(q_1,...,q_r)` of functions constructed from (2) satisfies

\[
h=H(q),\qquad L_iq=Y_i(q),
\tag{3}
\]

where the displayed reduced fields are locally Lipschitz. Then `q` supports unique predictive evolution `dot q=Y_0(q)+sum u_iY_i(q)`. If equality of `q` is equivalent to equality of all of (2), its fibres are exactly the predictive classes. At a point where `Dq` has rank `r`, any other smooth predictive state through which these functions factor smoothly has at least `r` coordinates. Thus a rank-`r` certificate is locally dimension-minimal among regular smooth reductions.

**Proof.** The chain rule and ODE uniqueness prove sufficiency. Any predictive reduction must retain each successive output derivative, hence every `q_i`. Differentiating their factorization through that reduction proves the rank lower bound. This proves minimality from separating functions, rather than assuming that the number of current readouts is the state dimension. Local sections justify smooth factorization for a regular submersion. At a singular quotient, the displayed regularity and global fibre checks are additional requirements. ∎

Energy is an output if its exact value must also survive the reduction. Add it to `h` before calculating (2). Input/output prediction alone does not assert preservation of hidden storage.

The [selected-skeleton implementation](</home/williaml/Cella Framework/engine/src/cella/continuation/selected_skeleton.py:58>) now identifies arrows by `(source,target,monodromy_parity)`; `route_class` and `morphism_digest` retain provenance without affecting equality or hashing. Those evidence fields do not automatically become state coordinates. Material history enters `z` only when the response law makes it affect future outputs. Arrow identity and predictive state equivalence are separate, explicitly typed comparisons.

### A terminating polynomial calculation, including discrete branches

Assume the fields and outputs are polynomial over a characteristic-zero field, with their coefficients fixed. On two copies of state space use paired derivatives

\[
\mathscr D_i=X_i(z)\cdot\partial_z+X_i(z')\cdot\partial_{z'},\qquad
I_0=\langle h_j(z)-h_j(z')\rangle,
\]
\[
I_{n+1}=I_n+\sum_i\mathscr D_i I_n.
\tag{4}
\]

For a finite generating set of `I_n`, it suffices to add its derivatives: differentiating polynomial multiples adds only existing ideal elements and multiples of these derivatives. The polynomial ring is Noetherian, so this increasing sequence stabilizes after finitely many strict increases. A Gröbner remainder test certifies stabilization; merely reaching a chosen iteration limit does not.

**Theorem.** On the admissible real domain,

\[
z\sim z'\quad\Longleftrightarrow\quad (z,z')\in V_{\mathbb R}(I_\infty).
\tag{5}
\]

**Proof.** The stable ideal is exactly the ideal generated by the differences of all words (2). Necessity follows from the duration-derivative argument. For sufficiency, let `g` be a finite generating vector of the stable ideal. Stability gives `mathscr D_i g=A_i(z,z')g` with polynomial matrix coefficients. Along any common-input pair of trajectories, `g` therefore solves a homogeneous linear ODE with zero initial value. It remains zero, including the original output differences. This avoids a generic-rank or Taylor-truncation assumption. ∎

The zero set can have several components: these encode actual alternative states with identical future outputs. A local Jacobian rank alone would miss them. Restrict the final relation to the declared real admissibility domain. For rational fields, use the corresponding localized ring and keep its denominators nonzero; the finite polynomial runner implements (4) directly, not an unverified extension through poles. When parameters are symbolic coefficients, specialize degenerate parameter strata separately.

### Smooth laws need a finite certificate or actual flows

There is no replacement of smooth future evolution by its entire jet at one state. Set

\[
f(z)=\begin{cases}e^{-1/z^2},&z>0,\\0,&z\le0,\end{cases}\qquad
\binom{\dot y}{\dot z}=
\begin{pmatrix}f(z)&0\\0&1\end{pmatrix}\binom{u_1}{u_2},\quad h=y.
\tag{6}
\]

The response is smooth and positive semidefinite. States `(y,z)=(0,0)` and `(0,-1)` agree on every function (2), because all derivatives of `f` vanish at both points. Under the same input `(1,1)`, the first has `y(t)=integral_0^t f(s)ds>0`, while the second has `y(t)=0` for `0<t<1`. Thus smooth flat material dependence can defeat even infinite local-jet equality. Equations (3) and (5) are sufficient certificates in their stated classes; (6) explains why one cannot silently extend the polynomial conclusion to arbitrary smooth laws.

## 2. Response descent, branch transport and regularity

Let `pi:S -> O` be a candidate observable state map and put `C=Dpi`. For the power-paired response

\[
\dot z=b(z)+\widehat M(z)C_z^Te,\qquad
\dot y=C_zb(z)+C_z\widehat M(z)C_z^Te,
\tag{7}
\]

an instantaneous law on `y` exists exactly when both coefficients

\[
c(z)=C_zb(z),\qquad B(z)=C_z\widehat M(z)C_z^T
\tag{8}
\]

factor through `pi`, for all admissible efforts. If the factors are locally Lipschitz, this gives a unique reduced evolution for these same ports. Positive semidefiniteness of `sym Mhat` transfers to `sym B` by congruence; the power pairing is unchanged. For a nonlinear response, require the whole pushed response to factor, including its zero-effort offset. The existing constitutive divisor factors can then be imposed before or after a valid descent in an adapted port frame.

This extends the [earlier submersion criterion](2026-09-05-cella-dependency-reconciliation.md#exact-reduction-target-dynamics-and-work-must-respect-the-state-map) to an explicit test at singular observations. Factoring coefficients and establishing a unique reduced ODE are separate checks. For example, the smooth flow `dot z=1` under `y=z^3` projects to `dot y=3|y|^(2/3)`. Its coefficient is single-valued, but starting at zero this equation allows arbitrary waiting times before `y=(t-T)^3`. The lifted flow selects `T=0`. Retaining the regular coordinate `z` adds no information, but restores a unique smooth evolution law.

For a finite deck group `Gamma` and the quotient `pi`, equivariance of each field under the **same input**,

\[
X_i(gz)=Dg_zX_i(z),
\tag{9}
\]

makes deck partners predictively equivalent for deck-invariant outputs. On a dense regular locus of a finite covering, it is also necessary for a vector field to descend through that full quotient: differentiate `pi g=pi` and invert `Dpi` there. Continuity extends (9) upstairs to the branch locus. For mobility and drift the corresponding sufficient rules are

\[
\widehat M(gz)=Dg_z\widehat M(z)Dg_z^T,\qquad
b(gz)=Dg_zb(z).
\tag{10}
\]

These rules test which algebraic sheet identifications are compatible with the dynamics. Kummer relations determine the available deck action; they do not impose (9). An input transformed by `g` defines a transformed experiment, which is a different comparison from holding the input fixed. The smoothness/uniqueness requirement in the chosen quotient coordinates remains as in (8).

## 3. The pinch retains a continuous direction

Use the existing finite-direction chart

\[
z=(\xi,d),\quad x=d\xi,\quad
F=(1-\xi^2)/p,\quad p=1+a^2>0,
\quad J_\pi=\begin{pmatrix}d&\xi\\0&1\end{pmatrix}.
\tag{11}
\]

Here `a` and the pinch chart label `sigma` are initially fixed. The old observables `(x,d)` collapse every finite `xi` at `d=0`.

**Theorem — transverse directions are distinguishable.** For any smooth resolved flow, at the exceptional fibre

\[
\dot x=\xi\dot d.
\tag{12}
\]

If a state there admits an input with `dot d != 0`, no state at another `xi` can have the same `(x,d)` future for every common input. Equality of the first rates for that input would give the same nonzero `dot d`, and (12) would force the same `xi`.

For (7), write `n(xi)=Mhat_dd(0,xi)`. The incremental response gives the constructive measurement

\[
B(0,\xi)=n(\xi)
\begin{pmatrix}\xi^2&\xi\\\xi&1\end{pmatrix},\qquad
\boxed{\xi=B_{xd}/B_{dd}\quad(n\ne0).}
\tag{13}
\]

Offsets cancel from this incremental measurement. If transverse drift supplies the crossing instead, (12) supplies its direction. This strengthens the [continuous-descent obstruction](2026-09-05-cella-constitutive-divisors.md#5-smooth-response-across-a-weighted-resolution): it identifies exactly which states the response distinguishes. A finite sheet label cannot replace an interval of distinguishable directions. The chart is still two-dimensional: `(xi,d)` replaces `(x,d)`.

## 4. Even the lapse requires both coordinates for the present load

Take the constant resolved family already constructed, including its positive-semidefinite boundary:

\[
\widehat M=\begin{pmatrix}\mu&\beta\\\beta&\nu\end{pmatrix},
\quad\mu>0,\ \nu>0,\ \mu\nu-\beta^2\ge0,
\quad b=0.
\tag{14}
\]

The controlled fields generated by the original efforts `(e_x,e_d)` are

\[
X_x=\binom{\mu d+\beta\xi}{\beta d+\nu\xi},\qquad
X_d=\binom{\beta}{\nu}.
\tag{15}
\]

**Theorem — exact minimum for the lapse.** With the two input channels independently available and the coefficients and `a` known, states have identical future `F` for all inputs **if and only if** their `(xi,d)` agree. Two smooth coordinates are necessary and sufficient, including at `xi=0` and `d=0`.

**Proof by reconstruction.** Put `c=mu nu+beta^2>0`. Direct differentiation gives

\[
L_dF=-2\beta\xi/p,\quad
L_xF=-2\xi(\mu d+\beta\xi)/p,
\]
\[
(L_dL_x-L_xL_d)F=-2c\xi/p.
\]

Consequently the following functions of the input/output derivative data recover the entire state:

\[
\boxed{
\xi=-\frac{p}{2c}(L_dL_x-L_xL_d)F,\qquad
d=\frac{L_x\xi-\beta\xi}{\mu}.
}
\tag{16}
\]

These formulas use derivatives of order at most three and have no state-dependent denominator. They hold also when `beta=0`. The Jacobian of the reconstructed pair with respect to `(xi,d)` is the identity. Section 1 therefore proves both global separation on the chart and the two-dimensional lower bound. The polynomial fields (15) supply sufficiency and regular local evolution. ∎

This is stronger than computing a Kummer sheet or declaring `F` incomplete: it supplies a finite input/output certificate and a minimum attained by the existing chart. No new material term has been introduced.

### A decisive horizon experiment

At `d=0, xi=+/-1`, both states have `F=0`. If `beta!=0`, the same radial effort already gives opposite lapse derivatives `-2 beta xi/p`.

Even when `beta=0`, use the same two-step input: first `e_d=1,e_x=0` for a duration `t`, then switch to `e_x=1,e_d=0`. After the first step, `d=nu t` and `xi` is unchanged. Immediately after switching,

\[
\dot F=-2\mu\nu t\xi/p.
\tag{17}
\]

The two null directions therefore produce opposite signs of `dot F`. Inputs can be arbitrarily small-duration so both trajectories stay in the native real chart. The experiment predicts departure to opposite signs of the rest-frame readout; gravitational trapping additionally uses the independently derived cut-expansion law.

Geometrically the deck is `tau(xi,d)=(-xi,-d)`. In (15), `X_x(tau z)=Dtau X_x(z)` but `X_d(tau z)=-Dtau X_d(z)`. Deck covariance transforms `e_d` to `-e_d`. Thus equal `F` under the deck does not identify the two **fixed-input** experiments above.

## 5. Exact changes in the count

**Only the radial input is available.** With `e_x=0`, predicting `F` needs no `d`: `dot xi=beta e_d`. If `beta!=0`, `xi` is a one-coordinate minimum, recovered from `L_dF`. If `beta=0`, `F` is constant and `F` itself is the minimum output state; its value must be retained but no evolving coordinate is needed beyond it. These reductions are statements about restricted experiments.

**The remaining constant reciprocal passive cases.** A symmetric positive-semidefinite response has `mu,nu>=0` and `beta^2<=mu nu`. This gives a complete fixed-coefficient classification for the lapse:

- `mu>0,nu>0`: (16) proves the two-coordinate minimum even on the rank-one boundary `beta^2=mu nu`. A one-dimensional instantaneous flow range need not imply a one-dimensional predictive state.
- `mu>0,nu=0`: positivity forces `beta=0`. The only field is `X_x=(mu d,0)`. The complete observation family is generated by `xi^2,xi d,d^2`, so predictive equivalence identifies exactly `(xi,d)` with `(-xi,-d)`. A regular polynomial presentation is the cone
  \[
  (U,V,W)=(\xi^2,\xi d,d^2),\quad U,W\ge0,\quad UW=V^2,
  \quad (\dot U,\dot V,\dot W)=(2\mu V,\mu W,0)e_x.
  \tag{18}
  \]
  This is a two-dimensional quotient with a singular vertex, represented by three polynomial coordinates and one relation. Two global Lipschitz coordinates also exist: `(A,B)=(U-W,2V)`, with `U=(sqrt(A^2+B^2)+A)/2`, `W=(sqrt(A^2+B^2)-A)/2` and `(dot A,dot B)=mu(B,sqrt(A^2+B^2)-A)e_x`. These give unique evolution, while the cone presentation retains smooth polynomial formulas.
- `mu=0`: positivity forces `beta=0`. The direction is constant, so `F` alone is the complete output state, for every `nu>=0`. This includes the zero response.

The second case follows from `L_xF=-2mu xi d/p`, `L_x^2F=-2mu^2d^2/p` and higher derivatives zero. Equality of the three quadratic products means equality of the rank-one matrix `(xi,d)(xi,d)^T`, hence equality up to overall sign, including the zero vector. Equation (18) proves sufficiency on the whole quotient. Thus the branch reduction is derived from the available dynamics.

**Material constants are unknown preparation variables.** Let `(mu,beta,nu)` be constant in time but variable across preparations in the strictly positive-definite interior `mu nu>beta^2` of (14). Under the same two input channels, even the output `d` has the exact minimum

\[
(\xi,d,\mu,\beta,\nu).
\tag{19}
\]

To prove this, regard the coefficient rates as zero. Output derivatives recover

\[
\nu=L_dd,\qquad
\beta=\frac{L_dL_xd}{2\nu},\qquad
\mu=\frac{L_dL_x^2d/\nu-3\beta^2}{\nu},\qquad
\xi=\frac{L_xd-\beta d}{\nu}.
\tag{20}
\]

Here `L_dL_xd=2 beta nu` and `L_dL_x^2d=nu(mu nu+3 beta^2)`. Together with the output `d`, (20) reconstructs five independent coordinates on the open positive-definite parameter domain. They close the evolution, so the minimum is exactly five. Calibrated constants are part of the law and contribute no unknown state. A constrained family such as `(mu,beta,nu)=(1+k^2,k,1)` contributes just one unknown material coordinate: `(xi,d,k)` is the minimum for `d`.

**The seat variable `a`.** If `a` is a fixed known parameter, the counts above stand. If it varies between preparations and is among the outputs to retain, append it: for `dot a=0` the fixed-material family has minimum `(a,xi,d)`. The local discrete label `sigma=+/-1` identifies which pinch neighbourhood reconstructs `gamma=sigma sqrt(1-Fd^2)` and `b=a gamma+xi d`; it is fixed in each chart. An evolving `a`, coframe, cut or material state requires its supplied evolution in (1). Their minimum is determined by the same closure test, rather than by adding an arbitrary number of variables.

## 6. Completing the material-memory question

Cella's [DIS chapter 10](</home/williaml/Cella Framework/research/derivations/DIS_chapters/10_bipolar_axes_and_constitutive_resistance.md:444>) explicitly admits path-dependent constitution. A current geometric state is therefore sufficient only if the material memory has been included or proved unobservable under the chosen outputs and inputs.

### Exact finite-state quotient for linear material variables

Let

\[
\dot z=(A_0+\sum_i u_iA_i)z+Bu,\qquad y=Cz.
\]

Start with `N_0=ker C` and iterate

\[
N_{k+1}=N_k\cap\bigcap_{i=0}^m A_i^{-1}N_k.
\tag{21}
\]

This descending sequence stabilizes after at most `dim z` strict dimension drops. Its limit `N_*` is the largest common invariant subspace inside `ker C`. Two preparations have identical outputs under every common input exactly when their difference lies in `N_*`; the affine forcing cancels in that difference equation. Thus the minimum dimension for arbitrary initial preparations is

\[
\boxed{\dim z-\dim N_*.}
\tag{22}
\]

Invariance proves sufficiency. Differentiating outputs under successive inputs proves necessity. Equivalently `N_*=intersection ker(C A_word)`. This is the dynamic version of the selected null-module quotient: only an observation-null subspace preserved by every allowed evolution can be removed. For one `A`, ordinary matrix observability computes (22).

For the time-invariant additive-input case `dot z=Az+Bu`, preparations reachable from rest occupy `R=span{im B,...,im A^(n-1)B}`. Their minimum linear realization dimension is `dim R-dim(R intersect N_*)`. For a bilinear law, a restricted preparation set need not be a linear subspace; restrict the predictive relation to the actual admissible set before counting its dimension.

### Passive memory can require any number of coordinates

For distinct positive rates `lambda_j` and positive weights `c_j`, consider

\[
\dot z_j=-\lambda_j z_j+u,\quad y=\sum_{j=1}^n c_jz_j,
\quad H=\tfrac12\sum_jc_jz_j^2.
\]
\[
\dot H=uy-\sum_jc_j\lambda_jz_j^2\le uy.
\tag{23}
\]

The observation matrix has entries `c_j(-lambda_j)^k`, `0<=k<n`; its determinant is a nonzero weighted Vandermonde determinant. Hence precisely `n` material coordinates are needed. Equal rates combine into one observable weighted mode; zero output weights remove a mode. The kernel for preparations generated from rest is `K(t)=sum c_j exp(-lambda_j t)`, and distinct rates also give full reachability.

The continuum version `lambda in [1,2]`, with

\[
\dot z(\lambda)=-\lambda z(\lambda)+u,\quad
y=\int_1^2z(\lambda)d\lambda,
\quad H=\tfrac12\int_1^2z(\lambda)^2d\lambda,
\tag{24}
\]

is likewise passive. Its kernel is

\[
K(t)=\int_1^2e^{-\lambda t}d\lambda
=\frac{e^{-t}-e^{-2t}}{t},\quad K(0)=1.
\]

Every moment matrix `H_N=(integral_1^2 lambda^(i+j)d lambda)_(0<=i,j<N)` is positive definite: for a nonzero polynomial `p`, its quadratic form is `integral_1^2 p(lambda)^2d lambda>0`. Thus the associated Hankel rank is infinite. A finite-dimensional linear realization would have moments `C A^k B` and Hankel rank at most its state dimension, a contradiction. Nor can a finite smooth realization around a regular rest equilibrium reproduce this linear input/output operator: linearizing it at rest would give just such a finite realization.

Consequently passivity and smooth response alone impose **no finite universal upper bound** on the material-state requirement. A finite relaxation model has the exact count (22); a distributed memory law such as (24) requires distributed state or an explicitly approximate reduction. This is an extension target determined by the selected material law, not a reason to prescribe a memoryless law.

## 7. Construction delivered and provenance

The [runner](../cella_predictive_state.py) implements the terminating polynomial pair-ideal construction and checks its invariance certificate. It also verifies the native reconstructions, the zero-coupling case, the horizon experiment, unknown-material counts, power and branch transport, delayed distinguishability, a nontrivial removable state, and finite material-memory reductions. The [receipt](cella-predictive-state-checks.json) records its exact checks; it supplements the proofs above.

Inputs read for this extension:

- [Selected Quotient Groupoids foundation](</home/williaml/Cella Framework/Papers_Library/05_expository_companions_and_research_maps/dbp_role_channel_and_orbit_geometry/SELECTED_QUOTIENT_GROUPOIDS_FOUNDATION_v1.0.md>), DAG `DBP:thm:sqg_foundation`, live hash `7a0bb32bf1642c3cd772ae150bffb5297bb92773c5e50743018a748e792f9043`: native morphisms preserve the selected null submodule. Equations (21)–(22) compute the largest removable subspace compatible with the specified dynamics; they do not identify coefficient modules with state manifolds.
- [General conjugate-Kummer criterion](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/galois_horizon_and_kummer_covers/GENERIC_SYMMETRIC_MONODROMY_OF_WEIGHTED_MULTIQUADRATIC_SUMS_v1.0.md:1326>) and [incidence workflow](</home/williaml/Cella Framework/Papers_Library/05_expository_companions_and_research_maps/galois_horizon_and_kummer_covers/MACAULAY2_REALIZATION_POSET_WORKFLOW_2026-07-10.md>): keep the actual cover and its relations when checking a proposed branch quotient. The polynomial pair ideal adds dynamic equivalence rather than imposing every algebraic deck identification.
- [DIS chapters 10](</home/williaml/Cella Framework/research/derivations/DIS_chapters/10_bipolar_axes_and_constitutive_resistance.md>) and [12](</home/williaml/Cella Framework/research/derivations/DIS_chapters/12_internal_capacity_energy_and_final_hypotheses.md>): constitutive dependence, internal state/history and work/storage are native research inputs; their candidate equations do not select (14), (23) or (24).
- Existing PINCH-1/2/3, [constitutive divisors](2026-09-05-cella-constitutive-divisors.md) and [dependency reconciliation](2026-09-05-cella-dependency-reconciliation.md): the resolved chart, the response family and power-preserving projection. Their existing proofs are retained.

Standard chain-rule/ODE uniqueness, polynomial ideal finiteness, Gröbner reduction and linear observability tools are retained. No online search or new external source was used. The general quotient machinery is not claimed as new control theory; the native minimum, uniform reconstruction through the pinch, parameter counts and constitutive continuation are the project-specific results.

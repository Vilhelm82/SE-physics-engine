**Seated Root × Cella: a programme for mathematical advancement**

5 September 2026. Replacement assessment based on the active Lorentzian construction.

**The highest-value next task is to derive a seat/load response that determines both state evolution and propagation.** That would connect the model's geometry to quantities it can predict. The strongest Cella contributions are its recovery of information from multiple presentations, its valuation calculus, and its algebraic treatment of covers. Their useful extensions concern the active response and its degenerations.

The Euclidean cell's scheduled removal is accepted. Its cocycles, round metrics, paravector dynamics and six-dimensional lift do not determine this programme's priorities. The calculations accompanying this assessment start with the current seat-derived Lorentzian form. They establish several supporting constructions; the proposed dynamics remains a research target.

The [fresh symbolic runner](</home/williaml/seated-root/docs/archive/lorentzian-upgrade-derivations.py>) and its [results](</home/williaml/seated-root/docs/receipts/lorentzian-upgrade-derivations.json>) contain the calculations developed here.

**1. Derive one response law with enough structure to produce predictions**

The current foundations already supply the relevant starting point: P11/P12, the T7 seat construction, the T8 operation connecting the two Clifford blocks, and THM-K's charged-load target. The next substantial advance would be a response operator whose coefficients come from a specified interaction between the seat and load. Its low-frequency response, dissipation, ruler exchange and propagation could then be calculated together. [Current primitives and THM-K target](</home/williaml/seated-root/docs/results/2026-09-04/2026-09-04-PRIMITIVES-v0.md:89>)

**Cella source update, 5 September:** the [DIS material/environment response discussion](</home/williaml/Cella Framework/research/derivations/DIS_chapters/10_bipolar_axes_and_constitutive_resistance.md:337>) belongs to legitimate DIS physics research in the Cella corpus, as William has confirmed after checking the source paths. It supplies an existing material/environment response programme to develop into the native seat/load law. Its state variables, stored quantity and response coefficients should be derived together; the reduction criterion applies where internal variables are eliminated. The [source reconciliation and derived extensions](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-cella-dependency-reconciliation.md) also identify Cella's algebraic indefinite coupling form, signed curvature channels and cut-area trace decomposition as inputs to the horizon continuation.

A concrete construction is to retain the load's internal variables while deriving the coupled equations. If their linearization has the form

\[
\begin{pmatrix}A&B\\ C&L\end{pmatrix}
\begin{pmatrix}f\\y\end{pmatrix}
=
\begin{pmatrix}j\\0\end{pmatrix},
\]

then eliminating the load gives

\[
A_{\mathrm{eff}}=A-BL_{\mathrm{ret}}^{-1}C.
\]

Here \(f\) denotes the retained perturbations, \(y\) the load variables, and the inverse uses the causal boundary condition of the proposed dynamics. In the static reciprocal quadratic case, \(C=B^{T}\); completing the square gives the same Schur complement. The runner checks that elimination explicitly.

**The mathematical work is to derive these blocks from the interaction.** This construction makes that work more focused: the same blocks determine the measurable response, stored energy, dissipation and possible memory effects. A state-dependent impedance becomes a result of elimination. Additional load variables can be retained whenever eliminating them would conceal relevant history.

Keep the two ruler channels unnamed during this derivation. For their exchange \(S\), use the algebraic projectors \((1+S)/2\) and \((1-S)/2\) to extract the even and odd responses. This directly continues the existing impedance work. Coupling to the T8 block-exchange operation can then determine which response actually involves the two sheets. [Impedance construction](</home/williaml/seated-root/hunch_z0_impedance.py:1>), [fourth-generator construction](</home/williaml/seated-root/prim_t8a_fourth_direction.py:1>)

There is an exact preliminary calculation that sharpens the propagation target. At the orthogonal reference state, write tangent Gram coordinates as

\[
x=(\delta a,\delta b,\delta\gamma).
\]

Independent ruler-pole reversals act by
\(\operatorname{diag}(-1,1,-1)\) and
\(\operatorname{diag}(1,-1,-1)\); ruler exchange swaps the first two coordinates. A symmetric quadratic response invariant under these operations has precisely the form

\[
\operatorname{diag}(A,A,D).
\]

This follows by solving the full invariance equations, which the runner does. Thus a neutral response has two independent quadratic coefficients at this state.

For a **candidate** local, reciprocal, second-order theory in which all three coordinates are dynamical, with isotropic spatial coupling and positive kinetic tensor, write

\[
K=\operatorname{diag}(k_d,k_d,k_t),\qquad
C_0=\operatorname{diag}(c_d,c_d,c_t).
\]

Its principal polynomial is

\[
\det(-\omega^2K+|\mathbf k|^2C_0)
=(c_d|\mathbf k|^2-k_d\omega^2)^2
 (c_t|\mathbf k|^2-k_t\omega^2).
\]

The two speed classes are therefore

\[
v_d^2=c_d/k_d,\qquad v_t^2=c_t/k_t.
\]

A common cone in this sector is equivalent to \(c_d/k_d=c_t/k_t\). Agreement with the seat's cone additionally fixes the common ratio in seat-clock and ruler units.

**Stronger target:** derive the kinetic and spatial response from the same seat/load interaction, identify its propagating and constrained modes, and determine whether it forces this equality. The response may instead produce distinct modes with calculable speeds. Either outcome advances the model. The displayed second-order sector is a tractable probe; the full kernel can reveal first-order terms, memory, constraints or different symmetries under a charged load.

For the electromagnetic sector, derive the general constitutive operator and propagation equations from the native/Cella interaction, then calculate their characteristics before specializing to a Hodge operator. The existing \(Z^{-1}\star_q\) result supplies a comparison point. The previously suggested premetric reconstruction framework remains a possible external comparison, registered as **EXT-006** in the [dependency ledger](/home/williaml/seated-root/docs/EXTERNAL-MATHEMATICS-DEBT.md). Its reconstruction theorem has not been adopted in the pinch proof; establish what the internal response construction supplies before importing it. [Existing native constitutive runner](</home/williaml/seated-root/thm_k_response_map.py:1>), [original external source: Hehl](https://arxiv.org/abs/1601.00320)

The thermal programme can use the same interaction: derive both its retarded response and its state-dependent correlations, carry the native complex-pivot cycle through that response, and calculate the resulting spectrum on the dynamically obtained trajectory. This continues T4's period calculation toward emission and flux. The load's statistical state belongs in that construction. [Native Wick-period calculation](</home/williaml/seated-root/prim_t4_hawking_period.py:1>)

**2. Recover the hidden depths using controlled presentations**

This is the most immediate overlap with Cella's observable-invariant work. T7 constructs a visible angle and two depths invisible to the static presentation. A useful upgrade is an explicit recovery theorem under a specified family of probes. [T7 realization and static readout](</home/williaml/seated-root/prim_t7_seat_form.py:84>)

Use the current Gram matrix

\[
G=
\begin{pmatrix}
-1&a&b\\a&1&\gamma\\b&\gamma&1
\end{pmatrix},\qquad
\delta=-\det G
=1-\gamma^2+a^2+b^2-2ab\gamma.
\]

Let \(C=\cos t\), \(s=\sqrt{1-C^2}>0\), \(A=\sqrt{1+a^2}\), \(B=\sqrt{1+b^2}\). In a selected-pole chart,

\[
q=\operatorname{diag}(1,1,-1),\quad
c=(0,0,1),\quad
h=(A,0,-a),\quad g=(BC,Bs,-b).
\]

Then \(\gamma=ABC-ab\).

**Probe hypothesis.** Read the projected ruler angle in a fixed reference detector while applying two calibrated common Lorentzian pivots to the frame. Equivalently, keep the frame fixed and vary the detector normal. No ruler is moved independently. The relation of this reference-detector protocol to the definitive P6 implementation is part of the proposed operational construction; a detector carried along with the same common isometry would give different data.

For a unit negative detector normal \(n\), the observed cosine is

\[
\mathcal C(n)=
\frac{\gamma+q(h,n)q(g,n)}
{\sqrt{(1+q(h,n)^2)(1+q(g,n)^2)}}.
\]

Probe along
\(n_x(z)=(z,0,\sqrt{1+z^2})\) and
\(n_y(z)=(0,z,\sqrt{1+z^2})\).
At \(z=0\), the parameter has the same first derivative as rapidity. Direct differentiation gives

\[
m_x=s^2\frac bB,\qquad
m_y=s\left(\frac aA-C\frac bB\right).
\]

Consequently,

\[
v_b=\frac{m_x}{s^2},\quad
v_a=\frac{m_y}{s}+Cv_b,\qquad
a=\frac{v_a}{\sqrt{1-v_a^2}},\quad
b=\frac{v_b}{\sqrt{1-v_b^2}}.
\]

The angle and two slopes recover all three continuous Gram parameters. The Jacobian is

\[
\det\frac{\partial(C,m_x,m_y)}{\partial(a,b,C)}
=-\frac{(1-C^2)^{3/2}}
{(1+a^2)^{3/2}(1+b^2)^{3/2}},
\]

so recovery is locally nonsingular at finite depths for \(0<t<\pi\). The exact fixture
\((a,b,C)=(3/4,4/3,3/5)\) gives
\((m_x,m_y)=(64/125,12/125)\) and recovers both depths.

**Continuation target:** derive the admissible finite-pivot protocol from P6, then construct an atlas of recovery procedures including the branch strata. At \(s=0\), this first-derivative protocol loses rank; higher derivatives or another admissible protocol become concrete development targets. Pole conventions and orientation history should be carried explicitly through the atlas. The result above recovers continuous Gram data on its selected chart.

Cella's \(n=3\) RoleChSpec theorem is directly useful as a proof strategy: determine the actual observational equivalence, construct separating data, and exhibit a recovery minor. That theorem already proves completeness for its same-gradient Hessian quotient. The Seated Root application needs its own probe maps, as the calculation above illustrates. [RoleChSpec theorem and proof](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/dbp_role_channel_and_orbit_geometry/dbp_orbit_calculus.tex:381>)

**3. Give the dynamics variables that continue across the ruler-plane boundary**

Set \(\varepsilon=1-\gamma^2\). The native lapse satisfies \(N^2=\varepsilon/\delta\). A useful regular variable is the dual normal

\[
\ell=(1,0,0)^T,\qquad
\nu=G^{-1}\ell
=\frac1\delta
\begin{pmatrix}
-\varepsilon\\a-b\gamma\\b-a\gamma
\end{pmatrix}.
\]

It satisfies

\[
q(\nu,h)=q(\nu,g)=0,\quad q(\nu,c)=1,\qquad
q(\nu,\nu)=-\varepsilon/\delta.
\]

Thus \(\nu\) stays finite at \(\gamma=\pm1\) wherever \(\delta\ne0\), and becomes null there. At \((a,b,\gamma)=(-1,0,1)\), \(\det G=-1\) and \(\nu=(0,-1,1)\).

T7g already identifies the null ruler separation. This construction supplies a normalization that extends that direction through neighbouring states without dividing by the vanishing ruler-plane minor. **Use it to formulate the response and continuation equations across that boundary.** At full frame rank loss, \(\delta=0\), homogeneous adjugate data and a stratified description provide the next construction to develop. [T7g null separations](</home/williaml/seated-root/prim_t7g_two_horizons.py:1>)

A second supporting construction is an explicit positive motion cost on the native state variables. It makes a possible dynamics calculation concrete while exposing the cost assumption.

For the unit negative seat \(c\), define

\[
h_c(u,v)=q(u,v)+2q(u,c)q(v,c).
\]

Writing \(u=\alpha c+u_\perp\) gives
\(h_c(u,u)=\alpha^2+q(u_\perp,u_\perp)>0\) for \(u\ne0\).
**Assume equal ruler-velocity costs measured by this seated form**, fix the seat, and eliminate common azimuthal velocity by minimization. With depths \(l_1,l_2\) and relative angle \(t\), the resulting metric is

\[
g_{\rm kin}=
\cosh(2l_1)\,dl_1^2+
\cosh(2l_2)\,dl_2^2+
\frac{\cosh^2l_1\cosh^2l_2}
{\cosh^2l_1+\cosh^2l_2}\,dt^2.
\]

The angular coefficient follows from minimizing
\(C_1\Omega^2+C_2(\Omega+\dot t)^2\), with \(C_i=\cosh^2l_i\).

At \(l_1=\operatorname{arsinh}1,\ l_2=0,\ t=\pi/4\), the ruler plane reaches \(\gamma=1\), while this candidate metric is
\(\operatorname{diag}(3,1,2/3)\): finite and positive. At the orthogonal state it is \(\operatorname{diag}(1,1,1/2)\), giving a concrete candidate kinetic tensor for the propagation calculation above.

**Stronger target:** obtain the kinetic cost from the derived load response and compare it with this native geometric candidate. Agreement would connect the geometry to the interaction law. Any additional terms would specify how the load changes that law. These regular state variables support a future dynamical horizon calculation; spacetime continuation requires the resulting evolution and field reconstruction.

**4. Extend Cella's valuation calculus to full tensors and measured response**

The completed local-curvature manuscript already supplies directional curvature channels, arbitrary normal-crossing codimension, weighted initial forms and the hierarchy through higher-unit-jet cancellation. Its weighted valuation selects the first nonzero combined coefficient along the actual approach. That is substantial machinery available for continuation work. [Completed calculus, directional identity](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/local_curvature_and_black_hole_metrics/LOCAL_CURVATURE_CALCULUS_COMPLETE_v1.0.md:177>), [weighted valuations and cancellation](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/local_curvature_and_black_hole_metrics/LOCAL_CURVATURE_CALCULUS_COMPLETE_v1.0.md:789>)

**The strongest extension is a calculus for full Lorentzian tensors, moving frames and source-dependent response near degeneracy.** This lets the corpus address rotating/off-diagonal geometries, changes of resolution and singular propagation using the same local methods.

For geometry, begin with a coframe \(e^a\) and signature matrix \(\eta\):

\[
g=\eta_{ab}e^a e^b,\qquad
de^a+\omega^a{}_b\wedge e^b=0,\qquad
\Omega^a{}_b=d\omega^a{}_b+\omega^a{}_c\wedge\omega^c{}_b.
\]

The proposed theorem should track weighted expansions of the coframe and its exterior derivatives, determine curvature/transport orders, and prove covariance under the stated class of regular frame changes. This retains the twisting and mixing data needed beyond coordinate-diagonal metrics. Begin with analytic or formal expansions and then extend to the smooth/polyhomogeneous classes appropriate to the application.

For response, retain the source and readout maps:

\[
T(s)=B(s)L(s)^{-1}C(s).
\]

Over a one-parameter formal power-series ring, a generically invertible operator admits Smith factorization

\[
L=U\,\operatorname{diag}(s^{\alpha_i})\,V,
\]

with \(U,V\) invertible over that ring. Therefore

\[
T=B V^{-1}\operatorname{diag}(s^{-\alpha_i})U^{-1}C.
\]

This expression identifies exactly where Cella's coefficient-sensitive valuation machinery enters: include the transformed source and readout, combine terms of equal order, and retain the first nonzero coefficient.

For example, take \(L=\operatorname{diag}(s,s^2)\) and drive \(C=(1,1)^T\). Readouts \((0,1)\), \((1,0)\), and \((1,-s)\) give respectively \(s^{-2}\), \(s^{-1}\), and \(0\). The runner verifies these exact responses.

**Research deliverable:** an invariant classification that distinguishes operator degeneration, excited singular modes, observed singular modes and regular continuation. For black-hole work, feed it the operator derived in item 1 and the native variables in item 3. The extension would also apply to constrained mechanical systems and Cella's coupling diagnostics. Multiple divisor parameters require the corpus's weighted methods and suitable local resolutions; the one-parameter Smith form is the initial tractable case.

**5. Transfer the Kummer method into the active normalization fields**

The reusable Cella result is its control of square classes through divisor valuations. Its all-\(k\) closure manuscript already gives a uniform norm-derived second radicand and the generic closure theorem. For the present project, apply that method to the functions the Lorentzian construction actually uses. [All-\(k\) construction and closure](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/galois_horizon_and_kummer_covers/ALL_K_TWO_RADICAL_KUMMER_CLOSURE_v1.0.tex:58>)

Over \(K=\mathbb Q(a,b,\gamma)\), consider

\[
\delta,\quad \varepsilon=1-\gamma^2,\quad
1+a^2,\quad 1+b^2.
\]

The prime divisors
\(\delta,\ 1-\gamma,\ 1+a^2,\ 1+b^2\)
give the identity parity matrix for these four radicands. Irreducibility of \(\delta\), viewed as a quadratic in \(\gamma\), follows from its discriminant \(4(1+a^2)(1+b^2)\), which is not a square in \(\mathbb Q(a,b)\). Hence

\[
[K(\sqrt\delta,\sqrt\varepsilon,\sqrt{1+a^2},\sqrt{1+b^2}):K]=16,
\qquad \operatorname{Gal}\cong(C_2)^4.
\]

This is the generic cover of these normalization functions. Its observable characters already separate several jobs:

- \(N^2=\varepsilon/\delta\) and the dual normal \(\nu\) are rational.
- \(N=\sqrt\varepsilon/\sqrt\delta\) reads the product of two radical signs.
- The visible cosine \((\gamma+ab)/\sqrt{(1+a^2)(1+b^2)}\) reads another sign product.
- Signed frame volume involves \(\sqrt\delta\).

**Stronger target:** derive the minimal subcover required by the actual response and transport, then classify its real continuation across the native strata. Some positive normalizations have canonical real branches; transport and orientation can retain different information. The algebraic degree above counts the specified field extension, not physical states. This calculation does not identify a click group.

This gives the corpus's algebraic-cover methods an immediate Lorentzian application with explicit functions and observables.

**6. A framework that survives removal of the physics labels**

The useful common structure can be expressed through:

- a state space, including its rank-changing strata;
- admissible transformations and controlled probes;
- observations and the equivalence they induce between states;
- load variables and interactions;
- evolution and measurable response.

These objects can describe a seated frame, an implicit engineering constraint or a family of algebraic solutions. Particular metrics, spin representations and physical units enter where the application derives or declares them.

A precise generalization target is **closure of observable dynamics**. Let \(\Phi\) collect the chosen observations and let \(V_L\) be a proposed state evolution under a specified load. On a regular quotient, a reduced evolution \(W_L\) exists exactly when

\[
D\Phi_x V_L(x)=W_L(\Phi(x))
\]

is well-defined: observationally equivalent states must have the same observable velocity under corresponding loads. This is a direct way to develop Cella's coupling-evolution programme. If the chosen curvature or channel data do not close the dynamics, enlarge the observable set, its derivative history, or the retained load state until closure can be established. The recovery calculation in item 2 supplies one local route for Seated Root.

This suggests the following assessment of the mathematical tools:

- **Lie groups and relational Gram geometry:** retain them for the seat-derived local transformations and continuous state. Advance them through admissible probes, state-dependent transport and strata.
- **Clifford and spin structures:** retain them for the established sheet and adjoint algebra. Their next contribution is an actual interaction using the block-exchange operator.
- **Holonomy and complex continuation:** retain the native boost-composition and period calculations. Derive the transport selected by the evolving state and evaluate its action on the actual response and correlation functions.
- **Exterior and constitutive calculus:** use them to formulate field response and calculate its characteristics. Derive the constitutive coefficients from the load.
- **Differential invariants and role-channel methods:** use them for observable recovery and dynamical closure, with the actual equivalence relation of each application.
- **Valuations and singularity theory:** extend the existing calculus to full tensors, propagators and coupled observations.
- **Kummer and monodromy methods:** apply them to the native observable fields and transport families.

**An entirely new algebra is not presently the most productive first investment.** The strongest opportunities for original mathematics lie in theorems joining these structures: recovery across changing rank, minimal observable states with closed dynamics, and response/transport classification through singular strata. Each has a concrete first case here. Their originality would need comparison with the relevant literature after a precise theorem is obtained.

**7. Stronger targets for the ROI outlook**

For the goal of advancing Seated Root while growing the corpus, I would reorder and sharpen the outlook as follows. These are proposed research updates; the source outlook has not been edited. [ROI outlook](</home/williaml/Cella Framework/docs/ROI_OUTLOOK_2026-07-17.md>)

1. **Promote geometry selection, D8, into response-derived geometry.** Derive the seat/load operator, its kinetic response, its interaction between cells and its propagating modes. This also gives Cella's coupling-field programme a concrete evolution problem.
2. **Use extractability, X3, for native controlled observations.** Complete the Lorentzian probe protocol, prove recovery on regular charts, and extend across the branch strata. This supplies an immediate mathematical case while the wider observational programme continues.
3. **Advance D10/D15 through full tensors and observed singularities.** Build on the completed cancellation calculus; add coframe transport and source/readout dependence.
4. **Direct the Kummer continuation toward observable subcovers and real dynamics.** Use the active normalization calculation as the first case. Broader monodromy questions become directly applicable when the model supplies its wave or transport family.

The first substantial dynamics deliverable should specify one seat/load interaction and compute its effective response. The first coupled deliverable should use that interaction to derive propagation between two cells. With the propagating sector established, extend to a spatial family and solve a controlled loading problem through \(\gamma=\pm1\), tracking whether its path encounters \(\delta=0\).

That would produce the ingredients for a serious black-hole prediction calculation: a generated state trajectory, a propagation operator, a continuation rule and a defined observation. Scattering, mode frequencies and fluxes can then be computed from the model's own evolution.

The present runner verifies the displayed algebraic constructions and conditional calculations. The seat/load law, operational realization of the probe, and nonlinear field evolution remain the next research work.

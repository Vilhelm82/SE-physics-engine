# RH-1: restoration as the natural direction of motion

6 September 2026. User clarification: “rest is the rest state and drift from rest takes work so it wants to return.” This changes the build target from measured angle correction to a physically attracting control configuration. The earlier finite-count servo is a comparison model; it does not establish the requested restoring mechanism.

**Next user-directed construction:** [the trine as DSP on `F=0`](../../../../results/2026-09-06/2026-09-06-RH-1-cella-constraint-surface.md) now supplies an explicit angle-chart closure, Cella's three signed curvature channels, and a restoring constitutive law along the common-angle actuator. This directs further work toward the constraint/material coupling. The operator-square proposal below remains a separate confinement observation.

## Exact energy target

Let `s` be a physical orientation coordinate, `lambda` the environmental state, and `V(s,lambda)` an actual stored energy. For fixed `lambda`, a local rest configuration requires

\[
V_s(s_*,\lambda)=0,\qquad V_{ss}(s_*,\lambda)=\kappa>0.
\]

The native DIS/Cella constitutive construction already permits

\[
I\ddot s+\gamma\dot s+V_s=f_{\rm ext},\qquad I>0,\quad\gamma>0,
\]
\[
E=\tfrac12I\dot s^2+V,
\qquad
\dot E=f_{\rm ext}\dot s-\gamma\dot s^2+V_\lambda\dot\lambda.
\]

**[derived | declared storage and response]** At fixed environment and zero external forcing, the energy decreases by damping. Linearized about a strict minimum, `I q_ddot+gamma q_dot+kappa q=0` has roots with negative real part. Displacing the stationary system raises its potential by `kappa*q^2/2+O(q^3)`; that energy must come from an external source or a fluctuation. With sufficient damping it settles instead of ringing. In the overdamped limit, `q_dot=-(kappa/gamma)q` and `V_dot=-(kappa^2/gamma)q^2`. These statements follow from the existing constitutive energy balance, not from a numerical success criterion.

For changing environment, `V_lambda*lambda_dot` accounts for work entering as the minimum moves. A simple local overdamped tracker obeys `q_dot=-(kappa/gamma)q-s_star_dot`; sustained drift therefore leaves lag about `-gamma*s_star_dot/kappa`. Finite temperature, a driven auxiliary mode, and implementation losses require their own noise/power terms. Autonomous operation does not imply zero supplied power.

**Construction still required:** derive a physical coupling whose minimum `s_*(lambda)` coincides with the good gate orientation. Substituting the optimum from a fidelity calculation into `V=kappa*(s-s_opt)^2/2` merely assumes this construction. It is an illustrative stability model, not native autotuning.

## Preserve every logical state at rest

The desired rest is a control configuration supporting the whole encoded qubit. It must not select a single logical state. A sufficient ideal separation would be a logical factor and a restoring factor, with `H_rest=I_L tensor h_rest` and dissipative operators `L_mu=I_L tensor ell_mu`. The restoring factor can relax while the logical density matrix is unchanged by that relaxation. Actual gate couplings and bath back-action must be tested against this separation; factorization is a declared target here.

In a code-subspace formulation, requiring every operator `P rho P` to be stationary under the restoring generator tests populations **and coherences**. Merely checking that the two basis populations survive is insufficient. A unique attractive pure state would initialize the qubit rather than preserve arbitrary encoded information. Dark-space invariance also does not, by itself, prove correction of errors originating outside it.

## What the current native operator does and does not supply

**[native spectral fact]** REALFIBER fixes the spectrum of the embedded rank-two operator to `(-a,0,0,+a)`, with `a>0`. Its instantaneous dark plane is not a ground subspace. For a frozen operator and a bath allowing downward bright transitions, ordinary energy relaxation can favor the `-a` bright state. This observation is not a no-go theorem for engineered reservoirs or driven stabilization. The finite-duration transported code is also not generally the instantaneous kernel throughout a loop.

**[native algebraic candidate, not a control implementation]** The polynomial

\[
B=H_0^2/a^2
\]

is exactly the bright-space projector: `B^2=B`, `B>=0`, and `ker B=ker H_0`. Thus `H_pen=Delta B`, `Delta>0`, puts both bright branches above the dark plane. If added to the original frozen operator, `H_0+Delta B` has spectrum `(Delta-a,0,0,Delta+a)` and the dark plane is lowest for `Delta>a`.

This is a concrete candidate for an energy penalty derived from the native operator. It introduces diagonal terms and potentially other couplings outside the original zero-diagonal family, changes resource bounds and finite-return dynamics, and supplies no dissipative return channel by itself. Moreover, it follows the chosen `H_0(s)`: it does not choose the correct `s` under unknown detuning. It is confinement algebra, not an accomplished angle autotune or a correction theorem for logical errors within the dark plane.

**[exact nominal back-action check | treating s as a reciprocal physical coordinate]** The whole trine gives

\[
K(s)=R(2s)G,\qquad R(u)=e^{-iuY},\quad GYG=-Y.
\]

For a frozen commanded angle during the nominal word, Duhamel's identity gives

\[
iU^\dagger\partial_sU=\int U(t)^\dagger(\partial_sH(t))U(t)\,dt.
\]

If the physical generalized force is `-partial_s H`, the word's code-restricted force impulse, with `hbar=1`, is therefore

\[
\mathcal T=-iK^\dagger\partial_sK=2Y.
\]

Every nominal intermediate dump returns fully to the code, so no nominal loss branch contributes to this identity. The impulse on an input is `2*<Y>`, independent of `s`, rather than a restoring term proportional to `-(s-s_*)`. The two Y eigenstates push oppositely. For a maximally mixed input its mean is zero at every angle; that does not give positive stiffness. Thus simply allowing the existing reciprocal angle coordinate to move does not supply a universal centering force. This conclusion is restricted to this native drive and conjugate-force model; an additional reference, countercoupling or dissipative construction can change it.

**[bounded numerical cross-check]** Central differences of the existing exact nominal trine endpoints at `s=0,+0.1,-0.1`, step `1e-6`, give `||T-2Y||` respectively `7.83e-10`, `3.64e-10`, `7.31e-11`; impulse eigenvalues are `(-2,+2)` within `6e-10`. Native Hamiltonians at `t=0.37*TAU` on all three arcs have spectrum `(-1,0,0,+1)` within `5e-16`. These checks corroborate algebra; they do not simulate a physical rotor or reservoir. Reproduce using `K(s)=compose_discard([trine_loop(float(a+s)) for a in TRINE])[0]` and `T=-1j*K(s).conj().T@(K(s+h)-K(s-h))/(2*h)` from `rlq/rh1_common.py`.

## Next construction and falsifiers

The first physical candidate should couple a reference/restoring coordinate to the native angle, with explicit storage and a dissipative outlet. Derive its force, fixed point and noise from that joint model. Then require:

- A displacement raises the declared physical storage, and the released system returns with a positive local restoring rate.
- The fixed point follows the low-error operating geometry when the actual environment changes, without a hidden optimum or true-error input.
- The force and damping preserve unknown logical data; trace preservation, logical coherences, leakage and erasure are retained in the test.
- The return rate exceeds the declared drift rate at counted drive power, bath noise and control cost. The moving-angle connection term and finite-word cancellation must be included if the frame moves during a gate.

The current rank-one measured dump cannot supply full-qubit restoration; RH-1's obstruction remains. A new coherent or dissipative instrument requires its own correction analysis. No artificial spring simulation is counted as evidence that such an instrument has been built.

## Provenance

Native inputs: canonical REALFIBER O1–O7, the real reference rotation and nominal trine reflection, the [constitutive extension](../../../../results/2026-09-05/2026-09-05-cella-constitutive-divisors.md), section 4 equations (17)–(19), and its DIS chapter 10 material-response and chapter 12 stored-energy backing passages. Standard passivity and finite-dimensional operator calculus are retained through EXT-004/022.

External comparison only: [Mirrahimi et al., Dynamically protected cat-qubits](https://arxiv.org/abs/1312.2017), primary PDF abstract and introductory encoding discussion; [Gertler et al., Protecting a Bosonic Qubit with Autonomous Quantum Error Correction](https://arxiv.org/abs/2004.09322), primary abstract. They establish background examples of engineered dissipative protection of encoded information; neither construction or performance result is transplanted into this model. See EXT-024 and SR-012.

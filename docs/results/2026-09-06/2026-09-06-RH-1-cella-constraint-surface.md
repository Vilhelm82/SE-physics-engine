# RH-1: the trine as DSP on a Cella constraint surface

6 September 2026. User direction: look up the Cella three-channel decomposition, take `F=0` as the constraint surface, and identify the trine with DSP. This note applies that direction in a concrete **ordered-angle realization**, `(D,S,P)=(alpha_1,alpha_2,alpha_3)`. This scalar chart is an explicit modeling choice, not a claim that the words Directive/Substrate/Product uniquely select three control angles in every physical realization.

The three ordered settings define the DSP coordinates. The three curvature channels describe the **whole coupling surface**; they are not one curvature channel per loop. The metric used for the geometric diagnostic is Euclidean in these equally scaled angular coordinates. `H_F` below is the Hessian of the constraint, not the driven quantum Hamiltonian.

## 1. Native Cella object recovered

The canonical [Three-Channel Strong Spec](</home/williaml/Cella Framework/Papers_Library/05_expository_companions_and_research_maps/dbp_role_channel_and_orbit_geometry/Three_Channel_KG_Strong_Spec.md>), sections 2–6, defines at a regular point of `F=0`

\[
g=\nabla F\ne0,\quad q=g^Tg,\quad H_F=\nabla^2F=H_c+H_s,
\quad H_s=\operatorname{diag}(H_F),\quad H_c=H_F-H_s.
\]

Its exact signed account is

\[
K_G=\kappa_c+\kappa_s+\kappa_{\rm int}
=-\frac{\det\begin{pmatrix}0&g^T\\g&H_F\end{pmatrix}}{q^2}.
\]

Pure off-diagonal Hessian terms give `kappa_c`, pure diagonal terms give `kappa_s`, and their mixed terms give `kappa_int`. The [Canonical Invariant Reduction Theorem](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/dbp_role_channel_and_orbit_geometry/Canonical_Invariant_Reduction_Theorem.md>) proves that changes of defining-function gauge can redistribute the channel account while preserving its sum. The physical DSP chart and defining function must therefore be recorded. Relabelling axes is not the same operation as physically changing the loop controls.

## 2. Exact nominal constraint from the native word

Let `R(u)=exp(-i u Y)` and `G=diag(-1,1)` on the logical code. A loop with angle `u` gives

\[
G(u)=R(u)G R(u)^T=R(2u)G.
\]

The three chronological loop settings obey

\[
G(P)G(S)G(D)=R\bigl(2(D-S+P)\bigr)G.
\]

Thus the local branch for the fixed target `G` has the closure relation

\[
\boxed{F_0(D,S,P)=D-S+P=0.}
\]

The existing trine is `(pi/3,2pi/3,pi/3)`, which lies on this plane. Equivalent distant branches up to global gate phase are not part of this local chart. It has

\[
g_0=(1,-1,1),\quad H_{F_0}=0,\quad
(\kappa_c,\kappa_s,\kappa_{\rm int})=(0,0,0).
\]

These zeros are exact. A flat constraint surface can support attracting normal dynamics; its Gaussian curvature is not a potential-energy stiffness. The plane also contains non-trine settings, so we retain the admissible common-angle family separately:

\[
x=x_0+s a,\quad x=(D,S,P),\quad a=(1,1,1).
\]

This preserves all relative loop angles. Since `a dot g_0=1`, common-angle motion is transverse to the closure surface. It is not the Euclidean normal, but it has the normal component needed to restore closure. Allowing an unrestricted gradient correction in all three coordinates would change the relative geometry.

## 3. Constitutive restoration on the constraint

Use the existing Cella passive-response framework with the following **declared material law**:

\[
E(x,\lambda)=\tfrac12 kF(x,\lambda)^2,\quad k>0,
\qquad M=\eta aa^T,\quad\eta>0,
\qquad \dot x=-M\nabla_x E.
\]

Here `lambda` collects environmental parameters. This specifies a restoring energy and a rank-one mobility along the existing actuator; it does not infer an energy from a curvature sign. At fixed environment,

\[
\boxed{\dot F=-\eta k(a\cdot g)^2 F,\qquad
\dot E=-\eta k^2F^2(a\cdot g)^2\le0.}
\]

The nominal word gives `F_dot=-eta*k*F`. A common displacement `s` from rest stores `k*s^2/2`, and its released motion returns exponentially in this overdamped response model. At `F=0`, `Hess(E)=k*g*g^T`: positive normal stiffness and zero tangent stiffness. Along the actuator the rest stiffness is `k*(a dot g)^2`. Every relative angle stays fixed because `x_dot` is parallel to `a`.

More generally a declared tangent flow `b_T`, satisfying `g dot b_T=0`, can be added without changing the fixed-environment storage identity. Under external force or changing environment, their power must be retained. In particular,

\[
\dot F=-\eta k(a\cdot g)^2F+F_\lambda\dot\lambda.
\]

Where `|a dot g|` is bounded away from zero, the same implicit surface has a locally unique intersection with the actuator family, and a locally positive restoring rate. If `a dot g=0`, this actuator loses linear restoring authority even if `g` is nonzero. This supplies a concrete geometric falsifier. Rest is invariant for the specified unforced law; we are not imposing an impassable barrier against every external input.

**Keep the two objects separate:** compute curvature from the regular signed `F`, and storage from `E`. Using `F^2` or gate infidelity as the defining function can make the gradient vanish at rest, where the regular Cella curvature formula would be undefined. Squaring the stored mismatch does not license squaring the curvature diagnostic.

## 4. A perturbed surface from the actual native propagator

To continue the nominal signed closure under noise, compute the full surviving map `K(D,S,P;lambda)`, its unitary polar `V`, and `W=VG`. On the local nonsingular branch with nonzero trace, use

\[
\boxed{F=\tfrac12\arctan\operatorname{Re}
\left(\frac{i\operatorname{tr}(YW)}{\operatorname{tr}W}\right).}
\]

This ratio removes the global phase. At zero noise it reduces locally to `D-S+P`. With noise it defines closure of the residual Y component. It does not assert that the other logical error components or erasure vanish. This is an offline, model-derived defining function; its evaluation from the propagator is not an implemented physical force or a freely available controller observation.

The [runner](../../../rh1_cella_surface.py) uses the existing Cella `three_channel_referee` from `engine/src/cella/reference_lift.py`, rather than a substitute curvature scalar. The nominal jet is derived symbolically and evaluated exactly. Perturbed jets are central finite differences at steps `0.002,0.001,0.0005`; exact arithmetic on their decimal values does **not** make these numerical jets exact. Their supplied-jet determinant partitions agree exactly, and medium/fine channel values agree within 0.1%.

At gain error `0.01`, the common-angle intersections and fine-grid channel values are:

- d-site detuning `+1e-4`: `s_*=-9.54496883e-5`, `a dot g=0.999963683`; `(kappa_c,kappa_s,kappa_int) ≈ (5.8877e-13,5.4973e-9,-1.4553e-10)` and `K_G≈5.3524e-9`.
- d-site detuning `-1e-4`: `s_*=+9.61462273e-5`, `a dot g=1.000037448`; channels approximately `(6.8653e-13,5.3717e-9,+1.5206e-10)` and `K_G≈5.5244e-9`.
- d-site detuning `+1e-3`: `s_*=-9.20031516e-4`, `a dot g=0.999687000`; channels approximately `(7.7472e-9,6.1456e-7,-1.6704e-7)` and `K_G≈4.5527e-7`.

The interaction term changes sign across the first two detuning cases and partially offsets self curvature at positive detuning. These are local channel diagnostics in the declared gauge, not coordinate-free statements about which physical component dissipates power. The common actuator stays transverse in all three cases. The root solver is an offline surface diagnostic; its access to the actual detuning is not passed off as autonomous tuning.

Run `python experiments/2026-09-06/rh1_cella_surface.py --json docs/receipts/rh1-cella-surface.json`. The [receipt](../../../../receipts/rh1-cella-surface.json) retains all derivative levels, full conditional/erasure metrics, exact nominal values, source hashes and scope. The native keystone `(-1/49,+1/49,-3/49)` also passes as a signed-channel cross-check. This is a local research consumer of the Cella referee, not a new substrate-promotion claim.

## 5. What this closes and what remains to build

The ordered-angle DSP realization now has a native signed closure surface, the full three-channel diagnostic, an admissible restoring mobility, and a proved energy-decay identity conditional on the declared constitutive law. It preserves the relative trine geometry and expresses the desired work asymmetry at the control level.

The physical remaining task is to realize the normal effort `-kF*g` and its dissipation through a reference/material coupling that preserves unknown logical information. That implementation must include noise, supplied power and the moving-frame term if the angle changes within a word. The earlier `2Y` reciprocal-back-action observation still applies to the bare drive; the new constitutive law is not silently attributed to that bare Hamiltonian. The operator-square confinement proposal is not needed to establish this constraint-surface model.

Native source ledger: DAG revision `575f1e672b44c5303a308fa04fef07fa46d61f1f186a65fad1396507d850871b`; backing files for `U-0039` and `U-0600` read through the DAG source tool with matching live hashes; Three-Channel Extension sections 1–3; Curvature Orbit Correction sections 1–3; RE-2's actual trine product and active lines; BR-8's storage/face conditions; the already inspected DIS material/storage passages. No online source or new external mathematical framework was used. See EXT-025 and the internal ledger receipt.

## 6. Recommended simplification — one restoring coordinate

Subsequent recommendation, 6 September 2026; no physical implementation is claimed. Since the allowed motion already satisfies `x=x_0+s*a`, pull back the signed constraint exactly:

\[
f(s,\lambda)=F(x_0+sa,\lambda),\quad
f_s=a\cdot\nabla F,\quad
U(s,\lambda)=\tfrac12 kf^2.
\]

The same declared rank-one constitutive law is exactly

\[
\boxed{\dot s=-\eta U_s=-\eta kf f_s,\qquad
\dot U=-\eta(U_s)^2\le0}
\]

at fixed environment. No three-dimensional Hessian or curvature channel is needed to evolve this one-coordinate law. The existing noisy receipts give `f_s` from `0.9996870` to `1.0000374` at the three tested rest points; this is local evidence, not a global capture-range bound. Retain Cella's channel account for diagnosis if a later failure requires it.

The recommended first physical candidate is the fixed-relative-angle trine with a single restoring reference coupling and damping. Its force must respond to the relevant mismatch, preserve logical data, and make both signs of a small displacement relax toward correct operation. A spring tied only to the original command rejects command-angle disturbances but does not follow detuning that moves the correct operating angle. A simulator-provided optimum or full propagator is not an available physical reference. Return rate, drift, back-action and power are the decisive tests; the missing coupling has not been replaced by this exact reduction.

This reuses the standard dissipative energy balance already retained in EXT-004; its primary source was re-read for this recommendation. See SR-013. The earlier claim of no online research applies to the original construction above.

**First approved trial completed:** [a passive two-level reference and spring](../../../../results/2026-09-06/2026-09-06-RH-1-passive-reference.md) supplies a new specified coupling whose rest matches this constraint to first order. It does not realize the exact `F^2` well. The small reference fails the noise budget, a larger reference improves the settled conditional channel, and drift lag/reciprocal logical back-action remain explicit. Hardware realization and exact data preservation are open.

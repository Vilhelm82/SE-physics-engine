# RH-1: a physical reference for the one-angle trine

6 September 2026. First bounded trial of the approved [simplification](2026-09-06-RH-1-cella-constraint-surface.md#6-recommended-simplification--one-restoring-coordinate). [Plan](superpowers/plans/2026-09-06-rh1-passive-reference.md), [runner](../../../rh1_passive_reference.py), [receipt](../../../../receipts/rh1-passive-reference.json).

**Verdict:** a specified equilibrium reference coupling produces the desired *local* moving rest point without an online angle optimum or measurement of the logical data. At the declared temperature, the small reference fails because its angle fluctuates too much. A larger reference improves the calculated settled conditional channel, but tracking lag can temporarily make the gate worse. Exact `F=0` locking, exact preservation of every logical state, and a hardware/resource advantage are not established.

This is a new, conditional Hamiltonian-level candidate, not a derived component of the original zero-diagonal control family. Its parameters were fixed before the numerical trials in the linked plan. The original designer predictions remain unchanged; this continuation is exploratory and unblinded.

![Model tracking lag and reference-size noise budget](../../../../figures/rh1-passive-reference.png)

## 1. What supplies the force

Keep the relative trine angles fixed: `alpha_j=TRINE_j+s`. Add an independent reference ensemble of N two-level systems and a torsional spring. In units `hbar=1`, take

\[
H_{\rm ref}(s,\delta_r)
=\sum_{j=1}^N\frac{\Omega X_j+(\delta_r-gs)Z_j}{2},
\qquad H_{\rm spring}=\frac{N k_0s^2}{2}.
\tag{1}
\]

Here `Omega>0` is a static transverse splitting; the angle changes the longitudinal splitting with sensitivity `-g`. The actual environmental field enters the reference as `delta_r`, and enters the native gate as `diag(0,0,delta_d,0)`. The first trial declares `delta_r=delta_d`. It is a shared physical input, not a supplied estimate of the gate error. A static field gradient coupled to a torsional coordinate is a possible interpretation of the cross term. Engineering that gradient and having the same coordinate realize the native drive rotation remain implementation requirements.

The reference operators act as identity on the logical data. The *gate's* dependence on the same coordinate is reciprocal and remains relevant in section 4. Independence of the reference operator is not a proof that the complete data/reference interaction factorizes.

Assume the reference thermalizes rapidly at temperature T while s moves slowly. With `z=delta_r-gs`, `E=sqrt(Omega^2+z^2)`, the two reference eigenvalues are `+/-E/2`. The partition function gives the isothermal potential

\[
V(s,\delta_r)=N\left[\frac{k_0s^2}{2}
-T\log\left(2\cosh\frac{E}{2T}\right)\right],
\tag{2}
\]
\[
-V_s=-N\left[k_0s+\frac{gz}{2E}\tanh\frac{E}{2T}\right].
\tag{3}
\]

At zero temperature the reference contribution is its ground energy `-NE/2`. Equation (3) also follows directly from the expectation of `-partial_s H_ref` in the Gibbs state; it is not a force proportional to a simulator's F.

**[derived | (1), Gibbs equilibrium, slow reference elimination]** For every fixed field, the potential is globally strictly convex if

\[
k_0>\frac{g^2}{2\Omega}.
\tag{4}
\]

To check this, differentiate `z tanh(E/2T)/E`. The derivative is nonnegative and at most `1/Omega`: use `x sech^2(x)<=tanh(x)` for `x>=0`. Thus `V_ss>=N(k0-g^2/(2 Omega))>0`. The spring dominates at large |s|, so there is one rest point and the force points toward it on both sides. The spin contribution softens the spring; approaching the zero-stiffness boundary increases susceptibility and fluctuations together.

For overdamped mean motion `gamma s_dot=-V_s`,

\[
\frac{dV}{dt}=-\gamma\dot s^2+V_{\delta_r}\dot\delta_r,
\qquad
V_{\delta_r}=-\frac{Nz}{2E}\tanh\frac{E}{2T}.
\tag{5}
\]

This is a **free-energy** balance in a fixed-temperature reduction. It is not a complete microscopic internal-energy/heat balance. Moving the environmental field supplies or removes work. At finite temperature individual trajectories fluctuate; monotonic relaxation describes the deterministic mean equation, not every realization.

## 2. Match the native constraint once

Set `chi=tanh(Omega/(2T))/(2 Omega)`. Near zero,

\[
V_s/N=(k_0-g^2\chi)s+g\chi\delta_r+O((|s|+|\delta_r|)^3).
\tag{6}
\]

The native signed constraint has `f(s,delta_d)=s+c delta_d+O((|s|+|delta_d|)^2)` at nominal gain. An offline central-difference calibration gives `c=0.956495758259`; halving the step changes it by `2.50e-10` relative. Choose the constants

\[
\Omega=g=1,\quad T=0.01,\quad
k_0=g^2\chi+g\chi/c=1.022741471337.
\tag{7}
\]

Then `s_rest=-c delta_r+O(delta_r^3)` and `V_ss/N=0.522741471337` at zero. The calibration is frozen at gain zero. Later tests include gain error `0.01`; no new coefficient is fitted.

This realizes the first-order location of the earlier proposed `F^2` well. It does **not** realize that well exactly. The physical reference's rest point is odd in detuning, whereas the noisy native optimum already has even-order corrections. The actual nonlinear reference force is retained throughout the trial.

At gain error `0.01`, before angle fluctuations:

| Detuning | Physical rest angle | Conditional error at rest | Fixed zero-angle error |
|---:|---:|---:|---:|
| `-1e-4` | `+9.5649572e-5` | `1.385e-12` | `2.465e-8` |
| `+1e-4` | `-9.5649572e-5` | `8.457e-13` | `2.429e-8` |
| `+1e-3` | `-9.5649218e-4` | `4.337e-9` | `2.256e-6` |

The larger-detuning case exposes the nonlinear mismatch; it is not repaired with a fitted term. Giving the reference a 5% sensitivity error increases the two small-detuning errors to `4.97e-11` and `6.69e-11` before thermal fluctuations. Drifts that do not couple to the reference, including an independent second native detuning, are outside this one-field construction.

## 3. Rest has a noise and response budget

A thermal harmonic well of stiffness `K=V_ss` has local angle variance `sigma_s^2=T/K`. This follows from its canonical Gaussian distribution. An equilibrium overdamped Markov bath would require `ds=-(V_s/gamma)dt+sqrt(2T/gamma)dW`; its stationary Fokker–Planck density is proportional to `exp(-V/T)`. The runner uses this equilibrium distribution for gate averages; it does not simulate a microscopic bath or claim its white-noise approximation at all frequencies.

For the exact nominal trine, `K_gate(s)=R(2s)G`. Averaging a Gaussian angle with mean m gives

\[
\epsilon_{\rm angle}=\frac{1-e^{-8\sigma_s^2}\cos(4m)}3
\simeq\frac83(m^2+\sigma_s^2).
\tag{8}
\]

The relative channel contracts the X and Z Bloch components by `exp(-8 sigma_s^2)` at zero mean; it preserves Y. Thus zero mean error does not preserve an arbitrary unknown qubit exactly. This statement also applies to a frozen quantum coordinate initially uncorrelated with the data: tracing that coordinate averages over its position distribution.

The actual noisy-gate calculation averages the **surviving CP maps and every leakage map before conditioning**, including logical coherences and input-dependent survival. It does not insert the mean angle into a fidelity function. Gaussian quadrature uses 13 nodes, checked against 21; conditional fidelity uses the existing equal-input Haar quadrature.

At `delta_d=+1e-4`, gain error `0.01`, temperature `0.01`:

| Reference spins N | RMS angle | Conditional error with thermal angle | Assessment against fixed zero angle |
|---:|---:|---:|---|
| `1e6` | `1.3831e-4` | `5.1010e-8` | Worse |
| `1e8` | `1.3831e-5` | `5.1094e-10` | Improved |
| `1e10` | `1.3831e-6` | `5.9467e-12` | Improved; deterministic errors becoming relevant |

For the three settled fields `(+1e-4,-1e-4,+1.2e-4)`, N=`1e8` gives mean conditional error `5.1133e-10`. A constant angle given all three future fields and set to their linearized mean correction gives `2.4125e-8`: a **47.18-fold** improvement in this settled, initially uncorrelated model. The baseline is a strong fixed comparison, not a global optimizer certificate. N is a reusable ensemble size, not a shot count; this is not an equal-cost comparison with the previous measured servo.

Erasure remains approximately `0.001591–0.001593` per word. The unconditional metric remains dominated by erasure. This construction improves orientation, not the existing gain-leakage instrument.

For N=`1e8`, displacing the angle by `1e-4` costs locally `0.2614` energy units, or `26.14 k_B T`. For N=`1e6` the same displacement costs only `0.2614 k_B T`; thermal wandering is unsurprising. Energy units are `hbar omega_g` and time units `1/omega_g`, where the native coupling amplitude is normalized to one. No device temperature, inertia, field-gradient capability or refrigeration power is inferred from these dimensionless choices.

## 4. Back-action and a moving frame

The earlier native calculation gives the data-dependent force impulse `2Y` per nominal word when s is reciprocal. It remains present. For one fresh data word every `T_word=9 TAU=54.7530121`, a worst-sign *mean* load has magnitude at most `2/T_word`. In the harmonic static approximation its displacement is bounded by `2/(K T_word)`: about `6.99e-10` radians for N=`1e8` here. This is a nominal averaged-load estimate, not cancellation of quantum back-action, a bound on all noisy transients, or a multiword coherence theorem. Unknown data can become correlated with the reference and its bath.

The deterministic response time is declared as `gamma/K=100` words. A displacement falls to 1% in about 461 words; the numerical value after 500 words is `0.00673795` of its initial size. Reference equilibration and momentum relaxation must be faster than this for the reduction to apply. The chosen rate is a parameter, not a measured hardware bandwidth.

To test motion rather than silently hold the angle, the runner integrates the native **laboratory-frame** Hamiltonian with `s(t)` and `delta_d(t)` throughout each arc. This includes the effect represented by `-s_dot Y_P` in a co-moving frame. It retains all three intermediate dump maps. Frozen-angle recovery and time-step convergence are checked.

The field has 1000-word plateaus joined by smooth 100-word ramps. During the first sign reversal:

- At word 1050, the field crosses zero while the reference still has angle `-7.99e-5`. Conditional error is `1.752e-8`, versus `7.92e-11` for the fixed-zero angle during that same moving-field word.
- At word 1100, after reaching `-1e-4`, error is `3.768e-8`, versus `2.465e-8` fixed.
- At word 1500 it has fallen to `1.966e-11`, before thermal fluctuations.

Thus passivity does not give instantaneous tracking or a pointwise guarantee of improvement during drift. In the linear limit the displacement from the moving rest obeys `q_dot=-r q-s_rest_dot`; a faster moving target causes lag. The thermal-channel and moving-gate trials are separate controlled approximations: the former assumes a nearly frozen Gaussian angle, the latter integrates the deterministic mean motion. A stochastic joint multiword simulation is still required before combining them into a full running-device claim.

For the 3000-word deterministic trace with N=`1e8`, excess mechanical friction dissipates `1.90023` energy units, and the environmental free-energy work is `1.44589`; their difference equals the free-energy decrease, including the initially displaced state. The numerical balance defect per reference is below `2.7e-14`. These numbers exclude refrigeration, bias-field generation and the full microscopic reference-bath heat. The ideal static coupling has no pump term; that does not establish zero implementation cost.

## 5. Decision and remaining construction

Retain this as a **local passive-reference candidate**. It supplies a concrete joint energy and a restoring sign, and has a falsifiable size/temperature/bandwidth budget. Reject the stronger claims that the bare trine self-centers, that finite-temperature rest cannot fluctuate, that the physical rest equals the complete native `F=0` surface, or that a finite reciprocal reference exactly preserves arbitrary logical information.

The next useful decision is a physical platform and transducer: what actual field causes the relevant gate drift, how the same field shifts the reference splitting, and how its displacement rotates the native control. That choice must supply the available stiffness, temperature, inertia, reference relaxation and logical bath coupling. Until then increasing N in the model is a resource sensitivity calculation, not a completed build recommendation.

## Reproduction and provenance

Run `python3 -m unittest test_rh1_passive_reference -q`, then `python3 experiments/2026-09-06/rh1_passive_reference.py --json docs/receipts/rh1-passive-reference.json --plot docs/figures/rh1-passive-reference.png`. Four independent checks cover Gibbs force, potential derivatives, two-sided restoring stiffness, the non-Gaussian canonical variance limit and logical coherence loss; the runner records 24 further calibration, channel, energy and integration checks. Passing checks validate the stated models and numerical calculations, not their physical realization.

An independent read-only review found no actionable force, channel, integration or claim-boundary error. It also compared the implemented nominal thermal channel with equation (8), agreeing below `4e-17` in its check. The retained approximation boundaries above remain part of the result.

Native inputs: the actual RH-1 propagator/instrument, REALFIBER and trine backing proofs already checked for EXT-022/025; the one-angle pullback; [BR-8](2026-09-05-response-on-branched-resolved-state-spaces.md#8-storage-constraints-and-real-crossing); and DIS chapters 10/12 on material response, internal state and storage. These support conditional constitutive evolution but do not select (1).

Retained external background: [de Voogd, Wagenaar and Oosterkamp (2017)](https://www.nature.com/articles/srep42239), Basic principles, Susceptibility and generalization to two-level systems. Their reciprocal resonator/reference treatment motivates retaining force, relaxation and back-action together. Equations (1)–(8) and the native matching calculation are derived here for the declared candidate; no experimental parameters or performance numbers are transferred. The [Glauber-spin mechanical model](https://arxiv.org/abs/1004.3215) was read at abstract level as a comparison only. Full search and assumption accounting is EXT-026/SR-014 in the [dependency ledger](EXTERNAL-MATHEMATICS-DEBT.md).

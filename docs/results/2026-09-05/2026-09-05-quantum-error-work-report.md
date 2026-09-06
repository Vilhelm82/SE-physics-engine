# Tonight's quantum error control and general response theory

5 September 2026 · Seated Root · Results through commit `df70de0`, 23:03 AEST

**The main result is a five-loop heralded reflection gate that includes finite pulse durations, finite absorbing dumps and first-order logical detuning correction.** Its final error budget counts both recorded erasures and unobserved background absorption. A separate coherent construction cancels the complete second-order encoded gain/reference response within the earlier echo's duration budget.

These are derived control constructions and simulated operating points for the declared quantum model. The strongest immediate continuation is to optimize the complete loss instrument, then correct its remaining quadratic response.

**The other central result is a general theory of response on branched, resolved state spaces.** It specifies which responses remain regular, which state distinctions affect prediction, and when eliminating internal variables preserves the response. It supplies a mathematical framework usable beyond the original physics interpretation.

## 1. General response theory

The construction separates the retained state, its observable presentation, effort/flow ports, and the law connecting flows to state evolution. It covers supplied branched covers and smooth resolutions, including continuous directions retained over a collapsed observable state. Passivity is a specified sector of the theory.

Its principal results are:

- **Complete maximal linear passive relations in finite dimensions.** Ideal constraints and changes in the admissible effort space are included. Power-preserving port transformations compose exactly, and passive network elimination remains defined when a matrix inverse ceases to exist.
- **An exact smooth-lift criterion.** For a generically invertible square observation Jacobian `A`, put `Δ=det A` and `C=adj A`. A response matrix `M` and drift `c`, pulled back to the retained chart, lift smoothly exactly when
  \[
  CMC^T\in(\Delta^2),\qquad Cc\in(\Delta),
  \]
  entry by entry in the smooth-function ring. The unique lifts are the corresponding quotients. This tests the complete matrix, including cancellations between its entries.
- **Complete boundary divisor laws.** Prescribed zero-output faces force different vanishing factors in dissipative and skew response. The mixed-boundary extension includes both two-sided interfaces and one-sided boundaries. Whether a zero-output face prevents a trajectory crossing is decided separately by the actual state-rate law.
- **Singular internal-port elimination.** For a passive block response `[[A,B],[C,D]]`, eliminating zero-flow internal ports gives a full external graph exactly when `im C ⊆ im D`, equivalently `B ker D = 0`. Otherwise the external object is a constrained relation. Smooth continuation requires a smooth internal solution; eliminating first and then taking a singular limit can give a different answer.
- **Predictive descent and branch symmetry.** States may be identified only when the declared common future inputs give identical outputs and the reduced evolution is well posed. Under a locally Lipschitz equivariant drive, a trajectory outside a symmetry-fixed stratum cannot reach it in finite time. This supplies a concrete condition on branch crossing.

For the quantum programme, these results identify where to retain slow or hidden state, how to treat a singular probe map, and which response zeros survive allowed changes. The quantum realization remains supplied by a Hamiltonian or completely positive channel; classical passivity alone does not establish quantum-channel admissibility. The gate results below use their explicit unitary and absorbing evolution.

Source: [general theory, BR-1–9](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-response-on-branched-resolved-state-spaces.md), with **163 recorded checks** in [its receipt](/home/williaml/seated-root/docs/receipts/branched-resolved-response-checks.json).

## 2. Exact returns supplied the primitive

The starting point is Cella's real, symmetric, zero-diagonal Hamiltonian family with spectrum `(a, 0, …, 0, −a)`. Its real dark-space holonomy was extended to the full group `O(m−2)`, including an explicit orientation-reversing loop. In the four-level realization, the logical code is `P = span(d,r)`, the two remaining states form the bright sector `Q`, and the target gate is

\[
G=\operatorname{diag}(-1,1)|_P.
\]

The three-arc loop has **exact finite-time returns** at

\[
\tau_n=\frac{\pi}{2a}\sqrt{16n^2-1},\qquad T_n=3\tau_n,
\qquad n=1,2,\ldots
\]

At `n=1`, the bare loop takes `aT = 18.2510`. The transported code executes a nonadiabatic holonomy and returns with zero nominal bright population. It can occupy the bright sector during transport; this distinction makes an endpoint dump useful.

The active logical direction's integrated bright exposure is also exact:

\[
X_n=3\tau_n\left(2v_n^2-\frac32v_n^4\right),
\qquad v_n=\frac1{4n},\qquad
X_n\sim\frac{3\pi}{4an}.
\]

This establishes two independent control quantities: endpoint return error and accumulated exposure. Finite ramps, perturbation integrals and independent laboratory-frame evolution were subsequently incorporated into the runners.

Sources: [quantum construction](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-quantum-prediction-targets.md), [exact returns and perturbation response](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-reflection-loop-finite-duration.md).

## 3. Coherent correction: from gain cancellation to certified complete code correction

The first composite word used forward and inverse loops with finite holds to cancel the leading common-gain leakage. At 1% gain error, the independently reproduced average infidelity fell from `5.30776 × 10⁻⁴` for the bare loop to `2.49482 × 10⁻⁷` for the composite.

Adding physical access to the reference direction produced a finite echo correcting the first-order response to all six static generators in

\[
\operatorname{Herm}(Q)\oplus\mathbb R D_d\oplus\mathbb R D_r,
\qquad D_d=|d\rangle\langle d|,\quad D_r=|r\rangle\langle r|.
\]

It also removed relative gain phase through cubic order and canceled the mixed second-order gain–noise response. An eleven-block solution compressed the echo's stretch multiplier from `24` to `11.8761128`: **49.484% of its duration and exposure at the same return index**. The changed detuning coefficients were calculated alongside that saving.

The subsequent detuning work followed two routes. A raised word canceled the complete quadratic static-noise logarithm on the whole Hilbert space. The shorter encoded construction instead closed the full second-order gain/reference-detuning error, including leakage, with separately certified **17-, 19- and 21-block solutions at `n=1,2,3`**. Their stretch multipliers are `23.1836`, `23.3059` and `23.8993`, all below the echo budget of `24`.

Two reusable mathematical results came out of this work:

- **Leading-response dimension:** time-reflection symmetry reduces the first surviving encoded error to six real coordinates in the two-code/two-bright system. The odd-time version, relevant to affine drift, has five.
- **Gain-order doubling:** for an analytic unitary primitive fixing the reference exactly, two signed harmonic conditions raise an encoded error amplitude `O(εᵖ)` to `O(ε²ᵖ)` after common phase removal. The native primitive therefore reaches **eighth-order pure-gain infidelity**.

The complete leading joint gain/reference infidelity is sixth order. The control roots have interval existence certificates; they are more than floating-point optimization candidates. An arbitrary-order finite Fourier recurrence supplies the higher response coefficients.

The resource distinction remains useful: these short encoded words retain gain/reference temporal filter onsets `O(ω⁴)` and `O(ω²)`, while the longer raised construction achieves `O(ω⁶)` and `O(ω⁴)`. Affine reference drift is the next coherent-control target.

Sources: [composite control](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-reflection-loop-composite-control.md), [reference echo](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-reflection-loop-reference-echo.md), [compression](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-reflection-loop-compression.md), [raised response](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-reflection-loop-second-order-detuning.md), [complete encoded correction](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-reflection-loop-complete-code-correction.md).

## 4. Endpoint erasure enabled the short detuning-corrected sequence

SPLIT-1 introduced bright dumps **between exact-return loops**. At a nominal return, the bright sector is empty, so an absorbing measurement there flags error-generated population. Its three-loop trine symmetrizes the leading gain-induced erasure across logical inputs. The subsequent projected-map calculation established its eighth-order conditional gain infidelity analytically.

The five-loop extension adds an orthogonal inverse pair and corrects detuning. Its chronological controls are

\[
\frac{\alpha}{\pi}=\left(\frac16,\frac56,\frac76,\frac34,\frac14\right),
\quad \text{directions}=(+,+,+,-,-),
\]
\[
\text{stretches}=(1,1,1,s,s),\qquad s=\frac{3\sqrt3}{4}.
\]

Each loop includes finite entry/exit ramps and ends with a finite absorbing hold at zero controlled Hamiltonian. Uncontrolled detuning persists during the entire sequence. For dump duration `t_d` and rate `κ`, the no-click dump is

\[
D=e^{-\kappa t_d/2}Q+
e^{-it_d(\delta_dD_d+\delta_rD_r)}P.
\]

The first three loops' traceless detuning response and the stretched pair's response cancel exactly. Equal dump times contribute only common phase at first order. This proves first-order logical correction for the same six static generators at **every integer return index**, including finite dump strength.

Under ideal projection, the orthogonal pair multiplies the trine's surviving map by a scalar, up to a fixed code-basis conjugation. Consequently its conditional gain fidelity is exactly the trine's at every gain value with nonzero success probability. The gain-only mean erasure begins at

\[
\bar p_{\rm erase}=\frac{25}{2}b_n^2\epsilon^2+O(\epsilon^3),
\qquad b_n=\frac\pi2\left(1-\frac1{16n^2}\right).
\]

The implemented finite dump retains and counts surviving bright amplitude. At fixed finite dump strength there is a small quadratic leakage floor; the eighth-order gain regime is therefore not the strict asymptote all the way to zero error.

With `n=1`, ramp slew parameter `ν/a²=10`, `κ/a=10` and `κt_d=40`, the **complete duration is `aT=124.8735`**. At `ε=10⁻³`, `Δ_d/a=10⁻⁴`, `Δ_r=0` and zero background absorption:

- Heralded erasure: **0.00269845%**.
- Conditional average infidelity, including surviving bright leakage: **2.24985 × 10⁻¹³**.
- The matched three-loop trine, with the same ramps and finite dumps, gives **2.40807 × 10⁻⁸**.

That is approximately **107,000-fold lower conditional infidelity**, for a duration increase from `68.2015/a` to `124.8735/a`, at this common operating point. At 1% gain and zero detuning, the five-loop gate gives **0.265153% heralded erasure** and **2.49522 × 10⁻¹⁵ conditional infidelity**.

Sources: [SPLIT-1](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-CLAIRE-split-loop.md), [projected gain identity, SC-7](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-reflection-loop-complete-code-correction.md#sc-7-comparison-with-the-endpoint-erasure-construction), [finite-dump derivation](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-reflection-loop-finite-dumps.md).

## 5. The final result includes the unheralded loss budget

The latest extension separates monitored dump absorption from unmonitored terminal background absorption during transport. Background losses stay in the accepted ensemble. Let their positive effects be `E_d` and `E_b`, and let `K` be the surviving map:

\[
K^\dagger K+E_d+E_b=I.
\]

With monitoring efficiencies `η_d,η_b`, set

\[
E_h=\eta_dE_d+\eta_bE_b,\qquad
E_u=(1-\eta_d)E_d+(1-\eta_b)E_b,
\quad p_{h,u}=\tfrac12\operatorname{tr}_P E_{h,u}.
\]

For `A=GK_PP`, `B=K_QP` and `A₀=A−tr(A)I/2`, the requested operating pair is

\[
\boxed{(p_h,q_u),\qquad
q_u=\frac{p_u+\|B\|_F^2/2+\|A_0\|_F^2/3}{1-p_h}.}
\]

Here `q_u` is the average input-state infidelity weighted by acceptance, **per accepted run**. Total absorption is `p_h+p_u`; residual error per attempted run is `(1−p_h)q_u`. The receipt separately retains the Bell-pair infidelity convention and the full state-dependent code maps and loss effects for a chosen decoder calculation.

All rows below use `n=1`, `aT=124.8735`, `κ/a=10`, `κt_d=40`, `Δ_d/a=10⁻⁴`, `Δ_r=0`, perfect dump monitoring and unmonitored background rate `γ`:

| Gain error ε | Background γ/a | Heralded fraction p_h | Unheralded residual q_u |
|---:|---:|---:|---:|
| 0.1% | 0 | 0.00269845% | 2.24985 × 10⁻¹³ |
| 0.1% | 10⁻⁶ | 0.00269841% | 6.06981 × 10⁻⁶ |
| 1% | 0 | 0.265058% | 3.47408 × 10⁻¹² |
| 1% | 10⁻⁶ | **0.265055%** | **5.93456 × 10⁻⁶** |
| 1% | 10⁻⁴ | 0.264728% | 5.93144 × 10⁻⁴ |

The highlighted pair is the starred point on the figure; its total absorption is **0.265647%**. It includes the background loss needed for an operational comparison. A decoder threshold has not yet been calculated.

![Finite-dump gate: duration, exposure, total absorption, and paired heralded and unheralded operating points](/home/williaml/seated-root/docs/figures/reflection-loop-finite-dump-tradeoffs.png)

The curves contain 168 operating points across `n=1,2,3`, both gain values and background rates from zero through `10⁻³a`. Independent density-matrix evolution with six qubit input states and separate absorbing registers checks the probability partition and accepted-state error.

Source: [FD-7 and full operating-point receipt](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-reflection-loop-finite-dumps.md#fd-7-the-pair-when-background-absorption-is-unheralded).

## 6. What increasing n actually buys

Higher `n` reduces integrated transport exposure. For the five-loop gate, the nominal mean exposures `aX̄` are **6.08633, 3.23334 and 2.17930**, while durations `aT` are **124.8735, 232.0902 and 338.1656**.

Endpoint gain-induced flags are controlled by the return amplitude. At 1% gain and zero detuning their probabilities are **0.265153%, 0.288862% and 0.288808%**. They do not inherit the exposure's `1/n` decrease.

The decrease does appear for background absorption during transport. With 1% gain, zero detuning and `γ/a=10⁻³`, total absorption falls from **0.851286%** to **0.594991%** and **0.490018%**. Higher `n` is therefore useful when that loss dominates. The residual detuning response and added duration must be included when choosing the operating point.

## 7. Strongest next targets

1. **Optimize the complete operating pair.** Use the existing maps to compare `(p_h,q_u,T)` over `n`, dump strength and monitoring efficiency at fixed error conditions. At `γ/a=10⁻⁶`, unobserved loss already dominates the displayed residual. Reducing exposure or observing those losses has the immediate contribution in that regime.
2. **Correct the quadratic full instrument.** Cancel the remaining logical detuning/mixed-gain distortion and state dependence of the flag effect while retaining finite dumps and their duration. The existing arbitrary-order response matrices supply the starting equations. The target includes both the surviving map and the flag effect.
3. **Complete affine reference-drift correction on the coherent branch.** Add its five independent response conditions to the certified fifteen-equation construction and seek a solution within `W≤24`. This advances the alternative that does not use heralding.
4. **Make a matched quantum-error-correction comparison.** Supply the complete instrument to a specified code/decoder, retaining leakage and state-dependent flags, and compare at equal duration and noise assumptions. This is the next step from the reported physical error pair to a logical error-rate prediction.

The most developed candidates for a research contribution are the full real-fibre orientation construction, the general response/order-doubling results, and the short detuning-corrected finite-dump instrument. Their derivations establish the stated behavior; historical originality and an advantage over competing implementations require the matched comparison.

## 8. Evidence and saved work

The committed receipts distinguish exact identities, interval certificates and numerical evolution. The current records include:

- General branched/resolved response: **163 checks**.
- Exact-return primitive: **34 symbolic and 28 numerical checks**.
- Composite control: **16 symbolic and 36 numerical checks**; reference echo: **11 symbolic and 45 numerical checks**.
- Compression: **103 checks**; second-order detuning: **105 checks**.
- Complete encoded correction: **124 checks and three interval root certificates**.
- Finite-dump instrument and loss readout: **139 numerical checks and six symbolic identities**.

The independent [quantum-information harness](/home/williaml/seated-root/qi_error_harness.py) reproduced the exact-return and initial composite results and checked the compressed word's printed moment equations. Later finite-dump verification includes separate chronological and density-matrix integrations. These are complementary checks with the scopes recorded in their receipts.

The main executable endpoints are [coherent encoded correction](/home/williaml/seated-root/reflection_loop_short_correction.py) and [finite-dump instrument](/home/williaml/seated-root/reflection_loop_finite_dump.py). Their detailed outputs are [the coherent receipt](/home/williaml/seated-root/docs/receipts/reflection-loop-short-checks.json) and [the finite-dump receipt](/home/williaml/seated-root/docs/receipts/reflection-loop-finite-dump-checks.json).

Key saved commits are `be1ae74` (general response theory), `ab25706` (exact returns), `74fb2a3` (reference echo), `eeaee3a` (compression), `597450d` (raised response), `5ed475f` (complete encoded correction), `b72203c` (SPLIT-1), `a4a372f` (five finite-dump loops) and `df70de0` (paired loss budget). The [external-mathematics ledger](/home/williaml/seated-root/docs/EXTERNAL-MATHEMATICS-DEBT.md) retains the general response dependencies under EXT-004/007–009 and the quantum realization and comparison sources under EXT-010–018. This report adds no external mathematical dependency.

# The five-loop gate in distance-dependent network QEC

6 September 2026

**The five-loop gate supports decreasing logical failure rates as the surface-code distance increases in this declared network model.** At the main operating point, logical X-memory failure falls from **0.107% at distance 3** to **0.023% at distance 5**, **0.003% at distance 7** and **0.001% at distance 9**. The last value represents two failures in 200,000 trials and has correspondingly wide uncertainty. The qLDPC code also completes syndrome extraction and recovery on all twelve logical observables.

The experiment inserts the finite-dump gate on every data qubit before each noisy syndrome-extraction round, measures the actual syndrome network, decodes the complete record, and checks the final logical observables. Every trial is decoded, including trials with heralded erasures.

The two working code families are rotated surface codes at distances **3, 5, 7 and 9**, and a **`[[72,12,6]]` bivariate bicycle qLDPC code**. Separate X- and Z-memory experiments test the two error sectors. Each block contains `d` noisy rounds. The qLDPC failure event is an error in **any of its twelve measured logical observables**; the surface-code event concerns its single logical observable.

## Main operating-point results

The following X-memory results use `γ/a=10⁻⁶`, the 1% native gain error and the `p_stack=10⁻³` noisy syndrome stack specified below. Every row contains **200,000 trials**. The two decoder columns use the same sampled trials.

| Surface-code distance | Flags used: logical failure | Flags hidden: logical failure |
|---:|---:|---:|
| 3 | 0.107% — 214 failures | 0.1625% — 325 failures |
| 5 | 0.023% — 46 failures | 0.0465% — 93 failures |
| 7 | 0.003% — 6 failures | 0.008% — 16 failures |
| 9 | 0.001% — 2 failures | 0.001% — 2 failures |

Z-memory gives 229, 26, 7 and 0 failures with flags used, versus 363, 72, 15 and 0 with flags hidden. Its distance-9 result is a **95% upper bound of `1.498 × 10⁻⁵` per block**, not a measured zero rate. The benefit of retaining flags is resolved at the smaller distances; these distance-9 samples do not resolve a decoder difference.

At distance 9, **170,988 of the 200,000 X-memory trials contained flags**: 85.5% of the blocks. The decoder retained all of them. The two observed logical failures are therefore a result of erasure correction across the whole experiment, with no rejection of flagged blocks.

The ideal-native-gate baseline, retaining the same noisy extraction stack, has 170, 25, 3 and 1 X-memory failures. This measures the contribution already present in the CNOT/preparation/measurement stack. The finite-dump gate adds errors, and the decoder recovers part of that added cost using its flags.

For the qLDPC code, **10,000 trials per basis**, each containing six noisy rounds, give:

| qLDPC memory | Flags used: any-logical failure | Flags hidden: any-logical failure |
|---|---:|---:|
| X, twelve logical observables | 0.10% — 10 failures | 0.12% — 12 failures |
| Z, twelve logical observables | 0.08% — 8 failures | 0.11% — 11 failures |

The corresponding ideal-native baselines have 7 and 4 failures. These qLDPC counts establish a working network/decoder example; they do not yet resolve the small difference between the two decoder configurations. The receipt retains per-logical counts and confidence intervals. Comparing its any-of-twelve failure probability directly with a single surface-code logical would use different success criteria.

### Increased background loss

With `γ/a=10⁻³`, all other native and stack parameters fixed, the surface-code X-memory failures are **757, 229, 66 and 13 out of 200,000**, respectively. Hidden-flag decoding gives **920, 320, 114 and 37**. The Z-memory counts with flags used are **709, 233, 60 and 14**. Distance suppression remains visible in this more lossy case.

The qLDPC counts become **32/10,000 for X memory and 20/10,000 for Z memory**, versus 43 and 31 when the flags are hidden. These observations belong to this finite code and decoder configuration; the experiment does not locate an asymptotic threshold.

![Rotated surface-code logical failure rates, with binomial intervals and erasure flags used or hidden](/home/williaml/seated-root/docs/figures/qec-distance-stack.png)

In the figure, `five` is the main `γ/a=10⁻⁶` point, `background` is `γ/a=10⁻³`, and `baseline` is the ideal native gate in the noisy stack. Solid and dashed curves compare flags used and hidden. The final [receipt](/home/williaml/seated-root/docs/receipts/qec-distance-stack-checks.json) contains **30 cases and 4,860,000 Monte Carlo trials**. The [30 sampled network definitions](/home/williaml/seated-root/docs/receipts/qec-distance-circuits) are saved separately by scenario and memory basis.

## The physical interface being tested

The supplied native point is `n=1`, gain error `ε=0.01`, `Δ_d/a=10⁻⁴`, `Δ_r=0`, `κ/a=10`, dump action `κt_d=40` and gate duration `aT=124.8735`. The main background rate is `γ/a=10⁻⁶`; a second run uses `10⁻³`. The complete matrices come from the [finite-dump receipt](/home/williaml/seated-root/docs/receipts/reflection-loop-finite-dump-checks.json).

Two additional operations define the QEC interface:

1. An ideal boundary repair preserves states in the code subspace and replaces every non-code outcome by the maximally mixed qubit. Unobserved loss is repaired without acquiring a classical flag; recorded dump loss keeps its flag.
2. Independent ideal Pauli wrappers randomize the instrument at each data-gate location. The experiment averages over this specified protocol. It does not replace an unrandomized coherent channel by its average infidelity.

These operations are explicit realization assumptions. Their finite implementation, additional noise and cost remain to be supplied. The present test measures the resulting gate interface inside real syndrome circuits; it does not model a lost carrier persisting through later entangling gates.

### Exact instrument-to-decoder calculation

Let `K` be the four-dimensional surviving map from the two-dimensional input code, `E_h` the monitored loss effect and `E_u` the unmonitored loss effect. Remove the target reflection from its code block:

\[
A=GK_{PP},\qquad B=K_{QP},\qquad R=E_u+B^\dagger B.
\]

Boundary repair gives the two completely positive maps

\[
\Phi_h(\rho)=\operatorname{tr}(E_h\rho)\frac I2,\qquad
\Phi_0(\rho)=A\rho A^\dagger+\operatorname{tr}(R\rho)\frac I2.
\]

Their sum is trace preserving because `A†A+R+E_h=I`. Expanding `A=Σ_P a_P P`, with `a_P=tr(PA)/2`, and averaging its conjugations by the four Paulis gives

\[
h=\frac12\operatorname{tr}E_h,\quad r=\frac12\operatorname{tr}R,
\qquad
\overline\Phi_h(\rho)=h\operatorname{tr}(\rho)\frac I2,
\]
\[
\boxed{\overline\Phi_0(\rho)=\sum_{P\in\{I,X,Y,Z\}}w_P P\rho P,
\qquad w_P=|a_P|^2+\frac r4,
\qquad \sum_Pw_P=1-h.}
\]

Pauli orthogonality removes the cross terms; averaging `R` makes it `rI`. The conditional Pauli probabilities are `w_P/(1−h)`. This derives the full decoder input, including its bias and loss contribution.

The same identity extends to an input code of dimension `q` with an orthogonal Weyl operator basis: `a_μ=tr(W_μ†A)/q`, `r=tr(R)/q` and `w_μ=|a_μ|²+r/q²`. It uses the declared maximally mixed replacement and Weyl randomization. The runner implements and checks the qubit case against direct four-Pauli averaging on a complete operator basis.

At the main operating point, the resulting probabilities are

\[
h=0.00265054791,\qquad
(p_X,p_Y,p_Z)\approx
(1.48363896,1.48364415,1.48363897)\times10^{-6}.
\]

The total conditional nonidentity Pauli probability is `4.45092208 × 10⁻⁶`. The raw instrument's `5.93455930 × 10⁻⁶` residual counted a missing carrier as orthogonal to the target. Replacing it by `I/2` gives one-half overlap with a pure target; the repaired conditional average infidelity is therefore `2.96728138 × 10⁻⁶` here. This change belongs to the declared repair channel.

## Syndrome circuits and decoders

The non-native network noise is fixed at **`p_stack=10⁻³`**. Preparation and measurement suffer the specified flip probability; Clifford operations and idle locations use their listed depolarizing channels. The native gate has its own matrix-derived channel above. The baseline uses an ideal native-gate location with the same noisy extraction stack, so it measures the stack's existing contribution.

**Surface codes:** Stim's rotated-memory circuits contain ancilla preparation, Hadamards, four layers of syndrome CNOTs, measurement and reset. Native gate layers are inserted at round starts with all measurement-record references preserved. Data noise before each round is also included. The logical X and Z sectors are decoded separately using PyMatching's minimum-weight matching. [Maintainer documentation](https://github.com/oscarhiggott/PyMatching).

**qLDPC:** the check matrices are constructed from commuting cyclic shifts,

\[
A=x^3+y+y^2,\quad B=y^3+x+x^2,\quad x^6=y^6=I,
\qquad H_X=[A\ B],\quad H_Z=[B^T\ A^T].
\]

The runner independently computes ranks and logical quotient bases. An exhaustive pair/triple meet-in-the-middle calculation excludes every nontrivial logical of weight at most five and exhibits weight-six logicals in both sectors. The network uses the source's seven CNOT-layer ordering. Both ancilla families are reset at cycle start and read at its end; added waits receive idle noise. The initial syndrome cycle and final data readout are ideal boundary conditions. [Code construction and source schedule](https://arxiv.org/pdf/2308.07915), [authors' parameters](https://github.com/sbravyi/BivariateBicycleCodes/blob/main/decoder_setup.py).

The qLDPC decoder is BP+OSD: minimum-sum belief propagation, 30 iterations, scaling factor `0.625`, and combination-sweep order 2. It decodes the network's fault-signature matrix, including propagated CNOT errors. Its returned correction is checked against the syndrome equation. [Decoder documentation](https://software.roffe.eu/ldpc/quantum_decoder.html).

For both families, heralds set the corresponding error-column probabilities to `1/2` for that shot. Hidden-flag decoding uses the unconditional channel priors on **the same sampled trials**. Other priors use the stated independent-mechanism and separate-sector decoder approximations; the sampled network retains its Pauli correlations. The `HERALDED_ERASE` network channel itself is sampled directly, preserving its flag/error correlation. Its approximate detector-model expansion is not used to generate shots. [Stim's channel specification](https://github.com/quantumlib/Stim/blob/main/doc/gates.md#the-heralded_erase-instruction).

## Verification and replay

The [runner](/home/williaml/seated-root/qec_distance_stack.py) checks the native channel conversion at `n=1,2,3`, noiseless encoded evolution and complete-erasure randomization. It corrects every individual network-fault signature in both test families: **55 per surface-code sector at distance 3 and 2,232 per qLDPC sector**. Erasure-aware recovery is checked separately: **1,392 patterns per surface sector** cover all pairs of native erasure locations, and **4,096 patterns per qLDPC sector** cover all error assignments on 128 sampled sets of five erased data qubits. The qLDPC distance proof checks all **59,640 triples per sector**, not a random subset.

An [independent forced-event test](/home/williaml/seated-root/test_qec_distance_stack.py) then activates one actual heralded erasure at a time. All **918 locations and 14,688 network shots** pass: the recorded flag identifies the intended location, and its observed syndrome/logical effect belongs to that location's decoder columns. This checks the flag-to-decoder connection independently of decoding the model's own columns.

The same test file verifies exact Clifford-tableau equality between the qLDPC seven-layer extraction unitary and sequential measurement of its stated check operators. This establishes the extraction action on all input states, beyond the two prepared memory states.

The runtime dependencies are pinned in [qec_stack_requirements.txt](/home/williaml/seated-root/qec_stack_requirements.txt). The receipt records package versions, random seeds, input and runner hashes, every failure count, logical-bit counts, paired-decoder disagreements and 95% binomial intervals. Zero-failure rows give a one-sided 95% upper bound. They are not zero-error claims.

```bash
uv venv --python 3.13 /home/williaml/.cache/seated-root-qec-venv
uv pip install --python /home/williaml/.cache/seated-root-qec-venv/bin/python -r qec_stack_requirements.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /home/williaml/.cache/seated-root-qec-venv/bin/python rlq/qec_distance_stack.py \
  --check --shots 200000 --bb-shots 10000 --distances 3 5 7 9 \
  --scenarios baseline five background \
  --json docs/receipts/qec-distance-stack-checks.json \
  --plot docs/figures/qec-distance-stack.png --circuits docs/receipts/qec-distance-circuits
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /home/williaml/.cache/seated-root-qec-venv/bin/python -m unittest -v test_qec_distance_stack
```

External code constructions and decoder dependencies are retained in [EXT-019/020 and search receipt SR-007](/home/williaml/seated-root/docs/EXTERNAL-MATHEMATICS-DEBT.md). The simulations condition on fixed gain/detuning parameters, independent reservoirs and independent Pauli wrappers. Correlated drift, imperfect repair and a native entangling gate are concrete continuation targets.

## What this advances

The finite-dump instrument now has a working path through encoding, repeated noisy extraction, erasure-aware decoding and logical readout. Its flags have measurable value in the surface-code stack, while its remaining error contribution sits above the separately measured syndrome-stack baseline.

The immediate physical extension is a finite, noisy realization of the repair/randomization interface. The next native control target is an entangling-gate instrument for the syndrome CNOTs, whose present errors are supplied independently. Those extensions would replace the remaining ideal interface operations and external gate-noise assumptions with calculated native operations. A larger qLDPC distance series and more rare-event samples would then test the resulting stack more sharply.

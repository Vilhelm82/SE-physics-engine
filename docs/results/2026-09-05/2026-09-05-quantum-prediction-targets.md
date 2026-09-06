# Quantum many-body and quantum-information prediction targets

5 September 2026 — source comparison, derived feasibility results and next theorem targets

**Recommendation:** begin with Cella's fixed-spectrum coupling fibres and their dark-space holonomy; pursue critical quantum memory as the main many-body programme; develop operator-response lifting as the bridge to scrambling. The most immediate contribution is a concrete orientation-reversing dark-space loop. The largest mathematical opportunity is a quantum response theory uniform through changes of the effective state space.

These are predictions within a declared quantum realization of the mathematics. Their mathematical validity, originality in the literature, experimental implementation and robustness are separate questions. The elementary witnesses below establish feasibility; the cited literature sets the baseline the extensions must improve.

## 1. The three proposed correspondences, made operational

**Criticality — singular elimination, BR-7.** Split a finite quantum generator into retained and discarded operator coordinates:

\[
\dot x=Ax+By,\qquad \dot y=Cx+Dy.
\]

Exact elimination gives

\[
\dot x(t)=Ax(t)+Be^{Dt}y(0)+\int_0^t Be^{D(t-s)}Cx(s)\,ds,
\qquad \Sigma(z)=B(zI-D)^{-1}C.
\tag{1}
\]

The quasistatic approximation is `A-BD^{-1}C`. When a discarded relaxation rate approaches zero, its inverse can become singular and its memory can become long. Its effect on a selected observable depends on the residue through `B,C`, preparation and readout. A vanishing eigenvalue with zero observable residue is not a divergent response in that channel. The upgrade target is the complete retained quantum dynamics at this singularity, including cancellations and the order of the time, size and coupling limits.

Equation (1) is ordinary linear elimination. BR-7's classical passivity conclusion does not by itself make a quantum reduction completely positive. Use the quantum compatibility condition below. A finite-size closing relaxation rate and a thermodynamic phase transition also require distinct limit statements.

**Topological protection — stable divisor constraints and transport, BR-4/8.** The relevant questions are which responses vanish, which perturbations preserve that vanishing, and which discrete transport invariant cannot change under the allowed deformations. A divisor supplies a zero and its order. A topological protection theorem additionally identifies an invariant and the class of gap-preserving perturbations. Section 2 constructs precisely such a discrete orientation invariant on a native Cella family.

**Scrambling — lifting the actual operator response, BR-3.** Full finite-dimensional unitary evolution is a regular invertible linear map on density operators. Information becomes difficult to recover through a restricted set of local probes. A branched or rank-deficient map must therefore be identified at that observation/response level. Section 4 supplies a positive matrix of quantum commutator responses on which congruence and divisor calculus apply; quantum realizability is supplied by its Hamiltonian or generator.

The common extension is **prediction on a retained operator state, with specified interventions and readouts**. It allows these mechanisms to meet without presuming that every instance of one phenomenon has the same singularity.

### Quantum realization to retain as standard input

Use a finite Hilbert space, positive unit-trace states and completely positive trace-preserving (CPTP) evolution. A Markovian example is

\[
\mathcal L\rho=-i[H,\rho]+\sum_\alpha
\left(L_\alpha\rho L_\alpha^\dagger-\tfrac12\{L_\alpha^\dagger L_\alpha,\rho\}\right),
\tag{2}
\]

in units with `hbar=1`. This is a supplied quantum realization, not a derivation of quantum postulates from Seated Root.

Exact observable-dependent quantum reduction already exists, including controlled generators and reduced Lindblad form. The 2025 controlled extension explicitly closes under all admissible control sequences. That result is a strong reusable starting point. Our extension target is parameter-dependent rank changes, uniform finite-time accuracy and compatible transport across the resulting strata. [Exact reduction](https://arxiv.org/html/2412.05102v3), [controlled reduction](https://arxiv.org/html/2510.25546v1).

## 2. First target: fixed-gap dark-space gates from REALFIBER

**Native source inspected:** [REALFIBER theorem](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/geometric_fault_localization_and_decomposition/REALFIBER_THEOREM.md>), obligations O1–O7. It supplies the whole real symmetric zero-diagonal family

\[
\mathcal F_m(a)=\{H:\operatorname{spec}H=(a,0^{m-2},-a),\ \operatorname{diag}H=0\},\quad a>0,
\tag{3}
\]

its connectedness for `m>=3`, and the requirement that active coupling partitions change through zero-amplitude faces. Interpret `H` as an `m`-level Hamiltonian, or as the one-excitation sector of an exchange-coupled spin network. Its dark projector is

\[
P_0=I-H^2/a^2.
\tag{4}
\]

It is smooth even where the active coupling graph changes. The dark-to-bright gap remains exactly `a`. Thus a coupling-coordinate singularity need not close the spectral protection gap.

### Derived extension: the full real dark-space holonomy is O(m−2)

For the connection `P_0 d` on the real dark bundle over (3), allow piecewise smooth loops through the zero-amplitude strata. At the single-edge base point `H_*=a(e_1e_2^T+e_2e_1^T)`,

\[
\boxed{\operatorname{Hol}_{H_*}(\ker H)=O(m-2),\qquad m\ge3.}
\tag{5}
\]

**Proof.** Parallel transport preserves the real inner product, so the group is contained in `O(m-2)`. The fixed-star subfamily `H=a(e_1q^T+qe_1^T)`, with `q` a unit vector perpendicular to `e_1`, has dark bundle `T S^{m-2}`. On a two-sphere using any two dark coordinate directions, an orthonormal latitude frame has connection

\[
\begin{pmatrix}0&-\cos\theta\\\cos\theta&0\end{pmatrix}d\phi.
\]

A meridian–latitude–return loop gives any rotation angle modulo `2pi`, while all other dark directions stay fixed. Coordinate-plane rotations generate `SO(m-2)` for `m>=4`; for `m=3` this subgroup is trivial.

To obtain the other component, set `H=a(pq^T+qp^T)` and use three successive arcs, each parameterized by `s in [0,pi/2]`:

\[
\begin{aligned}
1:&\quad p=\cos s\,e_1+\sin s\,e_3,&&q=e_2,\quad
D=-\sin s\,e_1+\cos s\,e_3;\\
2:&\quad p=e_3,&&q=\cos s\,e_2+\sin s\,e_1,\quad
D=-\cos s\,e_1+\sin s\,e_2;\\
3:&\quad p=\cos s\,e_3+\sin s\,e_2,&&q=e_1,\quad
D=\cos s\,e_2-\sin s\,e_3.
\end{aligned}
\tag{6}
\]

Each arc has orthonormal, disjoint-support `p,q`, zero diagonal and spectrum (3). The normalized dark vector satisfies `D^T dD=0`; the unused basis vectors are constant dark spectators. The endpoints match and `H` returns to `H_*`, but `D:e_3 -> -e_3`. Thus the dark holonomy is `diag(-1,1,...,1)`. Together with the rotations it generates all of `O(m-2)`. Each corner can be smoothed in time by stopping the parameter derivatives at the arc endpoints; the path and gap are unchanged. ∎

This is a full mathematical extension of the inspected Cella family. Its originality relative to the complete literature remains to be established. Holonomic gates and real `SO(k)` dark-space operations are established: a 2025 integrated-photonics experiment implements them up to six dark modes. The specific comparison target is the orientation-reversing component accessed by changing the active partition in the full fixed-spectrum fibre. [Chen et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12006547/).

**Immediate prediction.** Add the spectator `e_4`. In the adiabatic limit the logical basis `(e_3,e_4)` undergoes `diag(-1,1)`. The superposition `(e_3+e_4)/sqrt(2)` becomes `(-e_3+e_4)/sqrt(2)`: populations agree, while an interference measurement distinguishes them perfectly. Finite-time evolution needs its adiabatic error evaluated against the actual pulse duration and gap. No exact finite-time gate or general immunity to noise is asserted.

The determinant of real dark-space holonomy is a discrete invariant: it cannot change continuously under deformations preserving a real isolated band and the closed-loop identification. This protects the **orientation parity**. It does not protect every matrix entry of a higher-dimensional gate or cancel arbitrary dynamical phases. Complex couplings, band mixing and uncontrolled splitting of the dark band require their own error analysis.

**Finite-duration continuation:** [RF-1–7](../../../../results/2026-09-05/2026-09-05-reflection-loop-finite-duration.md) gives exact returns, nonadiabatic holonomy and perturbation response. [CP-1–4](../../../../results/2026-09-05/2026-09-05-reflection-loop-composite-control.md) cancels common first-order gain error. [Reference-access control, RE-1–6](../../../../results/2026-09-05/2026-09-05-reflection-loop-reference-echo.md) removes relative gain phase through cubic order and cancels static detuning and mixed gain–detuning response, with a derived erasure law and explicit duration/loss comparisons. [Compression and return-index tradeoffs, CM-1–5](../../../../results/2026-09-05/2026-09-05-reflection-loop-compression.md) gives a certified eleven-block member at 49.484% of the echo's duration/exposure, its full changed error response, and linear gain-drift cancellation. [Second-order detuning, SD-1–6](../../../../results/2026-09-05/2026-09-05-reflection-loop-second-order-detuning.md) adds certified logical cancellation at \(n=1,2,3\), full quadratic noise-log cancellation with exact gain preservation, and higher temporal rejection orders. [Complete encoded correction, SC-1–6](../../../../results/2026-09-05/2026-09-05-reflection-loop-complete-code-correction.md) closes the remaining leakage target within the echo budget at all three returns and raises pure-gain infidelity to eighth order. Next cancel affine reference drift within that budget.

### Immediate many-body continuation: occupation changes the protection

For a four-spin XY star, let the three centre-to-leaf exchange amplitudes be `x,y,z`. Its one-excitation matrix has the fixed spectrum (3) when `x^2+y^2+z^2=a^2`. In the two-excitation basis `(12,13,14;23,24,34)`, however,

\[
H_2=\begin{pmatrix}0&B\\B^T&0\end{pmatrix},\qquad
B=\begin{pmatrix}y&z&0\\x&0&z\\0&x&y\end{pmatrix},
\]
\[
\boxed{\det B=-2xyz,\qquad \det H_2=-4x^2y^2z^2.}
\tag{7}
\]

Consequently the one-excitation dark band persists while the two-excitation zero modes disappear when all three couplings are nonzero. Near `z=0`, with `xy!=0`, the pair departing zero has magnitude

\[
|E_{\rm small}|=\frac{2|xy|}{x^2+y^2}|z|+O(z^2).
\tag{8}
\]

The probe obtains (7) by restricting the full 16-dimensional spin Hamiltonian, independently of the displayed block. This gives a concrete occupation-dependent spectroscopic prediction and identifies exactly where a one-excitation construction must be upgraded before claiming many-body protection. The result is elementary; its value here is a decisive boundary and a calculable continuation, not an asserted first discovery of the XY spectrum.

## 3. Main many-body target: critical memory and changing reduced quantum state

The native inputs are [BR-3/4/7](../../../../results/2026-09-05/2026-09-05-response-on-branched-resolved-state-spaces.md), [predictive-state closure](../../../../results/2026-09-05/2026-09-05-cella-predictive-state.md), and Cella's normalized-incidence and special-fibre constructions. Retain the slow states before taking the limit that removes their relaxation rates.

The baseline is substantial. Quantum metastability already uses separated Liouvillian modes. A March 2026 construction preserves CPTP dynamics and supplies gap-dependent error bounds and perturbative reduction on a fixed reference centre manifold. The stronger target is a family with competing closing gaps and a changing reduced algebra, including controlled passage between its strata. [Metastability](https://arxiv.org/abs/1512.05801), [2026 reduction](https://arxiv.org/html/2603.11982v1).

### Exact calibration model

In a three-level system take

\[
L_{r,\theta}=r|0\rangle\langle B_\theta|,\quad
|B_\theta\rangle=\cos\theta|1\rangle+\sin\theta|2\rangle,\quad H=0.
\]

Prepare `|1>`. The exact excitation survival is

\[
\boxed{p_{\rm exc}(T)=\sin^2\theta+\cos^2\theta\,e^{-r^2T}.}
\tag{9}
\]

At every finite time the zero-coupling limit is one. For every `r>0`, the long-time limit is `sin^2 theta`. Resolving the joint limit requires `tau=r^2 T` and the coupling direction. The discarded bright population is precisely the state that a uniform description must retain on this timescale. The master equation and a Kraus representation are checked exactly in the probe.

The direction is a known control/constitutive datum when the couplings are known. It is not automatically an extra hidden physical state at `r=0`; hidden state must be identified by distinct futures under the **same** subsequent inputs. Equation (9) is an existing-type dark-state crossover used as a calibration, not the proposed novel prediction by itself.

**Strong theorem target:** construct a CPTP reduced model, or a sharp lower bound on the state that must remain, uniformly for `r -> 0`, `T -> infinity` with prescribed scaled times, and then for growing spin-system size. The new outputs sought are direction-dependent memory lifetimes, transitions between coherent and classical retained memory, and preparation-dependent critical relaxation that a single equilibrium order parameter omits.

An intrinsic leakage operator connects this problem to the divisor theorem. For a fixed code projector `P`, `Q=I-P`,

\[
K_P=P\mathcal L^\dagger(Q)P
=\sum_\alpha P L_\alpha^\dagger QL_\alpha P\succeq0,
\qquad \dot p_Q(0)=\operatorname{tr}(K_P\rho),\quad \rho=P\rho P.
\tag{10}
\]

It is independent of the chosen jump representation. If smooth `K_P(f)` vanishes at a two-sided parameter face, positivity gives `K_P=f^2\bar K_P` with positive smooth `bar K_P`; a one-sided rate can have a first-order zero. This imports BR-4 into an actual quantum response. Logical lifetime still requires the full projected generator: coherent leakage and dephasing inside `P` are not measured by this initial escape rate.

**First decisive calculation:** a collective-loss spin cluster with two independent weak perturbations. Compute its slow generator, its code leakage and the observable memory kernel; compare the full dynamics, the 2026 CPTP reduction and a direction-resolved reduction in the joint limit. Advance to a chain only after the finite cluster supplies an actual new coefficient, uniform bound or state lower bound. If the established reduction already provides the requested uniform prediction, retain it and move to the next degeneration.

### What “minimum memory” must mean

Specify the preparations, interventions, observation time `T` and tolerated error `epsilon`. Minimize the memory Hilbert dimension among physical encodings and updates reproducing all those experimental probabilities within that error. Keep this separate from the dimension of a linear observable span. Quantum models can encode distinguishable predictive states nonorthogonally; a scalar state count does not by itself count qubits.

Process tensors already describe multi-time memory, and quantum Markov order depends on the probing instruments. Restricted unitary sequences can witness and bound temporal correlations without determining the entire process. These tools are retained; the opportunity is their uniform and quantitatively minimal realization near the changing state spaces above. [Operational memory](https://arxiv.org/abs/1801.09811), [instrument dependence](https://arxiv.org/abs/1805.11341), [restricted control](https://quantum-journal.org/papers/q-2025-04-08-1695/), [memory bounds](https://arxiv.org/abs/2412.12812).

For a many-body bath, temporal entanglement of its influence functional provides a relevant computational comparison. Low temporal entanglement and its tensor-network use are already established; the new target is a certified, observable-specific memory/error tradeoff across the singular regime. [Sonner–Lerose–Abanin](https://arxiv.org/abs/2103.13741).

## 4. Scrambling target: quantum response zeros, growth and recovery

Use actual evolved operators `W(t)=U(t)^dagger W U(t)` and fixed local probes `V_alpha`. At infinite temperature define the real positive response matrix for Hermitian probes

\[
C_\alpha(t)=[W(t),V_\alpha],\qquad
K_{\alpha\beta}(t)=\frac1{\dim\mathcal H}\operatorname{tr}\big(C_\alpha(t)^\dagger C_\beta(t)\big)\succeq0.
\tag{11}
\]

Real changes of probe basis act by `K -> A K A^T`. Thus BR-3 applies to an **actual specified singular probe map**, and BR-4 to response zeros. General restriction to local probes is often rectangular; the relation-level transformation and predictive closure then replace a unique inverse lift.

The diagonal entries are squared-commutator diagnostics. Their use for operator spreading, weak-link slowdowns and distinct entanglement/transport timescales is established. Our target is the complete matrix and its parameter-dependent zero/jet structure, including cancellations that scalar diagnostics miss. [Operator spreading](https://arxiv.org/abs/1705.08975), [weak links](https://arxiv.org/abs/1705.10364).

For a controlled analytic finite system define the first observable order of a logical sector by

\[
m_* = \min\{|w|:\ P\,\mathcal L_w^\dagger(O)P
\text{ is not a scalar multiple of }P\},
\tag{12}
\]

where `w` ranges over ordered words of the allowed generators and `O` over permitted readouts. If the set is empty, those readouts never distinguish the sector under the declared controls. Otherwise an appropriate pulse word reveals its first nonzero response coefficient. This is the operator version of the existing predictive-state construction, and a natural place to apply Cella's weighted-jet machinery to competing vanishing couplings.

For a distance-`d_code` code, a weight-`w_0` observable and Hamiltonian terms supported on at most `k>=2` sites, an order-`m` nested commutator is a sum of terms individually supported on at most `w_0+m(k-1)` sites. Apply the code's detection property term by term. Consequently

\[
w_0+m(k-1)<d_{\rm code}\quad\Longrightarrow\quad
P\mathcal L_w^\dagger(O)P\text{ is scalar}.
\tag{13}
\]

The proof is locality of a nonzero commutator followed by the code's detection property. Geometry, conservation and cancellations may delay detection further. This is a baseline bound; a new result would give the exact onset orders and coefficients for a nontrivial controlled family, then a finite-time recovery or distinguishability bound. Finite-order vanishing alone does not give a sharp causal time cutoff. Local Hamiltonians have Lieb–Robinson tails. [Locality baseline](https://arxiv.org/abs/quant-ph/0603121).

**First calculation:** a small interacting spin chain with two tunable links, followed by a five-qubit code under an explicitly listed local control set. Derive the complete first nonzero response jets, identify which distinctions require nonlocal probes, and test whether the resolved description changes any predicted onset or recovery error beyond existing operator-growth calculations. A positive Gram matrix alone is not a count of entangled qubits or a reconstruction of the quantum state.

The probe includes `X_1 -> X_1 X_2` under CNOT: a previously commuting local probe becomes sensitive, while the full density-operator Jacobian stays invertible. This fixes where the lift must live—on selected response data—without identifying every quantum evolution with a branched coordinate map.

### Connection to the orientation-reversing gates

A June 2026 preprint already finds different spreading behaviour for certain random-circuit ensembles drawn from the two components of the orthogonal group. That abstract distinction is therefore prior work. The native opportunity is a physical gate family obtained from (5)–(6), with its actual pulse cost, correlations and noise, and a prediction for the resulting structured circuit. Its gate distribution must be checked against the preprint's invariance hypotheses before using its velocity result. [Tan–Brouwer, 2026](https://arxiv.org/abs/2606.03956).

## 5. Priority and success criteria

1. **Dark-space orientation gates:** strongest immediate native construction. Equation (5) is proved for the full declared fibre; (6) gives a measurable ideal gate. The linked control extensions now supply certified complete second-order encoded gain/reference correction within the echo budget, a general gain-order doubling theorem, arbitrary joint response, and explicit temporal/resource costs. The immediate coherent-control target is affine reference-drift cancellation at that duration. [FD-1–6](../../../../results/2026-09-05/2026-09-05-reflection-loop-finite-dumps.md) supplies a separate five-loop heralded gate with first-order logical detuning correction and finite absorbing dumps; its next target is quadratic correction of the full no-click/flag instrument. The occupation-dependent spectrum (7) is the first many-body extension test.
2. **Critical quantum memory:** strongest broad many-body target. Deliver a uniform CPTP response or a state-retention lower bound at competing gap closures, with an observable relaxation/lifetime prediction. Include time, size, error and control access in every minimality statement.
3. **Scrambling and logical observability:** strongest bridge between the predictive-state and weighted-jet work. Deliver exact response-onset orders and coefficients, then recovery bounds for specified controls. Use the full operator response rather than inferring scrambling from a spectral branch alone.

Exceptional-point sensing remains a secondary application: the full matrix exponential is regular at finite time when the generator is analytic, and divergent eigenvector coordinates must cancel in the physical response. Quantum-noise and resource comparisons are already central to that literature. It becomes worth pursuing here when a resolved calculation supplies a new operational precision bound or protocol. [Sensing baseline](https://arxiv.org/abs/1805.11760).

## 6. Provenance and replay

The inspected native inputs include REALFIBER O1–O7; the Kummer criterion and selected quotient work; the constitutive, tensor-valuation and predictive-state extensions; and DIS chapters 5, 10 and 12. The DIS Hamiltonian is a programme input, not an already constructed quantum channel. The bilinear coupling graph theorem was inspected but is not used as a quantum entanglement criterion. DAG keyword misses were not treated as proof of corpus-wide absence.

[Exact probes](../../../quantum_prediction_probes.py) replay the dark loop and sphere connection, the independently assembled spin sectors, the CPTP crossover and the operator-response example. [Receipt](../../../../receipts/quantum-prediction-probes.json). General proof (5) uses its explicit rotation and reflection generators; computational checks are witnesses, not a novelty certificate. External imports, comparison sources and the search queries are recorded in [the dependency ledger](EXTERNAL-MATHEMATICS-DEBT.md).

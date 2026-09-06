# Second-order detuning cancellation for the finite reflection

The target left by [CM-1–5](../../../../results/2026-09-05/2026-09-05-reflection-loop-compression.md) is closed in two forms:

- **Thirteen blocks:** certified controls at return indices $n=1,2,3$ cancel the complete second-order logical reference-detuning rotation. Second-order leakage remains. These gates retain the previous static-noise, mixed-error, gain-phase and erasure cancellations.
- **Full second-order correction:** a forward/inverse/forward construction cancels the quadratic static-noise logarithm on the whole state space, preserves the parent's gain-only gate exactly, and raises its temporal rejection orders. Applied to the eleven-block parent, it uses 33 blocks.

The shorter construction improves the available duration/loss choices; the general construction provides stronger correction with greater cost. Both costs and finite-error comparisons are below.

Inputs: [finite returns](../../../../results/2026-09-05/2026-09-05-reflection-loop-finite-duration.md), [composite control](../../../../results/2026-09-05/2026-09-05-reflection-loop-composite-control.md), [reference echo](../../../../results/2026-09-05/2026-09-05-reflection-loop-reference-echo.md), CM-1–5, and the native [REALFIBER proof](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/geometric_fault_localization_and_decomposition/REALFIBER_THEOREM.md>). The live Cella DAG's GFL:thm:realfiber node and canonical proof O1–O7 were checked. No new external search or import was used; finite-control composition remains retained EXT-018.

## SD-1. Exact response at arbitrary order

Use angular-frequency units, $Q=\operatorname{span}(e_1,e_2)$, $P=\operatorname{span}(d,r)$, and the target $G=\operatorname{diag}(-1,1)|_P$. Write $C$ for the wrapped, zero-ended primitive from RE-1. Its nominal full endpoint is $K=B|_Q\oplus G$. The reference perturbation is $\Delta V_r$, $V_r=|r\rangle\langle r|$.

First set $a=1$, so time is measured in $a^{-1}$ and detuning is $\delta=\Delta/a$. At any integer return index $n\ge1$,

\[
v=(16n^2-1)^{-1/2},\qquad \tau=\frac{\pi}{2v},\qquad \omega=4nv.
\tag{1}
\]

On each native arc the moving-frame generator has eigenvalues $0,\pm\omega$. If $\widehat L=L/\omega$, its projectors are
$\Pi_0=I-\widehat L^2$, $\Pi_\pm=(\widehat L^2\pm\widehat L)/2$. The real frame has exact integer coefficient matrices
$F(t)=F_0+F_c\cos(vt)+F_s\sin(vt)$, supplied explicitly by the runner. Consequently

\[
U_0(t)=\sum_p A_pe^{ipvt},\qquad
p=k-4nm,\quad k,m\in\{-1,0,1\}.
\tag{2}
\]

Thus $U_0^\dagger VU_0=\sum_k M_ke^{ikvt}$ is a finite Fourier sum for every fixed matrix $V$. Put

\[
I_k=\begin{cases}\pi/2&k=0,\\(i^k-1)/(ik)&k\ne0,\end{cases}
\quad
J_{k,l}=\begin{cases}(I_{k+l}-I_k)/(il)&l\ne0,\\
\pi^2/8&k=l=0,\\
\pi i^k/(2ik)+(i^k-1)/k^2&l=0,\ k\ne0.
\end{cases}
\tag{3}
\]

For $U_\delta(\tau)=U_0+\delta U_1+\delta^2U_2+\cdots$, ordered integration gives exactly

\[
U_1=-\frac{i}{v}U_0(\tau)\sum_kM_kI_k,\qquad
U_2=-\frac1{v^2}U_0(\tau)\sum_{k,l}M_kM_lJ_{k,l}.
\tag{4}
\]

**Arbitrary-order closure.** In $u=vt$, represent each interaction-frame coefficient as a finite sum $\sum_{k,d}C_{k,d}u^de^{iku}$. The recurrence is
$E_q(u)=-(i/v)\int_0^u M(t)E_{q-1}(t)\,dt$, $E_0=I$, using

\[
\int_0^u t^de^{ikt}\,dt=
\begin{cases}
u^{d+1}/(d+1),&k=0,\\
\displaystyle e^{iku}\sum_{j=0}^d
\frac{(-1)^jd!u^{d-j}}{(d-j)!(ik)^{j+1}}
-\frac{(-1)^dd!}{(ik)^{d+1}},&k\ne0.
\end{cases}
\tag{5}
\]

This finite algebra includes every repeated-frequency resonance. It computes response coefficients without fitting finite differences or integrating an ODE. Replacing $V$ by a formal linear combination also gives all ordered multichannel coefficients.

For the reference problem, pulling back $V_r$ by a code rotation gives
$\sin^2\alpha D_d+\cos^2\alpha D_r+\sin\alpha\cos\alpha X_{dr}$. The runner computes the three first-order and nine ordered second-order primitive coefficients. Holds and caps commute with these matrices, so their response is an exact exponential over their actual duration. Their physical ramps remain in the waveform.

A negative backward waveform maps a single-parameter coefficient $U_q$ to $(-1)^qU_q^\dagger$; stretching multiplies it by $s^q$. For temporally ordered multichannel coefficients, reversal also reverses the noise labels: $U_{ij}\mapsto U_{ji}^\dagger$ at second order. Polynomial convolution then gives every composite coefficient. The specialized quadratic sums, arbitrary-order recurrence, high-precision evaluation and independent laboratory ODE agree in the receipt, including the individual ordered coefficients before symmetrizing channel labels. Exactness refers to these identities; floating-point evaluations retain their numerical errors.

## SD-2. A certified thirteen-block solution

Take $\sigma_j=(-1)^{j-1}$ and

\[
\begin{aligned}
\beta&=(b_1,\ldots,b_6,0,-b_6,\ldots,-b_1),\\
s&=(w_1,\ldots,w_6,w_0,w_6,\ldots,w_1),\\
\alpha_j&=2\sum_{k<j}\sigma_k\beta_k+\sigma_j\beta_j.
\end{aligned}
\tag{6}
\]

The physical blocks are rotated $C$ or its negative backward waveform, stretched by $s_j\ge1$. The target is automatically $K$, and $\alpha_{14-j}=\alpha_j$. Impose the CM conditions

\[
\sum_js_je^{i(k\beta_j+\ell\alpha_j)}=0,\quad
(k,\ell)\in\{(1,0),(2,0),(1,2),(1,-2),(2,2),(2,-2)\},
\quad
\sum_js_j\sin2\alpha_j=0,\quad
\sum_j\sigma_je^{2i\beta_j}=0.
\tag{7}
\]

Under (6), these are eight real equations. They correct all six Hermitian generators of
$\mathfrak A=\operatorname{End}(Q)\oplus\mathbb CD_d\oplus\mathbb CD_r$, retain first-order gain correction and the mixed gain–noise cancellation, remove relative gain phase through cubic order, and give scalar first-order erasure.

**Why precisely two more logical equations suffice.** Let $J_Q$ exchange $e_1,e_2$ and fix $d,r$. The native wrapped primitive obeys

\[
H_C(T_c-t)=J_QH_C(t)J_Q.
\tag{8}
\]

This follows directly from the three real arc frames; the paired holds, middle inverse and zero-area caps satisfy the same identity. The palindromic physical angles, signs and stretches in (6) preserve (8). Since $V_r$ commutes with $J_Q$, transposing the time-ordered product gives
$U_\delta^T=J_QU_\delta J_Q$. Therefore $U_{\delta,PP}$ is symmetric.

Write $K^\dagger U_\delta=I+\delta A_1+\delta^2A_2+\cdots$. Conditions (7) give $A_1P=cP$, with imaginary $c$. Unitarity implies
$P(A_2+A_2^\dagger)P=-|c|^2P$, so $\operatorname{tf}_P PA_2P$ is anti-Hermitian. Multiplication of a symmetric $U_{\delta,PP}$ by $G$ makes its off-diagonals opposite. Hence the full second-order logical term has exactly the form

\[
B_\delta=\operatorname{tf}_P PA_2P
=\begin{pmatrix}iz&y\\-y&-iz\end{pmatrix},
\qquad y,z\in\mathbb R.
\tag{9}
\]

Solving $y=z=0$ cancels every logical component. The missing Pauli component is excluded by (8), independently of numerical residuals.

Fix $w_1=w_2=1,\ w_6=1.07$. Equations (7), (9) then form ten equations for ten free variables. [The control file](../../../../receipts/reflection-loop-detuning-controls.json) specifies certified centers and fixed coordinates for $n=1,2,3$, at $\nu/a^2=10$. Their total stretch weights are respectively

\[
W_1=16.425794058480463,\quad
W_2=16.42940267056479,\quad
W_3=16.433298758504733.
\tag{10}
\]

**Existence certificate.** Normalize the eight circular equations by 13 and the two response equations by $T_c^2$. Let $f$ denote this system and $x_0$ its stored decimal center. With a fixed decimal matrix $A$ approximating $Df(x_0)^{-1}$, use $T(x)=x-Af(x)$ on the radius-$10^{-35}$ cube in the ten free coordinates.

All controls in this cube satisfy $1\le s_j<3$; free angular derivatives of $\alpha_j$ have magnitude at most 4. Each quadratic Dyson product has at most $2N+4=30$ rotation factors and two stretch factors. The ordered integration region, divided by $T_c^2$, has measure at most $(3N)^2/2$. Product differentiation therefore bounds each response second partial by

\[
(4\cdot30+2)^2(3\cdot13)^2/2<1.2\cdot10^7.
\tag{11}
\]

The circular bounds are smaller. The certificate uses the deliberately loose common bound $10^{18}$. Interval evaluation of the finite sums supplies the center residual. Interval centered differences with step $h=10^{-35}$, enlarged by $10^{18}h$, enclose the Jacobian; the Hessian bound also encloses its variation across the cube.

The resulting contraction bounds at $n=1,2,3$ are below $2.17,2.62,3.11$ times $10^{-13}$, respectively; center displacements are below $4.9\times10^{-65}$. Thus $T$ maps the cube into itself and contracts. Its iterates converge to a unique fixed point there. Since $ADf$ is nonsingular, $A$ is nonsingular and the fixed point satisfies $f=0$. The certificate establishes exact nearby controls at these three indices; continuation to other indices or slew ratios requires further work.

The 13 blocks contain **221 physical stages**. Orthogonal rotation, negative reversal and stretches $s\ge1$ preserve real zero-diagonal rank-two control, amplitude/slew bounds and nominal parallel transport. At zero-Hamiltonian joins, the attached encoded plane is retained.

## SD-3. Full quadratic noise-log cancellation

This construction applies to any finite analytic unitary waveform, without the real-frame or reflection assumptions. Let its endpoint be $U_{\epsilon,\delta}$, where $\epsilon$ is common multiplicative gain error and $\delta V$ any fixed additive static perturbation. Define

\[
K_\epsilon=U_{\epsilon,0},\quad E_\epsilon(\delta)=K_\epsilon^\dagger U_{\epsilon,\delta},
\quad \log E_\epsilon(\delta)=\delta A_\epsilon+\delta^2B_\epsilon+\delta^3C_\epsilon+O(\delta^4).
\tag{12}
\]

Run the waveform forward, its negative backward waveform stretched by $s$, then forward again. The actual additive perturbation retains its laboratory sign throughout. The endpoint is

\[
\mathcal R_s(U)_{\epsilon,\delta}
=U_{\epsilon,\delta}U_{\epsilon,-s\delta}^\dagger U_{\epsilon,\delta},
\quad
K_\epsilon^\dagger\mathcal R_s(U)
=E_\epsilon(\delta)E_\epsilon(-s\delta)^\dagger E_\epsilon(\delta).
\tag{13}
\]

**Order theorem.** Expanding the three noncommuting exponentials gives

\[
\boxed{\log(K_\epsilon^\dagger\mathcal R_s(U))
=(2+s)\delta A_\epsilon+(2-s^2)\delta^2B_\epsilon
+(2+s^3)\delta^3C_\epsilon+O(\delta^4).}
\tag{14}
\]

The possible $[A_\epsilon,B_\epsilon]$ term cancels between the two ends. A noncommutative symbolic expansion is included in the runner. At $s=\sqrt2$, the quadratic logarithm vanishes **on the entire Hilbert space, for every common gain error** in the analytic neighborhood. At zero detuning,
$\mathcal R_s(U)_{\epsilon,0}=K_\epsilon K_\epsilon^\dagger K_\epsilon=K_\epsilon$ exactly.

The argument also permits $\delta V=\delta\sum_av_aV_a$ for arbitrary real $v_a$. Polarization of (14) cancels all homogeneous quadratic terms among these static noise channels. First-order correction still requires the parent's first-order condition; for the CM parent it holds for $V\in\mathfrak A$.

If $A_0P=cP$, removal of the common phase leaves the entire encoded-input error amplitude $O(\delta^3)$, including leakage. Pure-detuning average infidelity therefore starts at $O(\delta^6)$. The parent's mixed cancellation through total degree two is retained. At nonzero gain, (14) does not assert that the first-order noise response is still scalar on the code.

There is a definite cubic price: $C_\epsilon$ is multiplied by $2+2\sqrt2$. Within the symmetric three-block construction with forward stretch $p\ge1$ and inverse stretch $q\ge1$, quadratic cancellation requires $q=\sqrt2p$. Its least duration multiplier is $2+\sqrt2$; the cubic multiplier $2p^3+q^3$ is positive. This is a minimum for this construction, not a global control bound.

Applied to the eleven-block parent, (13) uses **33 blocks, 561 physical stages**. Its duration and first-order bright exposure both multiply by $2+\sqrt2$; its gain-only gate is unchanged at every $\epsilon$. Every physical ramp remains included.

## SD-4. Temporal orders improve too

Use the error vectors from CM-5:
$k_V(t)=(\operatorname{vec}(QU_0^\dagger VU_0P)/\sqrt2,\operatorname{vec}(\operatorname{tf}_P PU_0^\dagger VU_0P)/\sqrt3)$.
For gain, replace $V$ by $H_0(t)$. Write $M_j=\int(t-T/2)^jk(t)\,dt$.

After transport into the common initial frame, a stretched inverse contributes
$(-1)^js^{j+1}M_j$ for additive noise and $(-1)^{j+1}s^jM_j$ for gain, provided lower bad moments vanish. The powers follow by changing variables in the reversed integral; time offsets multiply only the vanishing lower moments.

The CM parent has $M_{\Delta,0}=0$ and $M_{g,0}=M_{g,1}=0$. In (13), the two forward contributions plus the inverse give $2-s^2$ for the first remaining moment in each case. Therefore

\[
\boxed{M_{\Delta,0}=M_{\Delta,1}=0,\qquad
M_{g,0}=M_{g,1}=M_{g,2}=0,}
\quad
\boxed{F_\Delta(\omega)=O(\omega^4),\qquad F_g(\omega)=O(\omega^6).}
\tag{15}
\]

Thus the raised word rejects affine additive drift and quadratic gain drift to first order in their amplitudes. Equation (15) concerns temporal linear response; (14) separately controls nonlinear static response. The receipt supplies both frequency coefficients, the reference spectra, and moment checks at all three evaluated indices.

## SD-5. Complete leading errors and resource selection

For the short word, after removing common phase, write the second-order leakage as $\epsilon^2L_g+\delta^2L_\delta$. The gain result inherited from CM is

\[
L_g=\ell_n v_\beta,\quad
\ell_n=\frac{\pi^2r_n}{12}\binom{-2\sqrt3+3i}{\sqrt3-6i},\quad
v_\beta=\left(\sum_j\sigma_j\cos\beta_j,\sum_j\sigma_j\sin\beta_j\right),
\quad r_n=1-\frac1{16n^2}.
\tag{16}
\]

The logical second-order detuning and mixed terms vanish, so the **complete leading joint infidelity** is

\[
1-F_{\rm av}=\tfrac12\|\epsilon^2L_g+\delta^2L_\delta\|_F^2
+O((|\epsilon|+|\delta|)^5).
\tag{17}
\]

| Return $n$ | Coefficient of $\epsilon^4$ | Coefficient of $\epsilon^2\delta^2$ | Coefficient of $\delta^4$ |
|---:|---:|---:|---:|
| 1 | 364.885273 | 1227.532186 | 4390.941169 |
| 2 | 401.635113 | 1181.109230 | 9848.012758 |
| 3 | 408.371283 | 1153.202519 | 16069.746193 |

The eleven-block $n=1$ gain coefficient was 8.45994; its detuning coefficient was $1.43578\times10^6$. Logical detuning cancellation trades against gain sensitivity and leaves measurable leakage. For the raised eleven-block word, the pure $\delta^6$ coefficients at $n=1,2,3$ are $3.72872\times10^{10}$, $2.54260\times10^{12}$, $2.89246\times10^{13}$. These are pure-detuning coefficients; they are not a joint sixth-order expansion.

Restore the amplitude cap $a$ and slew cap $\nu$. With $\theta=\pi/3$,

\[
\tau_n=\frac{\pi}{2a}\sqrt{16n^2-1},\quad
T_c=9\tau_n+\frac{2\theta}{a}+\frac{4a}{\nu}+\frac{2(1+\sqrt2)a}{\nu},
\quad
X_c=9\tau_n\left(\frac2{16n^2}-\frac3{2(16n^2)^2}\right).
\tag{18}
\]

The thirteen-block resources are $T=WT_c,\ \bar X=WX_c/2$. The raised resources are the CM eleven-block values multiplied by $2+\sqrt2$. Under the RE absorbing-sink model, $p_{\rm loss}=\kappa\bar X+O(\kappa^2)$ is state-independent on encoded inputs to first order.

At $\nu/a^2=10$, the following finite coherent errors use laboratory ODEs and the full waveform. Both columns have $\epsilon=10^{-3}$.

| Construction | $n$ | $aT$ | $a\bar X$ | $1-F_{\rm av}$, $\delta=10^{-4}$ | $1-F_{\rm av}$, $\delta=10^{-3}$ |
|---|---:|---:|---:|---:|---:|
| Thirteen blocks | 1 | 948.265 | 53.575 | $3.73864\times10^{-10}$ | $6.70302\times10^{-9}$ |
| Thirteen blocks | 2 | 1892.462 | 28.468 | $4.10035\times10^{-10}$ | $5.73805\times10^{-8}$ |
| Thirteen blocks | 3 | 2827.072 | 19.192 | $4.17527\times10^{-10}$ | $5.01380\times10^{-7}$ |
| Full second order | 1 | 2340.822 | 132.252 | $8.40675\times10^{-12}$ | $3.71836\times10^{-8}$ |
| Full second order | 2 | 4670.576 | 70.258 | $1.16709\times10^{-11}$ | $2.52534\times10^{-6}$ |
| Full second order | 3 | 6975.530 | 47.355 | $3.80531\times10^{-11}$ | $2.84661\times10^{-5}$ |

For a unitary endpoint, fidelity is evaluated stably as
$\|QU_\delta P\|_F^2/2+\|\operatorname{tf}(GU_{\delta,PP})\|_F^2/3$.
This avoids subtracting almost-equal numbers when errors are tiny. Step-doubling and independent full-waveform integration check the numerical amplitudes. Doubling from 120 to 240 integration steps per time unit changes all six $\delta=10^{-4}$ infidelities by less than $1.6\times10^{-7}$ relatively.

**Joint choice of $n$ and construction.** Restrict the candidate set to these gates, the prior echo and eleven-block gates, each at $n=1,2,3$. Require $aT\le3000$, set $\epsilon=10^{-3},\delta=10^{-4}$, and compare $I_{\rm coherent}+(\kappa/a)(a\bar X)$. This is a weak-loss estimate using the nominal loss coefficient; mixed loss/error corrections are beyond this comparison.

| $\kappa/a$ | Selected gate | Estimated infidelity | Previous best | Improvement factor |
|---:|---|---:|---|---:|
| $10^{-13}$ | Full second order, $n=1$ | $2.16320\times10^{-11}$ | Echo, $n=1$ | 2.55 |
| $10^{-10}$ | Thirteen blocks, $n=3$ | $2.33675\times10^{-9}$ | Eleven blocks, $n=1$ | 1.72 |
| $10^{-9}$ | Thirteen blocks, $n=3$ | $1.96097\times10^{-8}$ | Eleven blocks, $n=2$ | 1.17 |

![Duration, exposure, pure detuning order, finite-error loss choices and reference drift response](../../../../figures/reflection-loop-detuning-tradeoffs.png)

## SD-6. Replay and continuation

[Response runner](../../../reflection_loop_detuning_order.py) · [interval certificate](../../../reflection_loop_detuning_certificate.py) · [controls](../../../../receipts/reflection-loop-detuning-controls.json) · [105-check receipt and full matrices](../../../../receipts/reflection-loop-detuning-checks.json) · [dependency ledger](EXTERNAL-MATHEMATICS-DEBT.md).

The checks cover all three root certificates, every static-noise generator, mixed response, the complete logical target, independent response evaluations, cubic logical scaling, sextic detuning infidelity, temporal moment orders, finite ramps, loss slope and replay of the prior finite-error values. The existing compression runner's 103-check gate also passes after its coefficient interface was generalized to accept the thirteen-block controls.

**Completed continuation:** [SC-1–6](../../../../results/2026-09-05/2026-09-05-reflection-loop-complete-code-correction.md) cancels the remaining $QP$ reference-detuning response with certified 17-, 19- and 21-block controls at $n=1,2,3$, all within the prior echo's $W=24$ budget. Its additional signed moment raises pure-gain infidelity to eighth order, and its bivariate recurrence supplies the complete leading joint error. The next coherent-control target is affine reference-drift cancellation within the same budget. No global duration optimum is claimed.

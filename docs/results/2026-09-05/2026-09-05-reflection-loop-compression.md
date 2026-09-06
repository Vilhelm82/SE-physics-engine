# Reflection compression and the return-index tradeoff

An eleven-block construction implements the encoded reflection with **49.484% of the nested echo's duration and bright exposure at the same return index**. It preserves first-order correction of all six static-noise generators, the mixed gain–noise cancellation, relative gain-phase cancellation through cubic order, and state-independent first-order erasure. Its gain coefficient improves; its residual reference-detuning coefficient increases. Both coefficients enter the comparison below.

Increasing the return index $n$ independently lowers bright exposure. The reference echo works at every $n$, but its residual detuning and temporal-noise response must be included when selecting $n$.

Inputs: [CP-1–4](../../../../results/2026-09-05/2026-09-05-reflection-loop-composite-control.md), [RE-1–6](../../../../results/2026-09-05/2026-09-05-reflection-loop-reference-echo.md), and the native [REALFIBER proof](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/geometric_fault_localization_and_decomposition/REALFIBER_THEOREM.md>). Controls and perturbations have angular-frequency units. $P=\operatorname{span}(d,r)$, $Q=\operatorname{span}(e_1,e_2)$, and $G=\operatorname{diag}(-1,1)|_P$.

## CM-1. A response theorem for alternating reflections

Let $C$ be any finite zero-ended waveform with propagator in

\[
U_C(t)\in U(Q\oplus\mathbb C d)\oplus I_r,
\qquad C_0=B|_Q\oplus G.
\tag{1}
\]

The native wrapped composite satisfies (1). Its gain-perturbed endpoint obeys $C_0^\dagger C_\epsilon P=P+O(\epsilon^2)$. Its reference state is fixed at every gain and every intermediate time.

Choose an odd number $N$ of blocks, signs $\sigma_j=(-1)^{j-1}$, stretches $s_j\ge1$, and angles $\beta_j$. Define

\[
\phi_0=0,\qquad
\alpha_j=\phi_{j-1}+\sigma_j\beta_j,\qquad
\phi_j=\phi_{j-1}+2\sigma_j\beta_j.
\tag{2}
\]

The physical block is $C$, rotated by $O_{\alpha_j}$ and stretched by $s_j$, when $\sigma_j=+1$. When $\sigma_j=-1$, use its **negative backward waveform**, with the same rotation and stretch. Rotations specify coupling amplitudes; no instantaneous gate is applied. If $\sum_j\sigma_j\beta_j\in\pi\mathbb Z$, the encoded endpoint is $G$. The bright endpoint is $B$.

Write the primitive response map as

\[
\mathcal M_C(V)=\int U_C(t)^\dagger VU_C(t)\,dt.
\]

At nominal gain the complete response of this alternating word is exactly

\[
\boxed{\mathcal M_S(V)=\sum_j s_j O_{\beta_j}\,
\mathcal M_C(O_{\alpha_j}^{T}VO_{\alpha_j})\,O_{\beta_j}^{T}.}
\tag{3}
\]

**Proof.** Before a forward block the accumulated bright endpoint is $I_Q$; before an inverse it is $B$. The inverse's response is $C_0\mathcal M_C(V)C_0^\dagger$, so its preceding $B$ cancels. On $P$, the corresponding matrices are $O_{\beta_j}$ by (2), including the additional $G$ in an inverse. Stretching multiplies additive-noise moments by $s_j$. This proves (3) on the full state space, including the ramps.

The following finite moment conditions suffice:

\[
\boxed{\sum_j s_j e^{i(k\beta_j+\ell\alpha_j)}=0,
\quad (k,\ell)\in\{(1,0),(2,0),(1,2),(1,-2),(2,2),(2,-2)\},}
\tag{4}
\]
\[
\boxed{\sum_j s_j\sin(2\alpha_j)=0.}
\tag{5}
\]

Then, for every Hermitian

\[
V\in\mathfrak A=\operatorname{End}(Q)\oplus
\mathbb C|d\rangle\langle d|\oplus\mathbb C|r\rangle\langle r|,
\qquad \boxed{\mathcal M_S(V)P=c_S(V)P.}
\tag{6}
\]

**Proof of sufficiency.** For a bright-only perturbation, the primitive's $QP$ response has only a $d$ column and its $PP$ response only a $dd$ entry. The $k=1,2,\ell=0$ conditions respectively cancel the rotated column and make the projected response scalar.

For $V_r=|r\rangle\langle r|$, pulling back by $O_\alpha$ gives

\[
\sin^2\alpha\,|d\rangle\langle d|+
\cos^2\alpha\,|r\rangle\langle r|+
\sin\alpha\cos\alpha\,(|d\rangle\langle r|+|r\rangle\langle d|).
\]

After propagation by (1), its two $QP$ columns are arbitrary fixed vectors multiplied by $\sin^2\alpha$ and $\sin\alpha\cos\alpha$. Conditions (4) with $k=1$ cancel both. Its projected diagonal coefficients have the same angular weights; its real off-diagonal coefficient is multiplied by $\sin\alpha\cos\alpha$. Conditions (4) with $k=2$ cancel their traceless parts. The imaginary off-diagonal part is invariant under a real plane rotation; (5) cancels it. Interchanging sine and cosine proves the $V_d$ case. Linearity proves (6). The argument permits arbitrary primitive overlap integrals and arbitrary bright dynamics.

This is a family of sufficient constructions, not a claimed classification of every correcting control.

## CM-2. A certified eleven-block member

Use the antisymmetric angle list and symmetric stretch list

\[
\begin{aligned}
\beta&=(b_1,b_2,b_3,b_4,b_5,0,-b_5,-b_4,-b_3,-b_2,-b_1),\\
s&=(1,w_2,w_3,w_4,1,1,1,w_4,w_3,w_2,1).
\end{aligned}
\tag{7}
\]

The target-angle condition is automatic, $\alpha_{12-j}=\alpha_j$, and (4)–(5) reduce to seven real equations. Add the gain-phase condition

\[
\boxed{\sum_j\sigma_j e^{2i\beta_j}=0,}
\tag{8}
\]

which is real under (7). Thus eight real equations determine the five $b$'s and three $w$'s. One solution has

\[
\begin{aligned}
(b_1,b_2,b_3,b_4,b_5)&\simeq
(-0.3815850024,-1.3390135686,-2.5280157828,2.6041829583,1.4660588421),\\
(w_2,w_3,w_4)&\simeq(1.2644794881,1.0639750241,1.1096018909),\\
W=\sum_j s_j&=11.876112806251932\ldots.
\end{aligned}
\tag{9}
\]

Angles are radians. The rounded values in (9) are for reading; the runner supplies 62-digit centers and computes the physical angles from (2).

**Existence certificate.** Let $f$ be the eight reduced equations and $x_0$ the stored decimal center. With a fixed decimal approximation $A$ to $Df(x_0)^{-1}$, interval arithmetic bounds the map $x\mapsto x-Af(x)$ on the radius-$10^{-40}$ cube around $x_0$. Its contraction factor is below $3\times10^{-37}$, and its center displacement below $5\times10^{-62}$. It therefore maps the cube into itself and its iterates converge to a unique fixed point there. Since $ADf$ is nonsingular, so is $A$, and the fixed point satisfies $f=0$. All stretches remain at least one. This establishes an exact solution near the numerical center; an ODE fit is not the existence argument.

The physical waveform has 187 stages, including every entry, exit and sign ramp. Orthogonal rotation preserves spectral norm and rank; the original $PP$ Hamiltonian block is zero, so mixing $d,r$ also preserves zero diagonal. A stretch $s\ge1$ reduces amplitude and slew by $s^{-1}$ and $s^{-2}$. The complete waveform retains $\|H\|\le a$, $\|\dot H\|\le\nu$, real zero-diagonal rank two away from $H=0$, and nominal parallel transport. At $H=0$, the attached encoded plane is retained as in RE-1.

## CM-3. Gain, mixed errors and loss

For the native composite, set $r_n=1-1/(16n^2)$, $b_n=\pi r_n/2$, and

\[
\ell_n=\frac{\pi^2r_n}{12}
\begin{pmatrix}-2\sqrt3+3i\\\sqrt3-6i\end{pmatrix},\qquad
v=\sum_j\sigma_j(\cos\beta_j,\sin\beta_j)
=(0.68870509679255\ldots,0).
\]

The complete second-order code column is

\[
\boxed{S_0^\dagger S_\epsilon P=P+\epsilon^2
\begin{bmatrix}\ell_n v\\ i\sqrt3 b_n^2I_P\end{bmatrix}+O(\epsilon^3).}
\tag{10}
\]

To obtain (10), insert CP-2's primitive column into (3)'s alternating frame, using a minus sign for inverse gain errors and no stretch factor. Equation (8) makes the logical phase scalar. Each primitive's logical phase through cubic order is supported on its active line, while a return through its $O(\epsilon^2)$ leakage first affects the logical block at fourth order. Consequently the relative gain phase is $O(\epsilon^4)$. The leading gain infidelity is

\[
1-F_{\rm av}=\frac{15b_n^2(\pi/3)^2}{2}|v|^2\epsilon^4+O(\epsilon^5),
\tag{11}
\]

with coefficient $8.45994228$ at $n=1$, versus $48.15757075$ for the nested echo.

**Mixed cancellation.** Through first order in gain, primitive endpoints have unchanged code columns; their bright block $B_\epsilon$ still cancels between a forward block and the next inverse. Equation (3) therefore remains valid on code columns through first order in $\epsilon$, with $\mathcal M_C$ replaced by $\mathcal M_{C,\epsilon}$. The arbitrary-coefficient proof of (6) applies to this derivative too. Thus the mixed $\epsilon V$ response has neither leakage nor relative logical action for every $V\in\mathfrak A$.

For constant reference detuning $\delta=\Delta/a$, after removing common phase the full leading joint error is

\[
\boxed{1-F_{\rm av}=\tfrac12\|\epsilon^2L_g+\delta^2L_\delta\|_F^2
+\tfrac13\|\delta^2B_\delta\|_F^2+O((|\epsilon|+|\delta|)^5).}
\tag{12}
\]

The receipt stores all matrices, including the raw converging mixed residuals. At $n=1$, (12) is

\[
8.45994228\epsilon^4+94.4246044\epsilon^2\delta^2
+1.43577729\times10^6\delta^4.
\tag{13}
\]

Under the same absorbing bright-loss model as RE-5, let $X_c=3X_n$ be the active-input exposure of one composite. The primitive loss moment has only one code column. Conditions (4) give

\[
\mathcal L_SP=\frac{WX_c}{2}P,\qquad
K_SP=(1-\kappa WX_c/4)S_0P+O(\kappa^2),\qquad
p_{\rm loss}=\kappa WX_c/2+O(\kappa^2).
\tag{14}
\]

With one flag for all sinks this is state-independent erasure to first order. The resulting mean exposure is $\bar X_S=WX_c/2$, compared with $12X_c$ for the echo. Both duration and exposure ratios are $W/24=0.4948380336\ldots$.

## CM-4. What increasing $n$ buys

For every integer $n\ge1$, put

\[
\tau_n=\frac{\pi}{2a}\sqrt{16n^2-1},\qquad
X_n=3\tau_n\left(\frac{2}{16n^2}-\frac{3}{2(16n^2)^2}\right).
\]

Exposure decreases strictly with $n$. Indeed, writing $v=1/(4n)$, the derivative of $v\sqrt{1-v^2}(2-3v^2/2)$ is positive for $0<v\le1/4$: its numerator is $2-17v^2/2+6v^4>0$. Also $X_n\sim3\pi/(4an)$. Duration grows as $n$, while the gain coefficient in (11) remains bounded. The same all-$n$ cancellation proofs apply to both words.

The following calculation uses $a=1$, $\nu=10$, and includes finite ramps. $A_\Delta$ multiplies $\delta^4$ in the leading joint infidelity.

| Gate | $n$ | $aT$ | $a\bar X$ | $A_\Delta$ | Finite infidelity at $\epsilon=10^{-3},\delta=10^{-4}$ |
|---|---:|---:|---:|---:|---:|
| Echo | 1 | 1385.526 | 78.280 | $7.42334\times10^3$ | $4.73488\times10^{-11}$ |
| Echo | 2 | 2764.501 | 41.586 | $1.16678\times10^5$ | $6.17170\times10^{-11}$ |
| Echo | 3 | 4128.796 | 28.029 | $6.39959\times10^5$ | $1.13494\times10^{-10}$ |
| Eleven blocks | 1 | 685.611 | 38.736 | $1.43578\times10^6$ | $1.53140\times10^{-10}$ |
| Eleven blocks | 2 | 1367.980 | 20.578 | $2.43112\times10^7$ | $2.44377\times10^{-9}$ |
| Eleven blocks | 3 | 2043.085 | 13.870 | $1.23409\times10^8$ | $1.23599\times10^{-8}$ |

Thus the echo's $n=2,3$ exposure is respectively **53.12% and 35.81%** of its $n=1$ exposure. The quartic detuning coefficient simultaneously increases by factors $15.72$ and $86.21$; that tradeoff was absent from a gain-only comparison.

For the finite witness above, allowing all three echo durations and using $I\simeq I_{\rm coherent}+\kappa\bar X$, echo $n=2$ overtakes echo $n=1$ at $\kappa/a\simeq3.92\times10^{-13}$; $n=3$ overtakes $n=2$ at $3.82\times10^{-12}$. These are pairwise weak-loss comparisons under constant errors.

At the original echo's $n=1$ duration budget, the compressed $n=2$ gate fits and has **3.804 times less exposure**. Its finite coherent error is higher, and the pairwise loss crossover is $\kappa/a\simeq4.15\times10^{-11}$. Selection must use the allowed duration and the noise model together.

![Duration, exposure, complete joint error, and finite-frequency reference response](../../../../figures/reflection-loop-compression-tradeoffs.png)

## CM-5. Correlation constraints from the full temporal response

For any error waveform operator $V_a(t)$, define $M_a(t)=U_0(t)^\dagger V_a(t)U_0(t)$ and the error vector

\[
k_a(t)=\left(\frac{\operatorname{vec}(QM_a(t)P)}{\sqrt2},
\frac{\operatorname{vec}(\operatorname{tf}_P PM_a(t)P)}{\sqrt3}\right),\qquad
\widehat k_a(\omega)=\int e^{i\omega(t-T/2)}k_a(t)\,dt.
\tag{15}
\]

For weak, zero-mean classical errors $\eta_a(t)$, direct expansion of the first-order error amplitude gives

\[
\mathbb E[1-F_{\rm av}]_{(2)}=
\sum_{a,b}\int\!\int C_{ab}(t,t')\,
k_a(t)^\dagger k_b(t')\,dt\,dt',\quad
C_{ab}(t,t')=\mathbb E[\eta_a(t)\eta_b(t')].
\tag{16}
\]

Equation (16) retains correlations between gain and reference noise. The spectral Gram matrix $\widehat k_a(\omega)^\dagger\widehat k_b(\omega)$ is stored in the receipt. For a single random-phase sinusoidal error of RMS amplitude $\sigma$, the leading infidelity is $\sigma^2\big(\|\widehat k(\omega)\|^2+\|\widehat k(-\omega)\|^2\big)/2$. This comparison requires no assumed correlation time or bath spectrum.

**Linear gain-drift theorem.** Both constructions also cancel $\eta_g(t)=g_0+g_1t$ to first order, up to common code phase. For a primitive gain kernel, its zeroth code moment vanishes. Time reversal with a minus sign makes the first time moment of an inverse block equal to $+s_j$ times the forward moment on code columns; the time-offset terms multiply the vanished zeroth moment. Equations (4) with $\ell=0$ then cancel leakage and scalarize the projected first moment. For the nested echo the same argument uses its four equal sign weights and the primitive's vanishing projected gain kernel. Hence

\[
\widehat k_g(0)=\widehat k'_g(0)=0,\qquad
F_g(\omega)=O(\omega^4),\qquad
\widehat k_\Delta(0)=0,\quad F_\Delta(\omega)=O(\omega^2).
\tag{17}
\]

At $n=1,a=1,\nu=10$, the leading gain coefficients $F_g/\omega^4$ are approximately $3.6073\times10^{10}$ for the echo and $4.6454\times10^9$ for the compressed word. The reference coefficients $F_\Delta/\omega^2$ are $4.40494\times10^6$ and $1.15243\times10^8$, respectively. Compression improves the low-frequency gain response here and worsens the reference response. Increasing $n$ is therefore conditioned on temporal noise as well as residual static detuning.

## Replay and next target

[Runner](../../../reflection_loop_compression.py) · [full matrices, certificate and checks](../../../../receipts/reflection-loop-compression-checks.json) · [dependency ledger](EXTERNAL-MATHEMATICS-DEBT.md).

The proof of (3)–(6) is independent of the primitive's detailed trajectory. Its exact eleven-block solution is certified by interval contraction. The 103 checks cover arbitrary primitive matrices, all six physical noise generators, the complete mixed jets, finite-waveform integration and step doubling, every ramp, loss slopes, frequency-response convergence, and a finite time-varying detuning witness. The existing 52 composite checks and 56 echo checks also pass. The echo runner now exposes $n$ in its finite-gate and jet interfaces.

Finite-pulse response composition remains the retained EXT-018 baseline. No new external mathematics was adopted. The circular conditions, certified compressed member, all-$n$ resource comparison and gain-drift theorem are native extensions; this document makes no global time-minimum or historical-originality claim.

**Target completed in [SD-1–6](../../../../results/2026-09-05/2026-09-05-reflection-loop-second-order-detuning.md):** certified thirteen-block controls cancel the complete second-order logical detuning rotation at \(n=1,2,3\). A general forward/inverse/forward construction cancels the entire quadratic static-noise logarithm and improves the temporal rejection orders. Duration, gain sensitivity, remaining leakage and exposure enter the joint comparison. The next target is full second-order leakage cancellation within the prior echo's duration budget.

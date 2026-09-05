# Reference-access reflection: phase cancellation and finite echo

**The target is closed for the specified encoded qubit and quasi-static errors.** A three-composite construction removes the relative gain phase through cubic order. A finite five-block echo also cancels the first-order reference-detuning response, including leakage, and the mixed gain–detuning term. Both retain first-order gain cancellation. The remaining leading error is an explicit quartic response.

The echo corrects a larger class than reference detuning: all static site detunings and arbitrary Hermitian perturbations confined to the two bright modes. Under the declared absorbing-loss model, its first-order loss becomes a state-independent erasure when sink outcomes share a loss flag.

The duration and loss cost are substantial. The shorter construction and the previous composite remain useful alternatives; section RE-6 gives the comparison.

Inputs: [CP-1–4](2026-09-05-reflection-loop-composite-control.md), the native reference rotation in [RF-6C](2026-09-05-reflection-loop-finite-duration.md), and the REALFIBER real zero-diagonal control family. Hamiltonians have angular-frequency units. Errors are constant during a complete gate unless stated otherwise. Logical gates are compared up to a common phase on \(P=\operatorname{span}(d,r)\), with \(d=e_3,r=e_4\), and \(Q=\operatorname{span}(e_1,e_2)\).

## RE-1. Finite ramps that preserve the entire endpoint operator

The previous composite \(C_\epsilon\) starts and ends at \(aJ\). Wrap it with a zero-area entry from \(0\) to \(aJ\) and the reversed exit back to zero. If the common fractional gain is \(\epsilon\), both ramp propagators equal the identity exactly because their Hamiltonians commute and their signed areas vanish. Thus

\[
\bar C_{\epsilon,0}=C_\epsilon
\tag{1}
\]

on the full four-mode state, at every gain.

With \(|q|\le a\), \(|\dot q|\le\nu\), the shortest such entry is

\[
0\ \xrightarrow[\ a/(\sqrt2\nu)\ ]{\dot q=-\nu}\ -a/\sqrt2
\ \xrightarrow[\ (a+a/\sqrt2)/\nu\ ]{\dot q=+\nu}\ a.
\tag{2}
\]

Its signed area is \(-a^2/(4\nu)+a^2/(4\nu)=0\), and its duration is
\(\eta=(1+\sqrt2)a/\nu\).

**Minimum-time proof.** For an entry of duration \(h\), any admissible profile obeys \(q(t)\ge\max(-\nu t,a-\nu(h-t))\). The two lines meet at depth \(b=(\nu h-a)/2\); integrating gives the area lower bound \((a^2-2b^2)/(2\nu)\). Zero area therefore requires \(h\ge(1+\sqrt2)a/\nu\). Equation (2) attains it within the amplitude cap.

The wrapped duration is

\[
\bar T_c=T_c+2\eta.
\tag{3}
\]

The encoded plane is stationary on these ramps. At zero Hamiltonian it is retained as attached state; it is not reselected from the enlarged instantaneous kernel. Every subsequent block joins at \(H=0\), so coupling signs and the active reference direction can change continuously without extra unmodelled pulses.

## RE-2. Three rotated composites remove the gain phase

Define

\[
D_d=\operatorname{diag}(1,1,-1,1),\quad
D_r=\operatorname{diag}(1,1,1,-1),\quad D_b=D_dD_r,
\]

and let \(O_\alpha\) rotate the \((d,r)\) plane while fixing \(Q\). Put

\[
A_\epsilon=O_{\pi/3}\bar C_\epsilon O_{\pi/3}^{T},\qquad
\boxed{T_\epsilon=A_\epsilon D_dA_\epsilon D_dA_\epsilon.}
\tag{4}
\]

The physical chronological blocks are \(A\), \(D_dAD_d\), \(A\). A conjugated block is implemented by conjugating its Hamiltonian; the displayed sign matrices are **not instantaneous gates**. Equivalently, their active directions have angles \(\pi/3,2\pi/3,\pi/3\).

Let \(G=\operatorname{diag}(-1,1)\) on \(P\), and \(R=J|_Q\oplus G\). At nominal gain,

\[
T_0=R,\qquad R^2=I.
\tag{5}
\]

Indeed, the bright endpoint of \(C_0\) is \(B=Je^{-2i\pi J/3}\), so \(B^3=J\). The three logical reflections multiply to \(G\).

In the interaction frame of their preceding logical gates, the three active lines have angles \(\pi/3,0,-\pi/3\). Their projectors satisfy

\[
\sum_{j=1}^{3}\widetilde P_j=\frac32 I_P.
\tag{6}
\]

The primitive logical gain error acts only on its active line. Equation (6) makes its leading phase common to both logical states. Endpoint leakage of each \(C_\epsilon\) is \(O(\epsilon^2)\); a contribution returning through that leakage requires two such factors and first enters the logical block at order four. Therefore

\[
QT_0^\dagger T_\epsilon P=O(\epsilon^2),\qquad
\operatorname{tf}_P(PT_0^\dagger T_\epsilon P)=O(\epsilon^4),
\tag{7}
\]

where \(\operatorname{tf}_P X=X-\tfrac12\operatorname{tr}_P(X)I_P\).
This removes the quadratic and cubic relative gain phase.

For every fixed return index \(n\), write \(r_n=1-1/(16n^2)\),
\(b_n=\pi r_n/2\), and \(k_n=\sqrt3\pi^2r_n/4\). The full second-order code column is

\[
R^\dagger T_\epsilon P=P+\epsilon^2
\begin{bmatrix}
L_g\\ i\,3\sqrt3 b_n^2I_P
\end{bmatrix}
+O(\epsilon^3),\qquad
L_g=k_n\begin{pmatrix}0&i\\-1&-2i\end{pmatrix}.
\tag{8}
\]

This follows by inserting CP-2's complete code column into the three factors in (4); their first-order bright blocks cannot contribute to a code column at this order. It gives

\[
1-F_{\rm av}(T_{\epsilon,0})
=3k_n^2\epsilon^4+O(\epsilon^5)
=\frac{9\pi^4r_n^2}{16}\epsilon^4+O(\epsilon^5).
\tag{9}
\]

Removing the relative phase does not by itself minimize leakage. The coefficient in (9) is \(48.15757\) at \(n=1\), versus \(27.24191\) for the previous composite's total gain-only infidelity.

## RE-3. What the trine does to a finite reference perturbation

For any static Hermitian perturbation \(V\), define the complete waveform response

\[
\mathcal M_T(V)=\int_0^{T_T}U_T(t)^\dagger VU_T(t)\,dt,\qquad T_T=3\bar T_c.
\tag{10}
\]

This includes all arcs, holds and entry/exit ramps. For reference detuning \(V_r=|r\rangle\langle r|\), let \(U_{\bar C}(t)\) be the unrotated wrapped composite and set

\[
\mathcal A=\int|\langle d|U_{\bar C}(t)|d\rangle|^2dt,\qquad
\ell=\int\overline{\langle d|U_{\bar C}(t)|d\rangle}\,dt.
\]

Direct rotation and composition give

\[
\boxed{P\mathcal M_T(V_r)P=
\frac{9\mathcal A+3\bar T_c}{8}I_P
-\frac{3\sqrt3}{4}\operatorname{Im}\ell\,Y,}
\quad Y=\begin{pmatrix}0&-i\\i&0\end{pmatrix}.
\tag{11}
\]

**Proof.** A rotated block has projected response
\[
O_\alpha
\begin{pmatrix}
\mathcal A\sin^2\alpha&\ell\sin\alpha\cos\alpha\\
\bar\ell\sin\alpha\cos\alpha&\bar T_c\cos^2\alpha
\end{pmatrix}O_\alpha^T.
\]
Conjugate by the preceding logical reflections and sum at the three angles. Their \(\cos2\alpha=-1/2\) are equal, while their toggled projectors obey (6). The real traceless terms cancel; the imaginary off-diagonal term is (11). The runner checks this identity with arbitrary \(\mathcal A,\bar T_c,\operatorname{Re}\ell,\operatorname{Im}\ell\).

Thus the diagonal phase imbalance is removed, but a finite \(Y\) response and \(Q\mathcal M_TP\) can remain. Both must be included in an echo.

More generally, every Hermitian element of

\[
\mathfrak A=\operatorname{End}(Q)\oplus\mathbb C|d\rangle\langle d|
\oplus\mathbb C|r\rangle\langle r|
\tag{12}
\]

has equal logical diagonal entries after \(T\). For a bright-only \(V\), each rotated block contributes a multiple of its active projector, so (6) makes its logical response scalar. For the two diagonal logical perturbations, use (11) and their sum. This is a six-real-dimensional guaranteed noise class.

## RE-4. A finite echo of the complete error response

Let \(T^{D,[s]}\) mean the trine waveform conjugated by \(D\) and stretched by \(s\):
\[
H^{D,[s]}(t)=D H_T(t/s)D/s.
\]
Let \(T^{-1}\) mean its physical inverse waveform \(-H_T(T_T-t)\). Stretching keeps the intended gate and fractional gain error unchanged; static detuning and loss act for \(s\) times longer.

Execute

\[
\boxed{
T^{I,[1]}\ \longrightarrow\
(T^{-1})^{D_d,[2]}\ \longrightarrow\
T^{D_b,[2]}\ \longrightarrow\
(T^{-1})^{D_r,[2]}\ \longrightarrow\
T^{I,[1]}.}
\tag{13}
\]

Call the endpoint \(E\). All five nominal endpoints are \(R\), so \(E_0=R\). All boundaries have zero Hamiltonian. The full waveform remains real and zero-diagonal, with rank two wherever the coupling is nonzero. No control exceeds the original amplitude or slew caps.

**Response theorem.** For any static Hermitian \(V\),

\[
\boxed{\mathcal M_E(V)=
2\sum_{D\in\{I,D_d,D_b,D_r\}}D\,\mathcal M_T(DVD)\,D.}
\tag{14}
\]

**Proof.** An inverse waveform has response \(R\mathcal M_T(V)R\); its preceding forward endpoint conjugates this back to \(\mathcal M_T(V)\). The cumulative nominal endpoints alternate between \(I\) and \(R\), and each \(D\) commutes with \(R\). Stretching multiplies each response by its duration factor. The repeated identity pattern contributes weight two, matching the other three patterns. This gives (14), including errors during every finite control.

For \(V\in\mathfrak A\), \(DVD=V\). The four sign patterns cancel every \(QP\) entry and every off-diagonal logical entry. RE-3 supplies equal remaining diagonals. Hence, with \(c(V)=\tfrac12\operatorname{tr}_P\mathcal M_T(V)\),

\[
\boxed{\mathcal M_E(V)P=8c(V)P,\qquad V\in\mathfrak A.}
\tag{15}
\]

The complete first-order error on encoded inputs is a common phase. This corrects reference detuning, all four static site detunings, and both Hermitian bright-hopping directions. Equation (14) also gives the response outside this class; such errors are not assumed canceled.

For reference offset \(\Delta\), an exact endpoint formula useful for finite-error evaluation is
\[
\begin{aligned}
E_{\epsilon,\Delta}={}&T_{\epsilon,\Delta}D_r
T_{\epsilon,-2\Delta}^\dagger D_rD_bT_{\epsilon,2\Delta}D_bD_d\\
&\hspace{18mm}\cdot T_{\epsilon,-2\Delta}^\dagger D_dT_{\epsilon,\Delta}.
\end{aligned}
\tag{16}
\]
The negative detuning argument is essential: the inverse control does not reverse the physical detuning. Equation (16) accounts for that fact before taking the adjoint.

## RE-5. Gain preservation, mixed response and the surviving error

Let the trine's second-order gain column be \([L_g;i\chi I_P]\). The alternating forward/inverse signs in (13) give
\[
2L-D_dLD_d|_P+D_bLD_b|_P-D_rLD_r|_P=L.
\tag{17}
\]
Thus the echo retains exactly the trine's full second-order gain column, including its leakage coefficient (9). The relative gain phase remains \(O(\epsilon^4)\).

**The mixed \(\epsilon\Delta\) error also cancels, up to common phase.** Each \(C_\epsilon\) preserves \(P\) to first order in \(\epsilon\). Its endpoint is therefore block diagonal to that order, with a possibly changing bright block and the fixed logical reflection. The trine's projected occupation calculation in RE-3 still applies to that order with gain-dependent \(\mathcal A,\ell\). Its diagonal balance survives differentiation in \(\epsilon\).

Furthermore \(T_\epsilon\) commutes with every node-sign matrix up to \(O(\epsilon^2)\). In the inverse/forward pairs, the actual gain-dependent bright blocks cancel to first order, so the four-pattern averaging proof remains valid through degree one in \(\epsilon\). Therefore
\[
Q\partial_\epsilon\partial_\Delta(R^\dagger E)P=0,\qquad
\operatorname{tf}_P\!\left(P\partial_\epsilon\partial_\Delta(R^\dagger E)P\right)=0
\tag{18}
\]
at zero error.

For numerical coefficients, put \(\delta=\Delta/a\). After removing common phase, the degree-two leakage is \(\epsilon^2L_g+\delta^2L_\delta\), and the degree-two traceless logical error is \(\delta^2B_\delta\). Consequently
\[
\boxed{
1-F_{\rm av}
=\frac12\|\epsilon^2L_g+\delta^2L_\delta\|_F^2
+\frac13\|\delta^2B_\delta\|_F^2
+O((|\epsilon|+|\delta|)^5).
}
\tag{19}
\]

This uses the exact identity \(1-F_{\rm av}=p_{\rm leak,av}+
\|\operatorname{tf}_P(G^\dagger U_{PP})\|_F^2/3\) for unitary four-mode evolution. It counts leakage and removes only common logical phase.

At \(n=1\), \(\nu=10a^2\), the response is
\[
1-F_{\rm av}\simeq
48.15757\,\epsilon^4-5.27269\,\epsilon^2\delta^2+7423.34\,\delta^4.
\tag{20}
\]
The gain coefficient is exact by (9). The detuning and interference coefficients are numerical evaluations of the derived second-order waveform response. Their full matrices and convergence checks are in the receipt. The squared-norm form (19) explains the negative interference term and keeps the complete quartic form nonnegative.

To reproduce all second-order coefficients, the runner integrates one rotated wrapped composite using
\[
\dot U_{pq}=-iH_0U_{pq}-iH_0U_{p-1,q}-iV_rU_{p,q-1},
\quad p+q\le2,
\tag{21}
\]
with absent negative-index terms omitted, \(U_{00}(0)=I\), and all other initial coefficients zero. Coefficients include their Taylor factorials. Matrix convolution then composes (4) and (16), including stretching and the detuning sign. Independent integration of the entire 255-stage waveform checks this construction.

## RE-6. Duration, erasure response and when to use each gate

The full echo takes
\[
T_E=8T_T=24\bar T_c.
\tag{22}
\]
Within this four-pattern averaging construction, this is the shortest duration using blocks of duration at least \(T_T\): a reflection requires an odd number of involution blocks, so four patterns require at least five blocks. One pattern must occur twice. Equal total weight in all four patterns then requires at least \(2T_T\) per pattern. Equation (13) attains \(8T_T\). This is not a minimum over all possible real control paths.

For nominal bright loss, let \(X_c=3X_n\) be CP-3's exposure for the active input. The trine and echo loss moments obey
\[
P\mathcal L_TP=\tfrac32X_cI_P,\qquad
\boxed{\mathcal L_EP=12X_cP},\qquad
\mathcal L=\int U_0^\dagger Q_b(t)U_0\,dt.
\tag{23}
\]
The first identity follows from (6). Inverse traversal and stretching transform the loss moment just as in RE-4; the four-pattern sum also removes its \(QP\) block.

For absorption at rate \(\kappa\) into orthogonal sinks,
\[
K_EP=(1-6\kappa X_c)\,RP+O(\kappa^2),\qquad
p_{\rm loss}=12\kappa X_c+O(\kappa^2).
\tag{24}
\]
To this order, survival preserves the logical state and loss probability is independent of that state. Reporting all sink outcomes as one erasure flag gives a state-independent erasure channel followed by the target gate. Equation (24) is nominal in the other errors; mixed loss corrections are retained by the finite no-jump integration when evaluated.

For \(n=1,\nu=10a^2\), the durations in units \(1/a\) are \(57.2474\) for the previous composite, \(173.1907\) for the trine, and \(1385.5260\) for the echo. Their mean bright exposures are \(3.26165/a\), \(9.78496/a\), and \(78.27970/a\). Peak bright occupation remains bounded by the original \(15/64\) for every encoded input.

**Finite example:** common gain \(\epsilon=0.001\), static reference offset \(\Delta/a=10^{-4}\), and no loss. The calculated average logical infidelities are:

- Previous composite: \(5.47642\times10^{-6}\).
- Trine: \(2.38077\times10^{-7}\).
- Full echo: \(4.73488\times10^{-11}\).

The echo is approximately \(115{,}661\) times better than the previous composite in this declared coherent-error example. It requires about \(24.20\) times its duration and \(24\) times its mean loss exposure.

For any specified coherent-error pair or distribution, choose using
\[
\mathcal I_j\simeq \mathcal I_j^{\rm coherent}
+\kappa\,\bar X_j,\qquad j\in\{C,T,E\},
\tag{25}
\]
with exact coherent propagators where available. For this example, the leading loss budget selects the echo below \(\kappa/a\simeq3.48\times10^{-9}\), the trine between that rate and \(8.03\times10^{-7}\), and the previous composite above it. These boundaries belong to the stated weak-loss model. All three resource lines matter; comparing the echo only against the previous composite misses the useful middle regime.

Quasi-static detuning cancellation does not establish protection against rapidly varying noise; that requires its time-dependent response.

![Gain response, detuning response and loss comparison](/home/williaml/seated-root/docs/reflection-loop-reference-echo-performance.png)

## Evidence, provenance and next target

[Runner](../reflection_loop_reference_echo.py) · [matrix response and checks](reflection-loop-reference-echo-checks.json) · [dependency ledger](EXTERNAL-MATHEMATICS-DEBT.md).

The proofs establish the all-\(n\) gain column, trine occupation identity, four-pattern response theorem, mixed cancellation, and erasure response. Symbolic checks replay their algebra. Numerical checks separately cover finite ramps, the six noise generators, full waveform integration and step doubling, second-order response convergence, finite-error witnesses and absorbing loss.

Finite-strength group averaging and dynamically corrected gates are established control methods. [Khodjasteh–Viola](https://arxiv.org/pdf/0810.0698) supplies the comparison baseline for finite-pulse error composition, group cancellation and stretched controls. Its general idea is retained; the native trine, compiled node-sign construction, complete encoded response and resource comparison above are derived for this family's actual controls.

The next useful target is **compressing this correction under amplitude, slew and noise-correlation constraints**, using (19) and (25) as the objective. The present result supplies an exact feasible construction and a measured cost; a shorter candidate must preserve the complete code response, including mixed errors and the connecting ramps.

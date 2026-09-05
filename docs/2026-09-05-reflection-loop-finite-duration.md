# Reflection loop: finite duration, perturbations and control class

5 September 2026 · continuation of [REALFIBER dark-space holonomy](2026-09-05-quantum-prediction-targets.md)

**Result.** The original real coupling loop implements the exact logical reflection at a discrete, completely classified set of equal arc durations. Its first return is \(aT=3\pi\sqrt{15}/2\). These returns are nonadiabatic holonomies. The calculations below give their transient bright population, complete first-order Hamiltonian error functional, explicit quadratic and quartic error modes, and the distinction between protected orientation and the chosen logical gate.

Use angular-frequency Hamiltonians, \(i\dot\psi=H\psi\). If \(a\) denotes an energy instead, replace \(aT\) by \(aT/\hbar\). All results concern the declared finite-level Hamiltonian; extension to interacting occupation sectors remains a separate calculation.

## RF-1. Exact propagator of the original loop

The active basis is \((e_1,e_2,e_3)\); \(e_4\) is an uncoupled zero-energy reference. The encoded basis at the endpoints is \((e_3,e_4)\). Put \(\vartheta=\pi/2\). The continuous frames \(F_j=(p_j,q_j,D_j)\), for \(0\le s\le\vartheta\), are

\[
F_1=\begin{pmatrix}\cos s&0&-\sin s\\0&1&0\\\sin s&0&\cos s\end{pmatrix},\quad
F_2=\begin{pmatrix}0&\sin s&-\cos s\\0&\cos s&\sin s\\1&0&0\end{pmatrix},\quad
F_3=\begin{pmatrix}0&1&0\\\sin s&0&\cos s\\\cos s&0&-\sin s\end{pmatrix}.
\tag{1}
\]

Their endpoint frames match. \(F_1(0)=I\) and

\[
F_3(\vartheta)=R=
\begin{pmatrix}0&1&0\\1&0&0\\0&0&-1\end{pmatrix}.
\]

On arc \(j\), use \(s=\omega_j t\), \(\omega_j=\vartheta/\tau_j\), and the original Hamiltonian \(H_j=a_j(p_jq_j^T+q_jp_j^T)\). In moving coordinates,

\[
K_j=F_j^TH_jF_j-iF_j^T\dot F_j=a_jJ+\omega_jG_j,
\]
\[
J=\begin{pmatrix}0&1&0\\1&0&0\\0&0&0\end{pmatrix},\quad
G_1=G_3=\begin{pmatrix}0&0&i\\0&0&0\\-i&0&0\end{pmatrix},\quad
G_2=\begin{pmatrix}0&0&0\\0&0&-i\\0&i&0\end{pmatrix}.
\tag{2}
\]

Direct multiplication gives \(K_j^3=\Omega_j^2K_j\), \(\Omega_j^2=a_j^2+\omega_j^2\). Thus, at every intermediate time,

\[
E_j(t)=e^{-itK_j}
=I-i\frac{\sin(\Omega_jt)}{\Omega_j}K_j
+\frac{\cos(\Omega_jt)-1}{\Omega_j^2}K_j^2.
\tag{3}
\]

The full active endpoint propagator and total duration are

\[
U(T)=R\,E_3(\tau_3)E_2(\tau_2)E_1(\tau_1),\qquad T=\sum_j\tau_j.
\tag{4}
\]

The reference amplitude is exactly unchanged. Equations (1)–(4) also cover unequal positive durations and gaps. If the full Hamiltonian must close at its endpoints, require matching endpoint gaps; the encoded projector closes independently.

## RF-2. Complete equal-duration return set

Let \(a_j=a>0\), \(\tau_j=\tau\), and define

\[
\lambda=a\tau,\quad
\phi=\sqrt{\lambda^2+\vartheta^2},\quad
u=\lambda/\phi,\quad v=\vartheta/\phi,\quad
z=e_3^\dagger E_1E_2E_1e_3.
\]

The desired gate is \(\operatorname{diag}(-1,1)\) on the encoded space. It is exact if and only if \(z=1\): unitarity then excludes all endpoint leakage.

**Theorem.** Within this positive-gap, equal-duration, constant-speed family, the complete exact return set is

\[
\boxed{
\phi=2\pi n,\quad
\tau_n=\frac{\pi}{2a}\sqrt{16n^2-1},\quad
T_n=\frac{3\pi}{2a}\sqrt{16n^2-1},
\qquad n=1,2,\ldots .
}
\tag{5}
\]

**Proof.** Write \(c=\cos\phi\), \(b=\sin\phi\). Multiplying (3) gives

\[
\operatorname{Im}z=4bu v^2(c-1)\,[1+v^2(c-1)].
\tag{6}
\]

The last bracket is positive for \(\phi>\vartheta\). Indeed,
\(f(\phi)=\phi^2/\vartheta^2-(1-\cos\phi)\) vanishes at \(\vartheta\), and
\(f'(\phi)=2\phi/\vartheta^2-\sin\phi\ge4/\pi-1>0\).
Since \(u,v>0\), \(z=1\) therefore requires \(\sin\phi=0\).
At an odd multiple of \(\pi\),

\[
z-1=-2v^2(4v^4-8v^2+5)
=-2v^2[4(v^2-1)^2+1]<0.
\]

At every even multiple, all three \(E_j=I\). This proves necessity and sufficiency. ∎

The first return is \(aT_1=18.2510040419\). It is the fastest member of this specified family; no optimum over arbitrary pulse shapes or unequal durations is asserted. More generally, imposing \(\Omega_j\tau_j=2\pi n_j\) separately gives exact returns for unequal arcs.

**The finite-time return is itself holonomic.** On an arc starting in \(D_j(0)\), the moving-frame amplitudes on the two bright-frame vectors have opposite reality types: one is real and the other imaginary. Consequently
\(\chi^\dagger aJ\chi=0\) at every time. At the returns (5), each next arc starts in its dark vector. Together with the isolated reference, the actual evolved code satisfies

\[
P_{\mathrm{evolved}}(t)H(t)P_{\mathrm{evolved}}(t)=0.
\tag{7}
\]

It is cyclic and parallel transported, hence gives a nonadiabatic holonomy. Its evolving code generally differs from the instantaneous kernel of \(H(t)\). The parallel-transport criterion is established background; the exact return classification is the calculation for this native loop. [Pinske–Scheel, section II](https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.4.023086).

## RF-3. Finite-time error and transient bright exposure

For arbitrary durations define \(V=E_3E_2E_1\) and \(z=V_{33}\). Exact operational quantities are

\[
P_{\rm leak}=1-|z|^2,\qquad
F_+=\frac{|1+z|^2}{4},\qquad
F_{\rm av}=\frac{1+\operatorname{Re}z+|z|^2}{3}.
\tag{8}
\]

Leakage is for input \(e_3\), and is the maximum over encoded inputs. \(F_+\) compares the evolved equal superposition with its target. \(F_{\rm av}\) is pure-state fidelity averaged over the encoded qubit, counting leakage as error. To derive the last expression, average
\(|1-p+pz|^2\) over uniform \(p=|\langle e_3|\psi\rangle|^2\in[0,1]\).

The operator norm error on the entire encoded input space is

\[
\delta_0=\|(U-U_{\rm target})P\|
=\sqrt{2(1-\operatorname{Re}z)}
\le\min\left\{2,\sum_{j=1}^3
2\frac{\omega_j}{\Omega_j}
\left|\sin\frac{\Omega_j\tau_j}{2}\right|\right\}.
\tag{9}
\]

This follows by telescoping the product \(V-I\) on \(e_3\), using unitarity and
\(\|(E_j-I)e_3\|^2=2(\omega_j/\Omega_j)^2[1-\cos(\Omega_j\tau_j)]\).
For equal arcs it gives the uniform envelope
\(\delta_0\le\min(2,6v)\), hence \(O((aT)^{-1})\) amplitude error away from or at the exact returns.

Zero final leakage does not mean zero transient excitation. During a resonant arc,

\[
p_{\rm bright}(t)
=1-\left[1-v_n^2(1-\cos\Omega t)\right]^2,\qquad
v_n=\frac1{4n}.
\]

Therefore

\[
\boxed{
p_{\rm bright}^{\max}=4v_n^2(1-v_n^2),\qquad
\int_0^{T_n}p_{\rm bright}(t)\,dt
=T_n\left(2v_n^2-\frac32v_n^4\right).
}
\tag{10}
\]

At the first return these are \(15/64=23.4375\%\) and
\(a\int p_{\rm bright}dt=2.174436028\) (rounded).

For a declared absorbing loss with \(\sum_\alpha L_\alpha^\dagger L_\alpha=\kappa Q_{\rm bright}(t)\), and jumps into sinks with no return, first-order lost probability is \(\kappa\int p_{\rm bright}dt\). Its remainder is \(O((\kappa T)^2)\) for this finite-dimensional bounded generator. This prediction assumes bright-selective loss; arbitrary site loss uses its actual operator in the integral. Longer return members reduce this bright exposure, while accumulated reference detuning below grows with duration.

## RF-4. Complete first-order Hamiltonian error functional

Let \(U_0\) be the exact nominal evolution and \(H_\epsilon=H_0+\epsilon V(t)\), with integrable Hermitian \(V\). Put

\[
M[V]=\int_0^T U_0(t)^\dagger V(t)U_0(t)\,dt,\qquad
\kappa_V=|\epsilon|\int_0^T\|V(t)\|\,dt.
\]

Duhamel's formula and one further iteration give

\[
\|U_\epsilon(T)-U_0(T)\|\le\min(2,\kappa_V),
\]
\[
U_0(T)^\dagger U_\epsilon(T)=I-i\epsilon M+\mathcal R,\qquad
\|\mathcal R\|\le\frac{\kappa_V^2}{2}.
\tag{11}
\]

The remainder bound uses unitarity inside the double integral, so requires no exponential prefactor. Combined with (9), the perturbed encoded error is at most \(\delta_0+\kappa_V\), capped by two.

At an exact return, use \(P=|e_3\rangle\langle e_3|+|e_4\rangle\langle e_4|\), \(Q=I-P\). Then

\[
P_{\rm leak}(\psi)
=\epsilon^2\|QMP\psi\|^2+O(\epsilon^3),\qquad \psi=P\psi.
\tag{12}
\]

First-order leakage vanishes precisely when \(QMP=0\). First-order error of the complete logical gate, allowing an irrelevant common phase, vanishes precisely when

\[
\boxed{QMP=0,\qquad PMP=cP\quad\text{for some real }c.}
\tag{13}
\]

This classifies arbitrary coherent perturbations by their effect on this gate. It includes time-dependent coupling errors, detunings, complex phases and reference coupling; instantaneous gap preservation alone does not replace (13).

For a covariance model \(V=\sum_\alpha \xi_\alpha V_\alpha\), the same linear functional produces the entire leading leakage quadratic form through
\(P M_\alpha Q M_\beta P\) and the covariance of the real coefficients \(\xi_\alpha\). No independent-noise assumption is needed.

**Smooth pulse changes.** The piecewise smooth path (1), viewed as a single path \(H(s)\), is \(a\)-Lipschitz in operator norm. A monotone reparameterization \(\widetilde s(t)\) with the same endpoints therefore obeys

\[
\|\widetilde U(T)-U_0(T)\|
\le a\int_0^T|\widetilde s(t)-s(t)|\,dt.
\tag{14}
\]

This controls smoothing the velocity changes at the joins while preserving the real fixed-gap path. The original resonance condition is exact for its specified schedule; a smoothed schedule has the explicit error bound (14).

## RF-5. Pulse-error jets: a quadratic form and its quartic kernel

Fix the \(n\)-th equal-duration return. Let each gap be multiplied by \(1+\epsilon_j\), keeping its angle schedule and duration fixed. Define

\[
r_n=1-\frac1{16n^2},\qquad b_n=\vartheta r_n.
\]

Differentiating the exact exponential at \(\phi=2\pi n\) gives

\[
(E_3E_2E_1-I)e_3
=b_n\big[(\epsilon_1+\epsilon_3)e_1-\epsilon_2e_2\big]
+O(\|\epsilon\|^2).
\]

Thus

\[
\boxed{
P_{\rm leak}
=b_n^2\big[(\epsilon_1+\epsilon_3)^2+\epsilon_2^2\big]
+O(\|\epsilon\|^3).
}
\tag{15}
\]

There is no first-order logical phase for these gain errors. The quadratic form has exactly the kernel \(\epsilon_2=0,\ \epsilon_3=-\epsilon_1\). For zero-mean stochastic gain errors, its mean is
\(b_n^2[\operatorname{var}\epsilon_1+\operatorname{var}\epsilon_2+
\operatorname{var}\epsilon_3+2\operatorname{cov}(\epsilon_1,\epsilon_3)]\)
to leading order.

Two explicit predictions are

\[
\begin{aligned}
\epsilon_1=\epsilon_2=\epsilon_3=\epsilon:\quad&
P_{\rm leak}=5b_n^2\epsilon^2+O(\epsilon^3),\\
(\epsilon_1,\epsilon_2,\epsilon_3)=(\epsilon,0,-\epsilon):\quad&
P_{\rm leak}=
\big[\vartheta r_n(1-3r_n)\big]^2\epsilon^4+O(\epsilon^5).
\end{aligned}
\tag{16}
\]

For the second line, the middle propagator is \(I\). Expanding
\(E_1(\lambda-\lambda\epsilon)E_1(\lambda+\lambda\epsilon)e_3\)
leaves \(\lambda^2(E_1''-E_1'^2)e_3\,\epsilon^2
=\vartheta r_n(1-3r_n)e_1\epsilon^2\).
This proves the quartic coefficient rather than inferring it from a fitted slope.

At \(n=1\), the coefficients are \(10.84307124\) and \(7.124236652\). The exact leakage at a common +1% gain error is 0.00106145816. The Taylor expansion holds at fixed \(n\); its phase error scales with \(n\epsilon\), so it is not uniform as \(n\to\infty\). Use (3) or (9) in that joint limit. The cancellation requires the stated error correlation; it is not immunity to arbitrary amplitude noise. Duration errors with endpoints reached by rescaling each arc's speed enter through the same product \(a_j\tau_j\), so their first-order jets coincide. An endpoint-angle error is a different perturbation.

## RF-6. What survives perturbation

### A. Real isolated line, fixed reference

Suppose the active real symmetric \(3\times3\) loop is deformed continuously with the middle eigenvalue isolated throughout, and \(e_4\) remains a fixed uncoupled reference. Its real eigenline keeps holonomy sign \(-1\). For perturbations of norm at most \(\eta<a/2\), the straight interpolation from the original active loop keeps its middle-to-bright gap at least \(a-2\eta>0\), by the Hermitian eigenvalue perturbation bound.

In the adiabatic limit, for encoding into that perturbed eigenline, the logical gate is

\[
\operatorname{diag}\left(-e^{-i\int_0^T E_D(t)\,dt},\,1\right).
\tag{17}
\]

The sign is stable; a shifted eigenvalue produces a dynamical phase. Preparation in the original unperturbed code additionally requires its endpoint encoding mismatch to be included.

For real zero-diagonal active matrices with edge weights \(x,y,z\), the exact zero-eigenvalue condition is \(xyz=0\). This is the relevant structural divisor for the triad. Preserving it and the nonzero gap preserves the dark zero energy. An unwanted third edge supplies a calculable phase channel.

### B. Detuning and an unwanted third edge

For static site detunings \(V=\operatorname{diag}(\delta_1,\delta_2,\delta_3,\delta_4)\), the exact nominal resonant orbit has equal integrated occupation \(T_n/3\) on each active site. The three arcs are cyclic permutations of one another and each restarts in its dark vector. Consequently

\[
M_{33}=\frac{T_n}{3}(\delta_1+\delta_2+\delta_3),\qquad
M_{44}=T_n\delta_4.
\tag{18}
\]

Their difference is the first-order relative-phase coefficient. Traceless active detuning cancels this coefficient, but can still produce leakage through \(QMP\).

Reference detuning alone is particularly simple: the gate is exactly
\(\operatorname{diag}(-1,e^{-i\Delta T_n})\), with no leakage and
\(F_+=\cos^2(\Delta T_n/2)\).

For an equal additive offset \(\eta\) on every active edge, define
\(S=\mathbf1\mathbf1^T-I_3\). Then \(V=S\), \(\epsilon=\eta\) in (11).
Set \(k=4n\). Direct integration of the exact resonant orbit yields

\[
\frac{M_{33}[S]}{T_n}
=\frac{I(k)}{\vartheta},\qquad
I(k)=-\frac{2k^4-23k^2+33}{2(k^2-4)(k^2-1)}.
\tag{19}
\]

To obtain this, on arc 1 put \(A=k^{-1}\sin(ks)\),
\(B=1-k^{-2}+k^{-2}\cos(ks)\).
The offset readout is
\(\sin(2s)(A^2-B^2)+2\cos(2s)AB\).
Its integral from zero to \(\vartheta\) is \(I(k)\); the other arcs give the same integral. At the first return, \(M_{33}/T_1=-59/(60\pi)\). As \(n\to\infty\), it approaches \(-2/\pi\), agreeing with the instantaneous dark-energy derivative \(-\sin(2s)\). This is a finite-duration phase prediction, including an error that leaves the active matrix real.

### C. Reference mixing preserves parity but changes the gate

Let \(O_\alpha\) rotate \((e_3,e_4)\) through angle \(\alpha\), fixing \(e_1,e_2\), and replace the whole loop by \(O_\alpha H(t)O_\alpha^T\). This remains real, zero-diagonal, rank two and fixed-gap; its endpoint Hamiltonian is unchanged. Both its adiabatic gate and its exact resonant gate are

\[
G_\alpha=
\begin{pmatrix}
-\cos2\alpha&-\sin2\alpha\\
-\sin2\alpha&\cos2\alpha
\end{pmatrix},
\quad
\det G_\alpha=-1,\quad
\|G_\alpha-G_0\|=2|\sin\alpha|.
\tag{20}
\]

This explicit perturbation preserves orientation parity while rotating the logical reflection axis. Full logical protection therefore needs the reference attachment or an equivalent axis constraint, in addition to the real-band invariant.

## RF-7. Complete exact-tracking control class

The original real Hamiltonian already supplies the finite returns above. If additional complex couplings are allowed, arbitrary smooth timing admits exact tracking of the instantaneous dark projector \(P(t)\). With \(Q=I-P\), every Hermitian correction producing the prescribed parallel transport and zero logical Hamiltonian has the form

\[
\boxed{
H_{\rm corr}=i[\dot P,P]+QZQ,\qquad Z=Z^\dagger.
}
\tag{21}
\]

A specified logical Hamiltonian adds its \(PBP\) block. Proof: differentiating projector transport fixes the off-diagonal block \(QH_{\rm corr}P=iQ\dot PP\); Hermiticity fixes the reverse block; the parallel-transport condition fixes \(PH_{\rm corr}P=0\). The remaining bright block is free.

The canonical choice \(Z=0\) minimizes Hilbert–Schmidt norm and attains the minimum possible operator norm
\(\|H_{\rm corr}\|=\|\dot P\|\), since any admissible correction has that fixed off-diagonal block. For the triad it is

\[
H_{\rm corr}=i(\dot DD^T-D\dot D^T),\qquad
\|H_{\rm corr}\|=|\dot s|.
\tag{22}
\]

It is an imaginary coupling between the two outer modes. More generally, no real symmetric Hamiltonian can exactly transport a nonconstant real projector: \(\dot P=-i[H,P]\) would equate a real matrix with a purely imaginary one, forcing both to vanish. This includes tracking a changing real normalized \(D(t)\), even allowing phase or internal logical evolution. The original real finite-time gate therefore traverses a different, generally complex projector path, with the transient excursions in RF-3.

For constant speed, the corrected Hamiltonian has norm \(\sqrt{a^2+\omega^2}\). For a monotone traversal,
\(\int\|H_{\rm corr}\|dt=3\pi/2\).
If the total Hamiltonian norm is capped at \(A>a\), this particular corrected path requires
\(T\ge3\pi/[2\sqrt{A^2-a^2}]\). This is its resource constraint, not a universal quantum speed limit.

Transitionless driving is retained standard mathematics. Equation (21) completes the projector-level control class and (22) identifies its exact required coupling for this loop. It enlarges the available control set beyond the original real family. [Berry, section 2](https://michaelberryphysics.wordpress.com/wp-content/uploads/2013/07/berry415.pdf).

## Verification and immediate continuation

[Runner](../reflection_loop_dynamics.py) · [receipt](reflection-loop-dynamics-checks.json) · [dependency ledger](EXTERNAL-MATHEMATICS-DEBT.md).

The runner verifies the frame generators, exponential polynomial, return classification identities, zero mean energy, pulse-error derivatives, transient exposure, phase filters and projector correction. Independent laboratory-frame RK4 evolution checks the ideal gate, an unequal off-resonant schedule, the rotated gate and the perturbation remainder; step doubling and unitarity are checked separately. Numerical witnesses do not replace the general proofs.

![Finite-duration and perturbation error curves](/home/williaml/seated-root/docs/reflection-loop-error-curves.png)

The next control-design target is to make the **common** amplitude-error mode satisfy (13), while retaining a specified reference axis and bounded transient loss. Equation (15) supplies the error directions to cancel; equations (10), (18) and (19) supply the competing costs. This is a defined optimization problem over admissible pulses and controls, rather than an unspecified robustness claim.

Native source: REALFIBER O1–O7 and the previous full-family holonomy proof. Standard quantum evolution, nonadiabatic parallel transport, Duhamel/Dyson estimates, Hermitian spectral stability and transitionless driving are retained with provenance. The results here establish this loop's mathematics; literature originality remains open.

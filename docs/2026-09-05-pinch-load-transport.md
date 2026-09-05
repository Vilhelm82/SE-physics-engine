# Pinch selection and area transport from a load

5 September 2026

**The next target now has explicit evolution equations and a complete finite-time example.** A native load's rates select the retained pinch direction; its area response decides whether the same crossing is gravitationally trapped. A smooth reciprocal passive candidate is constructed below that crosses the pinch, changes sheet and forms a future outer trapping horizon. Its metric is regular at the crossing and satisfies null convergence.

The candidate's response and applied effort are declared inputs. The current primitives do not select a unique response: section 3 proves this with equally passive alternatives. The new result is a derivation of the response-to-outcome map, the conditions on a general response, and a solved coupled example.

## 1. Native evolution selects the pinch direction

Use the previous directional resolution, with

\[
x=b-a\gamma,\qquad d=\det[c,h,g],\qquad
\xi=x/d,\qquad F=\frac{1-\xi^2}{1+a^2}.
\]

Near either pinch branch, \((a,x,d)\) are smooth coordinates on the oriented native realization:

\[
\gamma=\sigma\sqrt{1-\frac{d^2-x^2}{1+a^2}},\qquad
b=a\gamma+x,\qquad \sigma=\pm1.
\tag{1}
\]

Let a smooth load-induced state field have rates

\[
\dot a=A(a,x,d),\qquad \dot x=P(a,x,d),\qquad
\dot d=Q(a,x,d).
\]

Dots denote a specified evolution parameter. After soldering, it can be the native clock's duration. For \(d\ne0\), differentiation gives the exact equations

\[
\boxed{\dot\xi=\frac{P-\xi Q}{d},\qquad
\dot F=-\frac{2\xi\dot\xi+2aF\dot a}{1+a^2}.}
\tag{2}
\]

**Crossing theorem.** Suppose a twice differentiable trajectory reaches \(x=d=0\), with \(Q_0=\dot d_0\ne0\). Then its normal has the unique finite directional limit

\[
\boxed{\xi_0=\frac{P_0}{Q_0},\qquad
\dot\xi_0=\frac{\dot P_0Q_0-P_0\dot Q_0}{2Q_0^2},\qquad
F_0=\frac{Q_0^2-P_0^2}{(1+a_0^2)Q_0^2}.}
\tag{3}
\]

**Proof.** Taylor-expand \(x=P_0s+\dot P_0s^2/2+o(s^2)\) and \(d=Q_0s+\dot Q_0s^2/2+o(s^2)\), divide, and substitute into \(F\). The signed volume changes sign at a transverse crossing. ∎

Thus \(|P_0|<|Q_0|\), equality, and \(|P_0|>|Q_0|\) select timelike, null, and spacelike normal limits respectively. At a null crossing, the second jet in (3) determines which side of the horizon follows.

This has an expression entirely in the native seat variables. From

\[
a=-\sinh\ell_1,\quad b=-\sinh\ell_2,\quad
d=\cosh\ell_1\cosh\ell_2\sin t,
\]

at \(t=0,\ell_2=\ell_1\) or \(t=\pi,\ell_2=-\ell_1\),

\[
\boxed{\xi_0=
\frac{\dot\ell_1-\sigma\dot\ell_2}{\cosh\ell_1\,\dot t}}
\qquad(\dot t\ne0).
\tag{4}
\]

The load therefore selects the pinch normal through the relative depth rate measured against angular separation. Multiplying all rates by one positive scalar changes duration, while this ratio remains invariant.

## 2. A load preserving the pinch induces a Riccati equation

Suppose \(P(a,0,0)=Q(a,0,0)=0\). Write its transverse linearization as

\[
\binom{P}{Q}=
\begin{pmatrix}m_{11}&m_{12}\\m_{21}&m_{22}\end{pmatrix}
\binom{x}{d}+O((|x|+|d|)^2).
\]

Substitution of \(x=\xi d\) into (2) gives a smooth local lift for finite \(\xi\), whose restriction to the resolved pinch is

\[
\boxed{\dot\xi=m_{12}+(m_{11}-m_{22})\xi-m_{21}\xi^2,\qquad
\dot d=d(m_{21}\xi+m_{22})+O(d^2).}
\tag{5}
\]

The roots of the first equation are the finite eigendirections of the load's transverse linearization. Their stability follows from its derivative. This supplies a calculable replacement for selecting an approach direction by hand.

Conversely, a smooth lift of the original-time vector field over an open interval of the pinch fiber requires \(P_0-\xi Q_0=0\) for every \(\xi\) in that interval, hence \(P_0=Q_0=0\). A transverse crossing as in (3) instead selects one fiber point along its trajectory. These are distinct evolution classes. Under local ODE uniqueness, an invariant pinch cannot be reached from outside in finite regular evolution time: uniqueness would continue the trajectory within the invariant set backwards as well.

Both native null directions \(\xi=\pm1\) are invariant under (5) precisely when

\[
m_{11}=m_{22},\qquad m_{12}=m_{21}.
\tag{6}
\]

For example, the declared passive relaxation

\[
\dot x=-\kappa x+\beta d,\qquad
\dot d=\beta x-\kappa d,\qquad \kappa>|\beta|,
\]

has storage \(H=(x^2+d^2)/2\), with

\[
-\dot H=\frac{\kappa-\beta}{2}(x+d)^2+
\frac{\kappa+\beta}{2}(x-d)^2>0
\]

away from the origin. It yields \(\dot\xi=\beta(1-\xi^2)\). For \(\beta>0\), fixed \(a\), and \(-1<\xi(0)<1\),

\[
\xi(T)=\tanh\!\left(\beta T+\operatorname{artanh}\xi(0)\right),
\]

so the null direction \(+1\) is approached asymptotically as the frame approaches the pinch. Equation (6) identifies the load coefficients responsible for that behavior; passivity alone does not impose them.

## 3. What must be supplied by the native load

The live [primitives](/home/williaml/seated-root/docs/2026-09-04-PRIMITIVES-v0.md:113) specify the proposed effort/flow relation and a state-dependent impedance. The [existing response runner](/home/williaml/seated-root/thm_k_response_map.py) constructs a constitutive map for differential forms; it does not specify the state field \((A,P,Q)\) or transverse area transport.

There is an exact non-uniqueness certificate. For any real \(u\), the response

\[
M_u=\begin{pmatrix}1+u^2&u\\u&1\end{pmatrix}
\]

is symmetric positive definite, has determinant one, and obeys

\[
e^TM_ue=e_x^2+(ue_x+e_d)^2.
\]

The same effort \(e=(0,1)^T\) gives \((P,Q)=(u,1)\), with the same input power \(e^TM_ue=1\). Yet (3) gives \(\xi_0=u\): timelike, null and spacelike pinch outcomes are all available. Consequently the response's coupling coefficients, beyond a scalar work or impedance measurement, contain essential prediction data.

The effort/flow framework is registered as external dependency **EXT-004**. Its use here is conditional; the power calculation and the non-uniqueness result above are direct proofs.

**Cella source reconciliation.** The [existing DIS constitutive programme](</home/williaml/Cella Framework/research/derivations/DIS_chapters/10_bipolar_axes_and_constitutive_resistance.md:337>) already makes resistance depend on environment and material, with evolving internal state and possible history dependence discussed in the transcript. William has confirmed that these cited chapters are legitimate DIS physics research. They supply an existing constitutive direction to develop into the native response, with their state, work and propagation relations made explicit. The general crossing theorem above applies to its twice differentiable trajectories; the instantaneous passive matrix is one illustrative response class. The [Cella reconciliation](/home/williaml/seated-root/docs/2026-09-05-cella-dependency-reconciliation.md) gives the enlarged-state rate formula and a channel-based area calculation.

## 4. Deriving the required area evolution in the clock frame

Let \(C\) be the transported unit native clock and \(R_n\) its orthogonal unit radial normal near the horizon. The native normal and future null normals are

\[
K=C+wR_n,\qquad w=\sqrt{1-F},\qquad
L=\frac{1+w}{2}(C+R_n),\qquad
N=\frac{C-R_n}{1+w}.
\]

Let \(\alpha=\theta_C\) and \(\beta_A=\theta_{R_n}\) be the two normal area rates. Linearity of the surface expansion gives

\[
\boxed{D:=\theta_K=\alpha+w\beta_A,\qquad
\rho=-\theta_N=\frac{\beta_A-\alpha}{1+w}.}
\tag{7}
\]

Combining this with the regular coincidence theorem gives the required clock evolution:

\[
\boxed{\alpha=-w\beta_A+F\Xi,\qquad
1+\frac{2\Xi_H}{\rho_H}>0.}
\tag{8}
\]

This displays the spatial transport contribution explicitly. In local transported surface labels it is the area-form equation \(\mathcal L_K\mu=F\Xi\mu\). A dynamical coupling must specify its trace source; trace-free shape evolution remains available.

There is also a focusing compatibility law. Assume a cut-compatible, twist-free pregeodesic outgoing congruence, \(\nabla_LL=\varkappa L\). Its optical matrix satisfies

\[
\nabla_LB=\varkappa B-B^2-\mathcal R.
\]

This follows by commuting the two covariant derivatives of \(L\) along transported separation vectors and projecting onto the screen. Taking the trace, with \(B=\theta_LI/2+\sigma_L\), gives

\[
L\theta_L=\varkappa\theta_L-\theta_L^2/2
-\|\sigma_L\|^2-\operatorname{Ric}(L,L).
\]

On a coincident horizon \(\xi=\tau=\pm1\), equation (8) and the exact native derivative of \(F\) then imply

\[
\boxed{L\xi\big|_H=
\frac{\tau(1+a_H^2)}{\rho_H+2\Xi_H}
\left(\|\sigma_L\|^2+\operatorname{Ric}(L,L)\right)_H.}
\tag{9}
\]

This is a testable compatibility condition on the state response, area response and transport curvature. The connection and propagation assumptions are recorded in **EXT-005**; no gravitational field equation is used to turn curvature into a matter source here.

## 5. A general passive response condition for coincidence

For this section, declare a local linear load whose generalized flow contains the convected area rate \(D\), with a power-conjugate area effort. Let

\[
\mathcal I=M\mathcal E,\qquad M=R+J,\qquad
R=R^T\succeq0,\quad J=-J^T.
\]

Assume smooth coefficients on both sides of a regular \(F=0\), freely variable efforts, and \(\rho_H>0\). This power port for \(D\) is a declared coupling hypothesis; deriving it from the native load is part of EXT-004.

**Passive boundary theorem.** The condition \(D|_{F=0}=0\) for every effort is equivalent to the local factorization

\[
\boxed{M=S_F(\bar R+\bar J)S_F,\qquad
S_F=\operatorname{diag}(I,F),\quad
\bar R\succeq0,\quad\bar J=-\bar J^T,}
\tag{10}
\]

with smooth barred matrices. In particular, the area self-response vanishes at least quadratically:

\[
\boxed{M_{AA}=F^2r,\qquad r\ge0.}
\tag{11}
\]

**Proof.** All-effort matching makes the area row of \(M\) zero at the boundary. Thus \(R_{AA}=0\). Positive semidefiniteness gives \(R_{Ai}=0\) there, so also \(J_{Ai}=0\). Smooth division gives one factor of \(F\) in each cross entry. The smooth nonnegative scalar \(R_{AA}\), vanishing on a two-sided regular hypersurface, also has zero first normal derivative; it has two factors of \(F\). This gives (10). Off the boundary, congruence preserves positivity of \(\bar R=S_F^{-1}RS_F^{-1}\); continuity extends positivity to the boundary. The converse follows directly by multiplication. ∎

Write the area row as \(D=F(c\cdot\mathcal E_0+Fr\mathcal E_A)\). Regular oriented coincidence for a particular load additionally requires

\[
1+2(c_H\cdot\mathcal E_{0,H})/\rho_H>0.
\]

If that orientation must hold for all unbounded efforts at fixed geometric data, then \(c_H=0\) is necessary. Equivalently the whole area row vanishes to second order, giving \(D=F^2\Psi\) and boundary factor one. Restricted efforts require only the displayed inequality.

These conditions classify universal matching. A particular applied load can have \(D=0\) with a full-rank response, as the following example demonstrates. The quadratic conclusion specifically depends on two-sided smooth passivity and the all-effort requirement.

## 6. A solved finite-time pinch and trapping crossing

Use dimensionless evolution units and a declared constant \(k>0\). Choose \(a=0\), let \(u(d)=1+2kd\), and prescribe the reciprocal load

\[
\begin{pmatrix}\dot x\\\dot d\\D\end{pmatrix}
=
\begin{pmatrix}1+u(d)^2&u(d)&0\\u(d)&1&0\\0&0&m_A\end{pmatrix}
\begin{pmatrix}0\\1\\0\end{pmatrix},\qquad m_A>0.
\tag{12}
\]

Its power for arbitrary efforts is \(e_x^2+(ue_x+e_d)^2+m_Ae_A^2\), and its determinant is \(m_A>0\). Thus the response is strictly passive, smooth and invertible. The selected effort gives \(\dot d=1,\dot x=1+2kd,D=0\). Starting at the pinch,

\[
\boxed{d(T)=T,\quad x(T)=T+kT^2,\quad
\xi(T)=1+kT,\quad F(T)=1-(1+kT)^2.}
\tag{13}
\]

**Resolved constitutive continuation.** [The Cella weighted-response extension](/home/williaml/seated-root/docs/2026-09-05-cella-constitutive-divisors.md), Theorem 3, now generates this same path from the constant resolved matrix \(\widehat M=\left(\begin{smallmatrix}1+k^2&k\\k&1\end{smallmatrix}\right)\) on \((\xi,d)\). Its pushforward differs from (12) for other efforts and stays smooth as a full response on the resolved state. [The implicit-cut calculation](/home/williaml/seated-root/docs/2026-09-05-cella-normal-geometry.md) then computes both expansions and the full Ricci tensor below directly from the selected metric and cut jets.

Equation (1) gives a real smooth native frame near \(T=0\), with original Gram determinant \(-T^2\). The signed volume changes sheet there. The resolved frame remains nondegenerate.

Take the clock and radial fields \(C=\partial_T\), \(R_n=\partial_r\). Locally \(\xi>0\), so \(w=\xi\). In spherical surface geometry, equation (7) with \(D=0\) becomes

\[
\partial_T\mathcal R+\xi(T)\partial_r\mathcal R=0.
\]

With initial area radius \(\mathcal R(0,r)=r\), its solution and metric are

\[
\boxed{\mathcal R=r-T-kT^2/2,\qquad
g=-dT^2+dr^2+\mathcal R^2d\Omega^2.}
\tag{14}
\]

The runner constructs the soldering explicitly and checks that it gives (14). Work on a regular angular chart with \(\mathcal R>0\) and \(|kT|<1/2\). The metric and native clock extend smoothly through \(T=0\).

Direct surface variation gives

\[
\boxed{\theta_L=F/\mathcal R,\qquad
\theta_N=-2/\mathcal R,\qquad
N\theta_L\big|_{T=0}=-k/r<0.}
\tag{15}
\]

Consequently the cuts are untrapped before the crossing, future marginal at the crossing, and future trapped after it. The horizon tube \(T=0\) is spacelike and future outer. This occurs in finite clock time at the native pinch.

Independent computation of the metric's complete Ricci tensor gives

\[
\operatorname{Ric}_{TT}=2k/\mathcal R,\quad
\operatorname{Ric}_{rr}=0,\quad
\operatorname{Ric}_{\theta\theta}=\xi^2-k\mathcal R,\quad
\operatorname{Ric}_{\phi\phi}=(\xi^2-k\mathcal R)\sin^2\theta.
\]

For every null vector \(V\), writing \(|V_\Omega|^2\) for its squared norm in the unit-sphere metric,

\[
\operatorname{Ric}(V,V)
=\frac{2k}{\mathcal R}(V^r)^2+
(\xi^2+k\mathcal R)|V_\Omega|^2\ge0.
\tag{16}
\]

Thus null convergence holds throughout this patch. The outgoing normal has nonaffinity \(k/2\), zero shear, and satisfies the focusing equation. At the horizon, (9) reproduces \(L\xi=k\) exactly.

For \(k=1,r=2\), the outgoing expansions at \(T=-1/4,0,1/4\) are respectively \(14/71,0,-18/55\), while the incoming expansions remain negative. These values are exact results from the solution, not fitted data.

## What this closes and what determines the prediction

The native rate-selection theorem, lifted Riccati equation, clock-frame area law, passive boundary classification and finite-time coupled example are established under their stated hypotheses. The example provides a fully specified continuation to assess, including a curvature and trapping calculation.

To select a physical prediction from the model, the remaining native calculation is now concrete: construct its actual effort and response, extract \(P,Q\) and their derivatives at a crossing or the transverse matrix \(m\) at an invariant pinch, and derive the area response. Equations (3), (5), (8) and (9) then determine the outcome. Equation (12) is an existence witness and candidate for that calculation; it is not asserted to be uniquely implied by the primitives.

The Cella/DIS physics programme supplies an existing constitutive direction involving internal variables, stored energy and response reduction. Develop its response law in the native variables and obtain area variation from the transported cut metric. Where internal variables are eliminated, establish that the reduction preserves the required evolution and work. Cella's trace channels and its signed directional curvature extension provide the corresponding area and curvature calculations; the exact derivations and source pointers are in the [reconciliation](/home/williaml/seated-root/docs/2026-09-05-cella-dependency-reconciliation.md).

Verification: the [executable derivation](/home/williaml/seated-root/pinch_load_transport.py) passes **63 exact symbolic checks**, recorded in [the results](/home/williaml/seated-root/docs/pinch-load-transport-checks.json). These include the native angle-rate formulas, both flow classes, the passive response identities, an explicit soldering, the full Ricci tensor and the finite-time trapping transition. External inputs and their rederivation targets are retained in the [first-principles debt ledger](/home/williaml/seated-root/docs/EXTERNAL-MATHEMATICS-DEBT.md). The effort/flow convention is discussed in [Bartel et al.](https://arxiv.org/html/2301.02024v1); the trapping comparison uses [Hayward's definitions](https://arxiv.org/abs/gr-qc/9303006). All load-selection and response-factorization proofs used here are stated above.

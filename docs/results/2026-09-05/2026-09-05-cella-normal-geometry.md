# Cella implicit channels: full normal geometry and area classification

5 September 2026.

**Result.** Cella's implicit-Hessian calculus extends to every regular, nondegenerate codimension-\(k\) submanifold in a fixed-signature ambient metric. It constructs the second fundamental form, every normal shape operator and its elementary polarizations, the area covector, and the normal connection along the submanifold. The total objects are invariant under \(\phi\mapsto A(x)\phi\). The original diagonal/off-diagonal channels have an explicit, generally nonzero gauge transfer. For a spacelike codimension-two cut, the area covector gives the complete pointwise trapping classification, including nonzero null mean curvature and zero mean curvature.

The construction is conditional on a specified ambient metric and implicit cut. It provides an actual derivative readout of area transport once the native cut and soldering are supplied. Unit-determinant normalization of a seat Gram matrix does not supply those derivatives or an area-conservation law.

## 1. Source, conventions and domain

The internal sources are [Mean Curvature Decomposition](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/geometric_fault_localization_and_decomposition/Mean Curvature Decomposition.md>) and [Three-Channel Mathematical Extension Notes](</home/williaml/Cella Framework/Papers_Library/05_expository_companions_and_research_maps/dbp_role_channel_and_orbit_geometry/Three_Channel_KG_New_Math_Extension.md>), especially sections 1–9. The latter already supplies every elementary curvature order for a Euclidean hypersurface, noncommutative power traces, and defining-function dependence of the raw channels. The extension here changes codimension and signature, carries the complete normal bundle, and computes the native area readout from relation jets.

Their odd-order sign conventions differ. The first source's mean curvature is the positive divergence of the unit gradient divided by the surface dimension. The extension notes use the shape convention \(S_n=-P\nabla n\), whose averaged trace is its negative. This note consistently uses

\[
B(X,Y)=(\nabla_XY)^\perp,\qquad
h(S_nX,Y)=g(B(X,Y),n),\qquad
\theta_n=-\operatorname{tr}S_n.
\tag{1}
\]

Thus an outward normal to a Euclidean round sphere has positive \(\theta_n\) and negative \(S_n\). No odd-order quantities from the two sources are identified without this sign conversion.

Let \((M^n,g)\) have a smooth nondegenerate metric of fixed signature and its Levi-Civita connection. Let \(\phi:U\to\mathbb R^k\) be smooth, \(\Sigma=\{\phi=0\}\), and \(m=n-k>0\). Along \(\Sigma\), assume

\[
\operatorname{rank}D=k,\qquad
D=d\phi,\qquad C=Dg^{-1}D^T\quad\hbox{is invertible}.
\tag{2}
\]

This is equivalent to nondegeneracy of the induced metric \(h=g|_{T\Sigma}\). Positive definiteness of \(h\) is needed only for the spacelike-cut interpretation below. A null hypersurface with \(C=0\) is outside this construction; the Lorentzian normal two-plane of a spacelike cut remains nondegenerate even when one chosen normal becomes null.

Standard linear algebra, metric volume density and Levi-Civita calculus are accepted with these hypotheses. Their definitions are not additional native physical laws. No new online source is used.

## 2. Reconstruction from the relation jet

Define the gradient frame, its conormal-dual normal frame, and the tangent projector by

\[
J=g^{-1}D^T,\qquad W=JC^{-1},\qquad P=I-WD.
\tag{3}
\]

Then \(DW=I_k\), \(P^2=P\), \(P^Tg=gP\), and

\[
T\Sigma=\ker D=\operatorname{im}P,\quad
N\Sigma=\operatorname{im}W,\quad
W^TgW=C^{-1}.
\]

Write \(\mathcal H^a=\nabla^2\phi^a\) and \(E^a=\mathcal H^a|_{T\Sigma\times T\Sigma}\). For tangent vectors \(X,Y\), the complete second fundamental form is

\[
\boxed{B(X,Y)=-W E(X,Y).}
\tag{4}
\]

**Proof.** Since \(Y\phi=0\) along the cut,

\[
0=X(Y\phi)=\nabla^2\phi(X,Y)+D(\nabla_XY).
\]

The normal component is \(WD\nabla_XY\), giving (4). This fixes the entire normal-valued bilinear form, including shear and mixed tangent entries; taking its trace is only one contraction.

For any tangent frame represented by an \(n\times m\) matrix \(T\), let \(h=T^TgT\) and \(E^a=T^T\mathcal H^aT\). If \(n=Wv=J\lambda\), then \(v=Dn\) and \(\lambda=C^{-1}Dn\). Therefore

\[
\boxed{S_n=-h^{-1}\sum_a\lambda_aE^a.}
\tag{5}
\]

Equivalently, the ambient endomorphism

\[
\mathscr S_n=-Pg^{-1}\left(\sum_a\lambda_a\mathcal H^a\right)P
\tag{6}
\]

equals \(S_n\) on tangents and vanishes on normals. Hence

\[
\det(I_n+z\mathscr S_n)=\det(I_m+zS_n),\qquad
\operatorname{tr}(\mathscr S_{n_1}\cdots\mathscr S_{n_r})
=\operatorname{tr}(S_{n_1}\cdots S_{n_r}).
\tag{7}
\]

These formulas do not require a simultaneous eigenbasis. For indefinite \(h\), self-adjoint operators need not have a real diagonalization; the polynomial identities still apply.

Set

\[
t_a=\operatorname{tr}_h E^a, \qquad
H=\operatorname{tr}_h B=-Wt.
\]

The volume/area covector and its metric dual are

\[
\boxed{\alpha(n)=\theta_n=t^TC^{-1}Dn,\qquad
H_\theta=\alpha^\sharp=Wt=-H.}
\tag{8}
\]

Here \(\frac12\mathcal L_ng|_{T\Sigma}=-g(B,n)\), so \(\theta_n\,d\mu_h\) is the first variation of the induced density \(d\mu_h=\sqrt{|\det h|}\,d^my\). A tangential deformation additionally contributes its tangential divergence. Equation (8) uses actual derivatives of the cut and ambient metric, not a determinant normalization of an unrelated frame.

The normalized conormal volume form is

\[
\omega_\perp=
\frac{d\phi^1\wedge\cdots\wedge d\phi^k}{\sqrt{|\det C|}}.
\tag{9}
\]

It has squared norm \(\operatorname{sgn}\det C\) and fixes the oriented normal volume once an orientation is chosen. In adapted coordinates \((\phi,y)\), the determinant identity

\[
|\det g_{(\phi,y)}|=\frac{|\det h|}{|\det C|}
\]

recovers the same induced density. This is the coarea factor for a nondegenerate level set; it does not impose any evolution of that density.

## 3. Exact defining-function and channel gauge laws

Let \(\widetilde\phi=A(x)\phi\) with \(A(x)\in GL(k,\mathbb R)\). On \(\Sigma\),

\[
\begin{aligned}
\widetilde D&=AD,& \widetilde C&=ACA^T,&
\widetilde W&=WA^{-1},& \widetilde P&=P,\\
\widetilde{\mathcal H}^{a}&=A^a{}_b\mathcal H^b+G^a,&
G^a&=dA^a{}_b\otimes d\phi^b+d\phi^b\otimes dA^a{}_b.
\end{aligned}
\tag{10}
\]

The omitted \(\phi^b\nabla^2A^a{}_b\) term vanishes on the cut. Since \(G|_{TT}=0\),

\[
\widetilde E=AE,\quad \widetilde t=At,\quad
\widetilde B=B,\quad\widetilde S_n=S_n,\quad
\widetilde\alpha=\alpha,\quad
\widetilde\omega_\perp=\operatorname{sgn}(\det A)\omega_\perp.
\tag{11}
\]

The frame coefficients of a fixed normal satisfy \(\widetilde v=Av\), \(\widetilde\lambda=A^{-T}\lambda\). These transformation rules prove invariance without choosing a preferred defining expression.

Equations (11) hold for a fixed geometric normal. If instead the chosen normal is the normalized defining gradient of a hypersurface, a negative scalar gauge reverses that normal and hence every odd-order shape scalar. The unoriented submanifold and \(B\) still remain fixed.

To retain Cella's original channels, fix an ambient role frame and split every covariant Hessian into its diagonal part \(\mathcal H_s\) and off-diagonal part \(\mathcal H_c\). This is a declared frame-dependent linear operation. It commutes with the relation-index matrix \(A\), but it need not commute with restriction to tangents. Define \(B_b=-W E_b\), \(b\in\{c,s\}\), and

\[
Q=-WA^{-1}(G_s|_{TT}).
\]

Then the exact channel law is

\[
\boxed{\widetilde B_s=B_s+Q,\qquad
\widetilde B_c=B_c-Q.}
\tag{12}
\]

For a fixed geometric normal \(n\), let \(R_n\) satisfy \(h(R_nX,Y)=g(Q(X,Y),n)\). Thus

\[
\widetilde S_{n,s}=S_{n,s}+R_n,\quad
\widetilde S_{n,c}=S_{n,c}-R_n,
\quad
\widetilde\theta_{n,s}=\theta_{n,s}-\operatorname{tr}R_n,
\quad
\widetilde\theta_{n,c}=\theta_{n,c}+\operatorname{tr}R_n.
\tag{13}
\]

In particular, absence of a mean-curvature interaction term does **not** make each raw mean channel a defining-function invariant. For the complete channel polynomial,

\[
\boxed{
\widetilde\chi(z;t,u;n)
=\det\!\left[I+z\{tS_{n,c}+uS_{n,s}+(u-t)R_n\}\right].}
\tag{14}
\]

Setting \(t=u\) removes the transfer at every curvature order. For a pure gauge deformation of a flat hypersurface, the two quadratic pure channels coincide and the mixed channel is minus twice either one, recovering the source's gauge-null direction. For a general initial shape, the changes contain its mixed pairings with \(R_n\); they are not confined to that one-dimensional pure-gauge pattern.

A channel split defined directly on the geometric tensor \(B\) by specified tensorial role projectors can be gauge invariant. That is an additional channel definition and must not be substituted silently for the raw Hessian split.

## 4. All elementary normal polarizations

For \(0\le r\le m\), define

\[
\det(I+zS_n)=\sum_{r=0}^m z^r\Phi_r(n),
\qquad \Phi_r(n)=\sigma_r(S_n).
\tag{15}
\]

Because \(n\mapsto S_n\) is linear, \(\Phi_r\) is a homogeneous degree-\(r\) polynomial on the normal fiber. Its full symmetric polarization is

\[
\boxed{
\Phi_r(n_1,\ldots,n_r)
=\frac1{r!}[z_1\cdots z_r]
\sigma_r\!\left(\sum_{i=1}^r z_iS_{n_i}\right).}
\tag{16}
\]

In particular,

\[
\Phi_1(n)=-\alpha(n),\qquad
\Phi_2(n_1,n_2)
=\frac12\big[
\operatorname{tr}S_{n_1}\operatorname{tr}S_{n_2}
-\operatorname{tr}(S_{n_1}S_{n_2})\big].
\tag{17}
\]

For a normal frame \(n_a\) and any declared finite channel split, every normal and channel order occurs as a coefficient of

\[
\det\!\left(I+z\sum_{a,b}u_{ab}S_{n_a,b}\right).
\tag{18}
\]

There is also a direct codimension-\(k\) extension of Cella's bordered determinant. For \(\mathcal H_\lambda=\sum_a\lambda_a\mathcal H^a\),

\[
\boxed{\det(I+zS_n)=
\frac{\det\begin{pmatrix}0_k&D\\D^T&g-z\mathcal H_\lambda\end{pmatrix}}
{(-1)^k\det g\det C}.}
\tag{18a}
\]

To prove it, use the tangent/normal basis \((T,W)\), in which \(D=(0,I_k)\). Expanding the bordered determinant along its constraint rows and columns leaves \((-1)^k\det(h-zE_\lambda)\), with the same basis determinant factor at \(z=0\). Their ratio is the stated characteristic polynomial. No inverse of \(g-z\mathcal H_\lambda\) is needed; the identity remains valid where that pencil is singular. Replacing \(\mathcal H_\lambda\) by a formal normal/channel sum gives all of (18) in one ambient formula.

Multilinearity of the determinant proves (15)–(18); no commutativity assumption is present. Equation (14) gives their raw-channel gauge transformation. With unnormalized rational normals and rational metric/relation jets, all entries and coefficients are rational. Square roots enter through chosen unit-normal or volume normalizations, rather than through the reconstruction itself.

The family \(\{\Phi_r\}\) contains every elementary invariant in every normal direction. It is not asserted to classify a tuple of shape operators up to simultaneous conjugation. Ordered traces from (7) and their commutators retain additional information. The geometric carrier is \((N\Sigma,g_N,B,\nabla^\perp)\), not merely its determinant fingerprint.

## 5. Normal transport and curvature from the same data

For a tangent vector field \(X\), define

\[
\nabla_X^\perp(Wv)=W(Xv+\omega_Xv).
\]

Differentiating \(DW=I\) gives the connection directly:

\[
\boxed{\omega_X=D\nabla_XW=-\mathcal H(X,W).}
\tag{19}
\]

The \(ab\) entry on the right is \(-\nabla^2\phi^a(X,W_b)\). In particular, mixed tangent-normal Hessian entries carry the coefficients of normal-frame transport; tangential Hessian entries alone carry \(B\).

Under (10),

\[
\boxed{\widetilde\omega_X=A\omega_XA^{-1}-(XA)A^{-1}.}
\tag{20}
\]

This follows either by differentiating \(\widetilde W=WA^{-1}\), or by substituting the full Hessian law into (19). Normal metric compatibility is

\[
X(C^{-1})=\omega_X^TC^{-1}+C^{-1}\omega_X.
\tag{21}
\]

The curvature matrix

\[
\Omega(X,Y)=X\omega_Y-Y\omega_X
+[\omega_X,\omega_Y]-\omega_{[X,Y]}
\tag{22}
\]

transforms as \(\widetilde\Omega=A\Omega A^{-1}\). Its geometric content is

\[
\boxed{g(R^\perp(X,Y)n_1,n_2)
=g(R^g(X,Y)n_1,n_2)
+h([S_{n_1},S_{n_2}]X,Y).}
\tag{23}
\]

Here \(R(X,Y)=\nabla_X\nabla_Y-\nabla_Y\nabla_X-\nabla_{[X,Y]}\). To obtain (23), substitute \(\nabla_Xn=-S_nX+\nabla_X^\perp n\) in the ambient curvature expression and take its normal pairing. The two second-fundamental terms combine into the displayed commutator. In a flat ambient space, noncommuting shapes can therefore give nonzero normal curvature. Vanishing connection coefficients at one point do not imply flatness.

The same tangent/normal splitting gives the useful complete compatibility identities

\[
\begin{aligned}
h(R^\Sigma(X,Y)Z,U)
&=g(R^g(X,Y)Z,U)
+g(B(X,U),B(Y,Z))-g(B(Y,U),B(X,Z)),\\
(\nabla_XB)(Y,Z)-(\nabla_YB)(X,Z)
&=(R^g(X,Y)Z)^\perp.
\end{aligned}
\tag{24}
\]

In the first line each extrinsic pairing is explicitly

\[
g(B(X,U),B(Y,Z))=E(X,U)^TC^{-1}E(Y,Z).
\]

These are applications of the declared connection, not a proposed native propagation law. At a point, the second jet of \(\phi\) and first jet of \(g\) determine (3)–(21), including all shape and area data. The second jet of \(g\) additionally determines ambient curvature; (23)–(24) then determine normal and intrinsic curvature without needing third derivatives of \(\phi\). Covariant variation of \(B\) generally needs its third jet. Transport along the cut requires these fields along the chosen path; a pointwise jet does not specify their future values. Transport between different physical cuts remains a separate selection/evolution datum.

## 6. Complete Lorentzian normal-plane classification

Now let \(\Sigma\) be spacelike of codimension two in a Lorentzian manifold with signature \((-+\cdots+)\). Choose future null normals \(L,N\) with \(g(L,N)=-1\), and set

\[
\ell=\alpha(L),\qquad \nu=\alpha(N).
\]

With the normal orientation fixed by this ordered pair,

\[
H_\theta=-\nu L-\ell N,\qquad
Z=-\nu L+\ell N,
\]

\[
\boxed{\alpha(Z)=0,\qquad
g(H_\theta,H_\theta)=-2\ell\nu,\qquad
g(Z,Z)=2\ell\nu.}
\tag{25}
\]

These algebraic identities classify every case:

| Expansions | Area dual \(H_\theta\) | Area-preserving normals | Cut classification |
|---|---|---|---|
| \(\ell<0,\nu<0\) | future timelike | one spacelike line | future trapped |
| \(\ell>0,\nu>0\) | past timelike | one spacelike line | past trapped |
| \(\ell\nu<0\) | spacelike | one timelike line | opposite expansions |
| \(\ell\nu=0\), exactly one nonzero | nonzero null | the same null line as \(H_\theta\) | marginal; future if the nonzero expansion is negative, past if positive |
| \(\ell=\nu=0\) | zero | the entire normal plane | extremal/zero mean vector |

For a normal \(V=pL+qN\), \(\alpha(V)=p\ell+q\nu\). Thus the first two rows also state that every future timelike normal decreases/increases area, respectively. In the opposite-sign row there is a future timelike zero-area direction. A positive null boost \(L\mapsto bL,N\mapsto N/b\) rescales \((\ell,\nu)\mapsto(b\ell,\nu/b)\), leaving \(H_\theta,Z\) and the classification fixed. Reversing the normal orientation changes the sign of \(Z\), not its line.

When \(\alpha=0\), the zero vector \(Z\) does not select a normal direction. Zero mean curvature also does not force zero second fundamental form. In an orthonormal tangent basis, normal-valued components

\[
B^0=\begin{pmatrix}1&0\\0&-1\end{pmatrix},\qquad
B^1=\begin{pmatrix}0&2\\2&0\end{pmatrix}
\tag{26}
\]

have zero trace and nonzero shear. They are realized at the origin by the flat-space graph \((t,w,u,v)=((u^2-v^2)/2,2uv,u,v)\). Consequently an extremal cut can have nontrivial extrinsic and normal curvature.

The classification is pointwise. A trapped or marginal cut must satisfy its corresponding conditions over the entire cut. An outer marginal condition also uses the derivative of an expansion. A global event horizon requires global propagation data and is not classified by (25).

For the full shape content behind these scalars, write

\[
Q_L=-hS_L=\frac{\ell}{m}h+\sigma_L,\qquad
Q_N=-hS_N=\frac{\nu}{m}h+\sigma_N,
\qquad\operatorname{tr}_h\sigma_L=\operatorname{tr}_h\sigma_N=0.
\]

Then

\[
2\Phi_2(L,N)=\left(1-\frac1m\right)\ell\nu
-\operatorname{tr}(\sigma_L^\sharp\sigma_N^\sharp),
\]

\[
R_\Sigma=R^g_T-2\left(1-\frac1m\right)\ell\nu
+2\operatorname{tr}(\sigma_L^\sharp\sigma_N^\sharp).
\tag{27}
\]

Here \(R^g_T\) is the full scalar contraction of ambient curvature on tangent directions. The scalar classification alone does not fix the shear term or the commutator in (23).

## 7. Native application: compute the area derivative from a selected cut

The native normal block is

\[
g(K,K)=-F,\qquad g(K,J)=1,\qquad g(J,J)=0,
\qquad F=\frac{1-\gamma^2}{1-\gamma^2+a^2+b^2-2ab\gamma}.
\tag{28}
\]

Its determinant is \(-1\), including \(F=0\). Given a smooth soldered metric and two independent implicit cut functions \(\phi^1,\phi^2\) whose normal plane is \(\operatorname{span}(K,J)\), equations (5) and (8) calculate

\[
\boxed{
D_{\rm area}=t^TC^{-1}DK,\qquad
\rho=t^TC^{-1}DJ,\qquad
\theta_L=t^TC^{-1}D\left(K+\frac F2J\right)
=D_{\rm area}+\frac{\rho F}{2}.}
\tag{29}
\]

Here \(D=d\phi\) is the constraint derivative; \(D_{\rm area}\) is the scalar denoted \(D\) in the existing trapping comparison. They are different objects. Equation (29) supplies the derivative calculation behind that scalar; the last equality alone was already established in the [comparison theorem](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-rest-frame-trapping-equivalence.md).

If \(\rho>0\), then \(N=-J\) has inward expansion \(-\rho\), and (25) specializes to the existing \(F_{\rm trap}=F+2D_{\rm area}/\rho\). The jet-level matching condition at a regular native boundary is now the concrete contraction

\[
\boxed{t^TC^{-1}DK=F\Xi,\qquad
1+\frac{2\Xi_H}{\rho_H}>0.}
\tag{30}
\]

This is necessary and sufficient for the regular sign-preserving coincidence stated there. Its left-hand side can be evaluated from the selected relation jets. A physically distinguished raw channel representation further resolves it by (13), with the gauge dependence kept explicit.

The two seat-screen Grams should remain separate:

\[
S=\begin{pmatrix}1+a^2&\gamma+ab\\\gamma+ab&1+b^2\end{pmatrix},
\quad\det S=\delta=d^2;
\qquad
\widehat S=\begin{pmatrix}1+a^2&a\xi\\a\xi&F+\xi^2\end{pmatrix},
\quad\det\widehat S=1.
\tag{31}
\]

The second is the resolved ruler frame's Gram matrix on its declared resolved branch. A physical cut metric is the pullback through the tangent soldering, including its scale and derivatives. If that pullback uses a square tangent map \(E_T\), then \(h=E_T^T\widehat S E_T\) and \(\det h=(\det E_T)^2\). Thus \(\det\widehat S=1\) leaves the entire \(E_T\) area factor and its deformation to be calculated. It proves no area conservation.

### Forward calculation from an implicit cut and a resolved load

Supply the metric and cut foliation

\[
g=-dT^2+dr^2+\mathcal R(T,r)^2d\Omega^2,\qquad
\phi=(T-T_0,r-r_0),\qquad\mathcal R>0.
\tag{32}
\]

Work in a regular angular chart. The tangent metric is \(h=\mathcal R^2d\Omega^2\), while the constraint derivative directly gives \(C=\operatorname{diag}(-1,1)\) and \(W=(\partial_T,\partial_r)\). The covariant Hessians are computed from \(\nabla^2T=-\Gamma^T\), \(\nabla^2r=-\Gamma^r\); their tangent restrictions and traces are

\[
E^T=-\frac{\mathcal R_T}{\mathcal R}h,\qquad
E^r=\frac{\mathcal R_r}{\mathcal R}h,\qquad
t=\frac2{\mathcal R}\begin{pmatrix}-\mathcal R_T\\\mathcal R_r\end{pmatrix}.
\tag{33}
\]

These quantities come from differentiating (32). No expansion value enters their construction. In the native \(a=0\) branch, take

\[
\begin{aligned}
K&=\partial_T+\xi\partial_r,& F&=1-\xi^2,\\
J&=\frac{\partial_r-\partial_T}{1+\xi},&
L&=\frac{1+\xi}{2}(\partial_T+\partial_r),&
N&=\frac{\partial_T-\partial_r}{1+\xi}.
\end{aligned}
\tag{34}
\]

For \(\xi>-1\), \(L,N\) are future null normals with \(g(L,N)=-1\). Substituting (33) into the general reconstruction (8) yields

\[
\boxed{
\begin{aligned}
D_{\rm area}&=\frac{2(\mathcal R_T+\xi\mathcal R_r)}{\mathcal R},\\
\rho&=\frac{2(\mathcal R_r-\mathcal R_T)}{\mathcal R(1+\xi)},\\
\theta_L&=\frac{(1+\xi)(\mathcal R_T+\mathcal R_r)}{\mathcal R},\\
\theta_N&=\frac{2(\mathcal R_T-\mathcal R_r)}{\mathcal R(1+\xi)}=-\rho.
\end{aligned}}
\tag{35}
\]

The full rest shape operator is \(S_K=-(\mathcal R_T+\xi\mathcal R_r)I_2/\mathcal R\). Thus the computation determines the second fundamental form's readout before taking its trace.

The [resolved-response theorem](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-cella-constitutive-divisors.md), Theorem 3, supplies a compatible constant material response. Choose \(k>0\) and

\[
\widehat M=\begin{pmatrix}1+k^2&k\\k&1\end{pmatrix}\succ0,
\qquad(e_x,e_d)=(0,1),\qquad(d,\xi)|_{T=0}=(0,1).
\]

The pulled-back effort is \((0,1)\), so the resolved flow gives \(d=T\), \(\xi=1+kT\), and \(x=T+kT^2\). Supply in addition the existing metric/area profile

\[
\mathcal R=r-T-\tfrac12kT^2.
\tag{36}
\]

This profile is chosen geometric transport data; the constant material matrix does not determine it. Differentiating it in (33) gives \(\mathcal R_T=-\xi\), \(\mathcal R_r=1\), hence

\[
\boxed{D_{\rm area}=0,\qquad \rho=\frac2{\mathcal R},\qquad
\theta_L=\frac F{\mathcal R},\qquad
\theta_N=-\frac2{\mathcal R}.}
\tag{37}
\]

Across the family of implicit cuts, direct differentiation further gives

\[
N\theta_L\big|_{T=0}=-\frac{k}{r}<0,\qquad
\partial_TF\big|_{T=0}=-2k.
\]

Thus the selected metric gives regular future outer marginal cuts at \(T=0\); nearby \(T>0\) cuts are future trapped while \(\mathcal R>0\) and \(\xi>-1\). The result follows from the cut Hessians of (32), with the supplied profile (36).

The runner also imports the [full-tensor curvature construction](/home/williaml/seated-root/cella_tensor_valuation.py) and applies it to this same metric. It computes the complete Ricci tensor in the \((T,r,\theta,\varphi)\) coordinate frame:

\[
\operatorname{Ric}
=\operatorname{diag}\!\left(
\frac{2k}{\mathcal R},0,
\xi^2-k\mathcal R,
(\xi^2-k\mathcal R)\sin^2\theta\right),
\tag{38}
\]

\[
\operatorname{Ric}(L,L)=\frac{k(1+\xi)^2}{2\mathcal R},\qquad
R[g]=\frac{2\xi^2}{\mathcal R^2}-\frac{4k}{\mathcal R}.
\]

This connects the resolved constitutive trajectory, implicit normal geometry and tensor-curvature extensions in one forward calculation. It introduces no field equation selecting the metric or its area profile.

**Exact remaining datum.** The native theory must specify the physical cut/foliation and its tangent soldering, or equivalently a selected implicit cut pair \(\phi\) with its second jet together with the ambient metric's first jet. These data determine \(D_{\rm area}\), \(\rho\), all shears, all normal shape invariants and the normal connection along the cut. The seat Gram values at one point do not determine those jets. A load-induced area rate plus spatial/material drift can provide them through a native evolution law; an independently declared area power port is not required by this geometry. Neither the diagonal channel decomposition nor a normalized frame selects that law.

## 8. Verification and provenance

Run [suites/cella_normal_geometry.py](/home/williaml/seated-root/cella_normal_geometry.py). Its exact output is [cella-normal-geometry-checks.json](/home/williaml/seated-root/docs/receipts/cella-normal-geometry-checks.json). The replay compares implicit reconstruction with direct embedding derivatives for a mixed codimension-two graph in \(1+4\) dimensions, checks all cubic normal monomials, tests a nonconstant matrix defining-function gauge including nonzero channel transfer, and checks the normal-connection law from independently differentiated gauged relations. A second graph has

\[
\omega_u(0)=\omega_v(0)=0,\qquad
\Omega_{uv}(0)=\begin{pmatrix}0&-4\\-4&0\end{pmatrix},
\]

matching its noncommuting shape operators through (23). A varying, off-diagonal ambient metric verifies the covariant Hessian terms against direct Christoffel acceleration. Symbolic mixed-determinant checks do not assume diagonal shapes or a diagonal tangent metric. The zero and nonzero-null mean cases are retained explicitly.

The native fixture differentiates the actual implicit cuts (32), constructs \(t,C,W\), and evaluates their normal expansions before specializing the metric. It also verifies the constant resolved load path, the outer derivative and every Ricci component of (36) using the tensor-curvature module. The metric and area profile remain declared inputs throughout.

The live Cella DAG entries \(GFL:thm:mean_curvature_decomp\) and \(DBP:thm:three_channel_kg_ext\) were inspected as provenance indexes and their source files read. The respective source SHA-256 values match the registered values:

- \(\texttt{ace98fbfd85a92cca4da4f348c410c81fe38a0409d57840915d98fa3359f9033}\).
- \(\texttt{dcb8e175e3c43f9b708afda98e8cbf02c55f910e2d67e28b300edf676ab7c2fd}\).

The original mean-curvature source declares the standard mean-curvature formula. That provenance remains. Equations (4), (10)–(23) give the full conditional geometric extension; the finite replay checks its implementation and examples. Native physical cut selection, soldering and evolution remain the specific open choices described above. No external search or new physical framework is adopted.

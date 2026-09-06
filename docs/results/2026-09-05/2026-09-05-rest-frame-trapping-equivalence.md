# Native rest-frame loss and gravitational trapping: the comparison theorem

5 September 2026

**Result.** The native rest criterion and gravitational trapping coincide under a precise area-transport condition. Arbitrary invertible soldering does not enforce that condition. Their exact separation is

\[
\boxed{F_{\rm trap}=F+\frac{2\theta_K}{\rho},\qquad
F=\frac{1-\gamma^2}{1-\gamma^2+a^2+b^2-2ab\gamma},\qquad
\rho=-\theta_N>0.}
\tag{1}
\]

Here \(K\) is the transported native rest normal, \(\theta_K\) measures its transverse area transport, and \(N\) is the future inward null normal. This gives both a general comparison theorem and a concrete target for the model's dynamics. The strong matching condition \(\theta_K=0\) admits changing transverse shape and angular drift; spherical symmetry is unnecessary.

**Dependency ledger:** derived from the current nondegenerate Lorentzian seat form, its fourth positive direction, a smooth invertible soldering, a compatible foliation by spacelike two-surfaces, and inward contraction. No gravitational field equation or native area-transport law is assumed. “Trapping” below concerns these surfaces and their null expansions; identifying a global event horizon would require global causal data.

## 1. Native normal and admissible surfaces

Use the current T7 form, extended by the positive T8 direction:

\[
G=\begin{pmatrix}-1&a&b\\a&1&\gamma\\b&\gamma&1\end{pmatrix},
\qquad \delta=-\det G>0,
\qquad k=-G^{-1}(1,0,0)^T.
\]

With \(c\) the native clock vector, \(q(k,c)=-1\) and \(q(k,k)=-F\). The resolved ruler plane is \(\mathcal P=\operatorname{span}(h,g,e_4)=k^\perp\). A unit timelike normal to this plane exists exactly when \(F>0\); its limiting normal is null at \(F=0\), then spacelike for \(F<0\). This is the native rest-frame criterion being compared.

For any smooth invertible soldering \(\mathsf E:TM\to\mathcal V\), set

\[
g=\mathsf E^*q,\qquad K=\mathsf E^{-1}k.
\]

Then \(g(K,K)=-F\) for every such map. Near either native horizon, the previously derived regular null partner \(\ell=(e_r-c)/(1+\sqrt{1-F})\) satisfies \(q(k,\ell)=1\). Its image \(J=\mathsf E^{-1}\ell\) gives

\[
g(K,K)=-F,\quad g(K,J)=1,\quad g(J,J)=0.
\tag{2}
\]

Assume the positive screen distribution \(\mathscr S=\operatorname{span}(K,J)^\perp\) integrates to the selected two-surfaces. Explicitly, Frobenius requires \(dK^\flat|_{\mathscr S}=dJ^\flat|_{\mathscr S}=0\). This concerns the two-dimensional screen, not integrability of the entire ruler plane.

The future outgoing and inward null normals are

\[
L=K+\frac F2J,\qquad N=-J,\qquad g(L,N)=-1.
\tag{3}
\]

For a normal vector \(V\), define \(\theta_V=\tfrac12 h^{AB}(\mathcal L_Vg)_{AB}\), with \(h\) the induced surface metric. Write \(D=\theta_K\) and \(\rho=\theta_J\), restricting to \(\rho>0\). Future trapping means both future null expansions are negative; future marginality means \(\theta_L=0,\theta_N<0\). Conditions describing a trapped or marginal surface must hold over the entire cut. These are the local definitions used in [Hayward's trapping-horizon framework](https://arxiv.org/abs/gr-qc/9303006).

## 2. Exact bridge and the geometric rest normal

**Theorem 1.** For every admissible soldering and selected surface,

\[
\boxed{\theta_L=D+\frac{\rho F}{2},\qquad \theta_N=-\rho.}
\tag{4}
\]

**Proof.** Expansion is linear over smooth coefficients for normal fields: in the screen projection of \(\mathcal L_{fV}g\), derivatives of \(f\) multiply components of \(V^\flat\) that vanish on the screen. Apply this to (3). Thus \(2\theta_L/\rho=F_{\rm trap}\), proving (1). Its positive, zero and negative signs respectively describe an expanding outgoing direction, future marginality and future trapping under the stated inward-contraction condition. ∎

There is a useful stronger geometric interpretation. Among normal vectors satisfying \(g(V,J)=1\), the unique vector transporting zero area is

\[
\boxed{K_{\rm area}=K-\frac D\rho J,\qquad
\theta_{K_{\rm area}}=0,\qquad
g(K_{\rm area},K_{\rm area})=-F_{\rm trap}.}
\tag{5}
\]

Indeed, the normalization forces \(V=K+tJ\); zero expansion forces \(t=-D/\rho\); (2) gives its norm. Consequently trapping has its own rest-normal interpretation. In this orientation the dual mean-curvature vector is

\[
Z=-\theta_NL+\theta_LN=\rho K-DJ=\rho K_{\rm area},
\qquad g(Z,Z)=-\rho^2F_{\rm trap}.
\]

The area-preserving normal and its relation to marginal surfaces belong to established extrinsic geometry; see [Anco's treatment of dual mean curvature](https://arxiv.org/abs/gr-qc/0402057). This supplies a general mathematical framework for the comparison. The native research question is whether the seat/load law makes the independently defined \(K\) agree with this area normal, or makes their causal transitions agree. Equation (5) does not replace the native normal by definition.

Under the positive rescaling \(K\mapsto bK,J\mapsto J/b\), the quantities \(F,D,\rho,F_{\rm trap}\) transform as \(b^2F,bD,\rho/b,b^2F_{\rm trap}\). The classifications and coincidence of their zero sets are invariant.

## 3. Necessary and sufficient conditions for regular coincidence

**Theorem 2.** Let \(H_0=\{F=0\}\) be regular, with \(dF\ne0\), and let \(\rho>0\). Locally at a point of \(H_0\), the following are equivalent:

1. \(\theta_L\) has the same zero hypersurface, is regular there, and has the same sign as \(F\) on both sides.
2. There is a smooth function \(\Xi\) such that

\[
\boxed{D=F\Xi,\qquad
\Lambda_H:=1+\frac{2\Xi_H}{\rho_H}>0.}
\tag{6}
\]

Then, after restricting the neighborhood if necessary,

\[
F_{\rm trap}=\Lambda F,\qquad
\theta_L=\frac{\rho\Lambda}{2}F,\qquad
\Lambda=1+\frac{2\Xi}{\rho}>0.
\tag{7}
\]

**Proof.** Coincidence implies \(D|_{H_0}=0\) by (4). Since \(dF\ne0\), use \(F\) as a local coordinate. The smooth division identity

\[
D(F,y)=F\int_0^1(\partial_FD)(tF,y)\,dt
\]

gives \(D=F\Xi\). On the boundary, \(d\theta_L=(\rho/2+\Xi)dF\). Regularity makes this coefficient nonzero; equal orientation makes it positive. Conversely, (6), continuity and (4) give (7), including regularity, equal zero sets and equal signs. ∎

For a marginal tube, additionally require \(H_0\) to be a union of the selected cuts. The strong sufficient law \(D=0\) gives exact equality \(F_{\rm trap}=F\) throughout a neighborhood. Boundary equality needs only the weaker (6). Cases with \(\Lambda_H=0\) require higher-order analysis and lie outside the regular theorem; (1) still applies. The existing Cella weighted initial-form calculus is a candidate continuation tool there.

For the usual outer condition,

\[
N\theta_L|_{H_0}=\frac{\rho_H\Lambda_H}{2}\,NF|_{H_0}.
\]

At the native branch \(\gamma=\sigma\), \(\sigma=\pm1\),

\[
dF|_{H_0}=-\frac{2\sigma}{(a-\sigma b)^2}d\gamma,
\qquad
N\theta_L|_{H_0}=-\frac{\sigma\rho_H\Lambda_H}{(a-\sigma b)^2}N\gamma.
\tag{8}
\]

Thus future outer marginality translates directly into \(\sigma N\gamma>0\) under (6).

## 4. What the previous soldering supplied, and a larger matching class

For the previously derived family

\[
g=-A^2F\,dv^2+2AB\,dv\,dr+
h_{AB}(dx^A-U^Adv)(dx^B-U^Bdv),
\]

with \(A,B>0\), the native normals and area density are

\[
K=A^{-1}(\partial_v+U^A\partial_A),\qquad
J=B^{-1}\partial_r,\qquad \mu=\sqrt{\det h}.
\]

Tracing the full spacetime Lie derivative gives

\[
\boxed{D=\frac{\partial_v\mu+\partial_A(\mu U^A)}{A\mu},
\qquad \rho=\frac{\partial_r\mu}{B\mu}.}
\tag{9}
\]

The spherical choice \(h=r^2d\Omega^2,U=0\) has \(D=0\). That is the actual mechanism behind the earlier coincidence.

A larger class, on regular angular charts, is

\[
h=R(r)^2\left[e^{2s(v,r,\theta)}d\theta^2+
e^{-2s(v,r,\theta)}\sin^2\theta\,d\phi^2\right],
\qquad U^\theta=0,\quad U^\phi=\Omega(v,r,\theta).
\]

Here \(\mu=R^2\sin\theta\), so \(D=0\) and \(\rho=2R'/(BR)>0\) whenever \(R>0,R'>0\). Shape can evolve and angular drift can be nonzero while the two criteria coincide. Global sphere constructions require smooth pole-compatible choices of \(s,\Omega\). This is a proved geometric class; selecting a member remains dynamical work.

The weaker transport target can be written as the trace law

\[
(\partial_v+\mathcal L_U)h_{AB}
=2\Sigma_{AB}+AF\Xi h_{AB},\qquad h^{AB}\Sigma_{AB}=0.
\tag{10}
\]

Its trace gives exactly \(D=F\Xi\). It leaves trace-free deformation available to the dynamics.

## 5. Why general soldering needs an additional transport law

There is a direct independence proof. Keep the native state, normal metric, \(K,J\) and screen distribution fixed, but change the soldering on the screen so that \(h\mapsto e^{2\chi}h\). Then

\[
D\mapsto D+2K\chi,\qquad
\rho\mapsto\rho+2J\chi.
\tag{11}
\]

At a point choose \(\chi=0,J\chi=0\), and choose \(K\chi\) freely. The pointwise metric, native scalar and inward expansion are unchanged while \(D\) is arbitrary. Therefore a condition on the native form and pointwise soldering alone cannot force coincidence. The missing information is precisely first-derivative area transport.

Two explicit realizations locate this freedom:

**Native transition in flat spacetime.** Choose \(a=-1,b=0,\gamma=1-(r-r_0)/L_0\), with \(r_0,L_0>0\). Near \(r=r_0\), \(\delta=2-\gamma^2>0\). In native coordinates use

\[
c=(1,0,0),\quad e_r=(-1,1,-\gamma)/\sqrt\delta,\quad
e_\perp=(0,0,1).
\]

These, with \(e_4\), form an orthonormal frame. The map
\(\mathsf E\partial_t=c,\mathsf E\partial_r=e_r,
\mathsf E\partial_\theta=r e_\perp,
\mathsf E\partial_\phi=r\sin\theta e_4\)
gives exactly \(g=-dt^2+dr^2+r^2d\Omega^2\). With \(w=\sqrt{1-F}=1/\sqrt\delta\),

\[
K=\partial_t+w\partial_r,\quad
J=\frac{\partial_r-\partial_t}{1+w},\quad
\theta_L=\frac{1+w}{r},\quad
\theta_N=-\frac{2}{r(1+w)},\quad F_{\rm trap}=(1+w)^2>0.
\]

The native \(F\) crosses zero while these spheres remain untrapped. The construction is smooth and full rank at the transition. It is an admissible arbitrary soldering; a native coupling law would need to determine whether it is admissible physically.

**A null native boundary with changing area.** In (9), set \(A=B=1,U=0,F=F(r)\), with a simple zero, and choose \(h=e^{2\lambda v}r^2d\Omega^2\). The boundary is null and generated by \(K=\partial_v\), but

\[
D=2\lambda,\quad \rho=2/r,\quad
\theta_L|_{H_0}=2\lambda,\quad F_{\rm trap}=F+2\lambda r.
\]

Thus even null-generator alignment leaves an independent area-transport requirement. This example asserts geometry, without imposing a gravitational field equation.

## 6. Null alignment and the next native theorem

For regular \(H_0\), the native \(K\) generates it as a null hypersurface exactly when

\[
\boxed{dF|_{H_0}=\eta K^\flat|_{H_0},\qquad \eta\ne0.}
\tag{12}
\]

The forward implication follows from the unique null normal line of a null hypersurface; the reverse follows because the annihilator of \(K^\flat\) is its null tangent hyperplane. This is a separate geometric condition from the expansion condition (6). A marginal tube can also be spacelike: if it is cut-adapted and \(NF\ne0\), its orthogonal tangent is \(T=L-(LF/NF)N\), with \(g(T,T)=2LF/NF\).

**The highest-value continuation is now explicit:** derive the screen/cut transport and equation (10), with \(1+2\Xi_H/\rho_H>0\), from the seat/load dynamics. This would promote conditional coincidence to a native theorem while retaining transverse deformation. The particularly strong result would be the area conservation law \(\partial_v\mu+\partial_A(\mu U^A)=0\).

If the independently derived dynamics instead gives \(D_H\ne0\), equation (1) supplies a computable separation between the native rest boundary and the marginal surface. That becomes a prediction target once the same dynamics fixes its magnitude and the physical surface selection. Imposing an area law after selecting a desired horizon would not establish this result.

The general comparison is settled here. The native dynamical transport theorem remains open, with its required equation and admissibility conditions specified.

## Verification and project files

The [executable derivation](/home/williaml/seated-root/rest_frame_trapping_equivalence.py) passes **40 exact symbolic checks**, recorded in [the results](/home/williaml/seated-root/docs/receipts/rest-frame-trapping-equivalence-checks.json). They include the full arbitrary-metric expansion calculation, native flat soldering, both separating examples, the anisotropic matching class and both native horizon branches. The smooth local equivalence proof is given above; symbolic checks support its identities without claiming to prove the unresolved transport law.

This extends [the horizon-crossing metric derivation](/home/williaml/seated-root/docs/results/2026-09-05/2026-09-05-horizon-crossing-metric.md). All three new deliverables are stored in `/home/williaml/seated-root`.

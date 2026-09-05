**A horizon-crossing metric from the native Lorentzian frame**

5 September 2026

The strongest result supported by the present construction is a **canonical regular radial block with a general spacetime completion**:

\[
\boxed{
ds^2=-A^2F\,dv^2+2AB\,dv\,dr
+h_{AB}(dx^A-U^A dv)(dx^B-U^B dv),
\qquad
F=\frac{1-\gamma^2}
{1-\gamma^2+a^2+b^2-2ab\gamma}.
}
\tag{1}
\]

The coefficients can depend on all coordinates. The transverse metric \(h_{AB}\) is positive definite, and \(A,B>0\). Its determinant is

\[
\boxed{\det g=-A^2B^2\det h,}
\tag{2}
\]

which stays nonzero through \(F=0\). This includes time-dependent and transversely drifting configurations, with arbitrary transverse geometry.

The native derivation fixes \(F\) and the regular radial frame. The map from coordinate displacements into that frame is the explicit spacetime-completion input below. Its coefficients \(A,B,h_{AB},U^A\) remain available to the seat/load law. Equation (1) is therefore a derived metric family under that input, rather than a uniquely selected dynamical solution.

The [symbolic derivation runner](</home/williaml/seated-root/horizon_crossing_metric.py>) verifies the frame identities, determinant, inverse, horizon expansions, radial geodesics, surface-gravity relation, a curvature component and the scalar wave operator. Its [results](</home/williaml/seated-root/docs/horizon-crossing-metric-checks.json>) accompany this document.

**1. Native data and the regular coframe**

Use T7's seat-derived form, with selected pole representatives:

\[
G=
\begin{pmatrix}
-1&a&b\\
a&1&\gamma\\
b&\gamma&1
\end{pmatrix},\qquad
\delta=-\det G
=1-\gamma^2+a^2+b^2-2ab\gamma>0.
\]

Write

\[
\varepsilon=1-\gamma^2,\qquad
F=\varepsilon/\delta,\qquad
W=\delta-\varepsilon=a^2+b^2-2ab\gamma.
\]

On the compact-ruler side, \(F=N^2\). The rational function \(F\) continues as a signed invariant across the ruler-plane boundary. The construction below never requires taking \(\sqrt F\).

For \(\gamma=\sigma\), with \(\sigma=\pm1\),

\[
\delta_H=(a-\sigma b)^2.
\]

Thus the two boundaries are inside the domain of this derivation wherever \(a-\sigma b\ne0\). These are the existing T7 horizons away from full frame rank loss. [Native state construction](</home/williaml/seated-root/prim_t7_seat_form.py:84>), [region structure](</home/williaml/seated-root/prim_t7f_regions.py:1>)

There is already a regular factorization of the entire native metric. For the coframe dual to \((c,h,g)\), set

\[
\begin{aligned}
\vartheta^0&=\theta^c-a\theta^h-b\theta^g,\\
\vartheta^1&=\sqrt{1+a^2}\,\theta^h
 +\frac{\gamma+ab}{\sqrt{1+a^2}}\,\theta^g,\\
\vartheta^2&=\frac{\sqrt\delta}{\sqrt{1+a^2}}\,\theta^g.
\end{aligned}
\]

Direct expansion gives

\[
q=-(\vartheta^0)^2+(\vartheta^1)^2+(\vartheta^2)^2.
\tag{3}
\]

Every coefficient remains real and finite at \(\gamma=\pm1\) for finite native data and \(\delta>0\). This factorization also covers zero tilt.

**2. Derive the radial frame from the seat's dual normal**

Let \(c,h,g\) now denote the frame vectors and define

\[
\nu=G^{-1}(1,0,0)^T
=\frac1\delta
\begin{pmatrix}
-\varepsilon\\a-b\gamma\\b-a\gamma
\end{pmatrix}.
\]

Then

\[
q(\nu,c)=1,\qquad q(\nu,h)=q(\nu,g)=0,\qquad q(\nu,\nu)=-F.
\]

Choose the time orientation containing \(c\), and put \(k=-\nu\). This fixes

\[
q(k,c)=-1,\qquad q(k,k)=-F.
\tag{4}
\]

Near either horizon, \(W>0\). Define

\[
w=\sqrt{W/\delta}=\sqrt{1-F},\qquad
e_r=\frac{k-c}{w},\qquad
e_\perp=\frac{b h-a g}{\sqrt W}.
\tag{5}
\]

The rational identities behind these normalizations are

\[
\begin{aligned}
q(k-c,k-c)&=W/\delta, & q(c,k-c)&=0,\\
q(bh-ag,bh-ag)&=W, &
q(c,bh-ag)&=q(k,bh-ag)=0.
\end{aligned}
\]

Consequently \((c,e_r,e_\perp)\) is an orthonormal frame and

\[
k=c+w e_r.
\tag{6}
\]

In the basis \((k,e_r,e_\perp)\), the native form is

\[
\begin{pmatrix}
-F&w&0\\w&1&0\\0&0&1
\end{pmatrix},
\qquad \det=-1.
\tag{7}
\]

This is a regular flow form obtained from the native data. At both horizons, \(w=1\). The radial direction in (5) is undefined at zero tilt, where the seat supplies no preferred radial direction; the full coframe (3) remains available there.

Now define the null partner

\[
\boxed{\ell=\frac{e_r-c}{1+w}.}
\tag{8}
\]

It obeys

\[
q(\ell,\ell)=0,\qquad q(k,\ell)=1,\qquad
q(\ell,e_\perp)=0.
\]

Therefore, in the basis \((k,\ell,e_\perp)\),

\[
\boxed{
[q]=
\begin{pmatrix}
-F&1&0\\1&0&0\\0&0&1
\end{pmatrix}.
}
\tag{9}
\]

At the horizon, \(\ell=(e_r-c)/2\), which is finite.

There is a useful uniqueness statement here. Within \(\operatorname{span}(c,k)\), write a null partner as \(x c+y k\), impose \(q(k,\ell)=1\), and solve its null equation. The two solutions have

\[
y=\frac1{w(1+w)}
\quad\hbox{or}\quad
y=-\frac1{w(1-w)}.
\]

Only the first stays finite at \(w=1\). Thus (8) is the unique finite branch within the seat/normal plane, with the prescribed pairing. The regular cross term in (9) has been derived by this normalization.

**3. Construct the spacetime metric explicitly**

To convert a form on the cell into a spacetime metric, specify a map \(\mathsf E\) taking coordinate displacements into the cell frame. This is often called a coframe or soldering map; it is additional structure beyond a state-dependent Gram matrix.

Use the current positive fourth-generator branch for the four-dimensional version. It adds a second positive direction transverse to \((k,\ell)\). The native \(2+1\) construction works identically with one transverse direction. [Fourth-generator construction](</home/williaml/seated-root/prim_t8a_fourth_direction.py:1>)

On a local null coordinate patch, take

\[
\mathsf E
=k\,A\,dv+\ell\,B\,dr
+\sum_{I=2}^{3}e_I E^I{}_A(dx^A-U^A dv).
\tag{10}
\]

Here \(e_I\) are orthonormal positive transverse frame vectors, and

\[
h_{AB}=\sum_I E^I{}_A E^I{}_B.
\]

Define \(g(X,Y)=q(\mathsf E X,\mathsf E Y)\). Substituting (9) immediately yields (1).

This completion allows arbitrary clock and radial coefficients, transverse drift and positive transverse geometry. It does not identify the state directions that leave \(F\) fixed as gauge: their response can enter these coefficients and the transport of the native rulers.

The metric's inverse is

\[
\begin{aligned}
g^{vv}&=0,&g^{vr}&=\frac1{AB},&g^{vA}&=0,\\
g^{rr}&=\frac{F}{B^2},&
g^{rA}&=\frac{U^A}{AB},&
g^{AB}&=h^{AB}.
\end{aligned}
\tag{11}
\]

Together, (2) and (11) prove regularity through \(F=0\) when the state and coframe coefficients extend smoothly with \(\delta>0\), \(AB\ne0\) and \(\det h>0\). With twice continuously differentiable coefficients, the connection and curvature are finite in a regular coordinate patch.

Null-coordinate metric forms are established mathematics; Gaussian null coordinates provide the corresponding local description near smooth Killing horizons. The specific result here is the construction of the radial frame and its coefficient \(F\) from the native Gram data. [Fontanella, equations (1.3)–(1.5)](https://arxiv.org/pdf/2211.03861)

**4. Identify the horizon using expansions**

For a spherical spacetime, set

\[
U^A=0,\qquad h_{AB}dx^A dx^B=r^2d\Omega^2,
\]

with \(r\) the areal radius. Then

\[
\boxed{ds^2=-A(v,r)^2F(v,r)\,dv^2
+2A(v,r)B(v,r)\,dv\,dr+r^2d\Omega^2.}
\tag{12}
\]

The future null vectors

\[
L=\partial_v+\frac{AF}{2B}\partial_r,\qquad
N=-\frac1{AB}\partial_r
\]

satisfy \(g(L,N)=-1\). Their spherical expansions are

\[
\boxed{\theta_L=\frac{AF}{Br},\qquad
\theta_N=-\frac2{ABr}.}
\tag{13}
\]

For \(A,B,r>0\):

- \(F>0\): outgoing expansion is positive, ingoing expansion is negative.
- \(F=0\): the sphere is future marginally trapped.
- \(F<0\): both expansions are negative.

The outer-horizon condition follows directly:

\[
N(\theta_L)\big|_H
=-\frac{F_r}{B_H^2r_H}.
\tag{14}
\]

Thus \(F_r>0\) gives the future outer branch in this orientation.

Time-dependent marginal surfaces are included. If \(F(v,r_H(v))=0\) and \(F_r\ne0\),

\[
\frac{dr_H}{dv}=-\frac{F_v}{F_r}
=-\frac{\gamma_v}{\gamma_r}\bigg|_H.
\tag{15}
\]

This last equality follows from the native horizon derivative in the next section.

For the nonspherical family, the zero set of \(F\) is null precisely when

\[
\left.
\left[
\frac{2F_r}{AB}(\partial_v+U^A\partial_A)F
+h^{AB}F_A F_B
\right]\right|_{F=0}=0,
\tag{16}
\]

assuming \(dF\ne0\). This is the explicit compatibility condition for a moving or distorted boundary. In coordinates adapted to a fixed null boundary \(r=r_H\), its tangential derivatives vanish and (16) holds. In spherical dynamics, (13) defines the marginal spheres irrespective of whether their evolving hypersurface is null. Global event-horizon identification uses the completed spacetime evolution.

**5. Surface gravity in native variables**

For the stationary spherical sector, normalize the Killing field as \(K=\partial_v\). On \(F=0\),

\[
K^\flat=AB\,dr,\qquad K^2=-A^2F.
\]

Using \(d(K^2)|_H=-2\kappa_H K^\flat|_H\) gives

\[
\kappa_H=\frac{A_H}{2B_H}F_r|_H.
\tag{17}
\]

The native quotient has a particularly simple horizon derivative:

\[
\boxed{
dF|_{\gamma=\sigma}
=-\frac{2\sigma}{(a-\sigma b)^2}\,d\gamma.
}
\tag{18}
\]

All terms involving derivatives of \(a,b\) vanish at this order because the numerator \(1-\gamma^2\) vanishes. Their horizon values supply the finite depth gap in the denominator.

Combining (17) and (18) gives

\[
\boxed{
\kappa_H
=-\frac{\sigma A_H\,\gamma_r|_H}
{B_H(a_H-\sigma b_H)^2},
\qquad \gamma_H=\sigma=\pm1.
}
\tag{19}
\]

This is a concrete connection between the native state trajectory and horizon geometry. Once the load law determines the crossing slope and the clock/radial response, it determines the surface gravity in the chosen clock normalization. A change in the normalization of \(K\) rescales \(\kappa\) accordingly.

Higher-order horizons also fit the construction. If \(\delta_H>0\) and
\(\gamma=\sigma(1-c\rho^p+o(\rho^p))\), \(c\ne0\), then

\[
F=\frac{2c}{(a_H-\sigma b_H)^2}\rho^p+o(\rho^p).
\]

Odd \(p\) changes the sign of \(F\); even \(p\) touches the boundary. For stationary data, \(p>1\) gives \(\kappa_H=0\). The metric's nondegeneracy is independent of that order.

**6. An infalling geodesic crosses in finite proper time**

In the stationary spherical sector, let a radial timelike geodesic have positive conserved energy

\[
E=A^2F\dot v-AB\dot r,
\]

where a dot denotes differentiation with respect to proper time and the norm is \(-1\). Solving the two equations for the ingoing branch gives

\[
\boxed{
\dot r=-\frac{\sqrt{E^2-A^2F}}{AB},\qquad
\dot v=\frac1{E+\sqrt{E^2-A^2F}}.
}
\tag{20}
\]

The second expression is the regularized exact solution; it contains no division by \(F\). At the horizon,

\[
\dot r_H=-\frac{E}{A_HB_H},\qquad
\dot v_H=\frac1{2E}.
\]

Both are finite, and the radial velocity is nonzero. Smooth coefficients therefore give a local continuation through the horizon in finite proper time.

**7. Differential quantities ready for a load law**

The geometric scalar wave operator of (12) is

\[
\boxed{
\square_g\psi=
\frac{
\partial_v(r^2\partial_r\psi)
+\partial_r\!\left(r^2\partial_v\psi+
\frac AB r^2F\,\partial_r\psi\right)}
{ABr^2}
+\frac1{r^2}\Delta_\Omega\psi.
}
\tag{21}
\]

At \(F=0\), its mixed principal term remains \(2(AB)^{-1}\partial_v\partial_r\). This supplies a regular geometric operator for minimally coupled scalar propagation experiments on the derived family.

A direct four-dimensional Ricci calculation also gives the horizon component of the Einstein tensor:

\[
\boxed{
G_{vv}|_H
=-\frac{A_H}{B_Hr_H}F_v|_H
=\frac{2\sigma A_H\,\gamma_v|_H}
{B_Hr_H(a_H-\sigma b_H)^2}.
}
\tag{22}
\]

Equation (22) is a geometric identity. A proposed source-response equation can be matched to it to determine horizon motion; no source equation was used to obtain it.

**8. A concrete two-boundary member**

Take a dimensionless radial variable \(\rho=(r-r_0)/L\), \(L>0\), and choose the explicit native state path

\[
a=-1,\qquad b=0,\qquad \gamma=1-\rho.
\]

It gives

\[
\delta=1+2\rho-\rho^2,\qquad
F=\frac{\rho(2-\rho)}{1+2\rho-\rho^2}.
\tag{23}
\]

The native frame has full rank on
\(1-\sqrt2<\rho<1+\sqrt2\). Both \(\rho=0\) and \(\rho=2\) lie inside this interval and satisfy \(\delta=1\).

For the declared unit clock and radial-response specialization \(A=B=1\), the spherical metric is explicitly

\[
ds^2=
-\frac{\rho(2-\rho)}{1+2\rho-\rho^2}\,dv^2
+2\,dv\,dr+r^2d\Omega^2.
\tag{24}
\]

Choose \(r_0\) so the patch has \(r>0\). The determinant stays \(-r^4\sin^2\theta\) on a regular angular chart through both boundaries. Their stationary surface gravities in this normalization are \(+1/L\) and \(-1/L\). This is an exact crossing example generated by a native state path, with its path and response coefficients explicitly chosen.

**What this derivation fixes next**

The regular radial structure is now explicit, and the metric family supports horizon crossing, trapped-surface calculations and wave propagation. The remaining dynamical calculation is to obtain \(A,B,h_{AB},U^A\) and the native state trajectory from the seat/load interaction.

In the spherical sector, that calculation has immediate quantitative targets: \(A/B\), the depth gap \(a-\sigma b\), and the derivatives \(\gamma_r,\gamma_v\). Equations (19), (15) and (22) convert them directly into surface gravity, horizon motion and a geometric flux component.

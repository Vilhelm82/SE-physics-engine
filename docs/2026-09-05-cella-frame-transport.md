# Cella accounts, regular native transport, and all local solderings

5 September 2026

**The finite-direction native frame admits an explicit regular presentation transport, a defining-function-gauge-invariant Cella jet/edge map, and a complete local classification of compatible connections and solderings.** Put

\[
p=1+a^2,\qquad F=\frac{1-\xi^2}{p},\qquad
\widehat G=
\begin{pmatrix}-1&a&\xi\\a&1&0\\\xi&0&F\end{pmatrix},
\qquad \eta=\operatorname{diag}(-1,1,1).
\]

The constructive result is

\[
\boxed{
R=\begin{pmatrix}
1&-a&-\xi\\
0&\sqrt p&a\xi/\sqrt p\\
0&0&1/\sqrt p
\end{pmatrix},\quad
\widehat G=R^T\eta R,\quad \det R=1,
}
\tag{1}
\]

\[
\boxed{
\Omega=R^{-1}\omega R+R^{-1}\mathrm dR,\qquad
\omega^T\eta+\eta\omega=0,\qquad
\mathsf E=R^{-1}\vartheta.
}
\tag{2}
\]

Here \(\omega\) is any smooth Lorentz-algebra-valued one-form and
\(\vartheta\) any invertible coframe. Equation (2) parameterizes **every**
metric-compatible native connection and **every** local soldering for this
resolved state. It gives

\[
g=\mathsf E^T\widehat G\mathsf E=\vartheta^T\eta\vartheta,
\quad
\mathcal R_\Omega=R^{-1}(\mathrm d\omega+\omega\wedge\omega)R,
\quad
R T_{\mathsf E}=\mathrm d\vartheta+\omega\wedge\vartheta.
\tag{3}
\]

Thus regularity and the admissible freedom are explicit. A chosen coframe
and zero torsion determine the physical Levi-Civita connection uniquely.
The native state and Cella presentation identities alone allow the whole
family (2).

The domain is a smooth local patch with finite smooth \(a,\xi\). In three
dimensions \(\mathsf E\) is a rank-three coframe. The existing positive
fourth-direction branch extends every formula by
\(\widehat G_4=\widehat G\oplus1\), \(R_4=R\oplus1\), and
\(\eta_4=\eta\oplus1\), with rank-four coframes. This retains that branch's
declared assumptions. The scalar \(d\) below denotes signed native volume;
\(\mathrm d\) is exterior differentiation.

## 1. What Cella's full presentation law supplies

The source results are the
[Gauge Channel Transport Law](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/dbp_role_channel_and_orbit_geometry/Gauge_Channel_Transport_Law.md>)
and [Canonical Invariant Reduction Theorem](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/dbp_role_channel_and_orbit_geometry/Canonical_Invariant_Reduction_Theorem.md>).
They concern regular relation jets in a fixed Euclidean presentation basis.
At a point, a defining-function gauge \(\widetilde f=\mu f\),
\(f=0\), \(\mu=1\), acts as

\[
g\mapsto g,\qquad H\mapsto H+g u^T+u g^T,
\qquad u=\nabla\log\mu.
\tag{4}
\]

There is a complete tensor quotient behind the invariant sums. For
\(q=g^Tg>0\) and \(P=I-gg^T/q\),

\[
\boxed{
\operatorname{Sym}^2(\mathbb R^n)^*/
\{g u^T+u g^T:u\in\mathbb R^n\}
\;\cong\;\operatorname{Sym}^2(g^\perp)^*,
\qquad [H]\mapsto PHP.
}
\tag{5}
\]

Indeed, the rank-two gauge term vanishes between tangent vectors.
Conversely the unique gauge taking a symmetric \(H\) to its tangent
representative is

\[
\boxed{
u_*=-\frac{Hg}{q}+\frac{g^THg}{2q^2}g,
\qquad H+g u_*^T+u_*g^T=PHP.
}
\tag{6}
\]

If \(g u^T+u g^T=0\), tangent-normal and normal-normal components give
\(u=0\); this proves uniqueness. The tangent quadratic form, with the
source's normal normalization, supplies all its elementary curvature
invariants. Equation (5) retains the whole tangent jet, not just one trace
or determinant.

Equation (6) is a **point-jet representative**. A field of such pointwise
choices need not be one scalar defining-function gauge. A prescribed
one-form \(u\cdot\mathrm dx\) integrates locally to \(\log\mu\) exactly
when it is closed. If \(\mu=1\) on an entire hypersurface, its tangential
derivative must also vanish there. For the source's keystone jet
\(g=(3,1,2)^T\),
\(H=\left(\begin{smallmatrix}2&1&0\\1&0&0\\0&0&2\end{smallmatrix}\right)\),

\[
P u_*=-PHg/q=(-1,-5,4)^T/98\ne0.
\tag{7}
\]

So independently normalizing every point while fixing every gradient is
not a single \(\mu=1\) surface gauge. Transport along a relation family
must respect its jet compatibility, in addition to the pointwise algebra.

### The exact native account embedding

[XLINK-1](/home/williaml/seated-root/xlink1_cella_coupling_form.py)
already realizes \(B=2I-\mathbf1\mathbf1^T\) as the triple-horizon Gram.
Use its trine centered on the native clock:

\[
J=\begin{pmatrix}
1/\sqrt3&1/\sqrt3&1/\sqrt3\\
2/\sqrt3&-1/\sqrt3&-1/\sqrt3\\
0&1&-1
\end{pmatrix},\qquad
J^T\eta J=B,\qquad J\mathbf1=\sqrt3 e_0.
\]

The resolved, state-dependent embedding is \(U=R^{-1}J\). With
\(c=e_0\) denoting the native clock coefficient vector,

\[
\boxed{
U^T\widehat G U=B,\qquad U\mathbf1=\sqrt3 c.
}
\tag{8}
\]

For the three-channel account \(C=(\kappa_c,\kappa_s,\kappa_{int})^T\),
write \(K=\Sigma C\), \(C_\perp=C-(K/3)\mathbf1\). Its complete native
encoding is

\[
\boxed{
V_C=UC=\frac{K}{\sqrt3}c+UC_\perp,\qquad
\widehat G(c,V_C)=-\frac K{\sqrt3},\qquad
\widehat G(UC_\perp,UC_\perp)=2\|C_\perp\|^2.
}
\tag{9}
\]

This identifies the invariant quotient with the clock component and the
whole zero-sum residue with its positive rest plane. The equal-channel
representative is the algebraic section \(K\mapsto(K/3)\mathbf1\);
realizability of a particular representative by an actual defining-function
gauge is governed by (4), not by the quotient notation alone. \(V_C\)
retains the units of the channel values; this encoding supplies no
position, effort, or source law.

The weighted **edge** vector
\(\nu=(g_1H_{23},g_2H_{13},g_3H_{12})^T\) is a different object. Its
encoding \(V_\nu=U\nu\) obeys
\(\widehat G(V_\nu,V_\nu)=\nu^TB\nu\), while its gradient-fixing gauge action (4) is

\[
\boxed{
\nu\mapsto\nu+M(g)u,\qquad
M(g)=\begin{pmatrix}
0&g_1g_3&g_1g_2\\
g_2g_3&0&g_1g_2\\
g_2g_3&g_1g_3&0
\end{pmatrix},\quad
\det M=2(g_1g_2g_3)^2.
}
\tag{10}
\]

When all \(g_i\ne0\), these point-gauge translations reach every edge
vector. They generally change its Lorentz norm: the keystone's opposite
vertex gauge gives \(\nu'=(3t,3t,2)^T\) and
\(\nu'^TB\nu'-\nu^TB\nu=-24t\). The other channels compensate in
\(\Sigma C\). Thus the Cella defining-function gauge is an exact
presentation action on jet/account values; it is not itself a Lorentz
change of native frame.

Here \(K\) is the invariant curvature of the represented Cella relation.
The curvature \(\mathcal R_\Omega\) of a spacetime transport in (3) is a
different typed object. Equations (8)–(13) encode the relation account;
identifying one of its readouts with a physical curvature requires the
corresponding native coupling rule.

### The complete nonzero gauge quotient and canonical edge map

The raw edge ambiguity has a canonical geometric resolution. At a regular
point of \(f=0\), now allow every nonzero real gauge value \(\mu\), and put
\(u=\nabla\mu\). Direct differentiation gives

\[
g'=\mu g,\qquad
H'=\mu H+g u^T+u g^T,\qquad
q'=\mu^2q,\qquad P'=P.
\tag{10a}
\]

In any Euclidean presentation dimension \(n\ge2\), define

\[
\boxed{
E=PHP,\qquad
T_{ijk}=\frac{g_iE_{jk}}q.
}
\tag{10b}
\]

Indices are lowered with the Euclidean presentation metric. Equivalently
\(T(X,Y)=g\,E(X,Y)/q\) is a vector-valued symmetric bilinear form. It
vanishes on normal input slots, takes values in the normal line, and is
exactly \(-B(PX,PY)\), the negative of the Euclidean normal-valued second
fundamental form extended tangentially. This is the standard second
fundamental form of the specified relation jet; it does not use a selected
physical spacetime metric.

Since \(E'=\mu E\), both scale factors cancel:

\[
\boxed{T'=T\quad\text{for every }\mu\ne0.}
\tag{10c}
\]

This includes negative rescaling: the normalized gradient and its scalar
shape operator reverse together, while the normal-valued tensor remains
fixed. The tensor is computed directly from the given relation. It does
not require integrating the pointwise gauge (6) into one defining
function over the whole surface.

**Complete quotient theorem.** The regular Euclidean second jets
\((g,H)\), \(g\ne0\), modulo (10a), are in bijection with

\[
\boxed{
(P,T),\qquad P^T=P^2=P,\quad \operatorname{rank}P=n-1,\qquad
T\in\ker P\otimes\operatorname{Sym}^2(\operatorname{im}P)^*,
}
\tag{10d}
\]

with the last two slots extended by \(P\). To construct a representative,
choose a unit normal \(n\in\ker P\), set \(A_{jk}=n^iT_{ijk}\), and choose
any \(s\ne0\) and any covector \(v_f\). Then
\(g=sn\), \(H=sA+g v_f^T+v_f g^T\) realize \((P,T)\).
Conversely, equal projectors imply \(g'=\mu g\) for a unique \(\mu\ne0\).
Equal tensors then imply \(P H'P=\mu PHP\). Their Hessian difference has
zero tangent restriction and hence is uniquely \(g u^T+u g^T\), by (5).
This proves both directions. Every such value/gradient pair
\((\mu,u)\) is locally the first jet of a smooth nonvanishing multiplier.

Retaining \(P\) is necessary at flat points: \(H=0\) with \(g=e_1\) or
\(g=e_2\) gives the same zero tensor but different tangent planes, and the
jets are not related by a defining-function rescaling. Thus \(T\) alone
is not a complete quotient on the zero-curvature stratum.

There is an exact reconstruction separating geometric and gauge data:

\[
\boxed{
v_f=\frac{Hg}{q}-\frac{g^THg}{2q^2}g=-u_*,
\qquad
H=E+g v_f^T+v_f g^T,\qquad
E_{jk}=\sum_i g_iT_{ijk}.
}
\tag{10e}
\]

Under (10a), \(v_f'=v_f+u/\mu\). This is a gauge remainder covector,
distinct from a native ruler or a physical load/transport law.

For the specified three role axes, extract the canonical edge readout

\[
\boxed{
\bar\nu=(T_{123},T_{213},T_{312})^T,\qquad
\nu=q\bar\nu+M(g)v_f,\qquad
\bar\nu'=\bar\nu.
}
\tag{10f}
\]

Indeed the three off-diagonal entries of (10e) give the decomposition
using exactly the existing \(M(g)\) in (10). Its full gauge transformation
is \(\nu'=\mu^2\nu+\mu M(g)u\), consistent with
\(M(g')=\mu^2M(g)\) and the shift in \(v_f\). For the source's Gaussian
coupling channel, applying its formula to the projected jet gives
\(\kappa_c^{\rm projected}=-\bar\nu^TB\bar\nu\). This projected readout is
well-defined under the full defining-function gauge and generally differs
from the raw coupling channel.

The native map is consequently

\[
\boxed{
\bar V=U\bar\nu,\qquad
\widehat G(\bar V,\bar V)=\bar\nu^TB\bar\nu,\qquad
D_0\bar V=U\,\mathrm d\bar\nu.
}
\tag{10g}
\]

The last identity applies where the relation jet field and the native
state field are jointly specified, or under a specified pullback relating
their bases. Its norm and vector are independent of the defining function
representing that relation. A passive native basis change acts on
\(\bar V\) by \(S^{-1}\), as in (14).

### Full axis covariance belongs to the tensor

For an orthogonal change of Euclidean presentation axes \(O\),

\[
\boxed{
P'=OPO^T,\qquad
T'_{ijk}=O_{ia}O_{jb}O_{kc}T_{abc}.
}
\tag{10h}
\]

This follows from \(g'=Og\), \(H'=OHO^T\), and preserves the complete
normal/tangent geometry. The three extracted edge components do not
carry this full tensor action. A specific counterexample is
\(f=z+(x^2-y^2)/2\) at the origin:

\[
g=e_3,\qquad H=\operatorname{diag}(1,-1,0),\qquad
\bar\nu=0,\qquad \|T\|^2=2.
\]

Rotate the two tangent axes by
\[
O=\begin{pmatrix}
1/\sqrt2&1/\sqrt2&0\\
-1/\sqrt2&1/\sqrt2&0\\
0&0&1
\end{pmatrix}.
\]
The same geometric tensor then has
\(\bar\nu'=(0,0,-1)^T\) and
\(\bar\nu'^TB\bar\nu'=1\), while its full squared norm remains \(2\).
No linear Lorentz transformation sends the initial zero edge vector to
this nonzero one. Thus (10g) gives a defining-function-gauge-invariant
readout in the specified role axes; its \(B\)-norm is not a scalar under
arbitrary rotations of those axes.

The full tensor supplies actual Euclidean invariants, for example
\[
\|T\|^2=\frac{\operatorname{tr}(E^2)}q,\qquad
\left\|\operatorname{tr}_{2,3}T\right\|^2
=\frac{(\operatorname{tr}E)^2}q.
\]
The trace vector is minus the mean-curvature vector. These contractions
and the full tensor retain geometry that edge extraction discards.

### Encoding the full shape operator in the native trine

There is an invertible bridge for the whole shape, with an exact tangent
frame transformation. Let a two-dimensional cut have positive metric \(h\),
and let \(S_n\) be its shape operator for a specified normal \(n\), using
the [implicit reconstruction](/home/williaml/seated-root/docs/2026-09-05-cella-normal-geometry.md).
In an orthonormal tangent basis write

\[
S_n=\begin{pmatrix}p&t\\t&r\end{pmatrix}.
\]

Take the tangent trine
\(e_1=(1,0)\), \(e_2=(-1,\sqrt3)/2\), \(e_3=(-1,-\sqrt3)/2\).
Its three directional normal curvatures are

\[
\boxed{
\chi(S_n)=
\begin{pmatrix}
p\\(p+3r-2\sqrt3t)/4\\(p+3r+2\sqrt3t)/4
\end{pmatrix},\qquad \chi_i=h(S_ne_i,e_i).
}
\tag{10i}
\]

These sample the shape operator; the earlier \(C\) denotes the source's
Gaussian-curvature channel account. Their complete linear inverse is

\[
\boxed{
p=\chi_1,\qquad
r=\frac{2(\chi_2+\chi_3)-\chi_1}{3},\qquad
t=\frac{\chi_3-\chi_2}{\sqrt3}.
}
\tag{10j}
\]

**Shape-trine theorem.** For every pair of self-adjoint shape operators
\(S,T\), this isomorphism obeys

\[
\boxed{
\Sigma\chi(S)=\frac32\operatorname{tr}S,\qquad
\chi(S)^TB\chi(T)
=-\frac32\bigl(\operatorname{tr}S\,\operatorname{tr}T
                    -\operatorname{tr}(ST)\bigr).
}
\tag{10k}
\]

In particular \(\chi(S)^TB\chi(S)=-3\det S\). Equations (10i)–(10j)
prove invertibility; substitution proves (10k), including all mixed terms.
The determinant identity is valid for every symmetric matrix, including
singular, indefinite and zero shapes.

The native image \(V_n=U\chi(S_n)\) therefore satisfies

\[
\boxed{
R V_n=\sqrt3
\begin{pmatrix}
(p+r)/2\\(p-r)/2\\-t
\end{pmatrix},\quad
\widehat G(c,V_n)=-\frac{\sqrt3}{2}\operatorname{tr}S_n,\quad
\widehat G(V_n,V_n)=-3\det S_n.
}
\tag{10l}
\]

The mean curvature occupies the clock component; the two shear components
occupy its rest plane. All three are retained by (10j).

For any orthogonal tangent change \(O\), define
\(\mathcal L_O=\chi\circ(S\mapsto OSO^T)\circ\chi^{-1}\).
Conjugation preserves traces and the polarized determinant, so

\[
\boxed{
\mathcal L_O^TB\mathcal L_O=B,\qquad
\mathcal L_O\mathbf1=\mathbf1,\qquad
\Sigma\mathcal L_O=\Sigma.
}
\tag{10m}
\]

This is a representation: \(\mathcal L_{O_1O_2}
=\mathcal L_{O_1}\mathcal L_{O_2}\). Tangent rotations act with twice
their angle on the two shear components. Reflections obey the same
quadratic-form law. Thus a tangent frame change induces a definite
clock-preserving Lorentz transformation of \(V_n\); its squared norm
remains invariant. A coorientation reversal gives \(S_n\mapsto-S_n\),
\(V_n\mapsto-V_n\), leaving that norm fixed.

This extends to every normal codimension: \(S_n\) and \(V_n\) depend
linearly on \(n\). Given the tangent metric and nondegenerate normal
pairing, the map \(n\mapsto V_n\) recovers the entire second fundamental
form by (10j) and
\(g(B(X,Y),n)=h(S_nX,Y)\). Equation (10k) simultaneously determines
all polarized normal determinant invariants. It retains nonzero
trace-free shapes, such as \(S=\operatorname{diag}(1,-1)\), whose encoded
norm is \(3\).

The encoded null condition is \(\det S_n=0\). A physical coupling must
determine its relation to \(F=0\) or \(\theta_L=0\). The construction is a
complete curvature-to-native-vector map for the stated cut data.

Equations (10b)–(10h) give the complete defining-function gauge quotient
and the canonical edge readout in fixed role axes. Equations (10i)–(10m)
add a full-shape encoding with covariant tangent-axis changes. Selecting
an actual load coefficient, response law, or spacetime attachment remains
physical information to supply to these maps.

## 2. The regular native connection and its full freedom

Use the convention \(D v=\mathrm dv+\Omega v\). A connection preserves
the native pairing exactly when

\[
\mathrm d\widehat G=\Omega^T\widehat G+\widehat G\Omega.
\]

Splitting \(\widehat G\Omega\) into symmetric and skew parts gives the
complete affine family

\[
\boxed{
\Omega=\tfrac12\widehat G^{-1}\mathrm d\widehat G
       +\widehat G^{-1}A,\qquad A^T=-A.
}
\tag{11}
\]

Conjugating by the explicit lift (1) gives the equivalent family (2).
There are three arbitrary one-forms in \(\omega\) in dimension three,
and six in the positive fourth-direction extension. A curve fixes
transport through \(\dot v+\Omega(\dot\gamma)v=0\); specifying the
Gram path alone leaves precisely this skew freedom.

The flat presentation representative \(\omega=0\) is particularly simple:

\[
\boxed{
\Omega_0=R^{-1}\mathrm dR
=\frac1p\begin{pmatrix}
0&-\mathrm da&-\mathrm d\xi\\
0&a\,\mathrm da&\xi\,\mathrm da+a\,\mathrm d\xi\\
0&0&-a\,\mathrm da
\end{pmatrix}.
}
\tag{12}
\]

It obeys
\(\mathrm d\Omega_0+\Omega_0\wedge\Omega_0=0\),
\(\operatorname{tr}\Omega_0=0\), and \(D_0c=0\). It transports the Cella
embedding exactly:

\[
\boxed{D_0U=0,\qquad D_0V_C=U\,\mathrm dC,\qquad
D_0V_\nu=U\,\mathrm d\nu.}
\tag{13}
\]

In particular an invariant-preserving Cella account path has a purely
rest-plane derivative in this representative. These equations construct
the covariant native presentation transport; a physical connection can
add the unrestricted \(R^{-1}\omega R\) in (2).

### Covariance and curvature distinguish the choices

Under a native basis change \(S\), coefficient vectors become
\(v'=S^{-1}v\), and

\[
\begin{aligned}
\widehat G'&=S^T\widehat G S,&
\Omega'&=S^{-1}\Omega S+S^{-1}\mathrm dS,\\
\mathsf E'&=S^{-1}\mathsf E,&
\mathcal R'&=S^{-1}\mathcal R S.
\end{aligned}
\tag{14}
\]

The skew parameter in (11) has the inhomogeneous transformation

\[
\boxed{A'=S^TAS+\tfrac12
(S^T\widehat G\,\mathrm dS-\mathrm dS^T\widehat G S).}
\tag{15}
\]

Consequently setting \(A=0\) independently in every basis is not a
covariant selector. It also does not give the flat presentation connection:
with \(Q=\widehat G^{-1}\mathrm d\widehat G\),

\[
\mathcal R_{Q/2}=-\tfrac14Q\wedge Q,
\qquad
\left.\mathcal R_{Q/2}(\partial_a,\partial_\xi)\right|_{a=\xi=0}
=\frac14\begin{pmatrix}0&0&0\\0&0&1\\0&-1&0\end{pmatrix}\ne0.
\tag{16}
\]

There is a direct obstruction to a selector depending only on the Gram
matrix and its derivatives: a varying Lorentz basis change leaves
\(\widehat G=\eta\) constant but requires the extra
\(S^{-1}\mathrm dS\) in the connection. The ordered seat lift (1)
supplies a definite presentation section. Recomputing an independent
triangular section after an arbitrary basis change may introduce an
additional varying Lorentz rotation; it must be transported as well.

## 3. Crossing the raw pinch

The [directional resolution](/home/williaml/seated-root/pinch_directional_resolution.py)
retains \(\xi=(b-a\gamma)/d\) and replaces the collapsing ruler by
\(v=(g-\gamma h)/d\). On a chosen real branch near \(d=0\),

\[
\gamma=\sigma\sqrt{1-Fd^2},\qquad
b=a\gamma+\xi d,\qquad \sigma=\pm1.
\]

The raw frame is the resolved frame times

\[
S_d=\begin{pmatrix}1&0&0\\0&1&\gamma\\0&0&d\end{pmatrix},
\quad \det S_d=d,\quad G_{\rm raw}=S_d^T\widehat G S_d.
\]

For \(d\ne0\), its connection has the basis-change form (14), so

\[
\boxed{
\operatorname{tr}\Omega_{\rm raw}=\frac{\mathrm dd}{d},
\qquad \operatorname{tr}\Omega=0.
}
\tag{17}
\]

A transverse zero of \(d\) therefore makes the raw frame's volume
connection singular. The resolved transport (12) contains no division by
\(d\) and is smooth for every finite smooth \(a,\xi\), including
\(d=0\) and both \(\xi=\pm1\). This is an explicit regular transport on
the resolved extension; it does not make the rank-deficient raw basis
invertible. The [native load calculation](/home/williaml/seated-root/docs/2026-09-05-pinch-load-transport.md)
already supplies the finite \(\xi\) and its rate at a transverse crossing.
Those rates can now be substituted directly into (12).

## 4. Soldering, torsion, and coordinate integration

Given (2), the spacetime connection induced by the native one has
coordinate matrices

\[
\boxed{
\Gamma_\mu=\mathsf E^{-1}
(\partial_\mu\mathsf E+\Omega_\mu\mathsf E).
}
\tag{18}
\]

Metric compatibility and all three identities (3) follow by substitution.
Conversely every native compatible connection determines the unique
\(\omega=R\Omega R^{-1}-\mathrm dR\,R^{-1}\), and every soldering
determines the unique \(\vartheta=R\mathsf E\). This proves completeness
of the parameterization.

For any chosen \(\vartheta\), zero torsion fixes \(\omega\) uniquely.
For example, if
\(\mathrm d\vartheta^a=-\tfrac12 C^a{}_{bc}\vartheta^b\wedge\vartheta^c\),
set \(C_{abc}=\eta_{ad}C^d{}_{bc}\) and
\(\omega_{ab}=\Gamma_{abc}\vartheta^c\). The explicit solution is

\[
\boxed{
\Gamma_{abc}=\tfrac12(C_{acb}-C_{cba}+C_{bac}).
}
\tag{19}
\]

This is the usual orthonormal-frame Levi-Civita formula. All other
compatible connections are \(\omega_{\rm LC}+\kappa\), with
\(\kappa_{ab}=-\kappa_{ba}\); their torsion is
\(\kappa\wedge\vartheta\). A difference with zero torsion would be
both skew in its first two lowered indices and symmetric in its last two,
forcing it to vanish. Thus there is no further connection ambiguity once
the full coframe and zero torsion are fixed.

An invertible coframe already defines a local soldering. It need not be a
coordinate differential. The stronger requirement of integrating a pair
\((\Omega,\mathsf E)\) into a **fixed flat vector space** has the exact
local conditions

\[
\boxed{\mathcal R_\Omega=0,\qquad T_{\mathsf E}=0.}
\tag{20}
\]

On a sufficiently small simply connected patch, flatness solves
\(\mathrm dH=H\Omega\). Then
\(\mathrm d(H\mathsf E)=H T_{\mathsf E}=0\), so
\(H\mathsf E=\mathrm dX\); invertibility makes \(X\) local coordinates.
The converse follows by differentiating. Equivalently, the affine
connection
\(\left(\begin{smallmatrix}\Omega&\mathsf E\\0&0\end{smallmatrix}\right)\)
has curvature
\(\left(\begin{smallmatrix}\mathcal R_\Omega&T_{\mathsf E}\\0&0\end{smallmatrix}\right)\).
Global integration additionally requires trivial holonomy and vanishing
periods. Curved spacetime coframes are fully allowed; they do not satisfy
the stronger flat integration condition (20).

For every smooth resolved state, an explicit regular flat completion is
\(\mathsf E_0=R^{-1}\mathrm dX\) with \(\Omega_0\). More generally any
local Lorentzian metric, expressed in an orthonormal coframe
\(\vartheta\), gives a completion \(R^{-1}\vartheta\). Native state
algebra therefore imposes no additional local metric restriction within
this class. The coframe also attaches the selected native clock and
rulers to spacetime directions; that attachment can matter even when two
coframes give the same metric. Passive basis changes act simultaneously
on all these data as in (14).

## 5. Integrating the native screen into surfaces

Near \(\xi=\pm1\), put
\(w=\sqrt{(a^2+\xi^2)/p}>0\). The native dual normal has orthonormal
components

\[
Rk=(1,-a/\sqrt p,-\xi/\sqrt p,0)^T,
\qquad k=c+w r_n.
\]

For the four-dimensional completion, the two-dimensional candidate screen
is the common kernel of the independent normal one-forms

\[
\alpha=\vartheta^0,\qquad
\beta=-\frac{a\vartheta^1+\xi\vartheta^2}{\sqrt{a^2+\xi^2}}.
\]

It integrates locally into surface cuts exactly when

\[
\boxed{
\alpha\wedge\beta\wedge\mathrm d\alpha=0,\qquad
\alpha\wedge\beta\wedge\mathrm d\beta=0.
}
\tag{21}
\]

For screen tangent \(X,Y\),
\(\alpha([X,Y])=-\mathrm d\alpha(X,Y)\), and similarly for \(\beta\).
These equations are precisely involutivity; local Frobenius integration
then gives the cuts. This condition is independent of merely preserving
the native Gram or even integrating the ambient coframe.

For an explicit obstruction, take the flat coframe
\((\mathrm dt,\mathrm dx,\mathrm dy,\mathrm dz)\),
\(a=z\), \(\xi=1\), \(d=t\). Then
\(\beta=-(z\,\mathrm dx+\mathrm dy)/\sqrt{1+z^2}\) and

\[
\boxed{
\alpha\wedge\beta\wedge\mathrm d\beta
=\frac{\mathrm dt\wedge\mathrm dx\wedge\mathrm dy\wedge\mathrm dz}
{1+z^2}\ne0.
}
\tag{22}
\]

Indeed \(X=\partial_x-z\partial_y\) and \(Y=\partial_z\) are screen
tangent, but \([X,Y]=\partial_y\) is not. All state and transport
coefficients are regular; the failure is specifically the absence of
surface cuts orthogonal to that prescribed native radial field. For an
actual cut and a torsion-free connection, the associated second
fundamental forms are then obtained by covariant differentiation. Their
curvature and optical propagation use (3); the congruence and source
conditions remain those declared for the physical completion.

## 6. A curved crossing with independently specified area

Use the previously solved native state path

\[
a=0,\qquad d=t,\qquad \xi=1+kt,\qquad k>0,
\]

and declare the smooth coframe, for \(\lambda>0\),

\[
\vartheta=(\mathrm dt,e^{\lambda t}\mathrm dx,
\mathrm dr,e^{\lambda t}\mathrm dz)^T.
\]

The resulting native soldering and metric are explicitly

\[
\boxed{
\mathsf E=(\mathrm dt+\xi\,\mathrm dr,
e^{\lambda t}\mathrm dx,\mathrm dr,e^{\lambda t}\mathrm dz)^T,
\quad
g=-\mathrm dt^2+\mathrm dr^2+e^{2\lambda t}
(\mathrm dx^2+\mathrm dz^2).
}
\tag{23}
\]

Its determinant is \(-e^{4\lambda t}\). Choose a neighborhood with
\(|kt|<1/2\) and \(1-Ft^2>0\), so the native raw branch is also real.
Both metric and resolved transport extend through the transverse signed
volume crossing at \(t=0\), where \(F=0\).

The torsion-free orthonormal connection has only

\[
\omega^0{}_i=\omega^i{}_0=\lambda\vartheta^i,
\qquad i=1,3.
\]

Equation (2) supplies its regular native connection. Independently
computing the complete coordinate Ricci tensor in the order \((t,x,r,z)\)
gives

\[
\boxed{
\operatorname{Ric}=\operatorname{diag}
(-2\lambda^2,2\lambda^2e^{2\lambda t},0,
2\lambda^2e^{2\lambda t}),\qquad
\operatorname{Scal}=6\lambda^2.
}
\tag{24}
\]

The native clock is \(C=\partial_t\), its radial unit field is
\(R_n=-\partial_r\), and \(K=C+\xi R_n\). The cuts \(t,r=\mathrm{constant}\)
have area density \(e^{2\lambda t}\). For the same normalization used in
the [load transport calculation](/home/williaml/seated-root/docs/2026-09-05-pinch-load-transport.md),

\[
L=\frac{1+\xi}{2}(C+R_n),\qquad
N=\frac{C-R_n}{1+\xi},\qquad g(L,N)=-1,
\]

\[
\boxed{
\theta_K=2\lambda,\qquad
\theta_L=\lambda(1+\xi),\qquad
\theta_N=\frac{2\lambda}{1+\xi}.
}
\tag{25}
\]

At the null native crossing these are \(2\lambda,2\lambda,\lambda\).
Thus the same finite-time native pinch path admits this expanding curved
completion as well as the previously specified trapping completion.
Equation (25) follows from its declared coframe; no condition forcing
\(\theta_K\) to vanish with \(F\) was inserted. Setting \(\lambda=0\)
also gives a flat completion of the same state. The difference is concrete
coframe and area response data.

## Closure and retained dependencies

For **EXT-001**, equations (1), (8), (12), and (18) construct the native
representatives, their regular transport and every local soldering. For
**EXT-005**, equations (2), (3), and (19) give the complete compatible
connection, curvature and torsion identities; choosing the coframe and
zero torsion determines the physical connection. For **EXT-007**,
equations (20)–(22) supply the exact local coordinate and cut integration
conditions, with a nonintegrability witness.

These are constructive closure statements within the declared smooth,
finite-direction class. The physical load still has to select its
coframe/attachment and state evolution; if torsion is allowed, its
contorsion is additional data. Cella point-gauge transport does not
silently supply those laws. Equations (23)–(25) exhibit the resulting
nonuniqueness quantitatively.

The \(\nu\)-map's defining-function ambiguity is now removed by the
complete quotient \((P,T)\), canonical edges \(\bar\nu\), and native image
\(U\bar\nu\) in (10b)–(10h), followed by the invertible full-shape trine
encoding and its tangent-frame action in (10i)–(10m). The remaining role-axis/native attachment is
specified data, rather than an unresolved defining-function gauge.

Standard exterior calculus, metric-compatible connections, Cartan
structure equations, the Levi-Civita construction, local flat-connection
integration, the Poincare lemma, and Frobenius are **retained mathematical
tools** with the domains stated above. Their native applications are
proved explicitly here; reproducing their general foundations is not a
remaining physical selection problem. No online source was consulted for
this extension. The Cella source files and native constructions above were
read directly, with provenance retained in the
[dependency ledger](/home/williaml/seated-root/docs/EXTERNAL-MATHEMATICS-DEBT.md).

The [exact runner](/home/williaml/seated-root/cella_frame_transport.py) records
its exact symbolic checks and their count in the
[verification receipt](/home/williaml/seated-root/docs/cella-frame-transport-checks.json).
They cover the full metric and connection identities, Cella tensor
quotient, full nonzero defining-function gauge and canonical edge map,
orthogonal tensor covariance, full-shape trine inverse and mixed determinant
form, account embedding, nonconstant basis changes, nonzero
curvature, raw-volume singularity, torsion, the independent full Ricci
calculation, and the surface integrability obstruction.

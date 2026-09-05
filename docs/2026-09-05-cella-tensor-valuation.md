# Cella tensor curvature and weighted finite-jet calculus

5 September 2026.

**Result.** Cella's diagonal scalar-curvature calculus extends to arbitrary coordinate metric matrices, every fixed signature, the entire Riemann tensor, and any finite contraction of curvature and its covariant derivatives. The extension uses polynomial curvature numerators and the metric determinant. For monomial/smooth metric germs these numerators have finite base support; a proved output-dependent Taylor depth computes every requested weighted coefficient. Leading terms are certified after numerator and determinant cancellations have both been resolved.

This is a mathematical extension conditional on a supplied metric and its Levi-Civita connection. It computes curvature from that metric. Native load dynamics, the spacetime coframe, and the physical selection of this connection remain model data.

The [exact replay](</home/williaml/seated-root/cella_tensor_valuation.py>) and its [JSON record](</home/williaml/seated-root/docs/cella-tensor-valuation-checks.json>) accompany the proofs.

## 1. Sources, retained tools, and the exact extension

The live Cella DAG identifies `DBP:thm:local_curvature_complete`, `DBP:thm:weighted_jet`, and `DBP:gap:II4`. The last explicitly requests a non-diagonal Laurent-jet ring, curvature cross terms and tensor pullback replay. The backing mathematics read for this extension is:

- [Local Curvature Calculus COMPLETE](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/local_curvature_and_black_hole_metrics/LOCAL_CURVATURE_CALCULUS_COMPLETE_v1.0.md>), especially the exact directional decomposition, finite-support theorem, master quadric, weighted fronts and cancellation hierarchy.
- [LEAD7 variable-transverse weighted-jet theorem](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/local_curvature_and_black_hole_metrics/LEAD7_VARIABLE_TRANSVERSE_WEIGHTED_JET_THEOREM_v1.0.md>) and its [replay](</home/williaml/Cella Framework/Papers_Library/07_certificates_data_and_reproducibility/local_curvature_and_black_hole_metrics/verify_lead7_variable_transverse_weighted_jet.py>).
- The current [signed diagonal extension](</home/williaml/seated-root/docs/2026-09-05-cella-dependency-reconciliation.md>), [native horizon metric](</home/williaml/seated-root/docs/2026-09-05-horizon-crossing-metric.md>), [pinch valuation](</home/williaml/seated-root/pinch3_cella_valuation.py>) and [coupling-form calculation](</home/williaml/seated-root/xlink1_cella_coupling_form.py>).

Standard coordinate differentiation, the Levi-Civita and curvature definitions, adjugate inversion, tensor pullback and Taylor's theorem are retained explicitly. There was no online search or newly adopted external source. SymPy supplies exact algebra in the replay; the propositions below are proved independently of its success.

The result is broader than extending a sign in the scalar formula: off-diagonal components participate in the determinant and in mixed first-derivative products. The old diagonal residues follow upon specialization, while genuinely new off-diagonal data can change them. A coordinate pullback of the same metric preserves its geometric curvature.

## 2. Universal polynomial curvature numerators

Let \(g=(g_{ij})\) be a symmetric \(C^2\) coordinate metric on an open set, of any constant nondegenerate signature. Put

\[
\Delta=\det g\ne0,\qquad A^{ij}=(\operatorname{adj}g)^{ij},
\qquad g^{ij}=A^{ij}/\Delta.
\tag{1}
\]

Repeated indices are summed. Use

\[
R^i{}_{jkl}
=\partial_k\Gamma^i_{lj}-\partial_l\Gamma^i_{kj}
+\Gamma^i_{ka}\Gamma^a_{lj}-\Gamma^i_{la}\Gamma^a_{kj},
\qquad R_{ijkl}=g_{ia}R^a{}_{jkl}.
\tag{2}
\]

This is Cella's convention: the round unit two-sphere has scalar curvature \(+2\). Define the first-kind Christoffel polynomial and the second-derivative polynomial

\[
C_{abc}=\tfrac12(g_{ab,c}+g_{ac,b}-g_{bc,a}),
\]

\[
T_{ijkl}=\tfrac12
(g_{il,jk}+g_{jk,il}-g_{ik,jl}-g_{jl,ik}).
\tag{3}
\]

**Theorem 1 — Determinant numerator identity.** Set

\[
\boxed{N_{ijkl}
=\Delta T_{ijkl}
+A^{ab}(C_{ail}C_{bjk}-C_{aik}C_{bjl}).}
\tag{4}
\]

Then

\[
\boxed{
R_{ijkl}=\frac{N_{ijkl}}\Delta,\quad
\operatorname{Ric}_{jl}=\frac{M_{jl}}{\Delta^2},\quad
R=\frac P{\Delta^3},}
\tag{5}
\]

where

\[
M_{jl}=A^{ik}N_{ijkl},\qquad P=A^{jl}M_{jl}.
\tag{6}
\]

Every numerator is a polynomial in metric entries and their first and second coordinate derivatives, with rational numerical coefficients. Their metric-factor degrees are respectively \(n+1,2n,3n-1\). Every term has total differential order two.

**Proof.** Substitute \(\Gamma^a_{bc}=g^{ad}C_{dbc}\) into (2), differentiate \(g^{-1}\) by \(\partial g^{-1}=-g^{-1}(\partial g)g^{-1}\), and lower the first index. The second-derivative terms are (3); the first-derivative terms combine to

\[
R_{ijkl}=T_{ijkl}
+g^{ab}(C_{ail}C_{bjk}-C_{aik}C_{bjl}).
\]

Equation (1) gives (4). Contracting the first and third slots gives Ricci, then contracting its two slots gives the scalar, yielding (5)–(6). A determinant has degree \(n\), an adjugate entry degree \(n-1\), and each differentiated metric counts as one metric factor. The stated degrees follow. No step uses positivity. \(\square\)

The determinant powers are universal bounds, not assertions that the resulting poles are sharp. Common polynomial factors and contracted tensor symmetries often cancel them. The algorithm retains these cancellations exactly.

**Corollary 1.1 — All algebraic curvature contractions.** A product of \(k\) lower Riemann tensors with \(p\) inverse-metric contractions or index raisings has the form

\[
\mathcal T=\frac{P_{k,p}}{\Delta^{k+p}},
\qquad
\deg_g P_{k,p}=k(n+1)+p(n-1).
\tag{7}
\]

A complete scalar contraction has \(p=2k\), hence

\[
\boxed{I_k=P_k/\Delta^{3k}.}
\tag{8}
\]

**Proof.** Replace each curvature factor by \(N/\Delta\) and each inverse metric by \(A/\Delta\), and carry out the finite index sums. This also supplies the numerator algorithm for every pairing of the tensor slots. \(\square\)

Polynomial combinations of these contractions follow termwise. In particular this covers scalar curvature, Ricci squared, Riemann squared and dimensionally defined Weyl contractions. A ratio of invariants needs a separate nonzero-denominator condition. An orientation/volume tensor, if requested, introduces its declared \(\sqrt{|\Delta|}\) factor; its valuation is obtained by the same nonzero-leading-coefficient rule. It is not silently included in (8).

## 3. Finite covariant derivatives are covered too

Let \(S_{i_1\ldots i_m}=P_{i_1\ldots i_m}/\Delta^q\) be any covariant tensor whose numerator is already known. Define

\[
\begin{aligned}
(\mathscr D_dP)_{i_1\ldots i_m}
={}&\Delta\partial_dP_{i_1\ldots i_m}
-q(\partial_d\Delta)P_{i_1\ldots i_m}\\
&-\sum_{\ell=1}^m A^{ab}C_{bdi_\ell}
P_{i_1\ldots i_{\ell-1}a i_{\ell+1}\ldots i_m}.
\end{aligned}
\tag{9}
\]

**Theorem 2 — Differential numerator recurrence.** For a sufficiently differentiable metric,

\[
\boxed{\nabla_dS=\frac{\mathscr D_dP}{\Delta^{q+1}}.}
\tag{10}
\]

Consequently \(\nabla^rR_{ijkl}\) has denominator \(\Delta^{r+1}\), numerator metric degree \((r+1)n+1\), and total differential order \(r+2\). Earlier derivative indices are included among the covariant slots when applying the next derivative.

A product of \(k\) curvature factors with a total of \(r\) covariant derivatives and \(p\) inverse metrics has

\[
\boxed{
\mathcal T_{k,r,p}=P_{k,r,p}/\Delta^{k+r+p},\qquad
\deg_gP_{k,r,p}=k(n+1)+rn+p(n-1).}
\tag{11}
\]

For a complete scalar contraction \(p=(4k+r)/2\), so the determinant bound is \(\Delta^{3k+3r/2}\); this scalar pairing requires \(r\) even.

**Proof.** Differentiate the quotient and subtract one Christoffel action for each lower tensor slot. Substitution of \(\Gamma^a_{di}=A^{ab}C_{bdi}/\Delta\) gives (9). Every recurrence increases the metric-factor degree by \(n\) and total differential order by one. Induction followed by the contraction count proves (11). \(\square\)

Thus every specified finite differential curvature invariant has a finite numerator construction. This theorem does not assert a single finite jet determines derivatives of arbitrarily high order.

## 4. Finite support for non-diagonal divisor germs

Use \(t=(z_1,\ldots,z_c,y_{c+1},\ldots,y_n)\), with \(z_a>0\), and write

\[
\boxed{g_{ij}(z,y)=\sum_{\alpha\in S_{ij}}z^\alpha h_{ij,\alpha}(z,y),}
\tag{12}
\]

where every \(S_{ij}\subset\mathbb R^c\) is finite and every coefficient extends smoothly to the corner. The coefficients need not individually be nonzero. Zero entries are allowed. Symmetry and interior nondegeneracy are the metric requirements. This includes arbitrary smooth coordinate metrics by taking \(S_{ij}=\{0\}\), as well as finite Laurent/real-monomial singularities. More rapidly growing arbitrary singular germs require a different asymptotic class.

Let \(e_i\) denote the \(i\)-th normal exponent vector for \(i\le c\), and zero for a tangential index. Finite support is propagated by the exact rules

\[
S(f+g)\subset S(f)\cup S(g),\quad
S(fg)\subset S(f)+S(g),\quad
S(\partial_i f)\subset S(f)-e_i.
\tag{13}
\]

Here \(+\) on supports is Minkowski addition, and the coefficients remain smooth. For a normal derivative the last rule is the exact identity

\[
\partial_i(z^\alpha h)
=z^{\alpha-e_i}(\alpha_i h+z_i\partial_i h).
\tag{14}
\]

If a displayed monomial coefficient vanishes, its zero contribution is retained until all equal exponents have been combined; (13) describes candidate support, not nonzero support.

**Theorem 3 — Finite numerator support.** Apply (13) to the determinant, adjugate, (3)–(4), and any finite recurrence/contraction above. Every numerator admits an exact finite representation

\[
P_{k,r,p}=\sum_{\beta\in B_{k,r,p}}z^\beta f_\beta(z,y).
\tag{15}
\]

The set \(B_{k,r,p}\) is finite and computable from the entry supports, index contractions and derivative indices alone. The coefficient germs depend on metric units and their derivatives through order \(r+2\). If

\[
\Delta=z^\delta d(z,y),\qquad d(0,y_0)\ne0,
\tag{16}
\]

then the tensor or invariant itself has the exact finite-base representation

\[
\boxed{\mathcal T_{k,r,p}
=\sum_{\beta\in B_{k,r,p}}
z^{\beta-(k+r+p)\delta}\,
\frac{f_\beta(z,y)}{d(z,y)^{k+r+p}}.}
\tag{17}
\]

**Proof.** Determinants and adjugates are finite sums of finite products. Equations (13)–(14) preserve the claimed class under every operation in (3)–(11). Equation (16) leaves a smooth nonzero denominator unit, so division introduces no additional base exponents. \(\square\)

For implementation, one can build the support sets before evaluating coefficient functions. For example, if \(B_\Delta,B_A^{ab},B_C^{abc},B_T^{ijkl}\) denote the supports from (3),

\[
\begin{aligned}
B_N^{ijkl}\subset {}&(B_\Delta+B_T^{ijkl})\\
&\cup\bigcup_{a,b}(B_A^{ab}+B_C^{ail}+B_C^{bjk})\\
&\cup\bigcup_{a,b}(B_A^{ab}+B_C^{aik}+B_C^{bjl}).
\end{aligned}
\tag{18}
\]

This replaces the diagonal theorem's at-most-\(n\) directional bases with a finite tensor support that includes every cross term.

**Useful full-matrix specialization.** Let

\[
g=DHD,\qquad D_{ii}=z^{d_i},\qquad
H=H^T,\quad\det H(0,y_0)\ne0.
\tag{19}
\]

Then

\[
g^{-1}=D^{-1}H^{-1}D^{-1},\qquad
\Delta=z^{2\sum_i d_i}\det H.
\tag{20}
\]

Equations (17)–(18) therefore apply with explicit determinant exponent \(\delta=2\sum_i d_i\). The signature is that of \(H\). In the lower Riemann formula, the linear base terms are

\[
d_i+d_l-e_j-e_k,\quad d_j+d_k-e_i-e_l,
\quad d_i+d_k-e_j-e_l,\quad d_j+d_l-e_i-e_k.
\tag{21}
\]

For a quadratic term add any one of

\[
\{d_a+d_i-e_l,\ d_a+d_l-e_i,\ d_i+d_l-e_a\}
\]

to any one of

\[
\{d_b+d_j-e_k,\ d_b+d_k-e_j,\ d_j+d_k-e_b\},
\]

then subtract \(d_a+d_b\), and include the other term in (4). All associated coefficients are smooth because \(H^{-1}\) is smooth.

The factors \(D\) are differentiated in the coordinate metric. If instead \(z^{d_i}dt^i\) is treated as a moving coframe, its exterior derivatives generally do not vanish; those frame commutators must enter a frame curvature calculation. Computing only the curvature of the matrix \(H\) in that coframe would discard real terms. A constant Gram matrix in a moving native frame likewise does not by itself imply zero curvature.

## 5. Weighted cancellation and a certified finite-jet algorithm

Fix a positive weight \(w\in\mathbb R_{>0}^c\), a positive approach ratio \(\rho\), and a transverse base point \(y_0\). Use the actual path

\[
z_a=\rho_a\epsilon^{w_a},\qquad y=y_0,\qquad\epsilon\downarrow0.
\tag{22}
\]

For a finite-base germ \(Q=\sum_\beta z^\beta f_\beta\), define grouped coefficients

\[
q_\lambda(\rho,y_0)
=\sum_{\substack{\beta\in B,\ \alpha\in\mathbb N^c\\
\langle w,\beta+\alpha\rangle=\lambda}}
\frac{\partial_z^\alpha f_\beta(0,y_0)}{\alpha!}\rho^{\beta+\alpha}.
\tag{23}
\]

In (23), \(\beta\) ranges over its finite base set, which may contain real or negative exponents; only \(\alpha\) is a nonnegative integer multi-index. There are finitely many contributing pairs below every finite weight because all \(w_a\) are positive.

**Theorem 4 — Quotient valuation after all cancellations.** Let \(\mathcal T=P/\Delta^q\) be any component or invariant from the preceding theorems. Suppose the grouped numerator and determinant expansions have first finite nonzero coefficients

\[
P(\rho\epsilon^w,y_0)=p_\nu\epsilon^\nu+o(\epsilon^\nu),
\qquad
\Delta(\rho\epsilon^w,y_0)=d_\delta\epsilon^\delta+o(\epsilon^\delta),
\]

with \(p_\nu d_\delta\ne0\). Then

\[
\boxed{
\operatorname{ord}_{w,\rho}\mathcal T=\nu-q\delta,
\qquad
\operatorname{lead}_{w,\rho}\mathcal T=p_\nu/d_\delta^q.}
\tag{24}
\]

The pole order is the negative of this valuation. Neither determinant monomialization nor a nonsingular entrywise initial matrix is required.

**Proof.** Expand each smooth coefficient to enough finite Taylor order for the first nonzero layer, substitute (22), and combine every equal weight. Divide the two displayed leading expansions. A nonzero \(d_\delta\) ensures the metric remains nondegenerate sufficiently far down the selected path. \(\square\)

The exhaustive cancellation mechanisms are Cella's structural, exponent-collision, ratio-front and higher-jet cancellations, now applied separately to numerator and determinant. Inverting an entrywise initial matrix before checking its determinant can miss the determinant's first surviving layer.

**Theorem 5 — An explicit safe jet depth.** Consider a numerator polynomial with metric-factor degree \(d\) and total differential order \(s\). Set

\[
\mu=\min_{i,j,\alpha\in S_{ij}}\langle w,\alpha\rangle,
\quad w_{\min}=\min_a w_a,\quad w_{\max}=\max_a w_a,
\quad b=d\mu-sw_{\max}.
\tag{25}
\]

Every numerator exponent has weight at least \(b\). To compute all coefficients through a requested finite weight \(\Lambda\), it suffices to replace each unit by its normal Taylor polynomial through degree

\[
\boxed{N=s+\max\left(0,\left\lfloor
\frac{\Lambda-b}{w_{\min}}\right\rfloor\right).}
\tag{26}
\]

Keep transverse coefficient derivatives through order \(s\). A full joint \((N+s)\)-jet is a sufficient input if those functions are represented by point jets. Equation (26) is conservative; the explicit support sets in (18) can reduce the depth.

**Proof.** A monomial term of the numerator uses \(d\) component monomials and \(s\) derivatives, so its weight is bounded below by (25). Let \(h_N\) be a Taylor polynomial in the normal variables. For all normal derivatives of order at most \(s\), Taylor remainder estimates give

\[
\partial_z^\alpha\partial_y^\gamma(h-h_N)
=O(|z|^{N+1-|\alpha|})
\]

for the finitely many transverse derivatives used. Expand the difference of the numerator built from \(h\) and from \(h_N\) as products containing at least one such remainder. Along (22), each difference term is

\[
O\!\left(\epsilon^{b+(N+1-s)w_{\min}}\right).
\]

The additional derivative losses in this bound are deliberately conservative. Equation (26) makes the exponent strictly greater than \(\Lambda\). Therefore every coefficient through \(\Lambda\) is unchanged. \(\square\)

For \(P_{k,r,p}\), use

\[
d=k(n+1)+rn+p(n-1),\qquad s=2k+r;
\]

for the determinant use \(d=n,s=0\). This specifies the input depth before any numerator is assembled. It replaces an unexplained fixed Laurent window such as `LO=-8, HI=8` with a proved cutoff derived from the actual question.

Once \(\nu,\delta\) are known, a quotient expansion through weight \(L\ge\nu-q\delta\) requires numerator layers through \(\nu+M\) and determinant layers through \(\delta+M\), where

\[
M=L-(\nu-q\delta).
\tag{27}
\]

After extracting the leading monomial, expand \((d_\delta+u)^{-q}\) by its ordinary binomial series. All positive exponent gaps are locally finite and have a smallest positive gap among the needed terms. Hence only finitely many products contribute through \(M\). This gives every requested inward coefficient, including after a front cancellation.

The executable interface supplies `numerator_jet` and `leading_through`. The latter separately searches numerator and determinant through the caller's stated bounds; it reports an unresolved later layer when either finite search is empty. Increasing the bounds terminates whenever Theorem 4's two finite-first-nonzero hypotheses hold. It does not mistake an empty finite jet for an identically zero smooth germ.

**No uniform cancellation depth exists.** For any integer \(m>2\), the off-diagonal metric

\[
g=dx^2+(1+x^m)^2(dy+dx)^2
\]

has

\[
R=-\frac{2m(m-1)x^{m-2}}{1+x^m}.
\tag{28}
\]

Its first nonzero curvature order can be delayed arbitrarily while all lower metric jets agree with a flat metric. The replay uses \(m=9\): a cutoff at weight four correctly returns no coefficient yet, while weight seven gives \(-144\).

For analytic coefficients the expansions converge locally; formal coefficients give a formal conclusion. For smooth coefficients, the theorem stops at a first finite nonzero layer. If all numerator layers vanish and the determinant has a finite leading layer, the resulting component is smaller than every power along the path, but need not be identically zero. If the determinant is also flat, no finite-jet quotient valuation follows. These are mathematical limits of finite jets, not missing implementation cases.

## 6. Coordinate covariance and what is actually preserved

**Theorem 6 — Full tensor pullback.** For a smooth diffeomorphism \(t=\phi(t')\) and \(J^i{}_a=\partial t^i/\partial t'^a\), set \(g'=J^T(g\circ\phi)J\). Then

\[
R[g']_{abcd}
=J^i{}_aJ^j{}_bJ^k{}_cJ^l{}_d(R[g]_{ijkl}\circ\phi),
\qquad I[g']=I[g]\circ\phi
\tag{29}
\]

for each complete contraction, including finite covariant derivatives. Thus scalar values and valuations along the transported geometric path agree exactly.

This retains the standard tensor naturality of the Levi-Civita curvature. In coordinates, substitution of \(g'\) into the connection formula produces the usual \(J^{-1}\partial^2\phi\) term; its nontensorial contributions cancel in (2). The numerator identities respect the same transformation. In particular,

\[
\Delta'=(\det J)^2(\Delta\circ\phi),\qquad
N'=(\det J)^2J^{\otimes4}(N\circ\phi),
\]

and a scalar numerator with denominator \(\Delta^q\) obeys

\[
P'=(\det J)^{2q}(P\circ\phi).
\tag{30}
\]

The determinant-density factors cancel in the quotient. Equation (30) also explains why a singular interior reparametrization can change numerator and determinant orders separately without changing the scalar along the mapped curve.

For a smooth nonsingular chart change at the limit point, the minimum valuation among components of a covariant tensor is invariant along the mapped path whenever a finite leading tensor exists: the leading tensor is multiplied by an invertible tensor power of \(J(0)\). The valuation of an individual component can change, and a scalar contraction can vanish while the tensor does not.

For strict boundary monomial changes

\[
z_a=c_a(y')\prod_b(z'_b)^{M_{ab}},\qquad y=\psi(y'),
\tag{31}
\]

the path rule is exact: \(w=Mw'\), \(\rho_a=c_a(y'_0)\prod_b(\rho'_b)^{M_{ab}}\). This includes divisor permutations and nonzero boundary rescalings, and interior monomial charts where applicable.

A general boundary-adapted change \(z_a=u_a(z',y')z'_{\sigma(a)}\) maps a straight monomial path to a path with inward corrections. Its generic first nonzero front transforms by the boundary units. After a special front cancellation, the inward corrections matter. Tangential mixing can matter even sooner. One must transport the actual path instead of setting the new transverse coordinates constant and silently identifying a different curve with the old one.

For example the scalar germ \(f=y+x^2\) has order two on \(y=0,x=\epsilon\). Under \(x=x',y=y'+x'\), it has order one on \(y'=0,x'=\epsilon\); those are different physical curves. Divisor orders as orders of divisibility and valuations along the same curve remain well defined. Literal higher Newton support at a frozen transverse slice is not invariant under arbitrary slice-changing maps.

These statements close the tensor-pullback calculation requested in `DBP:gap:II4` for transformations that actually act as metric diffeomorphisms. A Cella role operation must still be shown to induce the specified pullback before that operation can be called a physical gauge equivalence.

## 7. Recovering Cella and testing genuine cross terms

For a diagonal fixed-signature metric \(g_{ii}=\varepsilon_i e^{2u_i}\), substitute its zero off-diagonal entries into (3)–(6). Grouping by differentiation direction gives exactly

\[
R=-2\sum_j\frac{\mathcal Q_j(u)}{g_{jj}},
\]

\[
\mathcal Q_j
=\sum_{i\ne j}\left[u_{i,jj}+u_{i,j}^2-u_{j,j}u_{i,j}\right]
+\sum_{\substack{r<s\\r,s\ne j}}u_{r,j}u_{s,j}.
\tag{32}
\]

Consequently the original bases \(-P_a-2e_a,-P_\mu\), signed master quadrics, generic cubic pole and inverse-channel quartic pole are recovered. With arbitrary transverse coefficient functions, the two LEAD7 germs give

\[
R=\frac{P_1/P_0+Q_1/Q_0}{A_2}x^{-3}+O(x^{-2}),
\qquad
R=-\frac{14}{B}x^{-4}+O(x^{-3}),
\tag{33}
\]

under their stated nonzero leading-metric-coefficient hypotheses. The first displayed term is a pole of exactly that order only when its coefficient is nonzero. The exact replay recovers (33) while retaining arbitrary symbolic transverse functions and their derivatives.

For an actual coordinate pullback, take

\[
g=3x^2dx^2+2x^{-2}dy^2-5x^{-2}dz^2,
\qquad
(x,y,z)=(u(1+v),v+u,w+uv).
\tag{34}
\]

The Jacobian determinant is \(1+v-u\), nonzero near \(u=0,v>0\). The pulled-back metric has nonzero off-diagonal entries. All 81 lower curvature components agree with (29), and

\[
R[g']=-\frac{14}{3u^4(1+v)^4}.
\tag{35}
\]

At \(v=2\) its residue is \(-14/243\). The order-four law is unchanged along the corresponding divisor approach.

For independent cross-term data, however,

\[
g=3x^2dx^2+2c\,dx\,dy+2x^{-2}dy^2
\quad\Longrightarrow\quad
R=-\frac{12}{6-c^2}x^{-4},
\tag{36}
\]

where \(c^2\ne6\). Both the positive-definite and indefinite sectors are covered. Similarly

\[
g=3x^2dx^2+2cx\,dx\,dy+(2+5x)dy^2
\quad\Longrightarrow\quad
R=\frac5{6-c^2}x^{-3}+O(x^{-2}).
\tag{37}
\]

Thus diagonal entries alone do not determine the coefficients once independent off-diagonal entries are admitted.

Two exact cancellation controls make the weighted rule concrete. First,

\[
g=x^{-2}y^{-2}ds^2+x^2dx^2-y^2dy^2,
\qquad R=-6x^{-4}+6y^{-4}.
\tag{38}
\]

On \(x=y=\epsilon\), the scalar vanishes identically along the path. On \(x=\epsilon,y=\epsilon(1+\epsilon)\), the same leading face cancels but the next term is \(-24\epsilon^{-3}\). This is why a higher-order path correction cannot be discarded after a balanced-front cancellation.

Second, set

\[
g=\begin{pmatrix}1&1\\1&1+x-y+x^2\end{pmatrix},
\qquad
\Delta=x-y+x^2,\qquad
R=\frac{1+4y}{2(x-y+x^2)^2}.
\tag{39}
\]

On \(x=y=\epsilon\), the determinant's nominal weight-one face cancels, \(\Delta=\epsilon^2\), and

\[
R=\tfrac12\epsilon^{-4}+2\epsilon^{-3}.
\]

The theorem obtains this by finding the surviving determinant and numerator layers before dividing.

## 8. A full tensor can survive all scalar contractions

Consider the Lorentzian metric

\[
g=2du\,dv+dX^2+dY^2+(X^2-Y^2)du^2.
\tag{40}
\]

Direct use of (4) gives the nonzero independent curvature components

\[
R_{uXuX}=-1,\qquad R_{uYuY}=+1,
\tag{41}
\]

while \(\operatorname{Ric}=0\), \(R=0\), and \(R_{abcd}R^{abcd}=0\).

Indeed every algebraic complete contraction of its curvature tensors vanishes: each nonzero lower curvature factor contains two \(u\) indices and no \(v\) index, but the only nonzero inverse-metric pairing of a lower \(u\) is with \(v\). Some contraction factor must therefore be zero. The curvature tensor itself remains nonzero. This proves why a Lorentzian extension needs tensor valuation, rather than an inference of flatness from scalar-invariant cancellations.

The derivative extension is also checked on the non-diagonal metric \(dx^2+(1+x^3)^2(dy+dx)^2\). The replay verifies the two-dimensional metric-compatibility identity \(\nabla_mR_{ijkl}=\tfrac12(\partial_mR)(g_{ik}g_{jl}-g_{il}g_{jk})\) for its independent component, a second covariant derivative whose value requires the connection on the earlier derivative slot, the certified finite jet of the first derivative, and

\[
g^{ij}(\partial_iR)(\partial_jR)
=\frac{144(1-2x^3)^2}{(1+x^3)^4}.
\tag{42}
\]

The same expression is obtained by the general slot-pairing contraction of two covariant curvature derivatives, with its determinant denominator bound \(\Delta^9\).

## 9. Application to the native horizon-crossing metric

The native completion already gives

\[
g=-A^2F\,dv^2+2AB\,dv\,dr
+h_{AB}(dx^A-U^A dv)(dx^B-U^Bdv),
\qquad \det g=-A^2B^2\det h.
\tag{43}
\]

For smooth nonzero \(AB\) and positive \(h\), its determinant is a unit through \(F=0\). Theorem 1 therefore makes every curvature component finite in a regular coordinate chart there. Theorem 2 gives the same statement for each finite covariant derivative if the coefficients have the required additional differentiability. A vanishing lapse alone supplies no curvature pole.

For the genuinely off-diagonal spherical specialization

\[
g=-F(v,r)dv^2+2dv\,dr+r^2(d\theta^2+\sin^2\theta\,d\phi^2),
\tag{44}
\]

the exact full-tensor contraction yields

\[
\boxed{R=-F_{rr}-\frac4rF_r+\frac{2(1-F)}{r^2},}
\tag{45}
\]

\[
\boxed{R_{abcd}R^{abcd}
=F_{rr}^2+\frac{4F_r^2}{r^2}+\frac{4(1-F)^2}{r^4}.}
\tag{46}
\]

The full tensor retains time-dependent information that these two scalar contractions omit:

\[
R_{vrvr}=\tfrac12F_{rr},\qquad
R_{v\theta v\theta}=\tfrac r2(FF_r-F_v).
\tag{47}
\]

All these formulas permit either sign of \(F\) and use no division by \(F\). The coordinate patch has \(r>0\) and excludes angular-coordinate poles.

For the existing native path

\[
a=-1,\quad b=0,\quad\gamma=1-\rho,\quad
F=\frac{\rho(2-\rho)}{1+2\rho-\rho^2},\quad
r=r_0+L\rho,
\tag{48}
\]

with \(r_0,L>0\), the first horizon is \(\rho=0\). Here \(F_\rho=2\) and \(F_{\rho\rho}=-10\), giving

\[
\boxed{R_H=\frac{10}{L^2}-\frac8{Lr_0}+\frac2{r_0^2},}
\]

\[
\boxed{(R_{abcd}R^{abcd})_H
=\frac{100}{L^4}+\frac{16}{L^2r_0^2}+\frac4{r_0^4}.}
\tag{49}
\]

The fixed-signature tensor calculus therefore verifies an explicit smooth native horizon crossing with nonzero curvature. The calculation uses the declared unit clock/radial specialization and supplied state path; it does not select those functions dynamically.

## 10. Reproducibility and the native boundary

Run:

```bash
python -B /home/williaml/seated-root/cella_tensor_valuation.py
```

The replay independently constructs the differentiated-inverse Christoffel curvature at 24 rational two-jets, covering every signature in dimensions two, three and four. It compares 3,304 lower tensor components and checks their algebraic Bianchi identities. Further exact checks cover the sphere convention, signed Cella laws, arbitrary transverse functions, all 81 components of a nontrivial three-dimensional pullback, independent cross-term residues, numerator and determinant front cancellations, late finite jets, covariant derivatives and differential contractions, null curvature with vanishing scalar invariants, and the native horizon formulas above. The JSON records the current count and elapsed time.

The executable theorem interface supports arbitrary symbolic coordinate metrics and finite monomial-unit sums. Its weighted replay uses exact rational exponents and weights; the proofs allow real exponents and all positive real weights. Symbolic nonzero leading coefficients describe their stated nonvanishing parameter stratum. A coefficient's special zero set must be substituted and reclassified.

This completes the non-diagonal curvature, full tensor, finite-jet and coordinate-pullback computation. For EXT-005 it supplies all curvature terms required once the transport connection is specified. For the finite-jet/asymptotic part of EXT-007 it supplies a safe truncation and cancellation theorem. It does not settle unrelated integrability, ODE, or smooth-division claims by association.

The native theory must still provide the coframe that turns native Gram data into the spacetime metric, the evolution of its coefficients, and the reason its physical transport is the declared Levi-Civita connection. Once those are given, no diagonalization, exterior-metric guess or new curvature derivation is required: their complete tensor jets enter this calculator directly.

# Response on branched, resolved state spaces

5 September 2026 — general theory and exact constructions

**A response belongs to a retained state and a dual pair of effort/flow spaces.** Branch identification, resolution, port transformation and elimination of internal variables are distinct operations on that data. This theory gives their compatibility conditions, the complete finite-dimensional linear passive relation, smooth lifting criteria, and the divisor laws for both two-sided interfaces and one-sided boundaries. It then characterizes predictive reduction and finite-time passage through the strata.

No spacetime signature, metric, horizon interpretation or preferred group is required. Passivity is an optional sector of the theory. Algebraic statements concern the supplied finite cover; differential statements concern the supplied smooth resolution. Existence of a resolution is not inferred from the constitutive law. Differential equations are considered locally on their declared admissible domain, with explicit uniqueness hypotheses.

## 1. State, branch and response data

Use maps

\[
S\xrightarrow{\ r\ }X\xrightarrow{\ p\ }B,\qquad \pi=p\circ r.
\tag{1}
\]

`B` is an observable presentation, `X` a finite branched state cover, and `S` a supplied smooth resolution or retained-state extension. The resolution can retain continuous directions over a collapsed state. A physical real domain `S_adm` is specified separately from the algebraic cover. Ordinary smooth maps and unbranched state extensions also fit (1).

Over `S`, let `P` be a finite-rank real flow bundle, with effort bundle `P*` and canonical pairing

\[
\mathsf p(e,j)=\langle e,j\rangle.
\tag{2}
\]

The constitutive datum is a fibrewise relation `L_s subset P_s* direct-sum P_s`. A graph is the special case `j=Phi(s,e)`. State rates are attached by an anchor `V_s:P_s -> T_sS` and drift `b_s`:

\[
\dot s=b_s+V_sj.
\tag{3}
\]

This separates response channels from configuration coordinates. The native mobility case is `P=TS,V=id`. A relation need not determine a unique flow; an evolution claim requires the resulting admissible state rate to be single-valued and locally Lipschitz, or another separately established solution theory.

A relation is **passive** when `e(j)>=0` for all its pairs. It is **incrementally passive**, or monotone, when

\[
\langle e-e',j-j'\rangle\ge0
\quad\text{for every two pairs in the same fibre.}
\tag{4}
\]

For a linear relation the two conditions agree. An affine offset preserves (4) but contributes separately to absolute power. A nonlinear response with `Phi(s,0)=0` satisfying (4) is passive. These conditions concern the dual pairing, without selecting a positive metric on the state space.

## 2. Complete linear passive relations

**Theorem BR-1.** In an `n`-dimensional flow space `P`, every maximal linear monotone relation has the unique description

\[
\mathcal E\subseteq P^*,\qquad
K:\mathcal E\longrightarrow\mathcal E^*,\qquad
\operatorname{sym}K\succeq0,
\]
\[
\boxed{
L=\{(e,j):e\in\mathcal E,\ j|_{\mathcal E}=Ke\}.
}
\tag{5}
\]

Here `E` is its actual effort domain. Flow is determined modulo `E^ann`, the annihilator of `E`. The relation has dimension `n`. It is a full admittance graph exactly when `E=P*`. In that case

\[
j=(R+J)e,\qquad R=R^T\succeq0,\quad J=-J^T.
\tag{6}
\]

The cases `E=P*,K=0` and `E={0}` represent respectively zero flow for every effort and zero effort with arbitrary flow. The relation remains defined when a chosen impedance or admittance matrix representation ceases to exist.

**Proof.** For a linear passive relation put `E=pr_e L` and `N={j:(0,j) in L}`. Comparing `(e,j+tn)` for every real `t` forces `e(n)=0`, so `N subset E^ann`. Restriction of any permitted `j` to `E` therefore defines a single linear map `K:E -> E*`; its symmetric part is positive semidefinite. The relation (5) contains the original and remains passive. Maximality gives equality, including `N=E^ann`. Conversely, a monotone extension of (5) must annihilate every vertical vector in `E^ann`, so its effort lies in `E`. Testing against efforts `e+th` for both signs of small `t` forces its flow restriction to equal `Ke`. Thus it already belongs to (5). The dimension is `dim E+dim E^ann=n`. ∎

For a smooth family, (5) is a smooth bundle description on constant-domain-rank strata. The full `n`-plane relation can remain smooth while its effort projection changes rank. Parameter specialization must be performed on the relation, rather than on an inverse matrix formula with a vanishing denominator.

### A regular representation through changes of effort-domain rank

Choose an auxiliary positive identification of effort and flow coordinates; it is a port-unit convention, not a state metric or physical constitutive assumption. Put `w_+=e+j` and `w_-=e-j`. Every maximal linear passive relation has the unique form

\[
\boxed{
w_-=Cw_+,\quad C^TC\preceq I,\qquad
\binom ej=\frac12\binom{I+C}{I-C}w_+.
}
\tag{5a}
\]

Indeed `4e(j)=|w_+|^2-|w_-|^2`. The projection onto `w_+` is injective on a passive relation and, by its dimension, is an isomorphism. This proves (5a) and its converse. It also proves that smooth maximal passive relation families are exactly smooth contraction-matrix families in these coordinates. No effort-domain-rank hypothesis is needed. A graph exists when `I+C` is invertible and then `M=(I-C)(I+C)^{-1}`. Zero-flow and zero-effort relations are `C=I` and `C=-I` respectively.

For example `C(t)=(t^2-1)/(t^2+1)` describes the smooth relation spanned by `(e,j)=(t^2,1)`. Its admittance is `1/t^2` off zero; its value at zero is the ideal effort constraint `e=0`. The divergence belongs to the graph representation. This construction retains the standard linear scattering/Cayley representation and gives a regular coordinate for the relation-level continuation used here.

For a supplied effort in the domain of (5), the state rate is single-valued exactly when

\[
V_s\mathcal E_s^{\rm ann}=0.
\tag{5b}
\]

The possible flows differ by precisely that annihilator. Once (5b) holds, any section of the flow restriction gives the same velocity. If a prescribed effort is outside `E_s`, there is no constitutively admissible rate at that state. An evolution constrained to the compatible states also needs tangency and local uniqueness there; an extra algebraic selection is additional system data.

## 3. Power-preserving transformations and composition

For any linear map `A:P -> Q`, including a rank-deficient one, define

\[
\boxed{
A_\#L=\{(u,Aj):(A^Tu,j)\in L\}.
}
\tag{7}
\]

This is the ideal transformer: effort pulls back and flow pushes forward. Their power is identical. Passivity and incremental passivity are consequently preserved for every relation for which (7) is defined.

**Theorem BR-2.** Transformation (7) preserves maximal linear monotone relations and composes exactly:

\[
B_\#(A_\#L)=(BA)_\#L,\qquad (\operatorname{id})_\#L=L.
\tag{8}
\]

For (5), the transformed effort domain and bilinear response are

\[
\mathcal E_A=\{u:A^Tu\in\mathcal E\},\qquad
(K_Au)(v)=(K A^Tu)(A^Tv).
\tag{9}
\]

Its vertical freedom is exactly `A E^ann=E_A^ann`. This proves (5) for the transformed relation. Equation (8) follows directly by retaining the same intermediate flow witness. For a graph it gives

\[
M_Q=AM_PA^T,\qquad j_{0,Q}=Aj_{0,P}.
\tag{10}
\]

The equality `A E^ann=E_A^ann` follows because `E_A` is precisely the annihilator of `A E^ann`; finite-dimensional double annihilation is exact. Positivity follows by evaluating `K` on `A^Tu`. These arguments prove BR-2 without invertibility or a metric. ∎

### Changes of presentation and global coherence

For a state diffeomorphism `g` and invertible port change `T_s`, the local data transform as

\[
j'=Tj,\quad e'=T^{-T}e,\quad
L'_{gs}=\{(T^{-T}e,Tj):(e,j)\in L_s\},
\]
\[
V'_{gs}=Dg_sV_sT_s^{-1},\qquad b'_{gs}=Dg_sb_s.
\tag{11}
\]

Relations satisfying (11) on overlaps glue to a global relation when the state and port transitions obey their cocycle laws. Graphs glue by congruence (10). Nonlinear graphs obey `Phi'(gs,e')=T_s Phi(s,T_s^Te')`; their effort derivatives transform by the same congruence.

For a state map `phi:S_1 -> S_2`, a port map `A:P_1 -> phi*P_2` defines a response-compatible reduction when

\[
A_\#L_{1,s}=L_{2,\phi(s)},\qquad
D\phi_s b_{1,s}=b_{2,\phi(s)},\qquad
D\phi_s V_{1,s}=V_{2,\phi(s)}A_s.
\tag{12}
\]

Then full trajectories with pulled-back effort project to reduced trajectories. Composing two such pairs gives `(phi_2 phi_1,(phi_1* A_2)A_1)`; (8) and the chain rule preserve every condition in (12). If the reduced evolution is unique, this is an exact predictive reduction for those inputs and outputs. These are sufficient typed conditions on general port systems, not a claim that every output reduction must identify the entire port relation.

An action assigned to an arrow must agree for equal arrows. In the current selected skeleton, identity is `(source,target,monodromy_parity)`; route and provenance digests remain evidence. Material history enters the state when the response depends on it. It is not introduced through metadata-dependent equality.

## 4. Exact smooth lifting through a singular map

Specialize to state-rate ports and a smooth map between equal-dimensional charts, with Jacobian `A=Dpi`. Assume `det A` is nonzero on a dense open set. Let `M,c` be the lower response and drift pulled back to the upper chart; they may also be proposed lower-coordinate components over retained upper state.

The same algebraic criterion holds for any smooth generically invertible square port map. Interpreting such a matrix as a state-coordinate Jacobian additionally requires its actual integrability; no coordinate map is inferred from matrix invertibility alone.

Put `Delta=det A` and `C=adj A`.

**Theorem BR-3.** A smooth lifted graph and drift exist exactly when the entries satisfy smooth ideal divisibility

\[
CMC^T\in(\Delta^2)\operatorname{Mat}_n(C^\infty),\qquad
Cc\in(\Delta)(C^\infty)^n.
\tag{13}
\]

The lift is unique and is

\[
\boxed{
\widehat M=\frac{CMC^T}{\Delta^2},\qquad
\widehat c=\frac{Cc}{\Delta},\qquad
M=A\widehat M A^T,\quad c=A\widehat c.
}
\tag{14}
\]

**Proof.** Where `A` is invertible, power preservation forces (14). A smooth extension exists precisely when the displayed quotients are smooth. On the dense regular set, the reconstruction identities hold; continuity extends them everywhere. Any two smooth extensions agree on that dense set, hence agree everywhere. Symmetric positivity and skew symmetry transfer by congruence off the singular set and by continuity onto it. ∎

This criterion applies to a general Jacobian, including mixed minors. Orders of individual entries of `M` need not describe it in an unadapted frame. Divisibility in (13), not a finite collection of path limits, is the smooth condition.

If a factorization `A=L diag(q_i) R` by smooth invertible `L,R` and monomials `q_i` is supplied, put `N=L^{-1}ML^{-T}` and `v=L^{-1}c`. Then the criterion becomes

\[
N_{ij}\in(q_iq_j),\qquad v_i\in(q_i).
\tag{15}
\]

For a normal-crossing monomial, divisibility is equivalent to vanishing of the requisite coordinate-normal jets along every component, with smooth quotients. Such diagonal factorizations are not presumed for arbitrary multivariable Jacobians; (13) remains the general square-matrix criterion.

For arbitrary rank or rectangular `A`, (7) gives the exact transformed relation but an inverse lift need not exist or be unique. This is a change in the mathematical problem, resolved by retaining the relation and the missing internal directions rather than assigning an inverse.

### Refinements and branch-coordinate examples

For composable generically invertible smooth maps, lift first through one and then the other. If all the lifts are smooth, (8) and dense-set uniqueness show equality with lifting through the composite. Two resolutions with a supplied common refinement therefore agree there whenever their lifted tensors extend smoothly. Extra refinement does not automatically preserve smoothness: a constant positive mobility on a collapsed plane can acquire poles under a blow-up. Equation (13) determines the compatible class of resolutions.

For `x=s^k`, `k>=2`,

\[
\widehat M(s)=\frac{M(s^k)}{k^2s^{2k-2}},\qquad
\widehat c(s)=\frac{c(s^k)}{ks^{k-1}}.
\tag{16}
\]

Thus a scalar leading term `M(x)=x^m unit(x)` lifts smoothly exactly when `km>=2k-2`; this assertion concerns that nonzero leading monomial. For `k=2`, the one-sided base response `M(x)=x` lifts to the constant `1/4`. A constant positive base response does not have a smooth lift. A constant upper mobility pushes to `k^2s^{2k-2}`; its expression in `x` may have fractional regularity. Passing (13) and descending to the desired lower regularity are separate operations.

## 5. Divisor laws, including one-sided boundaries

Let a resolved chart have two-sided coordinates `f_a` and boundary coordinates `r_b>=0`, with independent coordinate differentials. For a graph `M=R+J`, prescribe which outputs vanish for every effort on which faces. An adapted port frame records two index sets for each channel: `I_i` for the `f` faces, and `H_i` for the `r` faces.

The output covectors are part of the input data. For arbitrary varying constraint covectors `C_a`, the intrinsic condition is `C_a M|_face=0`. Positivity separates it into `C_a R|_face=0` and `C_a J|_face=0`: the vanishing quadratic form `C_a R C_a^T` annihilates `R` in those directions. A common adapted coordinate incidence is an additional local property, not presumed for an arbitrary arrangement of subspaces.

### Two-sided interfaces

With `q_i=product_(a in I_i) f_a`, the existing constitutive theorem gives the complete class

\[
R=Q\bar RQ,\quad \bar R\succeq0,\qquad
J_{ij}=\Bigl(\prod_{a\in I_i\cup I_j}f_a\Bigr)\bar J_{ij},
\quad \bar J=-\bar J^T.
\tag{17}
\]

All barred functions are smooth. A common face contributes twice to dissipation and once to the skew response. The proof uses a nonnegative smooth diagonal entry's second-order zero and the positive-semidefinite Cauchy–Schwarz inequality, followed by smooth division. The full proof is retained in [constitutive divisors](2026-09-05-cella-constitutive-divisors.md#1-complete-linear-classification).

### Complete mixed-boundary extension

A one-sided passive scalar can vanish to first order: `R(r)=r` for `r>=0`. Its two-sided analogue `R(f)=f` is not passive on both sides. To obtain the whole mixed class, use the auxiliary real root cover

\[
r_b=t_b^2,\qquad
q_i(f,t)=\prod_{a\in I_i}f_a\prod_{b\in H_i}t_b,
\qquad Q=\operatorname{diag}(q_i).
\tag{18}
\]

This cover is a mathematical construction for the boundary coefficient law; its extra signs are not automatically physical state. Let `E_b` be diagonal, with entry `-1` for channels constrained on `r_b=0` and `+1` otherwise.

**Theorem BR-4.** Every smooth passive graph with the declared mixed face constraints, and only such a graph, has

\[
\boxed{
R(f,t^2)=Q(f,t)\bar R(f,t)Q(f,t),\quad
\bar R=\bar R^T\succeq0,
\quad \bar R(f,t^{(b)})=E_b\bar R(f,t)E_b,
}
\tag{19}
\]
\[
\boxed{
J_{ij}(f,r)=
\Bigl(\prod_{a\in I_i\cup I_j}f_a\Bigr)
\Bigl(\prod_{b\in H_i\cup H_j}r_b\Bigr)\bar J_{ij}(f,r),
\quad\bar J=-\bar J^T,
}
\tag{20}
\]

with smooth barred coefficients. Here `t^(b)` reverses only `t_b`. The parity condition in (19) is essential.

**Proof.** Pull a given `R` back by `r=t^2`. It is smooth and nonnegative for all signs of `t`, and its constrained rows vanish on the associated coordinate hypersurfaces. Apply (17) to obtain `Q bar R Q`. Root reversal leaves the original pullback fixed and sends `Q` to `E_b Q`; uniqueness off the divisor, then continuity, gives the parity rule. Conversely, that rule makes `Q bar R Q` even in every root coordinate. A smooth even function of `t` is a unique smooth function of `t^2` on the half-line, also with smooth parameters: its Taylor expansions have only even terms, and the remainder estimates after differentiation give all one-sided derivatives. Iterating this fact gives smooth descent to the corner. Positivity and the required row zeros follow directly. Skew symmetry transfers a row zero to its column; successive one-sided or two-sided smooth division gives exactly (20). ∎

In particular the complete diagonal factors are

\[
R_{ii}=\Bigl(\prod_{a\in I_i}f_a^2\Bigr)
        \Bigl(\prod_{b\in H_i}r_b\Bigr)\rho_i,\qquad\rho_i\ge0.
\tag{21}
\]

For off-diagonal entries, every common one-sided face needs at least one power, and a one-sided face constrained on only one row also needs at least one power by smoothness and that row's zero. These lower bounds alone are not sufficient for matrix positivity; (19) retains the full positive matrix and its parity compatibility.

Changing a defining function by a nonzero smooth unit rescales the factors and inversely rescales the barred matrix by congruence. Thus the statements describe sections of the corresponding divisor ideals; they do not require a global choice of defining functions. A one-sided boundary defining function changes by a positive unit. On overlaps, combine these rescalings with (11).

### Nonlinear responses

On a convex effort domain containing zero, write `j=j_0(s)+Phi(s,e)` with `Phi(s,0)=0` and `M=partial_e Phi`. The complete smooth incrementally passive class is obtained by imposing (17) or (19)–(20) on `M(s,e)` and the row integrability conditions

\[
\partial_{e_k}M_{ij}=\partial_{e_j}M_{ik}.
\]

Then

\[
\Phi(s,e)=\int_0^1M(s,te)e\,dt.
\tag{22}
\]

Closed row forms prove reconstruction; integrating the symmetric-positive Jacobian along effort segments proves incremental passivity. A universally vanishing total output also requires the corresponding `j_0` to vanish on its face. Potentials are a permitted reciprocal subclass, not a requirement on general response.

For a linear relation (5) with a restricted effort domain, a covector `ell` annihilates **every** permitted flow exactly when `ell in E` and `(Ke)(ell)=0` for all `e in E`. The first condition removes vertical reaction flow; the second removes the constitutive flow. Divisor graph formulas apply to a smooth adapted domain subbundle when their hypotheses hold.

## 6. Algebraic branching and dynamical descent

Cella's general conjugate-Kummer construction supplies, for a finite separable extension with normal closure `E`, radicands `r_(alpha,i)` and their square-relation module `R_rel`,

\[
W=\operatorname{span}_{\mathbb F_2}[r_{\alpha,i}]=V/R_{\rm rel},
\qquad \operatorname{Gal}(L/E)=W^\vee.
\tag{23}
\]

At a discrete valuation with residue characteristic different from two, tame inertia acts on a radical by the sign `(-1)^v(r)`. This classifies algebraic branch transport; the actual real admissible paths must still be specified. General finite covers can be used without a radical presentation or a symmetric base group. Relations and normalization are computed before projecting a defining polynomial: projected root collisions need not be branch points of the retained incidence cover.

Suppose a deck group acts on `S` and on ports by `T_g`. A constitutive law is compatible with that presentation action precisely when

\[
L_{gs}=\{(T_g^{-T}e,T_gj):(e,j)\in L_s\},\qquad
V_{gs}T_g=Dg_sV_s,\quad b_{gs}=Dg_sb_s.
\tag{24}
\]

For state-rate graphs these become `M(gs)=Dg M(s)Dg^T` and `b(gs)=Dg b(s)`. For a fixed-input comparison, the supplied input itself must be invariant; otherwise (24) relates two transformed experiments.

For a supplied finite group action these compatible passive graphs can also be constructed explicitly. Given any smooth passive seed `M_0`, put

\[
\mathcal A M_0(s)=\frac1{|G|}\sum_{g\in G}
T_g(g^{-1}s)M_0(g^{-1}s)T_g(g^{-1}s)^T.
\tag{24a}
\]

Congruence and summation preserve the positive symmetric part. Reindexing the sum proves (24), and every equivariant graph is fixed by `A`; hence the image is the entire compatible passive class. The same averaging of transformed graph maps preserves nonlinear incremental passivity. Declared output constraints are preserved when their family is group-invariant. This construction supplies a class of compatible laws; choosing a seed or identifying it with a physical material remains constitutive input.

**Theorem BR-5 — no finite passage through a fixed stratum under invariant drive.** If the full controlled vector field is locally Lipschitz and equivariant under a subgroup `H` for the applied input, `Fix(H)` is invariant. A trajectory starting outside that fixed set cannot reach it in finite time while those hypotheses hold.

**Proof.** For a fixed point `s`, transforming its trajectory by any `h in H` produces the same initial-value problem, hence the same trajectory. Thus it remains fixed. Local backward uniqueness at a proposed finite arrival would force the incoming trajectory to have been fixed as well. The same argument holds across finitely many input switches preserving the action. ∎

For example `x=s^2`, a smooth even mobility and finite base effort give `dot s=2s Mhat(s)e_x`. The root's fixed set is invariant. A finite upper-state effort driving `dot s=1` is a different input: its base effort would be `1/(2s Mhat)` where that expression exists. Branch-sensitive drive can distinguish or traverse states that invariant drive cannot identify dynamically.

### Descent and its regularity

For observations `y=pi(s)` and base efforts pulled back by `Dpi^T`, the exact graph coefficients are

\[
\dot y=c(s)+B(s)e,\qquad
c=D\pi\,b,\quad B=D\pi\,M\,D\pi^T.
\tag{25}
\]

They define an instantaneous response on `y` exactly when both are constant on the observation fibres and possess the requested regularity there. A locally Lipschitz descended field then gives unique reduced evolution. On a finite unbranched covering whose deck group is transitive on every fibre, deck covariance is equivalent to this tensor descent; at critical points regularity must also be checked. For a general observation with continuous fibres, response covariance alone does not establish predictive reduction.

The exact minimum is the quotient by equality of **all future outputs under the same admissible input**. The [predictive-state theorem](../../../../results/2026-09-05/2026-09-05-cella-predictive-state.md) gives its finite closure certificate and a terminating paired-state ideal construction for polynomial laws. It also proves why equality of all local jets is insufficient for arbitrary smooth flat laws. Continuous exceptional directions, discrete branch alternatives and internal memory are retained precisely when this test distinguishes their futures.

## 7. Internal-port elimination and specialization

**Theorem BR-6.** Direct sums, power-preserving transformations (7), and zero-flow elimination of internal ports preserve maximal linear monotone relations. For the last operation, swap effort and flow, apply (7) with the projection onto external ports, and swap back. The power pairing is unchanged by each swap; BR-2 proves the claim. Consequently any finite network assembled from these operations has an exact passive external relation, including ideal constraint ports.

For an admittance graph partitioned into external and internal blocks,

\[
\binom{j_e}{j_i}=
\begin{pmatrix}A&B\\C&D\end{pmatrix}
\binom{e_e}{e_i},\qquad j_i=0,
\tag{26}
\]

an invertible `D` gives the usual effective matrix

\[
M_{\rm eff}=A-BD^{-1}C.
\tag{27}
\]

Its symmetric part is positive semidefinite because the eliminated internal port has zero power. The singular case has an equally exact classification.

**Theorem BR-7.** If the full matrix in (26) is passive, then the external relation is a full graph exactly when

\[
\operatorname{im}C\subseteq\operatorname{im}D
\quad\Longleftrightarrow\quad B\ker D=0.
\tag{28}
\]

For any generalized inverse satisfying `DGD=D`, that graph is

\[
M_{\rm eff}=A-BGC,
\tag{29}
\]

independent of the choice of `G`. Otherwise let columns of `N` span `ker D`. The effort domain and free external reactions are

\[
\mathcal E_e=\ker(N^TC),\qquad
\mathcal E_e^{\rm ann}=\operatorname{im}(BN).
\tag{30}
\]

Equations (26) and (30) specify the complete relation, which may impose a constraint on efforts instead of a single-valued full admittance.

**Proof.** Write the full matrix as `R+J` with `R>=0`. If `Dk=0`, positivity gives `R_ii k=0`, and hence `J_ii k=0`; this also proves `ker D=ker D^T`. Positivity of the full block gives `R_ei k=0`, so `Bk=-C^Tk`. Solvability of `D e_i=-C e_e` for every `e_e` is therefore equivalent to `Bk=0` for every such `k`. These are exactly full-domain solvability and independence of external flow from the free internal solution. Generalized inverses differ only by such ineffective solutions. Without (28), compatibility is `N^TC e_e=0`, and the external freedom is `BN`, its exact annihilator by the same identity. The zero-power argument proves passivity. ∎

Offsets are retained in the affine equations: replace compatibility by `N^T(C e_e+j_0i)=0` and retain `j_0e` in the external flow. Incremental passivity is preserved. Absolute power must account for these offsets rather than absorbing them into the homogeneous passive claim.

### Elimination can fail to commute with a singular limit

Consider the smooth passive family

\[
M(t)=\begin{pmatrix}1&t\\t&t^2\end{pmatrix}
=\binom1t\begin{pmatrix}1&t\end{pmatrix}.
\tag{31}
\]

For `t!=0`, zero internal flow selects `e_i=-e_e/t` and the external admittance is zero. At `t=0`, the internal flow is automatically zero and the external admittance is one. Thus specializing the full relation and then eliminating differs from first eliminating and taking its limit. The missing assumption is visible: the selected internal effort diverges.

A smooth internal solution operator exists exactly when the columns of `C` lie in the image module `D C^infty^k`, that is, when a smooth matrix `K` solves `DK=C`. With the affine offset, also require a smooth solution of `Dv=j_0i`. If the external relation is a graph at each state, these conditions give a smooth effective response `A-BK` and offset `j_0e-Bv` through the singular stratum. For generically invertible `D`, smooth solvability is the adjugate/divisibility condition `adj(D)C in (det D)` columnwise. This is the internal-port analogue of BR-3.

Cella's split special fibres supply the same operational lesson for state branches: specialize the full retained object, then determine its components and compatible responses. A rank or degree change alone does not specify which component or internal solution evolves.

## 8. Storage, constraints and real crossing

If `H(s)` is a stored quantity and the constitutive effort is

\[
e=e_{\rm ext}-V_s^T dH_s,
\]

then every passive relation gives

\[
\boxed{
\dot H= dH(b)+\langle e_{\rm ext},j\rangle-\langle e,j\rangle
\le dH(b)+\langle e_{\rm ext},j\rangle.
}
\tag{32}
\]

External control, material exchange and imposed drift remain explicit. No conservative energy, potential or kinetic metric is selected by the abstract response law. When storage descends through a response-compatible state map, (32) descends as well. Under a zero-power interconnection, the internal power cancels and subsystem storage inequalities add.

Let `f=0` be a regular hypersurface in `S`, and write a selected controlled graph as `dot s=b+sum u_i X_i`. It is invariant for all inputs spanning the declared control space exactly when

\[
df(b)=0,\qquad df(X_i)=0\quad\text{on }f=0.
\tag{33}
\]

Tangency is necessary. In local coordinates with `f` one coordinate, the vanishing normal rate factors as `f times a locally bounded smooth coefficient`; uniqueness makes the face invariant and excludes finite arrival from outside. For a one-sided admissible face, inward or tangent drift is allowed. With all real control amplitudes, the boundary condition is `df(X_i)=0` and `df(b)>=0`. With a bounded control set, test the actual inequality `df(b+sum u_iX_i)>=0` over that set. These are the local inward-pointing conditions for unique smooth evolution in a coordinate half-space or corner.

**Theorem BR-8 — a zero-output face need not be a barrier.** Let the constrained flow covectors be `c_alpha in P*`. If `V^T df` belongs to their span on the face and `df(b)=0` there, their universal vanishing implies (33). A zero output in another direction need not do so. This follows by applying that linear combination to every allowed flow. The exact test is (33), using the supplied constraints and anchor.

For example, on state `(f,z)` take a two-sided passive mobility and finite input

\[
M=\operatorname{diag}(1,f^2),\qquad e=(1,1).
\]

Then `f(t)=t`, `z(t)=t^3/3` crosses the interface regularly while the `z` output vanishes there for every effort. If the universally vanishing output were instead the normal rate `dot f`, smooth unique dynamics would preserve the face. This distinction determines which constitutive channel an open-network interpretation constrains.

## 9. Native application and closure of the construction

The existing native finite-direction chart has

\[
x=d\xi,\quad F=\frac{1-\xi^2}{1+a^2},\quad
D\pi=\begin{pmatrix}d&\xi\\0&1\end{pmatrix}.
\tag{34}
\]

BR-3 recovers the weighted lift exactly. Over `d=0`, every smooth upper graph has pushed incremental response

\[
M_{\rm lower}(0,\xi)=\widehat M_{dd}(0,\xi)
\begin{pmatrix}\xi^2&\xi\\\xi&1\end{pmatrix}.
\tag{35}
\]

When the radial coefficient is nonzero, the ratio of its off-diagonal to radial entry recovers `xi`. Thus the response exists smoothly on the resolved space while its collapsed presentation omits predictive state. This is an exact instance of (25), not a failure of the constitutive law.

The constant native response remains an unconstrained member at `F=0`, as recorded by OC-1. When a particular load is required to have a vanishing output there, choose the actual output covector and apply BR-4 on the resolved chart. Near `xi=+/-1`, `dF` is nonzero, so the two-sided divisor law is available even at `d=0`. It concerns the retained state; raw `F=epsilon/delta` is still undefined at the collapsed pinch.

For an area output, the existing normal-geometry construction gives

\[
\theta_L=D+\rho F/2,\qquad \rho>0.
\tag{36}
\]

The constitutive area rate `D` and its drift come from the selected cut and work pullback. Universal vanishing of `D` on `F=0` does not mean `dot F=0`. BR-8 therefore permits finite horizon crossing with a constrained area port, while identifying the stronger normal-rate constraint that would prevent it. The actual coframe, cut, material coefficients and drives remain the model's constitutive data.

### Ruler exchange: two depth charts and the resolved action

[HG-1](../../../hg1_hbar_G_is_the_sheet.py) supplies the frame-level input: the actual ruler exchange is `(l_1,l_2,t)->(l_2,l_1,-t)`, reversing the signed frame volume `d`, while the lift with `+t` preserves it. Both induce the same swap of Gram coordinates. Its two depth-normalized directions are

\[
\xi_h=\frac{b-a\gamma}{d},\qquad
\xi_g=\frac{a-b\gamma}{d},\qquad
F=\frac{1-\xi_h^2}{1+a^2}=\frac{1-\xi_g^2}{1+b^2}.
\tag{37}
\]

The actual swap sends `xi_h` to `-xi_g`. This identifies its orientation character and provides the transition between the two normalizations. Transport of a response between these charts is (11), including the effort transformation. Equal numeric coefficient matrices in different depth charts need not represent the same law.

**Theorem BR-9 — resolved ruler exchange and its fixed-stratum consequence.** On either local chart with `gamma=sigma sqrt(1-Fd^2)`, the actual swap extends smoothly as

\[
\boxed{
S(a,\xi,d)=(a\gamma+\xi d,\ \gamma\xi-aFd,\ -d).
}
\tag{38}
\]

It preserves `F` and `gamma`, satisfies `S^2=id`, and commutes with the deck `tau(a,xi,d)=(a,-xi,-d)`.

**Proof.** Substitute `b=a gamma+xi d` and `1-gamma^2=Fd^2` into `(a-b gamma)/(-d)`. This gives the second component without dividing by `d`. Direct expansion gives

\[
1-(\gamma\xi-aFd)^2=F\bigl(1+(a\gamma+\xi d)^2\bigr).
\]

The same identities give `S^2=id`; substituting `(-xi,-d)` verifies commutation. They hold at `d=0` by smooth extension. ∎

On the regular part of the two horizons, `F=0,xi=+/-1`, this becomes

\[
S|_{H_\sigma}(a,\xi,d)=(\sigma a+\xi d,\ \sigma\xi,\ -d).
\tag{39}
\]

Thus HG-1's directional split is exact: on `sigma=-1`, the pair `(xi,d)` transforms by the deck; on `sigma=+1`, `xi` is fixed and `d` reverses. Equation (39) also records the simultaneous depth change. This matters when transporting a law previously studied at fixed `a`: the full ruler swap acts between those depth slices. Its orientation character and its action on all retained coordinates are both now specified.

Over the collapsed pinch,

\[
S(a,\xi,0)=(\sigma a,\sigma\xi,0).
\tag{40}
\]

On the positive branch the entire exceptional set `d=0` is fixed by `S`. A locally Lipschitz swap-equivariant vector field under swap-invariant drive therefore cannot reach or cross it in finite time, by BR-5. On the negative branch only `a=xi=d=0` is fixed in this chart; the horizon points `xi=+/-1` are not. There is no corresponding fixed-stratum prohibition there. This is a conditional distinction between the two resolved pinches, beyond the static equality of their lapse.

The condition on the drive is explicit. At the positive exceptional set,

\[
DS=\begin{pmatrix}1&0&\xi\\0&1&-aF\\0&0&-1\end{pmatrix},\qquad
S^*e=e\quad\Longleftrightarrow\quad
2e_d=\xi e_a-aF e_\xi.
\tag{41}
\]

An equivariant state-rate law with this invariant effort has zero `d` velocity there. An orientation-sensitive effort is not subject to that conclusion. A constant response in the previously fixed-`a` chart is not asserted to satisfy full swap covariance; (24) or the construction (24a) determines that separately.

HG-1's regular identity

\[
S\xi-\tau\xi=\frac{(b-a)(1+\gamma)}d
\tag{42}
\]

also has an exact resolved interpretation. Near `gamma=+1`, the closure of the regular equal-depth plane is

\[
\xi=\frac{aFd}{1+\gamma},
\tag{43}
\]

so it meets the exceptional set at `xi=0`. The equality of the two direction maps continues on that strict transform. The full exceptional fibre contains additional directions; (38) specifies their action rather than evaluating `0/0` in (42). This keeps the observed equality on equal-depth states through the resolution while preserving the extra directional state.

Finally the response makes the orientation dynamically distinguishable for the stated static readings. In the known diagonal native law `Mhat=diag(mu,nu)`, `mu,nu>0`, prepare `(xi,d)=(+1,0)` or `(-1,0)`. Apply the original `d` effort for time `T`, then the original `x` effort. The instantaneous lapse rates in the second stage are

\[
\dot F=-\frac{2\mu\nu T}{1+a^2}\,\xi,
\tag{44}
\]

with opposite signs. The same input distinguishes the deck partners despite equal initial `F`. This conclusion concerns the declared scalar readings and response experiment; adding a signed-volume or ordered-frame measurement changes the output set and therefore its predictive quotient. Identifying a ruler convention with an experimental port requires the transformation in (11), now concretely available from (38).

The construction is complete for the stated classes: maximal finite-dimensional linear passive relations, including smooth changes of effort-domain rank; smooth graph responses and their nonlinear integrable extensions; supplied branched covers and resolutions; mixed normal-crossing faces in an adapted port frame; graph lift, global covariance, passive network reduction, and unique controlled evolution with its predictive quotient. BR-9 gives the native ruler action and its controlled crossing consequence. Nonlinear relations with nonunique rates, nonsmooth dynamics, or a field theory require their own solution theory; the pointwise power and transformation identities remain applicable.

## 10. Proof sources, computation and retained background

The [exact runner](../../../branched_resolved_response.py) implements linear-relation transformation and internal-port elimination, graph lifting and polynomial checks of its divisibility conditions. Its [receipt](../../../../receipts/branched-resolved-response-checks.json) includes noninvertible maps, composition, mixed boundary parity, nonlinear effort integration, singular elimination before/after specialization, and native crossing examples. Symbolic checks verify these constructions and counterexamples; the general claims are supported by the proofs above.

The Cella inputs are:

- [Conjugate Kummer modules](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/galois_horizon_and_kummer_covers/KUMMER_MODULE_WREATH_LIFT_THEOREM_2026-07-10.md>) and its [general criterion](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/galois_horizon_and_kummer_covers/GENERIC_SYMMETRIC_MONODROMY_OF_WEIGHTED_MULTIQUADRATIC_SUMS_v1.0.md:1326>): the actual square-relation module and inertia action.
- [Normalized-incidence workflow](</home/williaml/Cella Framework/Papers_Library/05_expository_companions_and_research_maps/galois_horizon_and_kummer_covers/MACAULAY2_REALIZATION_POSET_WORKFLOW_2026-07-10.md>) and [split special fibre](</home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/galois_horizon_and_kummer_covers/ROTATING_KUMMER_RANK_JUMP_LEMMA_REPORT_2026-07-10.md:224>): preserve the full state before projecting or specializing.
- [Selected Quotient Groupoids](</home/williaml/Cella Framework/Papers_Library/05_expository_companions_and_research_maps/dbp_role_channel_and_orbit_geometry/SELECTED_QUOTIENT_GROUPOIDS_FOUNDATION_v1.0.md>), DAG `DBP:thm:sqg_foundation`, and the [selected-skeleton implementation](</home/williaml/Cella Framework/engine/src/cella/continuation/selected_skeleton.py:58>): typed transports, composition and evidence-independent arrow identity. Equations (11)–(12) give the response representation and its compatibility conditions.
- [Cella local divisor calculus](</home/williaml/Cella Framework/Papers_Library/01_completed_papers/local_curvature_and_black_hole_metrics/LOCAL_CURVATURE_CALCULUS_COMPLETE_v1.0.md>) and its [tensor extension](../../../../results/2026-09-05/2026-09-05-cella-tensor-valuation.md): full coefficient and divisor data, rather than only leading exponents, determine smooth continuation.
- DIS [constitution and history](</home/williaml/Cella Framework/research/derivations/DIS_chapters/10_bipolar_axes_and_constitutive_resistance.md>) and [stored-energy programme](</home/williaml/Cella Framework/research/derivations/DIS_chapters/12_internal_capacity_energy_and_final_hypotheses.md>), together with [constitutive divisors](../../../../results/2026-09-05/2026-09-05-cella-constitutive-divisors.md), [predictive state](../../../../results/2026-09-05/2026-09-05-cella-predictive-state.md), [frame transport](../../../../results/2026-09-05/2026-09-05-cella-frame-transport.md) and [normal geometry](../../../../results/2026-09-05/2026-09-05-cella-normal-geometry.md).

Standard dual linear algebra, passive relations, smooth division and ODE uniqueness are retained; so are the previously recorded Kummer and polynomial-closure tools. No online search or new external source was used. The relation calculus itself is not claimed as new passive-systems mathematics. The delivered extension is its complete integration with mixed boundary divisors, resolved lifting, specialization and native predictive/crossing data, with all comparison hypotheses explicit.

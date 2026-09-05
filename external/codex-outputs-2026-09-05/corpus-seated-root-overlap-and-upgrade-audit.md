# Cella corpus × Seated Root: overlap, continuation and upgrade audit

> **Library follow-up, 5 September:** the [mathematics-library addendum](mathematics-library-upgrade-addendum-2026-09-05.md) supplies later corrections to this assessment: the existing categorical implementation fails its equality laws; finite-order naturality and the four-charge rotating closure are already banked; the library already recognizes the squared-carrier and rounding-scope limitations. It also sharpens the proposed parity/period bridge.

5 September 2026. This extends the session’s mathematical-foundation audit to the Cella corpus. The objective is a model capable of earning new physical predictions, especially about singularities and black holes. Both bodies of work are candidates for improvement. A limitation identifies a continuation problem; a proved obstruction identifies which hypothesis or construction that continuation must change.

## 1. Verdict

**The corpus materially improves the available starting point.** Several ideas recommended in the first audit already have substantial development in Cella: active role-dependent observations, quotient/selection separation, typed boundaries, invariant recovery, curvature valuations, normalized covers and explicit continuation. They should be evaluated as existing mathematical assets before creating replacements.

The strongest route is a **general framework for observations, selections and evolution across changing strata**, with Seated Root as one realization. Cella supplies much of the observation and comparison language. Seated Root supplies a concrete signed-frame geometry on which to develop and test it. Neither currently supplies a uniquely selected physical evolution law connecting these structures to black-hole formation or singularity resolution.

There are three especially promising continuations:

1. **Singularity and observation valuations.** Combine the corpus’s asymptotic machinery with Seated Root’s exact Gram identities. A first classification is derived in §3 below.
2. **Selection and transport across boundaries.** Extend the selected-quotient framework to include actual state observations, changes of stratum and admissible evolution. Preserve the existing wall-isotropy obstruction rather than identifying different boundary types.
3. **Reconstruction of the quadratic response.** The corpus already produces a signature-(2,1) quadratic form from a curvature-channel calculation. Investigate whether an independently justified observation law connects that construction to Seated Root’s form. Equal signatures alone do not establish this connection.

The corpus also contains findings needing repair or sharper scope. I give exact examples and a continuation route for each. The result is not “use the corpus unchanged” or “discard the restricted parts.” It is a map of what can be reused, what needs a bridge, and what could become stronger mathematics.

## 2. What was inspected

The Cella DAG at revision `575f1e672b44c5303a308fa04fef07fa46d61f1f186a65fad1396507d850871b` contains 2,300 nodes and 4,335 edges. I screened the complete node index, searched mathematical topics across labels and summaries, inspected relevant dependency edges, and retrieved 17 distinct source files through `cella_dag_source`. All 17 matched their declared SHA-256 digests and were returned without truncation. Relevant definitions, propositions and proof sections were read; this is not independent verification of every proof in those files or all 2,300 nodes.

The node index references 167 distinct source pointers after removing line anchors. It is a discovery map, not an exhaustive inventory of every file or a guarantee that every claim is true. For example, `U-1918` remains open while a campaign has a `proves` edge to it whose basis concerns only the proved wall-isotropy fork. The underlying statement and scope determine what transfers.

On the Seated Root side I used the current primitives and runners, the earlier session audit, and the newer repository-local foundation audit. The live primitives now explicitly implement your ruling: **nothing is banned; every dependency must be derived or declared with its mathematical consequences.** That rule governs this comparison. The repository-local audit’s suggestions are evaluated as suggestions, including its claim that only one new mathematical development is needed and that the dynamics can already be identified with the Kähler dial. The corpus does not establish either of those stronger conclusions.

The output is an assessment and a proposed research direction. No corpus files, Seated Root files, DAG claims or statuses were changed.

## 3. A worked bridge: classification of lapse versus rank loss

This is an immediate mathematical gain from putting the two bodies of work together. Seated Root provides the exact identity; Cella’s valuation approach provides the right question to ask of it.

For

\[
G=\begin{pmatrix}-1&a&b\\a&1&\gamma\\b&\gamma&1\end{pmatrix},
\qquad \varepsilon=1-\gamma^2>0,
\qquad d=a-b\gamma,
\]

the compact-ruler sector satisfies

\[
-\det G=(1+b^2)\varepsilon+d^2,
\qquad
N^2=\frac{\varepsilon}{(1+b^2)\varepsilon+d^2}.
\]

These follow directly by expanding the determinant and completing the square in the Schur-complement numerator. [Seated Root identity][S1]

**Finite-endpoint valuation proposition.** Suppose, as \(s\downarrow0\),

\[
b\to b_0\in\mathbb R,\quad
\varepsilon\sim E s^p,\quad
d\sim D s^q,
\qquad E>0,\ D\ne0,\ p,q>0.
\]

Then all paths approach rank loss, but their lapse limits differ:

- If \(2q<p\), then
  \[
  N^2\sim\frac{E}{D^2}s^{p-2q}\to0.
  \]
- If \(2q=p\), then
  \[
  N^2\to\frac{E}{(1+b_0^2)E+D^2}.
  \]
- If \(2q>p\), then
  \[
  N^2\to\frac1{1+b_0^2}.
  \]

**Proof.** Compare the two positive denominator terms. Their orders are \(p\) and \(2q\); no cancellation is possible. Both vanish, so \(\det G\to0\). The limiting matrix has rank exactly two because its \((c,\hbar)\) principal minor is \(-1-a^2\), which never vanishes at a finite endpoint. Dividing numerator and denominator by the dominant power gives the three limits.

If instead \(d\to d_0\ne0\), the lapse vanishes while \(-\det G\to d_0^2\): this is the familiar off-rank-loss case. Paths with unbounded \(b\), non-power asymptotics or other observation rules need separate classifications.

**What this earns.** The previous counterexample becomes a whole classification of finite power-law approaches. It supplies exact targets for admissible evolution laws: which of these regimes can their solutions reach? It does not select a path, equate rank loss with a physical singularity, or establish a horizon. Those require dynamics and physical reconstruction.

**Continuation.** Extend the classification to unbounded paths, multiple vanishing minors, signed-volume observables and actual solution curves. Then express rates in a physically justified clock or invariant path parameter. The algebra above is elementary and makes no historical novelty claim; its systematic extension is a credible research target.

## 4. Active role jets and observation completeness

**Overlap: strong.** Cella’s active role calculus solves one implicit relation with different output coordinates and transports its derivatives. For \(P=f(D,S)\), with \(a=f_D,b=f_S\), one transposition gives

\[
(a,b)\mapsto(1/a,-b/a).
\]

Higher jets transform by the inverse-function theorem. The action is defined only where the required output derivative is nonzero. The orbit-quotient theorem and the separate requirement of orbit separation are appropriate foundations for reference-dependent readings. [Paper I, §§2–5][C1]

**Fit and transfer condition.** A Cella role change is an alternate graph chart of the same relation. A Seated Root pivot acts on its frame/reference data and may alter which information is accessible. To identify the two, construct a state-to-jet map and show that the relevant operations commute with it. Sharing three roles or an \(S_3\) label does not provide that map.

The corpus proves a useful same-gradient gauge-quotient result in dimension three: its full order-1/2 RoleChSpec distinguishes Hessians modulo \(H\mapsto H+gv^T+vg^T\), on the stated regular locus. This is different from recovering a whole physical state or proving completeness for arbitrary dimension. [Paper I, §8.2][C1]

**Upgrade: separate local recovery from global sign recovery.** The smaller squared coupling carrier has a nonzero Jacobian generically, but that proves local invertibility. At \(a=b=1\),

\[
O(A,B,C)=-\frac19\bigl(B^2,(A-B)^2,(C-B)^2\bigr).
\]

The jets \((2,1,3)\) and \((-2,-1,-3)\) give the same \(O=(-1,-1,-4)/9\). They are not related by an active role permutation preserving the first jet: only identity and input swap preserve \((a,b)=(1,1)\), and neither gives the negative jet.

Thus the paper’s local result is useful, but cannot be promoted to unrestricted global completeness. Generically, before further identifications, the three independent squarings forget eight sign choices. Retain a sign chamber, add appropriate signed readings, or use a richer verified fingerprint. [Paper I, §8.4][C1]

There is also a distinction between anisotropy and loss of rank. At \((A,B,C)=(2,1,2)\), all three channels agree, so the anisotropy scalar is zero, while the carrier Jacobian is \(8/729\ne0\). Equality of channels is not the same event as a channel vanishing or a recovery map becoming singular.

**Research opportunity.** Determine the minimal observable augmentation that globally separates the intended Seated Root states, then classify what happens when that augmentation fails at a boundary. Cella’s role-jet language makes this a precise recovery problem rather than a qualitative appeal to hidden information.

## 5. The corpus’s Lorentzian coupling form

**Overlap: unusually direct, but not yet a physical identification.** On a Euclidean implicit hypersurface, the pure coupling numerator can be written

\[
\Delta_c=\nu^TQ\nu,
\qquad Q=2I-\mathbf1\mathbf1^T,
\qquad
\nu=(g_1H_{23},g_2H_{13},g_3H_{12}).
\]

The eigenvalues of \(Q\) are \(2,2,-1\). The form arises from the declared curvature-channel construction, with a distinguished coordinate basis and the Euclidean gradient/Hessian machinery. It therefore supplies a concrete origin for an indefinite quadratic observable. [Canonical invariant reduction][C2]

This is worth investigating as an alternative route to Seated Root’s quadratic form. Its value is the explicit construction and dependency ledger, not the word “Lorentzian.” The vector \(\nu\) is a coupling-edge vector; Seated Root uses a frame of lines and a reference-dependent Gram form. Neither a shared signature nor an arbitrary linear isometry identifies their states, allowed operations, clocks or physical cones.

**Continuation.** First construct an observation-preserving map between these carriers, or prove that no such map exists under the proposed operations. Next determine which input assumptions force the curvature-channel observable itself. Finally classify the higher-dimensional channel forms: the corpus already has a general triangle decomposition and a named open signature question, rather than a reason to assume the three-dimensional form remains universal.

One prior corpus obstruction must survive this continuation. The additive gauge group \(\ker\Sigma\) over \(\mathbb Q\) is divisible. Every homomorphism to \(C_2\) is trivial because \(\phi(a)=2\phi(a/2)=0\). That particular gauge action cannot generate the desired parity. A discrete lattice, a disconnected symmetry group or an independently defined branched cover changes the mathematical problem, but each requires its own justification. [Wall-isotropy theorem, §1][C3]

## 6. Selected quotient groupoids and boundary transitions

**Overlap: the strongest existing organizational framework.** Cella already defines

\[
\mathbb X=(\mathcal G,\mathcal A,P,\ell,M,K,s),
\]

with a symmetry/transport groupoid, an admissible subgroupoid, stratum labels, a module representation, a null subrepresentation and a selected quotient class. Composition of native morphisms and flat scalar extension are proved. This is substantially more developed than a proposal to “use groupoids.” [SQG foundation, §§1–4][C4]

Two distinctions matter for an upgrade:

1. The selected natural transformation \(s:R\to M/K\) is a **coefficient class**, not automatically a section choosing a state from an observation fibre. The two selections must be related explicitly.
2. A functor from a groupoid to a poset is constant on each orbit. The current stratum label therefore does not itself describe evolution from one stratum to another. Likewise, invertible transport of a module cannot alone describe a rank-changing loss of accessible information.

These are design boundaries, not defects in the category-formation proof. Extend the framework with actual state/observation maps and explicit transition morphisms where the application needs them. A process category, correspondences or a stratified path construction are candidates. Exit-path theory is an existing model for organizing information over strata, under its own stratification hypotheses; its categorical arrows are not automatically physical time evolution. [Treumann][E1]

**The corpus has already tested a tempting over-unification.** A quadratic ordered-root cover has a fixed point where its roots meet. A two-candidate rounding output can have a free two-point orbit at its selection wall. Their stabilizers differ, so a wall-preserving equivariant isomorphism cannot identify them. The precision state cover and its candidate-output bundle are themselves distinct objects. [Wall-isotropy theorem, §§2–6][C3]

This suggests a stronger objective than insisting that all realizations be equivalent: a framework able to **classify and compare different boundary mechanisms**, preserving exactly what each reduction forgets. Regular-chamber comparison may still exist. Extending it across a wall requires additional structure and a suitable comparison notion, not merely more examples.

**Continuation for Seated Root.** Register its orientation-forgetting map, Clifford-block data and observation-critical sets as separate objects. Construct comparison maps, compute stabilizers and lost information, and identify which maps extend across each stratum. This could turn the current “sheet” language into a precise theorem family.

## 7. Lifted compatibility: a useful general construction with upgradeable limits

**Overlap: strong for diagnostics and state constraints; conditional for physics.** For a smooth constraint map \(F:\mathbb R^n\to\mathbb R^k\), the corpus studies

\[
L(x,X)=\sum_iF_i(x)^2-X^2.
\]

At a regular compatible point, \(F=0\) and \(\operatorname{rank}J=k\), its Hessian is

\[
\operatorname{Hess}L=
\begin{pmatrix}2J^TJ&0\\0&-2\end{pmatrix}.
\]

The rank, nullity and inertia statement follows directly. It applies outside physics without relabelling anything. [GR paper brief, §1][C5]

**A stronger formulation is available.** Since \(F\) is a submersion at a regular compatibility point, take \(y=F(x)\) as part of local coordinates and complete them by \(z\). The lifted zero set is then exactly

\[
\|y\|^2-X^2=0,
\]

times the free \(z\)-directions. This gives a local smooth-coordinate normal form; the Hessian alone only gives its quadratic approximation in the original coordinates. Singular-value cone aspect ratios refer to the declared input/residual metrics, not an invariant of the bare constraint set under arbitrary reparameterization.

**Generalize the residual geometry.** With a positive-definite residual weight \(W(x)\) and \(c(x)>0\),

\[
L_{W,c}=F^TWF-cX^2,
\qquad
\operatorname{Hess}L_{W,c}|_{F=X=0}
=\operatorname{diag}(2J^TWJ,-2c).
\]

The inertia remains \((k,1,n-k)\). The particular eigenvalue \(-2\) is a normalization, and the one-negative-direction choice comes from the scalar negative square in the lift. Replacing it by a plus sign changes the geometry. Thus the construction may support a conditional signature theorem, but does not independently select a physical time dimension.

**Repair the coupling-graph corollary.** The brief says an edge is present whenever some constraint has nonzero derivatives in both corresponding variables. Contributions can cancel: \(F_1=x+y,F_2=x-y\) have full-rank Jacobian and both use both variables, yet \(J^TJ=2I\), so the mixed Hessian entry is zero. Gradient support gives possible edges; the actual weighted inner product determines the edge.

The null space also gives tangents to the compatibility manifold; it is not automatically a physical symmetry group. At nonregular compatibility points the quadratic description becomes incomplete: \(F(x)=x^2\) gives the same lift Hessian at the origin as many higher-order constraints while its branches have different contact behaviour.

**Continuation.** Develop a metric-aware residual construction, exact graph cancellation rules and a higher-jet classification when \(J\) loses rank. For Seated Root, apply it only after defining which compatibility equations represent actual physical constraints. It can expose why a state ceases to satisfy them; it cannot choose those equations or resolve a spacetime singularity by introducing a conical residual coordinate.

## 8. Curvature decomposition, invariants and higher-dimensional growth

**Overlap: useful for interpreting what a chosen reference preserves and loses.** The corpus has an exact decomposition of elementary shape curvatures into channel coefficients, using a polynomial in the diagonal and off-diagonal Hessian parts. The total is invariant for the represented geometry; the split depends on the privileged basis and defining-function gauge. [Canonical reduction][C2]

This is a better foundation for partial readings than treating each channel as independently intrinsic. On a regular Euclidean hypersurface, the tangent restriction of \(H\) removes the gauge term \(gv^T+vg^T\). Extending to arbitrary coordinates requires the metric and covariant Hessian; extending to a null normal needs a different construction because the usual orthogonal projector is undefined there.

The shape-readout arm also contains an all-\(n\) triangle decomposition, with explicit nonzero-gradient and normalization scope. Its low-degree separation results are local or tested on specified strata, not global completeness results. [Arm closeout][C6]

One open item in that dated closeout is already superseded: the later F4 theorem proves that degree-\(d\) invariant monomial orbits of the pair module correspond to multigraphs with \(d\) edges and stabilize once \(n\ge2d\). This is an actual dimension-independent construction to reuse, not a new task to repeat. [F4 theorem][C7]

**Continuation.** Classify global fibres and degeneracies of these invariant maps; determine which extra moments recover discarded shape information; and compare that information with Seated Root’s orientation-sensitive readings. Keep the full role-labelled family where labels represent actual instruments. Scalar summation is a reduction, not permission to erase experimentally distinguishable channel information.

For the newest Seated Root curvature work, the first bridge must identify the geometry: Cella’s hypersurface shape curvature, an inverse-channel metric’s curvature and the frame quotient’s mechanical curvature are different tensors on different spaces. An explicit immersion, pullback or reduction is required before their formulas can be compared. A common divergence exponent does not identify them.

## 9. Curvature valuations and transfer laws: continuation after cancellation

**Overlap: immediately useful for the singularity programme.** The corpus’s pure normal form

\[
ds^2=Bx^2dx^2+x^{-2}\sum_{\alpha=1}^mA_\alpha dy_\alpha^2
\]

has scalar curvature

\[
R=-\frac{m(m+5)}B x^{-4}.
\]

The warped-product proof is short and sound. In proper normal distance \(r=\sqrt Bx^2/2\), it becomes \(R=-m(m+5)/(4r^2)\). This is a useful example of replacing coordinate exponents with an invariantly specified distance. [Normal-form paper, §2][C8]

Its restriction is the chosen metric germ, not the mathematics of curvature. Parity alone does not force the pole: \(x^2dx^2+dy^2+dz^2\) is even and flat on \(x>0\). The divergent transverse warp is essential. To transfer the result to Seated Root, derive its actual metric germ and identify which notion of distance it represents.

**Upgrade the corner theorem after cancellation.** The paper states that smooth unit factors leave the corner vertices and hence the order unchanged. This needs qualification when leading coefficients cancel. Compare

\[
g_0=\operatorname{diag}(x^2,1,1),\qquad
g_1=\operatorname{diag}(x^2,1+x,1).
\]

The positive factor \(1+x\) is a unit at the corner. The first metric is flat, but direct Ricci calculation gives

\[
R[g_1]=\frac{2+3x}{2x^3(1+x)^2},\qquad
x^3R[g_1]\to1.
\]

The pure leading-weight coefficient vanishes; the next jet of the unit creates an order-three pole. This preserves the pure monomial calculation while showing exactly what the extension must track. The later unified-spine draft already formulates its corner result more conditionally. [Corner claim][C8]; [unified spine, §11][C9]

**Continuation.** Build an iterative valuation ledger that records actual nonzero coefficients and unit jets after cancellation. Prove termination on a declared analytic/algebraic class, or characterize failures. Without a specified class, a uniform finite jet depth is not guaranteed. This is a stronger target than a universal pole exponent selected from leading weights alone.

**Transfer-function upgrade.** TFEF v6 correctly derives the amplitude asymptotic

\[
\frac{g(f+\delta f)-g(f)}{\delta f}
\sim\alpha c f^{\alpha-1}L(f)
\]

from its controlled branch assumptions and \(\eta=\delta f/f\to0\). But the further pointwise derivative claim \(D\log T\to\alpha-1\), with \(D=f\,d/df\), needs control of how the sampling step changes. [TFEF v6, Theorem 1][C10]

Take \(g(f)=f^2\) and

\[
\eta(f)=f\bigl(2+\sin(f^{-2})\bigr),\qquad
T(f)=f\bigl(2+\eta(f)\bigr).
\]

Then \(\eta\to0\) and \(T/(2f)\to1\), but along \(f_n=(2\pi n)^{-1/2}\), \(D\log T\to-\infty\). Differentiating a vanishing step error is not justified by its amplitude alone.

Retain the amplitude theorem. A sufficient repair for the derivative conclusion is a differentiable step schedule with \(D\eta\to0\), alongside the existing branch derivative controls. Alternatively formulate a finite-window slope statement with its own explicit sampling assumptions. This matters directly when Seated Root’s boundary exponents are estimated numerically.

The corpus’s precision-separation method is also useful, conditional on its amplitude model \(C_{p,k}=a+u_pb_k\). Vary arithmetic precision, approach path and geometric normalization independently. Precision stability can remove an arithmetic explanation; it cannot by itself distinguish an intrinsic curvature divergence from a coordinate-response divergence.

## 10. Galois covers, monodromy and real continuation

**Overlap: strong and reusable at theorem level.** The weighted multiquadratic paper contains a general decorated-cover construction. Adjoin every conjugate radicand to the normal closure of the base field, form their square-class span \(W\), and obtain

\[
1\longrightarrow W^\vee\longrightarrow\operatorname{Gal}(L/F)
\longrightarrow G\longrightarrow1.
\]

Its valuation-parity matrix tests square-class independence; full independence gives the maximal permutation/Kummer lift. Without maximality, the extension class must still be checked. [Weighted monodromy, §10][C11]

This is directly applicable methodology for Seated Root’s radical closure problem. The required work is to specify its actual acting group, include the orbit of its radicands, and compute their relations. The corpus’s all-\(k\) two-radical theorem is a substantial example of that method. Its physical entropy interpretation is explicitly restricted to the established four-charge case; the algebraic all-\(k\) extension does not supply new physical charge models automatically. [All-k closure, §§8–9][C12]

The horizon-cover paper also distinguishes projected sheet crossings, actual ramification and Kummer rank drop. It computes ramification on an incidence cover instead of using the eliminant discriminant as a universal oracle. That is closely aligned with the earlier Seated Root distinction between a Gram-coordinate fold and an intrinsic degeneration. [Horizon-cover paper, §9][C13]

**Continuation.** Apply the normalization and valuation machinery to the real Seated Root sectors and their intersections. Determine which algebraic loops correspond to admissible real paths, which require complex continuation, and how the allowed operations act on orientation and block labels. Extend over exceptional loci using their actual local geometry. A large generic group alone neither chooses the physical branch nor specifies its dynamical accessibility.

**Periods and continuation.** Cella’s CCE/Paper III has explicit route-dependent relative transport, coefficient-ring typing, compact corrections and controlled pole subtraction. These could strengthen Seated Root’s phase and continuation bookkeeping after a concrete cycle/period map is supplied. The source explicitly limits its native swept lattice to rank at most four inside a rank-twelve surface group, and leaves whole-plane transport and some comparisons open. [Paper III, §§7F–9][C14]

Those limitations are continuation targets: extend the route domain, calculate the native inclusion’s Smith factors, and prove the relevant maps commute. They are not grounds for assuming whole-surface equivalence. Nor does selecting a complex lateral continuation select the future direction of a black-hole interior.

## 11. Physical dynamics and constitutive reconstruction

The corpus’s coupling-field paper asks whether coupling geometry obeys a closed evolution law or action principle; it explicitly does not claim to have found one. Some introductory terminology predates the corrected active/passive role distinction, so the later role-calculus papers are stronger sources for that mathematics. [Coupling-field research paper, §§4–5][C15]

Accordingly, neither “the corpus already contains the dynamics” nor “the Kähler dial is the physical impedance” is established by this overlap. A Hermitian/Kähler structure supplies geometric tensors. It does not fix a generator, response kernel, time parameter, source law or dissipative environment. Such a structure remains an admissible declared input; its physical job must be demonstrated.

There is a constructive way forward. Given a candidate flow \(V\) and observation \(p\), form the successive observable derivatives

\[
p,\quad L_Vp,\quad L_V^2p,\ldots.
\]

Ask whether a finite collection separates the hidden state, or whether the flow descends to an autonomous observed law. Cella supplies finite-jet and recovery techniques; Seated Root supplies concrete hidden-depth examples. Established nonlinear observability theory provides a comparison framework, with local rank criteria under stated assumptions. It does not choose \(V\). [Hermann–Krener][E2]

**Continuation.** Derive or declare a bounded family of evolution/response laws from independently justified operations. Require preservation of admissibility, a defined clock, and a well-posed initial-value problem. Classify which laws reach the regimes in §3 and what their detector readings do there. If no reduced law exists, determine whether extra observations, state variables or memory close it.

For electromagnetic or causal-cone reconstruction, the first audit’s premetric constitutive route remains a candidate. The corpus search did not locate a completed Maxwell/Born/Clifford reconstruction theorem that closes those bridges. This is a scoped search result, not a claim that no relevant unindexed work exists.

The decisive physical comparison is then explicit: identify the spacetime reconstruction and which effective Schwarzschild/Kerr predictions it preserves, alters or ceases to support. A new observable must follow from the selected law and independently fixed inputs, rather than serving as the condition used to select that law.

## 12. Recommended continuation programme

These are proposed mathematical targets, not adopted changes to either project.

**First: complete the singularity/observation bridge.** Start with §3’s exact classification. Add signed volume, the projected-angle map and relevant principal minors. Classify finite and infinite approaches, distinguish coordinate rates from invariantly normalized rates, and prove which regimes an admissible flow reaches. Success is a theorem relating observation failure, rank change and dynamical accessibility. A result merely renaming a determinant zero would not complete it.

**Second: extend SQG with state observations and boundary transitions.** Keep coefficient quotients distinct from state fibres. Preserve the wall-isotropy comparison and allow different realizations to have different boundary types. Construct one explicit Seated Root realization and one non-physics realization with comparison maps. Success is useful common mathematics with precise obstruction statements, not universal equivalence by definition.

**Third: investigate the coupling-form reconstruction.** Determine whether the corpus’s curvature-derived quadratic form can be obtained from Seated Root’s independently motivated operations, and whether those operations justify the quadratic category. Compare nonquadratic alternatives where the same operational data permit them. An arbitrary isometry between two signature-(2,1) matrices is insufficient.

**Fourth: complete global observable recovery.** Extend the local RoleChSpec/shape-readout results by identifying sign information, degeneracy strata and sufficient extra measurements. Then test whether the augmented observations close a dynamical law. Both the corpus and Seated Root benefit from this result.

**Fifth: derive one physical consequence.** After the evolution and detector maps are fixed, target a transition-accessibility criterion, a finite/infinite tidal or detector response, or a source-matched rotational deviation. The choice should be made by what the derived law can distinguish, not by which familiar result is easiest to reproduce.

There is no reason to assume original mathematics can occur in only one place. The strongest opportunities are the connections between these structures, especially at singular strata. There is also no need to manufacture a new branch of mathematics before testing whether existing geometry, invariant theory and dynamics can supply those connections. A precise new theorem or a stronger construction is already meaningful progress toward the model’s goal.

## 13. Verification and remaining limits

An independent SymPy script checks six groups of exact examples: the lapse/rank identity and three valuation regimes; the coupling-form signature; cancellation in the lifted compatibility graph; squared-carrier ambiguity and anisotropy/rank separation; the corner-unit curvature counterexample; and the sampling-step counterexample to the transfer-slope claim. It completed successfully. The mathematical arguments precede those checks.

The 17-source hash record and verification files accompany this report:

- [Source manifest](</mnt/files/My Files/Documents/Codex/2026-09-05/ca/outputs/corpus-seated-root-source-manifest.json>)
- [Exact checks](</mnt/files/My Files/Documents/Codex/2026-09-05/ca/outputs/corpus-seated-root-checks.py>)
- [Check results](</mnt/files/My Files/Documents/Codex/2026-09-05/ca/outputs/corpus-seated-root-checks.json>)

No broad corpus test suite was run, no historical priority claim is established, and the proposed general framework has not been implemented or proved sufficient for new physics. The completed result is the source-backed overlap assessment, explicit upgrade routes and the first worked valuation bridge.

[S1]: /home/williaml/seated-root/prim_t7e_two_degeneracies.py:39
[C1]: </home/williaml/Cella Framework/Papers_Library/01_completed_papers/dbp_role_channel_and_orbit_geometry/dbp_orbit_calculus.tex>
[C2]: </home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/dbp_role_channel_and_orbit_geometry/Canonical_Invariant_Reduction_Theorem.md>
[C3]: </home/williaml/Cella Framework/Campaign_Library/02_residue_role_channel_and_coupling/CMP-0040__threeway_selected_skeleton_coherence/WALL_ISOTROPY_COMPARISON_THEOREM_2026-07-16.md>
[C4]: </home/williaml/Cella Framework/Papers_Library/05_expository_companions_and_research_maps/dbp_role_channel_and_orbit_geometry/SELECTED_QUOTIENT_GROUPOIDS_FOUNDATION_v1.0.md>
[C5]: </home/williaml/Cella Framework/Papers_Library/04_drafts_and_incomplete_papers/general_relativity_and_gravity/GR_PAPER_BRIEF.md>
[C6]: </home/williaml/Cella Framework/Papers_Library/04_drafts_and_incomplete_papers/geometric_fault_localization_and_decomposition/ARM_CLOSE_DRAFT.md>
[C7]: </home/williaml/Cella Framework/Papers_Library/02_theorems_and_lemmas/geometric_fault_localization_and_decomposition/F4_THEOREM.md>
[C8]: </home/williaml/Cella Framework/Papers_Library/01_completed_papers/local_curvature_and_black_hole_metrics/pfc_normal_forms.tex>
[C9]: </home/williaml/Cella Framework/Papers_Library/04_drafts_and_incomplete_papers/dbp_role_channel_and_orbit_geometry/DBP_UNIFIED_THEOREM_SPINE_DRAFT_v0_1.tex>
[C10]: </home/williaml/Cella Framework/Papers_Library/01_completed_papers/precision_flow_and_transfer_functions/transfer_function_exponent_family_v6.tex>
[C11]: </home/williaml/Cella Framework/Papers_Library/01_completed_papers/galois_horizon_and_kummer_covers/GENERIC_SYMMETRIC_MONODROMY_OF_WEIGHTED_MULTIQUADRATIC_SUMS_v1.0.md>
[C12]: </home/williaml/Cella Framework/Papers_Library/01_completed_papers/galois_horizon_and_kummer_covers/ALL_K_TWO_RADICAL_KUMMER_CLOSURE_v1.0.tex>
[C13]: </home/williaml/Cella Framework/Papers_Library/01_completed_papers/galois_horizon_and_kummer_covers/galois_horizon_cover_v1_0.tex>
[C14]: </home/williaml/Cella Framework/Papers_Library/05_expository_companions_and_research_maps/dbp_periods_landen_and_elliptic_structure/DBP_CURVATURE_PERIODS_OF_THE_DBP_QUADRIC_v1.0.md>
[C15]: </home/williaml/Cella Framework/Papers_Library/01_completed_papers/cella_residue_and_coupling_theory/coupling_field_theory_paper.md>
[E1]: https://arxiv.org/abs/0708.0659
[E2]: https://www.math.ucdavis.edu/~krener/1-25/10.IEEETAC77.pdf

# Mathematics-library review: corrections and stronger targets

5 September 2026. Addendum to this session’s corpus audit and proposed ROI update. This is a targeted source review, not a proof audit of every library item. Recommendations remain proposals; no corpus source, engine implementation or DAG record was edited.

## 1. What this directory adds

`/home/williaml/Cella Framework/01_MATHEMATICS_LIBRARY` is principally a topic navigation layer over the existing Papers_Library and campaign families. Following its links exposed 7,412 file paths resolving to 1,309 distinct files, including code, evidence and support material. Its 84 symlinks were unbroken. Those counts must not be interpreted as numbers of independent mathematical results.

The most consequential additional reading was:

- The native `Algebra/Groupoids/` reconciliation: README, canon, proofs and frontier.
- `CELLA_FRONTIER_RESEARCH_ROI_OUTLOOK_2026-07-17.md`, a different document from the previously requested `docs/ROI_OUTLOOK_2026-07-17.md`.
- The artificial-restriction reconciliation, finite-tower theorem, role-singularity valuation brief and rotating three-channel closure derivation.
- The actual selected-skeleton, AC-fold and CCE-5 implementations.

The Groupoids README labels its consolidation a draft awaiting acceptance. Physical presence does not make it accepted authority. Its reported categorical counterexamples were therefore checked directly against the implementation.

## 2. Corrections to the earlier assessment

### The common categorical implementation needs repair

The live `SelectedSkeletonMorphism` is a frozen dataclass whose equality includes its route string and digest. Composition embeds parenthesized route strings. Direct probes reproduce:

```text
reduce(0) == reduce(2)                          False
identity ; reduce(0) == reduce(0)               False
reduce(1+1) == reduce(1) ; reduce(1)            False
(a ; b) ; c == a ; (b ; c)                      False
```

The first is a failure of the intended parity reduction; the middle two include identity/functoriality failures; the last is an associativity failure under the implementation’s actual equality. The native integer-winding action remains a separate, valid construction. The abstract selected-quotient foundation is not refuted by this adapter defect.

**Immediate target:** define mathematical arrow equality and separate it from route provenance. For the AC fixture, a conventional action groupoid suffices: objects are the two sheets; arrows are `(source, parity)` with the target determined by the parity. Composition adds parity modulo two. The integer-winding reduction then has even-winding kernel on isotropy, and the expected functor laws can be proved. Retain the unreduced route as evidence attached to an arrow, rather than part of that arrow’s equality.

Do not erase winding globally: any observable sensitive to even winding cannot descend through this reduction. The stronger theorem classifies which observables and transport data survive the chosen quotient, and constructs a less destructive target when necessary.

CCE-5 independently rejects the empty route word. Its exact matrices remain useful, but a categorical encoding needs endpoint-specific identities and an explicit path-word equivalence. Equality of matrices must not silently become equality of paths: a representation may have a kernel.

### Finite-order role naturality is already established

The CCE-8 theorem explicitly covers every finite `N ≥ r ≥ 2` and every finite generator word:

\[
\tau_{N,r}(wf)=w\tau_{N,r}(f).
\]

Its triangular inverse-series argument supports this scope: a coefficient of degree `d` depends only on coefficients through degree `d`. Rechecking a few more finite orders is not a new frontier. Separate the next questions: more roles, extension across chart failure, global invariant separation, and analytic convergence. A formal inverse-limit construction and an analytic convergence theorem are different obligations.

### Two limitations were already caught within the library

The Groupoids reconciliation already rejects global completeness of the squared coupling carrier. Our earlier counterexample corroborates that internal finding; it is not a newly discovered corpus-wide omission. The full regular three-role RoleChSpec result has a different carrier and scope and must not be discarded with the squared carrier.

The same reconciliation corrects the rounding theorem: failure is contained in preimages of coarse midpoints; not every such point fails. With fine spacing 1, coarse spacing 2 and ties to even, `x = 3/4` rounds first to 1, then to 0; direct coarse rounding also gives 0. The frontier is the exact parity-dependent failure locus and its behaviour under composition.

### Four-charge rotating closure is also banked

The rotating source proves, relative to its stated static valuation and primitivity inputs,

\[
\operatorname{Gal}(\widetilde H_J/F(J))\cong C_2^2\wr S_5,
\qquad
\operatorname{Gal}(\widetilde L_J/F(J))\cong C_2^3\wr S_5.
\]

The live DAG marks this source proved/released, and its recorded digest matches the file. The source itself retains an older independent-counter-audit request; I have read the derivation, not rerun every upstream proof obligation.

Consequently, “rotating Kummer extension” must be narrowed as a proposed successor. Worthwhile remaining targets include exceptional parameter loci, broader canonically specified families, and R13b physical branch selection. The generic four-charge group is not an open target. At `J=0`, the source already identifies the returning square-class relation and augmented rank drop.

## 3. Strengthen the second outlook’s arithmetic programme

The Frontier ROI report’s parity/intersection and period-transport proposals are substantial research directions. Their initial construction needs to become more precise before comparing ranks.

### Start from square classes and the data valuations forget

For a declared field `K`, divisor set `D` and radicands `f_i`, define

\[
W=\langle[f_i]\rangle\subset K^*/K^{*2},
\qquad
v_D:W\longrightarrow\mathbf F_2^D.
\]

Classify the kernel and image, and determine their behaviour under the admitted base changes and monodromy. A full-column-rank valuation matrix proves independence on the chosen span; it does not establish that divisor parity classifies arbitrary square classes.

A useful general baseline is the Kummer exact sequence. On a scheme `U` with 2 invertible it gives

\[
0\to\mathcal O(U)^*/\mathcal O(U)^{*2}
\to H^1_{\mathrm{\acute et}}(U,\mu_2)
\to\operatorname{Pic}(U)[2]\to0.
\]

Thus units and line-bundle torsion require attention when moving from fields and valuations to a geometric family. This is established mathematics; the proposed contribution is the explicit Cella-family computation and comparison. [Stacks Project, Kummer theory](https://stacks.math.columbia.edu/tag/03PK).

Only then ask whether an independently constructed integral cycle lattice realizes the relevant module. Fix the lattice, comparison map and allowed subquotients before testing: “some subquotient” permits too much freedom to make rank matching informative. A mismatch establishes failure of that comparison, not automatically a canonical extra arithmetic summand.

### Distinguish a finite root cover from a family of curves

A finite separable root cover over its regular complex parameter locus has finite discrete fibres. Its sheet transport naturally acts on degree-zero data; those fibres do not supply nontrivial first-homology periods. This follows from the finite-cover structure, not a special defect of Cella. [Stacks Project, finite locally constant sheaves and finite étale covers](https://stacks.math.columbia.edu/tag/03RV).

The outlook’s proposed Gauss–Manin programme therefore needs an explicit additional construction if it intends first-homology periods: a family of curves or relative pairs, a differential and a cycle, together with a comparison to the root cover. Restricting the parameter base to one dimension does not alone change the fibres of the finite map into curves.

**Stronger target:** construct the smallest such family or relative pair naturally supplied by the existing equations. Prove which native period it represents and compare both monodromy representations on explicit generators. A differential equation satisfied by an algebraic root is not, on its own, the proposed bridge to the elliptic-period arm.

## 4. Give singularity continuation an independent first-principles target

The role-singularity brief is explicitly conjectural and already distinguishes carrier collapse from curvature valuation. It also warns that the reflection-fixed pole cannot be derived from its generic diagonal toy model. That is a useful starting point for growth.

Strengthen it to a classification of what local data determine an intrinsic response under declared changes of chart, frame and path parameter. Track cancellation and subleading unit jets; distinguish chart failure, observation ambiguity and metric rank loss.

For Seated Root the exact bridge remains

\[
N^2=\frac{\varepsilon}{(1+b^2)\varepsilon+d^2},
\quad \varepsilon=1-\gamma^2>0,\quad d=a-b\gamma.
\]

Its finite-endpoint approach regimes have already been classified in the earlier audit. The next mathematical question is which extra observations distinguish those regimes; the next physical question is which a declared evolution law reaches.

This programme need not wait for the parity/intersection programme. Its local geometric hypotheses can be investigated independently. Connecting its curvature response to a horizon discriminant is an additional theorem, not a prerequisite for doing the local mathematics.

Avoid promising one universal finite jet bound across arbitrary analytic germs. Arbitrarily late perturbations can agree to any prescribed finite order with a flat germ while changing its curvature. A stronger achievable target is a bound for a specified bounded family, or an adaptive procedure with clearly stated termination conditions.

## 5. Revised priority order and novelty target

1. **Repair the categorical target used by existing reductions.** Prove its laws, then reprove functoriality and characterize information loss. This is a bounded prerequisite, not a demand for a new foundational theory.
2. **Develop observable sufficiency across singular strata.** Continue the extractability and cancellation-sensitive valuation work together, using Seated Root and an independent Cella fixture. Couple it early to operational metric selection and admissible dynamics.
3. **Compute arithmetic information loss before proposing a cycle identification.** Extend the banked static/rotating results through exceptional fibres and an explicit valuation/cohomology comparison.
4. **Construct the actual family for the period bridge.** Then derive its differential system, integral transport and compatibility with the existing evaluator.
5. **Generalize the role calculus where it is genuinely unfinished.** Prioritize additional roles, global separation, boundary behaviour and analytic control over already-established finite-order naturality.

The most promising original contribution is an explicit theorem about which observations, reductions and transports remain compatible through degeneration, with minimal additional data when compatibility fails. Ordinary quotient categories, Kummer theory, local systems and valuation methods provide strong starting frameworks. Novel mathematics should emerge from a missing comparison, classification or extension theorem, rather than from replacing those names in advance.

## 6. Evidence and source locators

The replayable [implementation probes](mathematics-library-checks.py) produced the accompanying [results](mathematics-library-checks.json). They establish the reported counterexamples; they do not rerun the corpus’s historical assertion totals.

Principal sources, relative to `/home/williaml/Cella Framework/`:

```text
01_MATHEMATICS_LIBRARY/README.md
01_MATHEMATICS_LIBRARY/Algebra/Groupoids/{README,01_GROUPOIDS_CANON,02_PROOFS_AND_CERTIFICATION,03_FRONTIER}.md
Papers_Library/02_theorems_and_lemmas/cella_residue_and_coupling_theory/CCE_8_FINITE_TOWER_NATURALITY_THEOREM_v1.0.md
Papers_Library/04_drafts_and_incomplete_papers/dbp_role_channel_and_orbit_geometry/LEAD2_Role_Singularity_Valuation_Brief.md
Papers_Library/05_expository_companions_and_research_maps/cross_program_and_unclassified/CELLA_FRONTIER_RESEARCH_ROI_OUTLOOK_2026-07-17.md
Papers_Library/05_expository_companions_and_research_maps/dbp_role_channel_and_orbit_geometry/DBP_CCE_ARTIFICIAL_RESTRICTION_LIVE_RECONCILIATION_v1.0.md
Papers_Library/05_expository_companions_and_research_maps/galois_horizon_and_kummer_covers/ROTATING_THREE_CHANNEL_WREATH_CLOSURE_2026-07-10.md
engine/src/cella/continuation/{selected_skeleton,r3_ac_fold,cce5}.py
```

DAG source `DBP:doc:rotating_three_channel` matched SHA-256 `0a9d68ff15899a7143f49bd6e53fe46506a5fb92fe119372193b5465aa51d7e0` at graph revision `575f1e672b44c5303a308fa04fef07fa46d61f1f186a65fad1396507d850871b`. A text query for `selected_skeleton` returned no nodes; that limited query does not prove the reconciliation is wholly absent from the graph.

# Contents

Three results. Two corpora, no shared dependencies.

## A. Gram submersion curvature  — for the seated-root paper

**`GRAM_SUBMERSION_CURVATURE.md`** — the curvature deliverable.
**`INERTIA_NODE_DETECTOR.md`** — the boundary-classification deliverable (det I = 0 iff rank-one node).
Both are supported by the scripts below; `verify_abstract.py` (8/8) checks the detector claims.

Closed-form curvature of pi:(S^2)^3 -> elliptope, exponent law with proof,
C(b) closed form. Answers the KIN-2 curvature-tier debt.

| script | role | status |
|---|---|---|
| `prove_exponent.py` | kernel proof of exponents m=1, m=2; exact-rational pin | primary |
| `cb_final.py` | on-shell reduction to Gram entries | primary |
| `cb_sym5.py` | symmetric closed form (needs `Pred.srepr`) | primary |
| `cb_exact.py` | exact C^2 via independent 9x9 recipe | cross-check |
| `cb_wrap.py` | trine value 16/27, positivity, det-inertia identity (needs `Cb_symmetric.srepr`) | cross-check |
| `closed_form.py` | verifies F = 2 I^-1 sum(X_i x Y_i) vs holonomy | cross-check |
| `cb_symbolic2.py` | component-level assembly | intermediate |
| `isolate_exponent.py` | shared helpers; imported by several scripts | dependency |
| `Pred.srepr`, `Cb_symmetric.srepr` | cached exact polynomials | dependency |

Superseded by the closed form, kept only as the independent numerical route:
`holonomy_gram.py`, `compare_holonomies.py`, `node_test.py`, `pin2.py`.
These use finite-loop transport and need h <= 0.05*Delta. The closed form has
no h and supersedes them for all use.

## B. Groupoid cohomology — for the SQG / DBP corpus

**`GROUPOID_COHOMOLOGY_MATHEMATICS.md`** — the deliverable.

| script | role | result |
|---|---|---|
| `verify_transport_obstruction.py` | gauge quotient n=3,4,5 | 28/28 |
| `adjudicate_row6.py` | wall obstruction is pointwise; H^1(ker Sigma, C_2)=0 | 14/14 |
| `check_anchor_C_against_cce5.py` | CCE-5 arithmetic; transport group is Z | 21/21 |
| `repair_cpv.py` | CPV as conjugation-invariant lift | 12/12 |

## C. Scope-flagged

`wick_channels.py` (11/11) — computes on the DBP curvature-channels
(kappa_c, kappa_int, kappa_s), NOT the seated-root seats (c, hbar, G).
Per PAPER §10 these are unrelated referents sharing the word "channel".
Valid as DBP-side mathematics; does not apply to the seated-root paper.

## Verification

All scripts run standalone from this directory (python3 <file>).
Exact arithmetic throughout except the four superseded numerical scripts.

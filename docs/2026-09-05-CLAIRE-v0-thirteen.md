# A thirteen-block reflection word with the gain leakage column removed (v = 0)

Claire, 2026-09-05 evening.  Runner `claire_v0_thirteen.py` (30 seeded starts, ~6 min), receipt `docs/claire-v0-thirteen-checks.json`.
Judge: Will's own evaluators (`reflection_loop_detuning_order`, `_reference_echo`, `_compression`); the design is the input.
No existing document or runner is edited.  **Tier: DERIVED | Will's evaluators as stated.**  No interval certificate: the centre is a
double-precision root with residual 5.6e-12.  KILL: any printed check failing on re-run.

## 1. The idea, in one equation

CM-3 eq. (10): the composite word's eps^2 bright column is l_n v with

    v = sum_j sigma_j (cos beta_j, sin beta_j).

The eleven-block has v = (0.6887, 0); that is the whole of its gain infidelity, (15/2) b_n^2 theta^2 |v|^2 eps^4 = 8.46 eps^4,
and it is also the whole of the raised wrapper's error at eps = 1e-3, Delta/a = 1e-4 (8.41e-12 joint vs 8.38e-12 gain-only).
So impose **v = 0** as one more real equation (v_y = 0 is automatic for antisymmetric angles).  The eleven-block's system is
square (8 equations, 8 unknowns); v = 0 needs the thirteen's two extra blocks, where it competes with Will's second-order
logical detuning conditions (`logical_second` = 0).  In his 13-unknown parameterisation that is 8 + 3 + 1 = 12 residuals
(one logical component redundant by time-reflection symmetry) with two degrees of freedom left for total stretch W.

## 2. What was found

Two words, in order of discovery:

| word | conditions | W | gain (eps = 1e-3) | gain order | A_Delta | six generators | joint (1e-3, 1e-4) | duration (1/a) | exposure (1/a) |
|---|---|---:|---:|---:|---:|---|---:|---:|---:|
| v0-only thirteen | first-order (8) + v = 0 | 15.511 | 2.2e-20 | ~8 | 1.13e7 | cancelled (1e-13) | 1.13e-9 | 895 | 50.6 |
| **v0 + logical_second thirteen** | first-order + second-order logical Delta + v = 0 | **17.954** | **1.7e-20** | **7.95** | **1.85e4** | cancelled (1e-11) | **1.89e-12** | **1036** | **58.6** |

Reference points (Will, same errors): eleven-block 1.53e-10 (686/a, 38.7/a); his thirteen 3.74e-10 (948/a, 53.6/a, A_Delta 4391);
nested echo 4.73e-11 (1386/a, 78.3/a); raised wrapper 8.41e-12 (2341/a, 132/a).

So the second word is: 200x the error improvement of his thirteen for 9% more duration and exposure; 25x better than the
nested echo while 25% shorter and 25% less exposed; 4.4x better than the raised wrapper at 2.3x less duration and exposure.
One word, no wrapper.  Its residual error is entirely detuning: the gain term has left the budget (eps^8, 1e-20).

**Why eps^8 and not eps^6.** v = 0 removes the eps^2 leakage amplitude, so the next candidate is eps^3; measured order 7.95
says the eps^3 column vanishes too, presumably by the same palindromic symmetry that makes v_y automatic.  Not proved here.

## 3. What it says about his thirteen

`logical_second` = 0 does not fix A_Delta by itself: both words above satisfy it, with A_Delta = 1.1e7 and 1.9e4, against
his 4391.  His coefficient therefore also relies on a small second-order Delta **leakage** column (the bright block of the
Delta^2 response), which no stated condition constrains -- the solver happened to land where it is small.  Killing that column
explicitly is two more real conditions, which is a **fifteen-block** target: v = 0, logical_second = 0, bright Delta^2 column = 0,
with the slack spent on W.  Expected outcome: gain eps^8, A_Delta at or below 4391, at W near 18.

## 4. Centre (double precision; his parameterisation [b_1..b_6, w_1..w_6, w_centre])

See `docs/claire-v0-thirteen-checks.json`, key `centre`.  Three admissible solutions were found in 30 starts (W = 17.954,
19.215, 22.142); the residual of the best is 5.57e-12.  An interval-Newton certificate in the style of CM-2 is owed before
this is more than a numerical witness; the equations are his, so `reflection_loop_detuning_certificate.refine/certify` should
take the centre directly.

## 5. What is not claimed

Not certified; not compressed (W = 17.95 is the best of 30 starts, not a minimum); not evaluated at n = 2, 3 (the n-lever from
CM-4 should apply unchanged, since v = 0 is n-independent and A_Delta scales with T^4); temporal response (filter onsets,
linear drift) not computed for this word -- CM-5's argument applies to any word satisfying his first-order set, so the
omega^4 gain onset and omega^2 reference onset are expected, but they are not checked here.

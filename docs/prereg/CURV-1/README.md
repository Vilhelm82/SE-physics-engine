# CURV-1 — does the model own a rod-free curvature reading, and what does it do at the centre?
Opened 2026-09-02 by `thm_h2_d1.py` (commit `b4ac768`). Specs written 2026-09-02, same day.
Read `PROJECT-BRIEF.md` and `HANDOFF-2026-08-30.md` first.

## The question
THM-H parts 1 and 2 (`thm_h.py` 25/25 at `0ba60a4`, `thm_h2_d1.py` 30/30 at `b4ac768`) showed:
along the pinned family the seat's presented Gram is G'(r) = G + f(r) k k^T with
f = r_s/(r - r_s) and k_i = a_i . K, on both sheets; det G' = Delta r/(r - r_s); and no
rod-free click-invariant diverges at r = 0. Every divergence found there has a rod
(a presented length, r, or dr) in its denominator. That leaves one fork open:

- "no singularity" is a theorem of the model (the centre is a place where one rod
  reading vanishes and nothing else happens), or
- the model is too coarse to see one, because it has never computed a curvature.

CURV-1 closes the fork by building the model's own curvature readings — quantities
built from angles and rotors only — and evaluating them along the family, at the
centre and at the horizon.

## Frozen inputs (cite by check number; do not re-derive)
| Input | Content | Status |
|---|---|---|
| BARE-1 | rotor A = cosh(l/2) + sinh(l/2) K, two-sided sandwich X -> A X A | proved (`thm_e.py`) |
| H-1/H-2 | presented vector part v_i = a_i + (cosh l - 1) k_i K | proved (`thm_h.py`) |
| H-5/H-7 | G' = G + sinh^2(l) k k^T, det G' = Delta cosh^2(l) | proved (`thm_h.py`) |
| KIN-2a | pinning tanh(l) = sqrt(r_s/r) | declared |
| CONT-1 | interior sheet l = mu + i pi/2, r = r_s tanh^2(mu) | declared |
| D-1..D-4 | one coefficient f(r) = r_s/(r - r_s) on both sheets | proved (`thm_h2_d1.py`) |
| D-14..D-17 | two pinned seats: R(eps) x A(l2 - l1, K1); i pi/2 cancels inside; -> R at the centre; crossing = fixed quarter turn c = (1 + iK)/sqrt2 | proved (`thm_h2_d1.py`) |
| RULING-2 | K_C ~ a1+a2+a3 is the world layer (the hole); K_A = root is the seat layer (the infalling seat) | declared (Will, 09-02) |
| area-tangent | tan(Omega/2) = V/(1 + g12 + g13 + g23) for three unit vectors; B = S + iV, 2|arg B| = spherical excess | banked 08-30 (PRED-1, E-8-X) |

## The three paths
1. `PATH1-SPEC.md` — the presented seat triangle's spherical excess, profiled in r on
   both sheets. Uses only the presented Gram. Cheapest; first.
2. `PATH2-SPEC.md` — loop holonomy at the centre from part 2's relative rotors: a
   small loop of seats around r = 0, the product of their relative rotors, the limit
   as the loop shrinks. Written after Path 1's results are in.
3. Path 3 — held until 1 and 2 report. Not with a chosen constraint set: every
   candidate set that already carries a freeze date is run, and agreement or
   disagreement between them is the result.

Each path runs both K choices (A: K = root; C: K ~ a1+a2+a3).

## Classification (the same rule for every path, every endpoint)
At each endpoint (centre r -> 0+; horizon r -> r_s+ and r -> r_s-) the reading is
one of

| class | meaning | record |
|---|---|---|
| identity | 0 mod 2 pi | the sign of the half-angle (the spin lift, +1 or -1 in SU(2)) |
| finite-nontrivial | finite, not 0 mod 2 pi | the value |
| unbounded | Im Omega -> +-oo (e^{i Omega} -> 0 or oo), or no limit | the rate |
| indeterminate | a 0/0 whose limit along the family does not exist | the two expressions |

A path must be able to return every one of these. A check that can only print
"finite" is not a check.

## Hand-off
- The runner receives the spec file(s) and the repo. `PREDICTIONS.md` in this folder
  is for the designers and the comparison stage, not for the runner.
- Will's earlier, independent work on the same question (a different corpus) enters
  only at the comparison stage, after the model's own answers are in.
- Verification is re-running: the runner's suite must be green, self-contained
  (sympy/mpmath only), and committed together with its log and a results file.

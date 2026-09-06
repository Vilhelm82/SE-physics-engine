# RH-2P -- the eps^3 law is geometric: precision-scaling separation (Cella U-0784) on DB_stack

**[numerical, 50-digit; derived under Cella Theorem 2]** The tight-frame deviation of the stacked leak obeys
DB_stack = a0 eps^3 (1 - kappa eps + ...) with **a0 = 7.1409907** (second-order Richardson at 50 dps, converged to
seven digits; NOT sqrt(51) = 7.14143 -- the difference is 4e-4, far outside the convergence). The exponent is
exactly three: local log-log slope 2.9702 -> 2.9992 from below across eps = 2e-2 .. 3.1e-4, with a two-term
correction of exponent rho = 1.048 (U-0799). The amplitude is path-invariant across four arithmetic routes to
16+ digits at 50 dps, and eps-independent up to that O(eps) drift. The cubic law is physics, not the integrator.

**[provenance]** Run 2026-09-06 by Claire from HEAD `95ab648` (runner with predictions P1-P4 in its docstring
committed before the log). `rh2_precision.py`, receipt `rh2_precision-checks.json`, log `rh2_precision.log`.
20/20 numerical checks. The float64 closed form reproduces `rh1_common.trine_loop` to 0.0 (bit-identical
algebra). Precisions: complex64 (u = 2^-24), complex128 (2^-53), mpmath 30 dps (2.0e-31), 50 dps (2.7e-51).

**[Cella applied]** U-0784 (thm:precision-separation): C_{p,k} = a + u_p b_k, a geometric and path-invariant,
b_k conditioning and path-dependent, identifiable from >= 2 precisions. U-0786 (log-log transfer slope).
U-0799 (two-term finite-window bias). U-0776 (exponent-only non-identifiability -- the reason one precision
cannot settle this). Source: `Papers_Library/01_completed_papers/precision_flow_and_transfer_functions/
transfer_function_exponent_family_v6.tex`.

## Theorem 2 on the gate

| eps | a = C(mp50) | |a30-a50|/a | b_k1 (svd) | b_k2 (stable) | b_k3 (tr/det) | b_k4 (I-K^dagK) | f32 pred/meas k2 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2e-2   | 6.867658 | 7e-29 | -1.4e5 | -1.3e4 | 4.4e8   | -9.1e3 | 1.000 |
| 1e-2   | 7.011092 | 7e-28 | -1.7e6 | -7.3e4 | -4.2e10 | -2.4e6 | 0.999 |
| 5e-3   | 7.077808 | 2e-27 |  1.5e7 | -9.5e4 | -6.0e11 |  3.5e7 | 0.998 |
| 2.5e-3 | 7.109850 | 8e-28 |  7.9e6 | -1.4e6 | -1.5e14 | -1.8e8 | 0.886 |
| 1.25e-3| 7.125534 | 3e-26 |  3.9e8 | -9.4e5 | -5.7e15 | -1.9e9 | 0.356 |
| 6.25e-4| 7.133291 | 1e-25 |  2.6e9 | -4.7e6 |  4.0e17 |  4.5e9 | 0.034 |
| 3.1e-4 | 7.137148 | 6e-25 |  1.7e10| -1.9e7 | -6.4e16 |  3.6e14| 0.002 |

**[numerical]** a is identified: the two high precisions agree to 1e-25..1e-29, and all four paths coincide at
50 dps to 1e-42..1e-47. b_k spans ten orders across paths at the same eps -- the theorem's content, on this gate.
The stable 2x2 formula (k2) and the SVD (k1) hold in float64 down to eps = 3.1e-4 (DB = 2.2e-10); the tr/det
route (k3) returns 0.000 there and is 7x wrong at 6.25e-4; I - K^dag K (k4) degrades only to 0.6% at 3.1e-4.
float32 is predicted by a + u_32 b_k to 0.2% for k2 at eps >= 5e-3 and fails below 2.5e-3; the badly conditioned
paths are outside the linear regime at float32 everywhere. At eps = 6.25e-4 the float32 SVD returns DB = 0
exactly: the two singular values collide on the float32 grid -- the lattice-grain collapse (BACL, U-0766).

## Limits, by Richardson at 50 dps

| quantity | first-order sequence -> | second-order -> | value |
|---|---|---|---|
| a0 = lim DB/eps^3 | 7.15452, 7.14452, 7.14189, 7.14122, 7.14105, 7.14100 | 7.14119, 7.14101, 7.14099, 7.140991, 7.140991 | **7.1409907** |
| c0 = lim c/eps^2   | 16.3635, 16.2904, 16.2712, 16.2663, 16.2650, 16.2647 | 16.26600, 16.26476, 16.26463, 16.26461, 16.26461 | **16.264607** |

Against FD-5's endpoint-flag form (15/2) b_1^2 eps^2 for the original trine, c0 = 16.264607 gives b_1^2 = 2.16861.
[Comparison-stage reading; not re-derived here.]

**[derived + numerical]** Accumulator floor (RH-2) = c DB^2/6 = (c0 a0^2/6) eps^8 = **138.232 eps^8** asymptotically.
Local slope of the floor: 7.899 -> 7.9976 across the grid. At eps = 1e-2: 138.232e-16 + conditional 2.495e-15 =
1.63e-14 predicted against 1.553e-14 measured in RH-2 (5%; the O(c^2) remainder). So the unconditional gain error
of the accumulated gate is 138 eps^8 plus the survivor's own eighth-order conditional error -- same order, one
explicit coefficient, no discards.

## Predictions (docstring, pre-run) reconciled
P1 (a finite, nonzero, path-invariant, eps-independent): **hit**, with the eps-independence holding only up to the
O(eps) two-term drift, which U-0799 predicts and rho = 1.05 confirms.
P2 (b_k3, b_k4 >> b_k2; k4 fails outright by eps <= 1.25e-3): **half**. k3 fails as predicted; k4 is far more robust
than predicted (0.6% at the smallest eps, exact above it). The stable 2x2 formula is insensitive to the common-mode
part of the roundoff in forming I - K^dag K; my conditioning estimate ignored that.
P3 (slope -> 3.000, rho ~ 1): **hit**. Approach is from below; rho = 1.048.
P4 (float32 predicted within 2x only at eps >= 1e-2): **hit and better** -- 0.2% for the well-conditioned paths down
to 5e-3; the linear regime extends lower than I allowed for.

## Open
- **What is a0 = 7.1409907?** Not sqrt(51). A symbolic derivation from the three effective directions
  (a, 2a - b + pi, a - b) at (pi/3, 2pi/3, pi/3) with gain-perturbed reflections is the natural route; the cubic order
  says the tight-frame condition is exact at O(eps) AND O(eps^2) and first breaks at O(eps^3), which is a statement
  about the reflections, not the loops. Not attempted here.
- Five-loop word not run through this test (RH-2 measured the same eps^3 there at float64 only).
- rho = 1.048 rather than 1.000 over this grid: either a genuine non-integer correction or the next term's shadow;
  a finer grid at 50 dps would tell. Not pursued.

**Lineage, corrected:** RH-1's "unsearched" is now searched in the only corpus that mattered. The exact-cancellation
information loss (U-1558) is the flat zero of RH-1c; the early/late reassociation transition (U-0873, conjectured)
is RH-1b vs RH-2; the BACL invariant (U-0766) is the representability condition the fold-back needed. Cross-corpus
edges are candidates for the DAG when submissions reopen; none submitted.

Replay: `OPENBLAS_NUM_THREADS=1 python3 rh2_precision.py --json docs/rh2_precision-checks.json | tee docs/rh2_precision.log`

**Closed form found 2026-09-06 (RH-2A, `docs/2026-09-06-RH-2A-a0-closed-form.md`):** a0 = 3375 sqrt(5) pi^3/32768 exactly;
with x = (15 pi/32) eps, DB = sqrt(5) x^3 (1 - 55 eps/32 + ...), c = (15/2) x^2, floor = (25/4) x^8 = 138.2322 eps^8.
The "rho = 1.048" two-term fit is the -55/32 linear correction plus the next term's shadow; the exponent is exactly 1.

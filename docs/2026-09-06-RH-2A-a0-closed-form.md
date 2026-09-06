# RH-2A -- closed form for a0: the tight-frame deviation is sqrt(5) x^3, x = (15 pi/32) eps

**[proved | exact series, two independent routes, Cella-pinned]**

    a0 = 3375 sqrt(5) pi^3 / 32768 = 7.14099061748716977...
    c0 = 3375 pi^2 / 2048          = 16.2646068621467702...
    DB_stack = a0 eps^3 (1 - (55/32) eps + O(eps^2))
    floor    = c0 a0^2 / 6 = 64072265625 pi^8 / 4398046511104 = 138.232207895 eps^8

**[provenance]** `experiments/2026-09-06/rh2_a0_closed_form.py`, 12/12, receipt `rh2_a0-checks.json`, log `rh2_a0.log`. Route A: exact series
of the closed-form loop (same algebra as `rh1_common.trine_loop`) through the 2x2 frame reduction. Route B: 70-digit
third-order Richardson on `rh2_precision`'s closed form. Agreement 3.1e-16 (14 decimal digits certified by
`cella_two_route_compare`; constant pinned by `cella_arith_constant_pin`, digest `566bbaea565eabb4`). RH-2P's
Richardson value 7.1409907 and its measured a(3.1e-4) = 7.1371480 are reproduced (series predicts 7.1371551,
the residual being the O(eps^2) term). Derived 2026-09-06 by Claire.

## The reduction that made it a one-page derivation

At exact return every arc is a 2pi spin-1 rotation, so at eps = 0 each loop is the frame flip R = frame(2, pi/2):
a reflection I - 2 u u^T on P, a swap on Q, no mixing. At eps != 0 the loop is unitary on the 3-frame (p, q, D) and
fixes the spectator exactly, so on P it is

    S(alpha, eps) = I - (1 - sigma(eps)) u_alpha u_alpha^T,       E(alpha, eps) = ell(eps) u_alpha^T,

with ONE complex number sigma = <D|L|D>, ONE bright vector ell, and |sigma|^2 + |ell|^2 = 1 by unitarity of the D
column. (This is RH-1a's rank-one leak, now with its coefficient named.) Hence

    F^dag F = |ell|^2 sum_k conj(w_k) w_k^T,    w_1 = u_1,  w_2 = S_1 u_2,  w_3 = S_1 S_2 u_3,

and DB_stack is a function of sigma and the three angles alone. Everything else is series arithmetic.

## What sigma(eps) is, and why the exponent is three

    sigma = -1 + (1125 pi^2/2048) eps^2 + (1125 pi^2/32768)(-29 + 6 sqrt(15) i pi) eps^3 + O(eps^4)

The eps^2 term is REAL. A real shift of sigma turns the reflection into a partial reflection with the same axis, and
the augmented-triad frame stays tight under that: the traceless part of F^dag F vanishes at every order that only
sees Re sigma. The first IMAGINARY part -- the phase of sigma -- appears at eps^3, and only the phase tilts the
effective directions off the 60-degree frame. So t^2 + |w|^2 starts at eps^10, DB at eps^3, and a0 is proportional
to Im sigma_3 = 3375 sqrt(15) pi^3 / 16384. The sqrt(15) combines with the frame algebra's sqrt(3) to leave sqrt(5).

## The physical parameter

Let dphi = (15 pi/8) eps be the first-order over-rotation of each arc beyond 2pi, and q0 = 1/4 the D-to-bright mixing
ratio at eps = 0 (q = s/omega = (1/sqrt15)/(4/sqrt15)). Put x := q0 dphi = (15 pi/32) eps. Then

    c     = (15/2) x^2          [erasure; (3/2) from the tight-frame sum, 5 = |(2,-1)|^2 from the bright direction]
    DB    = sqrt(5) x^3         [tight-frame deviation]
    floor = (25/4) x^8          [accumulator unconditional gain error, RH-2]

and 15 pi/32 is exactly FD-5's b_1 (its "(15/2) b_n^2 eps^2" endpoint-flag law for the original trine): b_1 = q0 dphi/eps.
The bright direction (2, -1) is the arc pattern (G1, G2, G1) coupling D to p twice and to q once; sqrt(5) is its norm.

**[derived, n = 1 only; conjectured for general n]** With tau_n = (pi/2) sqrt(16 n^2 - 1): q0 = 1/(4n) and
dphi = 2 pi n eps (16 n^2 - 1)/(16 n^2), so b_n = pi (16 n^2 - 1)/(32 n^2). Whether c = (15/2) b_n^2 eps^2,
DB = sqrt(5) b_n^3 eps^3 and floor = (25/4) b_n^8 eps^8 hold at n = 2, 3 with the same rational prefactors is NOT
derived here: the bright-direction ratio 2:1 is n-independent, but the frame algebra was done at n = 1 only.
One rerun of this script with the n = 2 primitives decides it.

## Tiers
- a0, c0, kappa1, the floor coefficient, sigma to eps^6: **proved** (exact series; 12 symbolic identities checked).
- Route agreement 14 digits: **numerical**, certified by the router's two-route referee; assumes route independence,
  which holds here (series arithmetic vs. Richardson on a floating closed form share only the loop definition).
- The general-n formula for b_n: **derived** for the primitives; its consequences for c, DB, floor: **conjecture**.
- Nothing here touches the absorbing chi > 0 interface, detuning, or the five-loop word (whose own a0 is unknown
  and, by RH-2, of the same order: DB = 4.2e-6 at 1% against the trine's 7.0e-6).

Replay: `python3 experiments/2026-09-06/rh2_a0_closed_form.py | tee docs/rh2_a0.log`


---

## Addendum 2026-09-06: the five-loop word (RH-2P/5), `experiments/2026-09-06/rh2p_five.py`, 21/21

**[proved | same reduction, five angles, sigma -> conj(sigma) on the two inverse loops]**

    a0(five) = 2025 sqrt(5) pi^3 / 32768 = 4.28459437049230186 = (3/5) a0(trine) = (3 sqrt5 / 5) x^3
    c0(five) = 5625 pi^2 / 2048          = 27.1076781035779    = (25/2) x^2
    kappa1   = -55/32  (identical to the trine)
    floor    = (15/4) x^8 = 38443359375 pi^8 / 4398046511104 = 82.9393 eps^8

**[numerical]** Four precisions x four paths (Cella U-0784 again): a identified to 1e-25..1e-29, paths agree to
1e-42 at 50 dps, slope 2.9680 -> 2.9992 from below, floor slope -> 7.9976. Third-order Richardson at 70 dps:
a0 = 4.28459437049229802 (rel 9.0e-16 to the closed form), c0 rel 3.1e-15, kappa1 numerical -1.71899 (0.014%).
RH-2's float64 DB_stack(1e-2) = 4.204421e-6 reproduced to 2.8e-11. RH-2's measured self-calibrated
accumulator error 1.03e-14 at 1% against 82.94e-16 + conditional 2.5e-15 = 1.08e-14 (5%).
Conditioning as for the trine: tr/det (k3) returns 0 in float64 at eps <= 6.25e-4; I - K^dag K (k4) is 20% off
at 3.1e-4; the stable formula and the SVD hold to the bottom of the grid.

## Why 3/5, and what the five-loop word actually is

**[proved | leading order]** The leading coefficient of t^2 + |w|^2 is IDENTICAL for the two words:
sqrt(D2_10) = 11390625 sqrt(5) pi^5 / 67108864 for both. The absolute anisotropy of the stacked leak is the same;
the five-loop word only divides it by more erasure. Splitting the word (sandbox check, same series):

- loops 1-3 alone, angles (pi/6, 5pi/6, 7pi/6): D2_10 and c0 equal to the trine's exactly. Their effective
  directions are 30, 90, 150 degrees -- the augmented triad rotated by 30 degrees. Codex's word contains a
  tight-frame trine.
- loops 4-5 alone, angles (3pi/4, pi/4), both inverse: the traceless part of F^dag F vanishes at EVERY computed
  order (through eps^13). Their effective directions are 135 and 45 degrees -- orthogonal -- and for u_2 perp u_1
  the partial reflection S_1 fixes u_2 exactly, so the pair is a two-vector tight frame for all eps, not just at
  leading order.

So: **five-loop = (rotated tight trine) + (exactly tight orthogonal pair)**. The pair supplies the detuning
cancellation (FD-2) and (2/3 of the trine's) extra erasure, with no leading anisotropy. Hence
c0 scales by 5/3, a0 by 3/5, and the floor c0 a0^2/6 = |T|^2/(6 c) by 3/5.

**[derived, and worth saying out loud]** The accumulator floor is |T|^2/(6c): the SQUARE of the absolute
anisotropy divided by the total erasure. Adding loops whose leak is isotropic RAISES the erasure and LOWERS the
unconditional floor. More (isotropic) leakage is better for the accumulator, at the price of a register that must
capture more. That is not intuitive and it is exact at leading order.

**[not checked]** Equality of the full D2 series beyond the leading coefficient; general n; the absorbing interface.

Replay: `OPENBLAS_NUM_THREADS=1 python3 experiments/2026-09-06/rh2p_five.py | tee docs/rh2p_five.log`

# The split loop: return the data, flag the error

Claire, 2026-09-05 late.  Runner `experiments/2026-09-05/split1_split_loop.py`, receipt `docs/split1-checks.json`, log `split1_split_loop.log`.
No existing document or runner is edited.  **Tier: numerical witnesses on a declared model** (bright-selective absorbing
loss, static gain and detuning, ideal endpoint projection).  Kill conditions stated at the end.

## 1. The ask

Will, 2026-09-05: "find the loop that only returns the data and the error is split during the loop."  I.e. erasure
conversion: the gate returns the encoded plane exactly and every error goes somewhere flaggable, so the conditional
fidelity is set by what survives in the code, not by how much leaked.  The composite (CP) is the opposite philosophy --
it cancels the leak so no flag is needed.

## 2. Part A -- the obvious version fails, and RF-7 says why

Adiabatic dark loop (C^2-smoothed arcs) with bright-selective loss kappa ON during the loop, symmetrised by the trine
(pi/3, 0, -pi/3).  Premise: nominal bright population is zero, so anything bright is error and loss both flags it and
Zeno-suppresses it.

Result (a tau = 20, three loops): erasure 8.5% at kappa/a = 0.1, 46% at 1, 26% at 10; spread across states ~ the mean
(state-dependent); conditional infidelity 5e-4, 2.6e-2, 5.6e-3.  At a tau = 40, kappa/a = 1: erasure 27%, cond 6e-3.

The premise is wrong.  Adiabatic following IS an O(v/a) bright admixture, integrated over a loop of length 9 tau it is
order-unity exposure (~22/tau per trine), and loss attacks the transport itself -- Zeno freezing.  RF-7's theorem is the
reason: no real symmetric Hamiltonian transports a nonconstant real projector, so in this family transport always goes
through the bright sector.  **Bright loss during the loop cannot split error from data.**

## 3. Part B -- the split is at the return, and it is a point

At the exact return T_n (RF-2) the nominal bright population is exactly zero.  So any bright amplitude at the endpoint is
error by definition, and dumping it there -- loss or projection at the hold, where the encoded plane is stationary and
H P = 0 -- splits without back-action on the nominal state.

**Single exact-return loop + endpoint dump (18/a):** gain eps = 1%: heralded erasure 5.4e-4, conditional infidelity
**4.9e-8** (vs 5.3e-4 unheralded; vs the composite's 2.5e-7 at 3x the duration).  eps = 0.1%: 5.5e-6 and 5.1e-12.
But state-dependent: only d leaks, r never does; spread/mean = 1.9.

**Symmetrisation, and the angle correction.**  With dumps between loops the reflections reorient the effective active
directions: for loop angles (a, b, c) the effective directions are a, 2a - b + pi, a - b, and the gate needs c - b + a = 0.
The composite's (pi/3, 0, -pi/3) puts two effective directions on one line and leaves the spread at 1.1.  The tight-frame
solution -- three effective directions at 60 degree spacing, sum u u^T = (3/2) I -- is

    (a, b, c) = (pi/3, 2 pi/3, pi/3).

**Tight-frame trine + dumps (55/a):**

| errors | heralded erasure | spread / mean | conditional infidelity |
|---|---:|---:|---:|
| eps = 1% | 1.6e-3 | 1.4e-5 | **2.4e-15** |
| eps = 0.1% | 1.6e-5 | 1.4e-8 | ~1e-18 (numerical zero) |
| Delta/a = 1e-4 | 1.8e-9 | 1.4 | 2.4e-8 |
| Delta/a = 1e-3 | 1.8e-7 | 1.4 | 2.2e-6 |
| eps = 1%, Delta/a = 1e-4 | 1.6e-3 | 2.5e-3 | 2.4e-8 |

Gain error is converted entirely into a data-blind flag and the surviving gate is diag(-1, 1) to machine precision.  At
55/a this is the composite's duration (2.5e-7 unheralded) and 25x shorter than the reference echo (4.7e-11 at eps = 1e-3,
Delta = 1e-4, 1386/a).

## 4. What is not split, and what that changes

Detuning.  A site or reference offset is a logical rotation; nothing goes bright, so nothing is flagged, and the
conditional infidelity is first order in Delta (2.4e-8 at 1e-4, 2.2e-6 at 1e-3).  It is now the ONLY coherent error left.

That reframes the reference-echo programme: the echo no longer needs to cancel gain.  A word designed for detuning alone,
with a dump at every block boundary -- every block ends at H = 0 with nominal bright population zero, so the dump condition
holds throughout -- should be much shorter than the eleven's 686/a.  That is the next target, and it is Will's to take or
leave.

## 5. Tiers and kills

- Part A: DERIVED as a negative result on the declared loss model; the mechanism (RF-7) is Will's theorem.
- Part B: numerical witness; the tight-frame angles are derived (effective-direction algebra + tight-frame condition), the
  conditional numbers are ODE + projection at rtol 1e-11, 400 Haar states.  The dump is IDEAL (projection).  A physical dump
  (fast bright loss at the hold) needs its own finite-duration analysis; the hold can be as long as wanted since H P = 0.
- KILL: a state-dependent erasure surviving at first order for the tight-frame angles (spread/mean is 1e-5 at 1%, 1e-8 at
  0.1%: it does not); a conditional gain floor above 1e-12 (it is 1e-15 at 1%); the nominal block not equal to diag(-1, 1).
- NOT done: n = 2, 3; temporal response; a detuning-only word with dumps; the physical dump's finite time; the mixed
  gain-detuning term beyond the two rows above.

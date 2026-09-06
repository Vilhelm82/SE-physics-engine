# RH-4R -- FD's six readings with the accumulator on the same axes

**[numerical | FD curves from FD's receipt; accumulator from rh4_capture at the same (n, gain, delta, gamma) points]**
Figure `docs/figures/rh4-readings-vs-fd.png`, runner `experiments/2026-09-06/rh4_readings.py`, receipt `rh4_readings-checks.json`. Two register variants:
ideal (no decoherence) and realistic (gamma_r = gamma_phi = 1e-3, monitored). Capture t_c = 1/a per hold. n = 1, 2, 3.

## What the panels say

1. **Endpoint flags vs n.** FD trine/five as in FD-5. Accumulator with an ideal register: 0 at every n. Realistic register:
   0.035% -> 0.064% -> 0.087% for n = 1, 2, 3 -- the herald GROWS with n because the register waits longer.
2. **Accepted residual vs duration (eps = 1e-3, Delta_d = 1e-4, zero background).** FD five: 2e-13 / 3e-12 / 1.4e-11.
   Accumulator, ideal register: **1.5e-9 / 3e-10 / 1.5e-10** -- four orders ABOVE FD at n = 1, falling with n.
   Realistic register: 3e-7 / 5e-7 / 6e-7, rising with n.
3. **Register sweep** (n = 1, 1% gain, gamma = 1e-6): the accumulator's version of FD's dump-action tradeoff. Unmonitored
   loss crosses FD's 5.93e-6 residual immediately; monitored loss stays at it until 1e-3; dephasing crosses at ~1e-5.
4. **Total absorption vs gamma.** Ideal accumulator: total = FD's unobserved absorption exactly (heralded part zero);
   realistic register: a flat herald 3e-4..9e-4 replaces FD's 2.7e-3..2.9e-3.
5. **Paired operating points.** FD n = 1 at (0.265%, 5.93e-6). Accumulator, realistic register, same point:
   **(0.035%, 3.19e-5)** -- up and to the left: 7.6x less herald, 5.4x more unheralded residual. On these axes that is
   a TRADE, not a dominance. Only the ideal register (herald 0, residual 5.92e-6) dominates FD.
6. **Residual decomposition** (n = 1, eps = 1e-3, ideal register): above gamma ~ 1e-9 the accepted residual is unobserved
   absorption, as for FD; below it, the accumulator has a floor of **1.5e-9 from the recovered branch** that FD does not have.

## The floor, and why it is there  [derived + numerical]

recovered_error = c DB^2/6 to four digits (1.4937e-9 computed, 1.49e-9 measured), with DB_stack of the register map
**0.0182 at eps = 1e-3 and 0.0016 at eps = 1e-2: DB ~ delta/eps.** The gain leak is tight; the detuning-gain cross term
is not, and its relative weight is delta/eps. Hence the recovered-branch error c DB^2/6 ~ eps^2 (delta/eps)^2 = delta^2 --
**independent of gain**, ~1.2-1.5e-9 at delta_d/a = 1e-4, falling with n.
Mechanism: the leaked branch at dump k is E_k S_{k-1}...S_1 psi. It left the word after k loops and never received the
later loops' detuning cancellation. FD discards those branches; the five-loop word's detuning cancellation is a property
of the FULL survivor product only. The accumulator recovers partial-word states, and they carry the uncancelled
detuning of the partial products. A single 10 -> 2 isometry cannot absorb it: F = U H with H anisotropic is the
Knill-Laflamme obstruction again, at relative order delta/eps.

## Corrections to earlier statements today
- "The accumulator's operating point wants higher n" (RH-2A/RH-3 discussion): **only with an ideal register.** With a
  realistic register the herald and the residual both rise with n (panels 1, 2), because the register holds its shots
  longer. Which way n goes depends on register quality against background.
- "The accumulator's own floor is irrelevant": true at the operating point (gamma = 1e-6, where 1.5e-9 << 5.9e-6), false in
  a clean-background regime, where the delta^2 floor is the accumulator's limit and FD's conditional is 2e-13.

## What would remove the floor (open, not built)
Either a word whose PARTIAL products are detuning-cancelled (a stronger design constraint than FD-2's, so that every leaked
branch is as clean as the survivor), or a per-slot correction inside the register given delta-hat (the lock's estimate),
or accept it and note that at any background above ~1e-9 it is invisible. Kill for the first: if no five-or-more-loop word
with tight aggregate frame AND detuning-cancelled partial products exists, the floor is structural.

Replay: `~/.cache/seated-root-qec-venv/bin/python experiments/2026-09-06/rh4_readings.py`

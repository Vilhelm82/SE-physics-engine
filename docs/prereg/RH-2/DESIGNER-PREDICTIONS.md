# RH-2 designer predictions -- SEALED 2026-09-06 (Claire), written before rh2.py existed

Frozen against HEAD 866ff0e. Basis stated for each. Banked so they can be wrong on record.

RH-2a. Trine, gain-only: DB_stack in 5e-6 .. 2e-5 at all three eps (basis: compose_discard's
erasure_spread/erasure ~ 1.4e-5 is ~2 delta). Identity defect < 1e-12. c = 1.59e-3 at 1% (RH-1's number).
Five-loop, gain-only: DB_stack in 0.05 .. 0.5 -- the five-loop word was solved for detuning
cancellation, not for aggregate tightness, and its inverse stretched loops break the 60-degree frame.
Detuning-only DB_stack of order 0.5 .. 1 for both words (state-dependent leak, RH-1a).

RH-2b. Trine self-calibrated: 1 - F_avg in 3e-15 .. 1e-13 at eps_0, i.e. fork (a); the floor formula
c delta^2/6 gives ~3e-14 and the conditional is 2.5e-15, so the prediction is that the tight-frame
deviation, not the conditional error, sets the floor. Per-dump line reproduces RH-1b's 0.667 ratio.
Five-loop self-calibrated: "between" -- c delta^2/6 with delta ~ 0.1..0.5 gives 1e-6 .. 3e-5, well
below the discard 1.3e-3 but 10 orders above the conditional.
Fixed-calibration excess: exponent 2.0 +- 0.1 in (eps - eps_0), prefactor of order c.

RH-2c. Detuning passes through: the accumulator changes the 2.4e-8 detuning error by less than a
factor 2 in either direction.

Where I expect to be wrong, if anywhere: the five-loop DB_stack. The reflections between dumps
reorient effective directions (a, 2a-b+pi, a-b) for three loops; with five loops and two inverses I have
not computed the effective frame and the number could land anywhere from 1e-3 to 0.9.

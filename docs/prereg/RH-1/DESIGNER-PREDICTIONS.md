# RH-1 designer predictions -- SEALED 2026-09-06 (Claire)

Do not open before rh1a, rh1b, rh1c, rh1d have all run. Codex: this file is not part of your spec.
Frozen against HEAD 6cd6630. These are guesses with a stated basis, banked so they can be wrong on record.

RH-1a. Trine, gain-only: DB_k < 1e-4 at all three dumps (basis: total spread/mean 1.4e-5). FD, gain-only:
DB_k < 1e-3 at the three forward loops; the two stretched inverse loops LESS blind, DB_k in 1e-3..1e-1
(basis: they were solved for detuning cancellation, not for the tight-frame condition). Detuning-only
DB_k of order one everywhere (basis: spread/mean 1.4).

RH-1b. Trine: fork (a). Recovered unconditional within 10x of discard conditional at eps_0, i.e. 1e-14..1e-13.
Excess exponent in (eps - eps_0): 2.0 +- 0.1, prefactor of order c ~ 1e-3. FD: fork (a) only if RH-1a passes
at all five dumps; otherwise "between", dominated by the worst DB_k, and I expect that to be an inverse loop.
Replay variant indistinguishable from C_k at reported precision.

RH-1c. Fork (c), leaning (a): alignment 0.6..0.9 -- the window moves but a floor remains. Asym_N linear
through zero with nonzero slope. Quadratic exponent 2.00 +- 0.05 over seven points. Shot count for 1e-5
resolution at eps = 1e-2: 1e8..1e9 (basis: the whole delta effect is 1.8e-9 of the rate; Fisher info is small).
Dither's own erasure at delta = 0: data-blind (DB < 1e-4), because it is gain-channel bright population.

RH-1d(i). Trine, delta on D_d (half common): unitary fraction > 0.9 but Zc in 0.1..0.5 -- tilted, because the
common half is a gap change. FD pure differential: fork (a), Zc < 1e-2. FD pure common: fork (c) or (b),
non-unitary fraction > 10%.
RH-1d(ii). Gauge identity holds to 1e-12. Window returns to 0 under exact lock. Residual exponent 2.0.
RH-1d(iii). Common-mode lock does NOT return the window; erasure changes by less than 10%.

If RH-1a fork (b) fires on the trine, I was wrong about what 1.4e-5 on the total implies per dump and the
whole compressor reading needs re-deriving with the cumulative F_k in place of the per-loop E_k.

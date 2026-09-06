# RH-1 first passive-reference trial

Approved direction: [single-angle pullback](../../../../../2026-09-06-RH-1-cella-constraint-surface.md), section 6, and [physical-rest requirements](../../../../../2026-09-06-RH-1-restoring-dynamics.md).

This is a bounded research trial in the existing workspace. It does not change the gate implementation or the original RH-1 forecasts.

## Candidate fixed before numerical evaluation

Use a torsional spring and N independent reference two-level systems with a static transverse splitting and a longitudinal splitting sensitive to both the actual detuning and the angle:

`H_ref = sum_j [Omega X_j + (delta_ref - g s) Z_j]/2`, plus `N k0 s^2/2`.

Declare thermal equilibrium of the fast reference, a slow damped angle, and matched environmental coupling `delta_ref=delta_d`. Calibrate one small-signal coefficient at zero native gain/detuning; freeze it for all later trials. Neither a computed F nor a gate optimum enters the force law.

## Tasks and fixed checks

1. Derive the free energy, force, stiffness, small-signal calibration and environmental work. Require a stable minimum and return from both signs. Separate the free-energy Lyapunov law from microscopic internal energy.
2. Build a standalone runner. Evaluate both signs of detuning at `1e-4`, the wider `1e-3` case, gain `0` and `0.01`, and a 5% mismatch of reference sensitivity. No later fitted terms.
3. Test a moving minimum with relaxation time 100 native words, 1000-word drift plateaus, and native laboratory-frame integration during representative moving-angle gates. Compare the same fixed angle and record transient errors, erasure and completeness.
4. Include equilibrium angle fluctuations at dimensionless temperature `0.01`, `N=1e6,1e8,1e10`, plus the native reciprocal logical torque. Evaluate the full surviving CP map and its action on logical coherences; a mean-angle success is insufficient. Mechanical noise is initially slow relative to a word; label this approximation and its limits.
5. Retain as a candidate only if the declared model centers and improves the noisy conditional channel at a stated resource level. Reject any claim of exact arbitrary-data preservation at finite angle variance. Report power/work, return time, missing physical transduction and bath assumptions even if the numerical candidate passes.
6. Save the derivation and reproducible receipt; update the dependency ledger and current project handoff. Verify source hashes, targeted numerical checks and diff whitespace.

No fitted optimum, corrective postselection, new recovery channel or post-result nonlinear compensator is admissible in this trial.

## Completion, 6 September 2026

Tasks 1–6 completed for this bounded model trial. [Results](../../../../../2026-09-06-RH-1-passive-reference.md) retain the candidate's first-order matching, the failing small reference, thermal coherence loss, transient worsening, mechanical free-energy accounting and unresolved physical implementation. The constant comparison uses the linear correction to the mean of all future plateau fields; it is not labeled an exact optimizer. Four independent tests and 24 runner checks passed. An independent read-only review found no actionable mathematical, channel or implementation errors. The next physical-platform decision is documented, without promoting this conditional trial into a device claim.

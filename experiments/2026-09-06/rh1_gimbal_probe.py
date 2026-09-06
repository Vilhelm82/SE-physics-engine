"""Exploratory RH-1 angle servo; separate from the sealed RH-1 tests.

Plant: native trine, fixed lab detuning, instantaneous coherent P/Q dumps.
Controller receives only paired known-probe X counts (or exact expectations
in the explicitly ideal run). It does not receive plant parameters or fidelity.
One epoch is a whole calibration batch; drift is static within that batch.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))  # repo root on sys.path (reorg 2026-09-06)
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import minimize_scalar

from rlq.rh1_common import TRINE, TAU, compose_discard, trine_loop

X = np.array([[0., 1.], [1., 0.]])
GAIN = 0.5
SEED = 20260906
STATIC_SHIFT = -9.526823744383298e-5  # Earlier, supplied-estimate calibration.
DETUNINGS = [1e-4, -1e-4, 1.2e-4]


def plant(shift, dd, eps=.01):
    k, _, metrics = compose_discard(
        [trine_loop(float(a + shift), eps, dd) for a in TRINE])
    survive = np.sum(abs(k)**2, axis=0)
    means = np.real(np.diag(k.conj().T @ X @ k)) / survive
    return survive, means, metrics


def readout(survive, means, shots, rng):
    if shots is None:
        return float((means[0] - means[1]) / 2), None
    counts = []
    for p, m in zip(survive, means):
        accepted = int(rng.binomial(shots // 2, np.clip(p, 0, 1)))
        plus = int(rng.binomial(accepted, np.clip((1 + m) / 2, 0, 1)))
        counts.append(dict(attempted=shots // 2, accepted=accepted, plus=plus))
    if any(c['accepted'] == 0 for c in counts):
        return None, counts
    estimates = [2 * c['plus'] / c['accepted'] - 1 for c in counts]
    return float((estimates[0] - estimates[1]) / 2), counts


def update(shift, measured_contrast):
    # State and sensor only: no true detuning, Kraus matrix or error metric.
    if measured_contrast is None:
        return shift
    return float(shift - GAIN * measured_contrast / 4)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', type=Path)
    args = parser.parse_args()
    # Strong comparison: a fixed angle fitted using the whole future trace.
    # This information is never available to update().
    best_fixed = minimize_scalar(
        lambda s: np.mean([plant(float(s), d)[2]['conditional'] for d in DETUNINGS]),
        bounds=(-.001, .001), method='bounded', options=dict(xatol=1e-12))
    if not best_fixed.success:
        raise RuntimeError(best_fixed.message)
    checks, slopes = [], []
    for eps in (0., .01):
        for dd in (-1e-4, 1e-4):
            readings = [plant(s, dd, eps)[1] for s in (-1e-5, 1e-5)]
            contrasts = [(m[0] - m[1]) / 2 for m in readings]
            slope = float((contrasts[1] - contrasts[0]) / 2e-5)
            slopes.append(dict(eps=eps, dd=dd, slope=slope))
            checks.append(dict(name=f'local actuator slope {eps} {dd}',
                               passed=abs(slope - 4) < .001))
    runs = []
    for shots in (None, 10**6, 10**8):
        rng = np.random.default_rng(SEED)
        shift = 0.
        rows = []
        for n, dd in enumerate(np.repeat(DETUNINGS, 20)):
            survive, means, metrics = plant(shift, float(dd))
            measured, counts = readout(survive, means, shots, rng)
            baseline = plant(STATIC_SHIFT, float(dd))[2]
            optimized = plant(float(best_fixed.x), float(dd))[2]
            rows.append(dict(epoch=n, block=n // 20, dd=float(dd),
                             command=shift, measured_contrast=measured,
                             counts=counts, metrics=metrics,
                             static_conditional=baseline['conditional'],
                             optimized_fixed_conditional=optimized['conditional']))
            shift = update(shift, measured)
        settled = [r for r in rows if r['epoch'] % 20 >= 12]
        summary = dict(
            mean_settled_conditional=float(np.mean(
                [r['metrics']['conditional'] for r in settled])),
            mean_settled_static_conditional=float(np.mean(
                [r['static_conditional'] for r in settled])),
            mean_settled_optimized_fixed_conditional=float(np.mean(
                [r['optimized_fixed_conditional'] for r in settled])),
            final_conditional_by_block=[rows[i]['metrics']['conditional']
                                        for i in (19, 39, 59)],
            final_commands_by_block=[rows[i]['command'] for i in (19, 39, 59)],
            max_completeness_defect=max(
                r['metrics']['completeness_defect'] for r in rows),
            max_erasure=max(r['metrics']['erasure'] for r in rows),
            attempted_probe_shots=None if shots is None else 60 * shots,
            word_duration=9 * TAU,
            serial_word_time_per_epoch=None if shots is None else shots * 9 * TAU)
        checks.append(dict(name=f'instrument completeness shots={shots}',
                           passed=summary['max_completeness_defect'] < 1e-10))
        if shots is None:
            checks.append(dict(name='ideal readout reacquires after both steps',
                               passed=max(summary['final_conditional_by_block']) < 2e-12))
        runs.append(dict(attempted_shots_per_epoch=shots, summary=summary, rows=rows))
    result = dict(
        status='exploratory; one seeded trajectory per finite-shot setting',
        assumptions=[
            'Exact known d/r probe preparation and accepted-code X readout.',
            'Equal attempted shots per probe; every erasure is counted.',
            'Probes share the data-word plant parameters; no unknown data readout.',
            'Detuning and gain are constant within each full calibration batch.',
            'The next epoch receives the new command; no additional latency.',
            'No command slew, readout/preparation time or decoherence modeled.',
            'Expectation mode is an ideal sensor, not an implemented estimator.'],
        gain=GAIN, seed=SEED, fixed_plant_gain_error=.01,
        optimized_fixed_comparison=dict(command=float(best_fixed.x),
            mean_conditional=float(best_fixed.fun), bounds=[-.001, .001],
            status='Local bounded fit using all three future plateaus; not a global certificate'),
        slopes=slopes, runs=runs, checks=checks,
        provenance=dict(date='2026-09-06',
            head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
            sha256={p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                    for p in ('experiments/2026-09-06/rh1_gimbal_probe.py', 'rlq/rh1_common.py',
                              'rlq/reflection_loop_reference_echo.py',
                              'rlq/reflection_loop_dynamics.py')}))
    if args.json:
        args.json.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(slopes=slopes, checks=checks,
                          runs=[{k: v for k, v in r.items() if k != 'rows'}
                                for r in runs]), indent=2))
    if any(not c['passed'] for c in checks):
        raise SystemExit(1)


if __name__ == '__main__':
    main()

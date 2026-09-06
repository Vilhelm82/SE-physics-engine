"""Bounded Hamiltonian-level passive-reference trial, not a hardware controller.

The reference feels delta_ref through its physical Hamiltonian. An offline native
calibration fixes one constant; F and gate fidelity never enter its evolution.
Thermal gate averages assume an initially independent, nearly frozen angle.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2]))  # repo root on sys.path (reorg 2026-09-06)
import argparse
from dataclasses import dataclass, replace
from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from rh1_cella_surface import signed_constraint
from rlq.rh1_common import (G, I2, I4, TAU, TRINE, compose_discard,
                        haar_quadrature, trine_loop)
from rlq.reflection_loop_dynamics import hamiltonian
from rlq.reflection_loop_reference_echo import rotation

WORD_TIME = 9*TAU


@dataclass(frozen=True)
class Reference:
    k0: float
    gap: float = 1.
    coupling: float = 1.
    temperature: float = .01
    number: float = 1.
    rate: float = 1/(100*WORD_TIME)

    def __post_init__(self):
        if not (self.gap > 0 and self.number > 0 and self.temperature >= 0
                and self.rate > 0 and self.k0 > self.coupling**2/(2*self.gap)):
            raise ValueError('positive parameters and globally convex reference required')

    def spin(self, s, delta):
        detuning = delta-self.coupling*s
        energy = np.hypot(self.gap, detuning)
        polarization = 1. if self.temperature == 0 else np.tanh(energy/(2*self.temperature))
        return detuning, energy, polarization

    def potential(self, s, delta):
        """Isothermal free energy, not microscopic internal energy."""
        _, energy, _ = self.spin(s, delta)
        thermal = 0. if self.temperature == 0 else self.temperature*np.log1p(np.exp(-energy/self.temperature))
        return self.number*(self.k0*s*s/2-energy/2-thermal)

    def force(self, s, delta):
        detuning, energy, polarization = self.spin(s, delta)
        return -self.number*(self.k0*s+self.coupling*detuning*polarization/(2*energy))

    def field_derivative(self, s, delta):
        detuning, energy, polarization = self.spin(s, delta)
        return -self.number*detuning*polarization/(2*energy)

    def stiffness(self, s, delta):
        detuning, energy, polarization = self.spin(s, delta)
        susceptibility = self.gap**2*polarization/energy**3
        if self.temperature > 0:
            susceptibility += (detuning/energy)**2*(1-polarization**2)/(2*self.temperature)
        return self.number*(self.k0-self.coupling**2*susceptibility/2)

    @property
    def friction(self):
        return self.stiffness(0., 0.)/self.rate

    def rest(self, delta):
        bound = abs(self.coupling)/(2*self.k0)+1.
        return brentq(lambda s:self.force(s, delta), -bound, bound, xtol=5e-16)


def gaussian_nominal_error(mean, variance):
    # Stable form of [1-cos(4 mean) exp(-8 variance)]/3.
    return (-np.expm1(-8*variance)+2*np.sin(2*mean)**2*np.exp(-8*variance))/3


@lru_cache(None)
def native_word(shift, delta, gain=.01):
    return compose_discard([trine_loop(float(a+shift), gain, delta) for a in TRINE])


def calibrate():
    estimates = []
    for h in (1e-6, 5e-7):
        fs = (signed_constraint(TRINE+h, 0., 0.)[0]
              - signed_constraint(TRINE-h, 0., 0.)[0])/(2*h)
        fd = (signed_constraint(TRINE, 0., h)[0]
              - signed_constraint(TRINE, 0., -h)[0])/(2*h)
        estimates.append(fd/fs)
    c = estimates[-1]
    chi = np.tanh(1/(2*.01))/2
    return Reference(k0=chi*(1+1/c)), estimates


def thermal_channel(mean, variance, delta, nodes=13, gain=.01):
    """Average full native survival/leak maps before conditioning on success."""
    points, weights = np.polynomial.hermite.hermgauss(nodes)
    weights /= np.sqrt(np.pi)
    states, haar_weights = haar_quadrature()
    perpendicular = np.column_stack((-states[:, 1].conj(), states[:, 0].conj()))
    error, survival = np.zeros(len(states)), np.zeros(len(states))
    total, choi = np.zeros((2, 2), complex), np.zeros((4, 4), complex)
    erasure, mean_success, error_numerator = 0., 0., 0.
    for point, weight in zip(points, weights):
        k, leaks, metrics = native_word(float(mean+np.sqrt(2*variance)*point), delta, gain)
        a = G@k
        tf = a-np.trace(a)*I2/2
        error += weight*abs(np.sum(perpendicular.conj()*(states@tf.T), axis=1))**2
        survival += weight*np.sum(abs(states@a.T)**2, axis=1)
        total += weight*(k.conj().T@k+sum(f.conj().T@f for f in leaks))
        vector = a.ravel(order='F')
        choi += weight*np.outer(vector, vector.conj())
        erasure += weight*metrics['erasure']
        mean_success += weight*metrics['success']
        error_numerator += weight*np.linalg.norm(tf)**2/3
    return dict(conditional=float(haar_weights@(error/survival)),
        success_weighted=float(error_numerator/mean_success),
        erasure=float(erasure), success=float(mean_success),
        unconditional=float(error_numerator+erasure/2),
        completeness_defect=float(np.linalg.norm(total-I2)),
        survival_choi_min_eigenvalue=float(np.linalg.eigvalsh(choi)[0]),
        accepted_entanglement_fidelity=float((np.array([1,0,0,1])@choi@np.array([1,0,0,1])).real/(4*mean_success)))


def environment(t):
    """Actual field: plateaus joined by 100-word smooth ramps."""
    word = t/WORD_TIME
    for start, before, after in ((1000., 1e-4, -1e-4), (2000., -1e-4, 1.2e-4)):
        if word < start:
            return before, 0.
        if word < start+100:
            u = (word-start)/100
            return before+(after-before)*(3*u*u-2*u**3), (after-before)*6*u*(1-u)/(100*WORD_TIME)
    return 1.2e-4, 0.


def evolve(ref):
    def rhs(t, y):
        delta, velocity = environment(t)
        ds = ref.force(y[0], delta)/ref.friction
        return [ds, ref.friction*ds*ds, ref.field_derivative(y[0], delta)*velocity]
    sol = solve_ivp(rhs, (0, 3000*WORD_TIME), [0., 0., 0.],
        method='DOP853', rtol=1e-10, atol=1e-14, max_step=5*WORD_TIME, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def moving_native_word(start, angle, delta, gain=.01, max_step=.18):
    """Integrate the laboratory Hamiltonian; angle motion is not omitted."""
    loops = []
    for leg, base_angle in enumerate(TRINE):
        u = I4.copy()
        for arc in range(3):
            offset = start+(3*leg+arc)*TAU
            def rhs(local, flat):
                t = offset+local
                rot = rotation(float(base_angle+angle(t)))
                h = (1+gain)*rot@hamiltonian(arc, local, TAU)@rot.T
                h += np.diag([0., 0., delta(t), 0.])
                return (-1j*h@flat.reshape(4, 4)).ravel()
            sol = solve_ivp(rhs, (0, TAU), u.ravel(), method='DOP853',
                rtol=2e-12, atol=2e-14, max_step=max_step)
            if not sol.success:
                raise RuntimeError(sol.message)
            u = sol.y[:, -1].reshape(4, 4)
        loops.append(u)
    return compose_discard(loops)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', type=Path)
    parser.add_argument('--plot', type=Path)
    args = parser.parse_args()
    ref, slopes = calibrate()
    c = slopes[-1]
    checks = []
    def check(name, value, limit):
        checks.append(dict(name=name, value=float(value), limit=limit,
                           passed=bool(np.isfinite(value) and abs(value) <= limit)))
    check('calibration step convergence', slopes[0]/c-1, 1e-6)
    static = []
    for gain in (0., .01):
        for delta in (-1e-4, 1e-4, 1e-3):
            rest = ref.rest(delta)
            ideal_root = brentq(lambda s:signed_constraint(TRINE+s, gain, delta)[0], -.003, .003)
            static.append(dict(gain=gain, delta=delta, rest=rest, offline_F_root=ideal_root,
                F_at_rest=signed_constraint(TRINE+rest, gain, delta)[0],
                reference=native_word(rest, delta, gain)[2],
                fixed_zero=native_word(0., delta, gain)[2],
                mismatched_5_percent=native_word(ref.rest(1.05*delta), delta, gain)[2]))
    sol = evolve(ref)
    timeline = np.linspace(0, 3000*WORD_TIME, 601)
    states = sol.sol(timeline)
    balance = [ref.potential(s, environment(t)[0])-ref.potential(0., 1e-4)+heat-work
               for t, (s, heat, work) in zip(timeline, states.T)]
    check('isothermal free energy plus friction equals field work', max(abs(np.array(balance))), 2e-13)
    return_fraction = abs(sol.sol(500*WORD_TIME)[0]-ref.rest(1e-4))/abs(ref.rest(1e-4))
    check('settles within one percent after five relaxation times', return_fraction, .01)
    moving = []
    angle = lambda t:sol.sol(t)[0]
    field = lambda t:environment(t)[0]
    for word in (0, 50, 500, 1000, 1050, 1100, 1200, 1500, 2050, 2100, 2200, 2999):
        t = word*WORD_TIME
        k, _, m = moving_native_word(t, angle, field)
        moving.append(dict(start_word=word, angle=float(angle(t)), delta=float(field(t)),
            moving=m, fixed_zero=moving_native_word(t, lambda _:0., field)[2]))
        check(f'moving gate completeness at {word}', m['completeness_defect'], 2e-10)
    stationary = moving_native_word(0., lambda _:.0001, lambda _:1e-4)[0]
    check('moving integrator recovers frozen native word',
          np.linalg.norm(stationary-native_word(.0001, 1e-4)[0]), 2e-10)
    coarse = moving_native_word(1050*WORD_TIME, angle, field)[0]
    fine = moving_native_word(1050*WORD_TIME, angle, field, max_step=.09)[0]
    check('moving gate time-step convergence', np.linalg.norm(coarse-fine), 2e-10)
    resources = []
    for number in (1e6, 1e8, 1e10):
        local = replace(ref, number=number)
        center = local.rest(1e-4)
        stiffness = local.stiffness(center, 1e-4)
        variance = local.temperature/stiffness
        channel = thermal_channel(center, variance, 1e-4)
        check(f'thermal full instrument completeness N={number:g}', channel['completeness_defect'], 2e-10)
        check(f'thermal survival CP N={number:g}', min(0., channel['survival_choi_min_eigenvalue']), 2e-12)
        resources.append(dict(N=number, stiffness=stiffness, angle_rms=float(np.sqrt(variance)),
            variance=variance, pure_nominal_noise_error=gaussian_nominal_error(0., variance),
            nominal_XZ_visibility=float(np.exp(-8*variance)),
            worst_nominal_mean_data_load_bias=2/(WORD_TIME*stiffness),
            displacement_free_energy_1e_4=stiffness*1e-8/2,
            mechanical_friction_dissipation=number*float(sol.y[1,-1]),
            environmental_free_energy_work=number*float(sol.y[2,-1]), channel=channel))
    medium = resources[1]
    check('thermal quadrature convergence',
        thermal_channel(ref.rest(1e-4), medium['variance'], 1e-4, nodes=21)['conditional']
        - medium['channel']['conditional'], 2e-13)
    # A generous constant baseline gets all three future plateau detunings.
    plateau_deltas = [1e-4, -1e-4, 1.2e-4]
    fixed_shift = -c*float(np.mean(plateau_deltas))
    paired = [dict(delta=d, reference=thermal_channel(ref.rest(d), ref.temperature/(1e8*ref.stiffness(ref.rest(d),d)),d),
                   future_mean_constant=native_word(fixed_shift,d)[2]) for d in plateau_deltas]
    paths = [Path(__file__), Path('tests/test_rh1_passive_reference.py'), Path('experiments/2026-09-06/rh1_cella_surface.py'),
        Path('rlq/rh1_common.py'), Path('rlq/reflection_loop_dynamics.py'),
        Path('rlq/reflection_loop_reference_echo.py'), Path('rlq/reflection_loop_composite.py'),
        Path('docs/archive/superpowers/plans/2026-09-06-rh1-passive-reference.md')]
    output = dict(date='2026-09-06', model='Equilibrium independent two-level reference plus torsional spring; proposed coupling, no hardware realization.',
        calibration=dict(gain=0., delta=0., slope=c, steps=slopes, frozen_k0=ref.k0),
        parameters=dict(gap=ref.gap, coupling=ref.coupling, temperature=ref.temperature,
            rate=ref.rate, relaxation_words=100, word_time=WORD_TIME,
            reference_relaxation_assumption='fast compared with angle and 100-word environmental ramps; not a simulated microscopic bath'),
        static=static, moving=moving, resources=resources, plateau_comparison=paired,
        comparison_constant=fixed_shift,
        deterministic_balance=dict(max_defect=max(abs(np.array(balance))),
            heat_per_reference=float(sol.y[1,-1]), work_per_reference=float(sol.y[2,-1]),
            return_fraction_after_500_words=return_fraction),
        scope=['No gate optimum or F enters reference evolution.',
            'A matched environmental coupling, prescribed Hamiltonian, fast Gibbs equilibration and a slow classical angle are declared inputs.',
            'Thermal channels use quasi-static local Gaussian angles, initially independent of the data; moving gate tests are deterministic.',
            'Finite angle variance dephases unknown logical coherences. Data torque is bounded, not removed; repeated joint data-reference-bath dynamics remain untested.',
            'Dissipation is mechanical excess friction only. Refrigeration, bias-field generation and full microscopic reference-bath heat are not priced.'],
        checks=checks, sha256={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
    if args.json:
        args.json.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.3), layout='constrained')
        axes[0].plot(timeline/WORD_TIME, states[0]*1e6, label='Physical reference angle', color='#176e80')
        axes[0].plot(timeline/WORD_TIME, [-c*environment(t)[0]*1e6 for t in timeline],
                     '--', label='Linear gate correction', color='#b96e27')
        for left in (1000, 2000):
            axes[0].axvspan(left, left+100, color='#b96e27', alpha=.12)
        axes[0].set(xlabel='Elapsed native words', ylabel='Angle (microradians)',
                    title='The resting angle follows the field, with lag')
        axes[0].legend(frameon=False, fontsize=9)
        axes[1].loglog([r['N'] for r in resources], [r['channel']['conditional'] for r in resources],
                       'o-', color='#176e80', label='Settled thermal reference, delta=+1e-4')
        axes[1].axhline(native_word(0., 1e-4)[2]['conditional'], color='#b96e27',
                        linestyle='--', label='Fixed zero angle, same field')
        axes[1].set(xlabel='Number of reference spins', ylabel='Conditional gate infidelity',
                    title='A small reference is too noisy')
        axes[1].legend(frameon=False, fontsize=9)
        for ax in axes:
            ax.grid(alpha=.15)
            ax.spines[['top', 'right']].set_visible(False)
        fig.suptitle('Proposed equilibrium coupling | T=0.01 | response time=100 words', fontsize=12)
        fig.savefig(args.plot, dpi=160)
        plt.close(fig)
    print(json.dumps(dict(calibration=output['calibration'], checks_passed=sum(c['passed'] for c in checks),
        checks_total=len(checks), static=[dict(gain=s['gain'],delta=s['delta'],rest=s['rest'],
        corrected=s['reference']['conditional'],fixed=s['fixed_zero']['conditional']) for s in static],
        resources=[dict(N=r['N'],conditional=r['channel']['conditional'],rms=r['angle_rms']) for r in resources],
        plateau_mean=dict(reference=float(np.mean([p['reference']['conditional'] for p in paired])),
        future_constant=float(np.mean([p['future_mean_constant']['conditional'] for p in paired]))),
        failed=[c for c in checks if not c['passed']]), indent=2))
    if not all(c['passed'] for c in checks):
        raise SystemExit(1)


if __name__ == '__main__':
    main()

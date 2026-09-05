#!/usr/bin/env python3
"""Finite reference-access reflection gates, signed echo, and joint response jets.

Extends the committed real composite. All sign conjugations compile Hamiltonians;
they are not instantaneous gates. The full physical waveform starts/ends at zero.
"""
import argparse
from dataclasses import dataclass
import json

import numpy as np
import sympy as sp

from reflection_loop_dynamics import R, embed
from reflection_loop_composite import (
    Stage, P, Q, REFERENCE, composite, schedule, coefficients, ode,
)

DD = np.diag([1, 1, -1, 1])
DR = np.diag([1, 1, 1, -1])
DB = DD @ DR
TARGET = embed(R)
LOGICAL = np.diag([-1, 1])
EYE = np.eye(4, dtype=complex)
POWERS = ((0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2))


def rotation(alpha):
    out = np.eye(4)
    c, s = np.cos(alpha), np.sin(alpha)
    out[2:, 2:] = [[c, -s], [s, c]]
    return out


O = rotation(np.pi/3)


def ended_composite(n=1, gap=1.0, slew=None):
    """Minimum-time zero-area entry/exit ramps at the given amplitude and slew."""
    slew = 10*gap*gap if slew is None else slew
    tau = coefficients(n)['primitive_time']/(3*gap)
    if slew < gap*np.pi/(2*tau):
        raise ValueError('Slew cap is below the original loop arc requirement')
    b = gap/np.sqrt(2)
    entry = [Stage(b/slew, q_start=0, q_end=-b),
             Stage((gap+b)/slew, q_start=-b, q_end=gap)]
    exit_stages = [Stage((gap+b)/slew, q_start=gap, q_end=-b),
                   Stage(b/slew, q_start=-b, q_end=0)]
    return entry + schedule(n, gap, ramp=2*gap/slew) + exit_stages


@dataclass(frozen=True)
class FramedStage:
    base: Stage
    basis: np.ndarray
    stretch: float = 1.0
    reverse: bool = False

    @property
    def duration(self):
        return self.stretch*self.base.duration

    def local(self, t):
        s = t/self.stretch
        return self.base.duration-s if self.reverse else s

    def h(self, t):
        sign = -1 if self.reverse else 1
        return sign*self.basis @ self.base.h(self.local(t)) @ self.basis.T/self.stretch

    def loss_projector(self, t):
        return self.basis @ self.base.loss_projector(self.local(t)) @ self.basis.T

    def u(self, t):
        u = self.base.u(self.local(t))
        if self.reverse:
            u = u @ self.base.u(self.base.duration).conj().T
        return self.basis @ u @ self.basis.T


def trine_stages(n=1, gap=1.0, slew=None):
    base = ended_composite(n, gap, slew)
    return [FramedStage(stage, basis) for basis in (O, DD @ O, O) for stage in base]


def echo_stages(n=1, gap=1.0, slew=None):
    tri = trine_stages(n, gap, slew)
    out = []
    for gauge, stretch, reverse in ((EYE.real, 1, False), (DD, 2, True),
                                     (DB, 2, False), (DR, 2, True),
                                     (EYE.real, 1, False)):
        for stage in (list(reversed(tri)) if reverse else tri):
            out.append(FramedStage(stage.base, gauge @ stage.basis, stretch, reverse))
    return out


def trine_product(a):
    return a @ DD @ a @ DD @ a


def echo_product(t1, t2, tminus2):
    inv = tminus2.conj().T
    return t1 @ DR @ inv @ DR @ DB @ t2 @ DB @ DD @ inv @ DD @ t1


def exact_gain(error, kind='echo', n=1):
    a = O @ composite(error, n) @ O.T
    t = trine_product(a)
    return t if kind == 'trine' else echo_product(t, t, t)


def metrics(u, dissipative=False):
    """Haar average logical error, with leakage counted and global phase removed."""
    a = LOGICAL @ u[2:, 2:]
    b = a-np.trace(a)*np.eye(2)/2
    leak = float(np.linalg.norm(u[:2, 2:])**2/2)
    deficit = float(1-np.linalg.norm(a)**2/2) if dissipative else leak
    return dict(average_leakage=leak,
                average_infidelity=float(deficit+np.linalg.norm(b)**2/3),
                logical_traceless_norm=float(np.linalg.norm(b)),
                logical_global_phase=float(np.angle(np.trace(a))))


def integrate_block(error=0.0, delta=0.0, loss=0.0, steps=100, n=1):
    stages = [FramedStage(s, O) for s in ended_composite(n)]
    return ode(stages, error=error, detuning=delta, loss=loss, steps_per_time=steps)


def finite_gate(error=0.0, delta=0.0, loss=0.0, kind='echo', steps=100, n=1):
    t1 = trine_product(integrate_block(error, delta, loss, steps, n))
    if kind == 'trine':
        return t1
    t2 = trine_product(integrate_block(error, 2*delta, 2*loss, steps, n))
    tm2 = trine_product(integrate_block(error, -2*delta, 2*loss, steps, n))
    return echo_product(t1, t2, tm2)


def response_integrals(stages, perturbations, nodes=64):
    """Exact nominal propagators; quadrature of all supplied Hermitian errors."""
    points, weights = np.polynomial.legendre.leggauss(nodes)
    matrices = np.asarray(perturbations, complex)
    moments = np.zeros_like(matrices)
    loss = np.zeros((4, 4), complex)
    gain = np.zeros((4, 4), complex)
    prefix = EYE.copy()
    parallel = 0.0
    for stage in stages:
        for point, weight in zip(points, weights):
            t = stage.duration*(point+1)/2
            wt = weight*stage.duration/2
            u = stage.u(t) @ prefix
            moments += wt*(u.conj().T @ matrices @ u)
            loss += wt*(u.conj().T @ stage.loss_projector(t) @ u)
            term = u.conj().T @ stage.h(t) @ u
            gain += wt*term
            parallel = max(parallel, float(np.linalg.norm(P @ term @ P)))
        prefix = stage.u(stage.duration) @ prefix
    return dict(endpoint=prefix, moments=moments, loss=loss, gain=gain,
                parallel_error=parallel)


def block_jets(steps=100, n=1):
    """Independent ODE coefficients in gain epsilon and reference frequency delta.

    Taylor coefficients include factorials, i.e. U20 multiplies epsilon**2.
    Only one rotated, zero-ended C block is integrated; subsequent compositions
    propagate these coefficients algebraically, including delta's sign in inverses.
    """
    state = np.zeros((len(POWERS), 4, 4), complex)
    state[0] = EYE
    position = {power: j for j, power in enumerate(POWERS)}
    for stage in [FramedStage(s, O) for s in ended_composite(n)]:
        count = max(8, int(np.ceil(stage.duration*steps)))
        dt = stage.duration/count

        def rhs(t, values):
            h = stage.h(t)
            result = -1j*h @ values
            for j, (p, q) in enumerate(POWERS):
                if p:
                    result[j] -= 1j*h @ values[position[(p-1, q)]]
                if q:
                    result[j] -= 1j*REFERENCE @ values[position[(p, q-1)]]
            return result

        for j in range(count):
            t = j*dt
            k1 = rhs(t, state)
            k2 = rhs(t+dt/2, state+dt*k1/2)
            k3 = rhs(t+dt/2, state+dt*k2/2)
            k4 = rhs(t+dt, state+dt*k3)
            state += dt*(k1+2*k2+2*k3+k4)/6
    return {power: state[j] for j, power in enumerate(POWERS)}


def jet_multiply(a, b):
    out = {}
    for p, q in POWERS:
        out[p, q] = sum((a[i, j] @ b[p-i, q-j]
                        for i in range(p+1) for j in range(q+1)
                        if (i, j) in a and (p-i, q-j) in b), np.zeros((4, 4), complex))
    return out


def jet_constant(a):
    return {power: a if power == (0, 0) else np.zeros((4, 4), complex)
            for power in POWERS}


def jet_chain(items):
    out = jet_constant(EYE)
    for item in items:
        out = jet_multiply(out, item)
    return out


def combined_jets(a):
    d, b, r = map(jet_constant, (DD, DB, DR))
    t = jet_chain((a, d, a, d, a))
    t2 = {key: value*2**key[1] for key, value in t.items()}
    inv = {key: value.conj().T*(-2)**key[1] for key, value in t.items()}
    e = jet_chain((t, r, inv, r, b, t2, b, d, inv, d, t))
    return t, e


def quartic_response(jet, target=TARGET):
    leak, logical = {}, {}
    for power in ((2, 0), (1, 1), (0, 2)):
        w = target.conj().T @ jet[power]
        leak[power] = w[:2, 2:]
        block = w[2:, 2:]
        logical[power] = block-np.trace(block)*np.eye(2)/2
    terms = {}
    for first in leak:
        for second in leak:
            power = (first[0]+second[0], first[1]+second[1])
            terms[power] = terms.get(power, 0.0) + float(
                (np.vdot(leak[first], leak[second])/2
                 + np.vdot(logical[first], logical[second])/3).real)
    return dict(terms={f'{p},{q}': v for (p, q), v in sorted(terms.items())},
                leakage=leak, logical=logical)


def predict(response, error, delta):
    return sum(value*error**int(key.split(',')[0])*delta**int(key.split(',')[1])
               for key, value in response['terms'].items())


def symbolic_checks():
    checks = {}

    def exact(name, actual, expected=0):
        if isinstance(actual, sp.MatrixBase):
            expected = sp.zeros(*actual.shape) if expected == 0 else expected
            values = list(actual-expected)
        else:
            values = [actual-expected]
        assert all(sp.simplify(sp.expand(x)) == 0 for x in values), name
        checks[name] = True

    o = sp.eye(4)
    o[2:, 2:] = sp.Matrix([[sp.Rational(1, 2), -sp.sqrt(3)/2],
                          [sp.sqrt(3)/2, sp.Rational(1, 2)]])
    d, b, dr = map(lambda x: sp.Matrix(x), (DD, DB, DR))
    j = sp.Matrix([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    target = sp.Matrix(TARGET.real.astype(int))
    s = sp.eye(4)-sp.I*sp.sqrt(3)*j/2-j*j/2
    c0 = target*s*s
    a0 = o*c0*o.T
    t0 = a0*d*a0*d*a0
    exact('trine_full_endpoint_involution', t0, target)
    exact('trine_endpoint_squared', t0*t0, sp.eye(4))
    rn = sp.symbols('r_n', real=True)
    wc2 = sp.zeros(4)
    wc2[:, 2] = sp.Matrix([sp.pi**2*rn*(-2*sp.sqrt(3)+3*sp.I)/12,
                           sp.pi**2*rn*(sp.sqrt(3)-6*sp.I)/12,
                           sp.sqrt(3)*sp.I*sp.pi**2*rn*rn/2, 0])
    a2 = o*c0*wc2*o.T
    t2 = a2*d*a0*d*a0+a0*d*a2*d*a0+a0*d*a0*d*a2
    kt = sp.sqrt(3)*sp.pi**2*rn/4
    target_column = sp.zeros(4, 2)
    target_column[:2, :] = kt*sp.Matrix([[0, sp.I], [-1, -2*sp.I]])
    target_column[2:, :] = 3*sp.sqrt(3)*sp.I*sp.pi**2*rn*rn*sp.eye(2)/4
    exact('trine_complete_second_order_gain_column_all_n',
          (target*t2)[:, 2:], target_column)
    exact('trine_and_echo_quartic_gain_infidelity',
          (target_column[:2, :].H*target_column[:2, :]).trace()/2,
          9*sp.pi**4*rn*rn/16)

    alpha = (sp.pi/3, 2*sp.pi/3, sp.pi/3)
    prefix = sp.eye(2)
    projectors = sp.zeros(2)
    response = sp.zeros(2)
    aa, tt, kr, ki = sp.symbols('A T real_overlap imag_overlap', real=True)
    for angle in alpha:
        rot = sp.Matrix([[sp.cos(angle), -sp.sin(angle)],
                         [sp.sin(angle), sp.cos(angle)]])
        vector = rot[:, 0]
        pr = vector*vector.T
        projectors += prefix.T*pr*prefix
        sn, cs = sp.sin(angle), sp.cos(angle)
        m = rot*sp.Matrix([[aa*sn*sn, sn*cs*(kr+sp.I*ki)],
                           [sn*cs*(kr-sp.I*ki), tt*cs*cs]])*rot.T
        response += prefix.T*m*prefix
        prefix = (sp.eye(2)-2*pr)*prefix
    exact('toggled_trine_projectors_are_isotropic', projectors, sp.Rational(3, 2)*sp.eye(2))
    expected = (9*aa+3*tt)*sp.eye(2)/8
    expected += sp.Matrix([[0, 3*sp.sqrt(3)*sp.I*ki/4],
                          [-3*sp.sqrt(3)*sp.I*ki/4, 0]])
    exact('full_trine_projected_reference_response', response, expected)
    exact('trine_reference_diagonal_balance', response[0, 0]-response[1, 1])

    # Arbitrary QP entries and logical off-diagonal entries, equal code diagonals.
    entries = sp.symbols('x0:12', real=True)
    m = sp.zeros(4)
    m[:2, :2] = sp.Matrix([[entries[0], entries[1]+sp.I*entries[2]],
                          [entries[1]-sp.I*entries[2], entries[3]]])
    for row in range(2):
        for col in range(2):
            z = entries[4+2*(row*2+col)]+sp.I*entries[5+2*(row*2+col)]
            m[row, col+2] = z
            m[col+2, row] = sp.conjugate(z)
    common, x, y = sp.symbols('common x y', real=True)
    m[2:, 2:] = sp.Matrix([[common, x+sp.I*y], [x-sp.I*y, common]])
    averaged = sum((g*m*g for g in (sp.eye(4), d, b, dr)), sp.zeros(4))
    exact('four_sign_average_full_code_columns', averaged[:, 2:],
          (4*common*sp.eye(4))[:, 2:])
    l = sp.Matrix(4, 2, lambda i, j: sp.symbols(f'l{i}{j}'))
    l[2:, :] = sp.I*common*sp.eye(2)
    effect = (2*l-d*l*d[2:, 2:]+b*l*b[2:, 2:]-dr*l*dr[2:, 2:])
    exact('echo_preserves_second_order_gain_column', effect, l)
    amp, slew = sp.symbols('a nu', positive=True)
    low = amp/sp.sqrt(2)
    exact('zero_area_entry_ramp',
          -low**2/(2*slew)+(amp**2-low**2)/(2*slew))
    exact('zero_area_minimum_duration', low/slew+(amp+low)/slew,
          (1+sp.sqrt(2))*amp/slew)
    return checks


def run_checks():
    symbolic = symbolic_checks()
    numerical = {}

    def close(name, actual, expected, tol=1e-9):
        error = float(np.max(np.abs(np.asarray(actual)-np.asarray(expected))))
        assert error <= tol, (name, error, tol)
        numerical[name] = dict(error=error, tolerance=tol)

    def bounded(name, value, limit):
        assert value <= limit, (name, value, limit)
        numerical[name] = dict(value=float(value), upper_bound=float(limit))

    base, tri, ech = ended_composite(), trine_stages(), echo_stages()
    for name, stages in (('trine', tri), ('echo', ech)):
        close(name+'_starts_at_zero', stages[0].h(0), 0)
        close(name+'_ends_at_zero', stages[-1].h(stages[-1].duration), 0)
        close(name+'_continuous_joins', max(
            np.linalg.norm(left.h(left.duration)-right.h(0))
            for left, right in zip(stages, stages[1:])), 0)
    entry = base[1].u(base[1].duration) @ base[0].u(base[0].duration)
    close('entry_zero_area_identity_on_full_state', entry, EYE)
    close('exact_gain_nominal_trine', exact_gain(0, 'trine'), TARGET)
    close('exact_gain_nominal_echo', exact_gain(0), TARGET)
    for kind in ('trine', 'echo'):
        large = metrics(exact_gain(0.004, kind))['logical_traceless_norm']
        small = metrics(exact_gain(0.002, kind))['logical_traceless_norm']
        bounded(kind+'_relative_gain_phase_starts_at_fourth_order', small/large, 0.09)

    noise = [np.diag([int(i == j) for i in range(4)]) for j in range(4)]
    noise.extend([np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]),
                  np.array([[0, 1j, 0, 0], [-1j, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])])
    tdata = response_integrals(tri, noise)
    edata = response_integrals(ech, noise)
    for j, moment in enumerate(edata['moments']):
        common = np.trace(moment[2:, 2:])/2
        close(f'echo_cancels_noise_generator_{j}_on_entire_code',
              moment @ P, common*P, 2e-9)
    for name, data in (('trine', tdata), ('echo', edata)):
        close(name+'_propagator_integral_endpoint', data['endpoint'], TARGET)
        close(name+'_first_order_gain_on_code', data['gain'] @ P, 0)
        bounded(name+'_nominal_parallel_transport', data['parallel_error'], 2e-12)
    xc = 3*coefficients()['exposure']
    close('trine_isotropic_loss_exposure', tdata['loss'][2:, 2:], 1.5*xc*np.eye(2))
    close('echo_isotropic_loss_exposure', edata['loss'][2:, 2:], 12*xc*np.eye(2))
    close('echo_first_order_loss_is_state_independent_erasure',
          edata['loss'] @ P, 12*xc*P, 2e-9)

    coarse_a = block_jets(60)
    fine_a = block_jets(120)
    tc, ec = combined_jets(coarse_a)
    tf, ef = combined_jets(fine_a)
    first = TARGET.conj().T @ ef[0, 1]
    phase_slope = np.trace(first[2:, 2:])/2
    close('jet_reference_response_is_global_only', first @ P, phase_slope*P, 3e-7)
    close('jet_reference_matches_independent_quadrature',
          1j*first @ P, edata['moments'][3] @ P, 3e-7)
    for name, jet in (('trine', tf), ('echo', ef)):
        w = TARGET.conj().T @ jet[2, 0]
        close(name+'_quadratic_gain_phase_is_global',
              w[2:, 2:]-np.trace(w[2:, 2:])*np.eye(2)/2, 0, 2e-6)
    close('second_order_gain_column_preserved_by_echo',
          ef[2, 0][:, 2:], tf[2, 0][:, 2:], 2e-6)
    response = quartic_response(ef)
    response_coarse = quartic_response(ec)
    close('mixed_gain_detuning_leakage_coefficient_vanishes',
          response['leakage'][1, 1], 0, 4e-7)
    close('mixed_gain_detuning_logical_coefficient_is_global',
          response['logical'][1, 1], 0, 1e-7)
    for key, value in response['terms'].items():
        scale = max(1.0, abs(value))
        close('joint_response_convergence_'+key,
              response_coarse['terms'][key]/scale, value/scale, 2e-4)
    bounded('jet_ode_step_convergence',
            np.linalg.norm(fine_a[0, 0]-O@composite()@O.T)
            /np.linalg.norm(coarse_a[0, 0]-O@composite()@O.T), 0.08)

    # One direct integration covers all 255 physical stages, including inverse,
    # stretched, and sign-conjugated controls. Product evaluation only uses A.
    eps, delta = 0.002, 2e-5
    coarse_direct = ode(ech, error=eps, detuning=delta, steps_per_time=60)
    direct = ode(ech, error=eps, detuning=delta, steps_per_time=120)
    product = finite_gate(eps, delta, steps=240)
    close('full_waveform_ode_matches_block_construction', direct, product, 1e-8)
    bounded('full_waveform_step_doubling_ratio',
            np.linalg.norm(direct-product)/np.linalg.norm(coarse_direct-product), 0.08)
    bounded('full_waveform_ode_unitarity', np.linalg.norm(direct.conj().T @ direct-EYE), 2e-9)
    witnesses = []
    for eps, delta in ((0.001, 1e-5), (0.0005, 5e-6), (0.0, 1e-5)):
        actual = metrics(finite_gate(eps, delta, steps=120))['average_infidelity']
        leading = predict(response, eps, delta)
        witnesses.append(dict(gain=eps, reference_offset=delta,
                              actual_infidelity=actual, quartic_prediction=leading))
    # Leading-order comparison uses the smaller joint point; the finite point is
    # deliberately not asserted equal to a Taylor truncation.
    bounded('joint_quartic_relative_error_at_small_point',
            abs(witnesses[1]['actual_infidelity']/witnesses[1]['quartic_prediction']-1), 0.08)
    bounded('joint_quartic_relative_error_improves',
            abs(witnesses[1]['actual_infidelity']/witnesses[1]['quartic_prediction']-1)
            /abs(witnesses[0]['actual_infidelity']/witnesses[0]['quartic_prediction']-1), 0.7)

    rate = 1e-8
    damped = finite_gate(loss=rate, steps=120)
    no_loss = finite_gate(steps=120)
    close('finite_absorbing_loss_slope',
          (metrics(damped, True)['average_infidelity']
           -metrics(no_loss, True)['average_infidelity'])/rate, 12*xc, 0.02)
    # This is an explicit witness, not a presumed advantage for every noise mix.
    example_gain, example_delta = 0.001, 0.0001
    bare = composite(example_gain)
    bare[3, 3] = np.exp(-1j*example_delta*sum(s.duration for s in schedule()))
    short = finite_gate(example_gain, example_delta, kind='trine', steps=120)
    full = finite_gate(example_gain, example_delta, steps=120)
    cases = dict(composite=metrics(bare), trine=metrics(short), echo=metrics(full))
    echo_advantage = cases['composite']['average_infidelity']-cases['echo']['average_infidelity']
    crossovers = dict(
        echo_vs_composite=echo_advantage/(11.5*xc),
        echo_vs_trine=(cases['trine']['average_infidelity']
                      -cases['echo']['average_infidelity'])/(10.5*xc),
        trine_vs_composite=(cases['composite']['average_infidelity']
                           -cases['trine']['average_infidelity'])/xc)
    structured_quartic = response['terms'].copy()
    structured_quartic.update({'1,3': 0.0, '3,1': 0.0})
    report = dict(symbolic_checks=symbolic, numerical_checks=numerical,
                  symbolic_count=len(symbolic), numerical_count=len(numerical),
                  quartic_joint_infidelity=structured_quartic,
                  raw_numerical_quartic=response['terms'], witnesses=witnesses,
                  resource=dict(composite_duration=sum(s.duration for s in schedule()),
                                ended_composite_duration=sum(s.duration for s in base),
                                trine_duration=sum(s.duration for s in tri),
                                echo_duration=sum(s.duration for s in ech),
                                composite_average_exposure=xc/2,
                                trine_average_exposure=1.5*xc,
                                echo_average_exposure=12*xc,
                                reference_global_phase_slope=float(phase_slope.imag)),
                  example=dict(gain=example_gain, reference_offset=example_delta,
                               metrics=cases, leading_loss_crossovers=crossovers),
                  scope='All clocks, gains, detuning signs, finite ramps and absorbing '
                        'loss are explicit. Taylor coefficients are not finite-noise identities.')
    # Preserve the full matrix response for arbitrary correlated gain/detuning
    # distributions, not only a scalar curve or one chosen error direction.
    report['second_order_error_matrices'] = {
        name: {f'{p},{q}': {'real': matrix.real.tolist(), 'imag': matrix.imag.tolist()}
               for (p, q), matrix in response[name].items()}
        for name in ('leakage', 'logical')}
    return report


def plot(path, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), layout='constrained')
    errors = np.geomspace(2e-4, 0.02, 100)
    for kind, color in (('trine', '#146c94'), ('echo', '#bc5a3c')):
        values = [metrics(exact_gain(e, kind))['average_infidelity'] for e in errors]
        axes[0].loglog(errors, values, label=kind.capitalize(), color=color)
    axes[0].loglog(errors, [metrics(composite(e))['average_infidelity'] for e in errors],
                   label='Previous composite', color='#6a7485', linestyle='--')
    axes[0].set(title='Gain error: trine and echo overlap',
                xlabel='Common fractional gain error', ylabel='Average logical infidelity')
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.2, which='both')
    delta = np.geomspace(1e-7, 3e-4, 120)
    tc = report['resource']['composite_duration']
    d4 = report['quartic_joint_infidelity']['0,4']
    axes[1].loglog(delta, (tc*delta)**2/6, '--', color='#6a7485',
                   label='Previous composite: leading quadratic')
    axes[1].loglog(delta, d4*delta**4, color='#bc5a3c', label='Echo: leading quartic')
    witness = next(x for x in report['witnesses'] if x['gain'] == 0)
    axes[1].scatter([witness['reference_offset']], [witness['actual_infidelity']],
                    color='#bc5a3c', marker='o', label='Finite ODE witness')
    axes[1].set(title='Static reference detuning, zero gain error',
                xlabel=r'Reference offset $\Delta/a$', ylabel='Average logical infidelity')
    axes[1].legend(frameon=False, fontsize=9)
    axes[1].grid(alpha=0.2, which='both')
    rates = np.geomspace(1e-11, 1e-5, 140)
    colors = dict(composite='#6a7485', trine='#146c94', echo='#bc5a3c')
    for name, color in colors.items():
        offset = report['example']['metrics'][name]['average_infidelity']
        exposure = report['resource'][name+'_average_exposure']
        axes[2].loglog(rates, offset+rates*exposure, color=color, label=name.capitalize())
    axes[2].set(title='Choose using coherent error and loss',
                xlabel=r'Bright loss rate $\kappa/a$', ylabel='Average logical infidelity')
    axes[2].legend(frameon=False)
    axes[2].grid(alpha=0.2, which='both')
    axes[2].text(0.04, 0.08, 'Gain error: 0.1%\nReference offset: 0.0001 a\n'
                 'Leading absorbing-loss budget',
                 transform=axes[2].transAxes, fontsize=9)
    fig.suptitle('Reference-access reflection: phase cancellation with explicit duration and loss')
    fig.savefig(path, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json')
    parser.add_argument('--plot')
    args = parser.parse_args()
    report = run_checks()
    output = json.dumps(report, indent=2)
    if args.json:
        with open(args.json, 'w', encoding='utf-8') as target:
            target.write(output+'\n')
    else:
        print(output)
    if args.plot:
        plot(args.plot, report)

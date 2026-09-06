#!/usr/bin/env python3
"""Symmetry-reduced complete response for finite reflection words."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))  # repo root on sys.path (reorg 2026-09-06)
from pathlib import Path
import json
import numpy as np

from rlq.reflection_loop_detuning_order import (
    I4, primitive_noise_jet, expand_palindrome, reference_jet,
    arc_fourier, frame_coefficients, integrate_exponential_polynomial,
)
from rlq.reflection_loop_reference_echo import ended_composite, rotation
from rlq.reflection_loop_dynamics import J, GS
from math import factorial

CONTROL_PATH = Path(__file__).parent/'docs/receipts/reflection-loop-short-controls.json'


def parameters_batch(points):
    """Antisymmetric beta and symmetric stretches, one point per batch row."""
    points = np.atleast_2d(points)
    half = (points.shape[1]-1)//2
    beta = np.concatenate((points[:, :half], np.zeros((len(points), 1)),
                           -points[:, :half][:, ::-1]), axis=1)
    stretch = np.concatenate((points[:, half:2*half], points[:, -1:],
                              points[:, half:2*half][:, ::-1]), axis=1)
    signs = (-1.)**np.arange(2*half+1)
    alpha = 2*np.cumsum(signs*beta, axis=1)-signs*beta
    return beta, alpha, stretch, signs


def response_batch(points, n=1):
    """All nominal, linear and quadratic reference coefficients."""
    beta, alpha, stretch, signs = parameters_batch(points)
    co, sn = np.cos(alpha), np.sin(alpha)
    coeff = np.stack((sn*sn, co*co, sn*co), axis=-1)
    u, first, second = primitive_noise_jet(n)
    linear = np.einsum('bna,aij->bnij', coeff, first)
    quadratic = np.einsum('bna,bnc,acij->bnij', coeff, coeff, second)
    rot = np.broadcast_to(I4.real, (*alpha.shape, 4, 4)).copy()
    rot[..., 2, 2], rot[..., 2, 3] = co, -sn
    rot[..., 3, 2], rot[..., 3, 3] = sn, co
    result = np.zeros((len(beta), 3, 4, 4), complex)
    result[:, 0] = I4
    for j, sign in enumerate(signs):
        raw = np.stack((np.broadcast_to(u, linear[:, j].shape),
                        linear[:, j], quadratic[:, j]), axis=1)
        if sign < 0:
            raw = raw.conj().swapaxes(-1, -2)
        r = rot[:, j, None]
        raw = r@raw@r.swapaxes(-1, -2)
        raw *= (sign*stretch[:, j, None, None, None])**np.arange(3)[None, :, None, None]
        result = np.stack([sum(raw[:, k]@result[:, q-k] for k in range(q+1))
                           for q in range(3)], axis=1)
    return result


def real_leakage(quadratic):
    """Real coordinates when the lower code response has been corrected."""
    expj = np.array([[.5, 1j*np.sqrt(3)/2], [1j*np.sqrt(3)/2, .5]])
    return expj@quadratic[..., :2, 2:]@np.diag([1., -1j])


def equations_batch(points, n=1, gain_target=None, weight_target=None):
    beta, alpha, stretch, signs = parameters_batch(points)
    number = len(signs)
    c1, c2 = np.cos(beta), np.cos(2*beta)
    terms = [np.sum(stretch*c1, axis=1), np.sum(stretch*c2, axis=1)]
    for c in (c1, c2):
        for trig in (np.cos, np.sin):
            terms.append(np.sum(stretch*c*trig(2*alpha), axis=1))
    terms += [np.sum(stretch*np.sin(2*alpha), axis=1), c2@signs]
    circular = np.stack(terms, axis=1)/number
    second = response_batch(points, n)[:, 2]
    duration = sum(s.duration for s in ended_composite(n))
    logical = np.stack((-second[:, 2, 3].real,
                        (-second[:, 2, 2]-second[:, 3, 3]).imag/2), axis=1)
    leakage = real_leakage(second).real.reshape(-1, 4)
    result = np.concatenate((circular, logical/duration**2, leakage/duration**2), axis=1)
    if gain_target is not None:
        result = np.c_[result, (c1@signs-gain_target)/number]
    if weight_target is not None:
        result = np.c_[result, (stretch.sum(axis=1)-weight_target)/number]
    return result


def jacobian(function, point, step=1e-5):
    eye = np.eye(len(point))*step
    return ((function(point+eye)-function(point-eye))/(2*step)).T


def powers_through(order):
    return [(degree-q, q) for degree in range(order+1) for q in range(degree+1)]


def multiply_jets(left, right, powers):
    index = {power: j for j, power in enumerate(powers)}
    result = np.zeros_like(left)
    for k, (a, b) in enumerate(powers):
        for p in range(a+1):
            for q in range(b+1):
                result[..., k, :, :] += left[..., index[p, q], :, :]@right[..., index[a-p, b-q], :, :]
    return result


def arc_mixed_jet(arc, noises, n=1, order=4):
    """Exact finite sums for joint gain/reference coefficients at every order."""
    aa, speed = arc_fourier(arc, n)
    reference_modes = {}
    for p, ap in aa.items():
        for q, aq in aa.items():
            k = q-p
            reference_modes[k] = reference_modes.get(k, np.zeros_like(noises))+ap.conj().T@noises@aq
    jmat = np.zeros((4, 4), complex)
    jmat[:3, :3] = J
    khat = np.zeros((4, 4), complex)
    khat[:3, :3] = (J+speed*GS[arc])/(4*n*speed)
    projectors = {0: I4-khat@khat, 1: (khat@khat+khat)/2, -1: (khat@khat-khat)/2}
    f0, fc, _ = frame_coefficients(arc)
    initial_frame = f0+fc
    gain_modes = {}
    for p, pp in projectors.items():
        for q, pq in projectors.items():
            k = 4*n*(p-q)
            gain_modes[k] = gain_modes.get(k, np.zeros((4, 4), complex))+initial_frame@pp@jmat@pq@initial_frame.T
    modes = (gain_modes, reference_modes)
    endpoint = sum(a*(1j)**(k % 4) for k, a in aa.items())
    powers = powers_through(order)
    series = {(0, 0): {(0, 0): np.broadcast_to(I4, noises.shape).copy()}}
    result = np.zeros((len(noises), len(powers), 4, 4), complex)
    result[:, 0] = endpoint
    for idx, (p, q) in enumerate(powers[1:], 1):
        product = {}
        for channel, lower in ((0, (p-1, q)), (1, (p, q-1))):
            if min(lower) < 0:
                continue
            for k, mk in modes[channel].items():
                for (frequency, degree), value in series[lower].items():
                    key = (k+frequency, degree)
                    product[key] = product.get(key, np.zeros_like(noises))-1j*mk@value/speed
        terms = integrate_exponential_polynomial(product)
        series[p, q] = terms
        result[:, idx] = endpoint@sum(c*(np.pi/2)**degree*(1j)**(k % 4)
                                      for (k, degree), c in terms.items())
    return result


def mixed_jet(parameter_tuple, n=1, order=4, steps=None):
    """Fourier calculation, or independent laboratory ODE when steps is given."""
    _, angles, stretches, signs = parameter_tuple
    rotations = np.array([rotation(a) for a in angles])
    dr = np.diag([0., 0., 0., 1.])
    noises = rotations.transpose(0, 2, 1)@dr@rotations
    powers = powers_through(order)
    index = {p: j for j, p in enumerate(powers)}
    result = np.zeros((len(angles), len(powers), 4, 4), complex)
    result[:, 0] = I4
    cache = {}
    jmat = np.zeros((4, 4), complex)
    jmat[:3, :3] = J
    for stage in ended_composite(n):
        if steps is not None:
            count = max(8, int(np.ceil(stage.duration*steps)))
            dt = stage.duration/count

            def rhs(t, value):
                h = stage.h(t)
                out = -1j*h@value
                for k, (p, q) in enumerate(powers):
                    if p:
                        out[:, k] -= 1j*h@value[:, index[p-1, q]]
                    if q:
                        out[:, k] -= 1j*noises@value[:, index[p, q-1]]
                return out

            for j in range(count):
                t = j*dt
                k1 = rhs(t, result)
                k2 = rhs(t+dt/2, result+dt*k1/2)
                k3 = rhs(t+dt/2, result+dt*k2/2)
                k4 = rhs(t+dt, result+dt*k3)
                result += dt*(k1+2*k2+2*k3+k4)/6
            continue
        if stage.arc >= 0:
            if stage.arc not in cache:
                cache[stage.arc] = arc_mixed_jet(stage.arc, noises, n, order)
            block = cache[stage.arc]
            if stage.inverse:
                block = np.array([(-1.)**q for p, q in powers])[None, :, None, None]*block.conj().swapaxes(-1, -2)
        else:
            block = np.zeros_like(result)
            u = stage.u(stage.duration)
            area = (stage.q_start+stage.q_end)*stage.duration/2
            block[:, 0] = u
            for k, (p, q) in enumerate(powers[1:], 1):
                if p and q:
                    continue  # J and the rotated reference have orthogonal supports.
                if p:
                    block[:, k] = u@np.linalg.matrix_power(jmat, p)*(-1j*area)**p/factorial(p)
                else:
                    block[:, k] = u@noises*(-1j*stage.duration)**q/factorial(q)
        result = multiply_jets(block, result, powers)
    total = np.zeros_like(result[0])
    total[0] = I4
    for j, (r, stretch, sign) in enumerate(zip(rotations, stretches, signs)):
        block = result[j]
        if sign < 0:
            block = block.conj().swapaxes(-1, -2)
        block = r@block@r.T
        block *= np.array([(sign*stretch)**q for p, q in powers])[:, None, None]
        total = multiply_jets(block, total, powers)
    return {p: total[k] for k, p in enumerate(powers)}


def controls(n=1):
    record = json.loads(CONTROL_PATH.read_text())[str(n)]
    return expand_palindrome(np.array(record['center'], float))


def word_from_endpoint(u, parameter_tuple):
    out = I4.copy()
    for angle, sign in zip(parameter_tuple[1], parameter_tuple[3]):
        r = rotation(angle)
        out = r@(u if sign > 0 else u.conj().T)@r.T@out
    return out


def generic_gain_probe(parameter_tuple, seed):
    """Arbitrary bright dynamics and active-column error, independent of arcs."""
    from rlq.reflection_loop_compression import bad_response
    from rlq.reflection_loop_reference_echo import metrics
    rng = np.random.default_rng(seed)

    def exp_anti(a):
        values, vectors = np.linalg.eigh(1j*a)
        return (vectors*np.exp(-1j*values))@vectors.conj().T

    raw = rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2))
    bright = exp_anti((raw-raw.conj().T)/2)
    raw = rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2))
    moving = np.zeros((4, 4), complex)
    moving[:2, :2] = (raw-raw.conj().T)/4
    generator = np.zeros((4, 4), complex)
    column = rng.normal(size=2)+1j*rng.normal(size=2)
    generator[:2, 2], generator[2, :2] = column, -column.conj()
    generator[2, 2] = 1j*rng.normal()
    generator /= np.linalg.norm(generator)
    base = np.diag([1., 1., -1., 1.]).astype(complex)
    base[:2, :2] = bright
    derivative = moving@base+base@generator
    value, first = I4.copy(), np.zeros((4, 4), complex)
    for angle, sign in zip(parameter_tuple[1], parameter_tuple[3]):
        r = rotation(angle)
        u = r@(base if sign > 0 else base.conj().T)@r.T
        d = r@(derivative if sign > 0 else derivative.conj().T)@r.T
        first, value = d@value+u@first, u@value
    errors = []
    for rho in (.004, .002):
        u = exp_anti(rho*moving)@base@exp_anti(rho*generator)
        errors.append(metrics(word_from_endpoint(u, parameter_tuple))['average_infidelity'])
    return np.linalg.norm(bad_response(value.conj().T@first)), errors[1]/errors[0]


def run_checks():
    from rlq.reflection_loop_detuning_certificate import certify
    from rlq.reflection_loop_detuning_order import (
        higher_reference_jet, gain_word, finite_words, stages_from_parameters,
        time_moments,
    )
    from rlq.reflection_loop_compression import bad_response, packed_matrix, noise_generators, temporal_response
    from rlq.reflection_loop_reference_echo import metrics, response_integrals
    from rlq.reflection_loop_composite import coefficients
    records = json.loads(CONTROL_PATH.read_text())
    report = dict(checks={}, certificates={}, cases={})
    checks = report['checks']

    def close(name, value, target=0., tolerance=1e-9):
        error = float(np.max(abs(np.asarray(value)-target)))
        assert error <= tolerance, (name, error, tolerance)
        checks[name] = dict(error=error, tolerance=tolerance)

    def bound(name, value, maximum):
        assert value <= maximum, (name, value, maximum)
        checks[name] = dict(value=float(value), maximum=float(maximum))

    for seed in (7109, 9124, 543):
        first, ratio = generic_gain_probe(controls(1), seed)
        close(f'arbitrary_bright_dynamics_{seed}_first_response', first, tolerance=2e-12)
        close(f'arbitrary_bright_dynamics_{seed}_doubled_error_order', ratio, 1/16, .002)
    frequencies = np.r_[0., np.geomspace(1e-7, .002, 110)]
    for key, record in sorted(records.items(), key=lambda item: int(item[0])):
        n = int(key)
        report['certificates'][key] = certify(record)
        par = controls(n)
        x = np.array(record['center'], float)
        b, a, w, signs = par
        close(f'n{n}_complete_root_system', equations_batch(x, n, gain_target=0), tolerance=5e-12)
        for harmonic in (1, 2):
            close(f'n{n}_signed_moment_{harmonic}', signs@np.exp(1j*harmonic*b), tolerance=5e-13)
        stages = stages_from_parameters(par, n)
        response = response_integrals(stages, noise_generators())
        k = response['endpoint']
        for j, moment in enumerate(response['moments']):
            close(f'n{n}_static_noise_generator_{j}', bad_response(moment), tolerance=1e-8)
        exposure = w.sum()*3*coefficients(n)['exposure']/2
        close(f'n{n}_erasure_on_whole_code', response['loss'][:, 2:],
              exposure*I4[:, 2:], 3e-9)
        close(f'n{n}_zero_endpoints', max(np.linalg.norm(stages[0].h(0)),
              np.linalg.norm(stages[-1].h(stages[-1].duration))), tolerance=1e-12)
        close(f'n{n}_physical_joins', max(np.linalg.norm(l.h(l.duration)-r.h(0))
              for l, r in zip(stages, stages[1:])), tolerance=3e-12)
        bound(f'n{n}_echo_duration_budget', w.sum(), 24.)
        jet = mixed_jet(par, n, order=4)
        close(f'n{n}_full_endpoint', jet[0, 0], k, 5e-12)
        old = higher_reference_jet(par, n, order=4)
        for q in range(5):
            close(f'n{n}_independent_reference_order_{q}', (jet[0, q]-old[q])/(1+abs(old[q])), tolerance=2e-10)
        for power in ((1, 0), (0, 1), (2, 0), (1, 1), (0, 2), (3, 0)):
            close(f'n{n}_complete_code_coefficient_{power}', bad_response(k.conj().T@jet[power]), tolerance=3e-7)
        cubic_powers = [(2, 1), (1, 2), (0, 3)]
        vectors = np.array([bad_response(k.conj().T@jet[p]) for p in cubic_powers])
        gram = vectors.conj()@vectors.T
        close(f'n{n}_cubic_Gram_real_structure', gram.imag/(1+abs(gram)), tolerance=3e-9)
        real_form = real_leakage(jet[0, 3])
        close(f'n{n}_leading_leakage_real_form', real_form.imag/(1+abs(real_form)), tolerance=3e-8)
        gain8 = float(np.linalg.norm(bad_response(k.conj().T@jet[4, 0]))**2)
        delta6 = float(gram[2, 2].real)
        gain_size = .002
        gain_finite = [metrics(gain_word(s*gain_size, n, par))['average_infidelity'] for s in (1, -1)]
        close(f'n{n}_eighth_order_gain_finite_witness', np.mean(gain_finite)/(gain8*gain_size**8), 1., .015)
        finite = finite_words([1e-4, 1e-3], error=.001, parameter_tuple=par, n=n, steps=240)
        original, larger = [metrics(u) for u in finite]
        finer = finite_words([1e-4], error=.001, parameter_tuple=par, n=n, steps=480)[0]
        close(f'n{n}_finite_error_step_doubling', metrics(finer)['average_infidelity']/original['average_infidelity'],
              1., 5e-5)
        bound(f'n{n}_finite_unitarity', max(np.linalg.norm(u.conj().T@u-I4) for u in finite), 2e-9)
        negative = finite_words([-1e-4], error=-.001, parameter_tuple=par, n=n, steps=240)[0]
        monomials = np.array([.001**p*.0001**q for p, q in cubic_powers])
        prediction = float((monomials@gram@monomials).real)
        close(f'n{n}_joint_sextic_finite_witness', (original['average_infidelity']
              +metrics(negative)['average_infidelity'])/(2*prediction), 1., .015)
        tm = time_moments(stages)
        close(f'n{n}_constant_and_linear_gain_drift', tm[:2, 0], tolerance=1e-7)
        close(f'n{n}_constant_reference', tm[0, 1], tolerance=1e-8)
        tf = temporal_response(stages, frequencies)
        gain_frequency4 = float(np.linalg.norm(tm[2, 0])**2/4)
        reference_frequency2 = float(np.linalg.norm(tm[1, 1])**2)
        close(f'n{n}_reference_frequency_order', (tf['gram'][1, 1, 1].real
              +tf['negative_gram'][1, 1, 1].real)/(2*frequencies[1]**2*reference_frequency2), 1., .001)
        report['cases'][f'short_n{n}'] = dict(n=n, block_count=len(b), physical_stages=len(stages),
            weight=float(w.sum()), duration=sum(s.duration for s in stages), mean_exposure=float(exposure),
            gain_order=8, gain_coefficient=gain8, reference_order=6, reference_coefficient=delta6,
            cubic_powers=cubic_powers, cubic_gram=packed_matrix(gram),
            cubic_vectors=packed_matrix(vectors), cubic_gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
            response_coefficients={f'{p},{q}': packed_matrix(u) for (p, q), u in jet.items()},
            original=original, detuning_dominated=larger, joint_leading_prediction=prediction,
            gain_frequency_order=4, gain_frequency_coefficient=gain_frequency4,
            reference_frequency_order=2, reference_frequency_coefficient=reference_frequency2,
            temporal_frequency=frequencies.tolist(), temporal_gram=packed_matrix(tf['gram']),
            negative_temporal_gram=packed_matrix(tf['negative_gram']),
            temporal_reference=((tf['gram'][:, 1, 1].real+tf['negative_gram'][:, 1, 1].real)/2).tolist())
        if n == 1:
            numerical = mixed_jet(par, n, order=4, steps=120)
            for power in powers_through(4):
                close(f'independent_ODE_joint_coefficient_{power}', (numerical[power]-jet[power])/(1+abs(jet[power])), tolerance=5e-7)
            rate = 1e-8
            lossless, damped = finite_words([0., 0.], loss=[0., rate], parameter_tuple=par, n=n, steps=240)
            close('state_independent_loss_finite_slope',
                  (metrics(damped, True)['average_infidelity']-metrics(lossless, True)['average_infidelity'])/rate,
                  exposure, .01)
        print(f'n={n}: full joint correction, finite errors, temporal response and certificate checked', flush=True)
    old = json.loads((CONTROL_PATH.parent/'reflection-loop-detuning-checks.json').read_text())
    compression = json.loads((CONTROL_PATH.parent/'reflection-loop-compression-checks.json').read_text())
    for name, case in old['cases'].items():
        row = {key: case[key] for key in ('n', 'duration', 'mean_exposure', 'original', 'reference_order',
               'reference_coefficient', 'temporal_frequency', 'temporal_reference')}
        if name.startswith('thirteen'):
            gain = case['leading_joint']['gain4']
        else:
            parent = name.replace('raised_', 'compressed_')
            gain = compression['cases'][parent]['response']['terms']['4,0']
        row.update(gain_order=4, gain_coefficient=gain)
        report['cases'][name] = row
    report['budget_examples'] = []
    for rate in (1e-15, 1e-13, 1e-10, 1e-9):
        options = [(case['original']['average_infidelity']+rate*case['mean_exposure'], name)
                   for name, case in report['cases'].items() if case['duration'] <= 3000]
        best = min(options)
        previous = min(item for item in options if not item[1].startswith('short_'))
        report['budget_examples'].append(dict(duration_limit=3000, loss_rate=rate,
            gain=.001, delta=.0001, selected=best[1], predicted_infidelity=best[0],
            previous_best=previous[1], previous_infidelity=previous[0], improvement=previous[0]/best[0]))
    report['check_count'] = len(checks)
    report['scope'] = ('Certified full second-order encoded gain/reference-detuning correction at n=1,2,3; '
                       'pure gain begins at eighth-order infidelity and joint errors at sixth order. '
                       'Other generators retain first-order correction; affine reference drift remains.')
    return report


def plot(path, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = dict(compressed='#858c96', echo='#216d8a', thirteen='#b95e2a',
                  raised='#65549c', short='#087f5b')
    labels = dict(compressed='11 blocks', echo='Echo', thirteen='13 blocks',
                  raised='Raised word', short='Complete code correction')
    markers = dict(compressed='x', echo='o', thirteen='s', raised='D', short='^')
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.5), layout='constrained')
    for family, color in colors.items():
        for n in (1, 2, 3):
            row = report['cases'][f'{family}_n{n}']
            axes[0, 0].scatter(row['duration'], row['mean_exposure'], color=color,
                               marker=markers[family], s=45, label=labels[family] if n == 1 else None)
            offset = (4, -12) if family == 'short' else (4, 5)
            axes[0, 0].annotate(str(n), (row['duration'], row['mean_exposure']),
                               xytext=offset, textcoords='offset points', fontsize=8)
        row = report['cases'][family+'_n1']
        if family != 'raised':  # Same gain-only endpoint as the eleven-block parent.
            eps = np.geomspace(1e-4, .02, 120)
            axes[0, 1].loglog(eps, row['gain_coefficient']*eps**row['gain_order'],
                              color=color, label=labels[family])
        freq = np.asarray(row['temporal_frequency'])[1:]
        axes[1, 1].loglog(freq, row['temporal_reference'][1:], color=color, label=labels[family])
        rates = np.geomspace(1e-18, 1e-8, 160)
        rows = [r for name, r in report['cases'].items()
                if name.startswith(family+'_') and r['duration'] <= 3000]
        errors = np.min([r['original']['average_infidelity']+rates*r['mean_exposure'] for r in rows], axis=0)
        axes[1, 0].loglog(rates, errors, color=color, label=labels[family])
    axes[0, 0].set(xlabel=r'Duration $aT$', ylabel=r'Mean bright exposure $a\bar X$',
                   title='Resource costs; numbers denote return index n', xlim=(0, 8000), ylim=(0, 150))
    axes[0, 0].legend(frameon=False, fontsize=8, loc='upper right')
    axes[0, 1].set(xlabel=r'Common gain error $\epsilon$', ylabel='Leading average infidelity',
                   title='Pure gain at n=1: fourth versus eighth order')
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[1, 0].set(xlabel=r'Absorbing loss rate $\kappa/a$', ylabel='Finite coherent error + leading loss',
                   title=r'Choose among n = 1–3 with $aT\leq3000$')
    axes[1, 0].text(.04, .91, r'$\epsilon=10^{-3},\ \Delta/a=10^{-4}$',
                    transform=axes[1, 0].transAxes, fontsize=9)
    axes[1, 1].set(xlabel=r'Noise angular frequency $\omega/a$',
                   ylabel=r'Reference response $a^2F_\Delta(\omega)$',
                   title='Reference drift at n=1: the remaining tradeoff')
    for ax in axes.flat:
        ax.grid(which='both', alpha=.2)
    fig.suptitle('Complete encoded second-order correction within the previous echo budget')
    fig.savefig(path, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json')
    parser.add_argument('--plot')
    args = parser.parse_args()
    report = run_checks()
    output = json.dumps(report, indent=2)+'\n'
    if args.json:
        Path(args.json).write_text(output)
    else:
        print(output)
    if args.plot:
        plot(args.plot, report)

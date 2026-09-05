#!/usr/bin/env python3
"""Finite reflection compression by circular response moments.

Native inputs: REALFIBER, CP-1--4 and RE-1--6. Eleven alternating forward/inverse
composites replace the nested echo. All rotations compile real Hamiltonians;
every finite ramp is retained. No instantaneous control is applied.
"""
import argparse
import json

import mpmath as mp
import numpy as np
import sympy as sp

from reflection_loop_composite import P, Q, REFERENCE, composite, coefficients, ode
from reflection_loop_reference_echo import (
    EYE, POWERS, FramedStage, rotation, ended_composite, echo_stages,
    response_integrals, jet_constant, jet_multiply, quartic_response, metrics,
    block_jets, combined_jets, finite_gate, exact_gain, predict,
)

# A nearby exact root is certified below; these are not fitted ODE parameters.
# Five half-angles and three nontrivial stretches; the other stretches are one.
ROOT = (
    '-0.38158500240002994580545275152219225651079415259483318526603341',
    '-1.3390135686404914498955974417231346653992714035433985373306537',
    '-2.528015782758935599895692564796985622995847467963770734275491',
    '2.6041829582543396363327631476340590657423716411209034778798106',
    '1.4660588420813267016300986487435299973678576180910819006158223',
    '1.2644794881200301487515034671347902897769196615102302024761494',
    '1.0639750240984300163758432788794431941714227686806189337616642',
    '1.1096018909075060192151088261483118285402907038672034119411167',
)


def parameters():
    x = np.array(ROOT, dtype=float)
    beta = np.r_[x[:5], 0., -x[:5][::-1]]
    stretch = np.r_[1., x[5:], 1., 1., 1., x[5:][::-1], 1.]
    sign = (-1.)**np.arange(11)
    alpha = 2*np.cumsum(sign*beta)-sign*beta
    return beta, alpha, stretch, sign


def moment_residuals():
    b, a, w, sig = parameters()
    phases = (b, 2*b, b+2*a, b-2*a, 2*b+2*a, 2*b-2*a)
    return dict(circular=np.exp(1j*np.array(phases))@w,
                imaginary_overlap=w@np.sin(2*a),
                gain_phase=sig@np.exp(2j*b), target_angle=sig@b)


def root_expressions():
    x = sp.symbols('x0:8', real=True)
    b = list(x[:5])+[sp.S.Zero]+[-v for v in x[:5][::-1]]
    half = [sp.S.One]+list(x[5:])+[sp.S.One]
    w = half+[sp.S.One]+half[::-1]
    sig = [(-1)**j for j in range(11)]
    phi, a = sp.S.Zero, []
    for v, s in zip(b, sig):
        a.append(phi+s*v)
        phi += 2*s*v
    c1, c2 = ([sp.cos(k*v) for v in b] for k in (1, 2))
    f = [sum(v*c for v, c in zip(w, c1)),
         sum(v*c for v, c in zip(w, c2))]
    for cb in (c1, c2):
        for trig in (sp.cos, sp.sin):
            f.append(sum(v*c*trig(2*t) for v, c, t in zip(w, cb, a)))
    f.extend((sum(v*sp.sin(2*t) for v, t in zip(w, a)),
              sum(s*c for s, c in zip(sig, c2))))
    return x, sp.Matrix(f)


def certify_root():
    """Interval contraction certificate, not a floating residual assertion.

    T(x)=x-A f(x) maps an interval cube into itself and contracts there.
    Every arithmetic/trigonometric operation in the bounds uses mp.iv.
    """
    mp.mp.dps = 75
    mp.iv.dps = 70
    x, f = root_expressions()
    j = f.jacobian(x)
    fj = sp.lambdify([x], list(f), 'mpmath')
    jj = sp.lambdify([x], j.tolist(), 'mpmath')
    ivmod = [{'sin': mp.iv.sin, 'cos': mp.iv.cos}]
    fi = sp.lambdify([x], list(f), ivmod)
    ji = sp.lambdify([x], j.tolist(), ivmod)
    center = [mp.mpf(v) for v in ROOT]
    inv = mp.inverse(mp.matrix(jj(center)))
    radius = mp.mpf('1e-40')
    point = [mp.iv.mpf(v) for v in ROOT]
    cube = [v+mp.iv.mpf(['-1e-40', '1e-40']) for v in point]
    ai = [[mp.iv.mpf(mp.nstr(inv[r, c], 73)) for c in range(8)] for r in range(8)]
    jf, ff = ji(cube), fi(point)

    def upper_abs(v):
        return mp.mpf(abs(v).b)

    contraction = max(upper_abs(sum(abs((int(r == c))
        -sum(ai[r][k]*jf[k][c] for k in range(8))) for c in range(8))) for r in range(8))
    displacement = max(upper_abs(sum(ai[r][k]*ff[k] for k in range(8))) for r in range(8))
    assert contraction < 1
    assert displacement < (1-contraction)*radius
    assert min(center[5:])-radius > 1
    return dict(radius=mp.nstr(radius, 8), contraction_bound=mp.nstr(contraction, 12),
                newton_displacement_bound=mp.nstr(displacement, 12),
                root_exists_and_is_unique_in_cube=True,
                smallest_nontrivial_stretch_lower_bound=mp.nstr(min(center[5:])-radius, 18),
                center_residual=str(max(abs(v) for v in fj(center))))


def compressed_stages(n=1, gap=1., slew=None):
    _, a, w, signs = parameters()
    base = ended_composite(n, gap, slew)
    out = []
    for alpha, stretch, sign in zip(a, w, signs):
        reverse = sign < 0
        stages = list(reversed(base)) if reverse else base
        out.extend(FramedStage(s, rotation(alpha), stretch, reverse) for s in stages)
    return out


def compressed_gain(error=0., n=1):
    _, a, _, signs = parameters()
    c = composite(error, n)
    out = EYE.copy()
    for alpha, sign in zip(a, signs):
        rot = rotation(alpha)
        u = c if sign > 0 else c.conj().T
        out = rot @ u @ rot.T @ out
    return out


def compressed_jets(steps=100, n=1, parameter_tuple=None):
    """Integrate rotated laboratory perturbations together on one base waveform.

    Each physical reference projector is pulled back by its compiled rotation.
    Time stretching and the detuning sign inside inverse adjoints are explicit.
    """
    _, angles, stretches, signs = parameters() if parameter_tuple is None else parameter_tuple
    rotations = np.array([rotation(a) for a in angles])
    noises = rotations.transpose(0, 2, 1) @ REFERENCE @ rotations
    state = np.zeros((len(angles), len(POWERS), 4, 4), complex)
    state[:, 0] = EYE
    position = {power: j for j, power in enumerate(POWERS)}
    for stage in ended_composite(n):
        count = max(8, int(np.ceil(stage.duration*steps)))
        dt = stage.duration/count

        def rhs(t, values):
            h = stage.h(t)
            result = -1j*h @ values
            for j, (p, q) in enumerate(POWERS):
                if p:
                    result[:, j] -= 1j*h @ values[:, position[p-1, q]]
                if q:
                    result[:, j] -= 1j*noises @ values[:, position[p, q-1]]
            return result

        for j in range(count):
            t = j*dt
            k1 = rhs(t, state)
            k2 = rhs(t+dt/2, state+dt*k1/2)
            k3 = rhs(t+dt/2, state+dt*k2/2)
            k4 = rhs(t+dt, state+dt*k3)
            state += dt*(k1+2*k2+2*k3+k4)/6
    out = jet_constant(EYE)
    for j, (rot, stretch, sign) in enumerate(zip(rotations, stretches, signs)):
        block = {}
        for k, (p, q) in enumerate(POWERS):
            value = state[j, k]
            if sign < 0:
                value = value.conj().T
            block[p, q] = (stretch*sign)**q * (rot @ value @ rot.T)
        out = jet_multiply(block, out)
    return out


def noise_generators():
    vs = [np.diag([int(i == j) for i in range(4)]) for j in range(4)]
    for entry in (1., 1j):
        v = np.zeros((4, 4), complex)
        v[0, 1], v[1, 0] = entry, np.conjugate(entry)
        vs.append(v)
    return vs


def bad_response(m):
    """Linear error amplitude with the Haar-infidelity inner product."""
    p = m[..., 2:, 2:]
    trace = np.trace(p, axis1=-2, axis2=-1)
    tf = p-trace[..., None, None]*np.eye(2)/2
    return np.concatenate((m[..., :2, 2:].reshape(m.shape[:-2]+(4,))/np.sqrt(2),
                           tf.reshape(m.shape[:-2]+(4,))/np.sqrt(3)), axis=-1)


def temporal_response(stages, frequencies, nodes=48):
    """Gain/reference spectral Gram matrices and their first time moments.

    F_ab(omega)=<khat_a,khat_b>. The response is for real classical weak noise;
    no bath, Gaussian, or quasi-static assumption enters this second moment.
    """
    points, weights = np.polynomial.legendre.leggauss(nodes)
    moments = np.zeros((len(frequencies), 2, 8), complex)
    negative_moments = np.zeros_like(moments)
    first = np.zeros((2, 8), complex)
    second = np.zeros((2, 8), complex)
    dc = np.zeros((2, 8), complex)
    prefix, offset = EYE.copy(), 0.
    duration = sum(s.duration for s in stages)
    for stage in stages:
        for point, weight in zip(points, weights):
            local = stage.duration*(point+1)/2
            t = offset+local-duration/2
            wt = weight*stage.duration/2
            u = stage.u(local) @ prefix
            k = bad_response(u.conj().T @ np.array([stage.h(local), REFERENCE]) @ u)
            phase = np.exp(1j*np.asarray(frequencies)*t)[:, None, None]
            moments += wt*phase*k
            negative_moments += wt*phase.conj()*k
            first += wt*t*k
            second += wt*t*t*k
            dc += wt*k
        prefix = stage.u(stage.duration) @ prefix
        offset += stage.duration
    gram = np.einsum('wai,wbi->wab', moments.conj(), moments)
    negative_gram = np.einsum('wai,wbi->wab', negative_moments.conj(), negative_moments)
    first_gram = first.conj() @ first.T
    return dict(frequency=np.asarray(frequencies), gram=gram, negative_gram=negative_gram,
                amplitudes=moments, negative_amplitudes=negative_moments,
                low_frequency_gram=first_gram, second_time_moment=second,
                dc=dc, duration=duration)


def packed_matrix(m):
    return dict(real=np.asarray(m).real.tolist(), imag=np.asarray(m).imag.tolist())


def packed_response(r):
    return dict(terms=r['terms'], **{name: {f'{p},{q}': packed_matrix(v)
                for (p, q), v in r[name].items()} for name in ('leakage', 'logical')})


def run_checks():
    report = dict(root_certificate=certify_root(), checks={}, cases={})
    checks = report['checks']

    def close(name, value, expected=0., tolerance=1e-9):
        error = float(np.max(np.abs(np.asarray(value)-expected)))
        assert error <= tolerance, (name, error, tolerance)
        checks[name] = dict(error=error, tolerance=tolerance)

    def bounded(name, value, maximum):
        assert value <= maximum, (name, value, maximum)
        checks[name] = dict(value=float(value), maximum=float(maximum))

    beta, alpha, w, sig = parameters()
    for name, residual in moment_residuals().items():
        close('circular_conditions_'+name, residual, tolerance=2e-12)
    close('same_encoded_reflection', compressed_gain()[2:, 2:], np.diag([-1., 1.]))
    close('declared_bright_endpoint', compressed_gain()[:2, :2], composite()[:2, :2])
    # The root equations express cancellation for arbitrary local primitive
    # coefficients. Check them on independent arbitrary waveforms as well.
    rng = np.random.default_rng(922381)
    for trial in range(3):
        raw = rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3))
        local, _ = np.linalg.qr(raw)
        u = EYE.copy()
        u[:3, :3] = local
        summed = np.zeros((6, 4, 4), complex)
        for aa, bb, ww in zip(alpha, beta, w):
            a, b = rotation(aa), rotation(bb)
            vs = a.T @ np.array(noise_generators()) @ a
            summed += ww*(b @ u.conj().T @ vs @ u @ b.T)
        close(f'arbitrary_primitive_full_response_{trial}', bad_response(summed), tolerance=3e-12)

    frequencies = np.r_[0., np.geomspace(1e-7, .002, 121)]
    witness_gain, witness_delta = .001, .0001
    for n in (1, 2, 3):
        base_time = sum(s.duration for s in ended_composite(n))
        xc = 3*coefficients(n)['exposure']
        stages = compressed_stages(n)
        close(f'n{n}_zero_ends', max(np.linalg.norm(stages[0].h(0)),
              np.linalg.norm(stages[-1].h(stages[-1].duration))), tolerance=1e-12)
        close(f'n{n}_continuous_joins', max(np.linalg.norm(l.h(l.duration)-r.h(0))
              for l, r in zip(stages, stages[1:])), tolerance=1e-12)
        peak_h, peak_slew, spectral, diag = 0., 0., 0., 0.
        for stage in stages:
            for f in (.2, .5, .8):
                t = f*stage.duration
                h = stage.h(t)
                ev = np.linalg.eigvalsh(h)
                peak_h = max(peak_h, np.linalg.norm(h, 2))
                spectral = max(spectral, np.linalg.norm(ev[1:3]), abs(ev[0]+ev[3]))
                diag = max(diag, np.linalg.norm(np.diag(h)), np.linalg.norm(h.imag))
                dt = min(1e-5, stage.duration*1e-4)
                peak_slew = max(peak_slew, np.linalg.norm((stage.h(t+dt)-stage.h(t-dt))/(2*dt), 2))
        bounded(f'n{n}_amplitude_cap', peak_h, 1.+1e-12)
        bounded(f'n{n}_slew_cap', peak_slew, 10.+1e-8)
        close(f'n{n}_real_zero_diagonal_rank_two_family', max(spectral, diag), tolerance=1e-12)
        data = response_integrals(stages, noise_generators())
        close(f'n{n}_nominal_endpoint', data['endpoint'], composite(), 2e-12)
        bounded(f'n{n}_parallel_transport', data['parallel_error'], 1e-12)
        close(f'n{n}_gain_first_order', data['gain'] @ P, tolerance=3e-11)
        for j, moment in enumerate(data['moments']):
            close(f'n{n}_whole_code_static_noise_{j}', bad_response(moment), tolerance=2e-10)
        close(f'n{n}_state_independent_erasure', data['loss'] @ P, w.sum()*xc*P/2, 3e-10)
        # Independent Taylor integration versus exact moment quadrature.
        fine = compressed_jets(120, n)
        response = quartic_response(fine, target=composite())
        first = composite().conj().T @ fine[0, 1]
        close(f'n{n}_jet_matches_reference_quadrature',
              1j*first @ P, data['moments'][3] @ P, 5e-8)
        close(f'n{n}_mixed_gain_detuning_leakage', response['leakage'][1, 1], tolerance=8e-7)
        close(f'n{n}_mixed_gain_detuning_logical', response['logical'][1, 1], tolerance=3e-7)
        b = coefficients(n)['b']
        rn = 1-1/(16*n*n)
        v = np.array([sig@np.cos(beta), sig@np.sin(beta)])
        lc = np.pi**2*rn/12*np.array([-2*np.sqrt(3)+3j, np.sqrt(3)-6j])
        expected = np.vstack((np.outer(lc, v), 1j*np.sqrt(3)*b*b*np.eye(2)))
        close(f'n{n}_complete_quadratic_gain_column',
              (composite().conj().T @ fine[2, 0])[:, 2:], expected, 2e-7)
        gain_coefficient = coefficients(n)['composite_leakage']*np.dot(v, v)/2
        close(f'n{n}_quartic_gain_coefficient', response['terms']['4,0'], gain_coefficient, 1e-6)
        phase_big = metrics(compressed_gain(.004, n))['logical_traceless_norm']
        phase_small = metrics(compressed_gain(.002, n))['logical_traceless_norm']
        bounded(f'n{n}_relative_gain_phase_fourth_order', phase_small/phase_big, .11)
        actual = ode(stages, error=witness_gain, detuning=witness_delta, steps_per_time=120)
        bounded(f'n{n}_finite_ode_unitarity', np.linalg.norm(actual.conj().T @ actual-EYE), 5e-9)
        temporal = temporal_response(stages, frequencies)
        close(f'n{n}_dc_temporal_response', temporal['dc'], tolerance=3e-10)
        bounded(f'n{n}_linear_gain_drift_cancelled', abs(temporal['low_frequency_gram'][0, 0]), 1e-15)
        bounded(f'n{n}_spectral_gram_positive', -np.min(np.linalg.eigvalsh(temporal['gram'])), 1e-8)
        low_gain = float(np.vdot(temporal['second_time_moment'][0], temporal['second_time_moment'][0]).real/4)
        low_delta = float(temporal['low_frequency_gram'][1, 1].real)
        # The first nonzero asymptotic coefficients, not only a static zero.
        close(f'n{n}_gain_filter_fourth_order',
              temporal['gram'][1, 0, 0].real/(frequencies[1]**4*low_gain), 1., .003)
        close(f'n{n}_reference_filter_second_order',
              temporal['gram'][1, 1, 1].real/(frequencies[1]**2*low_delta), 1., .0001)
        row = dict(n=n, duration=w.sum()*base_time, mean_exposure=w.sum()*xc/2,
                   gain_metrics=metrics(compressed_gain(witness_gain, n)),
                   joint_metrics=metrics(actual), response=packed_response(response),
                   temporal=dict(frequency=frequencies.tolist(), gram=packed_matrix(temporal['gram']),
                                 negative_gram=packed_matrix(temporal['negative_gram']),
                                 gain_fourth_order_coefficient=low_gain,
                                 reference_second_order_coefficient=low_delta))
        report['cases'][f'compressed_n{n}'] = row
        # Replay the n parameter in the existing echo, without changing its word.
        old_a = block_jets(120, n)
        _, old_jet = combined_jets(old_a)
        old_response = quartic_response(old_jet)
        old_actual = finite_gate(witness_gain, witness_delta, steps=120, n=n)
        old_temporal = temporal_response(echo_stages(n), frequencies)
        bounded(f'n{n}_echo_linear_gain_drift_cancelled',
                abs(old_temporal['low_frequency_gram'][0, 0]), 1e-14)
        report['cases'][f'echo_n{n}'] = dict(
            n=n, duration=24*base_time, mean_exposure=12*xc,
            gain_metrics=metrics(exact_gain(witness_gain, n=n)), joint_metrics=metrics(old_actual),
            response=packed_response(old_response),
            temporal=dict(frequency=frequencies.tolist(), gram=packed_matrix(old_temporal['gram']),
              negative_gram=packed_matrix(old_temporal['negative_gram']),
              gain_fourth_order_coefficient=float(np.vdot(old_temporal['second_time_moment'][0],
                                                        old_temporal['second_time_moment'][0]).real/4),
              reference_second_order_coefficient=float(old_temporal['low_frequency_gram'][1, 1].real)))
        print(f'n={n}: full response, finite witness and temporal filters checked', flush=True)
        if n == 1:
            coarse = compressed_jets(60, n)
            cr = quartic_response(coarse, target=composite())
            for key in ('4,0', '2,2', '0,4'):
                close('joint_jet_convergence_'+key, cr['terms'][key]/response['terms'][key], 1., 2e-5)
            bounded('mixed_jet_step_doubling',
                    np.linalg.norm(response['leakage'][1, 1])/np.linalg.norm(cr['leakage'][1, 1]), .09)
            coarse_direct = ode(stages, error=witness_gain, detuning=witness_delta, steps_per_time=60)
            finer_direct = ode(stages, error=witness_gain, detuning=witness_delta, steps_per_time=240)
            bounded('laboratory_ode_step_doubling',
                    np.linalg.norm(actual-finer_direct)/np.linalg.norm(coarse_direct-finer_direct), .09)
            close('laboratory_ode_accuracy', actual, finer_direct, 3e-9)
            small = ode(stages, error=.0005, detuning=.00005, steps_per_time=120)
            relative = abs(metrics(small)['average_infidelity']/predict(response, .0005, .00005)-1)
            bounded('joint_quartic_finite_witness', relative, .035)
            rate = 1e-8
            damped = ode(stages, loss=rate, steps_per_time=120)
            clean = ode(stages, steps_per_time=120)
            slope = (metrics(damped, True)['average_infidelity']-metrics(clean, True)['average_infidelity'])/rate
            close('absorbing_loss_slope', slope, w.sum()*xc/2, .01)
            # Deliberately outside DC: a slow varying gain and detuning filter.
            tdense = temporal_response(stages, [1e-4, .01], nodes=96)
            tcoarse = temporal_response(stages, [1e-4, .01], nodes=48)
            close('spectral_quadrature_convergence', tcoarse['gram'], tdense['gram'], 1e-7)
            omega, amplitude = .001, 1e-6
            tdynamic = temporal_response(stages, [omega], nodes=64)
            sine = (tdynamic['amplitudes'][0, 1]-tdynamic['negative_amplitudes'][0, 1])/(2j)
            predicted_dynamic = amplitude**2*np.vdot(sine, sine).real
            duration = sum(s.duration for s in stages)
            dynamic = ode(stages, steps_per_time=120,
                          detuning_shape=lambda t: amplitude*np.sin(omega*(t-duration/2)))
            close('finite_time_varying_reference_witness',
                  metrics(dynamic)['average_infidelity']/predicted_dynamic, 1., .01)
    report['check_count'] = len(checks)
    report['parameters'] = dict(beta=beta.tolist(), alpha=alpha.tolist(), stretch=w.tolist(),
                                total_stretch=float(w.sum()), root=list(ROOT))
    report['witness'] = dict(gain=witness_gain, reference_offset=witness_delta)
    report['scope'] = ('Finite zero-ended real controls, common gain and quasi-static six-dimensional '
                       'noise class. Temporal Gram matrices describe weak classical noise. '
                       'No global optimality, arbitrary-noise correction or hardware advantage is assumed.')
    return report


def plot(path, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12.7, 8.6), layout='constrained')
    colors = ['#196c8b', '#b75b32', '#6b5797']
    rates = np.geomspace(1e-13, 1e-7, 180)
    deltas = np.geomspace(1e-6, 3e-4, 120)
    for n, color in zip((1, 2, 3), colors):
        for family, style, marker in (('echo', '--', 'o'), ('compressed', '-', 's')):
            row = report['cases'][f'{family}_n{n}']
            label = f'{"Echo" if family == "echo" else "11 blocks"}, n={n}'
            axes[0, 0].scatter(row['duration'], row['mean_exposure'], color=color, marker=marker,
                               facecolors='none' if family == 'echo' else color, s=65)
            axes[0, 0].annotate(label, (row['duration'], row['mean_exposure']),
                                xytext=(5, 7), textcoords='offset points', fontsize=8)
            terms = row['response']['terms']
            eps = report['witness']['gain']
            leading = terms['4,0']*eps**4+terms['2,2']*eps**2*deltas**2+terms['0,4']*deltas**4
            axes[0, 1].loglog(deltas, leading, style, color=color, label=label)
            axes[0, 1].scatter([report['witness']['reference_offset']],
                               [row['joint_metrics']['average_infidelity']], color=color,
                               marker=marker, s=22, facecolors='none' if family == 'echo' else color)
            axes[1, 0].loglog(rates, row['joint_metrics']['average_infidelity']
                              +rates*row['mean_exposure'], style, color=color, label=label)
            temporal = row['temporal']
            freq = np.array(temporal['frequency'])[1:]
            gram = (np.array(temporal['gram']['real'])+np.array(temporal['negative_gram']['real']))/2
            axes[1, 1].loglog(freq, gram[1:, 1, 1], style, color=color, label=label)
    axes[0, 0].set(xlabel=r'Duration $aT$', ylabel=r'Mean bright exposure $a\bar X$',
                    title='Both duration and exposure count', xlim=(0, 4700), ylim=(0, 90))
    axes[0, 1].set(xlabel=r'Static reference offset $\Delta/a$', ylabel='Average logical infidelity',
                    title='Joint quartic response; gain error 0.1%')
    axes[0, 1].legend(frameon=False, fontsize=8, ncol=2)
    axes[1, 0].set(xlabel=r'Absorbing loss rate $\kappa/a$', ylabel='Average logical infidelity',
                    title='Finite coherent witness + leading loss')
    axes[1, 0].text(.04, .9, r'$\epsilon=10^{-3},\;\Delta/a=10^{-4}$',
                    transform=axes[1, 0].transAxes, fontsize=9)
    axes[1, 1].set(xlabel=r'Noise angular frequency $\omega/a$',
                    ylabel=r'Reference response $a^2 F_\Delta(\omega)$',
                    title='Response to time-varying reference noise')
    for ax in axes.flat:
        ax.grid(alpha=.2, which='both')
    fig.suptitle('Reflection control: compression and the return-index tradeoff')
    fig.savefig(path, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json')
    parser.add_argument('--plot')
    parser.add_argument('--root-only', action='store_true')
    args = parser.parse_args()
    report = certify_root() if args.root_only else run_checks()
    output = json.dumps(report, indent=2)
    if args.json:
        with open(args.json, 'w', encoding='utf-8') as target:
            target.write(output+'\n')
    else:
        print(output)
    if args.plot:
        plot(args.plot, report)

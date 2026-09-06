#!/usr/bin/env python3
"""Native finite Fourier response and second-order detuning correction."""
import argparse
import json
from functools import lru_cache
from math import factorial
from pathlib import Path

import numpy as np

from reflection_loop_dynamics import J, GS, frame, embed
from reflection_loop_composite import composite, ode, coefficients
from reflection_loop_reference_echo import ended_composite, FramedStage, rotation, metrics
from reflection_loop_compression import (
    parameters, compressed_stages, bad_response, compressed_jets,
    temporal_response, packed_matrix, packed_response,
)
from reflection_loop_reference_echo import (
    response_integrals, quartic_response, echo_stages, finite_gate,
)
from reflection_loop_compression import noise_generators

I4 = np.eye(4, dtype=complex)
NOISE = np.zeros((3, 4, 4), complex)
NOISE[0, 2, 2] = NOISE[1, 3, 3] = 1
NOISE[2, 2, 3] = NOISE[2, 3, 2] = 1
CONTROL_PATH = Path(__file__).parent/'docs/reflection-loop-detuning-controls.json'


def integral_frequency(k):
    """Integral exp(i k u) du from 0 to pi/2; phases are exact fourth roots."""
    return np.pi/2 if k == 0 else ((1j)**(k % 4)-1)/(1j*k)


def ordered_frequency(k, l):
    """Integral exp(i k u) integral_0^u exp(i l v) dv du."""
    if l:
        return (integral_frequency(k+l)-integral_frequency(k))/(1j*l)
    if k == 0:
        return np.pi**2/8
    phase = (1j)**(k % 4)
    return np.pi*phase/(2j*k)+(phase-1)/(k*k)


def arc_fourier(arc, n=1):
    """U(t)=sum_p A_p exp(i p v t) at a native return, v=1/sqrt(16n²-1)."""
    constant, cosine, sine = frame_coefficients(arc)
    ff = {0: constant, 1: (cosine-1j*sine)/2, -1: (cosine+1j*sine)/2}
    v = 1/np.sqrt(16*n*n-1)
    khat = np.zeros((4, 4), complex)
    khat[:3, :3] = (J+v*GS[arc])/(4*n*v)
    projectors = {0: I4-khat@khat, 1: (khat@khat+khat)/2, -1: (khat@khat-khat)/2}
    out = {}
    f0 = constant+cosine
    for k, mat in ff.items():
        for m, proj in projectors.items():
            p = k-4*n*m
            out[p] = out.get(p, np.zeros((4, 4), complex))+mat@proj@f0.T
    return out, v


def frame_coefficients(arc):
    """Integer coefficients of the native frame F0+Fc cos(theta)+Fs sin(theta)."""
    constant = np.zeros((4, 4), complex)
    constant[3, 3] = 1
    cosine, sine = np.zeros_like(constant), np.zeros_like(constant)
    if arc == 0:
        constant[1, 1] = 1
        cosine[0, 0] = cosine[2, 2] = 1
        sine[0, 2], sine[2, 0] = -1, 1
    elif arc == 1:
        constant[2, 0] = 1
        cosine[0, 2], cosine[1, 1] = -1, 1
        sine[0, 1] = sine[1, 2] = 1
    else:
        constant[0, 1] = 1
        cosine[1, 2] = cosine[2, 0] = 1
        sine[1, 0], sine[2, 2] = 1, -1
    return constant, cosine, sine


@lru_cache(maxsize=None)
def arc_noise_jet(arc, n=1):
    aa, speed = arc_fourier(arc, n)
    modes = {}
    for p, ap in aa.items():
        for q, aq in aa.items():
            k = q-p
            modes[k] = modes.get(k, np.zeros_like(NOISE))+ap.conj().T@NOISE@aq
    endpoint = sum((a*(1j)**(p % 4) for p, a in aa.items()), np.zeros((4, 4), complex))
    one = -1j*sum((integral_frequency(k)*m for k, m in modes.items()), np.zeros_like(NOISE))/speed
    two = np.zeros((3, 3, 4, 4), complex)
    for k, mk in modes.items():
        for l, ml in modes.items():
            two -= ordered_frequency(k, l)*(mk[:, None]@ml[None, :])/(speed*speed)
    return endpoint, endpoint@one, endpoint@two


def multiply_noise_jets(left, right):
    a, b, c = left
    d, e, f = right
    return a@d, a@e+b@d, a@f+c@d+b[:, None]@e[None, :]


@lru_cache(maxsize=None)
def primitive_noise_jet(n=1):
    """All ordered quadratic coefficients for arbitrary real noise in span(Dd,Dr,Xdr)."""
    result = (I4.copy(), np.zeros_like(NOISE), np.zeros((3, 3, 4, 4), complex))
    for stage in ended_composite(n):
        if stage.arc >= 0:
            u, one, two = arc_noise_jet(stage.arc, n)
            if stage.inverse:
                # Reversal exchanges the later/earlier noise labels as well as
                # taking the matrix adjoint. Scalar contractions hide this swap.
                u, one, two = u.conj().T, -one.conj().transpose(0, 2, 1), two.conj().transpose(1, 0, 3, 2)
        else:
            u = stage.u(stage.duration)
            one = -1j*stage.duration*(u@NOISE)
            two = -stage.duration**2/2*(u@(NOISE[:, None]@NOISE[None, :]))
        result = multiply_noise_jets((u, one, two), result)
    return result


def ordered_primitive_ode(n=1, steps=120):
    """Independent time-ordered coefficients: noise i later than noise j."""
    value = np.zeros((13, 4, 4), complex)
    value[0] = I4
    for stage in ended_composite(n):
        count = max(8, int(np.ceil(stage.duration*steps)))
        dt = stage.duration/count

        def rhs(t, state):
            out = -1j*stage.h(t)@state
            out[1:4] -= 1j*NOISE@state[0]
            out[4:] -= (1j*NOISE[:, None]@state[None, 1:4]).reshape(9, 4, 4)
            return out

        for j in range(count):
            t = j*dt
            k1 = rhs(t, value)
            k2 = rhs(t+dt/2, value+dt*k1/2)
            k3 = rhs(t+dt/2, value+dt*k2/2)
            k4 = rhs(t+dt, value+dt*k3)
            value += dt*(k1+2*k2+2*k3+k4)/6
    return value[4:].reshape(3, 3, 4, 4)


def expand_palindrome(x):
    h = (len(x)-1)//2
    beta = np.r_[x[:h], 0., -x[:h][::-1]]
    stretches = np.r_[x[h:2*h], x[-1], x[h:2*h][::-1]]
    signs = (-1.)**np.arange(2*h+1)
    alpha = 2*np.cumsum(signs*beta)-signs*beta
    return beta, alpha, stretches, signs


def reduced_conditions(x):
    b, a, w, sig = expand_palindrome(x)
    c1, c2 = np.cos(b), np.cos(2*b)
    return np.array([w@c1, w@c2, (w*c1)@np.cos(2*a), (w*c1)@np.sin(2*a),
                     (w*c2)@np.cos(2*a), (w*c2)@np.sin(2*a), w@np.sin(2*a), sig@c2])


def reference_jet(parameter_tuple=None, n=1):
    _, angles, stretches, signs = parameters() if parameter_tuple is None else parameter_tuple
    u, first, second = primitive_noise_jet(n)
    result = [I4.copy(), np.zeros((4, 4), complex), np.zeros((4, 4), complex)]
    for angle, stretch, sign in zip(angles, stretches, signs):
        c, s = np.cos(angle), np.sin(angle)
        coeff = np.array([s*s, c*c, s*c])
        v1 = np.einsum('a,aij->ij', coeff, first)
        v2 = np.einsum('a,b,abij->ij', coeff, coeff, second)
        block = [u, v1, v2]
        if sign < 0:
            block = [z.conj().T for z in block]
        rot = rotation(angle)
        block = [(sign*stretch)**j*(rot@z@rot.T) for j, z in enumerate(block)]
        result = [sum((block[k]@result[j-k] for k in range(j+1)), np.zeros((4, 4), complex))
                  for j in range(3)]
    return result


def logical_second(parameter_tuple=None, n=1):
    jet = reference_jet(parameter_tuple, n)
    m = np.diag([-1., 1.])@jet[2][2:, 2:]
    m -= np.trace(m)*np.eye(2)/2
    return m


def stages_from_parameters(parameter_tuple, n=1):
    _, angles, stretches, signs = parameter_tuple
    base = ended_composite(n)
    return [FramedStage(stage, rotation(a), w, s < 0)
            for a, w, s in zip(angles, stretches, signs)
            for stage in (list(reversed(base)) if s < 0 else base)]


def order_raise(stages, stretch=np.sqrt(2)):
    """Forward, sqrt(2)-stretched physical inverse, forward."""
    middle = [FramedStage(s, I4.real, stretch, True) for s in reversed(stages)]
    return list(stages)+middle+list(stages)


def integrate_exponential_polynomial(terms):
    """Antiderivative from zero, including repeated/zero frequency resonances."""
    out = {}
    zero = np.zeros_like(next(iter(terms.values())))

    def add(key, value):
        out[key] = out.get(key, zero)+value

    for (k, degree), value in terms.items():
        if k == 0:
            add((0, degree+1), value/(degree+1))
        else:
            for j in range(degree+1):
                c = (-1)**j*factorial(degree)/factorial(degree-j)/(1j*k)**(j+1)
                add((k, degree-j), c*value)
            add((0, 0), -(-1)**degree*factorial(degree)/(1j*k)**(degree+1)*value)
    return out


def fourier_arc_jet(arc, noises, n=1, order=3):
    """Exact finite sums for arbitrary-order Dyson coefficients; no ODE or fit."""
    aa, speed = arc_fourier(arc, n)
    modes = {}
    for p, ap in aa.items():
        for q, aq in aa.items():
            k = q-p
            modes[k] = modes.get(k, np.zeros_like(noises))+ap.conj().T@noises@aq
    endpoint = sum((a*(1j)**(p % 4) for p, a in aa.items()), np.zeros((4, 4), complex))
    initial = np.broadcast_to(I4, noises.shape).copy()
    terms = {(0, 0): initial}
    result = [initial@ endpoint]
    for _ in range(order):
        product = {}
        for k, mk in modes.items():
            for (p, degree), value in terms.items():
                key = (k+p, degree)
                product[key] = product.get(key, np.zeros_like(noises))-1j*mk@value/speed
        terms = integrate_exponential_polynomial(product)
        value = sum((c*(np.pi/2)**degree*(1j)**(k % 4)
                     for (k, degree), c in terms.items()), np.zeros_like(noises))
        result.append(endpoint@value)
    return np.stack(result, axis=1)


def multiply_scalar_jets(left, right):
    return np.stack([sum((left[..., k, :, :]@right[..., j-k, :, :]
                          for k in range(j+1)), np.zeros_like(left[..., 0, :, :]))
                     for j in range(left.shape[-3])], axis=-3)


def compose_primitive_jets(jets, parameter_tuple):
    _, angles, stretches, signs = parameter_tuple
    result = np.zeros_like(jets[0])
    result[0] = I4
    for j, (angle, stretch, sign) in enumerate(zip(angles, stretches, signs)):
        rot = rotation(angle)
        block = jets[j]
        if sign < 0:
            block = block.conj().transpose(0, 2, 1)
        block = np.array([(sign*stretch)**q*(rot@value@rot.T) for q, value in enumerate(block)])
        result = multiply_scalar_jets(block, result)
    return result


def higher_reference_jet(parameter_tuple=None, n=1, order=3, steps=None):
    """Fourier recurrence, or an independent RK4 derivative integration."""
    par = parameters() if parameter_tuple is None else parameter_tuple
    angles = par[1]
    rotations = np.array([rotation(a) for a in angles])
    noises = rotations.transpose(0, 2, 1)@NOISE[1]@rotations
    result = np.zeros((len(angles), order+1, 4, 4), complex)
    result[:, 0] = I4
    arc_cache = {}
    for stage in ended_composite(n):
        if steps is not None:
            count = max(8, int(np.ceil(stage.duration*steps)))
            dt = stage.duration/count

            def rhs(t, value):
                out = -1j*stage.h(t)@value
                out[:, 1:] -= 1j*noises[:, None]@value[:, :-1]
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
            if stage.arc not in arc_cache:
                arc_cache[stage.arc] = fourier_arc_jet(stage.arc, noises, n, order)
            block = arc_cache[stage.arc]
            if stage.inverse:
                block = block.conj().transpose(0, 1, 3, 2)*((-1.)**np.arange(order+1))[None, :, None, None]
        else:
            u = stage.u(stage.duration)
            block = np.empty_like(result)
            block[:, 0] = u
            for q in range(1, order+1):
                block[:, q] = (-1j*stage.duration)**q/factorial(q)*(u@noises)
        result = multiply_scalar_jets(block, result)
    return compose_primitive_jets(result, par)


def raised_jet(jet, stretch=np.sqrt(2)):
    inv = jet.conj().transpose(0, 2, 1)*((-stretch)**np.arange(len(jet)))[:, None, None]
    return multiply_scalar_jets(jet, multiply_scalar_jets(inv, jet))


def finite_words(deltas, error=0., loss=0., parameter_tuple=None, n=1, steps=120):
    """Batched independent laboratory ODEs; one native C waveform per case/angle."""
    par = parameters() if parameter_tuple is None else parameter_tuple
    _, angles, stretches, signs = par
    rotations = np.array([rotation(a) for a in angles])
    noises = rotations.transpose(0, 2, 1)@NOISE[1]@rotations
    offsets = np.array(deltas)[:, None]*stretches[None, :]*signs[None, :]
    rates = np.broadcast_to(np.asarray(loss), (len(deltas),))[:, None]*stretches[None, :]
    value = np.broadcast_to(I4, (len(deltas), len(angles), 4, 4)).copy()
    for stage in ended_composite(n):
        count = max(8, int(np.ceil(stage.duration*steps)))
        dt = stage.duration/count

        def rhs(t, v):
            h = (1+error)*stage.h(t)+offsets[..., None, None]*noises[None, :]
            generator = -1j*h
            if np.any(rates):
                generator -= rates[..., None, None]*stage.loss_projector(t)/2
            return generator@v

        for j in range(count):
            t = j*dt
            k1 = rhs(t, value)
            k2 = rhs(t+dt/2, value+dt*k1/2)
            k3 = rhs(t+dt/2, value+dt*k2/2)
            k4 = rhs(t+dt, value+dt*k3)
            value += dt*(k1+2*k2+2*k3+k4)/6
    out = np.broadcast_to(I4, (len(deltas), 4, 4)).copy()
    for j, (rot, sign) in enumerate(zip(rotations, signs)):
        block = value[:, j]
        if sign < 0:
            block = block.conj().transpose(0, 2, 1)
        out = rot@block@rot.T@out
    return out


def control_parameters(n=1):
    record = json.loads(CONTROL_PATH.read_text())[str(n)]
    return expand_palindrome(np.array(record['center'], float))


def gain_word(error=0., n=1, parameter_tuple=None):
    par = parameters() if parameter_tuple is None else parameter_tuple
    u, out = composite(error, n), I4.copy()
    for a, sign in zip(par[1], par[3]):
        rot = rotation(a)
        out = rot@(u if sign > 0 else u.conj().T)@rot.T@out
    return out


def time_moments(stages, degree=3, nodes=64):
    points, weights = np.polynomial.legendre.leggauss(nodes)
    result = np.zeros((degree+1, 2, 8), complex)
    prefix, offset = I4.copy(), 0.
    duration = sum(s.duration for s in stages)
    for stage in stages:
        for point, wt in zip(points, weights):
            t = stage.duration*(point+1)/2
            u = stage.u(t)@prefix
            k = bad_response(u.conj().T@np.array([stage.h(t), NOISE[1]])@u)
            moment_time = offset+t-duration/2
            result += (moment_time**np.arange(degree+1))[:, None, None]*k*wt*stage.duration/2
        prefix = stage.u(stage.duration)@prefix
        offset += stage.duration
    return result


def symbolic_order_map():
    import sympy as sp
    a, b, c = sp.symbols('A B C', commutative=False)
    s = sp.symbols('s', real=True)
    forward = [1, a, b+a*a/2, c+(a*b+b*a)/2+a**3/6]
    inverse = [1, s*a, s*s*(-b+a*a/2), s**3*(c-(a*b+b*a)/2+a**3/6)]

    def product(x, y):
        return [sp.expand(sum(x[k]*y[j-k] for k in range(j+1))) for j in range(4)]

    word = product(forward, product(inverse, forward))
    logs = [word[1], word[2]-word[1]**2/2,
            word[3]-(word[1]*word[2]+word[2]*word[1])/2+word[1]**3/3]
    expected = [(2+s)*a, (2-s*s)*b, (2+s**3)*c]
    for actual, target in zip(logs, expected):
        assert sp.expand(actual-target) == 0
    assert sp.simplify((2-s*s).subs(s, sp.sqrt(2))) == 0
    return dict(first_log='(2+s) A', second_log='(2-s^2) B', third_log='(2+s^3) C',
                log_second_cancelled=True, gain_at_zero_detuning_preserved_exactly=True)


def run_checks():
    from reflection_loop_detuning_certificate import certify, primitive as precise_primitive
    import mpmath as mp
    report = dict(order_map=symbolic_order_map(), checks={}, cases={}, certificates={})
    controls = json.loads(CONTROL_PATH.read_text())
    checks = report['checks']

    def close(name, actual, expected=0., tolerance=1e-9):
        error = float(np.max(np.abs(np.asarray(actual)-expected)))
        assert error <= tolerance, (name, error, tolerance)
        checks[name] = dict(error=error, tolerance=tolerance)

    def bound(name, value, maximum):
        assert value <= maximum, (name, value, maximum)
        checks[name] = dict(value=float(value), maximum=float(maximum))

    frequencies = np.r_[0., np.geomspace(1e-6, .002, 100)]
    jq = I4.real.copy()
    jq[:2, :2] = [[0, 1], [1, 0]]
    for n in (1, 2, 3):
        report['certificates'][str(n)] = certify(controls[str(n)])
        par = control_parameters(n)
        close(f'n{n}_circular_equations', reduced_conditions(np.array(controls[str(n)]['center'], float)), tolerance=2e-12)
        st = stages_from_parameters(par, n)
        parent = compressed_stages(n)
        raised = order_raise(parent)
        moments = response_integrals(st, noise_generators())
        close(f'n{n}_full_nominal_endpoint', moments['endpoint'], composite(), 5e-12)
        for k, m in enumerate(moments['moments']):
            close(f'n{n}_six_noise_generator_{k}', bad_response(m), tolerance=3e-9)
        xc = 3*coefficients(n)['exposure']
        close(f'n{n}_isotropic_erasure', moments['loss'][:, 2:], par[2].sum()*xc/2*I4[:, 2:], 1e-9)
        for kind, stages in (('thirteen', st), ('raised', raised)):
            close(f'{kind}_n{n}_zero_endpoints',
                  max(np.linalg.norm(stages[0].h(0)), np.linalg.norm(stages[-1].h(stages[-1].duration))), tolerance=1e-12)
            close(f'{kind}_n{n}_continuous_joins',
                  max(np.linalg.norm(l.h(l.duration)-r.h(0)) for l, r in zip(stages, stages[1:])), tolerance=2e-12)
        # Exact native time-reflection symmetry supplies the missing logical
        # component; the certificate solves the other two real components.
        symmetry_error = 0.
        for left, right in zip(st, reversed(st)):
            for f in (.2, .7):
                symmetry_error = max(symmetry_error, np.linalg.norm(
                    left.h(f*left.duration)-jq@right.h((1-f)*right.duration)@jq))
        close(f'n{n}_complete_waveform_reciprocity', symmetry_error, tolerance=3e-12)
        analytic = higher_reference_jet(par, n, 3)
        fast = np.array(reference_jet(par, n))
        close(f'n{n}_independent_Fourier_recurrences', analytic[:3]/(1+np.abs(fast)), fast/(1+np.abs(fast)), 2e-10)
        mp.mp.dps = 70
        precise, _ = precise_primitive(mp.mp, n)
        pu, p1, p2 = primitive_noise_jet(n)
        close(f'n{n}_independent_precision_first_jet',
              np.array([a.tolist() for a in precise[1]], complex), p1, 2e-10)
        close(f'n{n}_independent_precision_second_jet',
              np.array([[a.tolist() for a in row] for row in precise[2]], complex), p2, 2e-8)
        relative = composite().conj().T@analytic
        pp2 = relative[2][2:, 2:]
        close(f'n{n}_complete_logical_detuning_second_order',
              pp2-np.trace(pp2)*np.eye(2)/2, tolerance=2e-8)
        numerical = compressed_jets(120, n, par)
        response = quartic_response(numerical, target=composite())
        close(f'n{n}_mixed_gain_detuning_leakage', response['leakage'][1, 1], tolerance=5e-6)
        close(f'n{n}_mixed_gain_detuning_logical', response['logical'][1, 1], tolerance=2e-6)
        close(f'n{n}_ODE_logical_detuning_second_order', response['logical'][0, 2], tolerance=2e-6)
        # Use the exact finite sums for the reported pure-reference coefficient.
        delta4 = float(np.linalg.norm(relative[2][:2, 2:])**2/2)
        b = par[0]
        v = np.array([par[3]@np.cos(b), par[3]@np.sin(b)])
        rn = 1-1/(16*n*n)
        lc = np.pi**2*rn/12*np.array([-2*np.sqrt(3)+3j, np.sqrt(3)-6j])
        lg = np.outer(lc, v)
        gain4 = float(np.linalg.norm(lg)**2/2)
        joint = dict(gain4=gain4, gain2_delta2=float(np.vdot(lg, relative[2][:2, 2:]).real), delta4=delta4)
        close(f'n{n}_analytic_gain_column', response['leakage'][2, 0], lg, 2e-6)
        originals = higher_reference_jet(n=n, order=3)
        rjet = raised_jet(originals)
        rr = composite().conj().T@rjet
        close(f'n{n}_raised_first_order_code', bad_response(rr[1]), tolerance=2e-9)
        close(f'n{n}_raised_complete_second_order_code', bad_response(rr[2]), tolerance=2e-8)
        delta6 = float(np.linalg.norm(bad_response(rr[3]))**2)
        tm = time_moments(raised)
        close(f'n{n}_raised_static_and_linear_reference', tm[:2, 1], tolerance=1e-7)
        close(f'n{n}_raised_constant_linear_quadratic_gain', tm[:3, 0], tolerance=1e-6)
        response_g6 = float(np.linalg.norm(tm[3, 0])**2/36)
        response_d4 = float(np.linalg.norm(tm[2, 1])**2/4)
        tf = temporal_response(raised, frequencies)
        close(f'n{n}_gain_spectral_sixth_order',
              (tf['gram'][1, 0, 0].real+tf['negative_gram'][1, 0, 0].real)/(2*frequencies[1]**6*response_g6), 1., .025)
        close(f'n{n}_reference_spectral_fourth_order',
              (tf['gram'][1, 1, 1].real+tf['negative_gram'][1, 1, 1].real)/(2*frequencies[1]**4*response_d4), 1., .001)
        common = dict(n=n, amplitude_cap=1., slew_cap=10.)
        rows = {
            'thirteen': dict(**common, duration=sum(s.duration for s in st), mean_exposure=par[2].sum()*xc/2,
                              leading_joint=joint, reference_order=4, reference_coefficient=delta4,
                              reference_jets=[packed_matrix(a) for a in analytic],
                              numerical_joint=packed_response(response)),
            'raised': dict(**common, duration=sum(s.duration for s in raised),
                           mean_exposure=(2+np.sqrt(2))*parameters()[2].sum()*xc/2,
                           reference_order=6, reference_coefficient=delta6,
                           reference_jets=[packed_matrix(a) for a in rjet],
                           gain_frequency_order=6, gain_frequency_coefficient=response_g6,
                           reference_frequency_order=4, reference_frequency_coefficient=response_d4,
                           temporal_frequency=frequencies.tolist(),
                           temporal_reference=((tf['gram'][:, 1, 1].real+tf['negative_gram'][:, 1, 1].real)/2).tolist()),
        }
        for delta, label in ((.0001, 'original'), (.001, 'detuning_dominated')):
            direct = finite_words([delta], error=.001, parameter_tuple=par, n=n)[0]
            uncorrected, reverse = finite_words([delta, -np.sqrt(2)*delta], error=.001, n=n)
            complete = uncorrected@reverse.conj().T@uncorrected
            rows['thirteen'][label] = metrics(direct)
            rows['raised'][label] = metrics(complete)
            bound(f'n{n}_{label}_finite_unitarity', max(np.linalg.norm(direct.conj().T@direct-I4),
                  np.linalg.norm(complete.conj().T@complete-I4)), 2e-8)
        pure13 = finite_words([1e-4, 5e-5], parameter_tuple=par, n=n)
        bound(f'n{n}_logical_detuning_starts_at_cubic_order',
              metrics(pure13[1])['logical_traceless_norm']/metrics(pure13[0])['logical_traceless_norm'], .14)
        for family, row in rows.items():
            report['cases'][f'{family}_n{n}'] = row
        print(f'n={n}: exact root, finite response, order raising and temporal cancellation checked', flush=True)
        if n == 1:
            close('ordered_multichannel_coefficients_against_ODE',
                  ordered_primitive_ode(n, 120), p2, 2e-7)
            coarser = compressed_jets(60, n, par)
            fine_error = np.linalg.norm(response['logical'][0, 2])
            coarse_error = np.linalg.norm(quartic_response(coarser, target=composite())['logical'][0, 2])
            bound('logical_cancellation_ODE_convergence', fine_error/coarse_error, .09)
            rk3 = higher_reference_jet(par, n, 3, steps=120)
            close('Fourier_third_order_against_independent_ODE',
                  rk3[3]/(1+np.abs(analytic[3])), analytic[3]/(1+np.abs(analytic[3])), 2e-7)
            x1, x2 = finite_words([1e-4, -np.sqrt(2)*1e-4], n=n, steps=240)
            pure = x1@x2.conj().T@x1
            close('sextic_reference_finite_witness', metrics(pure)['average_infidelity']/(delta6*1e-24), 1., .025)
            raw = ode(raised, error=.001, detuning=1e-4, steps_per_time=120)
            batch = finite_words([1e-4, -np.sqrt(2)*1e-4], error=.001, steps=240)
            close('whole_physical_waveform_against_batched_construction', raw,
                  batch[0]@batch[1].conj().T@batch[0], 2e-8)
            rate = 1e-8
            clean = finite_words([0., 0.])
            damped = finite_words([0., 0.], loss=[rate, np.sqrt(2)*rate])
            uk = damped[0]@damped[1].conj().T@damped[0]
            u0 = clean[0]@clean[1].conj().T@clean[0]
            close('raised_state_independent_loss_slope',
                  (metrics(uk, True)['average_infidelity']-metrics(u0, True)['average_infidelity'])/rate,
                  rows['raised']['mean_exposure'], .02)
    baseline = json.loads((CONTROL_PATH.parent/'reflection-loop-compression-checks.json').read_text())
    for name, old in baseline['cases'].items():
        n = old['n']
        u = (finite_gate(.001, .0001, n=n) if name.startswith('echo')
             else finite_words([.0001], error=.001, n=n)[0])
        close(name+'_baseline_replay', metrics(u)['average_infidelity']/old['joint_metrics']['average_infidelity'], 1., .0005)
        freq = old['temporal']['frequency']
        pg = np.array(old['temporal']['gram']['real'])[:, 1, 1]
        ng = np.array(old['temporal']['negative_gram']['real'])[:, 1, 1]
        report['cases'][name] = dict(n=n, duration=old['duration'], mean_exposure=old['mean_exposure'],
            original=metrics(u), reference_order=4,
            reference_coefficient=old['response']['terms']['0,4'],
            temporal_frequency=freq, temporal_reference=((pg+ng)/2).tolist())
    for n in (1, 2, 3):
        row = report['cases'][f'thirteen_n{n}']
        tr = temporal_response(stages_from_parameters(control_parameters(n), n), frequencies)
        row['temporal_frequency'] = frequencies.tolist()
        row['temporal_reference'] = ((tr['gram'][:, 1, 1].real+tr['negative_gram'][:, 1, 1].real)/2).tolist()
    report['budget_examples'] = []
    for rate in (1e-13, 1e-10, 1e-9):
        candidates = [(r['original']['average_infidelity']+rate*r['mean_exposure'], name)
                      for name, r in report['cases'].items() if r['duration'] <= 3000]
        old = [c for c in candidates if c[1].startswith(('echo_', 'compressed_'))]
        best, previous = min(candidates), min(old)
        report['budget_examples'].append(dict(duration_limit=3000, gain=.001, delta=.0001,
            loss_rate=rate, selected=best[1], predicted_infidelity=best[0],
            previous_best=previous[1], previous_infidelity=previous[0], improvement=previous[0]/best[0]))
    report['check_count'] = len(checks)
    report['scope'] = ('The thirteen-block roots cancel logical reference-detuning order two at n=1,2,3; '
                       'leakage remains at that order. The raised word cancels the complete quadratic '
                       'static-noise logarithm for arbitrary static noise vectors, and retains exactly '
                       'the parent gain-only gate. Resource comparisons use specified noise and duration.')
    return report


def plot(path, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.5), layout='constrained')
    colors = dict(compressed='#858c96', echo='#216d8a', thirteen='#b95e2a', raised='#65549c')
    labels = dict(compressed='11 blocks', echo='Echo', thirteen='13 blocks', raised='Full second order')
    markers = dict(compressed='x', echo='o', thirteen='s', raised='D')
    for family in colors:
        color = colors[family]
        for n in (1, 2, 3):
            row = report['cases'][f'{family}_n{n}']
            axes[0, 0].scatter(row['duration'], row['mean_exposure'], c=color, marker=markers[family],
                                s=45, label=labels[family] if n == 1 else None)
            axes[0, 0].annotate(str(n), (row['duration'], row['mean_exposure']),
                                xytext=(5, 5), textcoords='offset points', fontsize=9)
        row = report['cases'][family+'_n1']
        delta = np.geomspace(1e-6, 3e-4, 100)
        axes[0, 1].loglog(delta, row['reference_coefficient']*delta**row['reference_order'],
                          color=color, label=labels[family])
        freq = np.array(row['temporal_frequency'])[1:]
        axes[1, 1].loglog(freq, row['temporal_reference'][1:], color=color, label=labels[family])
    rates = np.geomspace(1e-14, 1e-7, 180)
    for family, color in colors.items():
        rows = [r for name, r in report['cases'].items()
                if name.startswith(family+'_') and r['duration'] <= 3000]
        values = np.min([r['original']['average_infidelity']+rates*r['mean_exposure'] for r in rows], axis=0)
        axes[1, 0].loglog(rates, values, color=color, label=labels[family])
    axes[0, 0].set(xlabel=r'Duration $aT$', ylabel=r'Mean bright exposure $a\bar X$',
                   title='Resource cost; labels give return index n', xlim=(0, 8000), ylim=(0, 150))
    axes[0, 0].legend(frameon=False, fontsize=8, loc='upper right')
    axes[0, 1].set(xlabel=r'Static reference offset $\Delta/a$', ylabel='Leading average infidelity',
                   title='Pure detuning, n=1: quartic versus sextic')
    axes[0, 1].legend(frameon=False, fontsize=8)
    axes[1, 0].set(xlabel=r'Absorbing loss rate $\kappa/a$', ylabel='Finite coherent error + leading loss',
                   title=r'Choose among n = 1–3 with $aT\leq3000$')
    axes[1, 0].text(.04, .91, r'$\epsilon=10^{-3},\ \Delta/a=10^{-4}$',
                    transform=axes[1, 0].transAxes, fontsize=9)
    axes[1, 0].legend(frameon=False, fontsize=8, loc='lower right')
    axes[1, 1].set(xlabel=r'Noise angular frequency $\omega/a$',
                   ylabel=r'Reference response $a^2 F_\Delta(\omega)$',
                   title='Reference drift: second versus fourth frequency order')
    for ax in axes.flat:
        ax.grid(alpha=.2, which='both')
    fig.suptitle('Second-order detuning cancellation with explicit duration, loss and temporal response')
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
        Path(args.json).write_text(output+'\n')
    else:
        print(output)
    if args.plot:
        plot(args.plot, report)

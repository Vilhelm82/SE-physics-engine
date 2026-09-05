#!/usr/bin/env python3
"""Real reflection-loop composite: exact response, finite ramps, independent ODE.

Uses the native frames from reflection_loop_dynamics.py. No external pulse
sequence or numerical code is imported. Hamiltonians have frequency units.
"""
import argparse
from dataclasses import dataclass
import json

import numpy as np
import sympy as sp

from reflection_loop_dynamics import (
    ANGLE, J, P, Q, R, embed, exact_loop, frame, hamiltonian, segment,
)

THETA = np.pi / 3
J4 = np.zeros((4, 4), complex)
J4[:3, :3] = J
REFERENCE = np.diag([0, 0, 0, 1])


def return_time(n=1, gap=1.0):
    if n < 1 or int(n) != n or gap <= 0:
        raise ValueError('n must be a positive integer and gap must be positive')
    return ANGLE * np.sqrt(16*n*n - 1) / gap


def hold(area):
    return embed(np.eye(3) - 1j*np.sin(area)*J
                 - 2*np.sin(area/2)**2*(J @ J))


def primitive(error=0.0, n=1, gap=1.0):
    tau = return_time(n, gap)
    return exact_loop((tau,)*3, (gap*(1+error),)*3)


def composite(error=0.0, n=1, gap=1.0, loop_errors=None, hold_errors=None):
    """Chronological F, S, F inverse, S, F; hold gains are imperfect too."""
    es = (error,)*3 if loop_errors is None else loop_errors
    hs = (error,)*2 if hold_errors is None else hold_errors
    f1, f2, f3 = (primitive(e, n, gap) for e in es)
    s1, s2 = (hold(THETA*(1+h)) for h in hs)
    return f3 @ s2 @ f2.conj().T @ s1 @ f1


def unitary_metrics(u):
    """Average logical infidelity includes leakage; reference is uncoupled."""
    z, reference = -u[2, 2], u[3, 3]
    leakage = float(np.linalg.norm(u[:2, 2])**2)
    return dict(leakage=leakage, relative_phase=float(np.angle(z*reference.conjugate())),
                average_infidelity=float(leakage/2 + abs(z-reference)**2/6))


def average_infidelity(k):
    """Also valid for no-jump K when all jumps exit to orthogonal sinks."""
    a = np.diag([-1, 1]) @ k[2:, 2:]
    return float(1 - (np.trace(a.conj().T @ a).real + abs(np.trace(a))**2)/6)


@dataclass(frozen=True)
class Stage:
    duration: float
    gap: float = 1.0
    arc: int = -1
    inverse: bool = False
    q_start: float = 0.0
    q_end: float = 0.0

    def h(self, t):
        if self.arc < 0:
            return (self.q_start + (self.q_end-self.q_start)*t/self.duration)*J4
        local = self.duration-t if self.inverse else t
        sign = -1 if self.inverse else 1
        return sign*hamiltonian(self.arc, local, self.duration, self.gap)

    def loss_projector(self, t):
        # The physical susceptibility is continued through the zero-gap hold.
        if self.arc < 0:
            return Q
        local = self.duration-t if self.inverse else t
        dark = frame(self.arc, ANGLE*local/self.duration)[:, 2]
        out = np.zeros((4, 4), complex)
        out[:3, :3] = np.eye(3)-np.outer(dark, dark)
        return out

    def u(self, t):
        """Nominal within-stage solution, independent of the composite product."""
        if self.arc < 0:
            area = self.q_start*t + (self.q_end-self.q_start)*t*t/(2*self.duration)
            return hold(area)

        def forward(local):
            f = frame(self.arc, ANGLE*local/self.duration)
            f0 = frame(self.arc, 0)
            e = segment(self.arc, self.duration, self.gap, local)
            return embed(f @ e @ f0.T)

        if self.inverse:
            return forward(self.duration-t) @ forward(self.duration).conj().T
        return forward(t)


def schedule(n=1, gap=1.0, ramp=None):
    """Continuous real waveform; two sign changes pass through H=0."""
    tau = return_time(n, gap)
    ramp = 0.2/gap if ramp is None else ramp
    if ramp <= 0:
        raise ValueError('a continuous sign switch requires positive ramp duration')
    f = [Stage(tau, gap, arc) for arc in range(3)]
    inv = [Stage(tau, gap, arc, True) for arc in (2, 1, 0)]
    plateau = Stage(THETA/gap, gap, q_start=gap, q_end=gap)
    down = Stage(ramp, gap, q_start=gap, q_end=-gap)
    up = Stage(ramp, gap, q_start=-gap, q_end=gap)
    return f + [plateau, down] + inv + [up, plateau] + f


def ode(stages, error=0.0, detuning=0.0, loss=0.0, gain_shape=None,
        steps_per_time=100):
    """Laboratory RK4, including the signed ramps and optional absorbing loss.

    gain_shape(global_time) specifies an additional dimensionless gain error.
    No moving-frame or matrix-product solution is used in the integration.
    """
    k = np.eye(4, dtype=complex)
    offset = 0.0
    for stage in stages:
        steps = max(8, int(np.ceil(stage.duration*steps_per_time)))
        dt = stage.duration/steps

        def rhs(t, state):
            gain = error + (0 if gain_shape is None else gain_shape(offset+t))
            h = (1+gain)*stage.h(t) + detuning*REFERENCE
            generator = -1j*h
            if loss:
                generator = generator - loss*stage.loss_projector(t)/2
            return generator @ state

        for j in range(steps):
            t = j*dt
            k1 = rhs(t, k)
            k2 = rhs(t+dt/2, k+dt*k1/2)
            k3 = rhs(t+dt/2, k+dt*k2/2)
            k4 = rhs(t+dt, k+dt*k3)
            k += dt*(k1+2*k2+2*k3+k4)/6
        offset += stage.duration
    return k


def nominal_integrals(stages, nodes=96):
    """DC/linear gain response, exposure, and a finite-variation error bound."""
    points, weights = np.polynomial.legendre.leggauss(nodes)
    total_time = sum(stage.duration for stage in stages)
    prefix = np.eye(4, dtype=complex)
    m0, m1 = np.zeros((4, 4), complex), np.zeros((4, 4), complex)
    exposure, variation_bound, parallel_error = 0.0, 0.0, 0.0
    offset = 0.0
    for stage in stages:
        for point, weight in zip(points, weights):
            t = stage.duration*(point+1)/2
            wt = weight*stage.duration/2
            u = stage.u(t) @ prefix
            b = u.conj().T @ stage.h(t) @ u
            normalized_time = (offset+t-total_time/2)/total_time
            m0 += wt*b
            m1 += wt*normalized_time*b
            psi = u[:, 2]
            exposure += wt*np.vdot(psi, stage.loss_projector(t) @ psi).real
            variation_bound += wt*abs(offset+t-total_time/2)*np.linalg.norm(b[:2, 2])
            parallel_error = max(parallel_error, float(np.linalg.norm(P @ b @ P)))
        prefix = stage.u(stage.duration) @ prefix
        offset += stage.duration
    return dict(endpoint=prefix, dc=m0, linear=m1, exposure=float(exposure),
                variation_bound=float(variation_bound), parallel_error=parallel_error)


def coefficients(n=1):
    r = 1-1/(16*n*n)
    b = ANGLE*r
    tau = return_time(n)
    v2 = 1/(16*n*n)
    return dict(b=b, primitive_leakage=5*b*b, composite_leakage=15*b*b*THETA**2,
                composite_phase=2*np.sqrt(3)*b*b,
                composite_infidelity=7.5*b*b*THETA**2+2*b**4,
                exposure=3*tau*(2*v2-1.5*v2*v2),
                peak_bright=4*v2*(1-v2), primitive_time=3*tau)


def run_checks():
    symbolic, numerical = {}, {}

    def exact(name, actual, expected=0):
        if isinstance(actual, sp.MatrixBase):
            expected = sp.zeros(*actual.shape) if expected == 0 else expected
            entries = list(actual-expected)
        else:
            entries = [actual-expected]
        assert all(sp.simplify(sp.expand(x)) == 0 for x in entries), name
        symbolic[name] = True

    def close(name, actual, expected, tolerance=1e-10):
        err = float(np.max(np.abs(np.asarray(actual)-np.asarray(expected))))
        assert err <= tolerance, (name, err, tolerance)
        numerical[name] = dict(error=err, tolerance=tolerance)

    def bounded(name, actual, bound):
        assert actual <= bound, (name, actual, bound)
        numerical[name] = dict(value=float(actual), upper_bound=float(bound))

    sj = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    sr = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, -1]])
    g1 = sp.Matrix([[0, 0, sp.I], [0, 0, 0], [-sp.I, 0, 0]])
    g2 = sp.Matrix([[0, 0, 0], [0, 0, -sp.I], [0, sp.I, 0]])
    r, lam = sp.symbols('r lambda', real=True)
    theta, q = sp.pi/3, sp.pi/2
    s = sp.eye(3)-sp.I*sp.sqrt(3)*sj/2-sj*sj/2
    bright = sp.diag(1, 1, 0)
    exact('hold_polynomial_on_full_bright_plane', (s*s-s+sp.eye(3))*bright)
    exact('hold_unitary', s.H*s, sp.eye(3))
    exact('hold_commutes_with_reflection', s*sr-sr*s)

    def mul(a, b):
        return [sum((a[k]*b[j-k] for k in range(j+1)), sp.zeros(3))
                .applyfunc(sp.simplify) for j in range(3)]

    v = [sp.eye(3), sp.zeros(3), sp.zeros(3)]
    for g in (g1, g2, g1):
        ell = lam*sj+q*g
        e = [sp.eye(3), -sp.I*r*ell,
             -sp.I*(r*(1-3*r)*ell/2+r*lam*sj)-r*r*ell*ell/2]
        v = mul(e, v)
    f = [sr*x for x in v]
    h = [s, -sp.I*theta*s*sj, -theta*theta*s*sj*sj/2]
    c = mul(mul(mul(mul(f, h), [x.H for x in f]), h), f)
    w1, w2 = c[0].H*c[1], c[0].H*c[2]
    expected = sp.Matrix([sp.pi**2*r*(-2*sp.sqrt(3)+3*sp.I)/12,
                          sp.pi**2*r*(sp.sqrt(3)-6*sp.I)/12,
                          sp.sqrt(3)*sp.I*sp.pi**2*r*r/2])
    exact('composite_nominal_full_operator', c[0], sr*s*s)
    exact('common_gain_first_order_on_code', w1[:, 2])
    exact('full_second_order_code_column_all_return_indices', w2[:, 2], expected)
    leak4 = sum(sp.conjugate(x)*x for x in w2[:2, 2])
    exact('quartic_leakage_coefficient', leak4, 5*sp.pi**4*r*r/12)
    exact('quadratic_relative_phase', sp.im(w2[2, 2]), sp.sqrt(3)*sp.pi**2*r*r/2)
    exact('primitive_has_no_quadratic_phase', sp.im(v[2][2, 2]))
    exact('n1_leakage_coefficient', leak4.subs(r, sp.Rational(15, 16)),
          375*sp.pi**4/1024)
    h1, h2 = sp.symbols('h_1 h_2', real=True)
    def hold_jet(gain):
        return [s, -sp.I*theta*gain*s*sj, -theta*theta*gain*gain*s*sj*sj/2]
    ch = mul(mul(mul(mul(f, hold_jet(h2)), [x.H for x in f]), hold_jet(h1)), f)
    wh2 = ch[0].H*ch[2]
    vec = sp.Matrix([2, -1, 0])
    joint_bright = -sp.I*q*r*theta*sj*(h1*sp.eye(3)+h2*s)*vec
    exact('separate_hold_gain_second_order_bright_response',
          bright*wh2[:, 2], joint_bright)
    exact('separate_hold_gain_phase_response', sp.im(wh2[2, 2]),
          sp.sqrt(3)*sp.pi**2*r*r/2)
    exact('separate_hold_gain_leakage_form', (joint_bright.H*joint_bright)[0],
          5*q*q*r*r*theta*theta*(h1*h1+h2*h2+h1*h2))
    e1, e2, e3 = sp.symbols('epsilon_1 epsilon_2 epsilon_3', real=True)
    mismatch = (e1*sp.eye(3)-e2*s.H+e3*s.H*s.H)*vec
    exact('independent_loop_gain_quadratic_form', (mismatch.H*mismatch)[0],
          sp.Rational(5, 2)*((e1-e2)**2+(e2-e3)**2+(e3-e1)**2))
    a, nu, time = sp.symbols('a nu time', positive=True)
    rho = 2*a/nu
    exact('linear_sign_ramp_zero_area', sp.integrate(a-2*a*time/rho, (time, 0, rho)))
    exact('sign_ramp_area_deficit', sp.integrate(2*a*time/rho, (time, 0, rho)),
          2*a*a/nu)

    for n in (1, 2, 4):
        coeff = coefficients(n)
        close(f'nominal_code_n{n}', composite(n=n) @ P, embed(R) @ P)
        eps = 2e-4/n
        samples = [unitary_metrics(composite(sign*eps, n)) for sign in (-1, 1)]
        leak = sum(x['leakage'] for x in samples)/(2*eps**4)
        phase = sum(x['relative_phase'] for x in samples)/(2*eps**2)
        close(f'quartic_leakage_n{n}', leak/coeff['composite_leakage'], 1, 4e-5)
        close(f'quadratic_phase_n{n}', phase/coeff['composite_phase'], 1, 4e-5)
        infid = sum(x['average_infidelity'] for x in samples)/(2*eps**4)
        close(f'quartic_average_infidelity_n{n}', infid/coeff['composite_infidelity'],
              1, 4e-5)

    eps = 1e-5
    coeff = coefficients()
    for index, direction in enumerate(((1, 0, 0), (0, 1, 0), (1, -2, 3)), 1):
        target = coeff['primitive_leakage']/2*sum(
            (direction[j]-direction[(j+1) % 3])**2 for j in range(3))
        observed = sum(unitary_metrics(composite(
            loop_errors=[sign*eps*x for x in direction], hold_errors=(0, 0)))['leakage']
            for sign in (-1, 1))/(2*eps*eps)
        close(f'between_loop_drift_{index}', observed/target, 1, 3e-6)
    close('hold_errors_alone_leave_code_exact',
          composite(loop_errors=(0, 0, 0), hold_errors=(0.21, -0.17)) @ P,
          embed(R) @ P)
    eps = 1e-4
    for index, gains in enumerate(((0, 0), (0.3, -0.7)), 1):
        h1n, h2n = gains
        expected_bright = (-1j*coeff['b']*THETA*J
                           @ (h1n*np.eye(3)+h2n*hold(THETA)[:3, :3])
                           @ np.array([2, -1, 0]))
        columns = [composite().conj().T @ composite(sign*eps,
                   hold_errors=tuple(sign*eps*x for x in gains))[:, 2]
                   for sign in (-1, 1)]
        observed = (columns[0]+columns[1])/(2*eps*eps)
        close(f'separate_hold_gain_bright_response_{index}',
              observed[:2], expected_bright[:2], 3e-5)

    stages = schedule()
    total_time = sum(stage.duration for stage in stages)
    close('duration_includes_both_finite_ramps', total_time,
          3*coeff['primitive_time']+2*THETA+0.4)
    continuity = max(np.linalg.norm(left.h(left.duration)-right.h(0))
                     for left, right in zip(stages, stages[1:]))
    bounded('physical_hamiltonian_continuous_at_all_joins', continuity, 2e-14)
    for index in (4, 8):
        close(f'sign_ramp_{index}_passes_through_zero',
              stages[index].h(stages[index].duration/2), np.zeros((4, 4)))
        close(f'sign_ramp_{index}_preserves_attached_code',
              stages[index].u(stages[index].duration/2) @ P, P)

    integ = nominal_integrals(stages)
    close('physical_stage_product_matches_composite', integ['endpoint'], composite())
    close('common_gain_response_integral_on_code', integ['dc'] @ P, np.zeros((4, 4)))
    close('exposure_is_exactly_three_primitive_loops', integ['exposure'],
          3*coeff['exposure'])
    bounded('zero_projected_energy_throughout_nominal_composite',
            integ['parallel_error'], 2e-13)
    off_resonant = composite(error=0.023)
    coarse = ode(stages, error=0.023, steps_per_time=60)
    fine = ode(stages, error=0.023, steps_per_time=120)
    coarse_error = np.linalg.norm(coarse-off_resonant)
    fine_error = np.linalg.norm(fine-off_resonant)
    bounded('independent_physical_ode_error', fine_error, 5e-9)
    bounded('rk4_step_doubling_ratio', fine_error/coarse_error, 0.08)
    bounded('independent_ode_unitarity', np.linalg.norm(fine.conj().T @ fine-np.eye(4)),
            2e-10)
    delta = 0.007
    detuned = ode(stages, error=0.023, detuning=delta, steps_per_time=120)
    phase_target = off_resonant.copy()
    phase_target[3, 3] = np.exp(-1j*delta*total_time)
    close('reference_phase_factorization_with_finite_ramps', detuned, phase_target, 5e-9)

    # Independent check of the response to drift across the complete waveform.
    drift = 1e-4
    drifting = ode(stages, gain_shape=lambda t: drift*(t-total_time/2)/total_time,
                   steps_per_time=120)
    predicted = composite() @ (np.eye(4)-1j*drift*integ['linear'])
    bounded('linear_drift_response_remainder', np.linalg.norm((drifting-predicted) @ P),
            2e-6)
    expected_drift_leak = np.linalg.norm(integ['linear'][:2, 2])**2
    close('linear_drift_leakage_coefficient',
          unitary_metrics(drifting)['leakage']/drift**2, expected_drift_leak, 3e-4)

    no_loss = ode(stages, steps_per_time=120)
    loss_results = []
    for rate in (2e-6, 1e-6):
        damped = ode(stages, loss=rate, steps_per_time=120)
        slope = (average_infidelity(damped)-average_infidelity(no_loss))/rate
        loss_results.append(slope)
    close('absorbing_loss_average_infidelity_slope', loss_results[-1],
          3*coeff['exposure']/2, 4e-4)
    bounded('absorbing_loss_first_order_convergence',
            abs(loss_results[-1]-3*coeff['exposure']/2)
            /abs(loss_results[0]-3*coeff['exposure']/2), 0.6)

    bare = unitary_metrics(primitive(0.01))
    comp = unitary_metrics(composite(0.01))
    advantage = bare['average_infidelity']-comp['average_infidelity']
    example = dict(common_gain_error=0.01, return_index=1, dimensionless_ramp_time=0.2,
                   primitive=bare, composite=comp,
                   leakage_reduction_factor=bare['leakage']/comp['leakage'],
                   average_infidelity_reduction_factor=bare['average_infidelity']
                   /comp['average_infidelity'],
                   dimensionless_primitive_time=coeff['primitive_time'],
                   dimensionless_composite_time=total_time,
                   dimensionless_primitive_exposure=coeff['exposure'],
                   dimensionless_composite_exposure=integ['exposure'],
                   nominal_peak_bright=coeff['peak_bright'],
                   gain_drift_leakage_coefficient=expected_drift_leak,
                   leading_loss_only_crossover=advantage/coeff['exposure'],
                   leading_reference_noise_only_crossover=np.sqrt(
                       6*advantage/(total_time**2-coeff['primitive_time']**2)))
    return dict(symbolic_checks=symbolic, numerical_checks=numerical,
                symbolic_count=len(symbolic), numerical_count=len(numerical),
                coefficients_n1=coeff, example=example,
                scope='Exact symbolic identities and independent finite-waveform probes; '
                      'noise crossovers are leading weak-noise estimates, not hardware data.')


def plot(path, report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    errors = np.geomspace(2e-4, 0.04, 180)
    bare = [unitary_metrics(primitive(e)) for e in errors]
    comp = [unitary_metrics(composite(e)) for e in errors]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), layout='constrained')
    colors = ('#6a7485', '#146c94')
    for data, color, label in zip((bare, comp), colors, ('Original loop', 'Composite')):
        axes[0].loglog(errors, [x['leakage'] for x in data], color=color, label=label)
        axes[1].loglog(errors, [x['average_infidelity'] for x in data], color=color,
                       label=label)
    axes[0].set(ylabel='Leakage probability', title='Correlated gain error')
    axes[1].set(ylabel='Average logical infidelity', title='Phase and leakage included')
    for ax in axes[:2]:
        ax.set_xlabel('Common fractional gain error')
        ax.grid(True, which='both', alpha=0.17)
        ax.legend(frameon=False)
    ex = report['example']
    max_loss = ex['leading_loss_only_crossover']
    max_phase = ex['leading_reference_noise_only_crossover']
    loss = np.linspace(0, max_loss, 200)
    phase = max_phase*np.sqrt(np.maximum(0, 1-loss/max_loss))
    axes[2].fill_between(loss, 0, phase, color=colors[1], alpha=0.2)
    axes[2].plot(loss, phase, color=colors[1])
    axes[2].set(xlabel=r'Bright loss rate $\kappa/a$',
                ylabel=r'Residual reference noise $\sigma_\Delta/a$',
                title='Composite wins below the boundary',
                xlim=(0, 1.08*max_loss), ylim=(0, 1.08*max_phase))
    axes[2].ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
    axes[2].grid(alpha=0.17)
    axes[2].text(0.04, 0.07, '1% fixed gain error; first return\n'
                 'Each sign ramp: 0.2/a\nLeading independent weak-noise model',
                 transform=axes[2].transAxes, fontsize=9)
    fig.suptitle('Real reflection loop: exact composite response and its resource boundary',
                 fontsize=14)
    fig.savefig(path, dpi=170)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', help='Write the verification receipt to this path')
    parser.add_argument('--plot', help='Write a performance figure to this path')
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

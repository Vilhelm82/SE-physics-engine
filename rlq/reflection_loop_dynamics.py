#!/usr/bin/env python3
"""Finite-duration reflection loop: exact identities and independent ODE probes.

Run with Python, NumPy and SymPy. JSON is printed; --plot PATH optionally writes
an error-curve figure using Matplotlib. Hamiltonians use angular-frequency units.
"""
import argparse
import json
import numpy as np
import sympy as sp

ANGLE = np.pi / 2
J = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]], complex)
G1 = np.array([[0, 0, 1j], [0, 0, 0], [-1j, 0, 0]], complex)
G2 = np.array([[0, 0, 0], [0, 0, -1j], [0, 1j, 0]], complex)
GS = (G1, G2, G1)
R = np.array([[0, 1, 0], [1, 0, 0], [0, 0, -1]], complex)
P = np.diag([0, 0, 1, 1])
Q = np.eye(4) - P


def frame(arc, angle):
    c, s = np.cos(angle), np.sin(angle)
    if arc == 0:
        return np.array([[c, 0, -s], [0, 1, 0], [s, 0, c]])
    if arc == 1:
        return np.array([[0, s, -c], [0, c, s], [1, 0, 0]])
    return np.array([[0, 1, 0], [s, 0, c], [c, 0, -s]])


def embed(matrix):
    out = np.eye(4, dtype=complex)
    out[:3, :3] = matrix
    return out


def segment(arc, duration, gap=1.0, elapsed=None):
    """Propagator in the continuous (p,q,D) frame, including both endpoints."""
    t = duration if elapsed is None else elapsed
    speed = ANGLE / duration
    omega = np.hypot(gap, speed)
    khat = (gap * J + speed * GS[arc]) / omega
    phase = omega * t
    return (np.eye(3) - 1j * np.sin(phase) * khat
            - 2 * np.sin(phase / 2)**2 * (khat @ khat))


def moving_product(durations, gaps=(1.0, 1.0, 1.0)):
    out = np.eye(3, dtype=complex)
    for arc, (duration, gap) in enumerate(zip(durations, gaps)):
        out = segment(arc, duration, gap) @ out
    return out


def exact_loop(durations, gaps=(1.0, 1.0, 1.0)):
    return embed(R @ moving_product(durations, gaps))


def nominal_at(arc, local_time, duration):
    """Full laboratory propagator, valid also away from resonance."""
    before = np.eye(3, dtype=complex)
    for previous in range(arc):
        before = segment(previous, duration) @ before
    return embed(frame(arc, ANGLE * local_time / duration)
                 @ segment(arc, duration, elapsed=local_time) @ before)


def hamiltonian(arc, local_time, duration, gap=1.0):
    f = frame(arc, ANGLE * local_time / duration)
    h = np.zeros((4, 4), complex)
    h[:3, :3] = gap * f @ J @ f.T
    return h


def ode_loop(durations, gaps=(1.0, 1.0, 1.0), perturbation=None, steps=800):
    """Independent RK4 integration of laboratory Schrödinger evolution.

    Does not use the moving-frame propagator or the resonance condition.
    The unitary defect and convergence on step doubling are checked below.
    """
    u = np.eye(4, dtype=complex)
    offset = 0.0
    for arc, (duration, gap) in enumerate(zip(durations, gaps)):
        dt = duration / steps

        def rhs(t, state):
            h = hamiltonian(arc, t, duration, gap)
            if perturbation is not None:
                h = h + perturbation(offset + t, arc, t)
            return -1j * h @ state

        for j in range(steps):
            t = j * dt
            k1 = rhs(t, u)
            k2 = rhs(t + dt / 2, u + dt * k1 / 2)
            k3 = rhs(t + dt / 2, u + dt * k2 / 2)
            k4 = rhs(t + dt, u + dt * k3)
            u = u + dt * (k1 + 2*k2 + 2*k3 + k4) / 6
        offset += duration
    return u


def interaction_integral(duration, perturbation, nodes=96):
    """Gauss-Legendre integral of U0(t)^dagger V(t) U0(t)."""
    points, weights = np.polynomial.legendre.leggauss(nodes)
    result = np.zeros((4, 4), complex)
    for arc in range(3):
        for point, weight in zip(points, weights):
            t = duration * (point + 1) / 2
            u = nominal_at(arc, t, duration)
            v = perturbation(arc * duration + t, arc, t)
            result += weight * duration / 2 * (u.conj().T @ v @ u)
    return result


def run_checks():
    symbolic, numerical = {}, {}

    def exact(name, actual, expected=0):
        if isinstance(actual, sp.MatrixBase):
            expected = sp.zeros(*actual.shape) if expected == 0 else expected
            entries = list(actual - expected)
        else:
            entries = [actual - expected]
        for entry in entries:
            assert sp.simplify(sp.trigsimp(sp.expand(entry))) == 0, (name, entry)
        assert name not in symbolic
        symbolic[name] = True

    def close(name, actual, expected, tolerance=1e-10):
        error = float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))
        assert error <= tolerance, (name, error, tolerance)
        numerical[name] = {'error': error, 'tolerance': tolerance}

    def bounded(name, actual, bound):
        assert actual <= bound, (name, actual, bound)
        numerical[name] = {'value': float(actual), 'upper_bound': float(bound)}

    x, a, w, t = sp.symbols('x a w t', real=True)
    sj = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    sg1 = sp.Matrix([[0, 0, sp.I], [0, 0, 0], [-sp.I, 0, 0]])
    sg2 = sp.Matrix([[0, 0, 0], [0, 0, -sp.I], [0, sp.I, 0]])
    c, s = sp.cos(x), sp.sin(x)
    fs = [sp.Matrix([[c, 0, -s], [0, 1, 0], [s, 0, c]]),
          sp.Matrix([[0, s, -c], [0, c, s], [1, 0, 0]]),
          sp.Matrix([[0, 1, 0], [s, 0, c], [c, 0, -s]])]
    for arc, (f, g) in enumerate(zip(fs, (sg1, sg2, sg1)), 1):
        k = a * sj + w * g
        h = a * f * sj * f.T
        dark = f[:, 2]
        projector = dark * dark.T
        cd = sp.I * w * (projector.diff(x)*projector-projector*projector.diff(x))
        exact(f'frame_{arc}_generator', a*sj-sp.I*w*f.T*f.diff(x), k)
        exact(f'frame_{arc}_hermiticity', k.H, k)
        exact(f'frame_{arc}_minimal_polynomial', k**3, (a*a+w*w)*k)
        exact(f'frame_{arc}_counterdiabatic_tracking', (h+cd)*dark,
              sp.I*w*dark.diff(x))
        exact(f'frame_{arc}_counterdiabatic_frame_cancellation',
              f.T*cd*f, -w*g)
        exact(f'frame_{arc}_cd_square', cd**2,
              w*w*(dark*dark.T+dark.diff(x)*dark.diff(x).T))

    # Trigonometric polynomial proof for the complete equal-duration family.
    u, v, co, si = sp.symbols('u v co si', real=True)
    k1, k2 = u*sj+v*sg1, u*sj+v*sg2
    e1 = sp.eye(3)-sp.I*si*k1+(co-1)*k1*k1
    e2 = sp.eye(3)-sp.I*si*k2+(co-1)*k2*k2
    total = e1*e2*e1
    ideal = sp.groebner([u*u+v*v-1, si*si+co*co-1], u, si, v, co)
    imaginary_z = sp.im(sp.expand(total[2, 2]))
    expected_im = 4*si*u*v*v*(co-1)*(1+v*v*(co-1))
    exact('equal_family_imaginary_return_amplitude',
          ideal.reduce(sp.expand(imaginary_z-expected_im))[1])
    odd_z = sp.re(sp.expand(total[2, 2])).subs({co: -1, si: 0})
    exact('odd_half_cycle_exclusion',
          ideal.reduce(sp.expand(odd_z-1+2*v*v*(4*v**4-8*v*v+5)))[1])
    exact('even_full_cycle_return', total.subs({co: 1, si: 0}), sp.eye(3))
    exact('single_arc_distance_squared', 2*(1-e1[2, 2]), 2*v*v*(1-co))
    exact('single_arc_one_zero_mean_energy', (e1[:, 2].H*sj*e1[:, 2])[0])
    exact('single_arc_two_zero_mean_energy', (e2[:, 2].H*sj*e2[:, 2])[0])
    bright = 1-(1-v*v*(1-sp.cos(x)))**2
    exact('cycle_average_bright_population', sp.integrate(bright, (x, 0, 2*sp.pi))/(2*sp.pi),
          2*v*v-sp.Rational(3, 2)*v**4)
    chi_p = sp.sin(4*x)/4
    chi_d = sp.Rational(15, 16)+sp.cos(4*x)/16
    edge_readout = (sp.sin(2*x)*(chi_p**2-chi_d**2)
                    +2*sp.cos(2*x)*chi_p*chi_d)
    exact('first_resonance_uniform_edge_offset_integral',
          sp.integrate(sp.expand_trig(edge_readout), (x, 0, sp.pi/2)), -sp.Rational(59, 120))

    # Derivatives at a resonant lambda, with lambda/phi=u, theta/phi=v.
    phi, lam, th = sp.symbols('phi lam th', positive=True)
    la = lam*sj+th*sg1
    radius = sp.sqrt(lam*lam+th*th)
    general_e = (sp.eye(3)-sp.I*sp.sin(radius)/radius*la
                 +(sp.cos(radius)-1)/radius**2*la**2)
    derivative = -sp.I*(lam/phi**2)*la
    second = (-sp.I*((1-3*lam*lam/phi**2)*la/phi**2
                    +2*lam*sj/phi**2) - lam*lam*la**2/phi**4)
    d = sp.Matrix([0, 0, 1])
    resonance_subs = {sp.sin(radius): 0, sp.cos(radius): 1}
    exact('resonant_exponential_first_derivative', general_e.diff(lam).subs(resonance_subs),
          derivative.subs(phi, radius))
    exact('resonant_exponential_second_derivative', general_e.diff(lam, 2).subs(resonance_subs),
          second.subs(phi, radius))
    exact('opposite_outer_second_order_amplitude',
          lam**2*(second-derivative**2)*d,
          sp.Matrix([th*(lam/phi)**2*(1-3*(lam/phi)**2), 0, 0]))
    eps1, eps2, eps3, beta = sp.symbols('eps1 eps2 eps3 beta', real=True)
    leakage_vector = -sp.I*beta*((eps1+eps3)*sg1+eps2*sg2)*d
    exact('independent_gain_error_vector', leakage_vector,
          sp.Matrix([beta*(eps1+eps3), -beta*eps2, 0]))
    exact('gain_error_quadratic_form', (leakage_vector.H*leakage_vector)[0],
          beta**2*((eps1+eps3)**2+eps2**2))
    rot = sp.Matrix([[sp.cos(x), -sp.sin(x)], [sp.sin(x), sp.cos(x)]])
    zgate = sp.diag(-1, 1)
    rotated = rot*zgate*rot.T
    exact('rotated_gate_retains_parity', rotated.det(), -1)
    exact('rotated_gate_error_norm_squared', (rotated-zgate).T*(rotated-zgate),
          4*sp.sin(x)**2*sp.eye(2))
    exact('spectator_phase_interference', (2+2*sp.cos(x))/4, sp.cos(x/2)**2)

    tau = ANGLE*np.sqrt(15)
    target = embed(R)
    for n in (1, 2, 4):
        duration = ANGLE*np.sqrt(16*n*n-1)
        close(f'exact_resonance_n{n}', exact_loop([duration]*3), target)
    mixed_durations = [ANGLE*np.sqrt(16*n*n-1)/gap
                       for n, gap in zip((1, 2, 3), (1.0, 0.8, 1.0))]
    close('unequal_resonance_family', exact_loop(mixed_durations, (1, .8, 1)), target)
    ode_coarse = ode_loop([tau]*3, steps=400)
    ode_fine = ode_loop([tau]*3, steps=800)
    coarse_error = np.linalg.norm(ode_coarse-target, 2)
    fine_error = np.linalg.norm(ode_fine-target, 2)
    bounded('laboratory_ode_step_doubling', fine_error/coarse_error, .07)
    close('laboratory_ode_exact_gate', ode_fine, target, 1e-9)
    close('laboratory_ode_unitarity', ode_fine.conj().T@ode_fine, np.eye(4), 1e-10)
    close('laboratory_ode_unequal_off_resonance',
          ode_loop([2.1, 4.3, 3.2], (.9, 1.2, .7), steps=500),
          exact_loop([2.1, 4.3, 3.2], (.9, 1.2, .7)), 2e-9)
    expected_common = 5*(ANGLE*15/16)**2
    expected_opposite = (ANGLE*15/16*(1-45/16))**2
    small = 1e-4
    common = moving_product([tau]*3, (1+small,)*3)
    opposite = moving_product([tau]*3, (1+small, 1, 1-small))
    close('common_gain_quadratic_coefficient', np.sum(abs(common[:2, 2])**2)/small**2,
          expected_common, .01)
    close('opposite_outer_quartic_coefficient', np.sum(abs(opposite[:2, 2])**2)/small**4,
          expected_opposite, .01)
    eps = np.array([.4, -.7, 1.2])*1e-5
    gain_result = moving_product([tau]*3, 1+eps)
    close('independent_gain_linear_response', gain_result[:2, 2],
          ANGLE*15/16*np.array([eps[0]+eps[2], -eps[1]]), 2e-8)
    for duration in (.2, 1.3, 4.7, 10.2, 23.0):
        product = moving_product([duration]*3)
        phase = np.hypot(duration, ANGLE)
        vnum = ANGLE/phase
        delta = np.linalg.norm(product[:, 2]-np.array([0, 0, 1]))
        bounded(f'finite_time_bound_{duration}', delta,
                min(2, 6*vnum*abs(np.sin(phase/2)))+1e-12)
    points, weights = np.polynomial.legendre.leggauss(96)
    exposure = 0.0
    occupations = np.zeros(3)
    for arc in range(3):
        for point, weight in zip(points, weights):
            local_time = tau*(point+1)/2
            state = nominal_at(arc, local_time, tau)[:, 2]
            dark = frame(arc, ANGLE*local_time/tau)[:, 2]
            exposure += weight*tau/2*(1-abs(dark.conj()@state[:3])**2)
            occupations += weight*tau/2*abs(state[:3])**2
    close('transient_bright_exposure', exposure, 3*tau*(2/16-1.5/256))
    close('equal_integrated_site_occupations', occupations, np.full(3, tau))
    detuning = np.diag([.2, -.7, .5, .3])
    perturbation = lambda total_time, arc, local_time: detuning
    m = interaction_integral(tau, perturbation)
    close('interaction_integral_hermiticity', m.conj().T, m)
    close('traceless_active_detuning_phase_cancellation', m[2, 2], 0)
    close('spectator_detuning_filter', m[3, 3], .3*3*tau)
    uniform_edges = np.zeros((4, 4), complex)
    uniform_edges[:3, :3] = np.ones((3, 3))-np.eye(3)
    edge_filter = interaction_integral(tau, lambda *args: uniform_edges)
    close('uniform_edge_offset_finite_time_phase', edge_filter[2, 2]/(3*tau),
          -59/(60*np.pi))
    for n in (2, 3):
        duration = ANGLE*np.sqrt(16*n*n-1)
        k = 4*n
        coefficient = -(2*k**4-23*k*k+33)/(2*(k*k-4)*(k*k-1))/ANGLE
        edge_filter_n = interaction_integral(duration, lambda *args: uniform_edges)
        close(f'uniform_edge_offset_phase_n{n}', edge_filter_n[2, 2]/(3*duration), coefficient)
    strength = 1e-3
    perturbed = ode_loop([tau]*3, perturbation=lambda *args: strength*detuning, steps=800)
    kappa = strength*3*tau*np.linalg.norm(detuning, 2)
    bounded('duhamel_perturbation_bound', np.linalg.norm(perturbed-target, 2), kappa+1e-9)
    bounded('interaction_expansion_remainder',
            np.linalg.norm(target.conj().T@perturbed-np.eye(4)+1j*strength*m, 2),
            kappa*kappa/2+1e-9)
    alpha = .13
    g = np.eye(4)
    g[2:, 2:] = [[np.cos(alpha), -np.sin(alpha)], [np.sin(alpha), np.cos(alpha)]]
    altered = ode_loop([tau]*3,
                       perturbation=lambda total, arc, local:
                       g@hamiltonian(arc, local, tau)@g.T-hamiltonian(arc, local, tau),
                       steps=800)
    close('laboratory_ode_rotated_reflection', altered, g@target@g.T, 1e-9)
    close('rotated_reflection_axis_error', np.linalg.norm((altered-target)@P, 2),
          2*abs(np.sin(alpha)), 1e-9)

    rows = []
    for n in range(1, 5):
        vn2 = 1/(16*n*n)
        duration = ANGLE*np.sqrt(16*n*n-1)
        product = moving_product([duration]*3, (1.01,)*3)
        rows.append({'n': n, 'gap_times_total_duration': 3*duration,
                     'peak_bright_population': 4*vn2*(1-vn2),
                     'gap_times_integrated_bright_population':
                     3*duration*(2*vn2-1.5*vn2*vn2),
                     'endpoint_leakage_at_one_percent_common_gain_error':
                     float(np.sum(abs(product[:2, 2])**2))})
    return {'symbolic_check_count': len(symbolic), 'symbolic_checks': symbolic,
            'numerical_check_count': len(numerical), 'numerical_checks': numerical,
            'return_family': rows,
            'scope': 'Exact closed-system equal-speed loop and declared perturbations; numerical ODE probes separately identified.',
            'novelty': 'Mathematical continuation of the native family; historical originality not established.'}


def make_plot(path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    tau = ANGLE*np.sqrt(15)
    times = np.linspace(3, 60, 2500)
    errors = []
    for total in times:
        z = moving_product([total/3]*3)[2, 2]
        errors.append(max(1e-16, (2-z.real-abs(z)**2)/3))
    gain_errors = np.geomspace(1e-4, .05, 200)
    common, opposite = [], []
    for eps in gain_errors:
        common.append(np.sum(abs(moving_product([tau]*3, (1+eps,)*3)[:2, 2])**2))
        opposite.append(np.sum(abs(moving_product([tau]*3, (1+eps, 1, 1-eps))[:2, 2])**2))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout='constrained')
    axes[0].semilogy(times, errors, color='#195c9b', lw=1.4)
    for n in (1, 2, 3):
        total = 3*ANGLE*np.sqrt(16*n*n-1)
        axes[0].axvline(total, color='#65717b', ls='--', lw=.9)
        axes[0].text(total+.5, .35, f'n={n}', color='#37434d')
    axes[0].set(xlabel=r'Total duration $aT$', ylabel='Average logical infidelity',
                title='Finite-time response; dashed lines are exact returns', ylim=(1e-9, 1))
    axes[1].loglog(gain_errors, common, color='#a54229', label='Common gain error: quadratic')
    axes[1].loglog(gain_errors, opposite, color='#196b50', label='Opposite outer errors: quartic')
    axes[1].loglog(gain_errors, 5*(ANGLE*15/16)**2*gain_errors**2, ':', color='#a54229', alpha=.6)
    axes[1].loglog(gain_errors, (ANGLE*15/16*(1-45/16))**2*gain_errors**4, ':', color='#196b50', alpha=.6)
    axes[1].set(xlabel=r'Fractional amplitude error $\epsilon$', ylabel='Endpoint leakage probability',
                title='First exact return; dotted lines are derived leading terms')
    axes[1].legend(frameon=False, loc='upper left', fontsize=9)
    for axis in axes:
        axis.grid(True, alpha=.15, which='both')
        axis.spines[['top', 'right']].set_visible(False)
    fig.savefig(path, dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot', help='Optional output image path')
    args = parser.parse_args()
    result = run_checks()
    if args.plot:
        make_plot(args.plot)
    print(json.dumps(result, indent=2))

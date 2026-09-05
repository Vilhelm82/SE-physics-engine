#!/usr/bin/env python3
"""Exact feasibility witnesses for the quantum prediction-target evaluation.

These are mathematical probes, not a simulator or a novelty certificate.
Run directly to print the JSON receipt; no files are modified.
"""
import json
import sympy as s


def run_checks():
    checks = {}

    def check(name, actual, expected=0):
        if isinstance(actual, s.MatrixBase):
            expected = s.zeros(*actual.shape) if expected == 0 else expected
            entries = list(actual - expected)
        else:
            entries = [s.sympify(actual) - expected]
        assert name not in checks, name
        for entry in entries:
            value = s.simplify(s.trigsimp(s.expand(entry)))
            if value != 0:
                raise AssertionError((name, value))
        checks[name] = True

    # Three real, fixed-gap arcs; an extra zero mode is the phase reference.
    q = s.symbols('q', real=True)
    e = [s.eye(4)[:, i] for i in range(4)]
    c, sn = s.cos(q), s.sin(q)
    paths = [(c*e[0]+sn*e[2], e[1], -sn*e[0]+c*e[2]),
             (e[2], c*e[1]+sn*e[0], -c*e[0]+sn*e[1]),
             (c*e[2]+sn*e[1], e[0], c*e[1]-sn*e[2])]
    initial_h = e[0]*e[1].T+e[1]*e[0].T
    last_h, last_dark = initial_h, e[2]
    for i, (p, b, dark) in enumerate(paths, 1):
        H = p*b.T+b*p.T
        P = s.eye(4)-H*H
        check(f'arc_{i}_bright_frame_orthonormal', p.row_join(b).T*p.row_join(b), s.eye(2))
        check(f'arc_{i}_zero_diagonal', s.Matrix(H.diagonal()))
        check(f'arc_{i}_fixed_spectrum_polynomial', H**3, H)
        check(f'arc_{i}_dark_projector', P*P, P)
        check(f'arc_{i}_dark_eigenvector', H*dark)
        check(f'arc_{i}_spectator_dark', H*e[3])
        check(f'arc_{i}_parallel_dark_frame', dark.row_join(e[3]).T*dark.row_join(e[3]).diff(q))
        check(f'arc_{i}_hamiltonian_join', H.subs(q, 0), last_h)
        check(f'arc_{i}_dark_join', dark.subs(q, 0), last_dark)
        last_h = H.subs(q, s.pi/2)
        last_dark = dark.subs(q, s.pi/2)
    check('loop_returns_hamiltonian', last_h, initial_h)
    check('loop_reverses_one_dark_vector', last_dark, -e[2])
    initial_frame = e[2].row_join(e[3])
    gate = initial_frame.T*last_dark.row_join(e[3])
    check('dark_gate_is_reflection', gate, s.diag(-1, 1))
    check('dark_gate_orientation', gate.det(), -1)
    check('superposition_interference_flips',
          ((s.Matrix([1, 1]).T*gate*s.Matrix([1, 1]))[0]/2)**2, 0)
    theta, phi = s.symbols('theta phi', real=True)
    n = s.Matrix([s.cos(theta), s.sin(theta)*s.cos(phi), s.sin(theta)*s.sin(phi)])
    sphere_frame = n.diff(theta).row_join(s.Matrix([0, -s.sin(phi), s.cos(phi)]))
    check('sphere_frame_orthonormal', sphere_frame.T*sphere_frame, s.eye(2))
    check('sphere_latitude_connection', sphere_frame.T*sphere_frame.diff(phi),
          s.Matrix([[0, -s.cos(theta)], [s.cos(theta), 0]]))

    # An actual four-spin XY Hamiltonian, independently restricted by occupation.
    x, y, z = s.symbols('x y z', real=True)
    raising = s.Matrix([[0, 0], [1, 0]])

    def on_site(operator, site):
        factors = [operator if j == site else s.eye(2) for j in range(4)]
        return s.kronecker_product(*factors)

    plus = [on_site(raising, j) for j in range(4)]
    Hfull = sum((g*(plus[0]*plus[j].T+plus[0].T*plus[j])
                 for j, g in enumerate((x, y, z), 1)), s.zeros(16))
    H1 = Hfull.extract([8, 4, 2, 1], [8, 4, 2, 1])
    H2 = Hfull.extract([12, 10, 9, 6, 5, 3], [12, 10, 9, 6, 5, 3])
    B = s.Matrix([[y, z, 0], [x, 0, z], [0, x, y]])
    check('one_excitation_rank_two_identity', H1**3, (x*x+y*y+z*z)*H1)
    check('two_excitation_block_from_spin_hamiltonian', H2,
          s.zeros(3).row_join(B).col_join(B.T.row_join(s.zeros(3))))
    check('two_excitation_off_block_determinant', B.det(), -2*x*y*z)
    check('two_excitation_full_determinant', H2.det(), -4*x*x*y*y*z*z)
    left, right = s.Matrix([x, -y, 0]), s.Matrix([0, y, -x])
    check('closing_mode_left_null', left.T*B.subs(z, 0))
    check('closing_mode_right_null', B.subs(z, 0)*right)
    check('closing_mode_linear_splitting', (left.T*B.diff(z)*right)[0]/(x*x+y*y),
          2*x*y/(x*x+y*y))

    # A CPTP critical crossover with approach-dependent dark state.
    r, time = s.symbols('r time', nonnegative=True)
    vac, one, two = [s.eye(3)[:, i] for i in range(3)]
    bright = s.cos(theta)*one+s.sin(theta)*two
    dark = -s.sin(theta)*one+s.cos(theta)*two
    PB, PD, P0 = bright*bright.T, dark*dark.T, vac*vac.T
    jump = r*vac*bright.T
    E = s.exp(-r*r*time)
    rho = (s.cos(theta)**2*(E*PB+(1-E)*P0)+s.sin(theta)**2*PD
           -s.cos(theta)*s.sin(theta)*s.exp(-r*r*time/2)*(bright*dark.T+dark*bright.T))
    generator = jump*rho*jump.T-(jump.T*jump*rho+rho*jump.T*jump)/2
    check('critical_crossover_initial_state', rho.subs(time, 0), one*one.T)
    check('critical_crossover_master_equation', rho.diff(time), generator)
    check('critical_crossover_trace', s.trace(rho), 1)
    K0 = P0+PD+s.exp(-r*r*time/2)*PB
    V = vac*bright.T
    check('critical_crossover_kraus_completeness', K0.T*K0+(1-E)*V.T*V, s.eye(3))
    check('critical_crossover_kraus_state', K0*one*one.T*K0.T+(1-E)*V*one*one.T*V.T, rho)
    check('critical_crossover_excitation_readout', s.trace((s.eye(3)-P0)*rho),
          s.sin(theta)**2+s.cos(theta)**2*E)
    check('critical_zero_generator_retains_initial_state', rho.subs(r, 0), one*one.T)

    # Operator growth need not introduce a singular full-state Jacobian.
    X, Z = s.Matrix([[0, 1], [1, 0]]), s.diag(1, -1)
    CNOT = s.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]])
    W, probe = s.kronecker_product(X, s.eye(2)), s.kronecker_product(s.eye(2), Z)
    check('initial_local_probe_commutes', W*probe-probe*W)
    evolved = CNOT.T*W*CNOT
    check('operator_support_spreads', evolved, s.kronecker_product(X, X))
    commutator = evolved*probe-probe*evolved
    check('commutator_response_appears', s.trace(commutator.H*commutator)/4, 4)
    full_jacobian = s.kronecker_product(CNOT, CNOT)
    check('full_density_evolution_jacobian_regular', full_jacobian.T*full_jacobian, s.eye(16))
    return {'check_count': len(checks), 'checks': checks,
            'scope': 'Exact feasibility witnesses; ideal adiabatic holonomy, finite-spin sectors, CPTP crossover, operator response',
            'novelty': 'Not established by these checks; see the prior-art comparison in the evaluation.'}


if __name__ == '__main__':
    print(json.dumps(run_checks(), indent=2))

#!/usr/bin/env python3
"""Exact predictive-state closure and native constitutive receipts.

Proof: docs/2026-09-05-cella-predictive-state.md.
Running this file prints JSON; it does not modify project or corpus files.
"""

import json

import sympy as s


def lie(function, field, coordinates):
    """Directional derivative, with the input held fixed."""
    return s.expand(sum(s.diff(function, q) * v
                        for q, v in zip(coordinates, field)))


def predictive_pair_ideal(coordinates, fields, outputs, max_rounds=12):
    """Compute the exact polynomial indistinguishability ideal over QQ.

    Parameters must be specialized rational coefficients or explicit state
    coordinates. No generic-parameter inversion or radicalization is hidden.
    A budget exhaustion raises; only zero derivative remainders certify closure.
    The returned real zero set must still be intersected with admissible states.
    """
    coordinates = tuple(coordinates)
    fields = tuple(tuple(field) for field in fields)
    if not coordinates or len(set(coordinates)) != len(coordinates):
        raise ValueError("state coordinates must be nonempty and distinct")
    if any(not isinstance(q, s.Symbol) for q in coordinates):
        raise ValueError("state coordinates must be symbols")
    if any(len(field) != len(coordinates) for field in fields):
        raise ValueError("field dimensions disagree with the state")
    if max_rounds < 0:
        raise ValueError("max_rounds must be nonnegative")
    twins = tuple(s.Symbol(f"{q.name}__paired", real=True) for q in coordinates)
    if set(twins) & set(coordinates) or len(set(twins)) != len(twins):
        raise ValueError("paired coordinate names collide")
    variables = coordinates + twins
    swap = dict(zip(coordinates, twins))
    fields = tuple(tuple(map(s.sympify, field)) for field in fields)
    outputs = tuple(map(s.sympify, outputs))
    for polynomial in outputs + tuple(v for field in fields for v in field):
        s.Poly(polynomial, *coordinates, domain=s.QQ)
    paired_fields = [field + tuple(v.xreplace(swap) for v in field)
                     for field in fields]
    generators = [h - h.xreplace(swap) for h in outputs]
    basis = s.groebner(generators, *variables, domain=s.QQ)
    history = [len(basis.polys)]
    rounds = 0
    while True:
        residuals = []
        for polynomial in basis.polys:
            for field in paired_fields:
                derivative = lie(polynomial.as_expr(), field, variables)
                remainder = s.expand(basis.reduce(derivative)[1])
                if remainder != 0 and remainder not in residuals:
                    residuals.append(remainder)
        if not residuals:
            return {
                "basis": basis,
                "twins": twins,
                "paired_fields": paired_fields,
                "rounds": rounds,
                "basis_sizes": history,
            }
        if rounds >= max_rounds:
            raise RuntimeError("iteration budget reached before an invariance certificate")
        basis = s.groebner([p.as_expr() for p in basis.polys] + residuals,
                           *variables, domain=s.QQ)
        rounds += 1
        history.append(len(basis.polys))


def linear_observation_rows(matrices, output):
    """Row span dual to the largest common invariant observation-null space."""
    output = s.Matrix(output)
    width = output.cols
    matrices = [s.Matrix(a) for a in matrices]
    if any(a.shape != (width, width) for a in matrices):
        raise ValueError("material matrix dimensions disagree")

    def row_basis(rows):
        rref, pivots = rows.rref()
        return rref[:len(pivots), :]

    rows = row_basis(output)
    while True:
        augmented = rows
        for a in matrices:
            augmented = augmented.col_join(rows * a)
        updated = row_basis(augmented)
        if updated.rows == rows.rows:
            return updated
        rows = updated


def run_checks():
    checks = {}
    certificates = {}

    def check(name, actual, expected=0):
        if isinstance(actual, s.MatrixBase):
            target = s.zeros(*actual.shape) if expected == 0 else expected
            residuals = list(actual - target)
        else:
            residuals = [s.sympify(actual) - s.sympify(expected)]
        for residual in residuals:
            reduced = s.simplify(s.cancel(s.expand(residual)))
            if reduced != 0:
                raise AssertionError((name, reduced))
        checks[name] = {"verified": True}

    def require(name, condition):
        if not condition:
            raise AssertionError(name)
        checks[name] = {"verified": True}

    def certify(name, result, expected_generators):
        basis = result["basis"]
        expected = s.groebner(expected_generators, *basis.gens, domain=s.QQ)
        require(f"{name}_exact_ideal",
                all(basis.reduce(p.as_expr())[1] == 0 for p in expected.polys)
                and all(expected.reduce(p.as_expr())[1] == 0 for p in basis.polys))
        for number, polynomial in enumerate(basis.polys):
            for index, field in enumerate(result["paired_fields"]):
                check(f"{name}_invariance_{number}_{index}",
                      basis.reduce(lie(polynomial.as_expr(), field, basis.gens))[1])
        certificates[name] = {
            "rounds": result["rounds"],
            "basis_sizes": result["basis_sizes"],
            "basis": [str(p.as_expr()) for p in basis.polys],
        }

    xi, d = s.symbols("xi d", real=True)
    mu, beta, nu = s.symbols("mu beta nu", real=True)
    a = s.Symbol("a", real=True)
    p = 1 + a ** 2
    coordinates = (xi, d)
    F = (1 - xi ** 2) / p
    x = d * xi
    Vx = s.Matrix([mu * d + beta * xi, beta * d + nu * xi])
    Vd = s.Matrix([beta, nu])
    Lx = lambda f: lie(f, Vx, coordinates)
    Ld = lambda f: lie(f, Vd, coordinates)
    mobility = s.Matrix([[mu, beta], [beta, nu]])
    observation = s.Matrix([x, d])
    jacobian = observation.jacobian(coordinates)
    ex, ed = s.symbols("ex ed", real=True)
    effort = s.Matrix([ex, ed])
    flow = mobility * jacobian.T * effort
    check("original_ports_generate_fields", flow, ex * Vx + ed * Vd)
    check("original_ports_power_pairing",
          (effort.T * jacobian * flow)[0],
          ((jacobian.T * effort).T * mobility * (jacobian.T * effort))[0])
    check("positive_response_schur_factor",
          mobility,
          s.Matrix([[1, 0], [beta / mu, 1]])
          * s.diag(mu, (mu * nu - beta ** 2) / mu)
          * s.Matrix([[1, beta / mu], [0, 1]]))
    check("native_lapse_first_x", Lx(F), -2 * xi * (mu * d + beta * xi) / p)
    check("native_lapse_first_d", Ld(F), -2 * beta * xi / p)
    check("native_input_commutator", Vx.jacobian(coordinates) * Vd,
          s.Matrix([mu * nu + beta ** 2, 2 * beta * nu]))
    mixed = Ld(Lx(F)) - Lx(Ld(F))
    check("native_lapse_mixed_difference", mixed,
          -2 * (mu * nu + beta ** 2) * xi / p)
    reconstructed_xi = s.cancel(-p * mixed / (2 * (mu * nu + beta ** 2)))
    reconstructed_d = s.cancel((Lx(reconstructed_xi) - beta * reconstructed_xi) / mu)
    check("reconstruct_direction_uniformly", reconstructed_xi, xi)
    check("reconstruct_volume_uniformly", reconstructed_d, d)
    reconstruction = s.Matrix([reconstructed_xi, reconstructed_d])
    check("minimal_two_state_jacobian", reconstruction.jacobian(coordinates), s.eye(2))
    check("reconstruction_at_bisector_and_pinch",
          reconstruction.jacobian(coordinates).subs({xi: 0, d: 0}), s.eye(2))
    check("zero_coupling_second_word", Ld(Lx(F)).subs(beta, 0), -2 * mu * nu * xi / p)
    check("zero_coupling_third_word", Lx(Ld(Lx(F))).subs(beta, 0),
          -2 * mu ** 2 * nu * d / p)

    # All directions with transverse response are distinguished at the collapse.
    r11, r12, r21, n = s.symbols("r11 r12 r21 n", real=True)
    general = s.Matrix([[r11, r12], [r21, n]])
    bxi, bd = s.symbols("bxi bd", real=True)
    pushed = jacobian * (s.Matrix([bxi, bd]) + general * jacobian.T * effort)
    check("general_transverse_direction_identity",
          (pushed[0] - xi * pushed[1]).subs(d, 0))
    boundary_response = (jacobian * general * jacobian.T).subs(d, 0)
    check("general_incremental_boundary_matrix", boundary_response,
          n * s.Matrix([[xi ** 2, xi], [xi, 1]]))
    check("boundary_increment_reconstructs_xi",
          s.cancel(boundary_response[0, 1] / boundary_response[1, 1]), xi)
    check("affine_drift_transverse_direction", pushed.subs({d: 0, ex: 0, ed: 0}),
          s.Matrix([xi * bd, bd]))

    # The closed ideal proves that neither a discrete sign nor a first-order
    # rank test suffices. Positive rational fixtures cover both beta strata.
    for name, values in (
        ("native_coupled", {mu: 2, beta: 1, nu: 1, a: 0}),
        ("native_uncoupled", {mu: 2, beta: 0, nu: 3, a: 0}),
        ("native_rank_one", {mu: 1, beta: 1, nu: 1, a: 0}),
    ):
        fields = [Vx.subs(values), Vd.subs(values)]
        result = predictive_pair_ideal(coordinates, fields, [F.subs(values)])
        xi2, d2 = result["twins"]
        certify(name, result, [xi - xi2, d - d2])
        require(f"{name}_initial_horizon_outputs_equal",
                F.subs({xi: 1, a: 0}) == F.subs({xi: -1, a: 0}))
        rejected = [p_.as_expr().subs({xi: 1, xi2: -1, d: 0, d2: 0})
                    for p_ in result["basis"].polys]
        require(f"{name}_horizon_partners_distinguished", any(v != 0 for v in rejected))
    try:
        predictive_pair_ideal(coordinates,
                              [Vx.subs({mu: 2, beta: 0, nu: 3}),
                               Vd.subs({mu: 2, beta: 0, nu: 3})], [1 - xi ** 2],
                              max_rounds=0)
    except RuntimeError:
        require("unfinished_ideal_does_not_self_certify", True)
    else:
        raise AssertionError("iteration budget incorrectly treated as completion")

    # Restricting experiments changes the exact quotient.
    radial = predictive_pair_ideal(coordinates, [s.Matrix([1, 1])], [1 - xi ** 2])
    xi2, d2 = radial["twins"]
    certify("radial_only", radial, [xi - xi2])
    radial_zero = predictive_pair_ideal(coordinates, [s.Matrix([0, 1])], [1 - xi ** 2])
    xi2, d2 = radial_zero["twins"]
    certify("radial_zero_coupling", radial_zero, [xi ** 2 - xi2 ** 2])
    folded = predictive_pair_ideal(coordinates, [(2 * d, 0)], [1 - xi ** 2])
    xi2, d2 = folded["twins"]
    certify("zero_radial_mobility_deck_quotient", folded,
            [xi ** 2 - xi2 ** 2, xi * d - xi2 * d2, d ** 2 - d2 ** 2])
    frozen = predictive_pair_ideal(coordinates, [(0, xi), (0, 1)], [1 - xi ** 2])
    xi2, d2 = frozen["twins"]
    certify("zero_directional_mobility_frozen_lapse", frozen, [xi ** 2 - xi2 ** 2])
    U, V, W = s.symbols("U V W", real=True)
    cone = s.Matrix([xi ** 2, xi * d, d ** 2])
    cone_rates = s.Matrix([lie(h, (mu * d, 0), coordinates) for h in cone])
    check("folded_cone_polynomial_dynamics", cone_rates,
          s.Matrix([2 * mu * xi * d, mu * d ** 2, 0]))
    check("folded_cone_relation_invariant",
          lie(U * W - V ** 2, (2 * mu * V, mu * W, 0), (U, V, W)))
    check("folded_two_coordinate_radius_squared",
          (xi ** 2 - d ** 2) ** 2 + (2 * xi * d) ** 2, (xi ** 2 + d ** 2) ** 2)
    check("folded_two_coordinate_A_dynamics",
          lie(xi ** 2 - d ** 2, (mu * d, 0), coordinates), mu * 2 * xi * d)
    check("folded_two_coordinate_B_dynamics",
          lie(2 * xi * d, (mu * d, 0), coordinates),
          mu * ((xi ** 2 + d ** 2) - (xi ** 2 - d ** 2)))
    fixed_materials = {mu: 2, beta: 1, nu: 1}
    resolved_controls = predictive_pair_ideal(
        coordinates, [mobility[:, i].subs(fixed_materials) for i in range(2)],
        [1 - xi ** 2])
    xi2, d2 = resolved_controls["twins"]
    certify("resolved_efforts_are_different_experiments", resolved_controls, [xi - xi2])

    # Exact two-step horizon witness when beta=0.
    t, elapsed = s.symbols("t elapsed", positive=True)
    after_radial = {d: nu * t, beta: 0}
    check("horizon_plus_after_switch", Lx(F).subs(after_radial).subs(xi, 1),
          -2 * mu * nu * t / p)
    check("horizon_minus_after_switch", Lx(F).subs(after_radial).subs(xi, -1),
          2 * mu * nu * t / p)
    # A complete rational-coefficient solution verifies that the rate witness
    # comes from actual trajectories, not just a formal tangent assignment.
    sx = s.cosh(elapsed) + t * s.sinh(elapsed)
    sd = s.sinh(elapsed) + t * s.cosh(elapsed)
    check("horizon_two_step_exact_xi_flow", s.diff(sx, elapsed), sd)
    check("horizon_two_step_exact_d_flow", s.diff(sd, elapsed), sx)
    check("horizon_two_step_lapse_slope", s.diff(1 - sx ** 2, elapsed).subs(elapsed, 0), -2 * t)
    check("horizon_two_step_matches_initial_direction", sx.subs(elapsed, 0), 1)
    check("horizon_two_step_matches_initial_volume", sd.subs(elapsed, 0), t)

    tau = {xi: -xi, d: -d}
    check("deck_fixes_lapse", F.xreplace(tau), F)
    check("deck_preserves_x_control", Vx.xreplace(tau), -Vx)
    check("deck_reverses_radial_control", Vd.xreplace(tau), Vd)
    check("deck_covariance_requires_input_transform",
          (ex * Vx + ed * Vd).xreplace(tau), -(ex * Vx - ed * Vd))

    # Every unknown constant in the positive material family is recoverable
    # from the single signed-volume output and a finite sequence of inputs.
    material_coordinates = (xi, d, mu, beta, nu)
    material_x = tuple(Vx) + (0, 0, 0)
    material_d = tuple(Vd) + (0, 0, 0)
    Lmx = lambda f: lie(f, material_x, material_coordinates)
    Lmd = lambda f: lie(f, material_d, material_coordinates)
    check("material_volume_first_probe", Lmd(d), nu)
    check("material_volume_second_probe", Lmd(Lmx(d)), 2 * beta * nu)
    check("material_volume_third_probe", Lmd(Lmx(Lmx(d))), nu * (mu * nu + 3 * beta ** 2))
    rn = Lmd(d)
    rb = s.cancel(Lmd(Lmx(d)) / (2 * rn))
    rm = s.cancel((Lmd(Lmx(Lmx(d))) / rn - 3 * rb ** 2) / rn)
    rx = s.cancel((Lmx(d) - rb * d) / rn)
    material_reconstruction = s.Matrix([rx, d, rm, rb, rn])
    check("five_material_state_reconstruction", material_reconstruction,
          s.Matrix(material_coordinates))
    check("five_material_state_minimum", material_reconstruction.jacobian(material_coordinates),
          s.eye(5))
    k = s.Symbol("k", real=True)
    correlated = {mu: 1 + k ** 2, beta: k, nu: 1}
    check("single_material_coordinate_is_recovered", rb.subs(correlated), k)
    check("single_material_parameter_does_not_add_three_states",
          s.Matrix([rx.subs(correlated), d, rb.subs(correlated)]).jacobian((xi, d, k)), s.eye(3))

    # Delayed dependence and a removable variable are independent controls
    # against merely reproducing the pinch-specific reconstruction.
    y, z, w = s.symbols("y z w", real=True)
    chain = predictive_pair_ideal((y, z, w), [(z, w, 0)], [y])
    y2, z2, w2 = chain["twins"]
    certify("delayed_state_chain", chain, [y - y2, z - z2, w - w2])
    check("delayed_states_same_immediate_output", y.subs(y, 0), 0)
    check("delayed_states_same_first_rate", lie(y, (z, w, 0), (y, z, w)).subs(z, 0), 0)
    check("delayed_states_second_rate_distinguishes", lie(z, (z, w, 0), (y, z, w)), w)
    removable = predictive_pair_ideal((y, z, w), [(z, -z, 0)], [y])
    y2, z2, w2 = removable["twins"]
    certify("removable_hidden_coordinate", removable, [y - y2, z - z2])
    constant = predictive_pair_ideal((y, z), [(1, z)], [s.Integer(1)])
    certify("constant_output_needs_no_state", constant, [])

    # Singular coordinates: fibre-constant coefficients can still lose ODE
    # uniqueness. Both an immediate trajectory and a waiting trajectory solve it.
    check("cube_projection_immediate_solution", s.diff(elapsed ** 3, elapsed), 3 * elapsed ** 2)
    wait = s.Symbol("wait", nonnegative=True)
    check("cube_projection_waiting_branch_equation",
          s.diff((elapsed - wait) ** 3, elapsed), 3 * (elapsed - wait) ** 2)
    check("cube_projection_waiting_branch_C1_join",
          s.diff((elapsed - wait) ** 3, elapsed).subs(elapsed, wait), 0)
    flat = s.exp(-1 / z ** 2)
    for order in range(6):
        check(f"flat_memory_right_jet_{order}",
              s.limit(s.diff(flat, z, order), z, 0, dir="+"))
    # Positivity on a subinterval is an exact lower bound for the distinguishing
    # integral; the proof covers every t in (0,1), not just these finite jets.
    check("flat_memory_positive_subinterval_bound",
          flat.subs(z, t / 2), s.exp(-4 / t ** 2))

    # Native selected-null quotient extended by common dynamical invariance.
    A = s.diag(-1, -1, -2)
    C = s.Matrix([[1, 2, 0]])
    obs = linear_observation_rows([A], C)
    check("linear_removable_material_count", obs.rows, 1)
    hidden = s.Matrix([-2, 1, 0])
    check("linear_hidden_mode_is_output_null", C * hidden)
    check("linear_hidden_mode_is_invariant", A * hidden, -hidden)
    A0 = s.Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
    A1 = s.Matrix([[0, 0, 0], [0, 0, 1], [0, 0, 0]])
    switched = linear_observation_rows([A0, A1], s.Matrix([[1, 0, 0]]))
    check("switched_material_needs_all_three_modes", switched, s.eye(3))
    rates = [s.Integer(1), s.Integer(2), s.Integer(4)]
    weights = [s.Integer(2), s.Integer(3), s.Integer(5)]
    decay = -s.diag(*rates)
    material_output = s.Matrix([weights])
    finite_rows = linear_observation_rows([decay], material_output)
    check("three_relaxations_exact_material_count", finite_rows, s.eye(3))
    obs_matrix = s.Matrix.vstack(*(material_output * decay ** j for j in range(3)))
    vandermonde = s.prod(weights) * s.prod(-rates[j] + rates[i]
                                           for i in range(3) for j in range(i + 1, 3))
    check("relaxation_observability_vandermonde", obs_matrix.det(), vandermonde)
    B = s.ones(3, 1)
    reachability = s.Matrix.hstack(*(decay ** j * B for j in range(3)))
    require("all_three_relaxations_reachable", reachability.det() != 0)
    zvec = s.Matrix([y, z, w])
    drive = s.Symbol("drive", real=True)
    H = sum(c * zz ** 2 / 2 for c, zz in zip(weights, zvec))
    check("relaxation_passivity_balance", lie(H, decay * zvec + B * drive, (y, z, w)),
          drive * (material_output * zvec)[0]
          - sum(c * rate * zz ** 2 for c, rate, zz in zip(weights, rates, zvec)))
    decay_merged = -s.diag(1, 1, 2)
    merged = linear_observation_rows([decay_merged], s.Matrix([[1, 2, 3]]))
    check("equal_decay_rates_merge_to_two_modes", merged.rows, 2)
    for size in range(1, 6):
        moment = s.Matrix(size, size, lambda i, j:
                          s.Rational(2 ** (i + j + 1) - 1, i + j + 1))
        require(f"continuum_memory_moment_minor_{size}", moment.det() > 0)
    return {
        "scope": "Minimum predictive state for declared controls/outputs; polynomial pair closure and native load",
        "proof": "docs/2026-09-05-cella-predictive-state.md",
        "check_count": len(checks),
        "checks": checks,
        "polynomial_certificates": certificates,
        "native_state_counts": {
            "known_coefficients_fixed_a_original_two_ports_output_F": 2,
            "known_coefficients_radial_only_beta_nonzero_output_F": 1,
            "unknown_constant_mu_beta_nu_original_two_ports_output_d": 5,
            "unknown_constant_k_correlated_family_original_two_ports_output_d": 3,
        },
        "boundaries": [
            "Counts refer to the stated law, inputs, outputs and admissible real chart.",
            "Finite checked flat jets and moment minors accompany all-order manuscript proofs.",
            "General smooth laws require a closure certificate or actual future-flow comparison.",
            "No physical coframe, cut or constitutive coefficients are selected by these checks.",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2))

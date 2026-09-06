#!/usr/bin/env python3
"""Exact constructions for response on branched, resolved state spaces.

Proof: docs/results/2026-09-05/2026-09-05-response-on-branched-resolved-state-spaces.md.
Prints a JSON receipt; does not write files. Symbolic matrix ranks are generic
over the coefficient field: specialize the full relation before a rank jump.
Smooth divisibility is decided here only for explicitly polynomial quotients.
"""

import json

import sympy as s


def simplify_matrix(matrix):
    return s.Matrix(matrix).applyfunc(lambda v: s.factor(s.cancel(v)))


def columns(vectors, rows):
    return s.Matrix.hstack(*vectors) if vectors else s.zeros(rows, 0)


def image_basis(matrix):
    matrix = simplify_matrix(matrix)
    return columns(matrix.columnspace(), matrix.rows)


def kernel_basis(matrix):
    matrix = s.Matrix(matrix)
    return columns(matrix.nullspace(), matrix.cols)


def graph_relation(matrix):
    matrix = s.Matrix(matrix)
    if matrix.rows != matrix.cols:
        raise ValueError("a graph matrix must be square")
    return s.eye(matrix.rows).col_join(matrix)


def passive_relation(effort_basis, response):
    """BR-1 in a basis of the effort domain; no positivity is presumed.

    The result is passive exactly when sym(response) is positive semidefinite.
    D(D^T D)^-1 is only a coordinate section for representing the restriction
    of flow to the effort domain, not a physical metric.
    """
    domain, response = s.Matrix(effort_basis), s.Matrix(response)
    n, k = domain.shape
    if domain.rank() != k or response.shape != (k, k):
        raise ValueError("independent effort columns and a matching K are required")
    vertical = kernel_basis(domain.T)
    section = domain * (domain.T * domain).inv() if k else s.zeros(n, 0)
    effort = domain.row_join(s.zeros(n, vertical.cols))
    flow = (section * response).row_join(vertical)
    return simplify_matrix(effort.col_join(flow))


def split_relation(relation):
    relation = s.Matrix(relation)
    if relation.rows % 2:
        raise ValueError("a relation needs equally sized effort and flow blocks")
    n = relation.rows // 2
    return relation[:n, :], relation[n:, :]


def push_relation(transform, relation):
    """Exact BR-2 image by solving A^T u = E t and retaining j'=A J t."""
    transform = s.Matrix(transform)
    effort, flow = split_relation(relation)
    m, n = transform.shape
    if effort.rows != n:
        raise ValueError("port transformation and relation dimensions disagree")
    k = effort.cols
    admissible = kernel_basis(transform.T.row_join(-effort))
    output = s.eye(m).row_join(s.zeros(m, k)).col_join(
        s.zeros(m, m).row_join(transform * flow))
    return image_basis(output * admissible)


def dual_relation(relation):
    effort, flow = split_relation(relation)
    return flow.col_join(effort)


def scattering_relation(contraction):
    """Regular BR-1 coordinate; passivity is equivalent to I-C^T C >= 0."""
    contraction = s.Matrix(contraction)
    if contraction.rows != contraction.cols:
        raise ValueError("the scattering coordinate must be square")
    identity = s.eye(contraction.rows)
    return (identity + contraction).col_join(identity - contraction) / 2


def relation_scattering(relation):
    """Inverse scattering coordinate for a maximal passive relation basis."""
    effort, flow = split_relation(image_basis(relation))
    if effort.rows != effort.cols or (effort + flow).det() == 0:
        raise ValueError("relation has no full scattering graph in these coordinates")
    return simplify_matrix((effort - flow) * (effort + flow).inv())


def clamp_internal(relation, external_count):
    """Set internal flows to zero, then project to external effort/flow."""
    effort, flow = split_relation(relation)
    n = effort.rows
    if not 0 <= external_count <= n:
        raise ValueError("external port count is out of range")
    allowed = kernel_basis(flow[external_count:, :])
    external = effort[:external_count, :].col_join(flow[:external_count, :])
    return image_basis(external * allowed)


def direct_sum_relations(*relations):
    blocks = [split_relation(relation) for relation in relations]
    return s.diag(*(e for e, _ in blocks)).col_join(
        s.diag(*(j for _, j in blocks)))


def same_relation(left, right):
    left, right = s.Matrix(left), s.Matrix(right)
    return (left.rows == right.rows and left.rank() == right.rank()
            == left.row_join(right).rank())


def lift_graph(jacobian, response, drift=None):
    """Unique rational candidate; smooth extension still requires BR-3.

    A zero determinant is rejected. A nonzero symbolic determinant can vanish
    on a stratum; cancellation here does not establish arbitrary C-infinity
    divisibility. The returned numerators retain its exact certificate target.
    """
    jacobian, response = s.Matrix(jacobian), s.Matrix(response)
    n = jacobian.rows
    if jacobian.cols != n or response.shape != (n, n):
        raise ValueError("lifting needs square matching matrices")
    determinant = s.factor(jacobian.det())
    if determinant == 0:
        raise ValueError("generically singular maps need relation-level lifting")
    adjugate = jacobian.adjugate()
    numerator = simplify_matrix(adjugate * response * adjugate.T)
    result = {"determinant": determinant, "numerator": numerator,
              "response": simplify_matrix(numerator / determinant**2)}
    if drift is not None:
        drift = s.Matrix(drift)
        if drift.shape != (n, 1):
            raise ValueError("drift dimension disagrees")
        result["drift_numerator"] = simplify_matrix(adjugate * drift)
        result["drift"] = simplify_matrix(adjugate * drift / determinant)
    return result


def polynomial_quotient(matrix, coordinates):
    """Certificate for polynomial fixtures, not a general smoothness oracle."""
    return all(s.cancel(v).is_polynomial(*coordinates) for v in matrix)


def descend_even_polynomial(expression, roots, boundaries):
    """Descend a polynomial in independent roots t_b through r_b=t_b^2."""
    if len(roots) != len(boundaries):
        raise ValueError("root and boundary coordinate counts disagree")
    result = s.Integer(0)
    for powers, coefficient in s.Poly(expression, *roots).terms():
        if any(power % 2 for power in powers):
            raise ValueError("an odd root term does not descend smoothly as a polynomial")
        result += coefficient * s.prod(r ** (power // 2)
                                       for r, power in zip(boundaries, powers))
    return s.expand(result)


def run_checks():
    checks = {}
    certificates = {}

    def check(name, actual, expected=0):
        if name in checks:
            raise AssertionError("duplicate check name: " + name)
        if isinstance(actual, s.MatrixBase):
            target = s.zeros(*actual.shape) if expected == 0 else s.Matrix(expected)
            residuals = list(actual - target)
        else:
            residuals = [s.sympify(actual) - s.sympify(expected)]
        for residual in residuals:
            residual = s.simplify(s.cancel(s.expand(residual)))
            if residual != 0:
                raise AssertionError((name, residual))
        checks[name] = {"verified": True}

    def require(name, condition):
        if name in checks or not condition:
            raise AssertionError(name)
        checks[name] = {"verified": True}

    def equal_relation(name, left, right):
        require(name, same_relation(left, right))

    # Linear relations, including constrained effort and vertical flow.
    M = s.Matrix([[3, 2, -1], [0, 4, 1], [1, -1, 2]])
    R = (M + M.T) / 2
    for size in range(1, 4):
        require(f"graph_positive_leading_minor_{size}", R[:size, :size].det() > 0)
    domain = s.Matrix([[1, 0], [0, 1], [1, 0]])
    K = s.Matrix([[2, 1], [-1, 3]])
    proper = passive_relation(domain, K)
    E, J = split_relation(proper)
    check("proper_domain_restriction", domain.T * J[:, :2], K)
    check("proper_vertical_annihilator", domain.T * J[:, 2:])
    check("proper_relation_dimension", proper.rank(), 3)
    check("proper_power_matrix", (E.T * J + J.T * E) / 2, s.diag(2, 3, 0))
    opened = passive_relation(s.eye(3), s.zeros(3))
    shorted = passive_relation(s.zeros(3, 0), s.zeros(0))
    equal_relation("open_is_zero_flow", opened, graph_relation(s.zeros(3)))
    equal_relation("short_is_zero_effort", shorted, s.zeros(3).col_join(s.eye(3)))
    A = s.Matrix([[1, 0, 1], [0, 1, 1]])
    B = s.Matrix([[2, -1]])
    maps = {"rectangular_full": A, "rectangular_rank_one": s.Matrix([[1, 0, 0], [2, 0, 0]]),
            "square_singular": s.diag(1, 0, 1), "zero": s.zeros(2, 3)}
    relations = {"graph": graph_relation(M), "proper": proper, "open": opened, "short": shorted}
    for rname, relation in relations.items():
        scattering = relation_scattering(relation)
        equal_relation(f"scattering_round_trip_{rname}", scattering_relation(scattering), relation)
        basis = scattering_relation(scattering)
        es, js = split_relation(basis)
        check(f"scattering_power_identity_{rname}",
              (es.T*js+js.T*es)/2, (s.eye(3)-scattering.T*scattering)/4)
        for aname, transform in maps.items():
            pushed = push_relation(transform, relation)
            check(f"maximal_dimension_{rname}_{aname}", pushed.rank(), transform.rows)
        equal_relation(f"composition_{rname}", push_relation(B, push_relation(A, relation)),
                       push_relation(B * A, relation))
        equal_relation(f"identity_{rname}", push_relation(s.eye(3), relation), relation)
    equal_relation("graph_push_is_congruence", push_relation(A, graph_relation(M)),
                   graph_relation(A * M * A.T))
    T = s.Matrix([[1, 2, 0], [0, 1, 1], [0, 0, 1]])
    equal_relation("port_inverse_recovers_proper_relation",
                   push_relation(T.inv(), push_relation(T, proper)), proper)
    u = s.Matrix(s.symbols("u0:2", real=True))
    j = s.Matrix(s.symbols("j0:3", real=True))
    check("rectangular_transform_preserves_power", (u.T * A * j)[0], ((A.T * u).T * j)[0])
    network = direct_sum_relations(proper, graph_relation(s.Matrix([[2]])))
    check("direct_sum_maximal_dimension", network.rank(), 4)
    en, jn = split_relation(network)
    check("direct_sum_power_matrix", (en.T * jn + jn.T * en) / 2, s.diag(2, 3, 0, 2))
    check("open_scattering_is_identity", relation_scattering(opened), s.eye(3))
    check("short_scattering_is_negative_identity", relation_scattering(shorted), -s.eye(3))
    require("proper_effort_free_flow_changes_state_for_identity_anchor",
            s.eye(3)*J[:, 2:] != s.zeros(3, 1))
    check("proper_effort_anchor_removes_all_free_flow", domain.T * J[:, 2:])

    # Arbitrary mixed-minor and weighted smooth-lift certificates.
    t, z, d, xi = s.symbols("t z d xi", real=True)
    contraction = s.Matrix([[(t**2-1)/(t**2+1)]])
    regular_relation = scattering_relation(contraction)
    check("smooth_constraint_family_contraction_defect", (s.eye(1)-contraction.T*contraction)[0],
          4*t**2/(1+t**2)**2)
    equal_relation("smooth_constraint_family_basis", regular_relation, s.Matrix([t**2, 1]))
    equal_relation("smooth_constraint_family_generic_graph", regular_relation, graph_relation(s.Matrix([[1/t**2]])))
    equal_relation("smooth_constraint_family_special_short", regular_relation.subs(t, 0), s.Matrix([0, 1]))
    require("smooth_constraint_family_graph_does_not_extend", not (1/t**2).is_polynomial(t))
    observation = s.Matrix([t*z, t**2/2])
    A = observation.jacobian([t, z])
    upper = s.Matrix([[2, 1], [1, 3]])
    upper_drift = s.Matrix([1, z])
    lower = A * upper * A.T
    lifted = lift_graph(A, lower, A * upper_drift)
    check("mixed_minor_lift_determinant", lifted["determinant"], -t**2)
    check("mixed_minor_lift_numerator", lifted["numerator"], t**4 * upper)
    check("mixed_minor_lift_reconstructs", lifted["response"], upper)
    check("mixed_minor_drift_numerator", lifted["drift_numerator"], -t**2 * upper_drift)
    check("mixed_minor_drift_reconstructs", lifted["drift"], upper_drift)
    require("mixed_minor_polynomial_certificate", polynomial_quotient(lifted["response"], (t, z)))
    pole = lift_graph(A, s.eye(2))
    require("constant_lower_response_not_polynomial", not polynomial_quotient(pole["response"], (t, z)))
    check("constant_lower_response_has_order_two_pole", pole["response"][0, 0], 1/t**2)
    L = s.Matrix([[1, z, 0], [0, 1, t], [0, 0, 1]])
    Q = s.diag(t, t**2, 1)
    W = s.Matrix([[1, 0, z], [0, 1, 0], [0, t, 1]])
    weighted = L * Q * W
    weighted_lower = weighted * M * weighted.T
    check("weighted_general_factor_lift", lift_graph(weighted, weighted_lower)["response"], M)
    adapted = simplify_matrix(L.inv() * weighted_lower * L.inv().T)
    check("weighted_adapted_divisors", simplify_matrix(Q.inv() * adapted * Q.inv()), W * M * W.T)
    A2 = s.Matrix([[1, t], [0, t]])
    nested_lower = A * A2 * upper * A2.T * A.T
    check("iterated_lift_matches_composite",
          lift_graph(A2, lift_graph(A, nested_lower)["response"])["response"],
          lift_graph(A * A2, nested_lower)["response"])
    for k, order in ((2, 1), (3, 2), (4, 2), (3, 1), (4, 1)):
        root_lift = lift_graph(s.Matrix([[k * t**(k-1)]]), s.Matrix([[t**(k*order)]]))
        require(f"root_order_criterion_k{k}_m{order}",
                polynomial_quotient(root_lift["response"], (t,)) == (k * order >= 2*k - 2))
    check("one_sided_linear_law_square_root_lift",
          lift_graph(s.Matrix([[2*t]]), s.Matrix([[t**2]]))["response"], s.Matrix([[s.Rational(1, 4)]]))

    # Complete mixed-boundary class: PSD factor and parity, not entry orders alone.
    f, r1, r2, t1, t2 = s.symbols("f r1 r2 t1 t2", real=True)
    q = s.diag(f*t1, f*t1, t2, 1)
    L = s.Matrix([[1, 0, 0, 0], [1, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    v = s.Matrix([t1, t1, t2, 1])
    barR = L * L.T + v * v.T
    rootR = q * barR * q
    for number, (root, parity) in enumerate(((t1, s.diag(-1, -1, 1, 1)),
                                            (t2, s.diag(1, 1, -1, 1))), 1):
        check(f"root_coefficient_parity_{number}", barR.subs(root, -root), parity * barR * parity)
        check(f"root_response_even_{number}", rootR.subs(root, -root), rootR)
    realR = rootR.applyfunc(lambda v: descend_even_polynomial(v, (t1, t2), (r1, r2)))
    check("root_descent_round_trip", realR.subs({r1: t1**2, r2: t2**2}), rootR)
    for i in (0, 1):
        check(f"two_sided_face_row_{i}", realR[i, :].subs(f, 0))
        check(f"one_sided_shared_face_row_{i}", realR[i, :].subs(r1, 0))
    check("second_boundary_face_row", realR[2, :].subs(r2, 0))
    check("shared_face_diagonal_order", realR[0, 0], f**2 * r1 * (1+r1))
    check("shared_face_off_diagonal_order", realR[0, 1], f**2 * r1 * (1+r1))
    check("single_row_boundary_off_diagonal_order", realR[0, 3], f * r1)
    effort = s.Matrix(s.symbols("e0:4", real=True))
    squares = ((L.T * q * effort).T * (L.T * q * effort))[0] + (v.T * q * effort)[0]**2
    check("mixed_boundary_global_sum_of_squares", (effort.T * rootR * effort)[0], squares)
    faces = [{f, r1}, {f, r1}, {r2}, set()]
    skew = s.zeros(4)
    for i in range(4):
        for j in range(i+1, 4):
            skew[i, j] = (i+j+1) * s.prod(faces[i] | faces[j])
            skew[j, i] = -skew[i, j]
    check("mixed_skew_is_skew", skew + skew.T)
    check("mixed_skew_has_zero_power", (effort.T * skew * effort)[0])
    check("common_face_skew_has_single_power", skew[0, 1], 2*f*r1)
    for i, row_faces in enumerate(faces):
        for face in sorted(row_faces, key=str):
            check(f"skew_constraint_row_{i}_{face}", skew[i, :].subs(face, 0))
    bad_root = s.diag(t1, 1) * s.Matrix([[2, 1], [1, 2]]) * s.diag(t1, 1)
    try:
        descend_even_polynomial(bad_root[0, 1], (t1,), (r1,))
    except ValueError:
        checks["positive_root_coefficient_without_parity_rejects_descent"] = {"verified": True}
    else:
        raise AssertionError("odd-root response falsely descended")
    require("one_sided_first_order_is_not_two_sided_passive", r1.subs(r1, -1) < 0)
    certificates["mixed_boundary"] = {
        "dissipation": str(s.factor(squares)),
        "shared_dissipative_entry": str(realR[0, 1]),
        "shared_skew_entry": str(skew[0, 1]),
        "root_parities": ["diag(-1,-1,1,1)", "diag(1,1,-1,1)"],
    }

    # Nonlinear integration and a negative integrability control.
    e0, e1, lam = s.symbols("e0 e1 lam", real=True)
    efforts = s.Matrix([e0, e1])
    phi = s.Matrix([e0, r1 * (e1 + e1**3)])
    jac = phi.jacobian(efforts)
    recovered = (jac.subs({e0: lam*e0, e1: lam*e1}, simultaneous=True) * efforts).applyfunc(
        lambda v: s.integrate(v, (lam, 0, 1)))
    check("nonlinear_radial_integration", recovered, phi)
    check("nonlinear_boundary_output", phi[1].subs(r1, 0))
    check("nonlinear_positive_jacobian", jac, s.diag(1, r1*(1+3*e1**2)))
    check("nonlinear_row_curl", s.Matrix([s.diff(jac[i, 0], e1)-s.diff(jac[i, 1], e0)
                                          for i in range(2)]))
    nonintegrable = s.Matrix([[1, e0], [-e0, 1]])
    check("nonintegrable_positive_symmetric_part", (nonintegrable+nonintegrable.T)/2, s.eye(2))
    check("nonintegrable_nonzero_row_curl", s.diff(nonintegrable[0, 0], e1)-s.diff(nonintegrable[0, 1], e0), -1)
    nonlinear_A = s.Matrix([[1, t], [0, 1]])
    pulled = phi.subs(dict(zip(efforts, nonlinear_A.T * efforts)), simultaneous=True)
    check("nonlinear_derivative_congruence", (nonlinear_A * pulled).jacobian(efforts),
          nonlinear_A * jac.subs(dict(zip(efforts, nonlinear_A.T * efforts)), simultaneous=True) * nonlinear_A.T)

    # Singular passive elimination, complete relations, and specialization.
    full = s.Matrix([[4, 1, 2], [-1, 3, 1], [0, 1, 2]])
    fullR = (full + full.T)/2
    for size in range(1, 4):
        require(f"network_positive_leading_minor_{size}", fullR[:size, :size].det() > 0)
    D = full[1:, 1:]
    schur = full[:1, :1] - full[:1, 1:] * D.inv() * full[1:, :1]
    eliminated = clamp_internal(graph_relation(full), 1)
    equal_relation("invertible_elimination_is_schur", eliminated, graph_relation(schur))
    projection = s.Matrix([[1, 0, 0]])
    equal_relation("zero_flow_elimination_dual_transform_construction", eliminated,
                   dual_relation(push_relation(projection, dual_relation(graph_relation(full)))))
    require("schur_passive", schur[0, 0] > 0)
    safe = s.Matrix([[3, 1, 0], [1, 2, 0], [0, 0, 0]])
    ga, gb, gc = s.symbols("ga gb gc", real=True)
    generalized = s.Matrix([[s.Rational(1, 2), ga], [gb, gc]])
    D = safe[1:, 1:]
    check("singular_generalized_inverse_identity", D * generalized * D, D)
    safe_schur = safe[:1, :1] - safe[:1, 1:] * generalized * safe[1:, :1]
    check("singular_effective_independent_of_inverse", safe_schur, s.Matrix([[s.Rational(5, 2)]]))
    equal_relation("singular_safe_elimination_graph", clamp_internal(graph_relation(safe), 1), graph_relation(safe_schur))
    check("singular_internal_kernel_invisible", safe[:1, 1:] * kernel_basis(D))
    constrained = s.Matrix([[2, 1], [-1, 0]])
    constrained_external = clamp_internal(graph_relation(constrained), 1)
    equal_relation("singular_constraint_retains_short_relation", constrained_external, s.Matrix([[0], [1]]))
    check("singular_constraint_passivity", (constrained+constrained.T)/2, s.diag(2, 0))
    require("singular_constraint_range_condition_fails", constrained[1:, 1:].row_join(constrained[1:, :1]).rank()
            > constrained[1:, 1:].rank())
    family = s.Matrix([[1, t], [t, t**2]])
    check("specialization_family_psd_factor", family, s.Matrix([1, t]) * s.Matrix([[1, t]]))
    generic = clamp_internal(graph_relation(family), 1)
    special = clamp_internal(graph_relation(family.subs(t, 0)), 1)
    equal_relation("generic_elimination_is_open", generic, graph_relation(s.zeros(1)))
    equal_relation("special_fibre_elimination_is_unit", special, graph_relation(s.eye(1)))
    require("elimination_and_specialization_differ", not same_relation(generic.subs(t, 0), special))
    internal_effort = -e0/t
    check("divergent_internal_solution_zero_flow", family[1, 0]*e0 + family[1, 1]*internal_effort)
    require("internal_solution_has_a_pole", not internal_effort.is_polynomial(t))
    compatible = t**2 * s.Matrix([[2, 1], [1, 1]])
    check("smooth_internal_solver_module_identity", compatible[1:, 1:] * s.ones(1), compatible[1:, :1])
    effective = compatible[:1, :1] - compatible[:1, 1:] * s.ones(1)
    equal_relation("smooth_solver_effective_response", clamp_internal(graph_relation(compatible), 1), graph_relation(effective))
    equal_relation("smooth_solver_specialization_commutes", clamp_internal(graph_relation(compatible.subs(t, 0)), 1),
                   graph_relation(effective.subs(t, 0)))
    certificates["singular_network"] = {
        "full_family": str(family), "generic_effective": "0", "special_effective": "1",
        "internal_effort": str(internal_effort), "regular_solver_effective": str(effective[0, 0]),
    }

    # Storage and the distinction between a zero output and a normal barrier.
    state = s.Matrix([f, z])
    anchor = s.Matrix([[1, z], [f, 1]])
    drift = s.Matrix([z, -f])
    energy_gradient = state
    external_effort = s.Matrix([e0, e1])
    actual_effort = external_effort - anchor.T * energy_gradient
    constitutive = s.Matrix([[2, 1], [-1, 3]])
    flow = constitutive * actual_effort
    energy_rate = (energy_gradient.T * (drift + anchor * flow))[0]
    check("storage_exact_power_balance", energy_rate,
          (energy_gradient.T * drift)[0] + (external_effort.T * flow)[0] - (actual_effort.T * flow)[0])
    check("storage_dissipation_squares", (actual_effort.T * flow)[0], 2*actual_effort[0]**2+3*actual_effort[1]**2)
    crossing = s.diag(1, f**2)
    path = s.Matrix([t, t**3/3])
    check("zero_output_regular_crossing_solution", path.diff(t), (crossing*s.ones(2, 1)).subs(f, t))
    check("crossing_face_second_output_zero", crossing[1, :].subs(f, 0))
    check("crossing_normal_speed_nonzero", path[0].diff(t), 1)
    barrier = s.diag(f**2, 1)
    check("normal_output_constraint_is_tangent", (barrier*s.ones(2, 1))[0].subs(f, 0))
    initial = s.symbols("initial", real=True)
    barrier_solution = initial / (1-initial*t)
    check("barrier_flow_exact", s.diff(barrier_solution, t), barrier_solution**2)
    check("barrier_zero_initial_stays_zero", barrier_solution.subs(initial, 0))
    root_flow = 2*t*e0
    check("invariant_base_drive_root_tangent", root_flow.subs(t, 0))
    check("invariant_base_drive_root_equivariant", root_flow.subs(t, -t), -root_flow)
    require("constant_upper_drive_not_fixed_input_equivariant", s.Integer(1) != -s.Integer(1))

    # Native chart and HG-1 ruler swap, including their resolved extensions.
    a, b, gamma, mu, beta, nu = s.symbols("a b gamma mu beta nu", real=True)
    p = 1+a**2
    F = (1-xi**2)/p
    epsilon = 1-gamma**2
    delta = epsilon + a**2+b**2-2*a*b*gamma
    native_A = s.Matrix([[d, xi], [0, 1]])
    native_M = s.Matrix([[mu, beta], [beta, nu]])
    native_lower = native_A * native_M * native_A.T
    check("native_lift_recovered", lift_graph(native_A, native_lower)["response"], native_M)
    check("native_exceptional_direction_response", native_lower.subs(d, 0), nu*s.Matrix([[xi**2, xi], [xi, 1]]))
    check("native_direction_recovered_from_response", s.cancel(native_lower[0, 1].subs(d, 0)/native_lower[1, 1].subs(d, 0)), xi)
    xi_h, xi_g = (b-a*gamma)/d, (a-b*gamma)/d
    check("two_depth_norms_same_lapse",
          ((1-xi_h**2)/p - (1-xi_g**2)/(1+b**2)).subs(d**2, delta))
    swap_xi = -xi_g
    check("swap_minus_deck_exact", swap_xi + xi_h, (b-a)*(1+gamma)/d)
    check("swap_deck_equal_direction_equal_depths", (swap_xi+xi_h).subs(b, a))
    check("negative_horizon_swap_deck_same_direction", (swap_xi+xi_h).subs(gamma, -1))
    check("positive_horizon_swap_fixes_direction", (swap_xi-xi_h).subs(gamma, 1))
    # Polynomial certificates on gamma^2+F*d^2=1, with p*F=1-xi^2.
    Fsymbol = s.symbols("F", real=True)
    transformed_a = a*gamma+xi*d
    transformed_xi = gamma*xi-a*Fsymbol*d
    constraint = gamma**2+Fsymbol*d**2-1
    lapse_residual = s.expand(1-transformed_xi**2-Fsymbol*(1+transformed_a**2))
    lapse_residual = s.rem(lapse_residual, constraint, gamma).subs(Fsymbol, F)
    check("resolved_swap_lapse_invariance", lapse_residual)
    check("resolved_swap_from_unresolved_formula",
          s.expand(swap_xi.subs(b, transformed_a)-transformed_xi), a*constraint/d)
    second_a = transformed_a*gamma-transformed_xi*d
    second_xi = gamma*transformed_xi+transformed_a*Fsymbol*d
    check("resolved_swap_involution_a", s.rem(s.expand(second_a-a), constraint, gamma))
    check("resolved_swap_involution_xi", s.rem(s.expand(second_xi-xi), constraint, gamma))
    state3 = s.Matrix([a, xi, d])
    tau_state = s.Matrix([a, -xi, -d])
    S_state = s.Matrix([transformed_a, transformed_xi, -d])
    check("swap_commutes_with_deck", S_state.subs({xi: -xi, d: -d}, simultaneous=True),
          s.Matrix([transformed_a, -transformed_xi, d]))
    for sigma in (1, -1):
        check(f"resolved_swap_exceptional_sigma_{sigma}", S_state.subs({d: 0, gamma: sigma}),
              s.Matrix([sigma*a, sigma*xi, 0]))
        check(f"resolved_swap_horizon_sigma_{sigma}", S_state.subs({Fsymbol: 0, gamma: sigma}),
              s.Matrix([sigma*a+xi*d, sigma*xi, -d]))
    plus_derivative = s.Matrix([[1, 0, xi], [0, 1, -a*F], [0, 0, -1]])
    positive_gamma = s.sqrt(1-F*d**2)
    positive_swap = S_state.subs({gamma: positive_gamma, Fsymbol: F}, simultaneous=True)
    check("positive_pinch_swap_derivative_from_full_map", positive_swap.jacobian(state3).subs(d, 0), plus_derivative)
    check("positive_pinch_swap_derivative_involution", plus_derivative**2, s.eye(3))
    averaged = simplify_matrix((s.eye(3)+plus_derivative*plus_derivative.T)/2)
    check("finite_group_average_swap_covariant_at_fixed_stratum", plus_derivative*averaged*plus_derivative.T, averaged)
    ea, ex = s.symbols("ea ex", real=True)
    invariant_effort = s.Matrix([ea, ex, (xi*ea-a*F*ex)/2])
    check("positive_pinch_invariant_effort", plus_derivative.T*invariant_effort, invariant_effort)
    check("invariant_response_and_drive_have_zero_normal_speed", (averaged*invariant_effort)[2])
    check("orientation_sensitive_drive_has_normal_speed", (averaged*s.Matrix([0, 0, 1]))[2], 1)
    # The strict transform, not every direction above a collapsed equal-depth point.
    equal_depth_equation = a*(gamma-1)+xi*d
    strict_plus = xi-a*Fsymbol*d/(1+gamma)
    check("equal_depth_strict_transform_regular_identity",
          equal_depth_equation-d*strict_plus, a*constraint/(1+gamma))
    check("equal_depth_strict_transform_positive_pinch", strict_plus.subs({d: 0, gamma: 1}), xi)
    # Same declared fixed-a diagonal law, same two-phase experiment as predictive state.
    elapsed = s.symbols("elapsed", real=True)
    original_x_field = s.Matrix([mu*d, nu*xi])
    F_gradient = s.Matrix([s.diff(F, xi), s.diff(F, d)])
    F_rate = (F_gradient.T*original_x_field)[0]
    for sign in (1, -1):
        check(f"fixed_input_deck_partner_lapse_rate_{sign}", F_rate.subs({xi: sign, d: nu*elapsed}),
              -2*mu*nu*elapsed*sign/p)
    certificates["ruler_swap"] = {
        "resolved_action": [str(transformed_a), str(transformed_xi), "-d"],
        "constraints": ["gamma**2 + F*d**2 = 1", "(1+a**2)*F = 1-xi**2"],
        "positive_exceptional_action": "(a,xi,0) -> (a,xi,0)",
        "negative_exceptional_action": "(a,xi,0) -> (-a,-xi,0)",
        "positive_fixed_stratum_invariant_effort": str(invariant_effort),
        "equal_depth_positive_strict_transform": "xi=0 at d=0",
    }
    return {
        "scope": "Finite-dimensional response relations, smooth resolved graphs, boundary divisors and controlled evolution",
        "proof": "docs/results/2026-09-05/2026-09-05-response-on-branched-resolved-state-spaces.md",
        "check_count": len(checks), "checks": checks, "certificates": certificates,
        "boundaries": [
            "General theorems have manuscript proofs; these exact checks certify implementations and examples.",
            "Symbolic relation ranks are generic; specialize full data before singular elimination.",
            "Polynomial cancellation certifies the displayed fixtures, not arbitrary smooth ideal membership.",
            "Ruler swap moves a as well as xi,d; fixed-a dynamics are not automatically swap covariant.",
            "No physical constitutive coefficients, cut, coframe or input selection are inferred.",
        ],
    }


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2))

#!/usr/bin/env python3
"""Exact Cella/native presentation transport and soldering classification.

The runner checks the constructive identities in
docs/results/2026-09-05/2026-09-05-cella-frame-transport.md.  Standard connection, torsion,
curvature, and local integrability mathematics is retained explicitly.
No state response, physical connection, or spacetime metric is selected by
these checks.  The curved example deliberately has no trapping coincidence.
"""

import json

import sympy as s


def reduced(value):
    # Exact rational cancellation is sufficient for the generic jet identities.
    # Full multivariate factorization can spend minutes on an already-cancelled
    # nonzero gauge tensor. Keep common factors without factoring its polynomial.
    value = s.factor_terms(s.cancel(value))
    return s.simplify(value) if value.has(s.sin, s.cos, s.sinh, s.cosh, s.exp) else value


def native_gram(a, xi):
    return s.Matrix([[-1, a, xi], [a, 1, 0],
                     [xi, 0, (1-xi**2)/(1+a**2)]])


def native_lift(a, xi):
    root = s.sqrt(1+a**2)
    return s.Matrix([[1, -a, -xi], [0, root, a*xi/root],
                     [0, 0, 1/root]])


def invariant_hypersurface_jet(gradient, hessian):
    """Return the complete Euclidean defining-function second-jet quotient.

    The invariant carrier is (projector, tensor).  ``tensor[i]`` is the
    matrix T_ijk in its last two slots; T is minus the normal-valued second
    fundamental form, extended to vanish on normal input slots.  ``q``,
    ``tangent_hessian`` and ``gauge_covector`` are reconstruction data and
    are not separately invariant under arbitrary nonzero rescaling of f.

    Domain: a nonzero real gradient and a symmetric Hessian at f=0 in a
    fixed orthonormal Euclidean presentation frame.  A symbolic nonzero
    gradient is accepted with that condition retained by the caller.
    """
    gradient = s.Matrix(gradient)
    hessian = s.Matrix(hessian)
    n = gradient.rows
    if gradient.cols != 1 or hessian.shape != (n, n):
        raise ValueError("Expected an n-vector and an n by n Hessian.")
    if any(reduced(value) != 0 for value in hessian-hessian.T):
        raise ValueError("The Hessian must be symmetric.")
    q = reduced((gradient.T*gradient)[0])
    if q == 0:
        raise ValueError("The defining-function gradient must be nonzero.")
    projector = (s.eye(n)-gradient*gradient.T/q).applyfunc(reduced)
    tangent_hessian = (projector*hessian*projector).applyfunc(reduced)
    gauge_covector = (hessian*gradient/q
        -(gradient.T*hessian*gradient)[0]*gradient/(2*q*q)).applyfunc(reduced)
    tensor = tuple((gradient[i]*tangent_hessian/q).applyfunc(reduced)
                   for i in range(n))
    return {"q": q, "projector": projector, "tensor": tensor,
            "tangent_hessian": tangent_hessian,
            "gauge_covector": gauge_covector}


def canonical_edge_vector(tensor):
    """Extract (T_123,T_213,T_312) in the declared three presentation axes.

    This readout is defining-function-gauge invariant. It does not carry
    the full O(3) tensor action and is not a complete geometric invariant.
    """
    if len(tensor) != 3 or any(part.shape != (3, 3) for part in tensor):
        raise ValueError("Canonical edge extraction requires dimension three.")
    return s.Matrix([tensor[0][1, 2], tensor[1][0, 2], tensor[2][0, 1]])


def rotate_jet_tensor(tensor, rotation):
    """Apply the full orthogonal tensor action, with the first slot sliced."""
    rotation = s.Matrix(rotation)
    n = len(tensor)
    return tuple(sum((rotation[i, j]*rotation*tensor[j]*rotation.T
                      for j in range(n)), s.zeros(n)).applyfunc(reduced)
                 for i in range(n))


def shape_trine(shape):
    """Encode a self-adjoint shape in an orthonormal two-dimensional cut frame."""
    shape = s.Matrix(shape)
    if shape.shape != (2, 2) or any(reduced(v) != 0 for v in shape-shape.T):
        raise ValueError("Expected a symmetric two-dimensional shape matrix.")
    p, t, r = shape[0, 0], shape[0, 1], shape[1, 1]
    return s.Matrix([p, (p+3*r-2*s.sqrt(3)*t)/4,
                    (p+3*r+2*s.sqrt(3)*t)/4])


def shape_from_trine(directions):
    """Recover the entire shape from the three trine directional curvatures."""
    v = s.Matrix(directions)
    if v.shape != (3, 1):
        raise ValueError("Expected three directional curvatures.")
    mixed = (v[2]-v[1])/s.sqrt(3)
    return s.Matrix([[v[0], mixed],
                     [mixed, (2*(v[1]+v[2])-v[0])/3]])


def lift_connection(lift, coordinates, orthonormal_connection=None):
    """Return Omega_mu = R^-1 (partial_mu R + omega_mu R)."""
    if orthonormal_connection is None:
        orthonormal_connection = [s.zeros(lift.rows) for _ in coordinates]
    inverse = lift.inv()
    return [(inverse*(lift.diff(x)+w*lift)).applyfunc(reduced)
            for x, w in zip(coordinates, orthonormal_connection)]


def curvature(connection, coordinates, i, j):
    return (connection[j].diff(coordinates[i])
            - connection[i].diff(coordinates[j])
            + connection[i]*connection[j]
            - connection[j]*connection[i]).applyfunc(reduced)


def torsion(coframe, connection, coordinates, i, j):
    return (coframe[:, j].diff(coordinates[i])
            - coframe[:, i].diff(coordinates[j])
            + connection[i]*coframe[:, j]
            - connection[j]*coframe[:, i]).applyfunc(reduced)


def run_checks():
    checks = {}

    def check(name, actual, expected=0):
        if isinstance(actual, s.MatrixBase) and expected == 0:
            expected = s.zeros(*actual.shape)
        difference = actual-expected
        entries = list(difference) if isinstance(difference, s.MatrixBase) else [difference]
        residuals = [reduced(entry) for entry in entries]
        if any(entry != 0 for entry in residuals):
            raise AssertionError((name, [r for r in residuals if r != 0]))
        checks[name] = {"verified": True}

    a, xi, d, gamma = s.symbols("a xi d gamma", real=True)
    p = 1+a*a
    eta = s.diag(-1, 1, 1)
    G = native_gram(a, xi)
    R = native_lift(a, xi)
    coordinates = [a, xi]
    O = lift_connection(R, coordinates)
    expected_a = s.Matrix([[0, -1, 0], [0, a, xi], [0, 0, -a]])/p
    expected_xi = s.Matrix([[0, 0, -1], [0, 0, a], [0, 0, 0]])/p
    check("resolved_gram_realization", R.T*eta*R, G)
    check("resolved_lift_volume", R.det(), 1)
    check("resolved_gram_volume", G.det(), -1)
    check("explicit_connection_a", O[0], expected_a)
    check("explicit_connection_xi", O[1], expected_xi)
    for i, x in enumerate(coordinates):
        check(f"presentation_compatibility_{x}",
              G.diff(x)-O[i].T*G-G*O[i])
        check(f"presentation_volume_preservation_{x}", s.trace(O[i]))
        for tau in (-1, 1):
            check(f"regular_horizon_transport_{x}_{tau}",
                  O[i].subs(xi, tau),
                  [expected_a, expected_xi][i].subs(xi, tau))
    check("presentation_curvature", curvature(O, coordinates, 0, 1))

    # Reuse the XLINK trine, now centered on the native clock c.
    r3 = s.sqrt(3)
    J = s.Matrix([[1/r3, 1/r3, 1/r3],
                  [2/r3, -1/r3, -1/r3], [0, 1, -1]])
    B = 2*s.eye(3)-s.ones(3)
    one = s.ones(3, 1)
    c = s.Matrix([1, 0, 0])
    U = R.inv()*J
    check("xlink_trine_normalization", J.T*eta*J, B)
    check("native_cella_embedding", U.T*G*U, B)
    check("cella_trace_seat", U*one, r3*c)
    C = s.Matrix(s.symbols("C_c C_s C_int", real=True))
    total = sum(C)
    residue = C-total*one/3
    check("channel_clock_projection", (c.T*G*U*C)[0], -total/r3)
    check("channel_full_decomposition", U*C, total*c/r3+U*residue)
    check("channel_rest_residue", (c.T*G*U*residue)[0])
    check("channel_rest_positive_form", (residue.T*U.T*G*U*residue)[0],
          2*(residue.T*residue)[0])
    for i, x in enumerate(coordinates):
        check(f"cella_embedding_parallel_{x}", U.diff(x)+O[i]*U)

    # Defining-function gauges act on Hessian/edge values, not as Lorentz
    # changes of the native frame. The full tangent jet is their quotient.
    g1, g2, g3 = s.symbols("g1 g2 g3", real=True)
    g = s.Matrix([g1, g2, g3])
    q = (g.T*g)[0]
    h11, h22, h33, h12, h13, h23 = s.symbols(
        "h11 h22 h33 h12 h13 h23", real=True)
    H = s.Matrix([[h11, h12, h13], [h12, h22, h23], [h13, h23, h33]])
    u = s.Matrix(s.symbols("u1 u2 u3", real=True))
    Hprime = H+g*u.T+u*g.T
    projector = s.eye(3)-g*g.T/q
    normalization_gauge = -H*g/q+(g.T*H*g)[0]*g/(2*q*q)
    check("cella_tangent_jet_gauge_invariance", projector*Hprime*projector,
          projector*H*projector)
    check("cella_canonical_point_jet", H+g*normalization_gauge.T
          +normalization_gauge*g.T, projector*H*projector)
    check("cella_canonical_normal_annihilation", projector*H*projector*g)
    M = s.Matrix([[0, g1*g3, g1*g2], [g2*g3, 0, g1*g2],
                  [g2*g3, g1*g3, 0]])
    nu = s.Matrix([g1*h23, g2*h13, g3*h12])
    nuprime = s.Matrix([g1*Hprime[1, 2], g2*Hprime[0, 2],
                        g3*Hprime[0, 1]])
    check("cella_full_edge_translation", nuprime-nu, M*u)
    check("cella_edge_translation_rank", M.det(), 2*(g1*g2*g3)**2)
    gauge_t = s.symbols("gauge_t", real=True)
    keystone_nu = s.Matrix([0, 0, 2])
    shifted_nu = keystone_nu+s.Matrix([3*gauge_t, 3*gauge_t, 0])
    check("cella_gauge_is_not_lorentz_transport",
          (shifted_nu.T*B*shifted_nu-keystone_nu.T*B*keystone_nu)[0],
          -24*gauge_t)
    keystone_g = s.Matrix([3, 1, 2])
    keystone_H = s.Matrix([[2, 1, 0], [1, 0, 0], [0, 0, 2]])
    keystone_P = s.eye(3)-keystone_g*keystone_g.T/14
    check("canonical_point_gauges_do_not_fix_a_whole_surface",
          -keystone_P*keystone_H*keystone_g/14,
          s.Matrix([-1, -5, 4])/98)

    # Full nonzero defining-function gauge, including a negative scale.
    # The invariant object is (P,T), not PHP without its scaling factor.
    invariant = invariant_hypersurface_jet(g, H)
    tensor = invariant["tensor"]
    gauge_covector = invariant["gauge_covector"]
    mu = s.symbols("mu", real=True, nonzero=True)
    scaled = invariant_hypersurface_jet(mu*g, mu*H+g*u.T+u*g.T)
    check("full_defining_gauge_projector", scaled["projector"], projector)
    check("full_defining_gauge_gradient_norm", scaled["q"], mu**2*q)
    check("full_defining_gauge_tangent_hessian", scaled["tangent_hessian"],
          mu*invariant["tangent_hessian"])
    check("full_defining_gauge_remainder", scaled["gauge_covector"],
          gauge_covector+u/mu)
    for i in range(3):
        check(f"full_defining_gauge_tensor_{i}", scaled["tensor"][i], tensor[i])
        check(f"invariant_tensor_symmetric_{i}", tensor[i], tensor[i].T)
        check(f"invariant_tensor_tangent_inputs_{i}", tensor[i]*g)
        check(f"invariant_tensor_normal_output_{i}",
              sum((projector[i, j]*tensor[j] for j in range(3)), s.zeros(3)))
    reconstructed = sum((g[i]*tensor[i] for i in range(3)), s.zeros(3))
    reconstructed += g*gauge_covector.T+gauge_covector*g.T
    check("full_jet_reconstruction_from_quotient_and_gauge", reconstructed, H)
    nubar = canonical_edge_vector(tensor)
    nubar_scaled = canonical_edge_vector(scaled["tensor"])
    check("canonical_edges_full_gauge_invariance", nubar_scaled, nubar)
    check("raw_edges_exact_geometric_plus_gauge_split", nu, q*nubar+M*gauge_covector)
    scaled_raw = s.Matrix([
        mu*g1*(mu*H+g*u.T+u*g.T)[1, 2],
        mu*g2*(mu*H+g*u.T+u*g.T)[0, 2],
        mu*g3*(mu*H+g*u.T+u*g.T)[0, 1]])
    check("raw_edges_full_gauge_action", scaled_raw, mu**2*nu+mu*M*u)
    check("native_canonical_edge_gauge_invariance", U*(nubar_scaled-nubar))
    check("complete_tensor_euclidean_norm", sum(s.trace(part.T*part) for part in tensor),
          s.trace(invariant["tangent_hessian"]**2)/q)

    # Actual relation differentiation in n=4 verifies the II interpretation
    # at nonzero graph slopes. Negative, nonconstant mu is differentiated
    # before the quotient is taken; no target tensor is inserted as input.
    x1, x2, x3, x4 = s.symbols("x1 x2 x3 x4", real=True)
    ambient_coordinates = (x1, x2, x3, x4)
    graph_coordinates = (x1, x2, x3)
    graph_hessian = s.Matrix([[1, 2, 1], [2, 3, -1], [1, -1, 5]])
    graph_x = s.Matrix(graph_coordinates)
    graph_f = (graph_x.T*graph_hessian*graph_x)[0]/2+x1/2+x2/3+x3/4
    relation = x4-graph_f
    embedding = s.Matrix([x1, x2, x3, graph_f])
    origin = dict.fromkeys(ambient_coordinates, 0)
    graph_gradient = s.Matrix([s.diff(relation, x).subs(origin) for x in ambient_coordinates])
    graph_ambient_hessian = s.hessian(relation, ambient_coordinates).subs(origin)
    graph_jet = invariant_hypersurface_jet(graph_gradient, graph_ambient_hessian)
    graph_tangent = embedding.jacobian(graph_coordinates).subs(origin)
    for i in range(3):
        for j in range(i, 3):
            minus_B = s.Matrix([(graph_tangent[:, i].T*part*graph_tangent[:, j])[0]
                                for part in graph_jet["tensor"]])
            direct_B = (s.eye(4)-graph_jet["projector"])*embedding.diff(
                graph_coordinates[i], graph_coordinates[j]).subs(origin)
            check(f"general_dimension_tensor_is_minus_embedded_II_{i}_{j}", minus_B, -direct_B)
    negative_mu = -2+x1+2*x4+x1*x3
    scaled_relation = negative_mu*relation
    negative_gradient = s.Matrix([s.diff(scaled_relation, x).subs(origin)
                                  for x in ambient_coordinates])
    negative_hessian = s.hessian(scaled_relation, ambient_coordinates).subs(origin)
    negative_jet = invariant_hypersurface_jet(negative_gradient, negative_hessian)
    check("negative_nonconstant_gauge_projector", negative_jet["projector"], graph_jet["projector"])
    for i in range(4):
        check(f"negative_nonconstant_gauge_full_tensor_{i}",
              negative_jet["tensor"][i], graph_jet["tensor"][i])

    # An O(n) frame change acts on all tensor indices. Exact rational
    # rotations with nonzero initial graph slopes avoid an axis-only test.
    rotation4 = s.Matrix([[s.Rational(3, 5), 0, 0, s.Rational(4, 5)],
                          [0, 0, 1, 0], [0, -1, 0, 0],
                          [-s.Rational(4, 5), 0, 0, s.Rational(3, 5)]])
    rotated4 = invariant_hypersurface_jet(rotation4*graph_gradient,
                                         rotation4*graph_ambient_hessian*rotation4.T)
    transformed4 = rotate_jet_tensor(graph_jet["tensor"], rotation4)
    check("full_tensor_rotation_is_orthogonal", rotation4.T*rotation4, s.eye(4))
    check("full_tensor_projector_covariance", rotated4["projector"],
          rotation4*graph_jet["projector"]*rotation4.T)
    for i in range(4):
        check(f"full_tensor_orthogonal_covariance_{i}", rotated4["tensor"][i], transformed4[i])

    # Edge extraction alone loses geometric information and cannot furnish
    # a Lorentz representation of all presentation-axis rotations.
    saddle_g = s.Matrix([0, 0, 1])
    saddle_H = s.diag(1, -1, 0)
    saddle = invariant_hypersurface_jet(saddle_g, saddle_H)
    tangent_rotation = s.Matrix([[1/s.sqrt(2), 1/s.sqrt(2), 0],
                                 [-1/s.sqrt(2), 1/s.sqrt(2), 0], [0, 0, 1]])
    rotated_saddle = invariant_hypersurface_jet(tangent_rotation*saddle_g,
                                               tangent_rotation*saddle_H*tangent_rotation.T)
    saddle_edges = canonical_edge_vector(saddle["tensor"])
    rotated_edges = canonical_edge_vector(rotated_saddle["tensor"])
    check("edge_extraction_omits_nonzero_shape", saddle_edges, s.zeros(3, 1))
    check("same_shape_rotated_edges_nonzero", rotated_edges, s.Matrix([0, 0, -1]))
    check("edge_B_norm_is_not_full_axis_invariant", (rotated_edges.T*B*rotated_edges)[0]
          -(saddle_edges.T*B*saddle_edges)[0], 1)
    check("complete_tensor_retains_saddle_curvature", sum(s.trace(part.T*part)
          for part in saddle["tensor"]), 2)
    check("complete_tensor_saddle_norm_axis_invariant", sum(s.trace(part.T*part)
          for part in rotated_saddle["tensor"]), 2)

    # T=0 does not determine a flat point's tangent plane. The pair (P,T)
    # is therefore required for a complete quotient, including flat jets.
    flat1 = invariant_hypersurface_jet(s.Matrix([1, 0, 0]), s.zeros(3))
    flat2 = invariant_hypersurface_jet(s.Matrix([0, 1, 0]), s.zeros(3))
    check("zero_tensor_retains_no_tangent_plane", sum(flat1["tensor"]+flat2["tensor"], s.zeros(3)))
    check("flat_jet_projector_distinguishes_gauge_orbits", flat1["projector"]-flat2["projector"],
          s.diag(-1, 1, 0))
    keystone_jet = invariant_hypersurface_jet(keystone_g, keystone_H)
    keystone_canonical_edges = canonical_edge_vector(keystone_jet["tensor"])
    keystone_native_vector = U*keystone_canonical_edges
    check("canonical_native_edge_norm", (keystone_native_vector.T*G*keystone_native_vector)[0],
          (keystone_canonical_edges.T*B*keystone_canonical_edges)[0])

    # The whole shape, unlike three ambient edge components, has an invertible
    # trine encoding and an induced Lorentz action for every tangent rotation.
    shp, sht, shr, thp, tht, thr = s.symbols("shp sht shr thp tht thr", real=True)
    shape = s.Matrix([[shp, sht], [sht, shr]])
    other_shape = s.Matrix([[thp, tht], [tht, thr]])
    chi = shape_trine(shape)
    psi = shape_trine(other_shape)
    cut_trine = [s.Matrix([1, 0]), s.Matrix([-1, s.sqrt(3)])/2,
                 s.Matrix([-1, -s.sqrt(3)])/2]
    check("shape_trine_is_directional_curvature", chi,
          s.Matrix([(v.T*shape*v)[0] for v in cut_trine]))
    check("shape_trine_full_inverse", shape_from_trine(chi), shape)
    arbitrary_chi = s.Matrix(s.symbols("chi1:4", real=True))
    check("shape_trine_full_surjection", shape_trine(shape_from_trine(arbitrary_chi)), arbitrary_chi)
    check("shape_trine_trace", sum(chi), 3*s.trace(shape)/2)
    check("shape_trine_determinant_form", (chi.T*B*chi)[0], -3*shape.det())
    check("shape_trine_polarized_determinant", (chi.T*B*psi)[0],
          -3*(s.trace(shape)*s.trace(other_shape)-s.trace(shape*other_shape))/2)
    check("shape_trine_native_components", R*U*chi,
          s.sqrt(3)*s.Matrix([(shp+shr)/2, (shp-shr)/2, -sht]))
    check("shape_trine_native_mean_clock", (c.T*G*U*chi)[0],
          -s.sqrt(3)*s.trace(shape)/2)
    angle = s.symbols("tangent_angle", real=True)
    rotate2 = s.Matrix([[s.cos(angle), s.sin(angle)],
                        [-s.sin(angle), s.cos(angle)]])
    action = s.Matrix.hstack(*[
        shape_trine(rotate2*shape_from_trine(s.eye(3)[:, i])*rotate2.T)
        for i in range(3)])
    check("shape_trine_full_tangent_rotation", shape_trine(rotate2*shape*rotate2.T), action*chi)
    check("shape_trine_rotation_preserves_B", action.T*B*action, B)
    check("shape_trine_rotation_preserves_trace", one.T*action, one.T)
    check("shape_trine_rotation_preserves_mean_section", action*one, one)
    reflection = s.diag(1, -1)
    reflected_action = s.Matrix.hstack(*[
        shape_trine(reflection*shape_from_trine(s.eye(3)[:, i])*reflection.T)
        for i in range(3)])
    check("shape_trine_reflection_preserves_B", reflected_action.T*B*reflected_action, B)
    normal1, normal2 = s.symbols("normal1 normal2", real=True)
    check("shape_trine_complete_normal_pencil",
          shape_trine(normal1*shape+normal2*other_shape), normal1*chi+normal2*psi)
    check("shape_trine_retains_saddle",
          (shape_trine(s.diag(1, -1)).T*B*shape_trine(s.diag(1, -1)))[0], 3)

    # Complete metric-compatible connection family and its curvature.
    z01, z02, z12 = s.symbols("z01 z02 z12", real=True)
    A = s.Matrix([[0, z01, z02], [-z01, 0, z12], [-z02, -z12, 0]])
    for x in coordinates:
        general = G.inv()*(G.diff(x)/2+A)
        check(f"all_metric_compatible_connections_{x}",
              G.diff(x)-general.T*G-G*general)
    boost01 = s.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    rotation12 = s.Matrix([[0, 0, 0], [0, 0, 1], [0, -1, 0]])
    omega = [xi*boost01, a*rotation12]
    curved = lift_connection(R, coordinates, omega)
    for i, x in enumerate(coordinates):
        check(f"orthonormal_free_connection_{x}",
              omega[i].T*eta+eta*omega[i])
        check(f"curved_native_metric_compatibility_{x}",
              G.diff(x)-curved[i].T*G-G*curved[i])
    check("curvature_full_conjugacy", curvature(curved, coordinates, 0, 1),
          R.inv()*curvature(omega, coordinates, 0, 1)*R)
    symmetric = [G.inv()*G.diff(x)/2 for x in coordinates]
    symmetric_curvature = curvature(symmetric, coordinates, 0, 1)
    Maurer = [G.inv()*G.diff(x) for x in coordinates]
    check("half_metric_connection_curvature", symmetric_curvature,
          -(Maurer[0]*Maurer[1]-Maurer[1]*Maurer[0])/4)
    check("half_metric_nonflat_witness", symmetric_curvature.subs({a: 0, xi: 0}),
          rotation12/4)

    # A nonconstant, nonorthogonal native basis change checks the full
    # affine connection transformation, including the skew-form correction.
    v = s.symbols("v", real=True)
    S = s.Matrix([[1, v, 0], [0, 1, 0], [0, 0, 1]])
    Gprime = S.T*G*S
    gauge_coordinates = [a, xi, v]
    original_connection = curved+[s.zeros(3)]
    transformed = []
    for x, old in zip(gauge_coordinates, original_connection):
        new = S.inv()*old*S+S.inv()*S.diff(x)
        transformed.append(new)
        check(f"basis_covariance_metric_{x}",
              Gprime.diff(x)-new.T*Gprime-Gprime*new)
        old_A = G*old-G.diff(x)/2
        check(f"basis_covariance_skew_parameter_{x}",
              Gprime*new-Gprime.diff(x)/2,
              S.T*old_A*S+(S.T*G*S.diff(x)-S.diff(x).T*G*S)/2)
    check("basis_covariance_curvature", curvature(transformed, gauge_coordinates, 0, 1),
          S.inv()*curvature(curved, coordinates, 0, 1)*S)
    check("basis_covariance_mixed_pure_gauge", curvature(transformed, gauge_coordinates, 0, 2))

    raw_change = s.Matrix([[1, 0, 0], [0, 1, gamma], [0, 0, d]])
    raw_R = R*raw_change
    check("raw_native_volume", raw_R.det(), d)
    check("raw_frame_rank_loss", raw_R.subs(d, 0).det())
    check("raw_transport_volume_pole",
          s.trace(raw_change.inv()*raw_change.diff(d)), 1/d)
    raw_G = raw_change.T*G*raw_change
    check("raw_gram_recovery", raw_G,
          s.Matrix([[-1, a, a*gamma+xi*d], [a, 1, gamma],
                    [a*gamma+xi*d, gamma, gamma**2+d*d*(1-xi*xi)/p]]))

    # The full orthonormal-frame Levi-Civita formula, checked for arbitrary
    # structure coefficients antisymmetric in the bracket indices.
    structure = s.MutableDenseNDimArray.zeros(3, 3, 3)
    for i in range(3):
        for j in range(3):
            for l in range(j+1, 3):
                structure[i, j, l] = s.Symbol(f"C{i}{j}{l}", real=True)
                structure[i, l, j] = -structure[i, j, l]
    coefficients = s.MutableDenseNDimArray.zeros(3, 3, 3)
    for i in range(3):
        for j in range(3):
            for l in range(3):
                coefficients[i, j, l] = (structure[i, l, j]
                    -structure[l, j, i]+structure[j, i, l])/2
    check("general_levi_civita_formula_metric_compatibility", s.Matrix([
        coefficients[i, j, l]+coefficients[j, i, l]
        for i in range(3) for j in range(3) for l in range(3)]))
    check("general_levi_civita_formula_zero_torsion", s.Matrix([
        coefficients[i, l, j]-coefficients[i, j, l]-structure[i, j, l]
        for i in range(3) for j in range(3) for l in range(3)]))

    # Nontrivial regular four-dimensional completion, with the same native
    # transverse pinch path as the prior load example. Its area is free.
    t, x, r, z = s.symbols("t x r z", real=True)
    k, lam = s.symbols("k lambda", positive=True)
    spacetime_coordinates = [t, x, r, z]
    state_xi = 1+k*t
    R4 = s.diag(native_lift(s.Integer(0), state_xi), 1)
    eta4 = s.diag(-1, 1, 1, 1)
    G4 = s.diag(native_gram(s.Integer(0), state_xi), 1)
    scale = s.exp(lam*t)
    theta = s.diag(1, scale, 1, scale)
    E = R4.inv()*theta
    physical_metric = s.diag(-1, scale**2, 1, scale**2)
    check("curved_soldering_metric", E.T*G4*E, physical_metric)
    check("curved_soldering_volume", E.det(), scale**2)
    check("curved_soldering_pinch_regular", E.subs(t, 0).det(), 1)
    check("same_state_null_pinch", (1-state_xi**2).subs(t, 0))
    omega4 = [s.zeros(4) for _ in spacetime_coordinates]
    for coordinate_index, spatial_index in ((1, 1), (3, 3)):
        omega4[coordinate_index][0, spatial_index] = lam*scale
        omega4[coordinate_index][spatial_index, 0] = lam*scale
    physical_connection = lift_connection(R4, spacetime_coordinates, omega4)
    flat_connection = lift_connection(R4, spacetime_coordinates)
    flat_E = R4.inv()
    check("same_native_state_flat_soldering", flat_E.T*G4*flat_E, eta4)
    for i in range(4):
        check(f"physical_compatibility_{i}", G4.diff(spacetime_coordinates[i])
              -physical_connection[i].T*G4-G4*physical_connection[i])
        for j in range(i+1, 4):
            check(f"physical_torsion_{i}_{j}",
                  torsion(E, physical_connection, spacetime_coordinates, i, j))
            check(f"flat_soldering_integrability_{i}_{j}",
                  torsion(flat_E, flat_connection, spacetime_coordinates, i, j))
            check(f"physical_curvature_conjugacy_{i}_{j}",
                  curvature(physical_connection, spacetime_coordinates, i, j),
                  R4.inv()*curvature(omega4, spacetime_coordinates, i, j)*R4)

    # Independently obtain coordinate Christoffel symbols from the metric,
    # and compare them to the connection induced by the native soldering.
    inverse_metric = physical_metric.inv()
    Christoffel = []
    induced = []
    for mu in range(4):
        matrix = s.zeros(4)
        for upper in range(4):
            for lower in range(4):
                matrix[upper, lower] = sum(
                    inverse_metric[upper, j]*(
                        s.diff(physical_metric[j, lower], spacetime_coordinates[mu])
                        +s.diff(physical_metric[j, mu], spacetime_coordinates[lower])
                        -s.diff(physical_metric[mu, lower], spacetime_coordinates[j]))/2
                    for j in range(4))
        Christoffel.append(matrix.applyfunc(reduced))
        induced.append(E.inv()*(E.diff(spacetime_coordinates[mu])
                                +physical_connection[mu]*E))
        check(f"native_connection_equals_metric_levi_civita_{mu}", induced[mu], Christoffel[mu])
    coordinate_curvature = {(i, j): curvature(Christoffel, spacetime_coordinates, i, j)
                            for i in range(4) for j in range(4)}
    Ricci = s.Matrix(4, 4, lambda i, j: sum(
        coordinate_curvature[upper, j][upper, i] for upper in range(4)))
    check("curved_example_full_ricci", Ricci,
          s.diag(-2*lam**2, 2*lam**2*scale**2, 0, 2*lam**2*scale**2))
    check("curved_example_scalar", s.trace(inverse_metric*Ricci), 6*lam**2)
    # Native k=c+w e_r has e_r=-partial_r for a=0 and xi>0.
    L = s.Matrix([(1+state_xi)/2, 0, -(1+state_xi)/2, 0])
    N = s.Matrix([1/(1+state_xi), 0, 1/(1+state_xi), 0])
    K = s.Matrix([1, 0, -state_xi, 0])
    check("example_null_normalization", (L.T*physical_metric*N)[0], -1)
    check("example_outgoing_null", (L.T*physical_metric*L)[0])
    check("example_incoming_null", (N.T*physical_metric*N)[0])
    check("example_native_normal_norm", (K.T*physical_metric*K)[0], -(1-state_xi**2))
    area = scale**2
    outgoing = sum(L[i]*s.diff(area, spacetime_coordinates[i]) for i in range(4))/area
    incoming = sum(N[i]*s.diff(area, spacetime_coordinates[i]) for i in range(4))/area
    convected = sum(K[i]*s.diff(area, spacetime_coordinates[i]) for i in range(4))/area
    check("native_area_not_selected_by_transport", convected, 2*lam)
    check("null_native_locus_not_automatically_marginal", outgoing.subs(t, 0), 2*lam)
    check("incoming_expansion_at_native_locus", incoming.subs(t, 0), lam)

    # Even a flat integrable ambient coframe need not integrate the native
    # radial screen. Here a=z, xi=1, d=t, and beta=-(z dx+dy)/sqrt(1+z^2).
    y = s.symbols("y", real=True)
    screen_coordinates = [t, x, y, z]
    alpha = s.Matrix([1, 0, 0, 0])
    beta = s.Matrix([0, -z, -1, 0])/s.sqrt(1+z*z)
    X = s.Matrix([0, 1, -z, 0])
    Y = s.Matrix([0, 0, 0, 1])
    check("screen_first_normal_independence", (alpha.T*X)[0])
    check("screen_second_normal_X", (beta.T*X)[0])
    check("screen_second_normal_Y", (beta.T*Y)[0])
    bracket = s.Matrix([sum(X[i]*s.diff(Y[j], screen_coordinates[i])
                           -Y[i]*s.diff(X[j], screen_coordinates[i]) for i in range(4))
                       for j in range(4)])
    check("screen_bracket_obstruction", (beta.T*bracket)[0], -1/s.sqrt(1+z*z))
    frobenius_coefficient = beta[1]*(s.diff(beta[3], y)-s.diff(beta[2], z)) \
        +beta[2]*(s.diff(beta[1], z)-s.diff(beta[3], x)) \
        +beta[3]*(s.diff(beta[2], x)-s.diff(beta[1], y))
    check("screen_frobenius_obstruction", frobenius_coefficient, 1/(1+z*z))

    return {
        "scope": "Complete local metric-compatible transport and soldering classification for finite resolved native state fields; Cella channel/edge encodings and integrability obstructions. Physical response and metric selection remain explicit inputs.",
        "verified_count": len(checks),
        "checks": checks,
        "results": {
            "regular_lift": "R=[[1,-a,-xi],[0,sqrt(p),a*xi/sqrt(p)],[0,0,1/sqrt(p)]], p=1+a^2, det R=1",
            "presentation_transport": "Omega0=R^-1 dR; curvature zero; regular at d=0 and xi=+/-1",
            "complete_connection_family": "Omega=R^-1(omega R+dR), omega^T eta+eta omega=0",
            "complete_soldering_family": "E=R^-1 theta for every invertible orthonormal coframe theta; g=theta^T eta theta",
            "cella_account_embedding": "U=R^-1 J, U^T G U=2I-11^T, U1=sqrt(3)c; V_C=(sum C/sqrt(3))c+U(C-mean(C)1)",
            "cella_edge_gauge": "nu'=nu+M(g)u; det M=2(g1g2g3)^2; generally not Lorentz transport",
            "torsion": "R T_E=dtheta+omega wedge theta; zero torsion uniquely fixes omega for a chosen theta",
            "coordinate_integration": "A fixed flat-frame coordinate integration exists locally exactly when connection curvature and coframe torsion vanish",
            "screen_integration": "For independent normal forms alpha,beta: alpha wedge beta wedge dalpha=alpha wedge beta wedge dbeta=0",
            "nonflat_metric_only_connection": "(1/2)G^-1 dG has curvature -[G^-1 partial_a G,G^-1 partial_xi G]/4, nonzero already at a=xi=0",
            "curved_crossing_example": "Same a=0,xi=1+kt,d=t with ds^2=-dt^2+exp(2 lambda t)(dx^2+dz^2)+dr^2; scalar curvature 6 lambda^2; outgoing expansion at t=0 is 2 lambda, not zero"
        }
    }


if __name__ == "__main__":
    print(json.dumps(run_checks(), indent=2))

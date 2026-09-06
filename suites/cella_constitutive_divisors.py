#!/usr/bin/env python3
"""Exact replay for the constitutive divisor classification and native ports.

Proofs: docs/2026-09-05-cella-constitutive-divisors.md.
JSON goes to stdout; no corpus files or input data are modified.
"""

import json

import sympy as s


checks = {}


def check(name, actual, expected=0):
    difference = actual if isinstance(actual, s.MatrixBase) and expected == 0 else actual - expected
    entries = list(difference) if isinstance(difference, s.MatrixBase) else [difference]
    for entry in entries:
        residual = s.simplify(s.cancel(s.expand(entry)))
        if residual != 0:
            raise AssertionError((name, residual))
    checks[name] = {"verified": True}


def divisor_response(coordinates, incidence, positive_factor, skew):
    """Generate the theorem's complete adapted linear class from barred data."""
    n = len(incidence)
    if positive_factor.rows != n or skew.shape != (n, n):
        raise ValueError("channel dimensions disagree")
    if skew + skew.T != s.zeros(n):
        raise ValueError("barred energy-neutral part must be skew")
    weights = [s.prod(coordinates[a] for a in sorted(indices)) for indices in incidence]
    q = s.diag(*weights)
    r = q * positive_factor * positive_factor.T * q
    j = s.Matrix(n, n, lambda i, k:
                 s.prod(coordinates[a] for a in sorted(incidence[i] | incidence[k]))
                 * skew[i, k])
    return r + j, r, j, q


f, g, a, b = s.symbols("f g a b", real=True)
incidence = [{0}, {0, 1}, {1}, set()]
factor = s.Matrix([[1, 0, 0, 0], [1, 1, 0, 0],
                   [1, a, 1, 0], [1, b, a, 1]])
skew = s.Matrix([[0, 2, a, 1], [-2, 0, 3, b],
                 [-a, -3, 0, 1], [-1, -b, -1, 0]])
mobility, dissipative, neutral, q = divisor_response((f, g), incidence, factor, skew)
effort = s.Matrix(s.symbols("e0:4", real=True))
check("multi_face_dissipative_symmetry", dissipative - dissipative.T)
check("multi_face_neutral_skew", neutral + neutral.T)
check("multi_face_power_positive_factor", (effort.T * mobility * effort)[0],
      (factor.T * q * effort).dot(factor.T * q * effort))
for i, surfaces in enumerate(incidence):
    for surface in surfaces:
        coordinate = (f, g)[surface]
        check(f"row_{i}_vanishes_on_face_{surface}",
              mobility.row(i).subs(coordinate, 0), s.zeros(1, 4))
check("shared_face_dissipative_factor", dissipative[0, 1], f ** 2 * g)
check("shared_face_neutral_factor", neutral[0, 1], 2 * f * g)
check("separate_faces_neutral_factor", neutral[0, 2], a * f * g)
check("unconstrained_channel_dissipation", dissipative[3, 3], 1 + a * a + b * b + 1)

# A shared-face energy-neutral response genuinely needs only one factor.
shared = f * f * s.eye(2) + f * s.Matrix([[0, 1], [-1, 0]])
e_shared = effort[:2, :]
check("shared_face_power", (e_shared.T * shared * e_shared)[0],
      f * f * (effort[0] ** 2 + effort[1] ** 2))
check("shared_face_rows_vanish", shared.subs(f, 0), s.zeros(2))
check("shared_face_first_derivative_survives", s.diff(shared[0, 1], f).subs(f, 0), 1)
bad_dissipative = s.Matrix([[f * f, f], [f, f * f]])
check("insufficient_dissipative_order_negative_minor",
      bad_dissipative.det().subs(f, s.Rational(1, 2)), -s.Rational(3, 16))

t = s.symbols("t", positive=True)
check("weighted_dissipative_order_four",
      dissipative[0, 1].subs({f: t, g: t ** 2}) / t ** 4, 1)
check("weighted_neutral_order_three",
      neutral[0, 1].subs({f: t, g: t ** 2}) / t ** 3, 2)

# Nonlinear, integrable, incrementally passive response with a horizon channel.
e0, e1, u = s.symbols("e0 e1 u", real=True)
potential = (e0 + f * e1) ** 2 / 2 + f * f * e1 ** 4 / 4
variables = s.Matrix([e0, e1])
response = s.Matrix([s.diff(potential, e) for e in variables])
jacobian = response.jacobian(variables)
q2 = s.diag(1, f)
bar = s.Matrix([[1, 1], [1, 1 + 3 * e1 ** 2]])
bar_factor = s.Matrix([[1, 0], [1, s.sqrt(3) * e1]])
check("nonlinear_positive_hessian_factor", jacobian, q2 * bar_factor * bar_factor.T * q2)
check("nonlinear_barred_hessian", bar, bar_factor * bar_factor.T)
for i in range(2):
    check(f"nonlinear_row_{i}_integrable",
          s.diff(jacobian[i, 1], e0) - s.diff(jacobian[i, 0], e1))
recovered = s.integrate(jacobian.subs({e0: u * e0, e1: u * e1}, simultaneous=True)
                        * variables, (u, 0, 1))
check("nonlinear_radial_reconstruction", recovered, response)
check("nonlinear_horizon_output", response[1].subs(f, 0))
nonintegrable = s.eye(2) + e0 * s.Matrix([[0, 1], [-1, 0]])
check("positive_symmetric_part_without_integrability",
      (nonintegrable + nonintegrable.T) / 2, s.eye(2))
check("row_curl_obstruction",
      s.diff(nonintegrable[0, 1], e0) - s.diff(nonintegrable[0, 0], e1), 1)

# Geometric observable ports with cross-coupled area dependence.
x, d, radius = s.symbols("x d radius", real=True, nonzero=True)
state = s.Matrix([x, d, radius])
observations = s.Matrix([x, d, 2 * s.log(radius) + x * x + x * d])
c = observations.jacobian(state)
v_state = s.Matrix(s.symbols("v_x v_d v_R", real=True))
e_obs = s.Matrix(s.symbols("e_x e_d e_A", real=True))
drift = s.Matrix([0, 0, s.Symbol("area_drift", real=True)])
v_obs = c * v_state + drift
check("derived_observable_power",
      ((c.T * e_obs).T * v_state)[0], (e_obs.T * (v_obs - drift))[0])
state_r = s.Matrix([[2, 1, 0], [1, 2, 1], [0, 1, 2]])
state_j = s.Matrix([[0, 1, 2], [-1, 0, 3], [-2, -3, 0]])
effective = c * (state_r + state_j) * c.T
check("observable_symmetric_part", (effective + effective.T) / 2, c * state_r * c.T)
check("observable_neutral_part", (effective - effective.T) / 2, c * state_j * c.T)
w, radial_gradient, radius_rate, e_area = s.symbols("w radial_gradient radius_rate e_area")
area_rate = 2 * radius_rate / radius + 2 * w * radial_gradient / radius
check("spherical_area_work_with_drift", 2 * e_area * radius_rate / radius,
      e_area * (area_rate - 2 * w * radial_gradient / radius))

rho, beta, cross, rself = s.symbols("rho beta cross rself", nonzero=True, real=True)
D = f * (beta + cross * e0 + f * rself * e1)
theta = D + rho * f / 2
check("affine_horizon_orientation", s.cancel(2 * theta / (rho * f)).subs(f, 0),
      1 + 2 * (beta + cross * e0) / rho)
check("orientation_without_unbounded_effort_dependence",
      s.diff((1 + 2 * (beta + cross * e0) / rho).subs(cross, 0), e0))

# Variable material inertia: identify the momentum exchange of the DIS equation.
tau = s.symbols("tau", real=True)
alpha, material, environment = [s.Function(n)(tau) for n in ("alpha", "material", "environment")]
mass = s.Function("mass")(material, environment)
storage = s.Function("storage")(alpha, material, environment)
damping = s.Function("damping")(material, environment)
force = s.Function("force")(tau)
velocity = s.diff(alpha, tau)
acceleration = (force - damping * velocity - s.diff(storage, alpha)) / mass
energy = mass * velocity ** 2 / 2 + storage
energy_rate = s.diff(energy, tau).subs(s.diff(alpha, tau, 2), acceleration)
expected_rate = (force * velocity - damping * velocity ** 2
                 + s.diff(mass, tau) * velocity ** 2 / 2
                 + s.diff(storage, material) * s.diff(material, tau)
                 + s.diff(storage, environment) * s.diff(environment, tau))
check("DIS_variable_constitution_energy", energy_rate, expected_rate)
momentum_rate = s.diff(mass * velocity, tau).subs(s.diff(alpha, tau, 2), acceleration)
check("DIS_material_momentum_exchange", momentum_rate,
      force - damping * velocity - s.diff(storage, alpha) + s.diff(mass, tau) * velocity)

# Cella weighted resolution determines the full native response transform.
xi, mu, nu, beta_native = s.symbols("xi mu nu beta_native", real=True)
resolved = s.Matrix([[mu, beta_native], [beta_native, nu]])
pushforward = s.Matrix([[d, xi], [0, 1]])
native = pushforward * resolved * pushforward.T
expected_native = s.Matrix([
    [mu * d ** 2 + 2 * beta_native * d * xi + nu * xi ** 2,
     beta_native * d + nu * xi],
    [beta_native * d + nu * xi, nu],
])
check("native_weighted_response", native, expected_native)
check("native_weighted_response_determinant", native.det(), d ** 2 * (mu * nu - beta_native ** 2))
contrast = s.Matrix([[1, -xi], [0, 1]])
check("native_contrast_coframe", contrast * pushforward, s.diag(d, 1))
check("native_contrast_response", contrast * native * contrast.T,
      s.diag(d, 1) * resolved * s.diag(d, 1))
native_effort = s.Matrix(s.symbols("native_e_x native_e_d"))
check("native_power_pullback", (native_effort.T * native * native_effort)[0],
      ((pushforward.T * native_effort).T * resolved * (pushforward.T * native_effort))[0])
clock, xi0 = s.symbols("clock xi0", real=True)
dt = nu * clock
xit = xi0 + beta_native * clock
xt = dt * xit
native_flow = native * s.Matrix([0, 1])
check("resolved_constant_load_x_trajectory",
      native_flow[0].subs({d: dt, xi: xit}), s.diff(xt, clock))
check("resolved_constant_load_d_trajectory", native_flow[1], s.diff(dt, clock))
check("resolved_constant_load_direction",
      s.limit(xt / dt, clock, 0), xi0)
check("resolved_constant_load_direction_rate",
      s.diff(xit, clock), beta_native)
k_native = s.symbols("k_native", real=True)
check("previous_native_crossing_recovered",
      xt.subs({nu: 1, xi0: 1, beta_native: k_native}), clock + k_native * clock ** 2)
check("original_full_rank_response_lift_obstruction",
      (pushforward.inv() * s.Matrix([[2, 1], [1, 1]]) * pushforward.inv().T
       ).subs(xi, 1)[0, 0], 1 / d ** 2)
weights = s.diag(d, d ** 2, 1)
material_response = s.Matrix([[2, 1, 1], [1, 3, 1], [1, 1, 2]])
weighted_response = weights * material_response * weights
check("weighted_1_2_mixed_response_order", weighted_response[0, 1], d ** 3)
check("weighted_1_2_second_response_order", weighted_response[1, 1], 3 * d ** 4)
check("weighted_1_2_lift_recovery", weights.inv() * weighted_response * weights.inv(),
      material_response)

# Smooth descent through the collapsed fibre is stronger than a smooth lift.
a11, a12, a21 = s.symbols("a11 a12 a21", real=True)
nxi = s.Function("nxi")(xi)
general_response = s.Matrix([[a11, a12], [a21, nxi]])
collapsed = (pushforward * general_response * pushforward.T).subs(d, 0)
check("collapsed_fibre_general_response", collapsed,
      nxi * s.Matrix([[xi ** 2, xi], [xi, 1]]))
check("descent_constant_22_forces_n_constant", s.diff(collapsed[1, 1], xi), s.diff(nxi, xi))
check("descent_constant_12_forces_n_zero",
      s.diff(collapsed[0, 1], xi).subs(s.diff(nxi, xi), 0), nxi)
resolved_general_flow = general_response * pushforward.T * native_effort
check("finite_effort_radial_boundary_rate", resolved_general_flow[1].subs(d, 0),
      nxi * (xi * native_effort[0] + native_effort[1]))
check("continuous_descent_blocks_transverse_rate",
      resolved_general_flow[1].subs({d: 0, nxi: 0}), 0)

print(json.dumps({
    "scope": "Smooth passive divisor responses, nonlinear integrability and geometric power ports",
    "check_count": len(checks),
    "checks": checks,
}, indent=2))

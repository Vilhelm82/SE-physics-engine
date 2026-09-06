#!/usr/bin/env python3
"""Exact checks for docs/results/2026-09-05/2026-09-05-cella-dependency-reconciliation.md.

Cella source proofs are identified in that document. This checks the reused
coupling form, the signed diagonal-curvature extension against Christoffel
contraction, and the normal-plane trace extension. It imports no web code.
Run from any directory; JSON is written to stdout, not into the Cella corpus.
"""

import json
from itertools import combinations

import sympy as s


checks = {}


def check(name, actual, expected=0):
    difference = actual - expected
    entries = list(difference) if isinstance(difference, s.MatrixBase) else [difference]
    for entry in entries:
        residual = s.simplify(s.cancel(s.expand(entry)))
        if residual != 0:
            raise AssertionError((name, residual))
    checks[name] = {"verified": True}


def diagonal_scalar_direct(diagonal, coordinates):
    """Contract the coordinate Christoffel/Ricci definition, independently of Q_j."""
    n = len(coordinates)
    gamma = [[[s.S.Zero for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                numerator = s.S.Zero
                if k == i:
                    numerator += s.diff(diagonal[k], coordinates[j])
                if k == j:
                    numerator += s.diff(diagonal[k], coordinates[i])
                if i == j:
                    numerator -= s.diff(diagonal[i], coordinates[k])
                gamma[k][i][j] = numerator / (2 * diagonal[k])
    scalar = s.S.Zero
    for i in range(n):
        ricci = s.S.Zero
        for k in range(n):
            ricci += s.diff(gamma[k][i][i], coordinates[k])
            ricci -= s.diff(gamma[k][i][k], coordinates[i])
            for l in range(n):
                ricci += gamma[k][k][l] * gamma[l][i][i]
                ricci -= gamma[k][i][l] * gamma[l][i][k]
        scalar += ricci / diagonal[i]
    return scalar


def diagonal_channels(diagonal, coordinates):
    """Cella directional channels, using u_i=log(abs(g_i))/2 on a fixed-sign patch."""
    n = len(coordinates)
    channels = []
    for j, coordinate in enumerate(coordinates):
        first = [s.diff(g, coordinate) / (2 * g) for g in diagonal]
        qj = sum(
            s.diff(first[i], coordinate) + first[i] ** 2 - first[j] * first[i]
            for i in range(n) if i != j
        )
        qj += sum(first[r] * first[t] for r, t in combinations(range(n), 2)
                  if r != j and t != j)
        channels.append(-2 * qj / diagonal[j])
    return channels


# Cella's coupling determinant supplies the indefinite edge form algebraically.
g1, g2, g3, h12, h13, h23 = s.symbols("g1 g2 g3 h12 h13 h23", real=True)
gradient = s.Matrix([g1, g2, g3])
hc = s.Matrix([[0, h12, h13], [h12, 0, h23], [h13, h23, 0]])
bordered = s.zeros(4)
bordered[0, 1:] = gradient.T
bordered[1:, 0] = gradient
bordered[1:, 1:] = hc
nu = s.Matrix([g1 * h23, g2 * h13, g3 * h12])
edge_form = 2 * s.eye(3) - s.ones(3)
check("cella_coupling_determinant", bordered.det(), (nu.T * edge_form * nu)[0])
e0 = s.ones(3, 1) / s.sqrt(3)
e1 = s.Matrix([1, -1, 0]) / s.sqrt(2)
e2 = s.Matrix([1, 1, -2]) / s.sqrt(6)
isometry = s.Matrix.hstack(e0, e1 / s.sqrt(2), e2 / s.sqrt(2))
check("cella_edge_form_native_signature", isometry.T * edge_form * isometry,
      s.diag(-1, 1, 1))

# A generic defining-function gauge moves the edge vector. This identifies the
# exact presentation-selection target instead of silently quotienting it away.
b1, b2, b3 = s.symbols("b1 b2 b3", real=True)
gauge_jet = s.Matrix([b1, b2, b3])
dh = gradient * gauge_jet.T + gauge_jet * gradient.T
dnu = s.Matrix([g1 * dh[1, 2], g2 * dh[0, 2], g3 * dh[0, 1]])
check("edge_gauge_map_determinant", dnu.jacobian(gauge_jet).det(),
      2 * (g1 * g2 * g3) ** 2)

# Arbitrary nonzero diagonal components; these identities include every fixed
# signature and arbitrary mixed coordinate dependence in dimensions 2, 3 and 4.
for dimension in (2, 3, 4):
    coordinates = s.symbols(f"t0:{dimension}", real=True)
    diagonal = [s.Function(f"g{i}")(*coordinates) for i in range(dimension)]
    check(f"signed_curvature_arbitrary_functions_n{dimension}",
          diagonal_scalar_direct(diagonal, coordinates),
          sum(diagonal_channels(diagonal, coordinates)))

# The completed load witness, on R>0 and away from spherical coordinate poles.
T, r, theta, phi, k = s.symbols("T r theta phi k", real=True)
radius = r - T - k * T ** 2 / 2
xi = 1 + k * T
diagonal = [-s.S.One, s.S.One, radius ** 2, radius ** 2 * s.sin(theta) ** 2]
coordinates = (T, r, theta, phi)
witness_channels = diagonal_channels(diagonal, coordinates)
expected_scalar = 2 * xi ** 2 / radius ** 2 - 4 * k / radius
check("load_witness_cella_scalar", sum(witness_channels), expected_scalar)
check("load_witness_direct_scalar", diagonal_scalar_direct(diagonal, coordinates),
      expected_scalar)
check("load_witness_horizon_scalar", expected_scalar.subs(T, 0),
      2 / r ** 2 - 4 * k / r)

# Reuse the weighted-jet reflection coefficient with either sign of its units.
x, y, z = s.symbols("x y z", real=True)
B, C1, C2 = s.symbols("B C1 C2", nonzero=True, real=True)
reflection = (B * x ** 2, C1 / x ** 2, C2 / x ** 2)
check("signed_reflection_existing_lead7_coefficient",
      sum(diagonal_channels(reflection, (x, y, z))), -14 / (B * x ** 4))

# Area variation on an arbitrary positive two-metric, allowing shear and
# off-diagonal metric components. Differentiate its determinant directly.
t = s.symbols("t", real=True)
q11, q12, q22 = [s.Function(name)(t) for name in ("q11", "q12", "q22")]
q = s.Matrix([[q11, q12], [q12, q22]])
check("area_log_determinant_trace", s.diff(q.det(), t) / (2 * q.det()),
      s.trace(q.inv() * q.diff(t)) / 2)
c11, c12, c22, b11, b12, b22 = s.symbols("c11 c12 c22 b11 b12 b22")
bc = s.Matrix([[c11, c12], [c12, c22]])
bs = s.Matrix([[b11, b12], [b12, b22]])
check("area_channel_trace_additivity", s.trace(q.inv() * (bc + bs)),
      s.trace(q.inv() * bc) + s.trace(q.inv() * bs))

F, rho, D, rho1, rho2, D1, D2 = s.symbols("F rho D rho1 rho2 D1 D2", real=True)
normal_metric = s.Matrix([[-F, 1], [1, 0]])
mean_dual = s.Matrix([rho, D + rho * F])
area_normal = s.Matrix([rho, -D])
check("normal_area_covector", normal_metric * mean_dual, s.Matrix([D, rho]))
check("normal_area_preservation", (s.Matrix([[D, rho]]) * area_normal)[0])
check("dual_normal_orthogonality", (mean_dual.T * normal_metric * area_normal)[0])
check("dual_normal_causal_norm", (area_normal.T * normal_metric * area_normal)[0],
      -rho ** 2 * F - 2 * rho * D)
check("dual_normal_channel_additivity", area_normal.subs({rho: rho1 + rho2, D: D1 + D2}),
      s.Matrix([rho1, -D1]) + s.Matrix([rho2, -D2]))

print(json.dumps({
    "scope": "Cella algebra and conditional horizon extensions; no native physical closure asserted",
    "check_count": len(checks),
    "checks": checks,
    "witness_scalar": str(expected_scalar),
}, indent=2))

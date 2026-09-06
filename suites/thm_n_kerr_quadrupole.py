#!/usr/bin/env python3
"""THM-N: Kerr O(J^2) falsifier for the minimal THM-M continuation.

The only continuation tested here is the one already available in THM-M:

1. each source element carries the static first-order presented field
   ``h_00 = rho`` and ``h_ij = rho n_i n_j``;
2. that rank-two field is Lorentz-boosted with the element;
3. source elements superpose linearly.

No Kerr profile, quadrupole, field equation, source-stress correction, or
second-order rescue term is admitted.  A uniformly rotating thin spherical
shell is used because THM-M used the same source and its first-order far field
depends only on J.

At second order in ``epsilon = Omega R / c``, the boosted shell develops a
direction-dependent ``P_2(cos Theta) / r`` term.  This happens both under the
literal component-only continuation of THM-M and when the field argument is
also transformed to each element's instantaneous rest frame.  Kerr has a
scalar mass monopole at order ``1/r`` in asymptotically Cartesian mass-centred
coordinates, so the minimal continuation is killed before its formal
``1/r^3`` quadrupole can count as a Kerr prediction.  If that forbidden term
is ignored, the formal quadrupole is only 3/10 of Kerr's invariant value.
"""

from __future__ import annotations

import sys

import sympy as sp


def _sphere_average(expr: sp.Expr, u: sp.Symbol, phi: sp.Symbol) -> sp.Expr:
    """Exact uniform average over a unit sphere using u = cos(theta)."""
    return sp.simplify(
        sp.integrate(sp.integrate(sp.expand(expr), (phi, 0, 2 * sp.pi)), (u, -1, 1))
        / (4 * sp.pi)
    )


def _legendre_decompose_quadratic(expr: sp.Expr, q: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    """Return c0, c2 for an axisymmetric quadratic c0 + c2 P2(q)."""
    c0, c2 = sp.symbols("c0 c2")
    p2 = (3 * q**2 - 1) / 2
    poly = sp.Poly(sp.expand(expr - c0 - c2 * p2), q)
    solution = sp.solve(poly.all_coeffs(), (c0, c2), dict=True)
    if len(solution) != 1:
        raise ValueError(f"expression does not have a unique l=0,2 decomposition: {expr}")
    return sp.simplify(solution[0][c0]), sp.simplify(solution[0][c2])


def derive_local_boost_weights() -> dict[str, sp.Expr]:
    """Derive the O(beta^2) component and field-argument weights exactly."""
    t, b_squared, b_dot_n_squared = sp.symbols("t B C", real=True)
    gamma_squared = 1 / (1 - t * b_squared)

    component_factor = gamma_squared * (1 + t * b_dot_n_squared)

    # At lab simultaneity, an event separation transforms to the element's
    # rest frame with s'^2/s^2 = 1 + gamma^2 (beta.n)^2.  The parallel unit
    # direction transforms by the same gamma before it is renormalised.
    rest_distance_squared = 1 + gamma_squared * t * b_dot_n_squared
    beta_dot_n_rest_squared = (
        gamma_squared * t * b_dot_n_squared / rest_distance_squared
    )
    covariant_factor = (
        gamma_squared
        * (1 + beta_dot_n_rest_squared)
        / sp.sqrt(rest_distance_squared)
    )

    component_t2 = sp.simplify(sp.diff(component_factor, t).subs(t, 0))
    covariant_t2 = sp.simplify(sp.diff(covariant_factor, t).subs(t, 0))
    return {
        "component_b_squared": sp.diff(component_t2, b_squared),
        "component_b_dot_n_squared": sp.diff(component_t2, b_dot_n_squared),
        "covariant_b_squared": sp.diff(covariant_t2, b_squared),
        "covariant_b_dot_n_squared": sp.diff(covariant_t2, b_dot_n_squared),
    }


def derive_second_order_coefficients() -> dict[str, sp.Expr]:
    """Derive the exact boosted-shell coefficients by direct sphere integration.

    Write the field direction as ``N_z = q`` and the source direction as ``e``.
    With ``a = N.e``, ``b = z-hat x e``, ``delta = R/r``, and
    ``S = r/|r N - R e|``, the Lorentz-transformed time component is

        F = <S> + epsilon^2 [<b^2 S> + <(b.N)^2 S^3>] + O(epsilon^4),

    where ``F = h'_00 / (2GM/(c^2 r))``.  The two required radial expansions are

        S   = 1 + delta a + delta^2(-1/2 + 3a^2/2) + ...,
        S^3 = 1 + 3delta a + delta^2(-3/2 + 15a^2/2) + ....

    Odd powers average to zero.  The returned ``r_inv`` coefficients multiply
    ``epsilon^2/r``; ``p2_r_inv3`` multiplies
    ``epsilon^2 R^2 P2(cos Theta)/r^3``.
    """
    u, phi, q = sp.symbols("u phi q", real=True)
    source_xy = sp.sqrt(1 - u**2)
    field_xy = sp.sqrt(1 - q**2)

    a = field_xy * source_xy * sp.cos(phi) + q * u
    b_squared = 1 - u**2
    b_dot_n_far_squared = (1 - q**2) * (1 - u**2) * sp.sin(phi) ** 2

    epsilon2_delta0 = _sphere_average(
        b_squared + b_dot_n_far_squared,
        u,
        phi,
    )
    epsilon2_delta2 = _sphere_average(
        b_squared * (-sp.Rational(1, 2) + sp.Rational(3, 2) * a**2)
        + b_dot_n_far_squared
        * (-sp.Rational(3, 2) + sp.Rational(15, 2) * a**2),
        u,
        phi,
    )

    r_inv_l0, r_inv_l2 = _legendre_decompose_quadratic(epsilon2_delta0, q)
    r_inv3_l0, r_inv3_l2 = _legendre_decompose_quadratic(epsilon2_delta2, q)

    return {
        "isotropic_r_inv": r_inv_l0,
        "p2_r_inv": r_inv_l2,
        "isotropic_r_inv3": r_inv3_l0,
        "p2_r_inv3": r_inv3_l2,
    }


def derive_covariant_argument_coefficients() -> dict[str, sp.Expr]:
    """Repeat the average after transforming the field argument as well.

    At lab simultaneity, the source-rest-frame distance is
    ``s' = s[1 + (beta.n)^2/2] + O(beta^4)``.  Combining this with the
    component transformation halves the ``(b.N)^2 S^3`` contribution.  This
    is the strongest minimal reading of BARE-1 covariance available without a
    retarded field or a new source-stress law.
    """
    u, phi, q = sp.symbols("u phi q", real=True)
    source_xy = sp.sqrt(1 - u**2)
    field_xy = sp.sqrt(1 - q**2)

    a = field_xy * source_xy * sp.cos(phi) + q * u
    b_squared = 1 - u**2
    b_dot_n_far_squared = (1 - q**2) * (1 - u**2) * sp.sin(phi) ** 2

    epsilon2_delta0 = _sphere_average(
        b_squared + sp.Rational(1, 2) * b_dot_n_far_squared,
        u,
        phi,
    )
    epsilon2_delta2 = _sphere_average(
        b_squared * (-sp.Rational(1, 2) + sp.Rational(3, 2) * a**2)
        + sp.Rational(1, 2)
        * b_dot_n_far_squared
        * (-sp.Rational(3, 2) + sp.Rational(15, 2) * a**2),
        u,
        phi,
    )

    r_inv_l0, r_inv_l2 = _legendre_decompose_quadratic(epsilon2_delta0, q)
    r_inv3_l0, r_inv3_l2 = _legendre_decompose_quadratic(epsilon2_delta2, q)
    return {
        "isotropic_r_inv": r_inv_l0,
        "p2_r_inv": r_inv_l2,
        "isotropic_r_inv3": r_inv3_l0,
        "p2_r_inv3": r_inv3_l2,
    }


def derive_second_order_coefficients_by_moments() -> dict[str, sp.Expr]:
    """Independent isotropic-tensor route for the same two coefficients."""
    q = sp.symbols("q", real=True)
    p2 = (3 * q**2 - 1) / 2

    # <b^2> = 2/3 and <(b.N)^2> = (1-q^2)/3.
    delta0 = sp.Rational(2, 3) + (1 - q**2) / 3

    # <b^2 a^2> = 2(2-q^2)/15.  For m = N x z,
    # <(m.e)^2> = (1-q^2)/3 and <(m.e)^2 (N.e)^2> = (1-q^2)/15.
    clock_delta2 = (
        -sp.Rational(1, 2) * sp.Rational(2, 3)
        + sp.Rational(3, 2) * (2 * (2 - q**2) / 15)
    )
    rods_delta2 = (
        -sp.Rational(3, 2) * ((1 - q**2) / 3)
        + sp.Rational(15, 2) * ((1 - q**2) / 15)
    )
    delta2 = sp.simplify(clock_delta2 + rods_delta2)

    c0_0, c2_0 = _legendre_decompose_quadratic(delta0, q)
    c0_2, c2_2 = _legendre_decompose_quadratic(delta2, q)
    if sp.simplify(delta2 + sp.Rational(2, 15) * p2) != 0:
        raise AssertionError("moment route did not reduce to the expected P2 sector")
    return {
        "isotropic_r_inv": c0_0,
        "p2_r_inv": c2_0,
        "isotropic_r_inv3": c0_2,
        "p2_r_inv3": c2_2,
    }


def evaluate_kerr_gate() -> dict[str, object]:
    """Evaluate the asymptotic and formal Kerr-quadrupole gates."""
    coeffs = derive_second_order_coefficients()
    covariant_argument_coeffs = derive_covariant_argument_coefficients()

    # Q_model = -(2/15) M Omega^2 R^4/c^2.
    # J_shell = (2/3) M R^2 Omega, while Q_Kerr = -J^2/(M c^2),
    # hence Q_model/Q_Kerr = (2/15)/(4/9) = 3/10.
    formal_ratio = sp.simplify(sp.Rational(2, 15) / sp.Rational(4, 9))
    asymptotic_pass = bool(coeffs["p2_r_inv"] == 0)
    covariant_argument_asymptotic_pass = bool(covariant_argument_coeffs["p2_r_inv"] == 0)
    quadrupole_pass = bool(formal_ratio == 1)
    return {
        "asymptotic_monopole_passes": asymptotic_pass,
        "covariant_argument_asymptotic_passes": covariant_argument_asymptotic_pass,
        "formal_quadrupole_passes": quadrupole_pass,
        "formal_q_over_kerr_q": formal_ratio,
        "status": (
            "PASSES_KERR"
            if asymptotic_pass and covariant_argument_asymptotic_pass and quadrupole_pass
            else "KILLED_BEFORE_QUADRUPOLE"
        ),
    }


def main() -> int:
    direct = derive_second_order_coefficients()
    covariant_argument = derive_covariant_argument_coefficients()
    boost_weights = derive_local_boost_weights()
    moments = derive_second_order_coefficients_by_moments()
    verdict = evaluate_kerr_gate()

    checks = [
        ("direct sphere integration agrees with the independent tensor-moment route", direct == moments),
        (
            "the exact spacetime boost makes the transformed-argument rod weight 1/2",
            boost_weights["covariant_b_dot_n_squared"] == sp.Rational(1, 2),
        ),
        ("isotropic O(epsilon^2/r) mass renormalisation is 8/9", direct["isotropic_r_inv"] == sp.Rational(8, 9)),
        ("forbidden P2/r coefficient is -2/9, not zero", direct["p2_r_inv"] == -sp.Rational(2, 9)),
        (
            "transforming the field argument leaves a forbidden P2/r coefficient -1/9",
            covariant_argument["p2_r_inv"] == -sp.Rational(1, 9),
        ),
        ("formal P2/r^3 coefficient is -2/15", direct["p2_r_inv3"] == -sp.Rational(2, 15)),
        ("minimal continuation fails the asymptotic Kerr gate", verdict["asymptotic_monopole_passes"] is False),
        ("formal quadrupole is 3/10 of Kerr and also fails", verdict["formal_q_over_kerr_q"] == sp.Rational(3, 10)),
    ]

    print("=" * 100)
    print("THM-N  Kerr O(J^2) kill test: minimal continuation of THM-M")
    print("=" * 100)
    for label, passed in checks:
        print(f"  [{'PASS' if passed else 'FAIL'}] {label}")

    print()
    print("NORMALISED FAR FIELD")
    print("  h'_00 / [2GM/(c^2 r)]")
    print("    = 1 + epsilon^2 [8/9 - (2/9) P2(cos Theta)")
    print("         - (2/15)(R/r)^2 P2(cos Theta) + O((R/r)^4)] + O(epsilon^4)")
    print("  epsilon = Omega R/c")
    print("  Full argument transform changes [8/9, -2/9] to [7/9, -1/9] at 1/r;")
    print("  it leaves the formal -2/15 coefficient at 1/r^3 unchanged.")
    print()
    print("KERR GATES")
    print("  Gate 1: Kerr permits no direction-dependent P2/r mass aspect in ACMC asymptotics.")
    print("          Component-only coefficient = -2/9; full-argument coefficient = -1/9.  KILL both ways.")
    print("  Gate 2 (formal only, because Gate 1 already failed):")
    print("          Q_model = -(2/15) M Omega^2 R^4/c^2")
    print("          J_shell = (2/3) M R^2 Omega")
    print("          Q_model / Q_Kerr = 3/10, with Q_Kerr = -J^2/(M c^2).  KILL.")
    print()
    print("VERDICT: KILLED_BEFORE_QUADRUPOLE")
    print("  THM-M's first-order clock+rods boost and linear superposition do not extend to Kerr at O(J^2).")
    print("  The failure is sharper than a wrong quadrupole: the rods generate an anisotropic 1/r term.")
    print("  This kills only the minimal second-order continuation, not the seated-root model.  A nonlinear")
    print("  source/field rule would be new declared or derived structure and must predict the cancellation")
    print("  before it may be tested against Kerr's invariant quadrupole.")

    return 0 if all(passed for _, passed in checks) else 1


if __name__ == "__main__":
    sys.exit(main())

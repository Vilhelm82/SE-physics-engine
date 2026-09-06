import sympy as sp

from thm_n_kerr_quadrupole import (
    derive_covariant_argument_coefficients,
    derive_local_boost_weights,
    derive_second_order_coefficients,
    evaluate_kerr_gate,
)


def test_second_order_coefficients_are_the_exact_shell_averages():
    """Catches a wrong angular average or a dropped rod-office term."""
    coeffs = derive_second_order_coefficients()

    assert coeffs["isotropic_r_inv"] == sp.Rational(8, 9)
    assert coeffs["p2_r_inv"] == -sp.Rational(2, 9)
    assert coeffs["p2_r_inv3"] == -sp.Rational(2, 15)


def test_transforming_the_field_argument_does_not_rescue_the_asymptotics():
    """Catches making the kill depend on a component-only boost convention."""
    coeffs = derive_covariant_argument_coefficients()

    assert coeffs["isotropic_r_inv"] == sp.Rational(7, 9)
    assert coeffs["p2_r_inv"] == -sp.Rational(1, 9)
    assert coeffs["p2_r_inv3"] == -sp.Rational(2, 15)


def test_covariant_argument_half_weight_follows_from_the_exact_boost():
    """Catches guessing the rest-frame distance correction or its sign."""
    weights = derive_local_boost_weights()

    assert weights == {
        "component_b_squared": 1,
        "component_b_dot_n_squared": 1,
        "covariant_b_squared": 1,
        "covariant_b_dot_n_squared": sp.Rational(1, 2),
    }


def test_minimal_extension_is_killed_before_the_kerr_quadrupole_gate():
    """Catches accidentally treating a direction-dependent 1/r term as Kerr-compatible."""
    verdict = evaluate_kerr_gate()

    assert verdict["asymptotic_monopole_passes"] is False
    assert verdict["covariant_argument_asymptotic_passes"] is False
    assert verdict["formal_quadrupole_passes"] is False
    assert verdict["formal_q_over_kerr_q"] == sp.Rational(3, 10)
    assert verdict["status"] == "KILLED_BEFORE_QUADRUPOLE"

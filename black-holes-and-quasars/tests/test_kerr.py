import pytest

from src.constants import M_SUN
from src.schwarzschild import schwarzschild_radius, isco_radius as schwarzschild_isco
from src.kerr import gravitational_radius, horizon_radius, isco_radius, photon_sphere_radius


def test_zero_spin_horizon_matches_schwarzschild():
    assert horizon_radius(M_SUN, 0.0) == pytest.approx(schwarzschild_radius(M_SUN))


def test_zero_spin_isco_matches_schwarzschild():
    assert isco_radius(M_SUN, 0.0) == pytest.approx(schwarzschild_isco(M_SUN))


def test_exactly_extremal_prograde_isco_equals_gravitational_radius():
    r_g = gravitational_radius(M_SUN)
    assert isco_radius(M_SUN, 1.0, prograde=True) == pytest.approx(r_g, rel=1e-6)


def test_near_extremal_prograde_isco_matches_bardeen_1972_table():
    """Bardeen, Press & Teukolsky 1972, ApJ 178, 347, table 1: r_isco = 1.237 M at a* = 0.998."""
    r_g = gravitational_radius(M_SUN)
    assert isco_radius(M_SUN, 0.998, prograde=True) == pytest.approx(1.237 * r_g, rel=1e-3)


def test_prograde_isco_smaller_than_retrograde():
    prograde = isco_radius(M_SUN, 0.9, prograde=True)
    retrograde = isco_radius(M_SUN, 0.9, prograde=False)
    assert prograde < retrograde


def test_prograde_photon_sphere_smaller_than_retrograde():
    prograde = photon_sphere_radius(M_SUN, 0.9, prograde=True)
    retrograde = photon_sphere_radius(M_SUN, 0.9, prograde=False)
    assert prograde < retrograde

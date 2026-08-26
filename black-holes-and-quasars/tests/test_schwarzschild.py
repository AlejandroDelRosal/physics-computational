import pytest

from src.constants import M_SUN
from src.schwarzschild import (
    schwarzschild_radius,
    photon_sphere_radius,
    isco_radius,
    gravitational_redshift,
)


def test_solar_schwarzschild_radius():
    assert schwarzschild_radius(M_SUN) == pytest.approx(2953, rel=1e-3)


def test_earth_schwarzschild_radius():
    m_earth = 5.972e24
    assert schwarzschild_radius(m_earth) == pytest.approx(8.87e-3, rel=1e-2)


def test_photon_sphere_is_1_5_schwarzschild_radii():
    r_s = schwarzschild_radius(M_SUN)
    assert photon_sphere_radius(M_SUN) == pytest.approx(1.5 * r_s)


def test_isco_is_3_schwarzschild_radii():
    r_s = schwarzschild_radius(M_SUN)
    assert isco_radius(M_SUN) == pytest.approx(3 * r_s)


def test_redshift_diverges_at_horizon():
    r_s = schwarzschild_radius(M_SUN)
    with pytest.raises(ValueError):
        gravitational_redshift(r_s, M_SUN)


def test_redshift_vanishes_far_away():
    r_far = schwarzschild_radius(M_SUN) * 1e12
    assert gravitational_redshift(r_far, M_SUN) == pytest.approx(0, abs=1e-6)

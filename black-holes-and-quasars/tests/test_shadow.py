import pytest

from src.constants import M_SUN, PC, MICROARCSEC
from src.shadow import shadow_angular_diameter


def test_m87_star_shadow_matches_eht_2019():
    """EHT Collaboration 2019, ApJL 875, L1: 42 +/- 3 microarcsec."""
    mass = 6.5e9 * M_SUN
    distance = 16.8e6 * PC
    predicted_uas = shadow_angular_diameter(mass, distance) / MICROARCSEC
    assert predicted_uas == pytest.approx(42, abs=6)


def test_sgr_a_star_shadow_matches_eht_2022():
    """EHT Collaboration 2022, ApJL 930, L12: 51.8 +/- 2.3 microarcsec."""
    mass = 4.15e6 * M_SUN
    distance = 8178 * PC
    predicted_uas = shadow_angular_diameter(mass, distance) / MICROARCSEC
    assert predicted_uas == pytest.approx(51.8, abs=5)

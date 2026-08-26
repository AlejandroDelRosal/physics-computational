import pytest

from src.constants import M_SUN
from src.accretion import eddington_luminosity, eddington_ratio, virial_mass_hbeta


def test_eddington_luminosity_one_solar_mass():
    assert eddington_luminosity(M_SUN) == pytest.approx(1.26e31, rel=0.02)


def test_eddington_ratio_at_eddington_limit():
    l_edd = eddington_luminosity(1e8 * M_SUN)
    assert eddington_ratio(l_edd, 1e8 * M_SUN) == pytest.approx(1.0)


def test_virial_mass_matches_shen2011_worked_example():
    mass = virial_mass_hbeta(fwhm_kms=1347.2, l5100_erg_s=10**44.435)
    assert mass / M_SUN == pytest.approx(10**7.39, rel=0.05)

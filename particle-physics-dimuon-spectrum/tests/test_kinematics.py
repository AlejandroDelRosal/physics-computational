import pytest

from src.kinematics import invariant_mass


def test_matches_published_mass_row_1():
    m = invariant_mass(8.38562, -1.76576, -0.10158, -8.19629, 7.24459, 2.02753, -0.2585, -6.94948)
    assert m == pytest.approx(3.83549, rel=1e-3)


def test_matches_published_mass_row_2():
    m = invariant_mass(5.90844, -0.578086, -4.61586, 3.64118, 5.93491, 1.91387, 4.5153, 3.34084)
    assert m == pytest.approx(9.47217, rel=1e-3)


def test_back_to_back_massless_pair_has_finite_mass():
    m = invariant_mass(10, 10, 0, 0, 10, -10, 0, 0)
    assert m == pytest.approx(20)

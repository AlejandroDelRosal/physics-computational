import pytest

from src.ginzburg_landau import build_grid, vortex_ansatz
from src.topology import winding_number


@pytest.mark.parametrize("winding", [-2, -1, 1, 2, 3])
def test_winding_number_matches_ansatz_charge(winding):
    x, y, xx, yy, _dx = build_grid(size=15, n_points=120)
    psi = vortex_ansatz(xx, yy, winding)
    assert winding_number(x, y, psi) == pytest.approx(winding, abs=1e-6)


def test_winding_number_is_additive_under_field_multiplication():
    """Fusing two vortices multiplies their order parameters and adds their
    charges: this is the group law of pi_1(U(1)) = Z acting on the field."""
    x, y, xx, yy, _dx = build_grid(size=15, n_points=120)
    psi_combined = vortex_ansatz(xx, yy, 1) * vortex_ansatz(xx, yy, 2)
    assert winding_number(x, y, psi_combined) == pytest.approx(3, abs=1e-6)


def test_uniform_field_has_zero_winding():
    x, y, xx, yy, _dx = build_grid(size=15, n_points=120)
    psi = 0 * xx + 1
    assert winding_number(x, y, psi) == pytest.approx(0, abs=1e-6)

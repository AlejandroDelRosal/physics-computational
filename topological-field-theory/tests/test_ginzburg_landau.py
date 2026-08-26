import numpy as np
import pytest

from src.ginzburg_landau import build_grid, vortex_ansatz, relax, free_energy
from src.topology import winding_number


def test_gradient_flow_decreases_free_energy():
    x, y, xx, yy, dx = build_grid(size=10, n_points=60)
    psi = vortex_ansatz(xx, yy, 1)
    energy_before = free_energy(psi, dx)
    psi = relax(psi, dx, steps=400)
    energy_after = free_energy(psi, dx)
    assert energy_after < energy_before


def test_relaxation_preserves_winding_number():
    x, y, xx, yy, dx = build_grid(size=10, n_points=60)
    psi = vortex_ansatz(xx, yy, 1)
    psi = relax(psi, dx, steps=400)
    assert winding_number(x, y, psi, radius=4) == pytest.approx(1, abs=0.05)


def test_field_amplitude_recovers_vacuum_far_from_core():
    x, y, xx, yy, dx = build_grid(size=10, n_points=60)
    psi = vortex_ansatz(xx, yy, 1)
    psi = relax(psi, dx, steps=400)
    edge_amplitude = np.abs(psi[30, -2])
    assert edge_amplitude == pytest.approx(1.0, abs=0.05)


def test_vortex_energy_scales_logarithmically_with_domain_size():
    """Kosterlitz & Thouless, J. Phys. C 6, 1181 (1973): E(R) ~ 2*pi*n^2*ln(R)."""
    sizes = [8, 14, 20]
    energies = []
    for size in sizes:
        x, y, xx, yy, dx = build_grid(size=size, n_points=70)
        psi = relax(vortex_ansatz(xx, yy, 1), dx, steps=1500)
        energies.append(free_energy(psi, dx))
    slope = np.polyfit(np.log(sizes), energies, 1)[0]
    assert slope == pytest.approx(2 * np.pi, rel=0.25)

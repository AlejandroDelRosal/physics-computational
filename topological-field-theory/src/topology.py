import numpy as np
from scipy.interpolate import RegularGridInterpolator


def winding_number(x, y, psi, center=(0.0, 0.0), radius: float = 5.0, n_samples: int = 400) -> float:
    """Topological charge: (1/2pi) times the total phase change around a loop.

    Classifies the field configuration as an element of the first homotopy
    group of the order parameter space, pi_1(U(1)) = Z (Mermin, Rev. Mod.
    Phys. 51, 591, 1979).
    """
    interpolate_real = RegularGridInterpolator((x, y), psi.real)
    interpolate_imag = RegularGridInterpolator((x, y), psi.imag)
    angles = np.linspace(0, 2 * np.pi, n_samples, endpoint=False)
    points = np.stack(
        [center[0] + radius * np.cos(angles), center[1] + radius * np.sin(angles)], axis=1
    )
    phase = np.arctan2(interpolate_imag(points), interpolate_real(points))
    delta_phase = np.diff(np.concatenate([phase, phase[:1]]))
    delta_phase = (delta_phase + np.pi) % (2 * np.pi) - np.pi
    return float(delta_phase.sum() / (2 * np.pi))

import numpy as np


def build_grid(size: float, n_points: int):
    x = np.linspace(-size, size, n_points)
    y = np.linspace(-size, size, n_points)
    dx = x[1] - x[0]
    xx, yy = np.meshgrid(x, y, indexing="ij")
    return x, y, xx, yy, dx


def vortex_ansatz(xx, yy, winding: int):
    r = np.sqrt(xx**2 + yy**2)
    theta = np.arctan2(yy, xx)
    return np.tanh(r) * np.exp(1j * winding * theta)


def multi_vortex_ansatz(xx, yy, cores: list[tuple[float, float, int]]):
    """cores: list of (x0, y0, winding). Charges combine additively: this is
    the group law of pi_1(U(1)) = Z applied to a multi-defect configuration."""
    psi = np.ones_like(xx, dtype=complex)
    for x0, y0, winding in cores:
        r = np.sqrt((xx - x0) ** 2 + (yy - y0) ** 2)
        theta = np.arctan2(yy - y0, xx - x0)
        psi = psi * np.tanh(r) * np.exp(1j * winding * theta)
    return psi


def laplacian(field, dx: float):
    lap = np.zeros_like(field)
    lap[1:-1, 1:-1] = (
        field[2:, 1:-1] + field[:-2, 1:-1] + field[1:-1, 2:] + field[1:-1, :-2] - 4 * field[1:-1, 1:-1]
    ) / dx**2
    return lap


def relax(psi0, dx: float, steps: int, dt: float | None = None):
    """Gradient flow on the Ginzburg-Landau free energy, dpsi/dt = -dF/dpsi*.

    Dirichlet boundary conditions pin the field to its initial value at
    the domain edge, which anchors the topological charge of the interior.
    """
    if dt is None:
        dt = 0.9 * dx**2 / 4
    psi = psi0.copy()
    boundary = psi0.copy()
    for _ in range(steps):
        psi = psi + dt * (laplacian(psi, dx) + psi * (1 - np.abs(psi) ** 2))
        psi[0, :], psi[-1, :] = boundary[0, :], boundary[-1, :]
        psi[:, 0], psi[:, -1] = boundary[:, 0], boundary[:, -1]
    return psi


def free_energy(psi, dx: float) -> float:
    dpsi_dx = np.zeros_like(psi)
    dpsi_dy = np.zeros_like(psi)
    dpsi_dx[1:-1, :] = (psi[2:, :] - psi[:-2, :]) / (2 * dx)
    dpsi_dy[:, 1:-1] = (psi[:, 2:] - psi[:, :-2]) / (2 * dx)
    gradient_term = np.abs(dpsi_dx) ** 2 + np.abs(dpsi_dy) ** 2
    potential_term = 0.5 * (1 - np.abs(psi) ** 2) ** 2
    interior = slice(1, -1)
    return float((gradient_term[interior, interior] + potential_term[interior, interior]).sum() * dx**2)

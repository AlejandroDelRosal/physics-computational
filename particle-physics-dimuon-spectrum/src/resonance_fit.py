import numpy as np
from scipy.optimize import curve_fit

PDG_Z_MASS_GEV = 91.1876
PDG_Z_WIDTH_GEV = 2.4952


def relativistic_breit_wigner(m, mass, width):
    m = np.asarray(m, dtype=float)
    return (mass * width) ** 2 / ((m**2 - mass**2) ** 2 + (mass * width) ** 2)


def voigt_resonance(m, mass, width, resolution, amplitude):
    """Relativistic Breit-Wigner convolved with a Gaussian detector resolution."""
    m = np.atleast_1d(np.asarray(m, dtype=float))
    half_range = 15 * max(width, resolution, 1.0)
    grid = np.linspace(mass - half_range, mass + half_range, 2000)
    bw = relativistic_breit_wigner(grid, mass, width)
    kernel = np.exp(-0.5 * ((m[:, None] - grid[None, :]) / resolution) ** 2)
    return amplitude * np.trapezoid(bw[None, :] * kernel, grid, axis=1)


def signal_plus_background(m, mass, width, resolution, amplitude, bg_amplitude, bg_slope, m_min):
    """Resonance plus an exponential Drell-Yan continuum under the peak."""
    m = np.asarray(m, dtype=float)
    return voigt_resonance(m, mass, width, resolution, amplitude) + bg_amplitude * np.exp(
        -bg_slope * (m - m_min)
    )


def fit_z_peak(masses, mass_range=(75, 107), bins=64):
    counts, edges = np.histogram(masses, bins=bins, range=mass_range)
    centers = (edges[:-1] + edges[1:]) / 2
    m_min = mass_range[0]

    def model(m, mass, width, resolution, amplitude, bg_amplitude, bg_slope):
        return signal_plus_background(m, mass, width, resolution, amplitude, bg_amplitude, bg_slope, m_min)

    p0 = [PDG_Z_MASS_GEV, PDG_Z_WIDTH_GEV, 2.0, counts.max(), counts.min() + 1, 0.05]
    bounds = ([80, 0.1, 0.1, 0, 0, 0], [100, 20, 20, np.inf, np.inf, 5])
    popt, pcov = curve_fit(model, centers, counts, p0=p0, bounds=bounds, sigma=np.sqrt(counts + 1))
    return popt, np.sqrt(np.diag(pcov)), centers, counts, m_min

import numpy as np

from .constants import G, C


def gravitational_radius(mass_kg: float) -> float:
    return G * mass_kg / C**2


def horizon_radius(mass_kg: float, spin: float) -> float:
    r_g = gravitational_radius(mass_kg)
    return r_g * (1 + (1 - spin**2) ** 0.5)


def isco_radius(mass_kg: float, spin: float, prograde: bool = True) -> float:
    """Bardeen, Press & Teukolsky 1972, ApJ 178, 347, eq. 2.21."""
    r_g = gravitational_radius(mass_kg)
    z1 = 1 + (1 - spin**2) ** (1 / 3) * (
        (1 + spin) ** (1 / 3) + (1 - spin) ** (1 / 3)
    )
    z2 = (3 * spin**2 + z1**2) ** 0.5
    sign = -1 if prograde else 1
    r_isco_over_rg = 3 + z2 + sign * ((3 - z1) * (3 + z1 + 2 * z2)) ** 0.5
    return r_isco_over_rg * r_g


def photon_sphere_radius(mass_kg: float, spin: float, prograde: bool = True) -> float:
    """Bardeen, Press & Teukolsky 1972, ApJ 178, 347, eq. 2.18."""
    r_g = gravitational_radius(mass_kg)
    sign = -1 if prograde else 1
    r_ph_over_rg = 2 * (1 + np.cos((2 / 3) * np.arccos(sign * spin)))
    return r_ph_over_rg * r_g

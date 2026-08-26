from .constants import G, C


def schwarzschild_radius(mass_kg: float) -> float:
    return 2 * G * mass_kg / C**2


def photon_sphere_radius(mass_kg: float) -> float:
    return 1.5 * schwarzschild_radius(mass_kg)


def isco_radius(mass_kg: float) -> float:
    return 3 * schwarzschild_radius(mass_kg)


def gravitational_redshift(r_m: float, mass_kg: float) -> float:
    r_s = schwarzschild_radius(mass_kg)
    if r_m <= r_s:
        raise ValueError("r must be outside the event horizon")
    return 1 / (1 - r_s / r_m) ** 0.5 - 1


def escape_velocity(r_m: float, mass_kg: float) -> float:
    return (2 * G * mass_kg / r_m) ** 0.5

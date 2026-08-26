from .kerr import gravitational_radius


def shadow_angular_diameter(mass_kg: float, distance_m: float) -> float:
    """Critical impact parameter b_crit = 3*sqrt(3)*r_g (EHT Collaboration 2019, ApJL 875, L1, eq. 2)."""
    r_g = gravitational_radius(mass_kg)
    b_crit = 3 * 3**0.5 * r_g
    return 2 * b_crit / distance_m

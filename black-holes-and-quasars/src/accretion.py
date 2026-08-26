import numpy as np

from .constants import G, C, M_PROTON, SIGMA_THOMSON, M_SUN


def eddington_luminosity(mass_kg: float) -> float:
    return 4 * np.pi * G * mass_kg * M_PROTON * C / SIGMA_THOMSON


def eddington_ratio(bolometric_luminosity_w: float, mass_kg: float) -> float:
    return bolometric_luminosity_w / eddington_luminosity(mass_kg)


def virial_mass_hbeta(fwhm_kms: float, l5100_erg_s: float) -> float:
    """Vestergaard & Peterson 2006, ApJ 641, 689, eq. 5."""
    log_m = 6.91 + 0.50 * np.log10(l5100_erg_s / 1e44) + 2 * np.log10(fwhm_kms / 1000)
    return 10**log_m * M_SUN

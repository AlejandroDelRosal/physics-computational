import numpy as np
import pytest

from src.data_loader import load_dimuon_events
from src.resonance_fit import fit_z_peak, relativistic_breit_wigner, PDG_Z_MASS_GEV, PDG_Z_WIDTH_GEV


def test_breit_wigner_peaks_at_resonance_mass():
    values = relativistic_breit_wigner(np.array([80, 91.1876, 100]), mass=91.1876, width=2.5)
    assert values[1] > values[0]
    assert values[1] > values[2]


def test_z_peak_fit_matches_pdg_on_real_data():
    """Particle Data Group world-average values: M_Z, Gamma_Z.

    The fit models the resonance as a relativistic Breit-Wigner
    convolved with a Gaussian detector resolution, plus an exponential
    Drell-Yan continuum background under the peak.
    """
    events = load_dimuon_events()
    (mass, width, _resolution, _amp, _bg_amp, _bg_slope), _errors, _c, _n, _m = fit_z_peak(events["M"])
    assert mass == pytest.approx(PDG_Z_MASS_GEV, abs=0.5)
    assert width == pytest.approx(PDG_Z_WIDTH_GEV, abs=0.5)

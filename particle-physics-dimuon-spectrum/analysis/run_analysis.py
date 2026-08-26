import pathlib

import numpy as np
import matplotlib.pyplot as plt

from src.data_loader import load_dimuon_events
from src.kinematics import invariant_mass
from src.resonance_fit import fit_z_peak, signal_plus_background, PDG_Z_MASS_GEV, PDG_Z_WIDTH_GEV

RESULTS_DIR = pathlib.Path(__file__).parent.parent / "results"

RESONANCES = [
    ("rho/omega", 0.78),
    ("phi", 1.02),
    ("J/psi", 3.10),
    ("psi(2S)", 3.69),
    ("Upsilon", 9.46),
    ("Z", 91.19),
]


def plot_full_spectrum(masses):
    fig, ax = plt.subplots(figsize=(9, 6))
    bins = np.logspace(np.log10(0.25), np.log10(300), 300)
    ax.hist(masses, bins=bins, histtype="stepfilled", color="#2b6cb0")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Dimuon invariant mass (GeV)")
    ax.set_ylabel("Events")
    ax.set_title(f"CMS Run2011A dimuon spectrum, {len(masses)} opposite-charge pairs")
    for name, mass in RESONANCES:
        ax.axvline(mass, color="gray", linewidth=0.5, linestyle="--")
        ax.text(mass, ax.get_ylim()[1] * 0.6, name, rotation=90, fontsize=8, ha="right", va="top")
    fig.savefig(RESULTS_DIR / "dimuon_spectrum.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_z_peak_fit(centers, counts, popt, m_min):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(centers, counts, yerr=np.sqrt(counts), fmt="o", markersize=3, color="#2b6cb0", label="Data")
    m_fine = np.linspace(centers.min(), centers.max(), 500)
    ax.plot(m_fine, signal_plus_background(m_fine, *popt, m_min), color="#c05621", label="Signal + Drell-Yan background")
    ax.set_xlabel("Dimuon invariant mass (GeV)")
    ax.set_ylabel("Events")
    ax.set_title("Z boson resonance")
    ax.legend()
    fig.savefig(RESULTS_DIR / "z_peak_fit.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    events = load_dimuon_events()

    computed_mass = invariant_mass(
        events["E1"], events["px1"], events["py1"], events["pz1"],
        events["E2"], events["px2"], events["py2"], events["pz2"],
    )
    residual = computed_mass - events["M"]
    print(f"Events: {len(events)}")
    print(f"Invariant mass reproduction: median residual {residual.median():.2e} GeV, max {residual.abs().max():.2e} GeV")

    plot_full_spectrum(events["M"])

    popt, errors, centers, counts, m_min = fit_z_peak(events["M"])
    mass, width = popt[0], popt[1]
    mass_err, width_err = errors[0], errors[1]
    plot_z_peak_fit(centers, counts, popt, m_min)

    print(f"\nZ boson resonance fit (Breit-Wigner convolved with detector resolution, plus Drell-Yan background)")
    print(f"  Fitted mass:  {mass:.3f} +/- {mass_err:.3f} GeV  (PDG: {PDG_Z_MASS_GEV} GeV)")
    print(f"  Fitted width: {width:.3f} +/- {width_err:.3f} GeV  (PDG: {PDG_Z_WIDTH_GEV} GeV)")


if __name__ == "__main__":
    main()

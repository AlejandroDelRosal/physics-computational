import pathlib

import numpy as np
import matplotlib.pyplot as plt

from src.constants import M_SUN, PC, MICROARCSEC
from src.data_loader import load_hbeta_subsample
from src.accretion import virial_mass_hbeta, eddington_luminosity, eddington_ratio
from src.schwarzschild import schwarzschild_radius
from src.shadow import shadow_angular_diameter

RESULTS_DIR = pathlib.Path(__file__).parent.parent / "results"


def compute_derived_quantities(df):
    df = df.copy()
    l5100_erg_s = 10 ** df["logL5100"]
    computed_mass = virial_mass_hbeta(df["W(BHb)"].to_numpy(), l5100_erg_s.to_numpy())
    df["log_mass_computed"] = np.log10(computed_mass / M_SUN)

    mass_kg = 10 ** df["logBH"] * M_SUN
    df["schwarzschild_radius_km"] = schwarzschild_radius(mass_kg) / 1000

    lbol_w = 10 ** df["logLbol"] * 1e-7
    df["eddington_ratio"] = eddington_ratio(lbol_w, mass_kg)
    return df


def plot_mass_validation(df):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(df["logBHHV"], df["log_mass_computed"], s=4, alpha=0.3)
    lims = [df["logBHHV"].min(), df["logBHHV"].max()]
    ax.plot(lims, lims, color="black", linewidth=1)
    ax.set_xlabel("log(M_BH / M_sun), Shen et al. 2011 catalog")
    ax.set_ylabel("log(M_BH / M_sun), reproduced from Vestergaard & Peterson 2006")
    ax.set_title("Virial black hole mass: reproduction vs published catalog")
    fig.savefig(RESULTS_DIR / "mass_validation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_mass_function(df):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(df["logBH"], bins=60, color="#2b6cb0")
    ax.set_xlabel("log(M_BH / M_sun)")
    ax.set_ylabel("Number of quasars")
    ax.set_title(f"Black hole mass distribution, {len(df)} SDSS DR7 quasars")
    fig.savefig(RESULTS_DIR / "mass_function.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_eddington_ratio(df):
    fig, ax = plt.subplots(figsize=(7, 5))
    log_ratio = np.log10(df["eddington_ratio"].clip(lower=1e-4))
    ax.hist(log_ratio, bins=60, color="#c05621")
    ax.axvline(0, color="black", linewidth=1, linestyle="--")
    ax.set_xlabel("log(L_bol / L_Eddington)")
    ax.set_ylabel("Number of quasars")
    ax.set_title("Eddington ratio distribution")
    fig.savefig(RESULTS_DIR / "eddington_ratio.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_schwarzschild_scaling(df):
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(df["logBH"], df["schwarzschild_radius_km"], s=4, alpha=0.3)
    ax.set_yscale("log")
    ax.set_xlabel("log(M_BH / M_sun)")
    ax.set_ylabel("Schwarzschild radius (km)")
    ax.set_title("Event horizon size across the SDSS DR7 quasar population")
    fig.savefig(RESULTS_DIR / "schwarzschild_scaling.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def report_eht_comparison():
    targets = {
        "M87*": dict(mass=6.5e9 * M_SUN, distance=16.8e6 * PC, measured_uas=42.0, ref="EHT Collaboration 2019, ApJL 875, L1"),
        "Sgr A*": dict(mass=4.15e6 * M_SUN, distance=8178 * PC, measured_uas=51.8, ref="EHT Collaboration 2022, ApJL 930, L12"),
    }
    print("\nBlack hole shadow diameter: theory vs Event Horizon Telescope")
    for name, t in targets.items():
        predicted_uas = shadow_angular_diameter(t["mass"], t["distance"]) / MICROARCSEC
        print(f"  {name}: predicted {predicted_uas:.1f} uas, measured {t['measured_uas']} uas ({t['ref']})")


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    df = load_hbeta_subsample()
    df = compute_derived_quantities(df)

    plot_mass_validation(df)
    plot_mass_function(df)
    plot_eddington_ratio(df)
    plot_schwarzschild_scaling(df)

    residual = df["log_mass_computed"] - df["logBHHV"]
    print(f"Sample size (H-beta, non-BAL): {len(df)}")
    print(f"Mass reproduction: median offset {residual.median():.3f} dex, scatter {residual.std():.3f} dex")
    print(f"Median black hole mass: {10**df['logBH'].median():.2e} M_sun")
    print(f"Median Eddington ratio: {df['eddington_ratio'].median():.3f}")

    report_eht_comparison()


if __name__ == "__main__":
    main()

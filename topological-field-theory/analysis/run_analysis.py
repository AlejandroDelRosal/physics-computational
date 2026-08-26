import pathlib

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb

from src.ginzburg_landau import build_grid, vortex_ansatz, multi_vortex_ansatz, relax, free_energy
from src.topology import winding_number

RESULTS_DIR = pathlib.Path(__file__).parent.parent / "results"


def domain_color(psi):
    magnitude = np.abs(psi)
    hue = (np.angle(psi) + np.pi) / (2 * np.pi)
    value = magnitude / magnitude.max()
    hsv = np.stack([hue, np.ones_like(hue), value], axis=-1)
    return hsv_to_rgb(hsv)


def plot_single_vortex(x, y, psi):
    xx, yy = np.meshgrid(x, y, indexing="ij")

    fig = plt.figure(figsize=(14, 6))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax1.plot_surface(xx, yy, np.abs(psi) ** 2, cmap="viridis", linewidth=0, antialiased=True)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_zlabel("|psi|^2")
    ax1.set_title("Order parameter density around an n=1 vortex")

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.imshow(domain_color(psi).transpose(1, 0, 2), origin="lower", extent=[x.min(), x.max(), y.min(), y.max()])
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_title("Domain coloring: hue = phase, brightness = |psi|")

    fig.savefig(RESULTS_DIR / "single_vortex.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_vortex_antivortex_pair(x, y, psi):
    xx, yy = np.meshgrid(x, y, indexing="ij")

    fig = plt.figure(figsize=(14, 6))
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")
    ax1.plot_surface(xx, yy, np.abs(psi) ** 2, cmap="magma", linewidth=0, antialiased=True)
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_zlabel("|psi|^2")
    ax1.set_title("Vortex-antivortex pair: order parameter density")

    ax2 = fig.add_subplot(1, 2, 2)
    ax2.imshow(domain_color(psi).transpose(1, 0, 2), origin="lower", extent=[x.min(), x.max(), y.min(), y.max()])
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_title("Domain coloring: opposite charges, opposite hue rotation")

    fig.savefig(RESULTS_DIR / "vortex_antivortex_pair.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_energy_scaling(sizes, energies, slope, intercept):
    fig, ax = plt.subplots(figsize=(7, 5))
    log_sizes = np.log(sizes)
    ax.plot(log_sizes, energies, "o", color="#2b6cb0", label="Relaxed vortex energy")
    fit_line = slope * log_sizes + intercept
    ax.plot(log_sizes, fit_line, color="#c05621", label=f"Fit: slope = {slope:.3f} (theory: 2*pi = {2*np.pi:.3f})")
    ax.set_xlabel("ln(domain size L)")
    ax.set_ylabel("Free energy")
    ax.set_title("Isolated vortex energy diverges logarithmically with system size")
    ax.legend()
    fig.savefig(RESULTS_DIR / "energy_scaling.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    RESULTS_DIR.mkdir(exist_ok=True)

    x, y, xx, yy, dx = build_grid(size=15, n_points=150)
    psi_single = relax(vortex_ansatz(xx, yy, 1), dx, steps=20000)
    plot_single_vortex(x, y, psi_single)
    n_single = winding_number(x, y, psi_single, radius=8)
    print(f"Single vortex: winding number = {n_single:.4f}")

    x2, y2, xx2, yy2, dx2 = build_grid(size=15, n_points=150)
    # A free vortex-antivortex pair attracts and annihilates under gradient
    # flow (real GL dynamics); a moderate step count relaxes the local core
    # profiles while keeping the pair separated for this snapshot.
    psi_pair = relax(multi_vortex_ansatz(xx2, yy2, [(5, 0, 1), (-5, 0, -1)]), dx2, steps=1000)
    plot_vortex_antivortex_pair(x2, y2, psi_pair)
    n_right = winding_number(x2, y2, psi_pair, center=(5, 0), radius=3)
    n_left = winding_number(x2, y2, psi_pair, center=(-5, 0), radius=3)
    n_total = winding_number(x2, y2, psi_pair, center=(0, 0), radius=13)
    print(f"\nVortex-antivortex pair: n_right(+1 core) = {n_right:.4f}, n_left(-1 core) = {n_left:.4f}, n_total (large loop) = {n_total:.4f}")

    sizes = [8, 12, 16, 20, 26, 32]
    energies = []
    for size in sizes:
        _x, _y, _xx, _yy, _dx = build_grid(size=size, n_points=140)
        psi = relax(vortex_ansatz(_xx, _yy, 1), _dx, steps=15000)
        energies.append(free_energy(psi, _dx))
    slope, intercept = np.polyfit(np.log(sizes), energies, 1)
    plot_energy_scaling(sizes, energies, slope, intercept)
    print(f"\nVortex energy scaling: fitted slope = {slope:.3f}, theory (Kosterlitz-Thouless 1973) = 2*pi*n^2 = {2*np.pi:.3f}")


if __name__ == "__main__":
    main()

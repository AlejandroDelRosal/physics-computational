import numpy as np
import matplotlib.pyplot as plt
import statistics

# Parameters
T_range = np.arange(1, 5, 0.25)
T_selected = [1.5, 2, 2.25, 4]
lattice_size = 10

# Generate initial spin array (all spins up)
def generate_spin_array(size):
    return np.ones((size, size), dtype=float)

# Monte Carlo using the Metropolis algorithm
def montecarlo_metropolis(spin_array, T, iterations=1000, kB=1):
    size = spin_array.shape[0]
    magnetization_list = []

    for _ in range(iterations):
        i, j = np.random.randint(0, size, 2)

        S_ij = spin_array[i, j]
        neighbors = (spin_array[(i+1) % size, j] +
                     spin_array[(i-1) % size, j] +
                     spin_array[i, (j+1) % size] +
                     spin_array[i, (j-1) % size])

        delta_E = 2 * S_ij * neighbors

        if delta_E < 0 or np.random.rand() < np.exp(-delta_E / (kB * T)):
            spin_array[i, j] *= -1  # Flip the spin

        magnetization_list.append(np.sum(spin_array))

    return magnetization_list, spin_array

# Main function
def main(T_range, T_selected, lattice_size, iterations=1000, kB=1):
    avg_magnetization = []
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))
    ax = np.ravel(ax)
    counter = 0
    time = range(iterations)

    for T in T_range:
        spin_array = generate_spin_array(lattice_size)
        magnetization, spin_array = montecarlo_metropolis(spin_array, T, iterations, kB)

        avg_magnetization.append(statistics.mean(np.array(magnetization) / (lattice_size * lattice_size)))

        # Plot time evolution of magnetization at selected temperatures
        if T in T_selected:
            ax[counter].plot(time, np.array(magnetization) / (lattice_size * lattice_size))
            ax[counter].set_ylabel(f"Magnetization at T = {T}")
            ax[counter].set_xlabel("Time")
            ax[counter].set_ylim(-1.1, 1.1)
            counter += 1

    plt.tight_layout()
    plt.show()

    # Plot average magnetization vs temperature
    plt.figure(figsize=(6, 4))
    plt.scatter(T_range, avg_magnetization, color="gray")
    plt.xlabel("Temperature (T)")
    plt.ylabel("Average Magnetization")
    plt.title("Magnetization vs Temperature")
    plt.show()

    return 0

# Run the simulation
main(T_range, T_selected, lattice_size, iterations=10000, kB=1)

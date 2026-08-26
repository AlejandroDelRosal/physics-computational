# Matrix method for the hydrogen atom

import numpy as np
import matplotlib.pyplot as plt

# Grid and system setup
dr = 0.1
R = 40
N = int((R - dr) / dr + 1)

# Create radial grid vector r
r = np.linspace(dr, R, N)  # N points from dr to R

# Build the second derivative matrix M (finite difference approximation)
m = np.zeros((N, N))
np.fill_diagonal(m, -2)
for i in range(1, N - 1):
    m[i, i - 1] = 1
    m[i, i + 1] = 1
m[0, 1] = 1
m[-1, -2] = 1
M = (-1 / (2 * dr**2)) * m

# Build the potential matrix V for the Coulomb potential (-1/r)
V = np.diag(-1 / r)

# Hamiltonian matrix H = T + V
H = M + V

# Compute eigenvalues and eigenvectors of the Hamiltonian
eigenvalues, eigenvectors = np.linalg.eigh(H)

# Sort eigenvalues and corresponding eigenvectors
sorted_indices = np.argsort(eigenvalues)
eigenvalues = eigenvalues[sorted_indices]
eigenvectors = eigenvectors[:, sorted_indices]

# Plot the square of the first three eigenfunctions
plt.figure(figsize=(10, 6))
for i in range(3):
    plt.plot(r, eigenvectors[:, i]**2, label=f'n={i+1}, E={eigenvalues[i]:.3f}')
plt.title('Hydrogen Atom: Radial Probability Densities')
plt.xlabel('r (atomic units)')
plt.ylabel('g(r)^2')
plt.legend()
plt.grid(True)
plt.show()

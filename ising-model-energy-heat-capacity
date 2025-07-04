import numpy as np
import matplotlib.pyplot as plt
from random import randrange, random

#constants
kB = 1; J = 1

T_i = 1; T_f = 5; dT = 0.25
Temperatures = np.arange(T_i, T_f + dT, dT)
# Size of the spin array (square NxN)
N = 10
steps = int(1e4)  

def energy(spins):
    """
    Calculate the energy of the spin system using periodic boundary conditions.
    """
    energy = 0.0
    for i in range(N):
        for j in range(N):
            S = spins[i, j]
            # Neighbors with periodic boundary conditions
            nb = spins[(i + 1) % N, j] + spins[i, (j + 1) % N]
            energy += -J * S * nb
    return energy

E_mean = []; Esqrd_mean = []; C = []
 # Initialization with all spins aligned
for n in range(len(Temperatures)):
    spins_0 = np.ones((N, N), dtype=int) 

    total_M = []
    E = []
    Esqrd = []

    T = Temperatures[n]
    beta = 1 / (kB * T)

    for k in range(steps):
        i = randrange(N)
        j = randrange(N)

        spins_new = spins_0.copy()
        spins_new[i, j] *= -1

        E_new = energy(spins_new)
        E_0 = energy(spins_0)
        dE = E_new - E_0
        #we accept or deny the energy shift
        if random() < np.exp(-beta * dE):
            spins_0 = spins_new
            E.append(E_new)
            Esqrd.append(E_new**2)
        else:
            E.append(E_0)
            Esqrd.append(E_0**2)
        total_M.append(np.sum(spins_0) / (N**2))

    E_mean_i = np.mean(E)
    Esqrd_mean_i = np.mean(Esqrd)
    sigma_Esqrd = Esqrd_mean_i - E_mean_i**2
    E_mean.append(E_mean_i)
    Esqrd_mean.append(Esqrd_mean_i)
    C.append(sigma_Esqrd / (kB * T**2))

E_mean = np.array(E_mean)
C = np.array(C)

# Graph of average energy per spin vs. temperature
plt.figure()
plt.scatter(Temperatures, E_mean / (N**2), label='Average energy per spin')
plt.title('Average energy per spin vs. temperature')
plt.xlabel('Temperature')
plt.ylabel('Average energy per spin')
plt.grid(True)
plt.legend()
plt.show()

# Plot of spin heat capacity vs. temperature
plt.figure()
plt.scatter(Temperatures, C / (N**2), label='Heat capacity per spin')
plt.title('Heat capacity per spin vs. temperature')
plt.xlabel('Tempeture')
plt.ylabel('Average heat capacity per spin')
plt.grid(True)
plt.legend()
plt.show()

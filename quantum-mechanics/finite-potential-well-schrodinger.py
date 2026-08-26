import numpy as np
import matplotlib.pyplot as plt

# Variables iniciales
L = 1
deltax = 0.001
N = int(L / deltax)  # Número de iteraciones dentro de la caja (x=0 a x=1)
E = 2  
deltaE = 0.1
limit = 2
deltaETol = 0.001
V_out = 50

# Define potential según el índice i del bucle
def V(i):
    return 0 if i < N else V_out  

# Bucle para ajustar E 
psi = [1, 1]  # Condiciones iniciales
i = 0
current_V = 0
# Primera iteración para determinar div_direction
while abs(psi[-1]) < limit:
    if i == N:
        current_V = V_out
    psi_next = 2 * psi[-1] - psi[-2] - 2 * (E - current_V) * psi[-1] * deltax**2
    psi.append(psi_next)
    i += 1
div_direction = np.sign(psi[-1])
if div_direction < 0:
    deltaE = -deltaE
E += deltaE

# Bucle principal de ajuste de E
while abs(deltaE) > deltaETol:
    psi = [1, 1]
    i = 0
    current_V = 0
    while abs(psi[-1]) < limit:
        if i == N:
            current_V = V_out
        psi_next = 2 * psi[-1] - psi[-2] - 2 * (E - current_V) * psi[-1] * deltax**2
        psi.append(psi_next)
        i += 1
    
    if np.sign(psi[-1]) * div_direction < 0:
        deltaE = -deltaE / 2
        div_direction *= -1
    E += deltaE

# Generar la gráfica
psi_mirror = np.flip(psi)
psi_full = np.concatenate((psi_mirror, psi))
L_total = deltax * len(psi)
x = np.linspace(-L_total, L_total, len(psi_full))

print(f'Valor de E: {E}')
plt.plot(x, psi_full)
plt.xlabel('x')
plt.ylabel('psi(x)')
plt.title('Solución de Pozo Finito para la Ecuación de Schrödinger')
plt.grid(True)
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Variables de la función de Lennard-Jones
epsilon = 10.0
sigma = 1.0

x = np.linspace(0.3, 3.5, 500)

# Función potencial de Lennard-Jones
def V(x):
    return 4 * epsilon * ((sigma / x)**12 - (sigma / x)**6)

# Parámetros iniciales
MP = 1.2
xl0 = 0.5
xr0 = 5.0
E = -1.5
deltaE = 0.1
deltax = 0.01

# Condiciones iniciales para las funciones de onda
psi_L = [-0.0001*deltax, 0.0]
psi_R = [-0.0001*deltax, 0.0]

# Funciones de integración
def integrate_left(x_start, psi_start, E):
    x = [x_start, x_start + deltax]
    psi = psi_start.copy()
    while x[-1] < MP:
        x_current = x[-1]
        psi_current = psi[-1]
        psi_double_prime = 2 * (V(x_current) - E) * psi_current  # Corregido: removido MP redundante
        psi_next = 2 * psi_current - psi[-2] + psi_double_prime * deltax**2
        psi.append(psi_next)
        x.append(x_current + deltax)
    return x, psi

def integrate_right(x_start, psi_start, E):
    x = [x_start, x_start - deltax]
    psi = psi_start.copy()
    while x[-1] > MP:
        x_current = x[-1]
        psi_current = psi[-1]
        psi_double_prime = 2 * (V(x_current) - E) * psi_current  # Corregido: removido MP redundante
        psi_next = 2 * psi_current - psi[-2] + psi_double_prime * deltax**2
        psi.append(psi_next)
        x.append(x_current - deltax)
    return x, psi

# Ciclo de búsqueda de energía
limit = 0.001
last = 0
condicion = False

while not condicion:
    # Integrar desde ambos lados
    x_left, psi_left = integrate_left(xl0, psi_L, E)
    x_right, psi_right = integrate_right(xr0, psi_R, E)
    
    # Re-escalar función derecha
    scaling_factor = psi_left[-1] / psi_right[-2]
    psi_right_scaled = [val * scaling_factor for val in psi_right]
    
    # Calcular derivadas
    dpsi_left = (psi_left[-1] - psi_left[-2]) / deltax
    dpsi_right = (psi_right_scaled[-1] - psi_right_scaled[-2]) / (-deltax)
    
    # Verificar convergencia
    if abs(dpsi_left - dpsi_right) > limit:
        now = 1 if (dpsi_left > dpsi_right) else -1
        
        if now != last:
            deltaE = -deltaE / 2
            last = now
        else:
            E += deltaE
            
        # Actualizar límite de convergencia
        limit = 0.001 * 0.5 * abs(dpsi_left + dpsi_right)
    else:
        condicion = True

# Resultados finales
plt.figure(figsize=(10, 6))
plt.plot(x_left, psi_left, label='Izquierda')
plt.plot(x_right, psi_right_scaled, label='Derecha (escalada)')
#plt.plot(x, V(x), 'k--', label='Potencial L-J')
plt.axvline(MP, color='gray', linestyle=':', label='MP')
plt.title(f'Solución para E = {E:.4f}')
plt.xlabel('Distancia (x)')
plt.ylabel('psi(x)')
plt.legend()
plt.grid(True)
plt.xlim(0.3, 3.5)
plt.show()

print(f'Energía encontrada: {E:.4f}')

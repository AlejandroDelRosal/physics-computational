
import numpy as np
import matplotlib.pyplot as plt

# Parámetros
omega_squared = 1  # omega^2 = k/m
m = 1              # masa
delta_t = 0.01     # paso de tiempo
y0 = 0.2           # posición inicial
v0 = 0             # velocidad inicial

# Tiempo total y número de pasos
t_total = 10
time_steps = int(t_total / delta_t)

# Inicialización de arrays
time = np.linspace(0, t_total, time_steps)
y = np.zeros(time_steps)
v = np.zeros(time_steps)
E = np.zeros(time_steps)

# Condiciones iniciales
y[0] = y0
v[0] = v0

# Cálculo de y usando el método de Euler y la fórmula iterativa
for i in range(1, time_steps):
    if i == 1:
        y[i] = y[0] + v[0] * delta_t
    else:
        y[i] = 2 * y[i-1] - y[i-2] - omega_squared * y[i-1] * delta_t**2
    v[i] = (y[i] - y[i-2]) / (2 * delta_t) if i > 1 else v[0]

    # Energía mecánica total (cinética + potencial)
    E[i] = 0.5 * m * v[i]**2 + 0.5 * m * omega_squared * y[i]**2

# Graficar resultados
plt.figure(figsize=(12, 6))

# Posición vs tiempo
plt.subplot(2, 1, 1)
plt.plot(time, y, label='Posición (y)')
plt.xlabel('Tiempo (s)')
plt.ylabel('Posición (m)')
plt.title('Oscilador armónico: Posición vs Tiempo')
plt.legend()

# Energía vs tiempo
plt.subplot(2, 1, 2)
plt.plot(time, E, label='Energía mecánica total (E)', color='orange')
plt.xlabel('Tiempo (s)')
plt.ylabel('Energía (J)')
plt.title('Energía mecánica total vs Tiempo')
plt.legend()

plt.tight_layout()
plt.show()

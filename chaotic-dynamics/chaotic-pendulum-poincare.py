# Using FD = 1.2 to observe chaotic behavior

import numpy as np
import matplotlib.pyplot as plt

# System parameters
L = g = 9.81       # Pendulum length and gravitational acceleration
q = 0.5            # Damping factor
FD = 1.2           # Amplitude of the external force
Omega_D = 2 / 3    # Driving force frequency
dt = 0.04          # Time step
t_max = 4000       # Total simulation time

# Initial conditions
theta1 = 0.2       # Initial angle (radians)
omega1 = 0.0       # Initial angular velocity (rad/s)

# Driving period
T_drive = 2 * np.pi / Omega_D

# Start recording after 100 drive periods
N_transients = 100
t_record_start = N_transients * T_drive

# Euler-Cromer integration step
def pendulum_step(theta, omega, t, dt):
    domega_dt = -g / L * np.sin(theta) - q * omega + FD * np.sin(Omega_D * t)
    omega_next = omega + domega_dt * dt
    theta_next = theta + omega_next * dt
    theta_next = (theta_next + np.pi) % (2 * np.pi) - np.pi  # Normalize angle
    return theta_next, omega_next

# Initialize variables
t = 0.0
theta = theta1
omega = omega1
poincare_theta = []
poincare_omega = []

# Simulation loop
while t < t_max:
    theta, omega = pendulum_step(theta, omega, t, dt)
    
    # Poincaré section: record only after transients and at multiples of drive period
    if t >= t_record_start:
        if abs(t % T_drive) < dt / 2:
            poincare_theta.append(theta)
            poincare_omega.append(omega)
    
    t += dt

# Plot Poincaré section
plt.figure(figsize=(8, 6))
plt.scatter(poincare_theta, poincare_omega, s=1, color='blue')
plt.title('Poincaré Section for FD = 1.2')
plt.xlabel(r'$\theta$ (radians)')
plt.ylabel(r'$\omega$ (rad/s)')
plt.grid(True)
plt.show()

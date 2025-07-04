
import numpy as np
import matplotlib.pyplot as plt

# Parameters of the system
g = L = 9.81       # Acceleration due to gravity in m/s^2 and Length of the pendulum in meters
q = 0.5        # Damping factor
FD1 = 1.2
FD2 = 0.5      # Amplitudes of the external force
Omega_D = 2.0 / 3.0  # Frequency of the external force
dt = 0.04      # Time step
t_max = 4000   # Total simulation time

# Initial conditions for FD=1.2
theta1 = 0.2   # Initial angle in radians
omega1 = 0.0   # Initial angular velocity

# Initial conditions for FD=0.5
theta2 = 0.201 # Slightly different initial angle
omega2 = 0.0   # Initial angular velocity

# Time array
t_values = np.arange(0, t_max, dt)

# Function for pendulum step with FD=1.2
def pendulum_step(theta, omega, t, dt, FD):
    dtheta_dt = omega
    domega_dt = -(g / L) * np.sin(theta) - q * omega + FD * np.sin(Omega_D * t)
    theta_next = theta + dtheta_dt * dt
    omega_next = omega + domega_dt * dt
    # Keep theta within [-pi, pi]
    theta_next = (theta_next + np.pi) % (2 * np.pi) - np.pi
    return theta_next, omega_next

# Simulation for FD=1.2
poincare_theta1 = []
poincare_omega1 = []
theta = theta1
omega = omega1

for t in t_values:
    # Poincaré section condition
    if abs(t % (2 * np.pi / Omega_D)) < dt / 2:
        poincare_theta1.append(theta)
        poincare_omega1.append(omega)
    # Advance the simulation
    theta, omega = pendulum_step(theta, omega, t, dt, FD1)

# Simulation for FD=0.5
poincare_theta2 = []
poincare_omega2 = []
theta = theta2
omega = omega2

for t in t_values:
    # Poincaré section condition
    if abs(t % (2 * np.pi / Omega_D)) < dt / 2:
        poincare_theta2.append(theta)
        poincare_omega2.append(omega)
    # Advance the simulation
    theta, omega = pendulum_step(theta, omega, t, dt, FD2)

# Plotting Poincaré sections
plt.figure(figsize=(12, 6))

# FD=1.2
plt.subplot(1, 2, 1)
plt.scatter(poincare_theta1, poincare_omega1, s=1, color='blue')
plt.title('Poincaré Section for FD = 1.2')
plt.xlabel(r'$\theta$ (radians)')
plt.ylabel(r'$\omega$ (rad/s)')
plt.grid(True)

# FD=0.5
plt.subplot(1, 2, 2)
plt.scatter(poincare_theta2, poincare_omega2, s=1, color='red')
plt.title('Poincaré Section for FD = 0.5')
plt.xlabel(r'$\theta$ (radians)')
plt.ylabel(r'$\omega$ (rad/s)')
plt.grid(True)

plt.tight_layout()
plt.show()

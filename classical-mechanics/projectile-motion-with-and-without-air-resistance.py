#Manuel Alejandro Del Rosal Saucedo
import numpy as np
import matplotlib.pyplot as plt

def simulate_projectile(v0, theta_deg, B_m=0, dt=0.01):
    # Constants
    g = 9.81  # gravitational acceleration (m/s^2)
    theta = np.radians(theta_deg)  # convert angle to radians
    
    # Initial conditions
    x0, y0 = 0, 0
    vx0 = v0 * np.cos(theta)
    vy0 = v0 * np.sin(theta)
    
    # Lists to store trajectory
    x_points = [x0]
    y_points = [y0]
    vx = vx0
    vy = vy0
    x = x0
    y = y0
    
    # Euler method iteration
    while y >= 0:
        # Calculate velocity magnitude for air resistance
        v = np.sqrt(vx**2 + vy**2)
        
        # Calculate accelerations
        ax = -B_m * v * vx if B_m > 0 else 0
        ay = -g - B_m * v * vy if B_m > 0 else -g
        
        # Update positions and velocities
        x = x + vx * dt
        y = y + vy * dt
        vx = vx + ax * dt
        vy = vy + ay * dt
        
        # Store positions
        if y >= 0:
            x_points.append(x)
            y_points.append(y)
    
    return np.array(x_points), np.array(y_points)

# Simulation parameters
v0 = 700  # initial velocity (m/s)
theta = 45  # angle (degrees)
B_m = 4e-5  # air resistance coefficient (m^-1)
dt = 0.01  # time step (s)

# Simulate both cases
x_no_resist, y_no_resist = simulate_projectile(v0, theta, B_m=0, dt=dt)
x_resist, y_resist = simulate_projectile(v0, theta, B_m=B_m, dt=dt)

# Create plot
plt.figure(figsize=(12, 6))
plt.plot(x_no_resist/1000, y_no_resist/1000, 'b-', label='Without air resistance')
plt.plot(x_resist/1000, y_resist/1000, 'r--', label='With air resistance')
plt.grid(True)
plt.xlabel('Distance (km)')
plt.ylabel('Height (km)')
plt.title(f'Projectile Motion (v₀={v0} m/s, θ={theta}°)')
plt.legend()

# Calculate and display impact points
print(f"Without air resistance:")
print(f"Maximum height: {np.max(y_no_resist)/1000:.2f} km")
print(f"Range: {x_no_resist[-1]/1000:.2f} km")
print(f"\nWith air resistance (B/m = {B_m}):")
print(f"Maximum height: {np.max(y_resist)/1000:.2f} km")
print(f"Range: {x_resist[-1]/1000:.2f} km")

plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Simulated Annealing parameters
T_max = 10
T_min = 0.001
tau = 10000
N = 25  # Number of cities
iterations_per_plot = 200  # How often to update the plot

# Generate random cities on a 10x10 grid
x = np.linspace(0, 10, 10)
y = np.linspace(0, 10, 10)
cities_x = np.random.choice(x, N)
cities_y = np.random.choice(y, N)

# Function to compute total route distance
def compute_distance(route):
    distance = 0.0
    for i in range(len(route)):
        current = route[i]
        next_city = route[(i + 1) % N]
        dx = cities_x[current] - cities_x[next_city]
        dy = cities_y[current] - cities_y[next_city]
        distance += np.sqrt(dx**2 + dy**2)
    return distance

# Initial configuration
current_route = list(range(N))
current_distance = compute_distance(current_route)
best_route = current_route.copy()
best_distance = current_distance

# Initial route plot
plt.figure(figsize=(10, 6))
plt.scatter(cities_x, cities_y, color='purple', marker='o', s=100, label='Cities')

# Connect cities in order
for i in range(N - 1):
    plt.plot([cities_x[current_route[i]], cities_x[current_route[i + 1]]],
             [cities_y[current_route[i]], cities_y[current_route[i + 1]]],
             'b-', alpha=0.7)

# Close the loop
plt.plot([cities_x[current_route[-1]], cities_x[current_route[0]]],
         [cities_y[current_route[-1]], cities_y[current_route[0]]],
         'b-', alpha=0.7)

# Label cities
for i in range(N):
    plt.text(cities_x[i] + 0.1, cities_y[i] + 0.1, str(i), fontsize=8)

plt.title(f'Initial Map with Route for {N} Cities')
plt.grid(True)
plt.legend()
plt.show()

# Optimization loop
T = T_max
t = 0
distance_history = []

while T > T_min:
    # Swap two cities
    i, j = np.random.choice(N, 2, replace=False)
    new_route = current_route.copy()
    new_route[i], new_route[j] = new_route[j], new_route[i]
    
    new_distance = compute_distance(new_route)
    delta = new_distance - current_distance
    
    # Acceptance criterion
    if delta < 0 or np.random.rand() < np.exp(-delta / T):
        current_route = new_route
        current_distance = new_distance
        
        if new_distance < best_distance:
            best_route = new_route.copy()
            best_distance = new_distance
    
    # Update temperature
    t += 1
    T = T_max * np.exp(-t / tau)
    distance_history.append(current_distance)
    
    # Show progress
    if t % iterations_per_plot == 0:
        plt.figure(figsize=(10, 6))
        plt.scatter(cities_x, cities_y, color='purple', marker='o', s=100)
        for k in range(N):
            start = current_route[k]
            end = current_route[(k + 1) % N]
            plt.plot([cities_x[start], cities_x[end]],
                     [cities_y[start], cities_y[end]], 'b-', alpha=0.5)
        for i in range(N):
            plt.text(cities_x[i] + 0.1, cities_y[i] + 0.1, str(i), fontsize=8)
        plt.title(f'Iteration: {t} - Temp: {T:.2f}\nCurrent Distance: {current_distance:.2f}')
        plt.grid(True)
        plt.show()

# Final results
plt.figure(figsize=(14, 6))

# Best route plot
plt.subplot(1, 2, 1)
plt.scatter(cities_x, cities_y, color='purple', marker='o', s=100)
for i in range(N):
    plt.text(cities_x[i] + 0.1, cities_y[i] + 0.1, str(i), fontsize=8)
for k in range(N):
    start = best_route[k]
    end = best_route[(k + 1) % N]
    plt.plot([cities_x[start], cities_x[end]],
             [cities_y[start], cities_y[end]], 'r-', linewidth=2)
plt.title(f'Best Route Found\nDistance: {best_distance:.2f}')
plt.grid(True)

# Convergence plot
plt.subplot(1, 2, 2)
plt.plot(distance_history)
plt.title('Algorithm Convergence')
plt.xlabel('Iterations')
plt.ylabel('Distance')
plt.grid(True)

plt.tight_layout()
plt.show()

# Print final route
print(f"\nBest distance found: {best_distance:.2f}")
print("Order of cities in the best route:")
print(' -> '.join(map(str, best_route)) + f' -> {best_route[0]}")

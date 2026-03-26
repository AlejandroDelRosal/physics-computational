import numpy as np
import matplotlib.pyplot as plt
from scipy.special import genlaguerre
from scipy.integrate import simpson

dr = 0.01
r = np.arange(0, 40 + dr, dr)
mu = 1  
hbar = 1  

for n in (1, 2, 3):
    l = 0  
    k = 2 * l + 1  # Parámetro de los polinomios de Laguerre
    
    # Energía del estado n (unidades atómicas)
    E = -1 / (2 * n**2)
    
    # Parámetro de escala
    epsilon = np.sqrt((-8 * mu * E) / hbar**2)
    x = r * epsilon  # Variable adimensional
    
    # Polinomio de Laguerre generalizado L_{n-l-1}^{2l+1}(x)
    L = genlaguerre(n - l - 1, k)
    g = np.exp(-x / 2) * x**((k + 1) / 2) * L(x)
    
    # Normalización (incluyendo Jacobiano x^2 / epsilon^3)
    g_squared = g**2
    g_squared_jacobian = g_squared * x**2 / epsilon**3
    norm = np.sqrt(simpson(g_squared_jacobian, x))
    g_normalized = g / norm
    
    plt.plot(r, g_normalized, label=f"n = {n}, l = {l}")

    

plt.xlabel("r (unidades atómicas")#Unidades que milagros nos dio igual a 1
plt.ylabel("g(r)")
plt.title("Funciones de onda radiales normalizadas ")
plt.legend()
plt.grid(True)
plt.show()

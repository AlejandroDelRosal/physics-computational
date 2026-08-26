import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad



#Potential Barrier 
def V(x, V1, V2, w, L):
    if abs(x) <= w:  # Barrier in the center
        return V1 
    elif abs(x) >= L: # Outside the wells
        return V2  
    else:
        return 0  


# Spatial range
x = np.linspace(-11, 11, 1000)
H = np.array([V(i, 100, 10**5, 0.1, 1) for i in x]) 

plt.plot(x, H, color='r')
plt.title("Potential Well")
plt.xlabel("x")
plt.ylabel("V(x)")
plt.show()

#Even Solution

# Initial conditions
dx = 0.001
E = 4.5
dE = 0.1  
MP = 0.3  # Matching point (at center of the barrier)
xl_0 = -1  
xr_0 = 1   
limit = 10**(-6)
last = 1

xl = np.arange(xl_0, MP+dx, dx)
xr = np.arange(xr_0, MP-dx, -dx)

psi_L = np.zeros(len(xl))
psi_R = np.zeros(len(xr))
psi_L[0] = -0.0001 * dx  
psi_R[0] = -0.0001 * dx

d_psi_L = -10
d_psi_R = 10


step=0
max_step=5000

#Shootong method, first we for the left side, then for the right one
while abs(d_psi_L-d_psi_R) > limit and step<max_step:
    step+=1
    for i in range(1, len(psi_L)-1):
        psi_L[i+1] = 2*psi_L[i]-psi_L[i-1]-2*(E-V(xl[i+1], 100, 10**5, 0.1, 1))*psi_L[i]*dx**2

    for i in range(1, len(psi_R)-1):
        psi_R[i+1] = 2*psi_R[i]-psi_R[i-1] -2*(E-V(xr[i+1], 100, 10**5, 0.1, 1))*psi_R[i]*dx**2

    # Scale psi_R to match psi_L at the boundary
    scale = psi_L[-1]/psi_R[-1]
    psi_R = scale*psi_R

    # Compute derivatives at the matching point
    d_psi_L = (psi_L[-1]-psi_L[-2]) /dx
    d_psi_R = -(psi_R[-1]-psi_R[-2]) /dx

    # Adjust energy based on derivative difference
    if d_psi_L > d_psi_R:
        E += dE
    else:
        E -=dE
    
    if E!= last:
        dE=-dE/2
        last *=-1

    E+=dE
    limit=0.5*((d_psi_L+d_psi_R)*0.001)

# Merge solutions
psi_R = np.flip(psi_R)
xr = np.flip(xr)
psi = np.concatenate((psi_L, psi_R))
x = np.concatenate((xl, xr))

# print("Energy:", E)
# plt.plot(x, psi)
# plt.title("Wavefunction for the Potential")
# plt.xlabel("x")
# plt.ylabel("Psi(x)")
# plt.grid(True)
# plt.show()


print("Energy:", E)
plt.plot(xl, psi_L, color='red', label='Psi_L')
plt.plot(xr, psi_R, color='green', label='Psi_R')
plt.title("Wavefunction (Even)")
plt.xlabel("x")
plt.ylabel("Psi(x)")
plt.legend()
plt.grid(True)
plt.show()

###################### ODD solution #################

# Initial conditions
dx = 0.001
E = 4.5
dE = 0.1  
MP = 0.3  # Matching point (at center of the barrier)
xl_0 = -1  
xr_0 = 1   
limit = 10**(-6)
last = 1

xl = np.arange(xl_0, MP+dx, dx)
xr = np.arange(xr_0, MP-dx, -dx)

psi_L = np.zeros(len(xl))
psi_R = np.zeros(len(xr))
psi_L[0] = 0.0001 * dx  
psi_R[0] = -0.0001 * dx

d_psi_L = -10
d_psi_R = 10


step=0
max_step=5000

#Shootong method, first we for the left side, then for the right one
while abs(d_psi_L-d_psi_R) > limit and step<max_step:
    step+=1
    for i in range(1, len(psi_L)-1):
        psi_L[i+1] = 2*psi_L[i]-psi_L[i-1]-2*(E-V(xl[i+1], 100, 10**5, 0.1, 1))*psi_L[i]*dx**2

    for i in range(1, len(psi_R)-1):
        psi_R[i+1] = 2*psi_R[i]-psi_R[i-1] -2*(E-V(xr[i+1], 100, 10**5, 0.1, 1))*psi_R[i]*dx**2

    # Scale psi_R to match psi_L at the boundary
    scale = psi_L[-1]/psi_R[-1]
    psi_R *= scale

    # Compute derivatives at the matching point
    d_psi_L = (psi_L[-1]-psi_L[-2]) /dx
    d_psi_R = -(psi_R[-1]-psi_R[-2]) /dx

    # Adjust energy based on derivative difference
    if d_psi_L > d_psi_R:
        E += dE
    else:
        E -=dE
    
    if E!= last:
        dE=-dE/2
        last *=-1

    E+=dE
    limit=0.5*((d_psi_L+d_psi_R)*0.001)

# Merge solutions
psi_R = np.flip(psi_R)
xr = np.flip(xr)
psi = np.concatenate((psi_L, psi_R))
x = np.concatenate((xl, xr))

print("Energy:", E)
plt.plot(x, psi)
plt.title("Wavefunction for the Potential")
plt.xlabel("x")
plt.ylabel("Psi(x)")
plt.grid(True)
plt.show()


# print("Energy:", E)

# plt.plot(xl, psi_L, color='red', label='Psi_L')
# plt.plot(xr, psi_R, color='green', label='Psi_R')
# plt.title("Wavefunction (Odd)")
# plt.xlabel("x")
# plt.ylabel("Psi(x)")
# plt.legend()
# plt.grid(True)
# plt.show()


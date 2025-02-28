import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Parameters
Q = 10.0       # Flow rate (L/h)
V = 100.0      # Reactor volume (L)
mu_max = 0.5   # Maximum cyanide degradation rate (h^-1)
Ks = 5.0       # Half-saturation constant for cyanide (mg/L)
Y = 0.2        # Yield coefficient for ammonia production
mu_nit = 0.3   # Maximum nitrification rate (h^-1)
K_nit = 2.0    # Half-saturation constant for ammonia (mg/L)
kd = 0.05      # Microbial decay rate (h^-1)
kd_nit = 0.05  # Nitrifying bacterial decay rate (h^-1)

# Inflow concentrations
C_in = 100.0   # Cyanide inflow concentration (mg/L)
N_in = 0.0     # Ammonia inflow concentration (mg/L)
X_in = 10.0    # Microbial population inflow (mg/L)
X_nit_in = 5.0 # Nitrifying bacterial population inflow (mg/L)

# Initial conditions
C0 = 0.0       # Initial cyanide concentration (mg/L)
N0 = 0.0       # Initial ammonia concentration (mg/L)
X0 = 1.0       # Initial microbial population (mg/L)
X_nit0 = 1.0   # Initial nitrifying bacterial population (mg/L)

# Time points
t = np.linspace(0, 100, 1000)  # Simulation time (0 to 100 hours)

# Model equations
def model(y, t):
    C, N, X, X_nit = y
    dCdt = (Q/V) * (C_in - C) - (mu_max * X * C) / (Ks + C)
    dNdt = (Q/V) * (N_in - N) + Y * (mu_max * X * C) / (Ks + C) - (mu_nit * X_nit * N) / (K_nit + N)
    dXdt = (Q/V) * (X_in - X) + mu_max * X * (C / (Ks + C)) - kd * X
    dX_nitdt = (Q/V) * (X_nit_in - X_nit) + mu_nit * X_nit * (N / (K_nit + N)) - kd_nit * X_nit
    return [dCdt, dNdt, dXdt, dX_nitdt]

# Solve ODE
y0 = [C0, N0, X0, X_nit0]
sol = odeint(model, y0, t)

# Extract results
C = sol[:, 0]
N = sol[:, 1]
X = sol[:, 2]
X_nit = sol[:, 3]

# Plot results
plt.figure(figsize=(12, 8))

plt.subplot(2, 2, 1)
plt.plot(t, C, label='Cyanide (C)')
plt.xlabel('Time (h)')
plt.ylabel('Concentration (mg/L)')
plt.title('Cyanide Concentration')
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(t, N, label='Ammonia (N)', color='orange')
plt.xlabel('Time (h)')
plt.ylabel('Concentration (mg/L)')
plt.title('Ammonia Concentration')
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(t, X, label='Microbial Population (X)', color='green')
plt.xlabel('Time (h)')
plt.ylabel('Population (mg/L)')
plt.title('Microbial Population')
plt.legend()

plt.subplot(2, 2, 4)
plt.plot(t, X_nit, label='Nitrifying Bacteria (X_nit)', color='red')
plt.xlabel('Time (h)')
plt.ylabel('Population (mg/L)')
plt.title('Nitrifying Bacterial Population')
plt.legend()

plt.tight_layout()
plt.show()
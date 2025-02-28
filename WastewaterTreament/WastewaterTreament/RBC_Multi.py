import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the system of differential equations for multiple bacterial species
def rbc_system_multispecies(t, y, Q_in, C_in, Q_out, mu_max, K_s, k_d, Y, Y_o, k_La, O_star):
    C = y[0]          # Cyanide concentration
    X = y[1:-1]       # Bacterial concentrations (one for each species)
    O = y[-1]         # Oxygen concentration

    # Calculate growth rates for each species
    mu = mu_max * (C / (K_s + C))  # Monod kinetics for each species

    # Cyanide mass balance
    dC_dt = (Q_in * C_in - Q_out * C) / V - np.sum((mu * X) / Y)

    # Bacterial growth dynamics for each species
    dX_dt = mu * X - k_d * X

    # Oxygen mass balance
    dO_dt = k_La * (O_star - O) - np.sum((mu * X) / Y_o)

    return [dC_dt, *dX_dt, dO_dt]

# Parameters (adjust based on actual system)
Q_in = 10.0          # Inflow rate (L/h)
C_in = 100.0         # Inflow cyanide concentration (mg/L)
Q_out = 10.0         # Outflow rate (L/h)
V = 1000.0           # Volume of the RBC system (L)
mu_max = np.array([0.5, 0.3])  # Maximum specific growth rates for each species (1/h)
K_s = np.array([10.0, 15.0])   # Half-saturation constants for each species (mg/L)
k_d = np.array([0.05, 0.03])   # Decay rates for each species (1/h)
Y = np.array([0.5, 0.4])       # Yield coefficients for each species (mg bacteria/mg cyanide)
Y_o = np.array([1.5, 1.2])     # Oxygen yield coefficients for each species (mg oxygen/mg cyanide)
k_La = 20.0          # Oxygen transfer coefficient (1/h)
O_star = 8.0         # Saturation oxygen concentration (mg/L)

# Initial conditions
C0 = 0.0             # Initial cyanide concentration (mg/L)
X0 = np.array([10.0, 5.0])  # Initial bacterial concentrations for each species (mg/L)
O0 = 8.0             # Initial oxygen concentration (mg/L)
y0 = [C0, *X0, O0]   # Combined initial conditions

# Time span for simulation (0 to 24 hours)
t_span = (0, 24)
t_eval = np.linspace(0, 24, 100)  # Time points for evaluation

# Solve the system of ODEs
sol = solve_ivp(rbc_system_multispecies, t_span, y0, args=(Q_in, C_in, Q_out, mu_max, K_s, k_d, Y, Y_o, k_La, O_star),
                t_eval=t_eval, method='RK45')

# Extract results
C = sol.y[0]         # Cyanide concentration over time
X = sol.y[1:-1]      # Bacterial concentrations over time (one array per species)
O = sol.y[-1]        # Oxygen concentration over time
t = sol.t            # Time points

# Plot results
plt.figure(figsize=(12, 10))

# Cyanide concentration
plt.subplot(3, 1, 1)
plt.plot(t, C, label='Cyanide Concentration (mg/L)', color='blue')
plt.xlabel('Time (h)')
plt.ylabel('Cyanide (mg/L)')
plt.title('Cyanide Concentration Over Time')
plt.legend()
plt.grid()

# Bacterial concentrations
plt.subplot(3, 1, 2)
for i, X_i in enumerate(X):
    plt.plot(t, X_i, label=f'Bacterial Species {i+1} (mg/L)')
plt.xlabel('Time (h)')
plt.ylabel('Bacteria (mg/L)')
plt.title('Bacterial Concentrations Over Time')
plt.legend()
plt.grid()

# Oxygen concentration
plt.subplot(3, 1, 3)
plt.plot(t, O, label='Oxygen Concentration (mg/L)', color='red')
plt.xlabel('Time (h)')
plt.ylabel('Oxygen (mg/L)')
plt.title('Oxygen Concentration Over Time')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()
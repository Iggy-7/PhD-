import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the system of differential equations
def rbc_system(t, y, Q_in, C_in, Q_out, mu_max, K_s, k_d, Y, Y_o, k_La, O_star):
    C, X, O = y  # Unpack variables: cyanide concentration, bacterial concentration, oxygen concentration

    # Monod kinetics for bacterial growth
    mu = mu_max * (C / (K_s + C))

    # Cyanide mass balance
    dC_dt = (Q_in * C_in - Q_out * C) / V - (mu * X) / Y

    # Bacterial growth dynamics
    dX_dt = mu * X - k_d * X

    # Oxygen mass balance
    dO_dt = k_La * (O_star - O) - (mu * X) / Y_o

    return [dC_dt, dX_dt, dO_dt]

# Parameters (adjust based on actual system)
Q_in = 10.0  # Inflow rate (L/h)
C_in = 100.0  # Inflow cyanide concentration (mg/L)
Q_out = 10.0  # Outflow rate (L/h)
V = 1000.0    # Volume of the RBC system (L)
mu_max = 0.5  # Maximum specific growth rate of bacteria (1/h)
K_s = 10.0    # Half-saturation constant for cyanide (mg/L)
k_d = 0.05    # Bacterial decay rate (1/h)
Y = 0.5       # Yield coefficient (mg bacteria/mg cyanide)
Y_o = 1.5     # Oxygen yield coefficient (mg oxygen/mg cyanide)
k_La = 20.0   # Oxygen transfer coefficient (1/h)
O_star = 8.0  # Saturation oxygen concentration (mg/L)

# Initial conditions
C0 = 0.0      # Initial cyanide concentration (mg/L)
X0 = 10.0     # Initial bacterial concentration (mg/L)
O0 = 8.0      # Initial oxygen concentration (mg/L)
y0 = [C0, X0, O0]

# Time span for simulation (0 to 24 hours)
t_span = (0, 24)
t_eval = np.linspace(0, 24, 100)  # Time points for evaluation

# Solve the system of ODEs
sol = solve_ivp(rbc_system, t_span, y0, args=(Q_in, C_in, Q_out, mu_max, K_s, k_d, Y, Y_o, k_La, O_star),
                t_eval=t_eval, method='RK45')

# Extract results
C = sol.y[0]  # Cyanide concentration over time
X = sol.y[1]  # Bacterial concentration over time
O = sol.y[2]  # Oxygen concentration over time
t = sol.t     # Time points

# Plot results
plt.figure(figsize=(12, 8))

# Cyanide concentration
plt.subplot(3, 1, 1)
plt.plot(t, C, label='Cyanide Concentration (mg/L)', color='blue')
plt.xlabel('Time (h)')
plt.ylabel('Cyanide (mg/L)')
plt.title('Cyanide Concentration Over Time')
plt.legend()
plt.grid()

# Bacterial concentration
plt.subplot(3, 1, 2)
plt.plot(t, X, label='Bacterial Concentration (mg/L)', color='green')
plt.xlabel('Time (h)')
plt.ylabel('Bacteria (mg/L)')
plt.title('Bacterial Concentration Over Time')
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
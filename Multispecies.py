import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib.pyplot as plt

data = pd.read_csv('mining_data.csv') 
time_data = data['time'].values
C_data = data['cyanide'].values
O_data = data['oxygen'].values

# Multispecies System Function
def multispecies_system(t, y, Q_in, C_in, O_in, V, D_C, D_O, mu_max, K_s, K_o, k_d, Y, Y_o, alpha, beta, dx, N, n_species):
    # Unpack variables
    C = y[:N]          
    O = y[N:2*N]       
    X = y[2*N:].reshape(n_species, N)  

    # Initialize derivatives
    dC_dt = np.zeros(N)
    dO_dt = np.zeros(N)
    dX_dt = np.zeros((n_species, N))

    # Diffusion terms 
    for i in range(1, N-1):
        dC_dt[i] = D_C * (C[i+1] - 2*C[i] + C[i-1]) / dx**2
        dO_dt[i] = D_O * (O[i+1] - 2*O[i] + O[i-1]) / dx**2
        for j in range(n_species):
            dX_dt[j, i] = D_X[j] * (X[j, i+1] - 2*X[j, i] + X[j, i-1]) / dx**2

    # Boundary conditions (Neumann no-flux)
    dC_dt[0] = dC_dt[1]    # dC/dx=0 at x=0
    dC_dt[-1] = dC_dt[-2]  # dC/dx=0 at x=L
    dO_dt[0] = dO_dt[1]    # dO/dx=0 at x=0
    dO_dt[-1] = dO_dt[-2]  # dO/dx=0 at x=L
    for j in range(n_species):
        dX_dt[j, 0] = dX_dt[j, 1]      # dX/dx=0 at x=0
        dX_dt[j, -1] = dX_dt[j, -2]    # dX/dx=0 at x=L

    # Biochemical reactions
    for i in range(n_species):
        competition = np.prod([(1 - alpha[i, j] * X[j]) for j in range(n_species) if j != i])
        cooperation = 1 + np.sum([beta[i, j] * X[j] for j in range(n_species) if j != i])
        mu = mu_max[i] * (C / (K_s[i] + C)) * (O / (K_o[i] + O)) * competition + cooperation
        dX_dt[i] += mu * X[i] - k_d[i] * X[i]
        dC_dt -= (mu * X[i]) / Y[i]
        dO_dt -= (mu * X[i]) / Y_o[i]

    return np.concatenate([dC_dt, dO_dt, dX_dt.flatten()])

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

# Parameters 
Q_in = 10.0
C_in = 100.0
O_in = 8.0
V = 1000.0
D_C = 1e-9
D_O = 2e-9
D_X = [1e-10, 1e-10]
K_o = [1.0, 2.0]
K_s = [1.0,2.0]
Y_o = [1.5, 1.2]
alpha = np.array([[0.0, 0.1], [0.05, 0.0]])  
beta = np.array([[0.0, 0.02], [0.01, 0.0]])   
dx = 0.01
N = 100
n_species = 2
t_span = (0, 24)
k_d = [0.01,0.1]
Y = [1.5,1.2]
mu_max =[1.0,2.0]
# Initial conditions
C0 = np.zeros(N)
C0[0] = C_in
O0 = np.full(N, O_in)
X0 = np.array([np.full(N, 10.0), np.full(N, 5.0)])
y0 = np.concatenate([C0, O0, X0.flatten()])

def multispecies_system(t, y, Q_in, C_in, O_in, V, D_C, D_O, mu_max, K_s, K_o, k_d, Y, Y_o, alpha, beta, dx, N, n_species):
    C = y[:N]
    O = y[N:2*N]
    X = y[2*N:].reshape(n_species, N)
    
    dC_dt = np.zeros(N)
    dO_dt = np.zeros(N)
    dX_dt = np.zeros((n_species, N))
    
    # Diffusion terms
    for i in range(1, N-1):
        dC_dt[i] = D_C * (C[i+1] - 2*C[i] + C[i-1]) / dx**2
        dO_dt[i] = D_O * (O[i+1] - 2*O[i] + O[i-1]) / dx**2
        for j in range(n_species):
            dX_dt[j,i] = D_X[j] * (X[j,i+1] - 2*X[j,i] + X[j,i-1]) / dx**2
    
    # Boundary conditions
    dC_dt[0] = 0
    dC_dt[-1] = D_C * (2*C[-2] - 2*C[-1]) / dx**2
    dO_dt[0] = 0
    dO_dt[-1] = D_O * (2*O[-2] - 2*O[-1]) / dx**2
    for j in range(n_species):
        dX_dt[j,0] = D_X[j] * (2*X[j,1] - 2*X[j,0]) / dx**2
        dX_dt[j,-1] = D_X[j] * (2*X[j,-2] - 2*X[j,-1]) / dx**2
    
    # Biochemical reactions
    for i in range(n_species):
        competition = np.prod([(1 - alpha[i,j]) * X[j] for j in range(n_species) if j != i])
        cooperation = 1 + np.sum([beta[i,j] * X[j] for j in range(n_species) if j != i])
        mu = mu_max[i] * (C / (K_s[i] + C)) * (O / (K_o[i] + O)) * competition + cooperation
        dX_dt[i] += mu * X[i] - k_d[i] * X[i]
        dC_dt -= (mu * X[i]) / Y[i]
        dO_dt -= (mu * X[i]) / Y_o[i]
    
    return np.concatenate([dC_dt, dO_dt, dX_dt.flatten()])

def objective_function(params):
    mu_max_new, K_s_new, k_d_new, Y_new = params
    mu_max_opt = np.array([mu_max_new, 0.4])
    K_s_opt = np.array([K_s_new, 15.0])
    k_d_opt = np.array([k_d_new, 0.03])
    Y_opt = np.array([Y_new, 0.4])
    
    sol = solve_ivp(multispecies_system, t_span, y0,
                   args=(Q_in, C_in, O_in, V, D_C, D_O, mu_max_opt, K_s_opt, K_o, k_d_opt, Y_opt, Y_o, alpha, beta, dx, N, n_species),
                   t_eval=time_data, method='RK45')
    
    C_model = sol.y[:N,:].mean(axis=0)
    O_model = sol.y[N:2*N,:].mean(axis=0)
    error = np.sum((C_model - C_data)**2) + np.sum((O_model - O_data)**2)
    return error

# Load experimental data
data = pd.read_csv('mining_data.csv')
time_data = data['time'].values
C_data = data['cyanide'].values
O_data = data['oxygen'].values

# Run optimization
initial_guess = [0.5, 10.0, 0.05, 0.5]
result = minimize(objective_function, initial_guess, method='L-BFGS-B',
                 bounds=[(0.1, 1.0), (5.0, 20.0), (0.01, 0.1), (0.1, 1.0)])
optimized_params = result.x
print("Optimized Parameters:", optimized_params)

# Update parameters with optimized values
mu_max_opt = np.array([optimized_params[0], mu_max[1]])
K_s_opt = np.array([optimized_params[1], K_s[1]])
k_d_opt = np.array([optimized_params[2], k_d[1]])
Y_opt = np.array([optimized_params[3], Y[1]])

# ReRun simulation with optimized parameters
sol = solve_ivp(multispecies_system, t_span, y0,
                args=(Q_in, C_in, O_in, V, D_C, D_O, mu_max_opt, K_s_opt, K_o, k_d_opt, Y_opt, Y_o, alpha, beta, dx, N, n_species),
                t_eval=np.linspace(0, 24, 100), method='RK45')

# Extract results
t = sol.t
C = sol.y[:N, :].mean(axis=0)  # Spatial average of cyanide
O = sol.y[N:2*N, :].mean(axis=0)
X = sol.y[2*N:, :].reshape(n_species, N, -1).mean(axis=1)  # Avg bacteria per species

# Plot
plt.figure(figsize=(12, 8))
plt.subplot(3, 1, 1)
plt.plot(t, C, 'b-', label='Model')
plt.scatter(time_data, C_data, color='r', label='Data')
plt.ylabel('Cyanide (mg/L)')
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(t, O, 'g-', label='Model')
plt.scatter(time_data, O_data, color='orange', label='Data')
plt.ylabel('Oxygen (mg/L)')
plt.legend()

plt.subplot(3, 1, 3)
for i in range(n_species):
    plt.plot(t, X[i], label=f'Species {i+1}')
plt.xlabel('Time (h)')
plt.ylabel('Bacteria (mg/L)')
plt.legend()

plt.tight_layout()
plt.show()
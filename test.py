import numpy as np
import matplotlib.pyplot as plt

from sklearn.utils import resample

from sse import sse

# Number of atoms in the chain
N = 2
Nb = N

# MC Cycles
therm_cycles = int(1e3)
mc_cycles = int(1e4)
bootstrap_cycles = mc_cycles

beta_vals = np.power(2.0, [-1, 0, 1, 2, 3, 4])
mean_n = np.zeros(beta_vals.shape)
std_n = np.zeros(beta_vals.shape)
mean_E = np.zeros(beta_vals.shape)
std_E = np.zeros(beta_vals.shape)

n_vals = np.zeros((len(beta_vals), mc_cycles))
for i, beta in enumerate(beta_vals):
    
    print(f"beta = {i+1}/{len(beta_vals)}", end='\r')
    
    rng = np.random.default_rng(np.random.MT19937())
    n_vals[i, :] = sse(beta, N, Nb, therm_cycles, mc_cycles, rng)
    
    mean_n_bs = np.zeros((bootstrap_cycles, mc_cycles))
    mean_E_bs = np.zeros((bootstrap_cycles, mc_cycles))
    
    for bs_cycle in range(bootstrap_cycles):
        mean_n_bs[bs_cycle, :] = resample(n_vals[i, :])
        mean_E_bs[bs_cycle, :] = - mean_n_bs[bs_cycle, :] / (beta * N)
    
    mean_E[i] = np.mean(mean_E_bs)
    std_E[i] = np.mean(np.std(mean_E_bs, axis=1))
    
    mean_n[i] = np.mean(mean_n_bs)
    std_n[i] = np.mean(np.std(mean_n_bs, axis=1))

for i, beta in enumerate(beta_vals):
    print(f"beta : {beta} | <E> = {mean_E[i]:.6f} | std_E = {std_E[i]:.6f} | <n> = {mean_n[i]:.2f} | std_n = {std_n[i]:.2f}")



import numpy as np
import matplotlib
import matplotlib.pyplot as plt

plt.style.use('seaborn')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Test temperatures
# T_vals = 6
# beta_vals = np.array([0.5, 1.0, 2.0, 4.0, 8.0, 16.0])
# T = 1.0 / beta_sim

# Custom temperatures
Ti = 0.05
Tf = 4.0
T_vals = 50
T = np.zeros(T_vals)
for i in range(T_vals):
    T[i] = (Tf - Ti) * (i + 1) / (T_vals)

beta = 1.0 / T

# Read the output file
FILENAME = "2D_L4_FM_Heisenberg_delta1.0_h0.0.csv" # "2D_L8_FM_Heisenberg_delta1.0_h0.0.csv"

E = np.zeros(T_vals)
E_std = np.zeros(T_vals)
C = np.zeros(T_vals)
C_std = np.zeros(T_vals)
m = np.zeros(T_vals)
m_std = np.zeros(T_vals)
m2 = np.zeros(T_vals)
m2_std = np.zeros(T_vals)
m_sus = np.zeros(T_vals)
m_sus_std = np.zeros(T_vals)
ms = np.zeros(T_vals)
ms_std = np.zeros(T_vals)
m2s = np.zeros(T_vals)
m2s_std = np.zeros(T_vals)
binder = np.zeros(T_vals)
binder_std = np.zeros(T_vals)

with open(FILENAME, "r") as f:
    header = f.readline()
    for j in range(T_vals):
        _, _, _, _, E[j], E_std[j], C[j], \
        C_std[j], m[j], m_std[j], m2[j], \
        m2_std[j], _, _, ms[j], ms[j], \
        m2s[j], m2s_std[j], _, _, m_sus[j], \
        m_sus_std[j], binder[j], binder_std[j] = [float(x) for x in f.readline().strip().split(",")]

# Plot the variables
fig, _ = plt.subplots(3, 3, sharex=True)
plt.figure(1)
fig.suptitle(FILENAME)

plt.subplot(3, 3, 1)
plt.errorbar(T, E, E_std, fmt=".--", label="SSE")
plt.xlabel(r"$T$")
plt.ylabel(r"$\langle E \rangle$")
plt.legend()

plt.subplot(3, 3, 2)
plt.errorbar(T, C, C_std, fmt=".--", label="SSE")
plt.xlabel(r"$T$")
plt.ylabel(r"$C$")
plt.legend()

plt.subplot(3, 3, 3)
plt.errorbar(T, binder, binder_std, fmt=".--", label="SSE")
plt.xlabel(r"$T$")
plt.ylabel(r"Binder Parameter")
plt.legend()

plt.subplot(3, 3, 4)
plt.errorbar(T, m2, ms_std, fmt=".--", label="SSE")
plt.xlabel(r"$T$")
plt.ylabel(r"$\langle m_s \rangle$")
plt.legend()

plt.subplot(3, 3, 5)
plt.errorbar(T, m2s, m2s_std, fmt=".--", label="SSE")
plt.xlabel(r"$T$")
plt.ylabel(r"$\langle m_s^2 \rangle$")
plt.legend()

plt.subplot(3, 3, 7)
plt.errorbar(T, m, m_std, fmt=".--", label="SSE")
plt.xlabel(r"$T$")
plt.ylabel(r"$\langle m \rangle$")
plt.legend()

plt.subplot(3, 3, 8)
plt.errorbar(T, m2, m2_std, fmt=".--", label="SSE")
plt.xlabel(r"$T$")
plt.ylabel(r"$\langle m^2 \rangle$")
plt.legend()

plt.subplot(3, 3, 9)
plt.errorbar(T, m_sus, m_sus_std, fmt=".--", label="SSE")
plt.xlabel(r"$T$")
plt.ylabel(r"$\chi$")
plt.legend()

plt.show()

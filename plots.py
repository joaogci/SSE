import numpy as np
import matplotlib.pyplot as plt
import matplotlib

plt.style.use('seaborn')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# filenames = ["1D_heisenberg_L8.csv", "1D_heisenberg_L16.csv", "1D_heisenberg_L32.csv", 
#              "1D_heisenberg_L64.csv", "1D_heisenberg_L128.csv"]
# filenames = ["2D_heisenberg_L4.csv", "2D_heisenberg_L8.csv"]
filenames = ["2D_heisenberg_L4_more_T.csv", "2D_heisenberg_L8_more_T.csv"]

T_vals = 10
L_vals = len(filenames)
L = np.array([4, 8])
# L = np.array([8, 16, 32, 64, 128])
T = np.zeros((L_vals, T_vals))
E = np.zeros((L_vals, T_vals))
C = np.zeros((L_vals, T_vals))
m = np.zeros((L_vals, T_vals))
m2 = np.zeros((L_vals, T_vals))
m_s = np.zeros((L_vals, T_vals))
m2_s = np.zeros((L_vals, T_vals))
n = np.zeros((L_vals, T_vals))

sus = np.zeros((L_vals, T_vals))

for i, filename in enumerate(filenames):
    with open("some_results/" + filename, "r") as file:
        header = file.readline()
        
        for j in range(T_vals):
            T[i, j], E[i, j], C[i, j], m[i, j], m2[i, j], m_s[i, j], m2_s[i, j], sus[i, j], n[i, j] = [float(x) for x in file.readline().strip().split(",")]
        T[i, :] = 1.0 / T[i, :]

plt.figure(1)
for i in range(L_vals):
    plt.plot(T[i, :], E[i, :], ".-", label=rf"$L=${L[i]}")
plt.xlabel(r"$T/J$", fontsize=16)
plt.ylabel(r"$\langle E \rangle$", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)

plt.figure(2)
for i in range(L_vals):
    plt.plot(T[i, :], C[i, :], ".-", label=rf"$L=${L[i]}")
plt.xlabel(r"$T/J$", fontsize=16)
plt.ylabel(r"$\langle C \rangle$", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)

plt.figure(3)
for i in range(L_vals):
    plt.plot(T[i, :], m2[i, :], ".-", label=rf"$L=${L[i]}")
plt.xlabel(r"$T/J$", fontsize=16)
plt.ylabel(r"$\langle m^2 \rangle$", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)

plt.figure(4)
for i in range(L_vals):
    plt.plot(T[i, :], m2_s[i, :], ".-", label=rf"$L=${L[i]}")
plt.xlabel(r"$T/J$", fontsize=16)
plt.ylabel(r"$\langle m^2 \rangle_s$", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)

plt.figure(5)
for i in range(L_vals):
    plt.plot(T[i, :], sus[i, :], ".-", label=rf"$L=${L[i]}")
plt.xlabel(r"$T/J$", fontsize=16)
plt.ylabel(r"$\langle \chi \rangle_s$", fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=16)

plt.show()


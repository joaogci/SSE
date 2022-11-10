import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from post_processing import read_sse_output

plt.style.use('seaborn')
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Read the output file
FILENAME = "output.csv"
sim_info, sampled = read_sse_output(FILENAME)

# Plot the variables
fig, _ = plt.subplots(3, 3, sharex=True)
plt.figure(1)
fig.suptitle(FILENAME)

plt.subplot(3, 3, 1)
plt.errorbar(sampled["T"], sampled["E"], sampled["E_std"], fmt=".--")
plt.xlabel(r"$T$")
plt.ylabel(r"$\langle E \rangle$")
plt.legend()

plt.subplot(3, 3, 2)
plt.errorbar(sampled["T"], sampled["C"], sampled["C_std"], fmt=".--")
plt.xlabel(r"$T$")
plt.ylabel(r"$C$")
plt.legend()

plt.subplot(3, 3, 4)
plt.errorbar(sampled["T"], sampled["m2"], sampled["m2_std"], fmt=".--")
plt.xlabel(r"$T$")
plt.ylabel(r"$\langle m_s \rangle$")
plt.legend()

plt.subplot(3, 3, 5)
plt.errorbar(sampled["T"], sampled["m2s"], sampled["m2s_std"], fmt=".--")
plt.xlabel(r"$T$")
plt.ylabel(r"$\langle m_s^2 \rangle$")
plt.legend()

plt.subplot(3, 3, 6)
plt.errorbar(sampled["T"], sampled["S_mean"], sampled["S_std"], fmt=".--")
plt.xlabel(r"$T$")
plt.ylabel(r"$S(\pi)$")
plt.legend()

plt.subplot(3, 3, 7)
plt.errorbar(sampled["T"], sampled["m"], sampled["m_std"], fmt=".--")
plt.xlabel(r"$T$")
plt.ylabel(r"$\langle m \rangle$")
plt.legend()

plt.subplot(3, 3, 8)
plt.errorbar(sampled["T"], sampled["m2"], sampled["m2_std"], fmt=".--")
plt.xlabel(r"$T$")
plt.ylabel(r"$\langle m^2 \rangle$")
plt.legend()

plt.subplot(3, 3, 9)
plt.errorbar(sampled["T"], sampled["m_sus"], sampled["m_sus_std"], fmt=".--")
plt.xlabel(r"$T$")
plt.ylabel(r"$\chi$")
plt.legend()

plt.figure(2)
for i in range(sim_info["L"]):    
    plt.errorbar(sampled["T"], sampled["corr_mean"][:, i], sampled["corr_std"][:, i], fmt=".--", label=f"i={i}")
    plt.xlabel(r"$T$")
    plt.ylabel(r"$\chi$")
    plt.legend()

plt.show()

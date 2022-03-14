import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

# -- Functions --
def P_ins(n):
    return np.min([1, Nb * beta / (2 * (M - n))])

def P_rem(n):
    return np.min([1, 2 * (M - n + 1) / (beta * Nb)])

# -- System Parameters --
# Number of spins
N = 2
# Dimensions
d = 1
# Number of Bonds
Nb = N * d
M = 100

# -- Simulation Parameters --
mc_cycles = 3000 # int(1e4)
rng = npr.default_rng(npr.MT19937(seed=1337))


def sse(beta):
    distrib = np.zeros(M)
    # -- Init Arrays -- 
    # M = np.max([4, N // 4])
    # len_op = 12
    opstring = np.zeros(M, dtype=np.int64)
    vertex_list = np.zeros(4 * M, dtype=np.int64) - 1
    vertex_init = np.zeros(N, dtype=np.int64) - 1
    vertex_last = np.zeros(N, dtype=np.int64) - 1

    bond_site = np.zeros((2, N), dtype=np.int64)
    spin = np.zeros(N)

    # -- Main program --

    # for M in range(1, M_max):
    n = 0

    # Init lattice 
    for b in range(N):
        bond_site[0, b] = b
        bond_site[1, b] = np.mod(b + 1, N) 

    # Init spins 
    for i in range(N):
        spin[i] = 1 if rng.random() <= 0.5 else -1

    n_vals = np.zeros(mc_cycles)

    for t in range(mc_cycles):
        n = np.count_nonzero(opstring)
        
        # diag_update()
        for p in range(M):
            if opstring[p] == 0:
                b = rng.integers(0, Nb)
                if spin[bond_site[0, b]] == spin[bond_site[1, b]]:
                    continue
                
                if rng.random() <= np.min([1, Nb * beta / (2 * (M - n))]):
                    opstring[p] = 2 * b
                    n += 1
            elif np.mod(opstring[p], 2) == 0:
                if rng.random() <= np.min([1, 2 * (M - n + 1) / (beta * Nb)]):
                    opstring[p] = 0
                    n += -1
            else:
                b = opstring[p] // 2
                spin[bond_site[0, b]] = - spin[bond_site[0, b]]
                spin[bond_site[1, b]] = - spin[bond_site[1, b]]

        vertex_list = np.zeros(4 * M, dtype=np.int64) - 1
        vertex_init = np.zeros(N, dtype=np.int64) - 1
        vertex_last = np.zeros(N, dtype=np.int64) - 1

        # create vertex_list
        for p in range(M):
            if opstring[p] == 0:
                # vertex_list[4 * p:4 * p + 4] = -1
                continue
            
            v0 = 4 * p
            bond = opstring[p] // 2
            i1 = bond_site[0, bond]
            i2 = bond_site[1, bond]
            
            v1 = vertex_last[i1]
            v2 = vertex_last[i2]
            
            if v1 != -1:
                vertex_list[v1] = v0
                vertex_list[v0] = v1
            else:
                vertex_init[i1] = v0
            
            if v2 != -1:
                vertex_list[v2] = v0 + 1
                vertex_list[v0 + 1] = v2
            else:
                vertex_init[i2] = v0 + 1
                
            vertex_last[i1] = v0 + 2
            vertex_last[i2] = v0 + 3

        for i in range(N):
            f = vertex_init[i]
            if f != -1:
                l = vertex_last[i]
                vertex_list[f] = l
                vertex_list[l] = f
                
        # loop_update()
        for v0 in range(0, 4 * M, 2):
            if vertex_list[v0] < 0:
                continue
            v_in = v0
            if rng.random() <= 0.5:
                # Transvese the loop, for all v in loop, set X[v]=-2
                # Flip the loop (change the operator types opstring[p=v/4] while loop is traversed)
                while True:
                    opstring[v_in//4] = opstring[v_in//4] ^ 1
                    vertex_list[v_in] = -2
                    v_out = v_in ^ 1
                    v_in = vertex_list[v_out]
                    vertex_list[v_out] = -2
                    
                    if v_in == v0:# or v_in < 0:
                        break
            else: 
                # Transverse the loop, for all v in loop, set X[v]=-1
                while True:
                    vertex_list[v_in] = -1
                    v_out = v_in ^ 1
                    v_in = vertex_list[v_out]
                    vertex_list[v_out] = -1
                    
                    if v_in == v0:# or v_in < 0:
                        break
        
        for i in range(N):
            if vertex_init[i] != -1:
                if vertex_list[vertex_init[i]] == -2:
                    spin[i] = - spin[i]
            else:
                if rng.random() <= 0.5:
                    spin[i] = - spin[i]

        n = np.count_nonzero(opstring)
        n_vals[t] = n
        # distrib[n] += 1
        print(f"{t=}", end='\r')
    print()
    return n_vals

beta_exp = np.array([-1, 0, 1, 2, 3, 4])
beta_vals = np.power(2.0, beta_exp)
n_vals = np.zeros((len(beta_vals), mc_cycles))
mean_n = np.zeros(beta_vals.shape)
std_n = np.zeros(beta_vals.shape)
mean_E = np.zeros(beta_vals.shape)
std_E = np.zeros(beta_vals.shape)

for i, beta in enumerate(beta_vals):
    print(f"Starting beta {i=}/{len(beta_vals) - 1}")
    n_vals[i, :] = sse(beta)

    mean_n[i] = np.mean(n_vals[i, 2500:])
    std_n[i] = np.std(n_vals[i, 2500:])
    
    E = - n_vals[i, 2500:] / (beta * N)
    mean_E[i] = np.mean(E)
    std_E[i] = np.std(E)

for i, beta in enumerate(beta_vals):
    print(f"beta : {beta} | <E> = {mean_E[i]:.6f} | std_E = {std_E[i]:.6f} | <n> = {mean_n[i]:.2f} | std_n = {std_n[i]:.2f}")


plt.figure(1)
plt.title(f"{mc_cycles=}")
for i in range(len(beta_vals)):
    plt.plot(np.arange(mc_cycles), n_vals[i, :], label=f"beta = {beta_vals[i]}")
    # print(f"beta = {beta_vals[i]} -> mean = {mean_n[i]}; std = {std_n[i]}")
plt.xlabel("MC Sweeps")
plt.ylabel(r"$n$")
plt.legend()

# plt.figure(3)
# plt.title(f"E")
# plt.plot(beta_vals, E, 'o')
# plt.xlabel(r"$\beta$")
# plt.ylabel(r"$E$")

# max_idx = np.where(distrib == np.max(distrib))[0][0]

# plt.figure(2)
# plt.title(f"{mc_cycles=}; {beta=}")
# plt.plot(np.arange(max_idx-130,max_idx+130), distrib[max_idx-130:max_idx+130]/np.sum(distrib[max_idx-130:max_idx+130]))
# plt.xlabel(r"$n$")
# plt.ylabel(r"$P(n)$")

plt.show()



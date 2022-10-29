import numpy as np
import sys
from itertools import product

SAVE_DIR = "tmp/"
FILENAME = str(sys.argv[1])
N_LEGS = 4

# Read the input file to get the system parameters
with open(FILENAME, "r") as f:
    f.readline()
    line = f.readline().strip().split(", ")

d = int(line[0])
L = int(line[1])
N = L ** d

S = float(line[2])
Sz = [m for m in np.arange(- S, S + 1, 1)]

delta = float(line[3])
h = float(line[4])
epsilon = float(line[5])

hb = h / (2 * d)
C = S**2 * np.abs(delta) + 2 * S * hb + epsilon

SAVE_FILENAME = "vtx_S" + str(S) + "_delta" + str(delta) + "_h" + str(h) + "_epsilon" + str(epsilon) + ".txt"
print(SAVE_FILENAME)

def Sp_op(state):
    if state[1] + 1 == state[0]:
        return np.sqrt((S - state[1]) * (S + state[1] + 1))
    return 0.0

def Sm_op(state):
    if state[1] - 1 == state[0]:
        return np.sqrt((S + state[1]) * (S - state[1] + 1))
    return 0.0
    
def Sz_op(state):
    if state[0] == state[1]:
        return state[0]
    return 0.0

def H(state1, state2):
    # < state1 | H | state2 >
    # Hb = 0.5 * (S+_i S-_j + S-_i S+_j) + C - delta Sz_i Sz_j + hb (Sz_i + Sz_j)
    # vtx type: 
    #   0 - diagonal
    #   1 - off_diagonal S+S-
    #   2 - off_diagonal S-S+
    
    term = 0.0
    allowed = False
    vtx_type = 0
    
    spin_i = (state1[0], state2[0])
    spin_j = (state1[1], state2[1])

    term1 = 0.5 * (Sp_op(spin_i) * Sm_op(spin_j))
    if term1 != 0.0:
        vtx_type = 1
        allowed = True
    term2 = 0.5 * (Sm_op(spin_i) * Sp_op(spin_j))
    if term2 != 0.0:
        if vtx_type == 1:
            print("error")
        vtx_type = 2
        allowed = True
        
    term3 = - delta * Sz_op(spin_i) * Sz_op(spin_j)
    if state1 == state2:
        allowed = True

    if vtx_type == 0:
        term += hb * (Sz_op(spin_i) + Sz_op(spin_j))
        term += C
    
    term += term1 + term2 + term3
    return term, allowed, vtx_type

# Generate the vertex information
# Generate all of the possible vertices
all_vertex_state = list(product(Sz, repeat=N_LEGS))
allowed_vertex_state = list()

# Possible updates
n_updates = 2 # +1 and -1 updates
updates = [-1, +1]
vtx = list()

# Select the allowed vertices
N_DIAGRAMS = 0
for vertex in all_vertex_state:
    term, allowed, vtx_type = H(vertex[2:], vertex[0:2])

    if allowed:
        allowed_vertex_state.append(list(vertex))
        
        vtx.append(dict())
        vtx[N_DIAGRAMS]["idx"] = N_DIAGRAMS
        vtx[N_DIAGRAMS]["type"] = vtx_type
        vtx[N_DIAGRAMS]["spin"] = list(vertex)
        vtx[N_DIAGRAMS]["H"] = np.around(term, decimals=10)
        
        vtx[N_DIAGRAMS]["new_vtx"] = [-1 + np.zeros((N_LEGS, N_LEGS), dtype=int), -1 + np.zeros((N_LEGS, N_LEGS), dtype=int)]
        vtx[N_DIAGRAMS]["prob"] = [np.zeros((N_LEGS, N_LEGS), dtype=float), np.zeros((N_LEGS, N_LEGS), dtype=float)]
        
        N_DIAGRAMS += 1

# Create transition table
for i in range(N_DIAGRAMS):
    for e in range(N_LEGS):
        for j, update in enumerate(updates):
            if np.abs(vtx[i]["spin"][e] + update) > S:
                continue

            for x in range(N_LEGS):
                new_state = vtx[i]["spin"].copy()
                new_state[e] = vtx[i]["spin"][e] + update
                
                if e == x:
                    new_state[x] = vtx[i]["spin"][x]
                elif e < 2 and x < 2:
                    new_state[x] = (new_state[2] + new_state[3]) - new_state[e]
                elif e >= 2 and x >= 2:
                    new_state[x] = (new_state[0] + new_state[1]) - new_state[e]
                elif e < 2 and x >= 2:
                    new_state[x] = (new_state[(e + 1) % 2] + new_state[e]) - new_state[(x - 2 + 1) % 2 + 2]
                elif e >= 2 and x < 2:
                    new_state[x] = (new_state[(e - 2 + 1) % 2 + 2] + new_state[e]) - new_state[(x + 1) % 2]
                        
                if new_state in allowed_vertex_state:
                    vtx[i]["new_vtx"][j][e, x] = allowed_vertex_state.index(new_state)

# Generate probablity table
# General method (Solution B)
for i in range(N_DIAGRAMS):
    for li in range(N_LEGS):
        for j in range(n_updates):
            A = np.zeros((N_LEGS, N_LEGS))
            b = np.zeros(N_LEGS)
            dim = N_LEGS
            
            vtx_idx = vtx[i]["new_vtx"][j][li, :]
            update_idx = np.zeros(N_LEGS, dtype=int)
            
            for l, v in enumerate(vtx_idx):
                if vtx[v]["spin"][l] - vtx[i]["spin"][l] > 0:
                    update_idx[l] = 0
                elif vtx[v]["spin"][l] - vtx[i]["spin"][l] == 0:
                    update_idx[l] = j
                else:
                    update_idx[l] = 1  

                b[l] = vtx[v]["H"] if v != -1 else 0.0
                dim += -1 if v == -1 else 0
            key = np.argsort(-b)
            
            if b[key[0]] <= np.sum(b[key[1:]]):
                if dim == N_LEGS:
                    A[key[2], key[3]] = 1.0
                    A[key[1], key[3]] = 0.0

                    A[key[3], key[2]] = 1.0
                    A[key[3], key[1]] = 0.0
                
                A[key[0], key[1]] = (b[key[0]] + b[key[1]] - b[key[2]] - b[key[3]]) / 2.0 + A[key[2], key[3]]
                A[key[0], key[2]] = (b[key[0]] - b[key[1]] + b[key[2]] - b[key[3]]) / 2.0 + A[key[1], key[3]]
                A[key[0], key[3]] = b[key[3]] - (A[key[2], key[3]] + A[key[1], key[3]])
                A[key[1], key[2]] = (- b[key[0]] + b[key[1]] + b[key[2]] + b[key[3]]) / 2.0 - (A[key[2], key[3]] + A[key[1], key[3]])

                A[key[1], key[0]] = (b[key[0]] + b[key[1]] - b[key[2]] - b[key[3]]) / 2.0 + A[key[2], key[3]]
                A[key[2], key[0]] = (b[key[0]] - b[key[1]] + b[key[2]] - b[key[3]]) / 2.0 + A[key[1], key[3]]
                A[key[2], key[1]] = (- b[key[0]] + b[key[1]] + b[key[2]] + b[key[3]]) / 2.0 - (A[key[2], key[3]] + A[key[1], key[3]])
                A[key[3], key[0]] = b[key[3]] - (A[key[2], key[3]] + A[key[1], key[3]])
            else:
                A[key[0], key[0]] = b[key[0]] - b[key[1]] - b[key[2]] - b[key[3]]

                A[key[0], key[1]] = b[key[1]]
                A[key[0], key[2]] = b[key[2]]
                A[key[0], key[3]] = b[key[3]]

                A[key[1], key[0]] = b[key[1]]
                A[key[2], key[0]] = b[key[2]]
                A[key[3], key[0]] = b[key[3]]

            for e in range(N_LEGS):
                for x in range(N_LEGS):
                    if b[e] > 0.0:
                        vtx[vtx_idx[e]]["prob"][update_idx[e]][e, x] = A[e, x] / b[e]

# Cumulative probabilities
for i in range(N_DIAGRAMS):
    for j in range(n_updates):
        for e in range(N_LEGS):
            cum_prob = 0.0
            for x in range(N_LEGS):
                if vtx[i]["prob"][j][e, x] > 0.0:
                    cum_prob += vtx[i]["prob"][j][e, x]
                    vtx[i]["prob"][j][e, x] = cum_prob

# Save the results to file
with open(SAVE_DIR + SAVE_FILENAME, "w") as f:
    f.write(str(N_DIAGRAMS) + " " + str(n_updates) + " " + str(N_LEGS) + "\n")
    f.write("\n")
    
    for i in range(N_DIAGRAMS):
        f.write(str(i) + "\n")
        f.write(str(vtx[i]["type"]) + "\n")
        f.write(str(vtx[i]["H"]) + "\n")
        
        for l in range(N_LEGS):
            f.write(str(int(2 * vtx[i]["spin"][l])) + " ")
        f.write("\n")
        
        for j in range(n_updates):
            for e in range(N_LEGS):
                for x in range(N_LEGS):
                    f.write(str(vtx[i]["new_vtx"][j][e, x]) + " ")
                f.write("\n")
                
            for e in range(N_LEGS):
                for x in range(N_LEGS):
                    f.write(str(vtx[i]["prob"][j][e, x]) + " ")
                f.write("\n")
        
        f.write("\n")

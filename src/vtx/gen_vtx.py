import numpy as np
import sys
from itertools import product

FILENAME = str(sys.argv[1])
N_LEGS = 4

# Read the input file to get the system parameters
with open(FILENAME, "r") as f:
    f.readline()
    line = f.readline().strip().split(", ")

d = int(line[0])
L = int(line[1])
N = L ** d
J = float(line[2])
delta = float(line[3])
h = float(line[4])
epsilon = float(line[5])

hb = h / (2 * d * J)
C = 0.25 * delta + hb + epsilon

vtx_name = "vtx_J" + str(J) + "_delta" + str(delta) + "_h" + str(h) + "_epsilon" + str(epsilon) + ".txt"
print(vtx_name)

def Sp(state):
    if state[1] < state[0]:
        return 1.0
    return 0.0

def Sm(state):
    if state[1] > state[0]:
        return 1.0
    return 0.0
    
def Sz(state):
    if state[0] == state[1]:
        return state[0] / 2.0
    return 0.0

# Hb = 0.5 * (S+_i S-_j + S-_i S+_j) + C - delta Sz_i Sz_j + hb (Sz_i + Sz_j)
def H(state1, state2):
    term = 0.0
    allowed = False
    diagonal = False
    
    spin_i = (state1[0], state2[0])
    spin_j = (state1[1], state2[1])

    term += 0.5 * (Sp(spin_i) * Sm(spin_j) + Sm(spin_i) * Sp(spin_j))
    if term != 0.0:
        diagonal = True
        
    term += - delta * Sz(spin_i) * Sz(spin_j)

    if term != 0.0:
        allowed = True

    if Sz(spin_i) != 0.0 and Sz(spin_j) != 0.0:
        term += hb * (Sz(spin_i) + Sz(spin_j))
        term += C
    
    return term, allowed, diagonal

def symmetry(vertex1, vertex2):
    if (vertex1[0] == vertex2[1] and vertex1[1] == vertex2[0]) and (vertex1[2] == vertex2[3] and vertex1[3] == vertex2[2]):
        return True
    elif (vertex1[0] == vertex2[2] and vertex1[1] == vertex2[3]) and (vertex1[2] == vertex2[0] and vertex1[3] == vertex2[1]):
        return True
    
    return False

# Generate the vertex information
# Generate all of the possible vertices
spin = [1, -1]
all_vertex_state = list(product(spin, repeat=N_LEGS))
allowed_vertex_state = list()

vtx = list()

# Select the valid vertices
N_DIAGRAMS = 0
for vertex in all_vertex_state:
    term, allowed, diagonal = H(vertex[0:2], vertex[2:])
    
    if allowed:
        allowed_vertex_state.append(list(vertex))
        
        vtx.append(dict())
        vtx[N_DIAGRAMS]["idx"] = N_DIAGRAMS
        vtx[N_DIAGRAMS]["type"] = 1 if diagonal else 0
        vtx[N_DIAGRAMS]["spin"] = list(vertex)
        vtx[N_DIAGRAMS]["H"] = term
        vtx[N_DIAGRAMS]["sym_to"] = list()
        
        vtx[N_DIAGRAMS]["new_vtx"] = np.zeros((N_LEGS, N_LEGS), dtype=int)
        vtx[N_DIAGRAMS]["prob"] = np.zeros((N_LEGS, N_LEGS), dtype=float)
        
        N_DIAGRAMS += 1

for i in range(N_DIAGRAMS):
    for j in range(N_DIAGRAMS):
        if (i != j) and (symmetry(vtx[i]["spin"], vtx[j]["spin"]) and j not in vtx[i]["sym_to"]):
            vtx[i]["sym_to"].append(j)

# Create transition table
for i in range(N_DIAGRAMS):
    state = vtx[i]["spin"]
    
    for e in range(N_LEGS):
        for x in range(N_LEGS):
            vtx[i]["new_vtx"][e, x] = -1
            
            new_state = state.copy()
            new_state[e] = - new_state[e]
            new_state[x] = - new_state[x]
            
            if new_state in allowed_vertex_state:
                new_i = allowed_vertex_state.index(new_state)
                vtx[i]["new_vtx"][e, x] = new_i

# Generate the transition probabilities
if int(sys.argv[2]) == 1:
    # Heat Bath method (Solution A)
    for i in range(N_DIAGRAMS):
        for e in range(N_LEGS):
            cum_prob = 0.0
            for x in range(N_LEGS):
                new_i = vtx[i]["new_vtx"][e, x]
                vtx[i]["prob"][e, x] = 0.0
                
                if new_i != -1:
                    denom = 0.0
                    
                    for x2 in range(N_LEGS):
                        new_i2 = vtx[i]["new_vtx"][e, x2]
                        if new_i2 != -1:
                            denom += vtx[new_i2]["H"]
                    if vtx[new_i]["H"] != 0:
                        cum_prob += vtx[new_i]["H"] / denom
                        vtx[i]["prob"][e, x] = cum_prob
else:
    # General method (Solution B)
    for i in range(N_DIAGRAMS):
        for ent in range(N_LEGS):
            A = np.zeros((N_LEGS, N_LEGS))
            b = np.zeros(N_LEGS)
            indx = vtx[i]["new_vtx"][ent, :]
            
            for e in range(N_LEGS):
                b[e] = vtx[indx[e]]["H"] if indx[e] != -1 else -1
            z = np.where(b == -1)[0][0]
            A[z, :] = -1
            A[:, z] = -1
            
            key = np.argsort(-b)
            if b[key[0]] <= b[key[1]] + b[key[2]]:
                A[key[0], key[1]] = (b[key[0]] + b[key[1]] - b[key[2]]) / 2
                A[key[0], key[2]] = (b[key[0]] - b[key[1]] + b[key[2]]) / 2
                A[key[1], key[2]] = (- b[key[0]] + b[key[1]] + b[key[2]]) / 2
                
                A[key[1], key[0]] = (b[key[0]] + b[key[1]] - b[key[2]]) / 2
                A[key[2], key[0]] = (b[key[0]] - b[key[1]] + b[key[2]]) / 2
                A[key[2], key[1]] = (- b[key[0]] + b[key[1]] + b[key[2]]) / 2
            else:
                A[key[0], key[0]] = b[key[0]] - b[key[1]] - b[key[2]]
                A[key[0], key[1]] = b[key[1]]
                A[key[0], key[2]] = b[key[2]]
                
                A[key[1], key[0]] = b[key[1]]
                A[key[2], key[0]] = b[key[2]]

            for e in range(N_LEGS):
                for x in range(N_LEGS):
                    if A[e, x] > 0.0:
                        vtx[indx[e]]["prob"][e, x] = A[e, x] / b[e]

    for i in range(N_DIAGRAMS):
        for e in range(N_LEGS):
            cum_prob = 0.0
            for x in range(N_LEGS):
                if vtx[i]["prob"][e, x] > 0.0:
                    cum_prob += vtx[i]["prob"][e, x]
                    vtx[i]["prob"][e, x] = cum_prob

with open("tmp/" + vtx_name, "w") as f:
    f.write(str(N_DIAGRAMS) + "\n")
    f.write("\n")
    
    for i in range(N_DIAGRAMS):
        f.write(str(i) + "\n")
        f.write(str(vtx[i]["type"]) + "\n")
        f.write(str(vtx[i]["H"]) + "\n")
        
        for l in range(N_LEGS):
            f.write(str(vtx[i]["spin"][l]) + " ")
        f.write("\n")
        
        for e in range(N_LEGS):
            for x in range(N_LEGS):
                f.write(str(vtx[i]["new_vtx"][e, x]) + " ")
            f.write("\n")
            
        for e in range(N_LEGS):
            for x in range(N_LEGS):
                f.write(str(vtx[i]["prob"][e, x]) + " ")
            f.write("\n")
        
        f.write("\n")

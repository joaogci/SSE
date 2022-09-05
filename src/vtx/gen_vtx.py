import numpy as np
import sys

FILENAME = str(sys.argv[1])
N_LEGS = 4
N_DIAGRAMS = 6

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
C = 0.25 * delta + h / (2 * d * J) + epsilon

vtx = list()

# Generate the vertex informatixn
for i in range(N_DIAGRAMS):
    vtx.append(dict())
    
    vtx[i]["type"] = 1 if (i == 3 or i == 4) else 0
    vtx[i]["transition_matx"] = np.zeros((N_LEGS, N_LEGS), dtype=int)
    vtx[i]["transition_prob"] = np.zeros((N_LEGS, N_LEGS), dtype=float)
    vtx[i]["spin"] = np.zeros(N_LEGS, dtype=int)
    
    if vtx[i]["type"] != 0:
        vtx[i]["H"] = 0.5
    else:
        if i == 5:
            vtx[i]["H"] = C - 0.25 * delta + h / (2 * d * J)
        elif i == 0:
            vtx[i]["H"] = C - 0.25 * delta - h / (2 * d * J)
        else:
            vtx[i]["H"] = C + 0.25 * delta
            
    # Generate transition matrix
    for e in range(N_LEGS):
        for x in range(N_LEGS):
            vtx[i]["transition_matx"][e, x] = -1
            
            if e == x:
                vtx[i]["transition_matx"][e, x] = i
            else:
                if i == 0:
                    if e == 0:
                        if x == 2:
                            vtx[i]["transition_matx"][e, x] = 2
                        elif x == 3:
                            vtx[i]["transition_matx"][e, x] = 4
                    elif e == 1:
                        if x == 3:
                            vtx[i]["transition_matx"][e, x] = 1
                        elif x == 2:
                            vtx[i]["transition_matx"][e, x] = 3
                    elif e == 2:
                        if x == 0:
                            vtx[i]["transition_matx"][e, x] = 2
                        elif x == 1:
                            vtx[i]["transition_matx"][e, x] = 3
                    elif e == 3:
                        if x == 1:
                            vtx[i]["transition_matx"][e, x] = 1
                        elif x == 0:
                            vtx[i]["transition_matx"][e, x] = 4
                elif i == 1:
                    if e == 0:
                        if x == 2:
                            vtx[i]["transition_matx"][e, x] = 5
                        elif x == 1:
                            vtx[i]["transition_matx"][e, x] = 4
                    elif e == 1:
                        if x == 3:
                            vtx[i]["transition_matx"][e, x] = 0
                        elif x == 0:
                            vtx[i]["transition_matx"][e, x] = 4
                    elif e == 2:
                        if x == 0:
                            vtx[i]["transition_matx"][e, x] = 5
                        elif x == 3:
                            vtx[i]["transition_matx"][e, x] = 3
                    elif e == 3:
                        if x == 1:
                            vtx[i]["transition_matx"][e, x] = 0
                        elif x == 2:
                            vtx[i]["transition_matx"][e, x] = 3
                elif i == 2:
                    if e == 0:
                        if x == 2:
                            vtx[i]["transition_matx"][e, x] = 0
                        elif x == 1:
                            vtx[i]["transition_matx"][e, x] = 3
                    elif e == 1:
                        if x == 3:
                            vtx[i]["transition_matx"][e, x] = 5
                        elif x == 0:
                            vtx[i]["transition_matx"][e, x] = 3
                    elif e == 2:
                        if x == 0:
                            vtx[i]["transition_matx"][e, x] = 0
                        elif x == 3:
                            vtx[i]["transition_matx"][e, x] = 4
                    elif e == 3:
                        if x == 1:
                            vtx[i]["transition_matx"][e, x] = 5
                        elif x == 2:
                            vtx[i]["transition_matx"][e, x] = 4
                elif i == 3:
                    if e == 0:
                        if x == 3:
                            vtx[i]["transition_matx"][e, x] = 5
                        elif x == 1:
                            vtx[i]["transition_matx"][e, x] = 2
                    elif e == 1:
                        if x == 2:
                            vtx[i]["transition_matx"][e, x] = 0
                        elif x == 0:
                            vtx[i]["transition_matx"][e, x] = 2
                    elif e == 2:
                        if x == 1:
                            vtx[i]["transition_matx"][e, x] = 0
                        elif x == 3:
                            vtx[i]["transition_matx"][e, x] = 1
                    elif e == 3:
                        if x == 0:
                            vtx[i]["transition_matx"][e, x] = 5
                        elif x == 2:
                            vtx[i]["transition_matx"][e, x] = 1
                elif i == 4:
                    if e == 0:
                        if x == 3:
                            vtx[i]["transition_matx"][e, x] = 0
                        elif x == 1:
                            vtx[i]["transition_matx"][e, x] = 1
                    elif e == 1:
                        if x == 2:
                            vtx[i]["transition_matx"][e, x] = 5
                        elif x == 0:
                            vtx[i]["transition_matx"][e, x] = 1
                    elif e == 2:
                        if x == 1:
                            vtx[i]["transition_matx"][e, x] = 5
                        elif x == 3:
                            vtx[i]["transition_matx"][e, x] = 2
                    elif e == 3:
                        if x == 0:
                            vtx[i]["transition_matx"][e, x] = 0
                        elif x == 2:
                            vtx[i]["transition_matx"][e, x] = 2
                elif i == 5:
                    if e == 0:
                        if x == 2:
                            vtx[i]["transition_matx"][e, x] = 1
                        elif x == 3:
                            vtx[i]["transition_matx"][e, x] = 3
                    elif e == 1:
                        if x == 3:
                            vtx[i]["transition_matx"][e, x] = 2
                        elif x == 2:
                            vtx[i]["transition_matx"][e, x] = 4
                    elif e == 2:
                        if x == 0:
                            vtx[i]["transition_matx"][e, x] = 1
                        elif x == 1:
                            vtx[i]["transition_matx"][e, x] = 4
                    elif e == 3:
                        if x == 1:
                            vtx[i]["transition_matx"][e, x] = 2
                        elif x == 0:
                            vtx[i]["transition_matx"][e, x] = 3
    
    for l in range(N_LEGS):
        if i == 0:
            vtx[i]["spin"][l] = -1
        elif i == 1:
            vtx[i]["spin"][l] = 1 if l % 2 != 0 else -1
        elif i == 2:
            vtx[i]["spin"][l] = -1 if l % 2 != 0 else 1
        elif i == 3:
            vtx[i]["spin"][l] = -1 if (l == 0 or l == 3) else 1
        elif i == 4:
            vtx[i]["spin"][l] = 1 if (l == 0 or l == 3) else -1
        elif i == 5:
            vtx[i]["spin"][l] = 1

for i in range(N_DIAGRAMS):
    # Heat Bath algorithm with cumulaative probabilites 
    for e in range(N_LEGS):
        cum_prob = 0.0
        for x in range(N_LEGS):
            new_i = vtx[i]["transition_matx"][e, x]
            vtx[i]["transition_prob"][e, x] = 0.0
            
            if new_i != -1:
                denom = 0.0
                
                for x2 in range(N_LEGS):
                    new_i2 = vtx[i]["transition_matx"][e, x2]
                    if new_i2 != -1:
                        denom += vtx[new_i2]["H"]
                if vtx[new_i]["H"] != 0:
                    cum_prob += vtx[new_i]["H"] / denom
                    vtx[i]["transition_prob"][e, x] = cum_prob

with open("tmp/vtx.txt", "w") as f:
    for i in range(N_DIAGRAMS):
        f.write(str(i) + "\n")
        f.write(str(vtx[i]["type"]) + "\n")
        f.write(str(vtx[i]["H"]) + "\n")
        
        for l in range(N_LEGS):
            f.write(str(vtx[i]["spin"][l]) + " ")
        f.write("\n")
        
        for e in range(N_LEGS):
            for x in range(N_LEGS):
                f.write(str(vtx[i]["transition_matx"][e, x]) + " ")
            f.write("\n")
            
        for e in range(N_LEGS):
            for x in range(N_LEGS):
                f.write(str(vtx[i]["transition_prob"][e, x]) + " ")
            f.write("\n")
        
        f.write("\n")


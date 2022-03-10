import numpy as np
import numpy.random as npr

# -- System Parameters --
# Number of spins
N = 8
# Number of Bonds
Nb = N
# Cut-off in the series expansion
M = 12

# -- Simulation Parameters --
T = 1.05
beta = 1 / T
rng = npr.default_rng(npr.MT19937(seed=0))


# spin = np.zeros(N)
# for i in range(N):
#     if rng.random() <= 0.5:
#         spin[i] = 1
#     else:
#         spin[i] = - 1

spin = np.array([1, 1, -1, -1, 1, -1, 1, -1])
print(f"Spin State = {spin}")

bond_site = np.zeros((2, N), dtype=np.int64)
for b in range(N):
    bond_site[0, b] = b
    bond_site[1, b] = np.mod(b + 1, N) 
# print(f"Bond Site = {bond_site[0, :]}")
# print(f"            {bond_site[1, :]}")

n = 8
opstring = np.array([4, 0, 9, 13, 6, 0, 0, 4, 13, 0, 9, 14], dtype=np.int64)
print(f"n = {n}")
print("p | a[p]  |  b[p]  |  opstring[p] ")
for p in range(M):
    a = 0
    b = 0
    if opstring[p] != 0:
        a = np.mod(opstring[p], 2) + 1
        b = opstring[p] // 2
    print(f"{p} |   {a}   |   {b}    |    {opstring[p]}")

vertex_list = np.zeros(4 * M, dtype=np.int64)
# v = 4*p + l
# l = np.mod(v, legs)
# p = 1 + (v - 1)//4

vertex_init = np.zeros(N, dtype=np.int64) - 1
vertex_last = np.zeros(N, dtype=np.int64) - 1

for p in range(M):
    if opstring[p] == 0:
        vertex_list[4 * p:4 * p + 4] = -1
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

print(f"Vertex list = ")
for p in range(M):
    print(f"[{vertex_list[4 * p]}, {vertex_list[4 * p + 1]}, {vertex_list[4 * p + 2]}, {vertex_list[4 * p + 3]}]")

def P_ins(n):
    return np.min([1, Nb * beta / (2 * (M - n))])

def P_rem(n):
    return np.min([1, 2 * (M - n + 1) / (beta * Nb)])


# diag_update()
# print()
# print("After diag_update()")
# print(f"Spin State = {spin}")
# print(f"n = {n}")
# print("p | a[p]  |  b[p]  |  opstring[p] ")
# for p in range(M):
#     a = 0
#     b = 0
#     if opstring[p] != 0:
#         a = np.mod(opstring[p], 2) + 1
#         b = opstring[p] // 2
#     print(f"{p} |   {a}   |   {b}    |    {opstring[p]}")
    
# diag_update()
# print()
# print("After second diag_update()")
# print(f"Spin State = {spin}")
# print(f"n = {n}")
# print("p | a[p]  |  b[p]  |  opstring[p] ")
# for p in range(M):
#     a = 0
#     b = 0
#     if opstring[p] != 0:
#         a = np.mod(opstring[p], 2) + 1
#         b = opstring[p] // 2
#     print(f"{p} |   {a}   |   {b}    |    {opstring[p]}")

def loop_update():
    for v0 in range(0, 4 * M, 2):
        # print(v0)
        if vertex_list[v0] < 0:
            continue
        v_in = v0
        if rng.random() < 0.5:
            # print("if", v0)
            # Transvese the loop, for all v in loop, set X[v]=-1
            while True:
                # print(v_in)
                opstring[v_in//4] = opstring[v_in//4] ^ 1
                vertex_list[v_in] = -2
                v_out = v_in ^ 1
                v_in = vertex_list[v_out]
                vertex_list[v_out] = -2
                
                if v_in == v0:
                    break
            
            # opstring[v_in//4] = opstring[v_in//4] ^ 1
            # vertex_list[v_in] = -2
            # v_out = v_in ^ 1
            # v_in = vertex_list[v_out]
            # vertex_list[v_out] = -2
            # while v_in != v0:
            #     opstring[v_in//4] = opstring[v_in//4] ^ 1
            #     vertex_list[v_in] = -2
            #     v_out = v_in ^ 1
            #     v_in = vertex_list[v_out]
            #     vertex_list[v_out] = -2
        else:
            # Transverse the loop, for all v in loop, set X[v]=-2
            # Flip the loop (change the operator types opstring[p=v/4] while loop is traversed)
            # print("else", v0)
            
            while True:
                vertex_list[v_in] = -1
                v_out = v_in ^ 1
                v_in = vertex_list[v_out]
                vertex_list[v_out] = -1
                
                if v_in == v0:
                    break
            
            # vertex_list[v_in] = -1
            # v_out = v_in ^ 1
            # v_in = vertex_list[v_out]
            # vertex_list[v_out] = -1
            # while v_in != v0:
            #     vertex_list[v_in] = -1
            #     v_out = v_in ^ 1
            #     v_in = vertex_list[v_out]
            #     vertex_list[v_out] = -1
    
    for i in range(N):
        if vertex_init[i] != -1:
            if vertex_list[vertex_init[i]] == -2:
                spin[i] = - spin[i]
        else:
            if rng.random() < 0.5:
                spin[i] = - spin[i]

loop_update()
print()
print("After loop_update()")
print(f"Spin State = {spin}")
print(f"n = {n}")
print("p | a[p]  |  b[p]  |  opstring[p] ")
for p in range(M):
    a = 0
    b = 0
    if opstring[p] != 0:
        a = np.mod(opstring[p], 2) + 1
        b = opstring[p] // 2
    print(f"{p} |   {a}   |   {b}    |    {opstring[p]}")

print(f"Vertex list = ")
for p in range(M):
    print(f"[{vertex_list[4 * p]}, {vertex_list[4 * p + 1]}, {vertex_list[4 * p + 2]}, {vertex_list[4 * p + 3]}]")




import numpy as np

def sse(beta, N, Nb, therm_cycles, mc_cycles, rng):
    # -- Init Arrays -- 
    M = np.max([4, N // 4])
    n = 0
    opstring = np.zeros(M, dtype=np.int64)

    bond_site = np.zeros((2, N), dtype=np.int64)
    spin = np.zeros(N)

    # -- Main program --
    # Init lattice 
    for b in range(N):
        bond_site[0, b] = b
        bond_site[1, b] = np.mod(b + 1, N) 

    # Init spins 
    for i in range(N):
        spin[i] = 1 if rng.random() <= 0.5 else -1

    for t in range(therm_cycles):
        n = np.count_nonzero(opstring)
        
        # diag_update()
        for p in range(M):
            if opstring[p] == 0:
                b = rng.integers(0, Nb)
                if spin[bond_site[0, b]] == spin[bond_site[1, b]]:
                    continue
                
                if rng.random() < np.min([1, Nb * beta / (2 * (M - n))]):
                    opstring[p] = 2 * b
                    n += 1
            elif np.mod(opstring[p], 2) == 0:
                if rng.random() < np.min([1, 2 * (M - n + 1) / (beta * Nb)]):
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
            if rng.random() < 0.5:
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
                if rng.random() < 0.5:
                    spin[i] = - spin[i]
        
        M_new = int(n * 4 / 3)
        if M_new > M:
            cp_opstring = opstring.copy()
            opstring = np.zeros(M_new, dtype=np.int64)
            opstring[:M] = cp_opstring
            
            M = M_new
            
    n_vals = np.zeros(mc_cycles)
    
    for t in range(mc_cycles):
        n = np.count_nonzero(opstring)
        
        # diag_update()
        for p in range(M):
            if opstring[p] == 0:
                b = rng.integers(0, Nb)
                if spin[bond_site[0, b]] == spin[bond_site[1, b]]:
                    continue
                
                if rng.random() < np.min([1, Nb * beta / (2 * (M - n))]):
                    opstring[p] = 2 * b
                    n += 1
            elif np.mod(opstring[p], 2) == 0:
                if rng.random() < np.min([1, 2 * (M - n + 1) / (beta * Nb)]):
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
            if rng.random() < 0.5:
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
                if rng.random() < 0.5:
                    spin[i] = - spin[i]
        
        n = np.count_nonzero(opstring)
        n_vals[t] = n
        
    return n_vals


import numpy as np

def read_sse_output(filename):
    """
        reads SSE simulation outputs
        
        parameters:
            filename - path for the output file
            
        output:
            dictionary with all of the sampled quantitites
            dictionary with simulation info
    """
    sim_info = dict()
    sampled = dict()
    
    with open(filename, "r") as file:
        # read simulation info
        file.readline() 
        line = file.readline().strip().split(',')

        sim_info["d"] = int(line[0])
        sim_info["L"] = int(line[1])
        sim_info["boundary_cond"] = line[2]
        sim_info["S"] = float(line[3])
        sim_info["delta"] = float(line[4])
        sim_info["h"] = float(line[5])
        sim_info["epsilon"] = float(line[6])
        
        file.readline()
        line = file.readline().strip().split(',')
        sim_info["therm_cycles"]= int(line[0])
        sim_info["mc_cycles"] = int(line[1])
        sim_info["n_bins"] = int(line[2])
        sim_info["cond"] = line[3]
        
        file.readline()
        line = file.readline().strip().split(',')
        sim_info["cpu_time"] = float(line[0])
        sim_info["n_threads"] = int(line[1])
        
        file.readline()
        line = file.readline().strip().split(',')
        sim_info["n_betas"] = int(line[0])
        sim_info["n_k"] = int(line[1])
        sim_info["x"] = int(line[2])
        sim_info["y"] = int(line[3])
        
        # read sampled quantities
        sampled["beta"] = np.zeros(sim_info["n_betas"])
        sampled["T"] = np.zeros(sim_info["n_betas"])
        sampled["n"] = np.zeros(sim_info["n_betas"])
        sampled["n_std"] = np.zeros(sim_info["n_betas"])
        sampled["E"] = np.zeros(sim_info["n_betas"])
        sampled["E_std"] = np.zeros(sim_info["n_betas"])
        sampled["C"] = np.zeros(sim_info["n_betas"])
        sampled["C_std"] = np.zeros(sim_info["n_betas"])
        sampled["m"] = np.zeros(sim_info["n_betas"])
        sampled["m_std"] = np.zeros(sim_info["n_betas"])
        sampled["m2"] = np.zeros(sim_info["n_betas"])
        sampled["m2_std"] = np.zeros(sim_info["n_betas"])
        sampled["m4"] = np.zeros(sim_info["n_betas"])
        sampled["m4_std"] = np.zeros(sim_info["n_betas"])
        sampled["ms"] = np.zeros(sim_info["n_betas"])
        sampled["ms_std"] = np.zeros(sim_info["n_betas"])
        sampled["m2s"] = np.zeros(sim_info["n_betas"])
        sampled["m2s_std"] = np.zeros(sim_info["n_betas"])
        sampled["m4s"] = np.zeros(sim_info["n_betas"])
        sampled["m4s_std"] = np.zeros(sim_info["n_betas"])
        sampled["m_sus"] = np.zeros(sim_info["n_betas"])
        sampled["m_sus_std"] = np.zeros(sim_info["n_betas"])
        sampled["S_mean"] = np.zeros(sim_info["n_betas"])
        sampled["S_std"] = np.zeros(sim_info["n_betas"])

        file.readline()
        for j in range(sim_info["n_betas"]):
            sampled["beta"][j], sampled["n"][j], _, sampled["n_std"][j], \
            sampled["E"][j], sampled["E_std"][j], sampled["C"][j], sampled["C_std"][j], \
            sampled["m"][j], sampled["m_std"][j], sampled["m2"][j], sampled["m2_std"][j], \
            sampled["m4"][j], sampled["m4_std"][j], sampled["ms"][j], sampled["ms_std"][j], \
            sampled["m2s"][j], sampled["m2s_std"][j], sampled["m4s"][j], sampled["m4s_std"][j], \
            sampled["m_sus"][j], sampled["m_sus_std"][j], sampled["S_mean"][j], sampled["S_std"][j] \
            = [float(x) for x in file.readline().strip().split(',')]
            
            sampled["T"][j] = 1.0 / sampled["beta"][j]
        
        # read equal time spin-spin correlation function
        sampled["corr_mean"] = np.zeros((sim_info["n_betas"], sim_info["L"]))
        sampled["corr_std"] = np.zeros((sim_info["n_betas"], sim_info["L"]))
        
        for j in range(sim_info["n_betas"]):
            file.readline()
            file.readline()
            file.readline()
            
            for i in range(sim_info["L"]):
                sampled["corr_mean"][j, i], sampled["corr_std"][j, i] = [float(x) for x in file.readline().strip().split(',')]

        # spin conductivity
        if sim_info["n_k"] != 0 and sim_info["cond"] != ".":
            sampled["w_k"] = np.zeros((sim_info["n_betas"], sim_info["n_k"]))
            sampled["g_spin_mean"] = np.zeros((sim_info["n_betas"], sim_info["n_k"]))
            sampled["g_spin_std"] = np.zeros((sim_info["n_betas"], sim_info["n_k"]))
            sampled["g_heat_mean"] = np.zeros((sim_info["n_betas"], sim_info["n_k"]))
            sampled["g_heat_std"] = np.zeros((sim_info["n_betas"], sim_info["n_k"]))

            if sim_info["cond"] == "both" or sim_info["cond"] == "spin":
                for j in range(sim_info["n_betas"]):
                    file.readline()
                    file.readline()
                    file.readline()
                    
                    for i in range(sim_info["n_k"]):
                        sampled["w_k"][j, i], sampled["g_spin_mean"][j, i], sampled["g_spin_std"][j, i] = [float(x) for x in file.readline().strip().split(',')]

            if sim_info["cond"] == "both" or sim_info["cond"] == "heat": 
                for j in range(sim_info["n_betas"]):
                    file.readline()
                    file.readline()
                    file.readline()
                    
                    for i in range(sim_info["n_k"]):
                        sampled["w_k"][j, i], sampled["g_heat_mean"][j, i], sampled["g_heat_std"][j, i] = [float(x) for x in file.readline().strip().split(',')]


    return sim_info, sampled

def read_exact_output(filename):
    """
        reads SSE simulation outputs
        
        parameters:
            filename - path for the output file
            
        output:
            dictionary with all of the sampled quantitites
            dictionary with simulation info
    """
    sim_info = dict()
    sampled = dict()
    
    with open(filename, "r") as file:
        # read simulation info
        file.readline() 
        line = file.readline().strip().split(',')

        sim_info["L"] = int(line[0])
        sim_info["boundary_cond"] = line[1]
        sim_info["S"] = float(line[2])
        sim_info["delta"] = float(line[3])
        sim_info["h"] = float(line[4])
        
        file.readline()
        line = file.readline().strip().split(',')
        sim_info["n_betas"] = int(line[0])
        sim_info["n_betas_k"] = int(line[1])
        sim_info["n_k"] = int(line[2])
        sim_info["x"] = int(line[3])
        sim_info["y"] = int(line[4])
        
        # read sampled quantities
        sampled["beta"] = np.zeros(sim_info["n_betas"])
        sampled["T"] = np.zeros(sim_info["n_betas"])
        sampled["E"] = np.zeros(sim_info["n_betas"])
        sampled["C"] = np.zeros(sim_info["n_betas"])
        sampled["m"] = np.zeros(sim_info["n_betas"])
        sampled["m2"] = np.zeros(sim_info["n_betas"])
        sampled["ms"] = np.zeros(sim_info["n_betas"])
        sampled["m2s"] = np.zeros(sim_info["n_betas"])
        sampled["m_sus"] = np.zeros(sim_info["n_betas"])
        
        file.readline()
        for j in range(sim_info["n_betas"]):
            sampled["beta"][j], sampled["E"][j], sampled["C"][j], \
            sampled["m"][j], sampled["m2"][j], sampled["ms"][j], \
            sampled["m2s"][j], sampled["m_sus"][j] \
            = [float(x) for x in file.readline().strip().split(',')]
            
            sampled["T"][j] = 1.0 / sampled["beta"][j]
        
        # spin conductivity
        if sim_info["n_k"] != 0:
            sampled["beta_k"] = np.zeros(sim_info["n_betas_k"])
            sampled["w_k"] = np.zeros((sim_info["n_betas_k"], sim_info["n_k"]))
            sampled["g_spin"] = np.zeros((sim_info["n_betas_k"], sim_info["n_k"]))
            sampled["g_heat"] = np.zeros((sim_info["n_betas_k"], sim_info["n_k"]))
        
            for j in range(sim_info["n_betas_k"]):
                file.readline()
                line = file.readline().strip().split()
                file.readline()
                sampled["beta_k"][j] = float(line[0])
                
                for i in range(sim_info["n_k"]):
                    sampled["w_k"][j, i], sampled["g_spin"][j, i] = [float(x) for x in file.readline().strip().split(',')]

            for j in range(sim_info["n_betas_k"]):
                file.readline()
                line = file.readline().strip().split()
                file.readline()
                sampled["beta_k"][j] = float(line[0])
                
                for i in range(sim_info["n_k"]):
                    sampled["w_k"][j, i], sampled["g_heat"][j, i] = [float(x) for x in file.readline().strip().split(',')]


    return sim_info, sampled

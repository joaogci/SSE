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
        sim_info["S"] = float(line[2])
        sim_info["delta"] = float(line[3])
        sim_info["h"] = float(line[4])
        sim_info["epsilon"] = float(line[5])
        
        file.readline()
        line = file.readline().strip().split(',')
        sim_info["therm_cycles"]= int(line[0])
        sim_info["mc_cycles"] = int(line[1])
        sim_info["n_bins"] = int(line[2])
        
        file.readline()
        line = file.readline().strip().split(',')
        sim_info["cpu_time"] = float(line[0])
        sim_info["n_threads"] = int(line[1])
        
        file.readline()
        line = file.readline().strip().split(',')
        sim_info["n_betas"] = int(line[0])
        
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
    
    return sim_info, sampled

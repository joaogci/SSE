# Stochastic Series Expansion - SSE

Implementation of Stochastic Series Expansion (SSE) Monte Carlo method for the spin-1/2 XXZ model. This version uses the Directed Loops [1] method for the loop update.
The implementation uses binnig for the estimation of the standard deviations of the sampled quantities. It is possible to run each bin in parallel mode, using openMP.

## Usage

To use the implementation, it is required that you have a C compiler (gcc-12, preferably), and a version of Python3 installed. You also need an uptodate version of NumPy. You can d´ so by running the following command on the terminal `pip3 install numpy`.

To run a simulation you will need to use the `run.sh` script and an input file. This input file contains all of the information about the simulation. This is the structure of the input file:
```
# d, L, J, delta, h, epsilon
1, 8, 1.0, 1.0, 0.0, 0.05

# therm_cycles, mc_cycles, n_bins
10000, 1000000, 10

# test_temps
# len_beta
10
# temp_range
0.05, 4.0
```
To use the test temperatures ($\beta = \{0.5, 1.0, 2.0, 4.0, 8.0, 16.0\}$), you will need to "uncomment" the line  `test_tmps`. Else, you will have to specify how many temperature point (`len_beta`) you want to generate, equaly spaced, between the `temp_range` variables.

After setting the simulation parameters, you will need to run it via the `run.sh` script. You can type 
```
$ ./run.sh -h
```
to get more information about the arguments. The most important ones are `-n <n_threads>` which specifies the number of threads used by openMP, `-i <input_name>` and `-o <output_name>` are used for specifing the input and output filenames. For ease of use, the extensions for the input and output files should be `.txt` and `.csv`, respectively.

The output file has the following structure:
```
beta,n,n2,n_std,E,E_std,C,C_std,m,m_std,m2,m2_std,m4,m4_std,ms,ms_std,m2s,m2s_std,m4s,m4s_std,sus,sus_std,binder,binder_std,binders,binders_std
...
```

[1] - "Directed loop updates for quantum lattice models", Olav F. Syljuåsen, 2003, Phys. Rev. E 67, 046701, https://doi.org/10.1103/PhysRevE.67.046701




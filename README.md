# Stochastic Series Expansion - SSE

Implementation of Stochastic Series Expansion (SSE) Monte Carlo method for the spin-1/2 XXZ model. This version uses the Directed Loops [1] method for the loop update.
The Hamiltonian of the simulated system is given by

$$ H = J \sum_{\langle i, j \rangle} \left[ \frac{1}{2} (S^+_i S^-_j + S^-_i S^+_j) + \Delta S^z_i S^z_j \right] - h \sum_i S^z_i $$

where $J$ is the coupling constant, $\Delta$ is the magnetic anisotropy along the $z$-direction and $h$ is an external magnetic field. In the code, $J = 1$, then if $\Delta > 0$ the system is antiferromagnetic and if $\Delta $ the system is ferromagnetic. 
The SSE method is the power series expansion of the partition function 

$$ Z = \text{Tr}\{e^{-\beta H}\} = \sum_{\alpha} \sum_{n=0}^{\infty} \frac{(-\beta)^n}{n!} \langle \alpha |H^n| \alpha \rangle $$

where $\beta = 1 / T$.

The implementation uses binnig for the estimation of the standard deviations of the sampled quantities. It is possible to run each bin in parallel mode, using openMP.

## Usage

To use and run the implementation, it is required that you have a C compiler (gcc-12, preferably), and a version of Python3 installed. To install GCC and Python3 you run the following commands if you have a Debian based OS (Ubuntu, Pop_OS, ...)
```bash
$ sudo apt install gcc python3
```
or if you have a MacOS based device with [Homebrew](https://brew.sh)
```bash
$ brew install gcc python3
```
You will also need an uptodate version of NumPy. You can do so by running the following command
```bash
$ pip3 install numpy
```

To run a simulation you will need to use the `run.sh` script and an input file. This input file contains all of the information about the simulation and simualted system. This is the structure of the input file:
```
# d, L, S, delta, h, epsilon
1, 8, 1/2, 1.0, 0.0, 0.05

# therm_cycles, mc_cycles, n_bins
10000, 1000000, 10

# test_temps
# len_beta
10
# temp_range
0.05, 4.0
```
The first line are the system parameters, `d` is the system dimension, `1` and `2`, `L` is the number of unit cells, `S`, `delta` and `h` are just the Hamiltonian parameters, and `epsilon` is a parameter for the Directed Loops method. It has to be a value larger or equal to 0. Simulations perform best for lower values of `epsilon`.

To use the test temperatures ( $\beta = \{0.5, 1.0, 2.0, 4.0, 8.0, 16.0\}$ ), you will need to "uncomment" the line  `test_tmps`. Else, you will have to specify how many temperature point (`len_beta`) you want to generate, equaly spaced, between the `temp_range` variables.

After setting the simulation parameters, you will need to run it via the `run.sh` script. You can type 
```bash
$ ./run.sh -h
```
to get more information about the arguments. The most important ones are `-n <n_threads>` which specifies the number of threads used by openMP, `-i <input_name>` and `-o <output_name>` are used for specifing the input and output filenames. For ease of use, the extensions for the input and output files should be `.txt` and `.csv`, respectively. For instance, to run a simulation where the input file is called `sim_input.txt` and we want an output in the directory `outputs/` named `sim_output.csv`, with 5 threads, you will have to write
```bash
$ ./run.sh -n 5 -i input.txt -o outputs/sim_output.csv
```


The output file has the following structure:
```
beta,n,n2,n_std,E,E_std,C,C_std,m,m_std,m2,m2_std,m4,m4_std,ms,ms_std,m2s,m2s_std,m4s,m4s_std,sus,sus_std,binder,binder_std,binders,binders_std
...
```

[1] - "Directed loop updates for quantum lattice models", Olav F. Syljuåsen, 2003, Phys. Rev. E 67, 046701, https://doi.org/10.1103/PhysRevE.67.046701




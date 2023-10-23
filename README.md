# Stochastic Series Expansion - SSE

Please cite the code uising the citation file provided in the repository.

Implementation of Stochastic Series Expansion (SSE) Monte Carlo method for the spin-S XXZ model. This version uses the Directed Loops [1] method for the loop update.
The Hamiltonian of the simulated system is given by

$$ H = J \sum_{\langle i, j \rangle} \left[ \frac{1}{2} (S^+_i S^-_j + S^-_i S^+_j) + \Delta S^z_i S^z_j \right] - h \sum_i S^z_i $$

where $J$ is the coupling constant, $\Delta$ is the magnetic anisotropy along the $z$-direction and $h$ is an external magnetic field. In the code, $J = 1$, then if $\Delta > 0$ the system is antiferromagnetic and if $\Delta < 0$ the system is ferromagnetic. 
The SSE method is the power series expansion of the partition function 

$$ Z = \text{Tr}\{e^{-\beta H}\} = \sum_{\alpha} \sum_{n=0}^{\infty} \frac{(-\beta)^n}{n!} \langle \alpha |H^n| \alpha \rangle $$

where $\beta \equiv 1 / T$.

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

To run a simulation you will need to use the `run.sh` script and three input files. The first input file `read.in` contains the information about the system and MC parameters for the simulation.
```
1, 8, 0.5, 1.0, 0.0, 0.05, PBC
10000, 1000000, 10
```
The first line is `d, L, S, delta, h, epsilon, boundary_cond` and the second is `therm_cycles, mc_cycles, n_bins`. `d` is the system dimension, `1` and `2`, `L` is the number of unit cells, `S` is the quantum spin number, `delta` and `h` are just the Hamiltonian parameters, `epsilon` is a parameter for the Directed Loops method and `boundary_cond` is the boundary condition for the lattice (can be `PBC` or `OBC`). `epsilon` has to be a value larger or equal to 0. Simulations perform best for low values of `epsilon`.

The second input file is `beta.in` which has the following structure
```
2
0.5
1.0
```
The first line is the number of beta values in the file and the next lines are the beta values for the MC simulation.

The final file `matsubara.in` is only available for OBC 1-dimensional problems. It contains the information about the measurement of spin conductivity.
```
5
4, 2
```
Here the first line is the number of matsubara frequencies to be computed by the system and the second is the x site for measuring the perturbation induced in the y site, respectively.

After setting the simulation parameters, you will need to run it via the `run.sh` script. You can type 
```bash
$ ./run.sh -h
```
to get more information about the arguments. The most important ones are `-n <n_threads>` which specifies the number of threads used by openMP.
```bash
$ ./run.sh -n 5 
```

[1] - "Directed loop updates for quantum lattice models", Olav F. Syljuåsen, 2003, Phys. Rev. E 67, 046701, https://doi.org/10.1103/PhysRevE.67.046701




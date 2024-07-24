# Stochastic Series Expansion - SSE

Please cite the code uising the citation file provided in the repository. [https://doi.org/10.5281/zenodo.10034200](https://doi.org/10.5281/zenodo.11172494)

Implementation of Stochastic Series Expansion (SSE) Monte Carlo method for the spin-S XXZ model. This version uses the Directed Loops [1] method for the loop update.
The Hamiltonian of the simulated system is given by

$$ H = J \sum_{\langle i, j \rangle} \left[ \frac{1}{2} (S^+_i S^-_j + S^-_i S^+_j) + \Delta S^z_i S^z_j \right] - h \sum_i S^z_i $$

where $J$ is the coupling constant, $\Delta$ is the magnetic anisotropy along the $z$-direction and $h$ is an external magnetic field. In the code, $J = 1$, then if $\Delta > 0$ the system is antiferromagnetic and if $\Delta < 0$ the system is ferromagnetic. 
The SSE method is the power series expansion of the partition function 

$$ Z = \text{Tr}\{e^{-\beta H}\} = \sum_{\alpha} \sum_{n=0}^{\infty} \frac{(-\beta)^n}{n!} \langle \alpha |H^n| \alpha \rangle $$

where $\beta \equiv 1 / T$.

The implementation uses binnig for the estimation of the standard deviations of the sampled quantities. It is possible to run each bin in parallel mode, using openMP.

## Usage

To use and run the implementation, it is required that you have a C and fortran90 compilers and a version of Python3 installed with numpy.

To compile the code, run 
```bash
$ source build.sh
```
in the main directory. This will compile the code and set the enviroment variable `$SSE_DIR`. To run the code, create a new directory and copy the parameter file (found in `scripts_and_parameters/Start/`) to the directory and type
```bash
$ $SSE_DIR/src/main n_threads $SSE_DIR
```
This will run the a SSE simulation with the paramertes specified in the parameter file using `n_threads` threads. To analyse the observables, 
```bash
$ $SSE_DIR/src/ana obs_1 obs_2 ...
```

[1] - "Directed loop updates for quantum lattice models", Olav F. Syljuåsen, 2003, Phys. Rev. E 67, 046701, https://doi.org/10.1103/PhysRevE.67.046701




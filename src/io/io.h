#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../sampling/sampling.h"
#include "../hamiltonian/hamiltonian.h"

#define BUFFER_SIZE 256

/*
 *  reads input file for the simulation - Fortran function 
 */
void read_parameters(int* Lx, int* Ly, char* boundary_condition, double* beta, double* epsilon, long* therm_cycles, int* n_bins, long* mc_sweeps, double* S, double* J_perp, double* J_par, double* h, double* D);

/*
 *  reads vertex information for the simulation 
 */
void read_vtx_info(Vertices** vtx, int* n_diagrams);

/*
 * writes observables
 */
void write_observables(Obs_scalar* obs, int n_scal);

#endif // IO_H

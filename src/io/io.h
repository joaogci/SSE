#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../sampling/sampling.h"
#include "../hamiltonian/hamiltonian.h"
#include "../sse/sse.h"

#define BUFFER_SIZE 256

/*
 *  reads input file for the simulation - Fortran function 
 */
void read_parameters(int* Lx, int* Ly, char* boundary_condition, double* beta, double* epsilon, long* therm_cycles, int* n_bins, long* mc_sweeps, double* S, double* J_perp, double* J_par, double* h, double* D);

/*
 *  reads input file for the simulation - Fortran function 
 */
void read_parameters_hyperbolic(int* p, int* q, int* nl, double* beta, double* epsilon, long* therm_cycles, int* n_bins, long* mc_sweeps, double* S, double* J_perp, double* J_par);

/*
 *  reads vertex information for the simulation 
 */
void read_vtx_info(Vertices** vtx, int* n_diagrams);

/*
 * reads the adjacency matrix file for hyperbolic lattices
 */
void read_adjacency_matrix(int*** adj_mat, int* N, int p, int q, int nl);

/*
 * writes observables
 */
void write_observables(Obs_scalar* obs_scal, int n_scal, Obs_latt* obs_eq, int n_eq);

/*
 * Writes transport observables
 */
void write_transport_obeservables(Obs_transport* obs, int n_transp);

/*
 * writes simulation statistics to file
 */
void write_sim_info(Sim_info sim);

/*
 * Returns the number of lines of a file
 */
int num_lines(FILE* fp);

#endif // IO_H

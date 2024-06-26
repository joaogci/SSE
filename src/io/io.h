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
 * reads parameters for analysis
 */
void read_parameters_analysis(int* n_rebin);

/*
 *  reads vertex information for the simulation 
 */
void read_vtx_info(Vertices** vtx, int* n_diagrams, char* sse_path);

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
 * Writes SSE configuration
 */
void write_configuration(int t_id, int N, SSE_config* conf);

/*
 * Reads SSE configuration
 */
void read_configuration(int t_id, int N, SSE_config* conf);

/*
 * Returns the number of lines of a file
 */
int num_lines(FILE* fp);

#endif // IO_H

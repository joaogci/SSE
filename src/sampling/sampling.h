#ifndef SAMPLING_H
#define SAMPLING_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../sse/sse.h"
#include "observables.h"
#include "../hamiltonian/lattice.h"
#include "../hamiltonian/hamiltonian.h"

/*
 * Sets the obervables to sample
 */
void set_observables(Obs_scalar* obs_scalar, int n_scal, Obs_latt* obs_eq, int n_eq, Lattice* latt);

/*
 * Reset observables
 */
void reset_observables(Obs_scalar* obs_scalar, int n_scal, Obs_latt* obs_eq, int n_eq);

/*
 * Sets the transport obervables to sample
 */
void set_observables_transport(Obs_transport* obs);

/*
 * Samples the set obsevables
 */
void sample(Obs_scalar* obs_scal, int n_scal, Obs_latt* obs_eq, int n_eq, XXZ_ham* ham, SSE_config* state);

/*
 * Sample scalar observables
 */ 
void sample_obs_scalar(Obs_scalar* obs, int n_scal, XXZ_ham* ham, SSE_config* state);

/*
 * Sample equal time observables
 */
void sample_obs_eq(Obs_latt* obs, int n_eq, XXZ_ham* ham, SSE_config* state);


// /* 
//  * function: sample 
//  *  samples the currect state in the simulation
//  * 
//  *  parameters:
//  *      (int) n: bin number
//  *      (int) t_idx: temperature index
//  *      (heisenberg_system *) system: system to sample from
//  *      (sse_state *) state: simulation state
//  *      (sampled_quantities *) samples: store the samples
//  */
// void sample(int n, int t_idx, heisenberg_system *system, sse_state *state, 
//     sampled_quantities *samples, pcg32_random_t* rng);

/*
 * function: compare
 * compares two numbers for the imaginary time assignments
 */
int compare( const void* num1, const void* num2);

#endif // SAMPLING_H

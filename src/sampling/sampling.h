#ifndef SAMPLING_H
#define SAMPLING_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

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

/* 
 * Free the observable structs
 */
void free_observables(Obs_scalar* obs_scalar, int n_scal, Obs_latt* obs_eq, int n_eq);

#endif // SAMPLING_H

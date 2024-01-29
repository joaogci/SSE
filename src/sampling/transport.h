#ifndef TRANSPORT_H
#define TRANSPORT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#include "../sse/sse.h"
#include "observables.h"
#include "../hamiltonian/lattice.h"
#include "../hamiltonian/hamiltonian.h"

/*
 * Sets the transport obervables to sample
 */
void set_observables_transport(Obs_transport* obs, int n_transp, double beta, XXZ_ham* ham);

/*
 * Resets transport obervables 
 */
void reset_observables_transport(Obs_transport* obs, int n_transp);

/*
 * Samples the set transport observables
 */
void sample_transport(Obs_transport* obs, int n_transp, XXZ_ham* ham, SSE_config* state, pcg32_random_t* rng);

/* 
 * Free the observable structs
 */
void free_observables_transport(Obs_transport* obs, int n_transp);

#endif // TRANSPORT_H

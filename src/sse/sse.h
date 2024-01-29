#ifndef SSE_H
#define SSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "../hamiltonian/hamiltonian.h"
#include "../hamiltonian/lattice.h"
#include "../rng/pcg_basic.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define MAX_LOOP_SIZE (100 * state->M)

/* 
 *  information about the state of the SSE simulation
 */
typedef struct SSE_config 
{
  double beta;
  int* spin_config;

  long M;
  long n;
  int* op_string;
  int* reduced_op_string;
  int* trans_op_string;
  double* op_tau;

  long n_loops;
  long loop_size;
  int* link;
  int* first;
  int* vertex_list;
} SSE_config;

/* 
 *  information about the SSE simulation
 */
typedef struct Sim_info
{
  long therm_cycles;
  long mc_sweeps;
  int n_bins;
  int n_threads;

  double wall_time;
  double avg_n_loops;
} Sim_info;


/* 
 *  initializes the sse_state struct
 */
void init_sse_config(double beta, int N, SSE_config *state);

/* 
 *  frees allocated memory in the sse_state struct
 */
void free_sse_config(SSE_config *state);

/* 
 *  resets the SSE state
 */
void reset_sse_config(int N, double Sz, SSE_config *state);

/* 
 *  diagonal update for the SSE MCS
 */
void diag_update(XXZ_ham *system, SSE_config *state, pcg32_random_t* rng);

/* 
 *  loop update for the SSE MCS
 */
void loop_update(XXZ_ham* ham, SSE_config* state, pcg32_random_t* rng);

/* 
 *  dinamically adjusts the expansion cutoff during the thermalization
 */
void ajust_cutoff(SSE_config* state, long t);

/* 
 *  tranlates the operator string in to a doubly linked list of vertecies
 */
void create_vtx_list(XXZ_ham* ham, SSE_config* state);

/* 
 *  computes probability for diagonal updates
 */
double prob(int state1, int state2, XXZ_ham* ham);

/*
 * randomly assigns times to the operators in the reduced operator string
 */
void assign_times(SSE_config* state, pcg32_random_t* rng);

/*
 * compares two numbers for the imaginary time assignments
 */
int compare( const void* num1, const void* num2);

#endif // SSE_H

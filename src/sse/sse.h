#ifndef SSE_H
#define SSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>

#include "../vtx/vtx_type.h"
#include "../rng/xorshiro256++.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define C ( 0.25 * delta + hb )
#define hb ( h / (2 * d * J) )

/* heisenberg_system
    Holds information about the simulated system.
    (int) d - dimension
    (int) L - number of unit cells
    (int) N - number of particles (pow(L, d))
    (int) Nb - number of bonds 
    (double) J - coupling constant
    (double) delta - z-axis anisotropy strength
    (double) h - z-axis magnetic field strength
    (int *) spin - spin state (length N)
    (int **) bond - lattice information (length Nb x 2) 
    (double *) prob - pre-computed probabilities for digonal update (length 3) */
typedef struct heisenberg_system 
{
    int d;
    int L;
    int N;
    int Nb;

    double J;
    double delta;
    double h;
    double epsilon;

    int *spin;
    int **bond;
    double prob[3];
} heisenberg_system;

/* sse_state
    Information about the state of the SSE simulation
    (int *) op_string - operator string (length M) 
    (int) M - maximum expansion of the partiton function
    (int) n - number of operators in the operator string
    (int) n_loops - number of loops to construct each MCS
    (int) loop_size - size of the loops
    (int *) link - linked list for the loop construction (length 4*n)
    (int *) first - first imaginary time which the spins appear (length N)
    (int *) vtx - vertex type at each imaginery time (length n)
    (vtx_element *) vtx_type - information about each vertex type (length 6) 
    (int *) red_op_string - reduced operator string (length n)
    (int *) trans_op_string - translation between the reduced and normal operator string (length n) */
typedef struct sse_state 
{
    int *op_string;
    int M;
    int n;

    int n_loops;
    int loop_size;
    int *link;
    int *first;
    int *vtx;
    vtx_element *vtx_type;

    int* red_op_string;
    int* trans_op_string;
} sse_state;

void init_heisenberg_system(int d, int L, double J, double delta, double h, double epsilon, heisenberg_system *system);
void init_sse_state(uint64_t seed, heisenberg_system *system, sse_state *state);
void reset_sse_state(heisenberg_system *system, sse_state *state);

void diag_update(double beta, heisenberg_system *system, sse_state *state);

void loop_update(heisenberg_system *system, sse_state *state);
void ajust_cutoff(sse_state *state, bool adjust_loop);

void create_vtx_list(heisenberg_system *system, sse_state *state);
double prob(int b, struct heisenberg_system *system);

void free_memory(heisenberg_system *system, sse_state *state);

#endif // SSE_H

#ifndef SSE_H
#define SSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>

#include "../vtx/vtx_type.h"
#include "../hamiltonian/hamiltonian.h"
#include "../rng/xorshiro256++.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

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
    double C;
    H_mat_term **H;

    int *spin;
    int **bond;
} heisenberg_system;

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

    uint64_t rng_state[4];
} sse_state;

void init_heisenberg_system(int d, int L, double J, double delta, double h, double epsilon, heisenberg_system *system);
void init_sse_state(uint64_t seed, heisenberg_system *system, sse_state *state);
void reset_sse_state(heisenberg_system *system, sse_state *state);

void diag_update(double beta, heisenberg_system *system, sse_state *state);

void loop_update(heisenberg_system *system, sse_state *state);
void ajust_cutoff(sse_state *state, bool adjust_loop);

void create_vtx_list(heisenberg_system *system, sse_state *state, int *red_op_string, int *trans_op_string);
double prob(int b, struct heisenberg_system *system);

void free_memory(heisenberg_system *system, sse_state *state);

#endif // SSE_H

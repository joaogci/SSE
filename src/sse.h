#ifndef SSE_H
#define SSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "functions/vtx_type.h"
#include "functions/hamiltonian.h"
#include "rng/xorshiro256++.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

struct heisenberg_system {
    int d;
    int N;
    int Nb;

    double J;
    double delta;
    double h;
    double epsilon;
    double C;
    struct H_term **H;

    int *spin;
    int **bond;
};

struct sse_state {
    int *op_string;
    int M;
    int n;

    int *link;
    int *first;
    int *vtx;
    struct vtx_type *vtx_type;
};

void simulate_sse(double *beta_vals, int len_beta, long mc_cycles, long therm_cycles, int n_bins, struct heisenberg_system *hberg_system, struct sse_state *sse_state);
void init_heisenberg_system(int d, int N, double J, double delta, double h, double epsilon, struct heisenberg_system *hberg_system);
void init_sse_state(uint64_t seed, struct heisenberg_system *hberg_system, struct sse_state *sse_state);
void reset_sse_state(struct heisenberg_system *hberg_system, struct sse_state *sse_state);
void sample(struct heisenberg_system *hberg_system, struct sse_state *sse_state);
void diag_update(double beta, struct heisenberg_system *hberg_system, struct sse_state *sse_state);
void loop_update(struct heisenberg_system *hberg_system, struct sse_state *sse_state);
void ajust_cutoff(struct sse_state *sse_state);
void create_vtx_list(struct heisenberg_system *hberg_system, struct sse_state *sse_state, int *red_op_string, int *trans_op_string);
void free_memory(struct heisenberg_system *hberg_system, struct sse_state *sse_state);

#endif // SSE_H

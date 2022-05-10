#ifndef SSE_H
#define SSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

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
    bool freed;
};

extern struct vtx_type *vtx_type;
extern struct heisenberg_system hberg_system;
extern struct sse_state sse_state;

void simulate_sse(double *beta_vals, int len_beta, long mc_cycles, long therm_cycles, int n_bins, int n_loops);
void init_heisenberg_system(int d, int N, double J, double delta, double h, double epsilon);
void init_sse_state(uint64_t seed);
void sample();
void diag_update(double beta);
void loop_update();
void ajust_cutoff();
void create_vtx_list();
void free_memory();

#endif // SSE_H

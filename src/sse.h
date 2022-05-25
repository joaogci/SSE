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
    int L;
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

struct sampled_quantities {
    int bins;
    int betas;
    double *beta_vals;

    double **n_bins;
    double **n2_bins;
    double **E_bins;
    double **C_bins;
    double **m_bins;
    double **m2_bins;
    double **ms_bins;
    double **m2s_bins;
    double **m_sus_bins;

    double *n_mean;
    double *n2_mean;
    double *n_std;

    double *E_mean;
    double *E_std;
    double *C_mean;
    double *C_std;

    double *m_mean;
    double *m_std;
    double *m2_mean;
    double *m2_std;
    double *ms_mean;
    double *ms_std;
    double *m2s_mean;
    double *m2s_std;
    double *m_sus_mean;
    double *m_sus_std;
};

void simulate_sse(double *beta_vals, int len_beta, long mc_cycles, long therm_cycles, int n_bins, struct heisenberg_system *hberg_system, struct sse_state *sse_state, struct sampled_quantities *samples);

void init_heisenberg_system(int d, int N, double J, double delta, double h, double epsilon, struct heisenberg_system *hberg_system);
void init_sse_state(uint64_t seed, struct heisenberg_system *hberg_system, struct sse_state *sse_state);
void init_samples(double *beta_vals, int len_beta, int n_bins, struct sampled_quantities *samples);
void reset_sse_state(struct heisenberg_system *hberg_system, struct sse_state *sse_state);

void sample(int n, int t_idx, struct heisenberg_system *hberg_system, struct sse_state *sse_state, struct sampled_quantities *samples);
void normalize(int t_idx, long mc_cycles, struct heisenberg_system *hberg_system, struct sampled_quantities *samples);

void diag_update(double beta, struct heisenberg_system *hberg_system, struct sse_state *sse_state);
double prob(int b, struct heisenberg_system *hberg_system);

void loop_update(struct heisenberg_system *hberg_system, struct sse_state *sse_state);
void ajust_cutoff(struct sse_state *sse_state);

void create_vtx_list(struct heisenberg_system *hberg_system, struct sse_state *sse_state, int *red_op_string, int *trans_op_string);

void write_to_file(char *filename, struct sampled_quantities *samples);
void free_memory(struct heisenberg_system *hberg_system, struct sse_state *sse_state, struct sampled_quantities *samples);

#endif // SSE_H

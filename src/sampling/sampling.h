#ifndef SAMPLING_H
#define SAMPLING_H

#include "../sse/sse.h"

typedef struct sampled_quantities 
{
    int bins;
    int betas;
    double *beta_vals;

    double **n_bins;
    double **n2_bins;
    double **E_bins;
    double **C_bins;
    double **m_bins;
    double **m2_bins;
    double **m4_bins;
    double **ms_bins;
    double **m2s_bins;
    double **m4s_bins;
    double **m_sus_bins;
    double **binder_bins;
    double **binders_bins;

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
    double *m4_mean;
    double *m4_std;
    double *ms_mean;
    double *ms_std;
    double *m2s_mean;
    double *m2s_std;
    double *m4s_mean;
    double *m4s_std;
    double *m_sus_mean;
    double *m_sus_std;

    double *binder_mean;
    double *binder_std;
    double *binders_mean;
    double *binders_std;
} sampled_quantities;

void init_samples(double *beta_vals, int len_beta, int n_bins, sampled_quantities *samples);

void sample(int n, int t_idx, heisenberg_system *system, sse_state *state, sampled_quantities *samples);
void normalize(long mc_cycles, sampled_quantities *samples, int N, int d, double J, double delta, double h, double epsilon);

void free_samples(sampled_quantities *samples);

#endif // SAMPLING_H

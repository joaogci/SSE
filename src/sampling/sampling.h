#ifndef SAMPLING_H
#define SAMPLING_H

#include "../sse/sse.h"

/*
 * struct: sampled_quantities
 *  contains all of the sampled quantities during the simulation
 *  
 *  quantities:
 *      energy
 *      specific heat
 *      magnetization
 *      staggered magnetization
 *      magnetic uniform susceptibility
 *      binder parameter
 */
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

/* 
 * function: init_samples
 *  initializes and allocates memory for the samples quantities struct
 * 
 *  parameters:
 *      (double *) beta_vals: array of the temperatures
 *      (int) len_beta: number of temperatures
 *      (int) n_bins: number of bins
 *      (sampled_quantities *) samples: struct to inilialize 
 */
void init_samples(double *beta_vals, int len_beta, int n_bins, 
    sampled_quantities *samples);

/* 
 * function: sample 
 *  samples the currect state in the simulation
 * 
 *  parameters:
 *      (int) n: bin number
 *      (int) t_idx: temperature index
 *      (heisenberg_system *) system: system to sample from
 *      (sse_state *) state: simulation state
 *      (sampled_quantities *) samples: store the samples
 */
void sample(int n, int t_idx, heisenberg_system *system, sse_state *state, 
    sampled_quantities *samples);

/* 
 * function: normalize 
 *  normalize the samples during the binning phase
 * 
 *  parameters:
 *      (long) mc_cycles: MCS for sampling
 *      (sampled_quantities *) samples: sampled quantities
 *      (int) N: number of particles
 *      (int) d: dimension of the system
 *      (double) S: quantum spin number
 *      (double) delta: anisotropy
 *      (double) h: applied magnetic field
 *      (double) epsilon: constant to the Hamiltonian
 */
void normalize(long mc_cycles, sampled_quantities *samples, int N, int d, 
    double S, double delta, double h, double epsilon);

/*
 * function: free_samples
 *  frees the allocated memory in the sampled_quantities struct
 * 
 *  parameters:
 *      (sampled_quantities *) samples: struct to free
 */
void free_samples(sampled_quantities *samples);

#endif // SAMPLING_H

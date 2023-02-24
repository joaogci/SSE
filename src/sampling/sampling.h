#ifndef SAMPLING_H
#define SAMPLING_H

#include "../sse/sse.h"
#include <time.h>

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
    int d;
    int L;
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

    double ***corr_bins;
    double **corr_mean;
    double **corr_std;

    double **S_bins;
    double *S_mean;
    double *S_std;

    int k_max;
    int x;
    int y;
    int max_samp;
    double **w_k;
    double ***g_spin_bins;
    double **g_spin_mean;
    double **g_spin_std;
    double ***g_heat_bins;
    double **g_heat_mean;
    double **g_heat_std;
} sampled_quantities;

/* 
 * function: init_samples
 *  initializes and allocates memory for the samples quantities struct
 * 
 *  parameters:
 *      (double *) beta_vals: array of the temperatures
 *      (int) len_beta: number of temperatures
 *      (int) n_bins: number of bins
 *      (int) d: number of dimensions of the system
 *      (int) L: number of unit lattices in the system
 *      (sampled_quantities *) samples: struct to inilialize 
 *      (int) max_samp: number of samples for the conductance
 */
void init_samples(double *beta_vals, int len_beta, int n_bins, int d, int L,
    sampled_quantities *samples, int max_samp);

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
    sampled_quantities *samples, pcg32_random_t* rng);

/* 
 * function: normalize 
 *  normalize the samples during the binning phase
 * 
 *  parameters:
 *      (long) mc_cycles: MCS for sampling
 *      (sampled_quantities *) samples: sampled quantities
 *      (int) N: number of particles
 *      (int) d: dimension of the system
 *      (int) boundary_cond: boundary condition of the lattice
 *      (double) S: quantum spin number
 *      (double) delta: anisotropy
 *      (double) h: applied magnetic field
 *      (double) epsilon: constant to the Hamiltonian
 */
void normalize(long mc_cycles, sampled_quantities *samples, int N, int d, int boundary_cond, 
    double S, double delta, double h, double epsilon);

/*
 * function: free_samples
 *  frees the allocated memory in the sampled_quantities struct
 * 
 *  parameters:
 *      (sampled_quantities *) samples: struct to free
 */
void free_samples(sampled_quantities *samples);

#if defined(SPIN_CONDUCTANCE) || defined(HEAT_CONDUCTANCE)

#include "acb_hypgeom.h"

/*
 * Functions to help computing the integral and factorial
 *  prefactors in the spin conductance formula. 
 */
double prefactor_spin_cond(int m, int n, int k);

/*
 * Functions to help computing the integral and factorial
 *  prefactors in the heat conductance formula. 
 */
double prefactor_heat_cond(int q, int n, int k);

/*
 * Computes the 1F1(a, b, x) hypergeometric function for an imaginary 
 *  argument x. 
 * Uses the ARB library for the calculation. (https://arblib.org)
 * parameters:
 *      (double *) re: real part of the result
 *      (double *) im: imaginary part of the result
 *      (double) a: parameter for F
 *      (double) b: parameter for F
 *      (double) x: imaginary part of the input
 */
void hyp1f1ix(double* re, double* im, double a, double b, double x);
#endif // CONDUCTANCE

#endif // SAMPLING_H

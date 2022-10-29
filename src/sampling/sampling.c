#include "sampling.h"

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
void sample(int n, int t_idx, heisenberg_system *system, sse_state *state, sampled_quantities *samples) 
{
    double m = 0.0;
    double m2 = 0.0;
    double m4 = 0.0;
    double ms = 0.0;
    double m2s = 0.0;
    double m4s = 0.0;

    samples->n_bins[t_idx][n] += state->n;
    samples->n2_bins[t_idx][n] += state->n * state->n;

    for (int i = 0; i < system->N; i++) {
        m += system->spin[i];
        ms += pow(- 1.0, i) * system->spin[i];
    }
    m *= 0.5;
    ms *= 0.5;
    m2 += m * m;
    m4 += m * m * m * m;
    m2s += ms * ms;
    m4s += ms * ms * ms * ms;

    for (int p = 0; p < state->M; p++) {
        if (state->op_string[p] % 3 != 0) {
            int b = (state->op_string[p] / 3) - 1;
            int a = state->op_string[p] % 3;

            if (a == 1) {
                system->spin[system->bond[b][0]] += 2;
                system->spin[system->bond[b][1]] += -2;
            } else if (a == 2) {
                system->spin[system->bond[b][0]] += -2;
                system->spin[system->bond[b][1]] += 2;
            }

            ms += 2 * pow(- 1.0, system->bond[b][0]) * system->spin[system->bond[b][0]];
        }

        if (state->op_string[p] != 0) {
            m2s += ms * ms;
            m4s += ms * ms * ms * ms;
        }
    }

    int norm = state->n > 0 ? state->n : 1;
    m /= system->N;
    m2 /= system->N;
    m4 /= system->N;
    ms /= (norm * system->N);
    m2s /= (norm * system->N);
    m4s /= (norm * system->N);

    samples->m_bins[t_idx][n] += m;
    samples->m2_bins[t_idx][n] += m2;
    samples->m4_bins[t_idx][n] += m4;
    samples->ms_bins[t_idx][n] += ms;
    samples->m2s_bins[t_idx][n] += m2s;
    samples->m4s_bins[t_idx][n] += m4s;
}

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
void normalize(long mc_cycles, sampled_quantities *samples, int N, int d, double S, double delta, double h, double epsilon) 
{
    for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
        for (int n = 0; n < samples->bins; n++) {
            samples->n_bins[t_idx][n] /= mc_cycles;
            samples->n2_bins[t_idx][n] /= mc_cycles;
            samples->E_bins[t_idx][n] = - samples->n_bins[t_idx][n] / (samples->beta_vals[t_idx] * N) + C + epsilon;
            samples->C_bins[t_idx][n] = (samples->n2_bins[t_idx][n] - samples->n_bins[t_idx][n] * samples->n_bins[t_idx][n] - samples->n_bins[t_idx][n]) / N;
            
            samples->m_bins[t_idx][n] /= mc_cycles;
            samples->m2_bins[t_idx][n] /= mc_cycles;
            samples->m4_bins[t_idx][n] /= mc_cycles;
            samples->ms_bins[t_idx][n] /= mc_cycles;
            samples->m2s_bins[t_idx][n] /= mc_cycles;
            samples->m4s_bins[t_idx][n] /= mc_cycles;
            samples->m_sus_bins[t_idx][n] = samples->beta_vals[t_idx] * (samples->m2_bins[t_idx][n] - samples->m_bins[t_idx][n] * samples->m_bins[t_idx][n]);
            samples->binder_bins[t_idx][n] = 1 - (samples->m4_bins[t_idx][n]) / (3 * samples->m2_bins[t_idx][n]);
            samples->binders_bins[t_idx][n] = 1 - (samples->m4s_bins[t_idx][n]) / (3 * samples->m2s_bins[t_idx][n]);

            samples->n_mean[t_idx] += samples->n_bins[t_idx][n];
            samples->n2_mean[t_idx] += samples->n2_bins[t_idx][n];
            samples->E_mean[t_idx] += samples->E_bins[t_idx][n];
            samples->C_mean[t_idx] += samples->C_bins[t_idx][n];

            samples->m_mean[t_idx] += samples->m_bins[t_idx][n];
            samples->m2_mean[t_idx] += samples->m2_bins[t_idx][n];
            samples->m4_mean[t_idx] += samples->m4_bins[t_idx][n];
            samples->ms_mean[t_idx] += samples->ms_bins[t_idx][n];
            samples->m2s_mean[t_idx] += samples->m2s_bins[t_idx][n];
            samples->m4s_mean[t_idx] += samples->m4s_bins[t_idx][n];
            samples->m_sus_mean[t_idx] += samples->m_sus_bins[t_idx][n];
            samples->binder_mean[t_idx] += samples->binder_bins[t_idx][n];
            samples->binders_mean[t_idx] += samples->binders_bins[t_idx][n];
        }
        samples->n_mean[t_idx] /= samples->bins;
        samples->n2_mean[t_idx] /= samples->bins;
        samples->E_mean[t_idx] /= samples->bins;
        samples->C_mean[t_idx] /= samples->bins;

        samples->m_mean[t_idx] /= samples->bins;
        samples->m2_mean[t_idx] /= samples->bins;
        samples->m4_mean[t_idx] /= samples->bins;
        samples->ms_mean[t_idx] /= samples->bins;
        samples->m2s_mean[t_idx] /= samples->bins;
        samples->m4s_mean[t_idx] /= samples->bins;
        samples->m_sus_mean[t_idx] /= samples->bins;
        samples->binder_mean[t_idx] /= samples->bins;
        samples->binders_mean[t_idx] /= samples->bins;

        for (int n = 0; n < samples->bins; n++) {
            samples->n_std[t_idx] += pow(samples->n_bins[t_idx][n] - samples->n_mean[t_idx], 2.0);
            samples->E_std[t_idx] += pow(samples->E_bins[t_idx][n] - samples->E_mean[t_idx], 2.0);
            samples->C_std[t_idx] += pow(samples->C_bins[t_idx][n] - samples->C_mean[t_idx], 2.0);

            samples->m_std[t_idx] += pow(samples->m_bins[t_idx][n] - samples->m_mean[t_idx], 2.0);
            samples->m2_std[t_idx] += pow(samples->m2_bins[t_idx][n] - samples->m2_mean[t_idx], 2.0);
            samples->m4_std[t_idx] += pow(samples->m4_bins[t_idx][n] - samples->m4_mean[t_idx], 2.0);
            samples->ms_std[t_idx] += pow(samples->ms_bins[t_idx][n] - samples->ms_mean[t_idx], 2.0);
            samples->m2s_std[t_idx] += pow(samples->m2s_bins[t_idx][n] - samples->m2s_mean[t_idx], 2.0);
            samples->m4s_std[t_idx] += pow(samples->m4s_bins[t_idx][n] - samples->m4s_mean[t_idx], 2.0);
            samples->m_sus_std[t_idx] += pow(samples->m_sus_bins[t_idx][n] - samples->m_sus_mean[t_idx], 2.0);
            samples->binder_std[t_idx] += pow(samples->binder_bins[t_idx][n] - samples->binder_mean[t_idx], 2.0);
            samples->binders_std[t_idx] += pow(samples->binders_bins[t_idx][n] - samples->binders_mean[t_idx], 2.0);
        }
        samples->n_std[t_idx] = sqrt(samples->n_std[t_idx] / samples->bins);
        samples->E_std[t_idx] = sqrt(samples->E_std[t_idx] / samples->bins);
        samples->C_std[t_idx] = sqrt(samples->C_std[t_idx] / samples->bins);

        samples->m_std[t_idx] = sqrt(samples->m_std[t_idx] / samples->bins);
        samples->m2_std[t_idx] = sqrt(samples->m2_std[t_idx] / samples->bins);
        samples->m4_std[t_idx] = sqrt(samples->m4_std[t_idx] / samples->bins);
        samples->ms_std[t_idx] = sqrt(samples->ms_std[t_idx] / samples->bins);
        samples->m2s_std[t_idx] = sqrt(samples->m2s_std[t_idx] / samples->bins);
        samples->m4s_std[t_idx] = sqrt(samples->m4s_std[t_idx] / samples->bins);
        samples->m_sus_std[t_idx] = sqrt(samples->m_sus_std[t_idx] / samples->bins);
        samples->binder_std[t_idx] = sqrt(samples->binder_std[t_idx] / samples->bins);
        samples->binders_std[t_idx] = sqrt(samples->binders_std[t_idx] / samples->bins);
    }
}

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
void init_samples(double *beta_vals, int len_beta, int n_bins, struct sampled_quantities *samples) 
{
    samples->bins = n_bins;
    samples->betas = len_beta;
    samples->beta_vals = beta_vals;

    samples->n_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->n2_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->E_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->C_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m2_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m4_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->ms_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m2s_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m4s_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m_sus_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->binder_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->binders_bins = (double **) malloc(len_beta * sizeof(double *));
    for (int i = 0; i < len_beta; i++) {
        samples->n_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->n2_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->E_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->C_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m2_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m4_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->ms_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m2s_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m4s_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m_sus_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->binder_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->binders_bins[i] = (double *) malloc(n_bins * sizeof(double));

        memset(samples->n_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->n2_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->E_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->C_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m2_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m4_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->ms_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m2s_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m4s_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m_sus_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->binder_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->binders_bins[i], 0.0, n_bins * sizeof(double));
    }

    samples->n_mean = (double *) malloc(len_beta * sizeof(double));
    samples->n_std = (double *) malloc(len_beta * sizeof(double));
    samples->n2_mean = (double *) malloc(len_beta * sizeof(double));

    samples->E_mean = (double *) malloc(len_beta * sizeof(double));
    samples->E_std = (double *) malloc(len_beta * sizeof(double));
    samples->C_mean = (double *) malloc(len_beta * sizeof(double));
    samples->C_std = (double *) malloc(len_beta * sizeof(double));

    samples->m_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m_std = (double *) malloc(len_beta * sizeof(double));
    samples->m2_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m2_std = (double *) malloc(len_beta * sizeof(double));
    samples->m4_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m4_std = (double *) malloc(len_beta * sizeof(double));
    samples->ms_mean = (double *) malloc(len_beta * sizeof(double));
    samples->ms_std = (double *) malloc(len_beta * sizeof(double));
    samples->m2s_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m2s_std = (double *) malloc(len_beta * sizeof(double));
    samples->m4s_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m4s_std = (double *) malloc(len_beta * sizeof(double));
    samples->m_sus_mean = (double *) malloc(len_beta * sizeof(double));
    samples->m_sus_std = (double *) malloc(len_beta * sizeof(double));

    samples->binder_mean = (double *) malloc(len_beta * sizeof(double));
    samples->binder_std = (double *) malloc(len_beta * sizeof(double));
    samples->binders_mean = (double *) malloc(len_beta * sizeof(double));
    samples->binders_std = (double *) malloc(len_beta * sizeof(double));

    memset(samples->n_mean, 0.0, len_beta * sizeof(double));
    memset(samples->n_std, 0.0, len_beta * sizeof(double));
    memset(samples->n2_mean, 0.0, len_beta * sizeof(double));
    memset(samples->E_mean, 0.0, len_beta * sizeof(double));
    memset(samples->E_std, 0.0, len_beta * sizeof(double));
    memset(samples->C_mean, 0.0, len_beta * sizeof(double));
    memset(samples->C_std, 0.0, len_beta * sizeof(double));
    memset(samples->m_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m_std, 0.0, len_beta * sizeof(double));
    memset(samples->m2_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m2_std, 0.0, len_beta * sizeof(double));
    memset(samples->m4_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m4_std, 0.0, len_beta * sizeof(double));
    memset(samples->ms_mean, 0.0, len_beta * sizeof(double));
    memset(samples->ms_std, 0.0, len_beta * sizeof(double));
    memset(samples->m2s_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m2s_std, 0.0, len_beta * sizeof(double));
    memset(samples->m4s_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m4s_std, 0.0, len_beta * sizeof(double));
    memset(samples->m_sus_mean, 0.0, len_beta * sizeof(double));
    memset(samples->m_sus_std, 0.0, len_beta * sizeof(double));
    memset(samples->binder_mean, 0.0, len_beta * sizeof(double));
    memset(samples->binder_std, 0.0, len_beta * sizeof(double));
    memset(samples->binders_mean, 0.0, len_beta * sizeof(double));
    memset(samples->binders_std, 0.0, len_beta * sizeof(double));
}

/*
 * function: free_samples
 *  frees the allocated memory in the sampled_quantities struct
 * 
 *  parameters:
 *      (sampled_quantities *) samples: struct to free
 */
void free_samples(sampled_quantities *samples) 
{
    for (int i = 0; i < samples->betas; i++) {
        free(samples->n_bins[i]);
        free(samples->n2_bins[i]);
        free(samples->E_bins[i]);
        free(samples->C_bins[i]);
        free(samples->m_bins[i]);
        free(samples->m2_bins[i]);
        free(samples->m4_bins[i]);
        free(samples->ms_bins[i]);
        free(samples->m2s_bins[i]);
        free(samples->m4s_bins[i]);
        free(samples->m_sus_bins[i]);
        free(samples->binder_bins[i]);
        free(samples->binders_bins[i]);
    }
    free(samples->n_bins);
    free(samples->n2_bins);
    free(samples->E_bins);
    free(samples->C_bins);
    free(samples->m_bins);
    free(samples->m2_bins);
    free(samples->m4_bins);
    free(samples->ms_bins);
    free(samples->m2s_bins);
    free(samples->m4s_bins);
    free(samples->m_sus_bins);
    free(samples->binder_bins);
    free(samples->binders_bins);
    free(samples->n_mean);
    free(samples->n_std);
    free(samples->n2_mean);
    free(samples->E_mean);
    free(samples->E_std);
    free(samples->C_mean);
    free(samples->C_std);
    free(samples->m_mean);
    free(samples->m_std);
    free(samples->m2_mean);
    free(samples->m2_std);
    free(samples->m4_mean);
    free(samples->m4_std);
    free(samples->ms_mean);
    free(samples->ms_std);
    free(samples->m2s_mean);
    free(samples->m2s_std);
    free(samples->m4s_mean);
    free(samples->m4s_std);
    free(samples->m_sus_mean);
    free(samples->m_sus_std);
    free(samples->binder_mean);
    free(samples->binder_std);
    free(samples->binders_mean);
    free(samples->binders_std);
}

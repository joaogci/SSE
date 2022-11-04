#include "sampling.h"


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
}

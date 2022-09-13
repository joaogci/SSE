#ifndef AUTOCORRELATION_H
#define AUTOCORRELATION_H

#include <stdlib.h>

#include "../sse/sse.h"

typedef struct autocorrelation
{
    long therm_cycles;
    long mc_cycles;

    int betas; 
    double *beta_vals;

    double **n_series;
    double **E_series;
    double **m_series;
    double **ms_series;
} autocorrelation;

void init_autocorrelation(double *beta_vals, int len_beta, long mc_cycles, long therm_cycles, autocorrelation *corr_series);
void measure_autocorrelation(int t_idx, int t, heisenberg_system *system, sse_state *state, autocorrelation *corr_series);
void free_autocorrelation(autocorrelation *corr_series);

#endif // AUTOCORRELATION_H

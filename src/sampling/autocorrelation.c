#include "autocorrelation.h"

void measure_autocorrelation(int t_idx, int t, heisenberg_system *system, sse_state *state, autocorrelation *corr_series) 
{
    double m = 0.0;
    double ms = 0.0;

    for (int i = 0; i < system->N; i++) {
        m += system->spin[i];
        ms += pow(-1.0, i) * system->spin[i];
    }
    m += 0.5;
    ms *= 0.5;

    for (int p = 0; p < state->M; p++) {
        if (state->op_string[p] % 2 == 1) {
            int b = (state->op_string[p] / 2) - 1;
            system->spin[system->bond[b][0]] = - system->spin[system->bond[b][0]];
            system->spin[system->bond[b][1]] = - system->spin[system->bond[b][1]];

            ms += 2 * pow(- 1.0, system->bond[b][0]) * system->spin[system->bond[b][0]];
        }
    }

    int norm = state->n > 0 ? state->n : 1;
    m /= system->N;
    ms /= (norm * system->N);

    corr_series->n_series[t_idx][t] = state->n;
    corr_series->E_series[t_idx][t] = - state->n / (corr_series->beta_vals[t_idx] * system->N) + system->J * (0.25 * system->delta + system->h / (2 * system->d * system->J)) + system->epsilon;
    corr_series->m_series[t_idx][t] = m;
    corr_series->ms_series[t_idx][t] = ms;
}

void init_autocorrelation(double *beta_vals, int len_beta, long mc_cycles, long therm_cycles, autocorrelation *corr_series)
{
    corr_series->mc_cycles = mc_cycles;
    corr_series->therm_cycles = therm_cycles;
    corr_series->beta_vals = beta_vals;
    corr_series->betas = len_beta;

    corr_series->n_series = (double **) malloc(len_beta * sizeof(double *));
    corr_series->E_series = (double **) malloc(len_beta * sizeof(double *));
    corr_series->m_series = (double **) malloc(len_beta * sizeof(double *));
    corr_series->ms_series = (double **) malloc(len_beta * sizeof(double *));
    for (int i = 0; i < len_beta; i++) {
        corr_series->n_series[i] = (double *) malloc((mc_cycles + therm_cycles) * sizeof(double));
        corr_series->E_series[i] = (double *) malloc((mc_cycles + therm_cycles) * sizeof(double));
        corr_series->m_series[i] = (double *) malloc((mc_cycles + therm_cycles) * sizeof(double));
        corr_series->ms_series[i] = (double *) malloc((mc_cycles + therm_cycles) * sizeof(double));

        memset(corr_series->n_series[i], 0.0, (mc_cycles + therm_cycles) * sizeof(double));
        memset(corr_series->E_series[i], 0.0, (mc_cycles + therm_cycles) * sizeof(double));
        memset(corr_series->m_series[i], 0.0, (mc_cycles + therm_cycles) * sizeof(double));
        memset(corr_series->ms_series[i], 0.0, (mc_cycles + therm_cycles) * sizeof(double));
    }
}

void free_autocorrelation(autocorrelation *corr_series)
{
    for (int i = 0; i < corr_series->betas; i++) {
        free(corr_series->n_series[i]);
        free(corr_series->E_series[i]);
        free(corr_series->m_series[i]);
        free(corr_series->ms_series[i]);
    }

    free(corr_series->n_series);
    free(corr_series->E_series);
    free(corr_series->m_series);
    free(corr_series->ms_series);
}


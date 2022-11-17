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
    for (int i = 0; i < 4; i++) {
        s[i] = rand() * (n + 1);
    }

    double m1 = 0.0;
    double m2 = 0.0;
    double m4 = 0.0;
    double ms = 0.0;
    double m2s = 0.0;
    double m4s = 0.0;
    double corr[system->L];
    double si[system->L];
    memset(corr, 0.0, system->L * sizeof(double));
    memset(si, 0.0, system->L * sizeof(double));

    // sample the first state 
    for (int i = 0; i < system->N; i++) {
        m1 += system->spin[i] * 0.5;
        ms += pow(- 1.0, i) * system->spin[i] * 0.5;

        corr[i] += system->spin[0] * system->spin[i] * 0.25;
        si[i] += system->spin[i] * 0.5;
    }
    m2 += m1 * m1;
    m4 += m1 * m1 * m1 * m1;
    m2s += ms * ms;
    m4s += ms * ms * ms * ms;

    // reduced operator string
    int red_op_string[state->n];
    int p_red = 0;
    for (int p = 0; p < state->M; p++) {
        if (state->op_string[p] != 0) {
            red_op_string[p_red] = state->op_string[p];
            p_red++;
        }
    }

    // propagate the system to sample the staggered magnetization
    for (int p = 0; p < state->n; p++) {
        if (red_op_string[p] % 3 != 0) {
            int b = (red_op_string[p] / 3) - 1;
            int a = red_op_string[p] % 3;
            int old_state = system->spin[system->bond[b][0]];

            if (a == 1) {
                system->spin[system->bond[b][0]] += 2;
                system->spin[system->bond[b][1]] += -2;
            } else if (a == 2) {
                system->spin[system->bond[b][0]] += -2;
                system->spin[system->bond[b][1]] += 2;
            }

            ms += pow(- 1.0, system->bond[b][0]) 
            * (system->spin[system->bond[b][0]] - old_state);
        }

        m2s += ms * ms;
        m4s += ms * ms * ms * ms;

        for (int i = 0; i < system->L; i++) {
            corr[i] += system->spin[0] * system->spin[i] * 0.25;
            si[i] += system->spin[i] * 0.5;
        }
    }

    // sample to the struct
    samples->n_bins[t_idx][n] += state->n;
    samples->n2_bins[t_idx][n] += state->n * state->n; 
    samples->m_bins[t_idx][n] += m1 / system->N;
    samples->m2_bins[t_idx][n] += m2 / (system->N * system->N);
    samples->m4_bins[t_idx][n] += m4 / (system->N * system->N * system->N * system->N);
    samples->ms_bins[t_idx][n] += ms / ((state->n + 1) * system->N);
    samples->m2s_bins[t_idx][n] += m2s / ((state->n + 1) * system->N * system->N);
    samples->m4s_bins[t_idx][n] += m4s / ((state->n + 1) * system->N * system->N * system->N * system->N);
    for (int i = 0; i < system->L; i++) {
        samples->corr_bins[t_idx][n][i] += corr[i] / (state->n + 1);
        samples->S_bins[t_idx][n] += pow(- 1.0, i) * corr[i] / (system->L * (state->n + 1));
    }

    // sample the spin conductance
#ifdef SPIN_COND
    if (state->n > 1) {
        int spinsum_x[state->n];
        int spinsum_y[state->n];
        int spin_prod[state->n + 1];
        memset(spinsum_x, 0, state->n * sizeof(int));
        memset(spinsum_y, 0, state->n * sizeof(int));
        memset(spin_prod, 0.0, (state->n + 1) * sizeof(int));

        for (int p = 0; p < state->n; p++) {
            for (int i = samples->x; i < system->L; i++) {
                spinsum_x[p] += system->spin[i];
            }
            for (int i = samples->y; i < system->L; i++) {
                spinsum_y[p] += system->spin[i];
            }

            if (red_op_string[p] % 3 != 0) {
                int b = (red_op_string[p] / 3) - 1;
                int a = red_op_string[p] % 3;

                if (a == 1) {
                    system->spin[system->bond[b][0]] += 2;
                    system->spin[system->bond[b][1]] += -2;
                } else if (a == 2) {
                    system->spin[system->bond[b][0]] += -2;
                    system->spin[system->bond[b][1]] += 2;
                }
            }
        }

        for (int p = 0; p < state->n; p++) {
            int m = 0;
            for (int p_prime = p; p_prime < p + state->n + 1; p_prime++) {
                spin_prod[m] += spinsum_y[p] * spinsum_x[p_prime % state->n];
                m++;
            }
        }

        for (int k = 0; k < samples->k_max; k++) {
            for (int m = 0; m < state->n + 1; m++) {
                samples->g_spin_bins[t_idx][n][k] += samples->w_k[t_idx][k] * samples->beta_vals[t_idx] * 
                    prefactor_spin_cond(m, state->n, samples->w_k[t_idx][k], samples->beta_vals[t_idx]) * 
                    spin_prod[m] * 0.25;
            }
        }
    }
#endif // SPIN_COND
}

#ifdef SPIN_COND
/*
 * Functions to help computing the integral and factorial
 *  prefactors in the spin conductance formula. 
 * The integral is computed using Monte-Carlo integration with 
 *  mc samples. 
 * The factorial prefactor is computed using the Stirling's approximation
 */
double prefactor_spin_cond(int m, int n, double w_k, double beta) 
{
    long mc = 100000;
    double val = 0.0;

    for (int i = 0; i < mc; i++) {
        double x = next_double();
        val += cos(w_k * beta * x) * pow(x, m) * pow(1 - x, n - m);
    }

    double f_n_1 = (n - 1) * log(n - 1) + 0.5 * 
        log(M_PI * (2 * (n - 1) + 1.0/3.0)) - (n - 1);

    double f_n_m = 0;
    if ((n - m) != 0) {
        f_n_m = (n - m) * log(n - m) + 0.5 * 
            log(M_PI * (2 * (n - m) + 1.0/3.0)) - (n - m);
    }
    
    double f_m = 0;
    if (m != 0) {
        f_m = m * log(m) + 0.5 * log(M_PI * (2 * (m) + 1.0/3.0)) - (m);
    }

    return exp(f_n_1 - f_n_m - f_m) * val / mc;
}
#endif // SPIN_COND

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
            samples->E_bins[t_idx][n] = - samples->n_bins[t_idx][n] 
                / (samples->beta_vals[t_idx] * N) + (N - 1) * (C_ + epsilon) / N; /* remove the added constant to the Hamiltonian */
            samples->C_bins[t_idx][n] = (samples->n2_bins[t_idx][n] 
                - samples->n_bins[t_idx][n] * samples->n_bins[t_idx][n] 
                - samples->n_bins[t_idx][n]) / N;
            
            samples->m_bins[t_idx][n] /= mc_cycles;
            samples->m2_bins[t_idx][n] /= mc_cycles;
            samples->m4_bins[t_idx][n] /= mc_cycles;
            samples->ms_bins[t_idx][n] /= mc_cycles;
            samples->m2s_bins[t_idx][n] /= mc_cycles;
            samples->m4s_bins[t_idx][n] /= mc_cycles;
            samples->m_sus_bins[t_idx][n] = samples->beta_vals[t_idx] 
                * (samples->m2_bins[t_idx][n] - samples->m_bins[t_idx][n] 
                * samples->m_bins[t_idx][n] * N);

            for (int i = 0; i < N; i++) {
                samples->corr_bins[t_idx][n][i] /= mc_cycles;
            }
            samples->S_bins[t_idx][n] /= mc_cycles;

            for (int k = 0; k < samples->k_max; k++) {
                samples->g_spin_bins[t_idx][n][k] /= mc_cycles;
            }

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

            for (int i = 0; i < N; i++) {
                samples->corr_mean[t_idx][i] += samples->corr_bins[t_idx][n][i];
            }
            samples->S_mean[t_idx] += samples->S_bins[t_idx][n];

            for (int k = 0; k < samples->k_max; k++) {
                samples->g_spin_mean[t_idx][k] += samples->g_spin_bins[t_idx][n][k];
            }
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

        for (int i = 0; i < N; i++) {
            samples->corr_mean[t_idx][i] /= samples->bins;
        }
        samples->S_mean[t_idx] /= samples->bins;

        for (int k = 0; k < samples->k_max; k++) {
            samples->g_spin_mean[t_idx][k] /= samples->bins;
        }

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

            for (int i = 0; i < N; i++) {
                samples->corr_std[t_idx][i] += pow(samples->corr_bins[t_idx][n][i] - samples->corr_mean[t_idx][i], 2.0);
            }
            samples->S_std[t_idx] += pow(samples->S_bins[t_idx][n] - samples->S_mean[t_idx], 2.0);

            for (int k = 0; k < samples->k_max; k++) {
                samples->g_spin_std[t_idx][k] += pow(samples->g_spin_bins[t_idx][n][k] - samples->g_spin_mean[t_idx][k], 2.0);
            }
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

        for (int i = 0; i < N; i++) {
            samples->corr_std[t_idx][i] = sqrt(samples->corr_std[t_idx][i] / samples->bins);
        }
        samples->S_std[t_idx] = sqrt(samples->S_std[t_idx] / samples->bins);

        for (int k = 0; k < samples->k_max; k++) {
            samples->g_spin_std[t_idx][k] = sqrt(samples->g_spin_std[t_idx][k] / samples->bins);
        }
    }
}

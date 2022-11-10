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
    double corr[system->L];


    // sample the first state 
    for (int i = 0; i < system->N; i++) {
        m += system->spin[i] * 0.5;
        ms += pow(- 1.0, i) * system->spin[i] * 0.5;

        corr[i] = 0.0;
        corr[i] += system->spin[0] * system->spin[i] * 0.25;
    }
    m2 += m * m;
    m4 += m * m * m * m;
    m2s += ms * ms;
    m4s += ms * ms * ms * ms;

    // propagate the system to sample the staggered magnetization
    for (int p = 0; p < state->M; p++) {
        if (state->op_string[p] % 3 != 0) {
            int b = (state->op_string[p] / 3) - 1;
            int a = state->op_string[p] % 3;
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

        if (state->op_string[p] != 0) {
            m2s += ms * ms;
            m4s += ms * ms * ms * ms;

            for (int i = 0; i < system->L; i++) {
                corr[i] += system->spin[0] * system->spin[i] * 0.25;
            }
        }
    }

    // sample to the struct
    samples->n_bins[t_idx][n] += state->n;
    samples->n2_bins[t_idx][n] += state->n * state->n; 
    samples->m_bins[t_idx][n] += m / system->N;
    samples->m2_bins[t_idx][n] += m2 / system->N;
    samples->m4_bins[t_idx][n] += m4 / system->N;
    samples->ms_bins[t_idx][n] += ms / ((state->n + 1) * system->N);
    samples->m2s_bins[t_idx][n] += m2s / ((state->n + 1) * system->N);
    samples->m4s_bins[t_idx][n] += m4s / ((state->n + 1) * system->N);
    for (int i = 0; i < system->L; i++) {
        samples->corr_bins[t_idx][n][i] += corr[i] / (state->n + 1);
        samples->S_bins[t_idx][n] += pow(- 1.0, i) * corr[i] / (system->L * (state->n + 1));
    }    
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
            samples->E_bins[t_idx][n] = - samples->n_bins[t_idx][n] 
                / (samples->beta_vals[t_idx] * N) + C + epsilon; /* remove the added constant to the Hamiltonian */
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
    }
}

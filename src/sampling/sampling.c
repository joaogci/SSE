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
void sample(int n, int t_idx, heisenberg_system *system, sse_state *state, sampled_quantities *samples, pcg32_random_t* rng) 
{
    double m1 = 0.0;
    double m2 = 0.0;
    double m4 = 0.0;
    double ms_tmp = 0.0;
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
        ms_tmp += pow(- 1.0, i) * system->spin[i] * 0.5;
        // corr[i] += system->spin[0] * system->spin[i] * 0.25;
        // si[i] += system->spin[i] * 0.5;
    }
    m2 += pow(m1, 2.0);
    m4 += pow(m1, 4.0);
    ms += ms_tmp;
    m2s += pow(ms_tmp, 2.0);
    m4s += pow(ms_tmp, 4.0);

    // propagate the system to sample the staggered magnetization
    for (int p = 0; p < state->n; p++) {
        int b = (state->red_op_string[p] / 3) - 1;
        int a = state->red_op_string[p] % 3;
        if (a == 1) {
            ms_tmp += pow(-1.0, system->bond[b][0]) - pow(-1.0, system->bond[b][1]);
        } else if (a == 2) {
            ms_tmp += pow(-1.0, system->bond[b][1]) - pow(-1.0, system->bond[b][0]);
        }
        ms += ms_tmp;
        m2s += pow(ms_tmp, 2.0);
        m4s += pow(ms_tmp, 4.0);

        // for (int i = 0; i < system->L; i++) {
        //     corr[i] += system->spin[0] * system->spin[i] * 0.25;
        //     si[i] += system->spin[i] * 0.5;
        // }
    }

    int norm = state->n + 1;
    // sample to the struct
    samples->n_bins[t_idx][n] += state->n;
    samples->n2_bins[t_idx][n] += state->n * state->n; 
    samples->m_bins[t_idx][n] += m1 / system->N;
    samples->m2_bins[t_idx][n] += m2 / pow(system->N, 2.0);
    samples->m4_bins[t_idx][n] += m4 / pow(system->N, 4.0);
    samples->ms_bins[t_idx][n] += ms / (system->N * norm);
    samples->m2s_bins[t_idx][n] += m2s / (pow(system->N, 2.0) * norm);
    samples->m4s_bins[t_idx][n] += m4s / (pow(system->N, 4.0) * norm);
    // for (int i = 0; i < system->L; i++) {
    //     samples->corr_bins[t_idx][n][i] += corr[i];
    //     samples->S_bins[t_idx][n] += pow(- 1.0, i) * corr[i] / system->L;
    // }

    // sample the spin and heat conductances
#ifdef L_SS
    // New sampling of Spin Conductance
    for (int samp = 0; samp < samples->max_samp; samp++) {
        // Assign random numbers to the imaginary times
        // state->n = 10000;
        double tau_spin[state->n + 2];
        tau_spin[0] = 0.0;
        for (int i = 1; i <= state->n; i++) {
            double tmp = pcg32_double_r(rng) * samples->beta_vals[t_idx];

            int j;
            for (j = i - 1; (j >= 1 && tau_spin[j] > tmp); j--)
                tau_spin[j + 1] = tau_spin[j];
        
            tau_spin[j + 1] = tmp;
        }
        tau_spin[state->n + 1] = samples->beta_vals[t_idx];

        int sum_x = 0.0;
        for (int a = samples->x + 1; a < system->L; a++) {
            sum_x += system->spin[a];
        }

        int sum_y = 0.0;
        for (int b = samples->y + 1; b < system->L; b++) {
            sum_y += system->spin[b];
        }

        for (int k = 0; k < samples->k_max; k++) {
            double Ka[2] = {};
            double Kb[2] = {};

            for (int p = 0; p <= state->n; p++) {
                Ka[0] += (sin(samples->w_k[t_idx][k] * tau_spin[p + 1]) - sin(samples->w_k[t_idx][k] * tau_spin[p])) * 0.5 * sum_x;
                Ka[1] += (cos(samples->w_k[t_idx][k] * tau_spin[p]) - cos(samples->w_k[t_idx][k] * tau_spin[p + 1])) * 0.5 * sum_x;

                Kb[0] += (sin(samples->w_k[t_idx][k] * tau_spin[p + 1]) - sin(samples->w_k[t_idx][k] * tau_spin[p])) * 0.5 * sum_y;
                Kb[1] += (cos(samples->w_k[t_idx][k] * tau_spin[p]) - cos(samples->w_k[t_idx][k] * tau_spin[p + 1])) * 0.5 * sum_y;

                if (p < state->n) {
                    int bond = (state->red_op_string[p] / 3) - 1;
                    int type = state->red_op_string[p] % 3;
                    if (bond == samples->x && type != 0) {
                        if (type == 1) {
                            sum_x += -2;
                        } else if (type == 2) {
                            sum_x += 2;
                        }
                    }
                    if (bond == samples->y && type != 0) {
                        if (type == 1) {
                            sum_y += -2;
                        } else if (type == 2) {
                            sum_y += 2;
                        }
                    }
                }
            }

            samples->L_SS_bins[t_idx][n][k] += (Ka[0] * Kb[0] + Ka[1] * Kb[1]) / (samples->w_k[t_idx][k] * samples->beta_vals[t_idx] * samples->max_samp);
        }
    }
#endif //L_SS

#ifdef L_HH
    // New Sampling of Heat Conductance
    if (state->n >= 2) {
        for (int samp = 0; samp < samples->max_samp; samp++) {
            // Assign random numbers to the imaginary times
            double tau_heat[state->n];
            for (int i = 0; i < state->n; i++) {
                double tmp = pcg32_double_r(rng) * samples->beta_vals[t_idx];

                int j;
                for (j = i - 1; (j >= 0 && tau_heat[j] > tmp); j--)
                    tau_heat[j + 1] = tau_heat[j];

                tau_heat[j + 1] = tmp;
            }

            int count_x[state->n];
            int count_y[state->n];
            int Cab = 0;
            for (int p = 0; p < state->n; p++) {
                int bond = (state->red_op_string[p] / 3) - 1;
                
                count_x[p] = 0;
                if (bond > samples->x) {
                    count_x[p] = 1;
                }
                
                count_y[p] = 0;
                if (bond > samples->y) {
                    count_y[p] = 1;
                }
                
                if (bond > samples->x && bond > samples->y) {
                    Cab++;
                }
            }

            for (int k = 0; k < samples->k_max; k++) {
                double Ka[2] = {};
                double Kb[2] = {};

                for (int p = 0; p < state->n; p++) {
                    Ka[0] += cos(samples->w_k[t_idx][k] * tau_heat[p]) * count_x[p];
                    Ka[1] += sin(samples->w_k[t_idx][k] * tau_heat[p]) * count_x[p];
                    
                    Kb[0] += cos(samples->w_k[t_idx][k] * tau_heat[p]) * count_y[p];
                    Kb[1] += sin(samples->w_k[t_idx][k] * tau_heat[p]) * count_y[p];
                }

                samples->L_HH_bins[t_idx][n][k] += samples->w_k[t_idx][k] * (Ka[0] * Kb[0] + Ka[1] * Kb[1] - Cab) / (samples->beta_vals[t_idx] * samples->max_samp);
            }
        }
    }
#endif // L_HH

#ifdef L_SH
    // Sampling of Spin-Seebeck Conductance
    if (state->n >= 1) {
        // Assign random numbers to the imaginary times
        double tau_heat[state->n + 2];
        tau_heat[0] = 0.0;
        for (int i = 1; i <= state->n; i++) {
            // double tmp = pcg32_double_r(rng) * samples->beta_vals[t_idx];
            double tmp = rand() * samples->beta_vals[t_idx] / RAND_MAX;

            int j;
            for (j = i - 1; (j >= 1 && tau_heat[j] > tmp); j--)
                tau_heat[j + 1] = tau_heat[j];
        
            tau_heat[j + 1] = tmp;
        }
        tau_heat[state->n + 1] = samples->beta_vals[t_idx];

        for (int k = 0; k < samples->k_max; k++) {
            double Ka[2] = {};
            double Kb[2] = {};

            double Ga[2] = {};
            double Gb[2] = {};

            int Ha[state->n];
            int Hb[state->n];
            double Sa[state->n + 1];
            double Sb[state->n + 1];

            for (int p = 0; p <= state->n; p++) {
                int bond = (state->red_op_string[p] / 3) - 1;

                Sa[p] = 0.0;
                for (int a = samples->x + 1; a < system->L; a++) {
                    Sa[p] += 0.5 * system->spin[a];
                }

                Sb[p] = 0.0;
                for (int a = samples->y + 1; a < system->L; a++) {
                    Sb[p] += 0.5 * system->spin[a];
                }

                Ha[p] = 0;
                if (p < state->n && bond > samples->x) {
                    Ha[p] = 1;
                }

                Hb[p] = 0;
                if (p < state->n && bond > samples->y) {
                    Hb[p] = 1;
                }

                if (p < state->n && state->red_op_string[p] % 3 != 0) {
                    int b_ = (state->red_op_string[p] / 3) - 1;
                    int a_ = state->red_op_string[p] % 3;

                    if (a_ == 1) {
                        system->spin[system->bond[b_][0]] += 2;
                        system->spin[system->bond[b_][1]] += -2;
                    } else if (a_ == 2) {
                        system->spin[system->bond[b_][0]] += -2;
                        system->spin[system->bond[b_][1]] += 2;
                    }
                }
            }

            for (int p = 0; p <= state->n; p++) {
                Ka[0] += (sin(samples->w_k[t_idx][k] * tau_heat[p + 1]) - sin(samples->w_k[t_idx][k] * tau_heat[p])) * Sa[p] / samples->w_k[t_idx][k];
                Ka[1] += (cos(samples->w_k[t_idx][k] * tau_heat[p]) - cos(samples->w_k[t_idx][k] * tau_heat[p + 1])) * Sa[p] / samples->w_k[t_idx][k];

                Kb[0] += (sin(samples->w_k[t_idx][k] * tau_heat[p + 1]) - sin(samples->w_k[t_idx][k] * tau_heat[p])) * Sb[p] / samples->w_k[t_idx][k];
                Kb[1] += (cos(samples->w_k[t_idx][k] * tau_heat[p]) - cos(samples->w_k[t_idx][k] * tau_heat[p + 1])) * Sb[p] / samples->w_k[t_idx][k];

                if (p < state->n) {
                    Ga[0] += cos(samples->w_k[t_idx][k] * tau_heat[p + 1]) * Ha[p];
                    Ga[1] += sin(samples->w_k[t_idx][k] * tau_heat[p + 1]) * Ha[p];

                    Gb[0] += cos(samples->w_k[t_idx][k] * tau_heat[p + 1]) * Hb[p];
                    Gb[1] += sin(samples->w_k[t_idx][k] * tau_heat[p + 1]) * Hb[p];
                }
            }
            
            samples->L_SH_bins[t_idx][n][k] += - samples->w_k[t_idx][k] * (Ka[0] * Gb[0] + Ka[1] * Gb[1]) / (samples->beta_vals[t_idx]);
            samples->L_HS_bins[t_idx][n][k] += - samples->w_k[t_idx][k] * (Ga[0] * Kb[0] + Ga[1] * Kb[1]) / (samples->beta_vals[t_idx]);
        }
    }
#endif // L_SH
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
 *      (int) boundary_cond: boundary condition of the lattice
 *      (double) S: quantum spin number
 *      (double) delta: anisotropy
 *      (double) h: applied magnetic field
 *      (double) epsilon: constant to the Hamiltonian
 */
void normalize(long mc_cycles, sampled_quantities *samples, int N, int d, int boundary_cond, double S, double delta, double h, double epsilon) 
{
    for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
        for (int n = 0; n < samples->bins; n++) {
            samples->n_bins[t_idx][n] /= mc_cycles;
            samples->n2_bins[t_idx][n] /= mc_cycles;
            samples->E_bins[t_idx][n] = - samples->n_bins[t_idx][n] 
                / (samples->beta_vals[t_idx] * N) + (N - boundary_cond) * (C_ + epsilon) / N; /* remove the added constant to the Hamiltonian */
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
                * N * (samples->m2_bins[t_idx][n] - samples->m_bins[t_idx][n] 
                * samples->m_bins[t_idx][n]);

            for (int i = 0; i < N; i++) {
                samples->corr_bins[t_idx][n][i] /= mc_cycles;
            }
            samples->S_bins[t_idx][n] /= mc_cycles;

            for (int k = 0; k < samples->k_max; k++) {
                samples->L_SS_bins[t_idx][n][k] /= mc_cycles;
                samples->L_HH_bins[t_idx][n][k] /= mc_cycles;
                samples->L_SH_bins[t_idx][n][k] /= mc_cycles;
                samples->L_HS_bins[t_idx][n][k] /= mc_cycles;
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
                samples->L_SS_mean[t_idx][k] += samples->L_SS_bins[t_idx][n][k];
                samples->L_HH_mean[t_idx][k] += samples->L_HH_bins[t_idx][n][k];
                samples->L_SH_mean[t_idx][k] += samples->L_SH_bins[t_idx][n][k];
                samples->L_HS_mean[t_idx][k] += samples->L_HS_bins[t_idx][n][k];
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
            samples->L_SS_mean[t_idx][k] /= samples->bins;
            samples->L_HH_mean[t_idx][k] /= samples->bins;
            samples->L_SH_mean[t_idx][k] /= samples->bins;
            samples->L_HS_mean[t_idx][k] /= samples->bins;
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
                samples->L_SS_std[t_idx][k] += pow(samples->L_SS_bins[t_idx][n][k] - samples->L_SS_mean[t_idx][k], 2.0);
                samples->L_HH_std[t_idx][k] += pow(samples->L_HH_bins[t_idx][n][k] - samples->L_HH_mean[t_idx][k], 2.0);
                samples->L_SH_std[t_idx][k] += pow(samples->L_SH_bins[t_idx][n][k] - samples->L_SH_mean[t_idx][k], 2.0);
                samples->L_HS_std[t_idx][k] += pow(samples->L_HS_bins[t_idx][n][k] - samples->L_HS_mean[t_idx][k], 2.0);
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
            samples->L_SS_std[t_idx][k] = sqrt(samples->L_SS_std[t_idx][k] / samples->bins);
            samples->L_HH_std[t_idx][k] = sqrt(samples->L_HH_std[t_idx][k] / samples->bins);
            samples->L_SH_std[t_idx][k] = sqrt(samples->L_SH_std[t_idx][k] / samples->bins);
            samples->L_HS_std[t_idx][k] = sqrt(samples->L_HS_std[t_idx][k] / samples->bins);
        }
    }
}

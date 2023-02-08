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

    // sample the spin and heat conductances
#ifdef CONDUCTANCE
    // New sampling of Spin Conductance
    for (int i = 0; i < 4; i++) {
        s[i] = rand() * (n * t_idx + 1);
    }

    // Assign random numbers to the imaginary times
    double res[samples->k_max];
    for (int k = 0; k < samples->k_max; k++) {
        res[k] = 0.0;
    }
    int max_samp = 10;

    for (int samp = 0; samp < max_samp; samp++) {
        double tau_spin[state->n + 2];
        tau_spin[0] = 0.0;
        for (int i = 1; i <= state->n; i++) {
            double tmp = next_double() * samples->beta_vals[t_idx];

            int j;
            for (j = i - 1; (j >= 1 && tau_spin[j] > tmp); j--)
                tau_spin[j + 1] = tau_spin[j];
        
            tau_spin[j + 1] = tmp;
        }
        tau_spin[state->n + 1] = samples->beta_vals[t_idx];
        
        for (int k = 0; k < samples->k_max; k++) {
            double Ka[2] = {};
            double Kb[2] = {};
            double Cab = 0.0;

            for (int p = 0; p <= state->n; p++) {
                double sum_a = 0.0;
                for (int a = samples->x; a < system->L; a++) {
                    sum_a += 0.5 * system->spin[a];
                }

                double sum_b = 0.0;
                for (int b = samples->y; b < system->L; b++) {
                    sum_b += 0.5 * system->spin[b];
                }

                Ka[0] += (cos(samples->w_k[t_idx][k] * tau_spin[p + 1] - M_PI_2) - cos(samples->w_k[t_idx][k] * tau_spin[p] - M_PI_2)) * sum_a / samples->w_k[t_idx][k];
                Ka[1] += (sin(samples->w_k[t_idx][k] * tau_spin[p + 1] - M_PI_2) - sin(samples->w_k[t_idx][k] * tau_spin[p] - M_PI_2)) * sum_a / samples->w_k[t_idx][k];

                Kb[0] += (cos(samples->w_k[t_idx][k] * tau_spin[p + 1] - M_PI_2) - cos(samples->w_k[t_idx][k] * tau_spin[p] - M_PI_2)) * sum_b / samples->w_k[t_idx][k];
                Kb[1] += (sin(samples->w_k[t_idx][k] * tau_spin[p + 1] - M_PI_2) - sin(samples->w_k[t_idx][k] * tau_spin[p] - M_PI_2)) * sum_b / samples->w_k[t_idx][k];
                
                Cab += - 2 * (1 - cos(samples->w_k[t_idx][k] * (tau_spin[p + 1] - tau_spin[p]))) * sum_b * sum_a / (samples->w_k[t_idx][k] * samples->w_k[t_idx][k]);

                if (p < state->n && red_op_string[p] % 3 != 0) {
                    int b_ = (red_op_string[p] / 3) - 1;
                    int a_ = red_op_string[p] % 3;

                    if (a_ == 1) {
                        system->spin[system->bond[b_][0]] += 2;
                        system->spin[system->bond[b_][1]] += -2;
                    } else if (a_ == 2) {
                        system->spin[system->bond[b_][0]] += -2;
                        system->spin[system->bond[b_][1]] += 2;
                    }
                }
            }
            
            res[k] += samples->w_k[t_idx][k] * (Ka[0] * Kb[0] + Ka[1] * Kb[1]) / samples->beta_vals[t_idx];
        }
    }
    for (int k = 0; k < samples->k_max; k++) {
        samples->g_spin_bins[t_idx][n][k] += res[k] / max_samp;
    }
    
    // New Sampling of Heat Conductance
    if (n >= 2) {
        double res[samples->k_max];
        for (int k = 0; k < samples->k_max; k++) {
            res[k] = 0.0;
        }
        int max_samp = 10;

        for (int samp = 0; samp < max_samp; samp++) {
            // Assign random numbers to the imaginary times
            double tau_heat[state->n];
            for (int i = 0; i < state->n; i++) {
                double tmp = next_double() * samples->beta_vals[t_idx];

                int j;
                for (j = i - 1; (j >= 0 && tau_heat[j] > tmp); j--)
                    tau_heat[j + 1] = tau_heat[j];

                tau_heat[j + 1] = tmp;
            }

            for (int k = 0; k < samples->k_max; k++) {
                double Ka[2] = {};
                double Kb[2] = {};
                double Cab = 0.0;

                for (int p = 0; p < state->n; p++) {
                    int bond = (red_op_string[p] / 3) - 1;

                    if (bond >= samples->x) {
                        Ka[0] += cos(samples->w_k[t_idx][k] * tau_heat[p]);
                        Ka[1] += sin(samples->w_k[t_idx][k] * tau_heat[p]);
                    }
                    
                    if (bond >= samples->y) {
                        Kb[0] += cos(samples->w_k[t_idx][k] * tau_heat[p]);
                        Kb[1] += sin(samples->w_k[t_idx][k] * tau_heat[p]);
                    }
                    
                    if (bond >= samples->x && bond >= samples->y) {
                        Cab++;
                    }
                }

                res[k] += samples->w_k[t_idx][k] * (Ka[0] * Kb[0] + Ka[1] * Kb[1] - Cab) / samples->beta_vals[t_idx];                    
            }
        }

        for (int k = 0; k < samples->k_max; k++) {
            samples->g_heat_bins[t_idx][n][k] += res[k] / max_samp;
        }
    }

#endif // CONDUCTANCE
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
                * (samples->m2_bins[t_idx][n] - samples->m_bins[t_idx][n] 
                * samples->m_bins[t_idx][n]);

            for (int i = 0; i < N; i++) {
                samples->corr_bins[t_idx][n][i] /= mc_cycles;
            }
            samples->S_bins[t_idx][n] /= mc_cycles;

            for (int k = 0; k < samples->k_max; k++) {
                samples->g_spin_bins[t_idx][n][k] /= mc_cycles;
                samples->g_heat_bins[t_idx][n][k] /= mc_cycles;
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
                samples->g_heat_mean[t_idx][k] += samples->g_heat_bins[t_idx][n][k];
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
            samples->g_heat_mean[t_idx][k] /= samples->bins;
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
                samples->g_heat_std[t_idx][k] += pow(samples->g_heat_bins[t_idx][n][k] - samples->g_heat_mean[t_idx][k], 2.0);
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
            samples->g_heat_std[t_idx][k] = sqrt(samples->g_heat_std[t_idx][k] / samples->bins);
        }
    }
}

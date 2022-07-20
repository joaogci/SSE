#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <omp.h>

#include "sse/sse.h"
#include "sampling/sampling.h"

int main(int argc, char **argv)
{
    int d = 1;
    int L = 2;
    double J = 1.0;
    double delta = 1.0;
    double h = 0.0;
    double epsilon = 0.05;

    uint64_t seed = (u_int64_t) time(NULL);
    long therm_cycles = 1e4;
    long mc_cycles = 1e4;
    int n_bins = 10;

    double beta_vals[] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};
    int len_beta = sizeof(beta_vals) / sizeof(beta_vals[0]);

    sampled_quantities *samples = (sampled_quantities *) malloc(sizeof(sampled_quantities));
    init_samples(beta_vals, len_beta, n_bins, samples);

    // for (int i = 0; i < 4; i++) {
    //     s[i] = seed * i;
    // }

    heisenberg_system *system = (heisenberg_system *) malloc(sizeof(heisenberg_system));
    sse_state *state = (sse_state *) malloc(sizeof(sse_state));

    init_heisenberg_system(d, L, J, delta, h, epsilon, system);
    init_sse_state(seed, system, state);

    for (int t_idx = 0; t_idx < len_beta; t_idx++) {
        double beta = beta_vals[t_idx];
        reset_sse_state(system, state);

        for (long t = 0; t < therm_cycles; t++) {
            diag_update(beta, system, state);

            for (int loop = 0; loop < state->n_loops; loop++) {
                loop_update(system, state);
            }

            ajust_cutoff(state, t % 1000 == 0);
        }

        for (int n = 0; n < n_bins; n++) {
            for (long t = 0; t < mc_cycles; t++) {
                diag_update(beta, system, state);

                // state->loop_size = 0;
                for (int loop = 0; loop < state->n_loops; loop++) {
                    loop_update(system, state);
                }

                sample(n, t_idx, system, state, samples);
            }
        }

        normalize(t_idx, mc_cycles, system, samples);

        printf("beta: %f; n_mean: %f +/- %f; E: %f +/- %f \n", beta, samples->n_mean[t_idx], samples->n_std[t_idx], samples->E_mean[t_idx], samples->E_std[t_idx]);
    }
    free_memory(system, state);
    free(state);
    free(system);

    write_to_file("1D_heisenberg_L2.csv", samples);

    free_samples(samples);
    free(samples);

    return 0;
}

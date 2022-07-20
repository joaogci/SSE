#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <omp.h>

#include "sse/sse.h"
#include "sampling/sampling.h"

#define SEED (u_int64_t) time(NULL)

int d = 1;
int L = 2;
double J = 1.0;
double delta = 1.0;
double h = 0.0;
double epsilon = 0.05;

long therm_cycles = 1e4;
long mc_cycles = 1e4;
int n_bins = 10;

double beta_vals[] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};
int len_beta = sizeof(beta_vals) / sizeof(beta_vals[0]);

sampled_quantities *samples;

void simulate(sampled_quantities *samples, int start_bin, int end_bin, int n_thread) 
{
    heisenberg_system *system = (heisenberg_system *) malloc(sizeof(heisenberg_system));
    sse_state *state = (sse_state *) malloc(sizeof(sse_state));

    init_heisenberg_system(d, L, J, delta, h, epsilon, system);
    init_sse_state(SEED * n_thread, system, state);

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

        for (int n = start_bin; n < end_bin; n++) {
            for (long t = 0; t < mc_cycles; t++) {
                diag_update(beta, system, state);

                for (int loop = 0; loop < state->n_loops; loop++) {
                    loop_update(system, state);
                }

                sample(n, t_idx, system, state, samples);
            }
        }

        #pragma omp master
        printf("beta: %f \n", beta);
    }

    free_memory(system, state);
    free(state);
    free(system);
}

int main(int argc, char **argv)
{
    samples = (sampled_quantities *) malloc(sizeof(sampled_quantities));
    init_samples(beta_vals, len_beta, n_bins, samples);

    int start, end;
    // clock_t start_clock, end_clock;
    time_t start_clock, end_clock;
    double cpu_time_used;

    omp_set_num_threads(10);

    // start_clock = clock();
    time(&start_clock);
    #pragma omp parallel shared(samples) private(start, end)
    {
        int T = omp_get_num_threads();
        int n = omp_get_thread_num();

        start = (n * n_bins) / T;
        end = ((n + 1) * n_bins / T); 

        simulate(samples, start, end, n + 1);
    }
    // end_clock = clock();
    // cpu_time_used = ((double) (end_clock - start_clock)) / CLOCKS_PER_SEC;
    time(&end_clock);
    cpu_time_used = difftime(end_clock, start_clock);
    printf("done in %.2lfs \n", cpu_time_used);

    normalize(mc_cycles, samples, pow(L, d), J, 0.25 * delta + 0.5 * h / J + epsilon);

    write_to_file("1D_heisenberg_L2.csv", samples);

    free_samples(samples);
    free(samples);

    return 0;
}

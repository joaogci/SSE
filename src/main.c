#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "sse/sse.h"
#include "sampling/sampling.h"
#include "io/io.h"

#ifdef M_AUTOCORRELATION
    #include "sampling/autocorrelation.h"
#endif

#ifndef TESTING
#define SEED (u_int64_t) time(NULL)
#else
#define SEED (u_int64_t) 2
#endif

int d;
int L;
double J;
double delta;
double h;
double epsilon;

long therm_cycles;
long mc_cycles;
int n_bins;
int n_threads;

double *beta_vals;
int len_beta;

sampled_quantities *samples;
#ifdef M_AUTOCORRELATION
    autocorrelation *corr_series;
#endif

void simulate(int start_bin, int end_bin, int t_id, char *vtx_file)
{
    heisenberg_system *system = (heisenberg_system *) malloc(sizeof(heisenberg_system));
    sse_state *state = (sse_state *) malloc(sizeof(sse_state));

    init_heisenberg_system(d, L, J, delta, h, epsilon, system);
    init_sse_state(SEED * t_id, system, state);
    read_vtx_info(vtx_file, &(state->vtx_type), &(state->n_diagrams));

    for (int t_idx = 0; t_idx < len_beta; t_idx++) {
        clock_t start_clock = clock();

        double beta = beta_vals[t_idx];
        reset_sse_state(system, state);

        for (long t = 0; t < therm_cycles; t++) {
            diag_update(beta, system, state);

            create_vtx_list(system, state);
            for (int loop = 0; loop < state->n_loops; loop++) {
                loop_update(system, state);
            }

            ajust_cutoff(state, t % 1000 == 0);
            #ifdef M_AUTOCORRELATION
                measure_autocorrelation(t_idx, t, system, state, corr_series);
            #endif
        }

        for (int n = start_bin; n < end_bin; n++) {
            for (long t = 0; t < mc_cycles / n_bins; t++) {
                diag_update(beta, system, state);

                create_vtx_list(system, state);
                for (int loop = 0; loop < state->n_loops; loop++) {
                    loop_update(system, state);
                }

                sample(n, t_idx, system, state, samples);
                #ifdef M_AUTOCORRELATION
                    measure_autocorrelation(t_idx, t + therm_cycles, system, state, corr_series);
                #endif
            }
        }

        #pragma omp barrier
        #pragma omp master 
        {
            time_t t = time(NULL);
            char *buff = ctime(&t);
            buff[strcspn(buff, "\n")] = 0;

            clock_t end_clock = clock();
            double cpu_time_used = ((double) (end_clock - start_clock)) / CLOCKS_PER_SEC;

            printf("%s | beta: %.4lf | time: %.5lfs \n", buff, beta, cpu_time_used / n_threads);
        }
    }

    free_memory(system, state);
    free(state);
    free(system);
}

int main(int argc, char **argv)
{
    #ifndef M_AUTOCORRELATION
        if (argc < 5) {
            printf("Please provide the input and outut file names for the program to work. \n");
            printf("Usage: %s n_threads input_name.txt vtx_name.txt output_name.csv", argv[0]);
            exit(1);
        }
    #else 
        if (argc < 6) {
            printf("Please provide the input, outut and correlation series file names for the program to work. \n");
            printf("Usage: %s n_threads input_name.txt vtx_name.txt output_name.csv corr_name.csv", argv[0]);
            exit(1);
        }
    #endif

    read_inputs(argv[2], &d, &L, &J, &delta, &h, &epsilon, &therm_cycles, &mc_cycles, &n_bins, &beta_vals, &len_beta);
    n_threads = atoi(argv[1]);

    if (n_bins < n_threads) {
        printf("The number of bins (%d) should be equal or larger than the number of threads (%d).\n", n_bins, n_threads);
        exit(1);
    }

    #ifdef M_AUTOCORRELATION
        if (n_threads != 1) {
            printf("Please only use one thread. \n");
            exit(1);
        }

        corr_series = (autocorrelation *) malloc(sizeof(autocorrelation));
        init_autocorrelation(beta_vals, len_beta, mc_cycles, therm_cycles, corr_series);
    #endif
    samples = (sampled_quantities *) malloc(sizeof(sampled_quantities));
    init_samples(beta_vals, len_beta, n_bins, samples);
    
    omp_set_num_threads(n_threads);

    time_t t = time(NULL);
    printf(" -- Starting SSE simulation of the Heisenberg model -- \n");
    printf("   d: %d | L: %d | J: %.2lf | delta: %.2lf | h: %.2lf | epsilon: %.2lf \n", d, L, J, delta, h, epsilon);
    printf("   n_threads: %d | therm_cycles: %ld | mc_cycles: %ld | n_bins: %d \n", n_threads, therm_cycles, mc_cycles, n_bins);
    #ifdef M_AUTOCORRELATION
        printf("    Measuring time series \n");
    #endif
    printf("   Simulation started at: %s ", ctime(&t));
    printf("\n");

    clock_t start_clock = clock();
    #pragma omp parallel shared(samples)
    {
        int team_size = omp_get_num_threads();
        int t_id = omp_get_thread_num();

        int start = (t_id * n_bins) / team_size;
        int end = ((t_id + 1) * n_bins) / team_size; 

        simulate(start, end, t_id + 1, argv[3]);
    }
    clock_t end_clock = clock();
    double cpu_time_used = ((double) (end_clock - start_clock)) / CLOCKS_PER_SEC;
    
    printf("\n");
    printf("Simulation finished in %.5lfs \n", cpu_time_used / n_threads);
    printf(" -- Writing simulation results to file -- \n");

    normalize(mc_cycles, samples, pow(L, d), d, J, delta, h, epsilon);
    write_outputs(argv[4], samples);
    
    printf(" -- Results written with success -- \n");

    #ifdef M_AUTOCORRELATION
        write_autocorrelation(argv[5], corr_series);
        printf(" -- Correlation series written with success -- \n");
    #endif    

    free_samples(samples);
    free(samples);
    #ifdef M_AUTOCORRELATION
        free_autocorrelation(corr_series);
        free(corr_series);
    #endif

    return 0;
}

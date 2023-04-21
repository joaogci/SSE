#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "sse/sse.h"
#include "sampling/sampling.h"
#include "io/io.h"

#ifndef TESTING
#define SEED (u_int64_t) time(NULL)
#else
#define SEED (u_int64_t) 2
#endif

// gobal variables for the system and simulation
int d;
int L;
int boundary_cond;
double S;
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

// variables for debugging
double *cpu_time;
int *loop_size;
int *n_loops;

void simulate(int start_bin, int end_bin, int t_id, char *vtx_file)
{
    heisenberg_system *system = (heisenberg_system *) malloc(sizeof(heisenberg_system));
    sse_state *state = (sse_state *) malloc(sizeof(sse_state));

    init_heisenberg_system(d, L, boundary_cond, S, delta, h, epsilon, system);
    init_sse_state(SEED * t_id, system, state);

    pcg32_random_t rng;
    pcg32_srandom_r(&rng, (SEED * t_id) ^ (intptr_t)&rng, (SEED * t_id));

    #pragma omp critical 
    read_vtx_info(vtx_file, &(state->vtx_type), &(state->n_diagrams));

    for (int t_idx = 0; t_idx < len_beta; t_idx++) {
        clock_t start_clock = clock();

        double beta = beta_vals[t_idx];
        if (t_idx > 0) { reset_sse_state(system, state); }

        for (long t = 0; t < therm_cycles; t++) {
            diag_update(beta, system, state, &rng);

            create_vtx_list(system, state);
            loop_update(system, state, &rng);

            ajust_cutoff(state, t);
        }

        clock_t end_clock = clock();
        printf("therm ended: %lu time: %lf \n", state->n, ((double) (end_clock - start_clock)) / CLOCKS_PER_SEC/ n_threads);

        start_clock = clock();
        for (int n = start_bin; n < end_bin; n++) {
            for (long t = 0; t < mc_cycles; t++) {
                diag_update(beta, system, state, &rng);
                
                state->loop_size = 0;
                create_vtx_list(system, state);
                loop_update(system, state, &rng);

                sample(n, t_idx, system, state, samples, &rng);
            }
        }
        end_clock = clock();
        printf("mc cycles ended: %lu time: %lf \n", state->n, ((double) (end_clock - start_clock)) / CLOCKS_PER_SEC/ n_threads);


        // clock_t end_clock = clock();

        time_t t = time(NULL);
        char *buff = ctime(&t);
        buff[strcspn(buff, "\n")] = 0;

        cpu_time[t_id - 1] = ((double) (end_clock - start_clock)) / CLOCKS_PER_SEC;
        n_loops[t_id - 1] = state->n_loops;

        #pragma omp barrier
        #pragma omp master 
        {
            double max = cpu_time[0];
            int max_id = 0;
            int avg_n_loops = 0;

            for (int i = 0; i < n_threads; i++) {
                if (max < cpu_time[i]) {
                    max = cpu_time[i];
                    max_id = i;
                }

                avg_n_loops += n_loops[i];
            }
            avg_n_loops /= n_threads;

            printf("%s | beta: %.4lf | n_loops: %d | max_time: %.5lfs (T_id: %d) \n", 
                buff, beta, avg_n_loops, max / n_threads, max_id);
            fflush(stdout);
        }
        #pragma omp barrier
    }

    free_memory(system, state);
    free(state);
    free(system);
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        printf("Please provide the output file names for the program to work. \n");
        printf("Usage: %s n_threads vtx_name.txt \n", argv[0]);
        exit(1);
    }

    read_inputs(&d, &L, &boundary_cond, &S, &delta, &h, &epsilon, &therm_cycles, 
        &mc_cycles, &n_bins, &beta_vals, &len_beta);
    n_threads = atoi(argv[1]);

    if (n_bins < n_threads) {
        printf("The number of bins (%d) should be equal or larger than" 
            " the number of threads (%d).\n", n_bins, n_threads);
        exit(1);
    }
    if (boundary_cond == 1 && d != 1) {
        printf("OBC can only be used in a one-dimensional system. \n");
        exit(1);
    }
#if KINETIC
    if (boundary_cond != 1 || d != 1) {
        printf("Conductances can only be computed in a one-dimensional system"
                " with OBC. \n");
        exit(1);
    }
#endif // KINETIC

    samples = (sampled_quantities *) malloc(sizeof(sampled_quantities));
    int max_samp = 1;
    if (argc == 4) {
        max_samp = atoi(argv[3]);
    }
    init_samples(beta_vals, len_beta, n_bins, d, L, samples, max_samp);
    
    omp_set_num_threads(n_threads);
    cpu_time = (double *) malloc(sizeof(double) * n_threads);
    loop_size = (int *) malloc(sizeof(int) * n_threads);
    n_loops = (int *) malloc(sizeof(int) * n_threads);

    time_t t = time(NULL);
    printf(" -- Starting SSE simulation of the spin-S XXZ model -- \n");
    printf("   d: %d | L: %d | %s | S: %.1lf | delta: %.2lf | h: %.2lf | epsilon: %.2lf \n", 
        d, L, (boundary_cond == 0) ? "PBC" : "OBC" , S, delta, h, epsilon);
    printf("   n_threads: %d | therm_cycles: %ld | mc_cycles: %ld | n_bins: %d \n", 
        n_threads, therm_cycles, mc_cycles, n_bins);
#ifdef KINETIC
    printf("   Computing kinetic coefficients | Max Matsubara Freq: %d | x: %d | y: %d | samples: %d \n", samples->k_max, samples->x, samples->y, samples->max_samp);
#endif // KINETIC
    printf("   Simulation started at: %s ", ctime(&t));
    printf("\n");
    fflush(stdout);

    clock_t start_clock = clock();
    #pragma omp parallel shared(samples, cpu_time, loop_size, n_loops)
    {
        int team_size = omp_get_num_threads();
        int t_id = omp_get_thread_num();

        int start = (t_id * n_bins) / team_size;
        int end = ((t_id + 1) * n_bins) / team_size; 

        simulate(start, end, t_id + 1, argv[2]);
    }
    clock_t end_clock = clock();
    double cpu_time_used = ((double) (end_clock - start_clock)) / (CLOCKS_PER_SEC * n_threads);
    free(cpu_time);
    free(loop_size);
    free(n_loops);
    
    printf("\n");
    printf("Simulation finished in %.5lfs \n", cpu_time_used);
    printf(" -- Writing simulation results to file -- \n");
    fflush(stdout);

    normalize(mc_cycles, samples, pow(L, d), d, boundary_cond, S, delta, h, epsilon);
    char *file_name = write_outputs(samples, d, L, boundary_cond, S, delta, h, epsilon,
        therm_cycles, mc_cycles, cpu_time_used, n_threads);
    
    printf(" -- Results written with success to file: %s -- \n", file_name);
    fflush(stdout);

    free_samples(samples);
    free(samples);
    free(file_name);

    return 0;
}

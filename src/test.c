#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "sse.h"

double mean(int *arr, int size) {
    double mean = 0;
    int i;
    for (i = 0; i < size; i++) {
        mean += arr[i];
    }
    return mean / size;
}

int main(int argc, char **argv) {
    int N = 2;
    double beta;
    double beta_vals[6] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};

    long therm_cycles = 1e5;
    long mc_cycles = 1e6;

    struct sse_state *sse_state = malloc(sizeof(struct sse_state));
    struct vertex *vertex = malloc(sizeof(struct vertex));

    int t;
    int n_vals[mc_cycles];
    int n2_vals[mc_cycles];

    for (int i = 0; i < 6; i++) {
        beta = beta_vals[i];

        init_system(N, beta, sse_state, vertex, (uint64_t) time(NULL));

        for (t = 0; t < therm_cycles; t++) {
            diag_update(sse_state);
            create_vertex_list(sse_state, vertex);
            loop_update(sse_state, vertex);
            ajust_cutoff(sse_state, vertex);
        }

        for (t = 0; t < mc_cycles; t++) {
            diag_update(sse_state);
            create_vertex_list(sse_state, vertex);
            loop_update(sse_state, vertex);

            n_vals[t] = sse_state->n;
            n2_vals[t] = sse_state->n * sse_state->n;
        }

        double n_mean = mean(n_vals, mc_cycles);
        double n2_mean = mean(n2_vals, mc_cycles);
        double E_mean = - n_mean / (N * beta) + 0.25;
        double C_mean = (n2_mean - n_mean * (n_mean - 1)) / N;

        printf("beta: %f | E_mean: %f | C_mean: %f | n_mean %f \n", beta, E_mean, C_mean, n_mean);
    }

    free_state(&sse_state);
    free_vertex(&vertex);

    return 0;
}
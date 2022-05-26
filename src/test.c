#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include "sse.h"

int main(int argc, char **argv) {
    int d = 1;
    int N = 2;
    double J = 1.0;
    double delta = 1.0;
    double h = 0.0;
    double epsilon = 0.5;

    uint64_t seed = (u_int64_t) time(NULL);
    long therm_cycles = 1e4;
    long mc_cycles = 1e4;
    int n_bins = 10;

    double beta_vals[] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};
    int len_beta = sizeof(beta_vals) / sizeof(beta_vals[0]);

    struct heisenberg_system *hberg_system = (struct heisenberg_system *) malloc(sizeof(struct heisenberg_system));
    struct sse_state *sse_state = (struct sse_state *) malloc(sizeof(struct sse_state));    
    struct sampled_quantities *samples = (struct sampled_quantities *) malloc(sizeof(struct sampled_quantities));

    init_heisenberg_system(d, N, J, delta, h, epsilon, hberg_system);
    init_sse_state(seed, hberg_system, sse_state);

    simulate_sse(beta_vals, len_beta, mc_cycles, therm_cycles, n_bins, hberg_system, sse_state, samples);

    write_to_file("1D_heisenberg_N2.csv", samples);

    free_memory(hberg_system, sse_state, samples);
    free(hberg_system);
    free(sse_state);
    free(samples);
}

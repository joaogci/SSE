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
    double epsilon = 0.0;

    uint64_t seed = (uint64_t) time(NULL);
    long therm_cycles = 1e4;
    long mc_cycles = 1e4;
    int n_bins = 1;
    int n_loops = 1;

    double beta_vals[] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};
    int len_beta = sizeof(beta_vals) / sizeof(beta_vals[0]);

    struct heisenberg_system *hberg_system = (struct heisenberg_system *) malloc(sizeof(struct heisenberg_system));
    struct sse_state *sse_state = (struct sse_state *) malloc(sizeof(struct sse_state));    

    init_heisenberg_system(d, N, J, delta, h, epsilon, hberg_system);
    init_sse_state(seed, hberg_system, sse_state);

    // for (int i = 0; i < 6; i++) {
    //     printf("vtx_type -> %d; H -> %f \n", i + 1, sse_state->vtx_type[i].H);
    //     for (int li = 0; li < 4; li++) {
    //         for (int le = 0; le < 4; le++) {
    //             printf("li %d; le %d -> %f; new_vtx_type -> %d\n", li, le, sse_state->vtx_type[i].prob_exit[li][le], sse_state->vtx_type[i].new_vtx_type[li][le]);
    //         }
    //     }
    //     printf("\n");
    // }
    // exit(1);

    simulate_sse(beta_vals, len_beta, therm_cycles, mc_cycles, n_bins, n_loops, hberg_system, sse_state);

    free_memory(hberg_system, sse_state);
    free(hberg_system);
    free(sse_state);
}

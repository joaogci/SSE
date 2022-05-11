#include <stdio.h>
#include <stdlib.h>

#include "sse.h"

int main(int argc, char **argv) {
    double beta_vals[6] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};

    struct heisenberg_system *hberg_system = (struct heisenberg_system *) malloc(sizeof(struct heisenberg_system));
    struct sse_state *sse_state = (struct sse_state *) malloc(sizeof(struct sse_state));    

    init_heisenberg_system(1, 2, 1.0, 1.0, 0.0, 0.0, hberg_system);
    init_sse_state(10, hberg_system, sse_state);

    simulate_sse(beta_vals, 6, 1e4, 1e4, 1, 1, hberg_system, sse_state);

    free_memory(hberg_system, sse_state);
    free(hberg_system);
    free(sse_state);
}

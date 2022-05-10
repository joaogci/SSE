#include <stdio.h>
#include <stdlib.h>

#include "sse.h"

int main(int argc, char **argv) {
    double beta_vals[6] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};

    init_heisenberg_system(1, 2, 1.0, 1.0, 0.0, 0.0);
    init_sse_state(10);

    simulate_sse(beta_vals, 6, 1e4, 1e4, 1, 1);

}

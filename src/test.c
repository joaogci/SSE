#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "sse.h"

int main(int argc, char **argv) {

    int a;
    printf("%d \n", a);

    struct sse_state *sse_state = malloc(sizeof(struct sse_state));
    struct vertex *vertex = malloc(sizeof(struct vertex));

    uint64_t seed = (u_int64_t) time(NULL);
    printf("seed: %u \n", seed);
    init_system(4, &sse_state, seed);

    printf("N = %d \n", sse_state->N);
    printf("spin:\n [ ");
    for (int i = 0; i < sse_state->N; i++) {
        printf("%d, ", sse_state->spin[i]);
    }
    printf("] \n");
    printf("bond_site : \n");
    for (int b = 0; b < sse_state->N; b++) {
        printf("bond %d: %d -> %d \n", b, sse_state->bond_site[b][0], sse_state->bond_site[b][1]);
    }
    printf("\n");
    printf("M init = %d \n", sse_state->M);
    printf("opstring: \n [ ");
    for (int i = 0; i < sse_state->M; i++) {
        printf("%d, ", sse_state->opstring[i]);
    }
    printf("]\n");

    create_vertex_list(sse_state, &vertex);

    free_memory(&sse_state);

    return 0;
}

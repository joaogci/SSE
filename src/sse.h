#ifndef SSE_H
#define SSE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "xorshiro256++.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))


struct sse_state {
    int N;
    int Nb;
    int *spin;
    int **bond_site;

    int *opstring;
    int M;
    int n;
};

struct vertex {
    int *vertex_list;
    int *vertex_last;
    int *vertex_init;
    bool freed;
};

/*
 * Initializes the system. Allocates arrays and starts the system on a 
 * random spin configuration.
 * 
 */
void init_system(int N, struct sse_state **sse_state, uint64_t seed) {
    int i;
    // Initialize RNG
    for (i = 0; i < 4; i++) {
        s[i] = seed * i;
    }

    (*sse_state)->N = N;
    (*sse_state)->Nb = N;
    (*sse_state)->n = 0;

    (*sse_state)->spin = malloc(N * sizeof(int));
    (*sse_state)->bond_site = malloc(N * sizeof(int*));
    /*
     * TODO: add more boundary conditions
     *       random spin initialization
     * Init spins and lattice with PBC. Can change BC later.
     * Start with all spins up. 
     */
    for (i = 0; i < N; i++) {
        (*sse_state)->bond_site[i] = malloc(2 * sizeof(int));

        (*sse_state)->bond_site[i][0] = i;
        (*sse_state)->bond_site[i][1] = (i + 1) % N;

        (*sse_state)->spin[i] = 1;
        if (next_double() < 0.25) {
            (*sse_state)->spin[i] = -1;
        }
    }

    (*sse_state)->M = MAX(4, N / 4);
    (*sse_state)->opstring = malloc((*sse_state)->M * sizeof(int));

    for (i = 0; i < (*sse_state)->M; i++) {
        (*sse_state)->opstring[i] = 0;
    }
}

void create_vertex_list(struct sse_state *sse_state, struct vertex **vertex) {
    if ((*vertex)->freed != true) {
        printf("hello \n");
    } else {
        printf("bye \n");
    }
}

void free_memory(struct sse_state **sse_state) {
    free((*sse_state)->spin);
    free((*sse_state)->opstring);

    for (int i = 0; i < (*sse_state)->N; i++) {
        free((*sse_state)->bond_site[i]);
    }
    free((*sse_state)->bond_site);
    free((*sse_state));
    (*sse_state) = NULL;
}


#endif // SSE_H

#ifndef SSE_H
#define SSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "../rng/xorshiro256++.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))


struct sse_state {
    int d;
    int L;
    int N;
    int Nb;
    double beta;
    double pre_computed_prob;

    int *spin;
    int **bond_site;

    int *opstring;
    int M;
    int n;
};

struct vertex {
    int *list;
    int *last;
    int *init;
    bool freed;
};

/*
 * Initializes the system. Allocates arrays and starts the system on a 
 * random spin configuration.
 * 
 */
void init_system(int d, int L, double beta, struct sse_state *sse_state, struct vertex *vertex, uint64_t seed) {
    int i, j;
    // Initialize RNG
    for (i = 0; i < 4; i++) {
        s[i] = seed * i;
    }

    sse_state->L = L;
    sse_state->N = pow(L, d);
    sse_state->d = d;
    sse_state->Nb = sse_state->N * d;
    sse_state->beta = beta;
    sse_state->pre_computed_prob = sse_state->Nb * beta / 2.0;
    sse_state->n = 0;

    sse_state->spin = (int*) malloc(sse_state->N * sizeof(int));
    sse_state->bond_site = (int**) malloc(sse_state->Nb * sizeof(int*));
    /*
     * TODO: add more boundary conditions
     *       random spin initialization
     * Init spins and lattice with PBC. Can change BC later.
     * Start with all spins up. 
     */
    if (d == 1) {
        for (i = 0; i < L; i++) {
            sse_state->bond_site[i] = (int*) malloc(2 * sizeof(int));

            sse_state->bond_site[i][0] = i;
            sse_state->bond_site[i][1] = (i + 1) % sse_state->N;
        }
    } else if (d == 2) {
        for (i = 0; i < L; i++) {
            for (j = 0; j < L; j++) {
                sse_state->bond_site[j * L + i] = (int*) malloc(2 * sizeof(int));
                
                sse_state->bond_site[j * L + i][0] = j * L + i;
                sse_state->bond_site[j * L + i][1] = j * L + (i + 1) % L;

                sse_state->bond_site[(j * L + i) + sse_state->N][0] = j * L + i;
                sse_state->bond_site[(j * L + i) + sse_state->N][1] = ((j + 1) % L) * L + i;
            }
        }
    }

    for (i = 0; i < sse_state->N; i++) {
        sse_state->spin[i] = 1;
        if (next_double() < 0.5) {
            sse_state->spin[i] = -1;
        }
    }

    sse_state->M = MAX(4, sse_state->N / 4);
    sse_state->opstring = (int*) malloc(sse_state->M * sizeof(int));
    memset(sse_state->opstring, 0, sse_state->M * sizeof(int));

    vertex->init = (int*) malloc(sse_state->N * sizeof(int));
    vertex->last = (int*) malloc(sse_state->N * sizeof(int));
    vertex->list = (int*) malloc(4 * sse_state->M * sizeof(int));

    vertex->freed = false;
}

void create_vertex_list(struct sse_state *sse_state, struct vertex *vertex) {
    int p, i, v0, bond;
    int i1, i2, v1, v2;

    if (vertex->freed) {
        vertex->list = (int*) malloc(4 * sse_state->M * sizeof(int));

        vertex->freed = false;
    }

    memset(vertex->init, -1, sse_state->N * sizeof(int));
    memset(vertex->last, -1, sse_state->N * sizeof(int));
    memset(vertex->list, -1, 4 * sse_state->M * sizeof(int));

    for (p = 0; p < sse_state->M; p++) {
        if (sse_state->opstring[p] == 0) {
            continue;
        }

        v0 = 4 * p;
        bond = (sse_state->opstring[p] / 2) - 1;

        i1 = sse_state->bond_site[bond][0];
        i2 = sse_state->bond_site[bond][1];

        v1 = vertex->last[i1];
        v2 = vertex->last[i2];

        if (v1 != -1) {
            vertex->list[v1] = v0;
            vertex->list[v0] = v1;
        } else {
            vertex->init[i1] = v0;
        }

        if (v2 != -1) {
            vertex->list[v2] = v0 + 1;
            vertex->list[v0 + 1] = v2;
        } else {
            vertex->init[i2] = v0 + 1;
        }

        vertex->last[i1] = v0 + 2;
        vertex->last[i2] = v0 + 3;
    }

    for (i = 0; i < sse_state->N; i++) {
        if (vertex->init[i] != -1) {
            vertex->list[vertex->init[i]] = vertex->last[i];
            vertex->list[vertex->last[i]] = vertex->init[i];
        }
    }
}

void diag_update(struct sse_state *sse_state) {
    int p, b;
    for (p = 0; p < sse_state->M; p++) {
        if (sse_state->opstring[p] == 0) {
            b = next() % sse_state->Nb;
            if (sse_state->spin[sse_state->bond_site[b][0]] != sse_state->spin[sse_state->bond_site[b][1]]) {
                if (next_double() <= sse_state->pre_computed_prob / (sse_state->M - sse_state->n)) {
                    sse_state->opstring[p] = 2 * (b + 1);
                    sse_state->n++;
                }
            }
        } else if (sse_state->opstring[p] % 2 == 0) { 
            if (next_double() <= (sse_state->M - sse_state->n + 1.0) / sse_state->pre_computed_prob) { 
                sse_state->opstring[p] = 0;
                sse_state->n--;
            }
        } else {
            b = (sse_state->opstring[p] / 2) - 1;
            sse_state->spin[sse_state->bond_site[b][0]] = - sse_state->spin[sse_state->bond_site[b][0]];
            sse_state->spin[sse_state->bond_site[b][1]] = - sse_state->spin[sse_state->bond_site[b][1]];
        }
    }
}

void loop_update(struct sse_state *sse_state, struct vertex *vertex) {
    int v0, v_in, v_out;

    for (v0 = 0; v0 < 4 * sse_state->M; v0+=2) {
        if (vertex->list[v0] < 0) {
            continue;
        }

        v_in = v0;
        if (next_double() <= 0.5) {
            while (true) {
                sse_state->opstring[v_in / 4] ^= 1;
                vertex->list[v_in] = -2;
                v_out = v_in ^ 1;
                v_in = vertex->list[v_out];
                vertex->list[v_out] = -2;

                if (v_in == v0) {
                    break;
                }
            }
        } else {
            while (true) {
                vertex->list[v_in] = -1;
                v_out = v_in ^ 1;
                v_in = vertex->list[v_out];
                vertex->list[v_out] = -1;

                if (v_in == v0) {
                    break;
                }
            }
        }
    }

    int i;
    for (i = 0; i < sse_state->N; i++) {
        if (vertex->init[i] != -1) {
            if (vertex->list[vertex->init[i]] == -2) {
                sse_state->spin[i] = - sse_state->spin[i];
            }
        } else if (next_double() <= 0.5) {
            sse_state->spin[i] = - sse_state->spin[i];
        }
    }
}

void ajust_cutoff(struct sse_state *sse_state, struct vertex *vertex) {
    int M_new = sse_state->n * 1.33;

    if (M_new > sse_state->M) {
        int opstring_cpy[sse_state->M];
        memcpy(opstring_cpy, sse_state->opstring, sse_state->M * sizeof(int));
        free(sse_state->opstring); 
        sse_state->opstring = (int*) malloc(M_new * sizeof(int));
        memset(sse_state->opstring, 0, M_new * sizeof(int));
        memcpy(sse_state->opstring, opstring_cpy, sse_state->M * sizeof(int));

        sse_state->M = M_new;

        free(vertex->list);
        vertex->freed = true;
    }
}

void free_vertex(struct vertex **vertex) {
    free((*vertex)->init);
    free((*vertex)->last);
    free((*vertex)->list);
    free((*vertex));
    *vertex = NULL;
}

void free_state(struct sse_state **sse_state) {
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

#include "sse.h"

/* 
 * function: diag_update
 *  diagonal update for the SSE MCS
 *  inserts or removes a diagonal operator from the operator string
 *   with some probablity
 * 
 *  paramters:
 *      (double) beta: temperature
 *      (heisenberg_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void diag_update(double beta, heisenberg_system *system, sse_state *state, pcg32_random_t* rng) 
{
    for (int p = 0; p < state->M; p++) {
        if (state->op_string[p] == 0) {
            // no operator -> insert

            int b = pcg32_boundedrand_r(rng, system->Nb);
            if (pcg32_double_r(rng) <= (system->Nb * beta * prob(b, system, state)) / (state->M - state->n)) {
                state->op_string[p] = 3 * (b + 1);
                state->n++;
            }
        } else if (state->op_string[p] % 3 == 0) {
            // diagonal operator -> remove

            int b = (state->op_string[p] / 3) - 1;
            if (pcg32_double_r(rng) <= (state->M - state->n + 1) / (beta * system->Nb * prob(b, system, state))) {
                state->op_string[p] = 0;
                state->n--;
            }
        } else {
            // off-diagonal operator -> propagate

            int b = (state->op_string[p] / 3) - 1;
            int a = state->op_string[p] % 3;

            if (a == 1) {
                system->spin[system->bond[b][0]] += 2;
                system->spin[system->bond[b][1]] += -2;
            } else if (a == 2) {
                system->spin[system->bond[b][0]] += -2;
                system->spin[system->bond[b][1]] += 2;
            }
        }
    }
}

/* 
 * function: loop_update
 *  loop update for the SSE MCS
 *  creates a loop within the vertex list and updates the configuration
 *  off-diagonal update
 * 
 *  parameters:
 *      (heisenberg_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void loop_update(heisenberg_system *system, sse_state *state, pcg32_random_t* rng) 
{
    // if no operators just randomize the spins
    if (state->n == 0) {
        for (int loop = 0; loop < state->n_loops; loop++) {
            for (int i = 0; i < system->N; i++) {
                // update always but choose an update randomly
                system->spin[i] = system->Sz[pcg32_boundedrand_r(rng, system->n_proj)];
            }
        }

        return;
    }

    for (int loop = 0; loop < state->n_loops; loop++) {
        int j0 = pcg32_boundedrand_r(rng, 4 * state->n); /* select the first in leg randomly */
        int j = j0;

        // select first update randomly
        int update = 2;
        int update_idx = 1;
        if (pcg32_double_r(rng) <= 0.5) { update = -2; update_idx = 0; }

        // if the update is not allowed on the selected vertex, quit the loop update
        if (abs(state->vtx_type[state->vtx[j / 4]].spin[j % 4] + update) > 2 * system->S) {
            continue; 
        }

        int loop_size = 0;
        // begin the loop 
        while (true) {
            int p = j / 4;
            int li = j % 4;
            int le;

            // select an exit leg
            double r = pcg32_double_r(rng);
            for (le = 0; le < 4; le++) {
                if (state->vtx_type[state->vtx[p]].prob_exit[update_idx][li][le] == 1.0 || 
                        r < state->vtx_type[state->vtx[p]].prob_exit[update_idx][li][le]) {
                    break;
                }
            } 

            // the new update is given by (state_new[le] - state_old[le])
            int old_state = state->vtx_type[state->vtx[p]].spin[le];
            state->vtx[p] = state->vtx_type[state->vtx[p]].new_vtx_type[update_idx][li][le];
            int new_state = state->vtx_type[state->vtx[p]].spin[le];

            if (new_state - old_state == -2) {
                update_idx = 0;
            } else if (new_state - old_state == 2) {
                update_idx = 1;
            } else if (new_state - old_state == 0) {
                if (update_idx == 1) { update_idx = 0; }
                else if (update_idx == 0) { update_idx = 1; }
            }

            j = 4 * p + le;

            // record the loop size if no bounce happens
            if (li != le) { loop_size++; }

            if (j == j0) { break; }
            
            // go to linking vertex
            j = state->link[j];
            if (j == j0) { break; }

            // end it loop is too large
            if (loop_size >= 60 * state->M) {
                printf("aborted loop update \n");
                return;
            }
        }

        state->loop_size += loop_size;
    }

    // translate the updated vertex list to the string operator
    for (int p = 0; p < state->n; p++) {
        state->red_op_string[p] = 3 * (state->red_op_string[p] / 3) + state->vtx_type[state->vtx[p]].type;
        state->op_string[state->trans_op_string[p]] = state->red_op_string[p];
    }

    // flip the update spins and flip randomly the non-updated ones
    for (int i = 0; i < system->N; i++) {
        if (state->first[i] != -1) {
            int p = state->first[i] / 4;
            int l = state->first[i] % 4;
            system->spin[i] = state->vtx_type[state->vtx[p]].spin[l];
        }
        else {
            system->spin[i] = system->Sz[pcg32_boundedrand_r(rng, system->n_proj)];
        }
    }
}

/* 
 * function: ajust_cutoff
 *  dinamically adjusts the expansion cutoff during the thermalization
 *   part of the simulation
 * 
 *  parameters:
 *      (sse_state *) state: SSE state
 *      (bool) adjust_loop
 */
void ajust_cutoff(sse_state *state, bool adjust_loop) 
{
    u_int64_t M_new = state->n * 1.33;

    if (M_new > state->M) {
        u_int64_t opstring_cpy[state->M];
        memcpy(opstring_cpy, state->op_string, state->M * sizeof(u_int64_t));
        free(state->op_string); 
        state->op_string = (u_int64_t*) malloc(M_new * sizeof(u_int64_t));
        memset(state->op_string, 0, M_new * sizeof(u_int64_t));
        memcpy(state->op_string, opstring_cpy, state->M * sizeof(u_int64_t));

        state->M = M_new;
    }

    if (adjust_loop) {
        state->n_loops = 2 * state->M * state->n_loops / state->loop_size;
        state->loop_size = 0;
    }
}

/* 
 * function: init_heisenberg_system 
 *  initializes the heisenberg_system struct
 * 
 *  parameters:
 *      (int) d: dimension
 *      (int) L: number of unit cells
 *      (int) boundary_cond: boundary condition of the lattice
 *      (double) S: spin quantum number
 *      (double) delta: z-axis anisotropy strength
 *      (double) h: z-axis magnetic field strength
 *      (double) epsilon: constant added to the Hamiltonian
 *      (heisenberg_system *) system: system to be initialized
 */
void init_heisenberg_system(int d, int L, int boundary_cond, double S, double delta, double h, double epsilon, heisenberg_system *system) 
{
    system->d = d;
    system->L = L;
    system->boundary_cond = boundary_cond;
    system->N = pow(L, d);
    system->Nb = system->N * d - system->boundary_cond;

    system->S = S;
    system->n_proj = 2.0 * S + 1.0;
    system->Sz = (int *) malloc(system->n_proj * sizeof(int));
    for (int i = 0; i < system->n_proj; i++) {
        system->Sz[i] = - 2 * (S - i);
    }
    
    system->h = h;
    system->delta = delta;
    system->epsilon = epsilon;

    system->spin = (int *) malloc(system->N * sizeof(int));
    system->bond = (int **) malloc(system->Nb * sizeof(int *));
    for (int i = 0; i < system->Nb; i++) { system->bond[i] = (int*) malloc(2 * sizeof(int)); }
    if (d == 1) {
        for (int i = 0; i < system->Nb; i++) {
            system->bond[i][0] = i;
            system->bond[i][1] = (i + 1) % system->N;
        }
    } else if (d == 2) {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                system->bond[j * L + i][0] = j * L + i;
                system->bond[j * L + i][1] = j * L + (i + 1) % L;

                system->bond[(j * L + i) + system->N][0] = j * L + i;
                system->bond[(j * L + i) + system->N][1] = ((j + 1) % L) * L + i;
            }
        }
    } else if (d == 3) {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                for (int k = 0; k < L; k++) {
                    system->bond[k * L * L + j * L + i][0] = k * L * L + j * L + i;
                    system->bond[k * L * L + j * L + i][1] = k * L * L + j * L + (i + 1) % L;

                    system->bond[(k * L * L + j * L + i) + system->N][0] = k * L * L + j * L + i;
                    system->bond[(k * L * L + j * L + i) + system->N][1] = k * L * L + ((j + 1) % L) * L + i;

                    system->bond[(k * L * L + j * L + i) + 2 * system->N][0] = k * L * L + j * L + i;
                    system->bond[(k * L * L + j * L + i) + 2 * system->N][1] = ((k + 1) % L) * L * L + ((j + 1) % L) * L + i;
                }
            }
        }
    }

    for (int i = 0; i < system->N; i++) {
        system->spin[i] = system->Sz[0];
    }
}

/* 
 * function: init_sse_state 
 *  initializes the sse_state struct
 * 
 *  parameters:
 *      (uint64_t) seed: seed for the RNG
 *      (heisenberg_system *) system: simulated system
 *      (sse_state *) state: sse_state to be initialized
 */
void init_sse_state(uint64_t seed, heisenberg_system *system, sse_state *state) 
{
    state->n = 0;
    state->M = MAX_(4, system->N / 4);
    state->op_string = (u_int64_t *) malloc(state->M * sizeof(u_int64_t));
    memset(state->op_string, 0, state->M * sizeof(u_int64_t));

    state->n_loops = MAX_(4, system->N / 4);
    state->loop_size = 0;
    state->first = (u_int64_t *) malloc(system->N * sizeof(u_int64_t));
    state->vtx = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
    state->link = (u_int64_t *) malloc(4 * state->n * sizeof(u_int64_t)); 
    state->red_op_string = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
    state->trans_op_string = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
}

/* 
 * function: reset_sse_state
 *  resets the SSE state
 *  deletes all of the operator string and spin states
 * 
 *  parameters:
 *      (heisenberg_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void reset_sse_state(heisenberg_system *system, sse_state *state) 
{
    for (int i = 0; i < system->N; i++) {
        system->spin[i] = system->Sz[0];
    }

    state->n = 0;
    state->M = MAX_(4, system->N / 4);
    state->loop_size = 0;
    state->n_loops = MAX_(4, system->N / 4);

    free(state->op_string);
    state->op_string = (u_int64_t *) malloc(state->M * sizeof(u_int64_t));
    memset(state->op_string, 0, state->M * sizeof(u_int64_t));
}

/* 
 * function: create_vtx_list
 *  tranlates the operator string in to a doubly linked list of vertecies
 * 
 *  parameters:
 *      (heisenberg_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void create_vtx_list(heisenberg_system *system, sse_state *state) 
{
    free(state->vtx);
    free(state->link);
    free(state->red_op_string);
    free(state->trans_op_string);
    state->vtx = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
    state->link = (u_int64_t *) malloc(4 * state->n * sizeof(u_int64_t)); 
    state->red_op_string = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
    state->trans_op_string = (u_int64_t *) malloc(state->n * sizeof(u_int64_t));
    
    u_int64_t last[system->N];
    memset(last, -1, system->N * sizeof(u_int64_t));
    memset(state->first, -1, system->N * sizeof(u_int64_t));
    memset(state->link, -1, 4 * state->n * sizeof(u_int64_t));
    memset(state->vtx, -1, state->n * sizeof(u_int64_t));
    memset(state->red_op_string, 0, state->n * sizeof(u_int64_t));
    memset(state->trans_op_string, 0, state->n * sizeof(u_int64_t));
    
    // reduce the operator string to the times where it has an operator
    // create a translatio between the full and reduced one
    int p_red = 0;
    for (int p = 0; p < state->M; p++) {
        if (state->op_string[p] != 0) {
            state->red_op_string[p_red] = state->op_string[p];
            state->trans_op_string[p_red] = p;
            p_red++;
        }
    }

    // create the vertex list
    int l[4];
    for (int p = 0; p < state->n; p++) {
        // propagate the system
        int b = (state->red_op_string[p] / 3) - 1;
        int a = state->red_op_string[p] % 3;
        int v0 = 4 * p;

        int i1 = system->bond[b][0];
        int i2 = system->bond[b][1];

        l[0] = system->spin[i1];
        l[1] = system->spin[i2];
        if (a == 1) {
            system->spin[i1] += 2;
            system->spin[i2] += -2;
        } else if (a == 2) {
            system->spin[i1] += -2;
            system->spin[i2] += 2;
        }
        l[2] = system->spin[i1];
        l[3] = system->spin[i2];

        // vertex at time p
        for (int i = 0; i < state->n_diagrams; i++) {
            if (l[0] == state->vtx_type[i].spin[0] && 
            l[1] == state->vtx_type[i].spin[1] && 
            l[2] == state->vtx_type[i].spin[2] && 
            l[3] == state->vtx_type[i].spin[3]) { state->vtx[p] = state->vtx_type[i].indx; }
        }
        
        // link the legs of the vertex
        int v1 = last[i1];
        int v2 = last[i2];

        if (v1 != -1) {
            state->link[v1] = v0;
            state->link[v0] = v1;
        } else {
            state->first[i1] = v0;
        }

        if (v2 != -1) {
            state->link[v2] = v0 + 1;
            state->link[v0 + 1] = v2;
        } else {
            state->first[i2] = v0 + 1;
        }

        last[i1] = v0 + 2;
        last[i2] = v0 + 3;
    }

    for (int i = 0; i < system->N; i++) {
        if (state->first[i] != -1) {
            state->link[state->first[i]] = last[i];
            state->link[last[i]] = state->first[i];
        }
    }
}

/* 
 * function: prob
 *  computes probability for diagonal updates
 *  
 *  parameters:
 *      (int) b: bond for the insertion/removal of operator
 *      (heisenberg_system *) system: system
 *      (sse_state* ) state: SSE state
 *  
 *  returns:
 *      (double) prob: probability for the update
 */
double prob(int b, heisenberg_system *system, sse_state *state) 
{
    for (int i = 0; i < state->n_diagrams; i++) {
        if (state->vtx_type[i].type == 0 && state->vtx_type[i].spin[0] == system->spin[system->bond[b][0]] && state->vtx_type[i].spin[1] == system->spin[system->bond[b][1]]) {
            return state->vtx_type[i].H;
        }
    }

    return 0.0;
}

/* 
 * function: free_memory
 *  frees allocated memory in the heisenberg_system and sse_state structs
 * 
 *  parameters:
 *      (heisenberg_system *) system: heisenberg_system struct
 *      (sse_state *) state: sse_state struct
 */
void free_memory(heisenberg_system *system, sse_state *state) 
{
    for (int i = 0; i < system->Nb; i++) { free(system->bond[i]); }
    free(system->bond);
    free(system->spin);
    free(system->Sz);

    free(state->first);
    free(state->op_string);

    free(state->vtx);
    free(state->link);
    free(state->red_op_string);
    free(state->trans_op_string);

    free(state->vtx_type);
}

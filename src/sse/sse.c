#include "sse.h"

void diag_update(double beta, heisenberg_system *system, sse_state *state) 
{
    for (int p = 0; p < state->M; p++) {
        if (state->op_string[p] == 0) {
            int b = next() % system->Nb;
            if (next_double() <= MIN(1.0, (system->Nb * beta * prob(b, system)) / (state->M - state->n))) {
                state->op_string[p] = 2 * (b + 1);
                state->n++;
            }
        } else if (state->op_string[p] % 2 == 0) {
            int b = (state->op_string[p] / 2) - 1;
            if (next_double() <= MIN(1.0, (state->M - state->n + 1) / (beta * system->Nb * prob(b, system)))) {
                state->op_string[p] = 0;
                state->n--;
            }
        } else {
            int b = (state->op_string[p] / 2) - 1;
            system->spin[system->bond[b][0]] = - system->spin[system->bond[b][0]];
            system->spin[system->bond[b][1]] = - system->spin[system->bond[b][1]];
        }
    }
}

void loop_update(heisenberg_system *system, sse_state *state) 
{
    if (state->n == 0) {
        for (int i = 0; i < system->N; i++) {
            if (next_double() <= 0.5) { system->spin[i] = - system->spin[i]; }
        }
        return;
    }

    state->vtx = (int *) malloc(state->n * sizeof(int));
    state->link = (int *) malloc(4 * state->n * sizeof(int)); 
    
    int red_op_string[state->n];
    int trans_op_string[state->n];
    create_vtx_list(system, state, red_op_string, trans_op_string);

    int j0 = next() % (4 * state->n);
    int j = j0;
    while (true) {
        int p = j / 4;
        int li = j % 4;
        state->loop_size++;

        double r = next_double();
        int vtx_type = state->vtx[p];
        int le;
        for (le = 0; le < 4; le++) {
            if (r <= state->vtx_type[vtx_type].prob_exit[li][le]) {
                state->vtx[p] = state->vtx_type[vtx_type].new_vtx_type[li][le];
                break;
            }
        }
        
        j = 4 * p + le;
        state->loop_size++;
        if (j == j0) { break; }
        
        j = state->link[j];
        if (j == j0) { break; }
    }

    for (int p = 0; p < state->n; p++) {
        state->op_string[trans_op_string[p]] = 2 * (red_op_string[p] / 2) + state->vtx_type[state->vtx[p]].type;
    }

    for (int i = 0; i < system->N; i++) {
        if (state->first[i] != -1) {
            int p = state->first[i] / 4;
            int l = state->first[i] % 4;
            system->spin[i] = state->vtx_type[state->vtx[p]].spin[l];
        }
        else if (next_double() <= 0.5){ system->spin[i] = - system->spin[i]; }
    }

    free(state->vtx);
    free(state->link);
}

void ajust_cutoff(sse_state *state, bool adjust_loop) 
{
    int M_new = state->n * 1.33;

    if (M_new > state->M) {
        int opstring_cpy[state->M];
        memcpy(opstring_cpy, state->op_string, state->M * sizeof(int));
        free(state->op_string); 
        state->op_string = (int*) malloc(M_new * sizeof(int));
        memset(state->op_string, 0, M_new * sizeof(int));
        memcpy(state->op_string, opstring_cpy, state->M * sizeof(int));

        state->M = M_new;
    }

    if (adjust_loop) {
        if (state->loop_size > 1000) { state->n_loops = 2000 * state->M * state->n_loops / state->loop_size; }
        state->loop_size = 0;
    }
}

void init_heisenberg_system(int d, int L, double J, double delta, double h, double epsilon, heisenberg_system *system) 
{
    system->d = d;
    system->L = L;
    system->N = pow(L, d);
    system->Nb = system->N * d; // For PBC
    
    system->J = J;
    system->h = h;
    system->delta = delta;
    system->epsilon = epsilon;

    system->prob[0] = C + epsilon + 0.25 * delta; 
    system->prob[1] = C + epsilon - 0.25 * delta + hb;
    system->prob[2] = C + epsilon - 0.25 * delta - hb;

    system->spin = (int *) malloc(system->N * sizeof(int));
    system->bond = (int **) malloc(system->Nb * sizeof(int *));
    for (int i = 0; i < system->Nb; i++) {
        system->bond[i] = (int*) malloc(2 * sizeof(int));
    }
    if (d == 1) {
        for (int i = 0; i < L; i++) {
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
    }

    for (int i = 0; i < system->N; i++) {
        system->spin[i] = 1;
        if (2 * i < system->N) { system->spin[i] = -1; }
    }
}

void init_sse_state(uint64_t seed, heisenberg_system *system, sse_state *state) 
{
    for (int i = 0; i < 4; i++) { s[i] = seed * (i + 1); }
    state->vtx_type = create_vtx_type_list(system->J, system->delta, system->h, system->epsilon);

    state->n = 0;
    state->M = MAX(4, system->N / 4);
    state->op_string = (int *) malloc(state->M * sizeof(int));
    memset(state->op_string, 0, state->M * sizeof(int));

    state->n_loops = MAX(4, system->N / 4);
    state->loop_size = 0;
    state->first = (int *) malloc(system->N * sizeof(int));
}

void reset_sse_state(heisenberg_system *system, sse_state *state) 
{
    for (int i = 0; i < system->N; i++) {
        system->spin[i] = 1;
        if (next_double() < 0.5) { system->spin[i] = -1; }
    }

    state->n = 0;
    state->M = MAX(4, system->N / 4);
    state->loop_size = 0;
    state->n_loops = MAX(4, system->N / 4);

    free(state->op_string);
    state->op_string = (int *) malloc(state->M * sizeof(int));
    memset(state->op_string, 0, state->M * sizeof(int));
}

void create_vtx_list(heisenberg_system *system, sse_state *state, int *red_op_string, int *trans_op_string) 
{
    int last[system->N];
    memset(last, -1, system->N * sizeof(int));
    memset(state->first, -1, system->N * sizeof(int));
    memset(state->link, -1, 4 * state->n * sizeof(int));
    memset(state->vtx, -1, state->n * sizeof(int));
    memset(red_op_string, 0, state->n * sizeof(int));
    memset(trans_op_string, 0, state->n * sizeof(int));
    
    int p_red = 0;
    for (int p = 0; p < state->M; p++) {
        if (state->op_string[p] != 0) {
            red_op_string[p_red] = state->op_string[p];
            trans_op_string[p_red] = p;
            p_red++;
        }
    }

    int l[4];
    for (int p = 0; p < state->n; p++) {
        int b = (red_op_string[p] / 2) - 1;
        int v0 = 4 * p;

        int i1 = system->bond[b][0];
        int i2 = system->bond[b][1];

        l[0] = system->spin[i1];
        l[1] = system->spin[i2];
        if (red_op_string[p] % 2 != 0) {
            system->spin[i1] = - system->spin[i1];
            system->spin[i2] = - system->spin[i2];
        }
        l[2] = system->spin[i1];
        l[3] = system->spin[i2];

        for (int i = 0; i < N_DIAGRAMS; i++) {
            if (l[0] == state->vtx_type[i].spin[0] && 
            l[1] == state->vtx_type[i].spin[1] && 
            l[2] == state->vtx_type[i].spin[2] && 
            l[3] == state->vtx_type[i].spin[3]) { state->vtx[p] = state->vtx_type[i].indx; }
        }

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

double prob(int b, heisenberg_system *system) 
{
    if (system->spin[system->bond[b][0]] != system->spin[system->bond[b][1]]) {
        return system->prob[0];
    } else if (system->spin[system->bond[b][0]] == system->spin[system->bond[b][1]]) {
        if (system->spin[system->bond[b][1]] == 1) {
            return system->prob[1];
        } else {
            return system->prob[2];
        }
    }

    return 0.0;
}

void free_memory(heisenberg_system *system, sse_state *state) 
{
    for (int i = 0; i < system->N; i++) { free(system->bond[i]); }
    free(system->bond);
    free(system->spin);

    free(state->first);
    free(state->vtx_type);
    free(state->op_string);
}

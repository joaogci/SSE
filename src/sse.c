#include "sse.h"

void simulate_sse(double *beta_vals, int len_beta, long mc_cycles, long therm_cycles, int n_bins, int n_loops, struct heisenberg_system *hberg_system, struct sse_state *sse_state) {
    int n, t_idx, loop;
    long t;
    double beta;
    double n_mean;

    for (t_idx = 0; t_idx < len_beta; t_idx++) {
        beta = beta_vals[t_idx];
        reset_sse_state(hberg_system, sse_state);

        for (t = 0; t < therm_cycles; t++) {
            diag_update(beta, hberg_system, sse_state);
            for (loop = 0; loop < n_loops; loop++) {
                loop_update(hberg_system, sse_state);
            }
            ajust_cutoff(sse_state);
        }

        n_mean = 0.0;

        for (n = 0; n < n_bins; n++) {
            for (t = 0; t < mc_cycles; t++) {
                diag_update(beta, hberg_system, sse_state);
                for (loop = 0; loop < n_loops; loop++) {
                    loop_update(hberg_system, sse_state);
                }
                sample(hberg_system, sse_state);
                n_mean += sse_state->n;
            }
        }

        n_mean = n_mean / (mc_cycles * n_bins);
        double E = - n_mean / (beta * hberg_system->N) + hberg_system->J * hberg_system->C;
        printf("beta: %f; n_mean: %f; E: %f \n", beta, n_mean, E);
    }
}

void sample(struct heisenberg_system *hberg_system, struct sse_state *sse_state) {

}

void init_heisenberg_system(int d, int N, double J, double delta, double h, double epsilon, struct heisenberg_system *hberg_system) {
    hberg_system->d = d;
    hberg_system->N = N;
    hberg_system->Nb = N * d;         // For PBC
    
    hberg_system->J = J;
    hberg_system->h = h;
    hberg_system->delta = delta;
    hberg_system->epsilon = epsilon;
    hberg_system->C = 0.25 * delta + 0.5 * h / J + epsilon;
    hberg_system->H = create_hamiltonian(J, delta, h, hberg_system->C);

    hberg_system->spin = (int *) malloc(N * sizeof(int));
    hberg_system->bond = (int **) malloc(N * sizeof(int *));
    for (int i = 0; i < N; i++) {
        hberg_system->bond[i] = (int*) malloc(2 * sizeof(int));

        hberg_system->bond[i][0] = i;
        hberg_system->bond[i][1] = (i + 1) % N;

        hberg_system->spin[i] = 1;
        if (next_double() < 0.5) {
            hberg_system->spin[i] = -1;
        }
    }
}

void init_sse_state(uint64_t seed, struct heisenberg_system *hberg_system, struct sse_state *sse_state) {
    for (int i = 0; i < 4; i++) {
        s[i] = seed * i;
    }
    sse_state->vtx_type = create_vtx_type_list(hberg_system->J, hberg_system->delta, hberg_system->h, hberg_system->C);

    sse_state->n = 0;
    sse_state->M = MAX(4, hberg_system->N / 4);
    sse_state->op_string = (int *) malloc(sse_state->M * sizeof(int));
    memset(sse_state->op_string, 0, sse_state->M * sizeof(int));

    sse_state->link = (int *) malloc(4 * sse_state->M * sizeof(int));
    sse_state->first = (int *) malloc(hberg_system->N * sizeof(int));
    sse_state->vtx = (int *) malloc(sse_state->M * sizeof(int));
}

void reset_sse_state(struct heisenberg_system *hberg_system, struct sse_state *sse_state) {
    for (int i = 0; i < hberg_system->N; i++) {
        hberg_system->spin[i] = 1;
        if (next_double() < 0.5) {
            hberg_system->spin[i] = -1;
        }
    }

    sse_state->n = 0;
    sse_state->M = MAX(4, hberg_system->N / 4);

    free(sse_state->op_string);
    sse_state->op_string = (int *) malloc(sse_state->M * sizeof(int));
    memset(sse_state->op_string, 0, sse_state->M * sizeof(int));

    if (!sse_state->freed) {
        free(sse_state->link);
        free(sse_state->vtx);
    }
    sse_state->link = (int *) malloc(4 * sse_state->M * sizeof(int));
    sse_state->vtx = (int *) malloc(sse_state->M * sizeof(int));
    sse_state->freed = false;    
}

void diag_update(double beta, struct heisenberg_system *hberg_system, struct sse_state *sse_state) {
    int p, b;
    double prob;

    for (p = 0; p < sse_state->M; p++) {
        prob = 0.0;

        if (sse_state->op_string[p] == 0) {
            b = next() % hberg_system->Nb;
            
            if (hberg_system->spin[hberg_system->bond[b][0]] != hberg_system->spin[hberg_system->bond[b][1]]) {
                prob = hberg_system->H[1][1].value;
            }
            else if (hberg_system->spin[hberg_system->bond[b][0]] == hberg_system->spin[hberg_system->bond[b][1]]) {
                if (hberg_system->spin[hberg_system->bond[b][1]] == 1) {
                    prob = hberg_system->H[0][0].value;
                }
                else {
                    prob = hberg_system->H[3][3].value;
                }
            }

            if (next_double() <= MIN(1.0, (hberg_system->Nb * beta * prob) / (sse_state->M - sse_state->n))) {
                sse_state->op_string[p] = 2 * (b + 1);
                sse_state->n++;
            }
        }
        else if (sse_state->op_string[p] % 2 == 0) {
            b = (sse_state->op_string[p] / 2) - 1;

            if (hberg_system->spin[hberg_system->bond[b][0]] != hberg_system->spin[hberg_system->bond[b][1]]) {
                prob = hberg_system->H[1][1].value;
            }
            else if (hberg_system->spin[hberg_system->bond[b][0]] == hberg_system->spin[hberg_system->bond[b][1]]) {
                if (hberg_system->spin[hberg_system->bond[b][1]] == 1) {
                    prob = hberg_system->H[0][0].value;
                }
                else {
                    prob = hberg_system->H[3][3].value;
                }
            }

            if (next_double() <= MIN(1.0, (sse_state->M - sse_state->n + 1) / (beta * hberg_system->Nb * prob))) {
                sse_state->op_string[p] = 0;
                sse_state->n--;
            }
        }
        else {
            b = (sse_state->op_string[p] / 2) - 1;
            hberg_system->spin[hberg_system->bond[b][0]] = - hberg_system->spin[hberg_system->bond[b][0]];
            hberg_system->spin[hberg_system->bond[b][1]] = - hberg_system->spin[hberg_system->bond[b][1]];
        }
    }
}

void loop_update(struct heisenberg_system *hberg_system, struct sse_state *sse_state) {
    int j0, j;
    int b, p, li, le, l;
    int i;
    double r, cum_prob;
    int vtx_type;

    create_vtx_list(hberg_system, sse_state);

    if (sse_state->n == 0) {
        for (i = 0; i < hberg_system->N; i++) {
            if (sse_state->first[i] == -1 && next_double() <= 0.5) {
                hberg_system->spin[i] = - hberg_system->spin[i];
            }
        }
        return;
    }

    do {
        j0 = next() % (4 * sse_state->M);
    } while (sse_state->link[j0] < 0 && sse_state->vtx[j0 / 4] == 0);

    j = j0;
    while (true) {
        p = j / 4;
        li = j % 4;

        r = next_double();
        cum_prob = 0.0;
        vtx_type = sse_state->vtx[p] - 1;
        for (le = 0; le < 4; le++) {
            cum_prob += sse_state->vtx_type[vtx_type].prob_exit[li][le];
            if (r <= cum_prob) {
                sse_state->vtx[p] = sse_state->vtx_type[vtx_type].new_vtx_type[li][le];
                break;
            }
        }
        
        j = 4 * p + le;
        if (j == j0) {
            break;
        }
        
        j = sse_state->link[j];
        if (j == j0) {
            break;
        }
    }

    for (p = 0; p < sse_state->M; p++) {
        if (sse_state->vtx[p] != 0) {
            b = (sse_state->op_string[p] / 2) - 1;
            sse_state->op_string[p] = 2 * (b + 1) + sse_state->vtx_type[sse_state->vtx[p] - 1].type;
        }
    }

    for (i = 0; i < hberg_system->N; i++) {
        if (sse_state->first[i] != -1) {
            p = sse_state->first[i] / 4;
            l = sse_state->first[i] % 4;
            hberg_system->spin[i] = sse_state->vtx_type[sse_state->vtx[p] - 1].spin[l];
        }
        else {
            if (next_double() <= 0.5) {
                hberg_system->spin[i] = - hberg_system->spin[i];
            }
        }
    }
}

/*
    // Write to file
    FILE *fp;
    fp = fopen("1D_heisenberg_L16.csv", "w");

    fprintf(fp, "beta,E,C,m,m2,m_s,m2_s,n\n");
    for (int i = 0; i < beta_len; i++) {
        fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f\n", beta_vals[i], E_mean[i], C_mean[i], m_mean[i], m2_mean[i], ms_mean[i], m2s_mean[i], n_mean[i]);
    }

    fclose(fp);
*/

void ajust_cutoff(struct sse_state *sse_state) {
    int M_new = sse_state->n * 1.33;

    if (M_new > sse_state->M) {
        int opstring_cpy[sse_state->M];
        memcpy(opstring_cpy, sse_state->op_string, sse_state->M * sizeof(int));
        free(sse_state->op_string); 
        sse_state->op_string = (int*) malloc(M_new * sizeof(int));
        memset(sse_state->op_string, 0, M_new * sizeof(int));
        memcpy(sse_state->op_string, opstring_cpy, sse_state->M * sizeof(int));

        sse_state->M = M_new;

        free(sse_state->link);
        free(sse_state->vtx);
        sse_state->freed = true;
    }
}

void create_vtx_list(struct heisenberg_system *hberg_system, struct sse_state *sse_state) {
    int p, i, v0, b;
    int i1, i2, v1, v2;
    int last[hberg_system->N];
    int l[4];

    if (sse_state->freed) {
        sse_state->link = (int *) malloc(4 * sse_state->M * sizeof(int));
        sse_state->vtx = (int *) malloc(sse_state->M * sizeof(int));

        sse_state->freed = false;
    }

    memset(last, -1, hberg_system->N * sizeof(int));
    memset(sse_state->first, -1, hberg_system->N * sizeof(int));
    memset(sse_state->link, -1, 4 * sse_state->M * sizeof(int));
    memset(sse_state->vtx, 0, sse_state->M * sizeof(int));

    for (p = 0; p < sse_state->M; p++) {
        if (sse_state->op_string[p] == 0) {
            continue;
        }

        b = (sse_state->op_string[p] / 2) - 1;
        v0 = 4 * p;

        if (sse_state->op_string[p] % 2 != 0) {
            hberg_system->spin[hberg_system->bond[b][0]] = - hberg_system->spin[hberg_system->bond[b][0]];
            hberg_system->spin[hberg_system->bond[b][1]] = - hberg_system->spin[hberg_system->bond[b][1]];

            l[0] = - hberg_system->spin[hberg_system->bond[b][1]];
            l[1] = - hberg_system->spin[hberg_system->bond[b][0]];
            l[2] = hberg_system->spin[hberg_system->bond[b][1]];
            l[3] = hberg_system->spin[hberg_system->bond[b][0]];
        }
        else if (sse_state->op_string[p] % 2 == 0) {
            l[0] = hberg_system->spin[hberg_system->bond[b][1]];
            l[1] = hberg_system->spin[hberg_system->bond[b][0]];
            l[2] = hberg_system->spin[hberg_system->bond[b][1]];
            l[3] = hberg_system->spin[hberg_system->bond[b][0]];
        }

        if (l[0] == l[1] && l[1] == l[2] && l[2] == l[3] && l[3] == -1) {
            sse_state->vtx[p] = 1;
        }
        else if (l[0] == -1 && l[1] == 1 && l[2] == -1 && l[3] == 1) {
            sse_state->vtx[p] = 2;
        }
        else if (l[0] == 1 && l[1] == -1 && l[2] == 1 && l[3] == -1) {
            sse_state->vtx[p] = 3;
        }
        else if (l[0] == -1 && l[1] == 1 && l[2] == 1 && l[3] == -1) {
            sse_state->vtx[p] = 4;
        }
        else if (l[0] == 1 && l[1] == -1 && l[2] == -1 && l[3] == 1) {
            sse_state->vtx[p] = 5;
        }
        else if (l[0] == l[1] && l[1] == l[2] && l[2] == l[3] && l[3] == 1) {
            sse_state->vtx[p] = 6;
        }

        i1 = hberg_system->bond[b][0];
        i2 = hberg_system->bond[b][1];

        v1 = last[i1];
        v2 = last[i2];

        if (v1 != -1) {
            sse_state->link[v1] = v0;
            sse_state->link[v0] = v1;
        } else {
            sse_state->first[i1] = v0;
        }

        if (v2 != -1) {
            sse_state->link[v2] = v0 + 1;
            sse_state->link[v0 + 1] = v2;
        } else {
            sse_state->first[i2] = v0 + 1;
        }

        last[i1] = v0 + 2;
        last[i2] = v0 + 3;
    }

    for (i = 0; i < hberg_system->N; i++) {
        if (sse_state->first[i] != -1) {
            sse_state->link[sse_state->first[i]] = last[i];
            sse_state->link[last[i]] = sse_state->first[i];
        }
    }
}

void free_memory(struct heisenberg_system *hberg_system, struct sse_state *sse_state) {

}


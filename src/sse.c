#include "sse.h"

void simulate_sse(double *beta_vals, int len_beta, long mc_cycles, long therm_cycles, int n_bins, struct heisenberg_system *hberg_system, struct sse_state *sse_state, struct sampled_quantities *samples) {
    int n, t_idx, loop;
    long t;
    double beta;

    init_samples(beta_vals, len_beta, n_bins, samples);

    for (t_idx = 0; t_idx < len_beta; t_idx++) {
        beta = beta_vals[t_idx];
        reset_sse_state(hberg_system, sse_state);

        for (t = 0; t < therm_cycles; t++) {
            diag_update(beta, hberg_system, sse_state);

            for (loop = 0; loop < sse_state->n_loops; loop++) {
                loop_update(hberg_system, sse_state);
            }

            ajust_cutoff(sse_state, t % 1000 == 0);
        }

        for (n = 0; n < n_bins; n++) {
            for (t = 0; t < mc_cycles; t++) {
                diag_update(beta, hberg_system, sse_state);

                sse_state->loop_size = 0;
                for (loop = 0; loop < sse_state->n_loops; loop++) {
                    loop_update(hberg_system, sse_state);
                }

                sample(n, t_idx, hberg_system, sse_state, samples);
            }
        }

        normalize(t_idx, mc_cycles, hberg_system, samples);
        
        printf("beta: %f; n_mean: %f +/- %f; E: %f +/- %f \n", beta, samples->n_mean[t_idx], samples->n_std[t_idx], samples->E_mean[t_idx], samples->E_std[t_idx]);
    }
}

void diag_update(double beta, struct heisenberg_system *hberg_system, struct sse_state *sse_state) {
    int p, b;

    for (p = 0; p < sse_state->M; p++) {
        if (sse_state->op_string[p] == 0) {
            b = next() % hberg_system->Nb;
            if (next_double() <= MIN(1.0, (hberg_system->Nb * beta * prob(b, hberg_system)) / (sse_state->M - sse_state->n))) {
                sse_state->op_string[p] = 2 * (b + 1);
                sse_state->n++;
            }
        }
        else if (sse_state->op_string[p] % 2 == 0) {
            b = (sse_state->op_string[p] / 2) - 1;
            if (next_double() <= MIN(1.0, (sse_state->M - sse_state->n + 1) / (beta * hberg_system->Nb * prob(b, hberg_system)))) {
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
    int j0, j, vtx_type;
    int p, li, le, l, i;
    double r;
    int red_op_string[sse_state->n];
    int trans_op_string[sse_state->n];

    if (sse_state->n == 0) {
        for (i = 0; i < hberg_system->N; i++) {
            if (next_double() <= 0.5) {
                hberg_system->spin[i] = - hberg_system->spin[i];
            }
        }
        return;
    }

    sse_state->vtx = (int *) malloc(sse_state->n * sizeof(int));
    sse_state->link = (int *) malloc(4 * sse_state->n * sizeof(int)); 
    
    create_vtx_list(hberg_system, sse_state, red_op_string, trans_op_string);

    j0 = next() % (4 * sse_state->n);
    j = j0;
    while (true) {
        p = j / 4;
        li = j % 4;
        sse_state->loop_size++;

        r = next_double();
        vtx_type = sse_state->vtx[p] - 1;
        for (le = 0; le < 4; le++) {
            if (r <= sse_state->vtx_type[vtx_type].prob_exit[li][le]) {
                sse_state->vtx[p] = sse_state->vtx_type[vtx_type].new_vtx_type[li][le];
                break;
            }
        }
        
        j = 4 * p + le;
        sse_state->loop_size++;
        if (j == j0) {
            break;
        }
        
        j = sse_state->link[j];
        if (j == j0) {
            break;
        }
    }

    for (p = 0; p < sse_state->n; p++) {
        sse_state->op_string[trans_op_string[p]] = 2 * (red_op_string[p] / 2) + sse_state->vtx_type[sse_state->vtx[p] - 1].type;
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

    free(sse_state->vtx);
    free(sse_state->link);
}

void sample(int n, int t_idx, struct heisenberg_system *hberg_system, struct sse_state *sse_state, struct sampled_quantities *samples) {
    int i;//, p, b;
    // int norm;
    double m = 0.0;
    double m2 = 0.0;
    double ms = 0.0;
    double m2s = 0.0;

    samples->n_bins[t_idx][n] += sse_state->n;
    samples->n2_bins[t_idx][n] += sse_state->n * sse_state->n;

    for (i = 0; i < hberg_system->N; i++) {
        m += hberg_system->spin[i];
        ms += pow(- 1.0, i) * hberg_system->spin[i];
    }
    m *= 0.5;
    ms *= 0.5;
    m2 += m * m;
    m2s += ms * ms;

    // for (p = 0; p < sse_state->M; p++) {
    //     if (sse_state->op_string[p] % 2 == 1) {
    //         b = (sse_state->op_string[p] / 2) - 1;
    //         hberg_system->spin[hberg_system->bond[b][0]] = - hberg_system->spin[hberg_system->bond[b][0]];
    //         hberg_system->spin[hberg_system->bond[b][1]] = - hberg_system->spin[hberg_system->bond[b][1]];

    //         ms += 2 * pow(- 1.0, hberg_system->bond[b][0]) * hberg_system->spin[hberg_system->bond[b][0]];
    //     }

    //     if (sse_state->op_string[p] != 0) {
    //         m2s += ms * ms;
    //     }
    // }

    // norm = sse_state->n > 0 ? sse_state->n : 1;
    m /= hberg_system->N;
    m2 /= hberg_system->N;
    ms /= hberg_system->N;
    m2s /= hberg_system->N;

    samples->m_bins[t_idx][n] += m;
    samples->m2_bins[t_idx][n] += m2;
    samples->ms_bins[t_idx][n] += ms;
    samples->m2s_bins[t_idx][n] += m2s;
}

void normalize(int t_idx, long mc_cycles, struct heisenberg_system *hberg_system, struct sampled_quantities *samples) {
    int n;

    for (n = 0; n < samples->bins; n++) {
        samples->n_bins[t_idx][n] /= mc_cycles;
        samples->n2_bins[t_idx][n] /= mc_cycles;
        samples->E_bins[t_idx][n] = - samples->n_bins[t_idx][n] / (samples->beta_vals[t_idx] * hberg_system->N) + hberg_system->J * hberg_system->C;
        samples->C_bins[t_idx][n] = (samples->n2_bins[t_idx][n] - samples->n_bins[t_idx][n] * samples->n_bins[t_idx][n] - samples->n_bins[t_idx][n]) / hberg_system->N;
        
        samples->m_bins[t_idx][n] /= mc_cycles;
        samples->m2_bins[t_idx][n] /= mc_cycles;
        samples->ms_bins[t_idx][n] /= mc_cycles;
        samples->m2s_bins[t_idx][n] /= mc_cycles;
        samples->m_sus_bins[t_idx][n] = samples->beta_vals[t_idx] * (samples->m2_bins[t_idx][n] - samples->m_bins[t_idx][n] * samples->m_bins[t_idx][n]);
        
        samples->n_mean[t_idx] += samples->n_bins[t_idx][n];
        samples->n2_mean[t_idx] += samples->n2_bins[t_idx][n];
        samples->E_mean[t_idx] += samples->E_bins[t_idx][n];
        samples->C_mean[t_idx] += samples->C_bins[t_idx][n];

        samples->m_mean[t_idx] += samples->m_bins[t_idx][n];
        samples->m2_mean[t_idx] += samples->m2_bins[t_idx][n];
        samples->ms_mean[t_idx] += samples->ms_bins[t_idx][n];
        samples->m2s_mean[t_idx] += samples->m2s_bins[t_idx][n];
        samples->m_sus_mean[t_idx] += samples->m_sus_bins[t_idx][n];
    }
    samples->n_mean[t_idx] /= samples->bins;
    samples->n2_mean[t_idx] /= samples->bins;
    samples->E_mean[t_idx] /= samples->bins;
    samples->C_mean[t_idx] /= samples->bins;

    samples->m_mean[t_idx] /= samples->bins;
    samples->m2_mean[t_idx] /= samples->bins;
    samples->ms_mean[t_idx] /= samples->bins;
    samples->m2s_mean[t_idx] /= samples->bins;
    samples->m_sus_mean[t_idx] /= samples->bins;

    for (n = 0; n < samples->bins; n++) {
        samples->n_std[t_idx] += pow(samples->n_bins[t_idx][n] - samples->n_mean[t_idx], 2.0);
        samples->E_std[t_idx] += pow(samples->E_bins[t_idx][n] - samples->E_mean[t_idx], 2.0);
        samples->C_std[t_idx] += pow(samples->C_bins[t_idx][n] - samples->C_mean[t_idx], 2.0);

        samples->m_std[t_idx] += pow(samples->m_bins[t_idx][n] - samples->m_mean[t_idx], 2.0);
        samples->m2_std[t_idx] += pow(samples->m2_bins[t_idx][n] - samples->m2_mean[t_idx], 2.0);
        samples->ms_std[t_idx] += pow(samples->ms_bins[t_idx][n] - samples->ms_mean[t_idx], 2.0);
        samples->m2s_std[t_idx] += pow(samples->m2s_bins[t_idx][n] - samples->m2s_mean[t_idx], 2.0);
        samples->m_sus_std[t_idx] += pow(samples->m_sus_bins[t_idx][n] - samples->m_sus_mean[t_idx], 2.0);
    }
    samples->n_std[t_idx] = sqrt(samples->n_std[t_idx] / samples->bins);
    samples->E_std[t_idx] = sqrt(samples->E_std[t_idx] / samples->bins);
    samples->C_std[t_idx] = sqrt(samples->C_std[t_idx] / samples->bins);

    samples->m_std[t_idx] = sqrt(samples->m_std[t_idx] / samples->bins);
    samples->m2_std[t_idx] = sqrt(samples->m2_std[t_idx] / samples->bins);
    samples->ms_std[t_idx] = sqrt(samples->ms_std[t_idx] / samples->bins);
    samples->m2s_std[t_idx] = sqrt(samples->m2s_std[t_idx] / samples->bins);
    samples->m_sus_std[t_idx] = sqrt(samples->m_sus_std[t_idx] / samples->bins);
}

void write_to_file(char *filename, struct sampled_quantities *samples) {
    int t_idx;

    FILE *fp;
    fp = fopen(filename, "w");

    fprintf(fp, "beta,n,n2,n_std,E,E_std,C,C_std,m,m_std,m2,m2_std,ms,ms_std,m2s,m2s_std,sus,sus_std\n");
    for (t_idx = 0; t_idx < samples->betas; t_idx++) {
        fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", 
        samples->beta_vals[t_idx], 
        samples->n_mean[t_idx],
        samples->n2_mean[t_idx], 
        samples->n_std[t_idx],
        samples->E_mean[t_idx],
        samples->E_std[t_idx],
        samples->C_mean[t_idx],
        samples->C_std[t_idx],
        samples->m_mean[t_idx],
        samples->m_std[t_idx],
        samples->m2_mean[t_idx],
        samples->m2_std[t_idx],
        samples->ms_mean[t_idx],
        samples->ms_std[t_idx],
        samples->m2s_mean[t_idx],
        samples->m2s_std[t_idx],
        samples->m_sus_mean[t_idx],
        samples->m_sus_std[t_idx]);
    }

    fclose(fp);
}

void ajust_cutoff(struct sse_state *sse_state, bool adjust_loop) {
    int M_new = sse_state->n * 1.33;

    if (M_new > sse_state->M) {
        int opstring_cpy[sse_state->M];
        memcpy(opstring_cpy, sse_state->op_string, sse_state->M * sizeof(int));
        free(sse_state->op_string); 
        sse_state->op_string = (int*) malloc(M_new * sizeof(int));
        memset(sse_state->op_string, 0, M_new * sizeof(int));
        memcpy(sse_state->op_string, opstring_cpy, sse_state->M * sizeof(int));

        sse_state->M = M_new;
    }

    if (adjust_loop) {
        if (sse_state->loop_size > 1000) {
            sse_state->n_loops = 2000 * sse_state->M * sse_state->n_loops / sse_state->loop_size;
        }
        sse_state->loop_size = 0;
    }
}

void create_vtx_list(struct heisenberg_system *hberg_system, struct sse_state *sse_state, int *red_op_string, int *trans_op_string) {
    int p, i, v0, b;
    int i1, i2, v1, v2;
    int last[hberg_system->N];
    int l[4];
    int p_red;

    memset(last, -1, hberg_system->N * sizeof(int));
    memset(sse_state->first, -1, hberg_system->N * sizeof(int));
    memset(sse_state->link, -1, 4 * sse_state->n * sizeof(int));
    memset(sse_state->vtx, 0, sse_state->n * sizeof(int));
    memset(red_op_string, 0, sse_state->n * sizeof(int));
    memset(trans_op_string, 0, sse_state->n * sizeof(int));
    
    p_red = 0;
    for (p = 0; p < sse_state->M; p++) {
        if (sse_state->op_string[p] != 0) {
            red_op_string[p_red] = sse_state->op_string[p];
            trans_op_string[p_red] = p;
            p_red++;
        }
    }

    for (p = 0; p < sse_state->n; p++) {
        b = (red_op_string[p] / 2) - 1;
        v0 = 4 * p;

        i1 = hberg_system->bond[b][0];
        i2 = hberg_system->bond[b][1];

        l[0] = hberg_system->spin[i1];
        l[1] = hberg_system->spin[i2];
        if (red_op_string[p] % 2 != 0) {
            hberg_system->spin[i1] = - hberg_system->spin[i1];
            hberg_system->spin[i2] = - hberg_system->spin[i2];
        }
        l[2] = hberg_system->spin[i1];
        l[3] = hberg_system->spin[i2];

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

double prob(int b, struct heisenberg_system *hberg_system) {
    if (hberg_system->spin[hberg_system->bond[b][0]] != hberg_system->spin[hberg_system->bond[b][1]]) {
        return hberg_system->H[1][1].value;
    }
    else if (hberg_system->spin[hberg_system->bond[b][0]] == hberg_system->spin[hberg_system->bond[b][1]]) {
        if (hberg_system->spin[hberg_system->bond[b][1]] == 1) {
            return hberg_system->H[0][0].value;
        }
        else {
            return hberg_system->H[3][3].value;
        }
    }

    return 0.0;
}

void init_heisenberg_system(int d, int L, double J, double delta, double h, double epsilon, struct heisenberg_system *hberg_system) {
    int i, j;

    hberg_system->d = d;
    hberg_system->L = L;
    hberg_system->N = pow(L, d);
    hberg_system->Nb = hberg_system->N * d;         // For PBC
    
    hberg_system->J = J;
    hberg_system->h = h;
    hberg_system->delta = delta;
    hberg_system->epsilon = epsilon;
    hberg_system->C = 0.25 * delta + 0.5 * h / J + epsilon;
    hberg_system->H = create_hamiltonian(J, delta, h, hberg_system->C);

    hberg_system->spin = (int *) malloc(hberg_system->N * sizeof(int));
    hberg_system->bond = (int **) malloc(hberg_system->Nb * sizeof(int *));
    for (i = 0; i < hberg_system->Nb; i++) {
        hberg_system->bond[i] = (int*) malloc(2 * sizeof(int));
    }
    if (d == 1) {
        for (i = 0; i < L; i++) {
            hberg_system->bond[i][0] = i;
            hberg_system->bond[i][1] = (i + 1) % hberg_system->N;
        }
    } else if (d == 2) {
        for (i = 0; i < L; i++) {
            for (j = 0; j < L; j++) {
                hberg_system->bond[j * L + i][0] = j * L + i;
                hberg_system->bond[j * L + i][1] = j * L + (i + 1) % L;

                hberg_system->bond[(j * L + i) + hberg_system->N][0] = j * L + i;
                hberg_system->bond[(j * L + i) + hberg_system->N][1] = ((j + 1) % L) * L + i;
            }
        }
    }

    for (i = 0; i < hberg_system->N; i++) {
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

    sse_state->n_loops = MAX(4, hberg_system->N / 4);
    sse_state->loop_size = 0;
    sse_state->first = (int *) malloc(hberg_system->N * sizeof(int));
}

void init_samples(double *beta_vals, int len_beta, int n_bins, struct sampled_quantities *samples) {
    int i;

    samples->bins = n_bins;
    samples->betas = len_beta;
    samples->beta_vals = beta_vals;

    samples->n_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->n2_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->E_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->C_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m2_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->ms_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m2s_bins = (double **) malloc(len_beta * sizeof(double *));
    samples->m_sus_bins = (double **) malloc(len_beta * sizeof(double *));
    for (i = 0; i < len_beta; i++) {
        samples->n_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->n2_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->E_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->C_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m2_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->ms_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m2s_bins[i] = (double *) malloc(n_bins * sizeof(double));
        samples->m_sus_bins[i] = (double *) malloc(n_bins * sizeof(double));

        memset(samples->n_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->n2_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->E_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->C_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m2_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->ms_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m2s_bins[i], 0.0, n_bins * sizeof(double));
        memset(samples->m_sus_bins[i], 0.0, n_bins * sizeof(double));
    }

    samples->n_mean = (double *) malloc(len_beta *sizeof(double));
    samples->n_std = (double *) malloc(len_beta *sizeof(double));
    samples->n2_mean = (double *) malloc(len_beta *sizeof(double));

    samples->E_mean = (double *) malloc(len_beta *sizeof(double));
    samples->E_std = (double *) malloc(len_beta *sizeof(double));
    samples->C_mean = (double *) malloc(len_beta *sizeof(double));
    samples->C_std = (double *) malloc(len_beta *sizeof(double));

    samples->m_mean = (double *) malloc(len_beta *sizeof(double));
    samples->m_std = (double *) malloc(len_beta *sizeof(double));
    samples->m2_mean = (double *) malloc(len_beta *sizeof(double));
    samples->m2_std = (double *) malloc(len_beta *sizeof(double));
    samples->ms_mean = (double *) malloc(len_beta *sizeof(double));
    samples->ms_std = (double *) malloc(len_beta *sizeof(double));
    samples->m2s_mean = (double *) malloc(len_beta *sizeof(double));
    samples->m2s_std = (double *) malloc(len_beta *sizeof(double));
    samples->m_sus_mean = (double *) malloc(len_beta *sizeof(double));
    samples->m_sus_std = (double *) malloc(len_beta *sizeof(double));

    memset(samples->n_mean, 0.0, len_beta *sizeof(double));
    memset(samples->n_std, 0.0, len_beta *sizeof(double));
    memset(samples->n2_mean, 0.0, len_beta *sizeof(double));
    memset(samples->E_mean, 0.0, len_beta *sizeof(double));
    memset(samples->E_std, 0.0, len_beta *sizeof(double));
    memset(samples->C_mean, 0.0, len_beta *sizeof(double));
    memset(samples->C_std, 0.0, len_beta *sizeof(double));
    memset(samples->m_mean, 0.0, len_beta *sizeof(double));
    memset(samples->m_std, 0.0, len_beta *sizeof(double));
    memset(samples->m2_mean, 0.0, len_beta *sizeof(double));
    memset(samples->m2_std, 0.0, len_beta *sizeof(double));
    memset(samples->ms_mean, 0.0, len_beta *sizeof(double));
    memset(samples->ms_std, 0.0, len_beta *sizeof(double));
    memset(samples->m2s_mean, 0.0, len_beta *sizeof(double));
    memset(samples->m2s_std, 0.0, len_beta *sizeof(double));
    memset(samples->m_sus_mean, 0.0, len_beta *sizeof(double));
    memset(samples->m_sus_std, 0.0, len_beta *sizeof(double));
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
    sse_state->loop_size = 0;
    sse_state->n_loops = MAX(4, hberg_system->N / 4);

    free(sse_state->op_string);
    sse_state->op_string = (int *) malloc(sse_state->M * sizeof(int));
    memset(sse_state->op_string, 0, sse_state->M * sizeof(int));
}

void free_memory(struct heisenberg_system *hberg_system, struct sse_state *sse_state, struct sampled_quantities *samples) {
    int i;

    for (i = 0; i < 4; i++) {
        free(hberg_system->H[i]);
    }
    for (i = 0; i < hberg_system->N; i++) {
        free(hberg_system->bond[i]);
    }
    free(hberg_system->H);
    free(hberg_system->bond);
    free(hberg_system->spin);

    free(sse_state->first);
    free(sse_state->vtx_type);
    free(sse_state->op_string);

    for (i = 0; i < samples->betas; i++) {
        free(samples->n_bins[i]);
        free(samples->n2_bins[i]);
        free(samples->E_bins[i]);
        free(samples->C_bins[i]);
        free(samples->m_bins[i]);
        free(samples->m2_bins[i]);
        free(samples->ms_bins[i]);
        free(samples->m2s_bins[i]);
        free(samples->m_sus_bins[i]);
    }
    free(samples->n_bins);
    free(samples->n2_bins);
    free(samples->E_bins);
    free(samples->C_bins);
    free(samples->m_bins);
    free(samples->m2_bins);
    free(samples->ms_bins);
    free(samples->m2s_bins);
    free(samples->m_sus_bins);
    free(samples->n_mean);
    free(samples->n_std);
    free(samples->n2_mean);
    free(samples->E_mean);
    free(samples->E_std);
    free(samples->C_mean);
    free(samples->C_std);
    free(samples->m_mean);
    free(samples->m_std);
    free(samples->m2_mean);
    free(samples->m2_std);
    free(samples->ms_mean);
    free(samples->ms_std);
    free(samples->m2s_mean);
    free(samples->m2s_std);
    free(samples->m_sus_mean);
    free(samples->m_sus_std);
}

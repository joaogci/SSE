#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <omp.h>

#include "sse.h"

double mean(int *arr, int size) {
    double mean = 0;
    int i;
    for (i = 0; i < size; i++) {
        mean += arr[i];
    }
    return mean / size;
}

int main(int argc, char **argv) {
    int d = 1;
    int L = 2;
    double beta;
    double beta_vals[6] = {0.5, 1.0, 2.0, 4.0, 8.0, 16.0};
    // double beta_vals[50];
    // int beta_len = sizeof(beta_vals) / sizeof(beta_vals[0]);
    int beta_len = 6;

    // for (int i = 0; i < beta_len; i++) {
    //    beta_vals[i] = 1.0 / ((2.0 - 0.0) * (i + 1) / beta_len);
    // }

    long therm_cycles = 1e4;
    long mc_cycles = 1e4;

    struct sse_state *sse_state = malloc(sizeof(struct sse_state));
    struct vertex *vertex = malloc(sizeof(struct vertex));

    long t;
    int n_vals[mc_cycles];
    int n2_vals[mc_cycles];
    double n_mean[beta_len];
    double n2_mean[beta_len];
    double E_mean[beta_len];
    double C_mean[beta_len];
    double ms_mean[beta_len];
    double m2s_mean[beta_len];
    double m_mean[beta_len];
    double m2_mean[beta_len];

    double sus_mean[beta_len];
    double suss_mean[beta_len];

    for (int i = 0; i < beta_len; i++) {
        beta = beta_vals[i];
        n_mean[i] = 0.0;
        n2_mean[i] = 0.0;
        E_mean[i] = 0.0;
        C_mean[i] = 0.0;
        ms_mean[i] = 0.0;
        m2s_mean[i] = 0.0;
        m_mean[i] = 0.0;
        m2_mean[i] = 0.0;

        sus_mean[i] = 0.0;
        suss_mean[i] = 0.0;

        init_system(d, L, beta, sse_state, vertex, (uint64_t) time(NULL));

        for (t = 0; t < therm_cycles; t++) {
            diag_update(sse_state);
            create_vertex_list(sse_state, vertex);
            loop_update(sse_state, vertex);
            ajust_cutoff(sse_state, vertex);
        }

        for (t = 0; t < mc_cycles; t++) {
            diag_update(sse_state);
            create_vertex_list(sse_state, vertex);
            loop_update(sse_state, vertex);

            n_vals[t] = sse_state->n;
            n2_vals[t] = sse_state->n * sse_state->n;
            
            // Sampling magnetization
            double m = 0.0;
            double m2 = 0.0;
            double ms = 0.0;
            double m2s = 0.0;
            for (int j = 0; j < sse_state->N; j++) {
                ms += pow(- 1.0, j) * sse_state->spin[j];
                m += sse_state->spin[j];
            }
            m *= 0.5;
            ms *= 0.5;

            int b;
            for (int p = 0; p < sse_state->M; p++) {
                if (sse_state->opstring[p] % 2 == 1) {
                    b = (sse_state->opstring[p] / 2) - 1;
                    sse_state->spin[sse_state->bond_site[b][0]] = - sse_state->spin[sse_state->bond_site[b][0]];
                    sse_state->spin[sse_state->bond_site[b][1]] = - sse_state->spin[sse_state->bond_site[b][1]];

                    ms += 2 * pow(- 1.0, sse_state->bond_site[b][0]) * sse_state->spin[sse_state->bond_site[b][0]];
                    m += 2 * sse_state->spin[sse_state->bond_site[b][0]];
                }

                if (sse_state->opstring[p] != 0) {
                    m2s += ms * ms;
                    m2 += m * m;
                }
            }
            if (sse_state->n != 0) {
                m2 /= (sse_state->n * sse_state->N * sse_state->N);
                m /= (sse_state->n * sse_state->N);
                m2s /= (sse_state->n * sse_state->N * sse_state->N);
                ms /= (sse_state->n * sse_state->N);
            }
            else {
                m2 /= (sse_state->N * sse_state->N);
                m /= (sse_state->N);
                m2s /= (sse_state->N * sse_state->N);
                ms /= (sse_state->N);
            }

            m_mean[i] += m;
            m2_mean[i] += m2;
            ms_mean[i] += ms;
            m2s_mean[i] += m2s;

            sus_mean[i] += m2;
            suss_mean[i] += m2s;
        }

        n_mean[i] = mean(n_vals, mc_cycles);
        n2_mean[i] = mean(n2_vals, mc_cycles);
        E_mean[i] = - n_mean[i] / (sse_state->N * beta) + 0.25;
        C_mean[i] = (n2_mean[i] - n_mean[i] * n_mean[i] - n_mean[i]) / sse_state->N;
        m_mean[i] /= mc_cycles;
        m2_mean[i] /= mc_cycles;
        ms_mean[i] /= mc_cycles;
        m2s_mean[i] /= mc_cycles;

        sus_mean[i] = beta * sus_mean[i] / mc_cycles;
        suss_mean[i] = beta * suss_mean[i] / mc_cycles;

        printf("beta: %f | E: %f | C: %f | m: %f | m2: %f | m_s: %f | m2_s: %f | sus: %f | suss: %f | n %f \n", beta, E_mean[i], C_mean[i], m_mean[i], m2_mean[i], ms_mean[i], m2s_mean[i], sus_mean[i], suss_mean[i], n_mean[i]);
    }

    free_state(&sse_state);
    free_vertex(&vertex);

    // Write to file
    FILE *fp;
    fp = fopen("2D_heisenberg_L16_more_T.csv", "w");

    fprintf(fp, "beta,E,C,m,m2,m_s,m2_s,sus,suss,n\n");
    for (int i = 0; i < beta_len; i++) {
        fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", beta_vals[i], E_mean[i], C_mean[i], m_mean[i], m2_mean[i], ms_mean[i], m2s_mean[i], sus_mean[i], suss_mean[i], n_mean[i]);
    }

    fclose(fp);

    return 0;
}

#include "io.h"

#define BUFFER_SIZE 256

void read_inputs(char *file_name, int *d, int *L, double *J, double *delta, double *h, double *epsilon, long *therm_cycles, long *mc_cycles, int *n_bins, double **beta_vals, int *len_beta) 
{
    FILE *input_file;
    char buffer[BUFFER_SIZE];
    int count = 0;
    double Ti, Tf;

    input_file = fopen(file_name, "r");

    if (input_file != NULL) {
        while (fgets(buffer, BUFFER_SIZE, input_file) != NULL) {
            if (strcmp(buffer, "test_temps\n") == 0) {
                (*len_beta) = 6;
                (*beta_vals) = (double *) malloc((*len_beta) * sizeof(double));
                for (int i = 0; i < (*len_beta); i++) {
                    (*beta_vals)[i] = pow(2.0, i - 1.0);
                }
                break;
            } else if (buffer[0] != '#' && buffer[0] != '\n') {
                switch (count) {
                case 0:
                    sscanf(buffer, "%d, %d, %lf, %lf, %lf, %lf ", d, L, J, delta, h, epsilon);
                    break;
                case 1:
                    sscanf(buffer, "%ld, %ld, %d ", therm_cycles, mc_cycles, n_bins);
                    break;
                case 2:
                    sscanf(buffer, "%d ", len_beta);
                    break;
                case 3:
                    sscanf(buffer, "%lf, %lf ", &Ti, &Tf);
                    break;
                }
                count++;
            }
        }

        if (count == 4) {
            (*beta_vals) = (double *) malloc((*len_beta) * sizeof(double));
            for (int i = 0; i < (*len_beta); i++) {
                (*beta_vals)[i] = 1.0 / ((Tf - Ti) * (i + 1) / (*len_beta));
            }
        }
    } else {
        printf("Error opening the %s file. Check if the file exists. \n", file_name);
        exit(1);
    }

    fclose(input_file);
}

void write_outputs(char *file_name, sampled_quantities *samples) 
{
    FILE *output_file;

    output_file = fopen(file_name, "w");

    if (output_file != NULL) {
        fprintf(output_file, "beta,n,n2,n_std,E,E_std,C,C_std,m,m_std,m2,m2_std,m4,m4_std,ms,ms_std,m2s,m2s_std,m4s,m4s_std,sus,sus_std,binder,binder_std,binders,binders_std\n");
        for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
            fprintf(output_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", 
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
            samples->m4_mean[t_idx],
            samples->m4_std[t_idx],
            samples->ms_mean[t_idx],
            samples->ms_std[t_idx],
            samples->m2s_mean[t_idx],
            samples->m2s_std[t_idx],
            samples->m4s_mean[t_idx],
            samples->m4s_std[t_idx],
            samples->m_sus_mean[t_idx],
            samples->m_sus_std[t_idx],
            samples->binder_mean[t_idx],
            samples->binder_std[t_idx],
            samples->binders_mean[t_idx],
            samples->binders_std[t_idx]);
        }
    } else {
        printf("Error in opening the %s file. \n", file_name);
        exit(1);
    }

    fclose(output_file);
}


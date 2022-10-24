#include "io.h"

#define BUFFER_SIZE 256

void read_inputs(char *file_name, int *d, int *L, double *S, double *delta, double *h, double *epsilon, long *therm_cycles, long *mc_cycles, int *n_bins, double **beta_vals, int *len_beta) 
{
    char buffer[BUFFER_SIZE];
    char S_eval[BUFFER_SIZE];
    int count = 0;
    double Ti, Tf;

    FILE *input_file;
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
                    sscanf(buffer, "%d, %d, %s, %lf, %lf, %lf ", d, L, S_eval, delta, h, epsilon);
                    for (int i = 0; i < sizeof(S_eval)/sizeof(S_eval[0]); i++) { 
                        if (S_eval[i] == ',') {
                            S_eval[i] = '\0'; 
                            break;
                        } 
                    }
                    (*S) = te_interp(S_eval, 0);
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

void read_vtx_info(char *file_name, vtx_element **vtx, int *n_diagrams, int *n_updates, int *n_legs) 
{
    FILE *vtx_file;
    vtx_file = fopen(file_name, "r");

    if (vtx_file != NULL) {
        fscanf(vtx_file, "%d %d %d \n", n_diagrams, n_updates, n_legs);
        (*vtx) = (vtx_element *) malloc((*n_diagrams) * sizeof(vtx_element));

        for (int i = 0; i < (*n_diagrams); i++) {
            (*vtx)[i].spin = (int *) malloc((*n_legs) * sizeof(int));
            (*vtx)[i].new_vtx_type_ = (int ***) malloc((*n_updates) * sizeof(int**));
            (*vtx)[i].prob_exit_ = (double ***) malloc((*n_updates) * sizeof(double**));
            for (int j = 0; j < (*n_updates); j++) {
                (*vtx)[i].new_vtx_type_[j] = (int **) malloc((*n_legs) * sizeof(int *));
                (*vtx)[i].prob_exit_[j] = (double **) malloc((*n_legs) * sizeof(double *));
                for (int l = 0; l < (*n_legs); l++) {
                    (*vtx)[i].new_vtx_type_[j][l] = (int *) malloc((*n_legs) * sizeof(int));
                    (*vtx)[i].prob_exit_[j][l] = (double *) malloc((*n_legs) * sizeof(double));
                }
            }

            fscanf(vtx_file, "%d \n", &((*vtx)[i].indx));
            fscanf(vtx_file, "%d \n", &((*vtx)[i].type));
            fscanf(vtx_file, "%lf \n", &((*vtx)[i].H));
            fscanf(vtx_file, "%d %d %d %d \n", &((*vtx)[i].spin[0]), &((*vtx)[i].spin[1]), &((*vtx)[i].spin[2]), &((*vtx)[i].spin[3]));
            
            for (int j = 0; j < (*n_updates); j++) {
                for (int l = 0; l < N_LEGS; l++) {
                    fscanf(vtx_file, "%d %d %d %d \n", &((*vtx)[i].new_vtx_type_[j][l][0]), &((*vtx)[i].new_vtx_type_[j][l][1]), &((*vtx)[i].new_vtx_type_[j][l][2]), &((*vtx)[i].new_vtx_type_[j][l][3]));
                }

                for (int l = 0; l < N_LEGS; l++) {
                    fscanf(vtx_file, "%lf %lf %lf %lf \n", &((*vtx)[i].prob_exit_[j][l][0]), &((*vtx)[i].prob_exit_[j][l][1]), &((*vtx)[i].prob_exit_[j][l][2]), &((*vtx)[i].prob_exit_[j][l][3]));
                }
            }
        }
    } else {
        printf("Error opening the %s file. \n", file_name);
        exit(1);
    }

    fclose(vtx_file);
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

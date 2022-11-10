#include "io.h"

#define BUFFER_SIZE 256

/*
 * function: read_inputs
 *  reads input file for the simulation and stores them in the 
 * parameters of the function
 *  
 *  parameters:
 *      (int *) d: dimension
 *      (int *) L: length of the system
 *      (double *) S: spin quantum number
 *      (double *) delta: anisotropy
 *      (double *) h: applied magnetic field
 *      (double *) epsilon: constant added to the Hamiltonian
 *      (long *) therm_cycles: MCS for thermalization
 *      (long *) mc_cycles: MCS for sampling
 *      (int *) n_bins: number of sampling bins
 *      (double **) beta_vals: array to store simulation temperatures
 *      (int *) len_betas: number of temperatures 
 */
void read_inputs(int *d, int *L, double *S, 
    double *delta, double *h, double *epsilon, long *therm_cycles, 
    long *mc_cycles, int *n_bins, double **beta_vals, int *len_beta) 
{
    char buffer[BUFFER_SIZE];
    int count = 0;
    double Ti, Tf;
    FILE *input_file;

    input_file = fopen("read.in", "r");
    if (input_file != NULL) {
        fgets(buffer, BUFFER_SIZE, input_file);
        sscanf(buffer, "%d, %d, %lf, %lf, %lf, %lf", d, L, S, delta, h, epsilon);
        
        fgets(buffer, BUFFER_SIZE, input_file);
        sscanf(buffer, "%ld, %ld, %d ", therm_cycles, mc_cycles, n_bins);
    } else {
        printf("Error opening the read.i file. Check if the file exists. \n");
        exit(1);
    }
    fclose(input_file);

    input_file = fopen("beta.in", "r");
    if (input_file != NULL) {
        fgets(buffer, BUFFER_SIZE, input_file);
        sscanf(buffer, "%d ", len_beta);
        (*beta_vals) = (double *) malloc((*len_beta) * sizeof(double));
        
        for (int i = 0; i < (*len_beta); i++) {
            fgets(buffer, BUFFER_SIZE, input_file);
            sscanf(buffer, "%lf ", &((*beta_vals)[i]));
        }
    } else {
        printf("Error opening the beta.in file. Check if the file exists. \n");
        exit(1);
    }
    fclose(input_file);
}

/*
 * function: read_vtx_info
 *  reads vertex information for the simulation and 
 * stores it in the parameters of the function
 *  
 *  parameters:
 *      (char *) file_name: file name
 *      (vtx_element **) d: dimension
 *      (int *) n_diagrams: number of diagrams
 */
void read_vtx_info(char *file_name, vtx_element **vtx, int *n_diagrams) 
{
    FILE *vtx_file;
    vtx_file = fopen(file_name, "r");

    if (vtx_file != NULL) {
        fscanf(vtx_file, "%d \n", n_diagrams);
        (*vtx) = (vtx_element *) malloc((*n_diagrams) * sizeof(vtx_element));

        for (int i = 0; i < (*n_diagrams); i++) {
            fscanf(vtx_file, "%d \n", &((*vtx)[i].indx));
            fscanf(vtx_file, "%d \n", &((*vtx)[i].type));
            fscanf(vtx_file, "%lf \n", &((*vtx)[i].H));
            fscanf(vtx_file, "%d %d %d %d \n", &((*vtx)[i].spin[0]), &((*vtx)[i].spin[1]), &((*vtx)[i].spin[2]), &((*vtx)[i].spin[3]));
            
            for (int j = 0; j < N_UPDATES; j++) {
                for (int l = 0; l < N_LEGS; l++) {
                    fscanf(vtx_file, "%d %d %d %d \n", &((*vtx)[i].new_vtx_type[j][l][0]), &((*vtx)[i].new_vtx_type[j][l][1]), &((*vtx)[i].new_vtx_type[j][l][2]), &((*vtx)[i].new_vtx_type[j][l][3]));
                }

                for (int l = 0; l < N_LEGS; l++) {
                    fscanf(vtx_file, "%lf %lf %lf %lf \n", &((*vtx)[i].prob_exit[j][l][0]), &((*vtx)[i].prob_exit[j][l][1]), &((*vtx)[i].prob_exit[j][l][2]), &((*vtx)[i].prob_exit[j][l][3]));
                }
            }
        }
    } else {
        printf("Error opening the %s file. \n", file_name);
        exit(1);
    }

    fclose(vtx_file);
}

/*
 * function: write_outputs
 *  writes the simulation outputs to file
 *  
 *  parameters:
 *      (char *) file_name: file name
 *      (sampled_quantities *) samples: sampled quantities
 * during the simulation
 */
void write_outputs(char *file_name, sampled_quantities *samples, 
    int d, int L, double S, double delta, double h, double epsilon,
    long therm_cycles, long mc_cycles, double cpu_time_used, int n_threads) 
{
    FILE *output_file;
    output_file = fopen(file_name, "w");

    if (output_file != NULL) {
        fprintf(output_file, "d,L,S,delta,h,epsilon\n");
        fprintf(output_file, "%d,%d,%lf,%lf,%lf,%lf\n", d, L, S, delta, h, epsilon);
        
        fprintf(output_file, "therm_cycles,mc_cycles,n_bins\n");
        fprintf(output_file, "%ld,%ld,%d \n", therm_cycles, mc_cycles, samples->bins);
        
        fprintf(output_file, "cpu_time,n_threads\n");
        fprintf(output_file, "%lf,%d\n", cpu_time_used, n_threads);

        fprintf(output_file, "n_betas\n");
        fprintf(output_file, "%d\n", samples->betas);

        fprintf(output_file, "beta,n,n2,n_std,E,E_std,C,C_std,m,m_std,m2,m2_std,m4,m4_std,ms,ms_std,m2s,m2s_std,m4s,m4s_std,sus,sus_std,S_mean,S_std\n");
        for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
            fprintf(output_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", 
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
            samples->S_mean[t_idx],
            samples->S_std[t_idx]);
        }

        
        for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
            fprintf(output_file, "beta\n");
            fprintf(output_file, "%lf\n", samples->beta_vals[t_idx]);
            fprintf(output_file, "corr_mean,corr_std\n");
            for (int i = 0; i < L; i++) {
                fprintf(output_file, "%lf,%lf\n", samples->corr_mean[t_idx][i], samples->corr_std[t_idx][i]);
            }
        }
    } else {
        printf("Error in opening the %s file. \n", file_name);
        exit(1);
    }

    fclose(output_file);
}

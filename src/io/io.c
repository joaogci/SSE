#include "io.h"

void read_vtx_info(Vertices** vtx, int* n_diagrams)
{
  FILE *vtx_file;
  vtx_file = fopen("vertices_info", "r");

  if (vtx_file != NULL) {
    fscanf(vtx_file, "%d \n", n_diagrams);
    (*vtx) = (Vertices*) malloc((*n_diagrams) * sizeof(Vertices));

    for (int i = 0; i < (*n_diagrams); i++) {
      fscanf(vtx_file, "%d \n", &((*vtx)[i].indx));
      fscanf(vtx_file, "%d \n", &((*vtx)[i].type));
      fscanf(vtx_file, "%lf \n", &((*vtx)[i].weight));
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
    printf("Error opening the vertices_info file. \n");
    exit(1);
  }

  fclose(vtx_file);
}

void write_observables(Obs_scalar* obs, int n_scal)
{
  int i;
  FILE* out;
  char filename[BUFFER_SIZE];

  for (i = 0; i < n_scal; i++) {
    strcpy(filename, obs[i].filename);
    // strcpy(filename, "../");
    // strcat(filename, obs[i].filename);
    strcat(filename, "_scal");

    out = fopen(filename, "a");
    write_obs_scalar(out, &(obs[i]));
    fclose(out);
  }
}



// void read_vtx_info(char *file_name, vtx_element **vtx, int *n_diagrams) 
// {
//     FILE *vtx_file;
//     vtx_file = fopen(file_name, "r");

//     if (vtx_file != NULL) {
//         fscanf(vtx_file, "%d \n", n_diagrams);
//         (*vtx) = (vtx_element *) malloc((*n_diagrams) * sizeof(vtx_element));

//         for (int i = 0; i < (*n_diagrams); i++) {
//             fscanf(vtx_file, "%d \n", &((*vtx)[i].indx));
//             fscanf(vtx_file, "%d \n", &((*vtx)[i].type));
//             fscanf(vtx_file, "%lf \n", &((*vtx)[i].H));
//             fscanf(vtx_file, "%d %d %d %d \n", &((*vtx)[i].spin[0]), &((*vtx)[i].spin[1]), &((*vtx)[i].spin[2]), &((*vtx)[i].spin[3]));
            
//             for (int j = 0; j < N_UPDATES; j++) {
//                 for (int l = 0; l < N_LEGS; l++) {
//                     fscanf(vtx_file, "%d %d %d %d \n", &((*vtx)[i].new_vtx_type[j][l][0]), &((*vtx)[i].new_vtx_type[j][l][1]), &((*vtx)[i].new_vtx_type[j][l][2]), &((*vtx)[i].new_vtx_type[j][l][3]));
//                 }

//                 for (int l = 0; l < N_LEGS; l++) {
//                     fscanf(vtx_file, "%lf %lf %lf %lf \n", &((*vtx)[i].prob_exit[j][l][0]), &((*vtx)[i].prob_exit[j][l][1]), &((*vtx)[i].prob_exit[j][l][2]), &((*vtx)[i].prob_exit[j][l][3]));
//                 }
//             }
//         }
//     } else {
//         printf("Error opening the %s file. \n", file_name);
//         exit(1);
//     }

//     fclose(vtx_file);
// }

// /*
//  * function: write_outputs
//  *  writes the simulation outputs to file
//  *  
//  *  parameters:
//  *      (sampled_quantities *) samples: sampled quantities
//  * during the simulation
//  *  returns:
//  *      (char *) file_name: name of save file
//  */
// char *write_outputs(sampled_quantities *samples, 
//     int d, int L, int boundary_cond, double S, double delta, double h, double epsilon,
//     long therm_cycles, long mc_cycles, double cpu_time_used, int n_threads) 
// {
//     char *file_name = (char *) malloc(BUFFER_SIZE * sizeof(char));
//     char *buffer = (boundary_cond == 0) ? "PBC" : "OBC";

// #ifndef KINETIC
//     if (delta > 0) {
//         sprintf(file_name, "%dD_L%d_%s_AFM_XXZ_S%g_delta%g_h%g_ep%g.csv", d, L, buffer, S, fabs(delta), h, epsilon);
//     } else if (delta < 0) {
//         sprintf(file_name, "%dD_L%d_%s_FM_XXZ_S%g_delta%g_h%g_ep%g.csv", d, L, buffer, S, fabs(delta), h, epsilon);
//     } else {
//         sprintf(file_name, "%dD_L%d_%s_XY_S%g_h%g_ep%g.csv", d, L, buffer, S, h, epsilon);
//     }
// #else
//     if (delta > 0) {
//         sprintf(file_name, "%dD_L%d_%s_AFM_XXZ_S%g_delta%g_h%g_ep%g_x%d_y%d.csv", d, L, buffer, S, fabs(delta), h, epsilon, samples->x, samples->y);
//     } else if (delta < 0) {
//         sprintf(file_name, "%dD_L%d_%s_FM_XXZ_S%g_delta%g_h%g_ep%g_x%d_y%d.csv", d, L, buffer, S, fabs(delta), h, epsilon, samples->x, samples->y);
//     } else {
//         sprintf(file_name, "%dD_L%d_%s_XY_S%g_h%g_ep%g_x%d_y%d.csv", d, L, buffer, S, h, epsilon, samples->x, samples->y);
//     }
// #endif // KINETIC

//     FILE *output_file;
//     output_file = fopen(file_name, "w");

//     if (output_file != NULL) {
//         fprintf(output_file, "d,L,boundary_cond,S,delta,h,epsilon\n");
//         fprintf(output_file, "%d,%d,%s,%lf,%lf,%lf,%lf\n", d, L, buffer, S, fabs(delta), h, epsilon);
        
//         fprintf(output_file, "therm_cycles,mc_cycles,n_bins\n");
//         fprintf(output_file, "%ld,%ld,%d \n", therm_cycles, mc_cycles, samples->bins);
        
//         fprintf(output_file, "cpu_time,n_threads\n");
//         fprintf(output_file, "%lf,%d\n", cpu_time_used, n_threads);

// #if defined(L_SS) && defined(L_HH) && defined(L_SH)
//         char cond[5] = "full";
// #else 
//     #if defined(L_SS) && defined(L_HH)
//         char cond[5] = "diag";
//     #elif defined(L_SH)
//         char cond[5] = "offd";
//     #elif defined(L_HH) && defined(L_SH)
//         char cond[5] = "hhsh"; 
//     #elif defined(L_SS) && defined(L_SH)
//         char cond[5] = "sssh";
//     #elif defined(L_HH)
//         char cond[5] = "hh";
//     #elif defined(L_SS)
//         char cond[5] = "ss";
//     #else 
//         char cond[5] = "."; 
//     #endif
// #endif

//         fprintf(output_file, "n_betas,n_k,x,y,kinetic\n");
//         fprintf(output_file, "%d,%d,%d,%d,%s\n", samples->betas, samples->k_max, samples->x, samples->y, cond);

//         fprintf(output_file, "beta,n,n2,n_std,E,E_std,C,C_std,m,m_std,m2,m2_std,m4,m4_std,ms,ms_std,m2s,m2s_std,m4s,m4s_std,sus,sus_std,S_mean,S_std\n");
//         for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
//             fprintf(output_file, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", 
//             samples->beta_vals[t_idx], 
//             samples->n_mean[t_idx],
//             samples->n2_mean[t_idx], 
//             samples->n_std[t_idx],
//             samples->E_mean[t_idx],
//             samples->E_std[t_idx],
//             samples->C_mean[t_idx],
//             samples->C_std[t_idx],
//             samples->m_mean[t_idx],
//             samples->m_std[t_idx],
//             samples->m2_mean[t_idx],
//             samples->m2_std[t_idx],
//             samples->m4_mean[t_idx],
//             samples->m4_std[t_idx],
//             samples->ms_mean[t_idx],
//             samples->ms_std[t_idx],
//             samples->m2s_mean[t_idx],
//             samples->m2s_std[t_idx],
//             samples->m4s_mean[t_idx],
//             samples->m4s_std[t_idx],
//             samples->m_sus_mean[t_idx],
//             samples->m_sus_std[t_idx],
//             samples->S_mean[t_idx],
//             samples->S_std[t_idx]);
//         }
        
//         for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
//             fprintf(output_file, "beta\n");
//             fprintf(output_file, "%lf\n", samples->beta_vals[t_idx]);
//             fprintf(output_file, "corr_mean,corr_std\n");
//             for (int i = 0; i < L; i++) {
//                 fprintf(output_file, "%lf,%lf\n", samples->corr_mean[t_idx][i], samples->corr_std[t_idx][i]);
//             }
//         }
// #ifdef L_SS
//         for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
//             fprintf(output_file, "beta\n");
//             fprintf(output_file, "%lf\n", samples->beta_vals[t_idx]);
//             fprintf(output_file, "w_k,L_SS_mean,L_SS_std\n");
//             for (int k = 0; k < samples->k_max; k++) {
//                 fprintf(output_file, "%lf,%lf,%lf\n", samples->w_k[t_idx][k], samples->L_SS_mean[t_idx][k], samples->L_SS_std[t_idx][k]);
//             }
//         }
// #endif // L_SS
// #ifdef L_HH
//         for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
//             fprintf(output_file, "beta\n");
//             fprintf(output_file, "%lf\n", samples->beta_vals[t_idx]);
//             fprintf(output_file, "w_k,L_HH_mean,L_HH_std\n");
//             for (int k = 0; k < samples->k_max; k++) {
//                 fprintf(output_file, "%lf,%lf,%lf\n", samples->w_k[t_idx][k], samples->L_HH_mean[t_idx][k], samples->L_HH_std[t_idx][k]);
//             }
//         }
// #endif // L_HH
// #ifdef L_SH
//         for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
//             fprintf(output_file, "beta\n");
//             fprintf(output_file, "%lf\n", samples->beta_vals[t_idx]);
//             fprintf(output_file, "w_k,L_SH_mean,L_SH_std\n");
//             for (int k = 0; k < samples->k_max; k++) {
//                 fprintf(output_file, "%lf,%lf,%lf\n", samples->w_k[t_idx][k], samples->L_SH_mean[t_idx][k], samples->L_SH_std[t_idx][k]);
//             }
//         }

//         for (int t_idx = 0; t_idx < samples->betas; t_idx++) {
//             fprintf(output_file, "beta\n");
//             fprintf(output_file, "%lf\n", samples->beta_vals[t_idx]);
//             fprintf(output_file, "w_k,L_HS_mean,L_HS_std\n");
//             for (int k = 0; k < samples->k_max; k++) {
//                 fprintf(output_file, "%lf,%lf,%lf\n", samples->w_k[t_idx][k], samples->L_HS_mean[t_idx][k], samples->L_HS_std[t_idx][k]);
//             }
//         }
// #endif // L_HH
//     } else {
//         printf("Error in opening the %s file. \n", file_name);
//         exit(1);
//     }

//     fclose(output_file);
//     return file_name;
// }

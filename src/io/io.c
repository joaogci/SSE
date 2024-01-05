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
    strcat(filename, "_scal");

    out = fopen(filename, "a");
    write_obs_scalar(out, &(obs[i]));
    fclose(out);
  }
}


#include "io.h"

void read_vtx_info(Vertices** vtx, int* n_diagrams)
{
  int status;
  char* sse_path;
  char command[BUFFER_SIZE];
  FILE* vtx_file;

  sse_path = getenv("SSE_DIR");
  if (sse_path == NULL) {
    printf("The SSE_DIR enviroment variable is not set. Please run the build.sh script in the SSE directory. \n");
    printf("source build.sh \n");
    exit(1);
  }

  strcpy(command, "");
  strcat(command, "python3 ");
  strcat(command, sse_path);
  strcat(command, "/src/hamiltonian/gen_vtx.py parameters");
  status = system(command);
  if (status != 0) {
    printf("Error executing %s.\n", command);
    printf("Error core: %d\n", status);
    exit(1);
  }

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

void write_observables(Obs_scalar* obs_scal, int n_scal, Obs_latt* obs_eq, int n_eq)
{
  int i;
  FILE* out;
  FILE* out_i;
  FILE* out_k;
  char filename[BUFFER_SIZE];

  for (i = 0; i < n_scal; i++) {
    strcpy(filename, obs_scal[i].filename);
    strcat(filename, "_scal");

    out = fopen(filename, "a");
    write_obs_scalar(out, &(obs_scal[i]));
    fclose(out);
  }

  for (i = 0; i < n_eq; i++) {
    strcpy(filename, obs_eq[i].filename);
    strcat(filename, "_eqR");
    out_i = fopen(filename, "a");

    strcpy(filename, obs_eq[i].filename);
    strcat(filename, "_eqK");
    out_k = fopen(filename, "a");

    write_obs_latt(out_i, out_k, &(obs_eq[i]));

    fclose(out_i);
    fclose(out_k);
  }
}


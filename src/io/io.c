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

void read_hyperbolic_lattice(double*** pos, int*** adj_mat, int** bulk, int** sublattice, int* N, int p, int q, int nl)
{
  int i, j;
  char filename[BUFFER_SIZE];
  char *sse_path;
  FILE* file;

  sse_path = getenv("SSE_DIR");
  if (sse_path == NULL) {
    printf("The SSE_DIR enviroment variable is not set. Please run the build.sh script in the SSE directory. \n");
    printf("source build.sh \n");
    exit(1);
  }

  sprintf(filename, "%s/src/hamiltonian/adjacency_matrices/hopping_%d_%d_%d.dat", sse_path, p, q, nl);
  file = fopen(filename, "r");

  if (file != NULL) {
    fscanf(file, "%d \n", N);
    (*adj_mat) = (int**) malloc((*N) * sizeof(int*));
    
    for (i = 0; i < (*N); i++) {
      (*adj_mat)[i] = (int*) malloc((*N) * sizeof(int));

      for (j = 0; j < (*N); j++) {
        fscanf(file, "%d.", &((*adj_mat)[i][j]));
      }
    }
  } else {
    printf("Error opening the adjacency matrix file. \n");
    printf("The values of p, q and nl do not have an adjacency matrix assigned. \n");
    exit(1);
  }

  fclose(file); 

  sprintf(filename, "%s/src/hamiltonian/adjacency_matrices/coord_%d_%d_%d.dat", sse_path, p, q, nl);
  file = fopen(filename, "r");

  if (file != NULL) {
    (*pos) = (double**) malloc((*N) * sizeof(double*));
    
    for (i = 0; i < (*N); i++) {
      (*pos)[i] = (double*) malloc(2 * sizeof(double));

      fscanf(file, "%lf %lf\n", &((*pos)[i][0]), &((*pos)[i][1]));
    }
  } else {
    printf("Error opening the coordinates file. \n");
    printf("The values of p, q and nl do not have coordinates assigned. \n");
    exit(1);
  }
  
  fclose(file);
  
  sprintf(filename, "%s/src/hamiltonian/adjacency_matrices/bulk_%d_%d_%d.dat", sse_path, p, q, nl);
  file = fopen(filename, "r");

  if (file != NULL) {
    (*bulk) = (int*) malloc((*N) * sizeof(int));
    
    for (i = 0; i < (*N); i++) {
      fscanf(file, " %d\n", &((*bulk)[i]));
    }
  } else {
    printf("Error opening the bulk file. \n");
    printf("The values of p, q and nl do not have bulk assigned. \n");
    exit(1);
  }
  
  fclose(file);

  sprintf(filename, "%s/src/hamiltonian/adjacency_matrices/sublattice_%d_%d_%d.dat", sse_path, p, q, nl);
  file = fopen(filename, "r");

  if (file != NULL) {
    (*sublattice) = (int*) malloc((*N) * sizeof(int));
    
    for (i = 0; i < (*N); i++) {
      fscanf(file, " %d\n", &((*sublattice)[i]));
    }
  } else {
    printf("Error opening the sublattice file. \n");
    printf("The values of p, q and nl do not have sublattice assigned. \n");
    exit(1);
  }
  
  fclose(file);
}

void write_observables(Obs_scalar* obs_scal, int n_scal, Obs_latt* obs_eq, int n_eq)
{
  int i;
  FILE* out,* out_i,* out_k,* info;
  char filename[BUFFER_SIZE];

  for (i = 0; i < n_scal; i++) {
    strcpy(filename, obs_scal[i].filename);
    strcat(filename, "_scal");

    out = fopen(filename, "a");
    write_obs_scalar(out, &(obs_scal[i]));
    fclose(out);

    strcat(filename, "_info");
    info = fopen(filename, "w");
    write_obs_scalar_info(info, &(obs_scal[i]));
    fclose(info);
  }

  for (i = 0; i < n_eq; i++) {
    strcpy(filename, obs_eq[i].filename);
    strcat(filename, "_eqR");
    out_i = fopen(filename, "a");

    // strcpy(filename, obs_eq[i].filename);
    // strcat(filename, "_eqK");
    // out_k = fopen(filename, "a");
    out_k = NULL;

    write_obs_latt(out_i, out_k, &(obs_eq[i]));

    fclose(out_i);
    // fclose(out_k);

    strcpy(filename, obs_eq[i].filename);
    strcat(filename, "_eqR_info");
    info = fopen(filename, "w");
    write_obs_eq_info(info, &(obs_eq[i]));
    fclose(info);

    // strcpy(filename, obs_eq[i].filename);
    // strcat(filename, "_eqK_info");
    // info = fopen(filename, "w");
    // write_obs_eq_info(info, &(obs_eq[i]));
    // fclose(info);
  }
}

void write_transport_obeservables(Obs_transport* obs, int n_transp)
{
  int i;
  FILE* out,* info;
  char filename[BUFFER_SIZE], tmp[BUFFER_SIZE];

  for (i = 0; i < n_transp; i++) {
    strcpy(filename, obs[i].filename);
    
    strcat(filename, "_");
    sprintf(tmp, "%d", obs[i].x);
    strcat(filename, tmp);
    
    strcat(filename, "_");
    sprintf(tmp, "%d", obs[i].y);
    strcat(filename, tmp);
    strcat(filename, "_transp");

    out = fopen(filename, "a");
    write_obs_transport(out, &(obs[i]));
    fclose(out);

    strcat(filename, "_info");
    info = fopen(filename, "w");
    write_obs_transport_info(info, &(obs[i]));
    fclose(info);
  }
}

void write_sim_info(Sim_info sim)
{
  FILE* out;

  out = fopen("info", "a");
  fprintf(out, "  -------------------------------  \n");
  fprintf(out, "Thermalization cycles:   %ld \n", sim.therm_cycles);
  fprintf(out, "Number of bins, sweeps:  %d, %ld\n", sim.n_bins, sim.mc_sweeps);
  fprintf(out, "Number of threads:       %d \n", sim.n_threads);
  fprintf(out, "Average number of loops: %lf \n", sim.avg_n_loops);
  fprintf(out, "Wall Time:               %lfs \n", sim.wall_time);
  fclose(out);
}

int num_lines(FILE* fp)
{
  int ch, lines;

  lines = 0;
  while(!feof(fp))
  {
    ch = fgetc(fp);
    if(ch == '\n')
    {
      lines++;
    }
  }

  return lines;
}

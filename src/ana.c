#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "sampling/observables.h"
#include "hamiltonian/lattice.h"
#include "io/io.h"
#include "analysis/analysis.h"

// Lattice
int Lx, Ly;
char boundary_condition;
double a_1[2], a_2[2];
Lattice latt;

// Transport
int x, y, n_max;
double beta;

// Observables
Ana_scalar obs_scal;
Ana_latt obs_eq;
Ana_transp obs_transp;

// Rebin
int n_rebin;


int main(int argc, char** argv)
{
  int i, n, n_bins;
  size_t len;
  char filename[BUFFER];
  char* line;
  FILE* out;
  FILE* info,* res;

  if (argc <= 1) {
    printf("Please provide the files to be analysed. \n");
    printf("Usage: %s file_1 file_2 ... \n ", argv[0]);
    exit(1);
  }

  read_parameters_analysis(&n_rebin);

  // Analysis
  for (i = 1; i < argc; i++) {
    if (strstr(argv[i], "_info") == NULL && (strstr(argv[i], "_scal") != NULL || strstr(argv[i], "_eq") != NULL || strstr(argv[i], "_transp") != NULL)) {
    // strcmp(argv[i], "vertices_info") != 0 && strcmp(argv[i], "parameters") != 0 && strcmp(argv[i], "info") != 0) {
      if (access(argv[i], F_OK) != 0) {
        printf("%s file does not exist. \n", argv[i]);
        continue;
      }
      
      if (strstr(argv[i], "J") != NULL || strstr(argv[i], "_eqJR") != NULL || strstr(argv[i], "_eqKJ") != NULL || strstr(argv[i], "_transpJ") != NULL) {
        continue;
      }

      strcpy(filename, argv[i]);
      strcat(filename, "_info");
      info = fopen(filename, "r");

      line = NULL;
      getline(&line, &len, info);
      getline(&line, &len, info);
      line[strlen(line) - 1] = '\0';

      // strcpy(filename, argv[i]);
      // filename[strlen(filename) - 5] = '\0';
      
      switch (line[0]) 
      {
        case 'i':
          obs_scal.obs_mean = 0.0;
          obs_scal.obs_std = 0.0;

          strcpy(filename, argv[i]);
          res = fopen(filename, "r");
          len = num_lines(res);
          n_bins = len;
          fclose(res);
          
          printf("Analysing %s \n", filename);
          printf("Number of bins: %d \n", n_bins);
          printf("Number of effective bins: %d \n", n_bins / n_rebin);

          res = fopen(filename, "r");
          analyse_scal(res, &obs_scal, n_bins, n_rebin);
          fclose(res);

          strcat(filename, "J");
          out = fopen(filename, "w");
          write_scal(out, &obs_scal);
          fclose(out);
          
          break;
        case 'e':
          getline(&line, &len, info);
          fscanf(info, "L1: %d\n", &Lx);
          fscanf(info, "L2: %d\n", &Ly);
          fscanf(info, "a1: (%lf, %lf)\n", &(a_1[0]), &(a_1[1]));
          fscanf(info, "a2: (%lf, %lf)\n", &(a_2[0]), &(a_2[1]));
          fscanf(info, "boundary condition: %c\n", &boundary_condition);

          make_lattice(Lx, Ly, a_1, a_2, boundary_condition, &latt);
          
          obs_eq.obs_mean = (double _Complex*) malloc(latt.N * sizeof(double _Complex));
          obs_eq.obs_std = (double _Complex*) malloc(latt.N * sizeof(double _Complex));
          for (n = 0; n < latt.N; n++) {
            obs_eq.obs_mean[n] = 0.0;
            obs_eq.obs_std[n] = 0.0;
          }

          strcpy(filename, argv[i]);
          // strcat(filename, "R");
          res = fopen(filename, "r");
          len = num_lines(res);
          n_bins = len / latt.N;
          fclose(res);

          printf("Analysing %s \n", filename);
          printf("Number of bins: %d \n", n_bins);
          printf("Number of effective bins: %d \n", n_bins / n_rebin);

          res = fopen(filename, "r");
          analyse_latt(res, &obs_eq, &latt, n_bins, n_rebin);
          fclose(res);

          filename[strlen(filename) - 1] = '\0';
          strcat(filename, "JR");
          out = fopen(filename, "w");
          write_latt_r(out, &obs_eq, &latt);
          fclose(out);

          for (n = 0; n < latt.N; n++) {
            obs_eq.obs_mean[n] = 0.0;
            obs_eq.obs_std[n] = 0.0;
          }

          filename[strlen(filename) - 2] = '\0';
          strcat(filename, "K");
          res = fopen(filename, "r");
          len = num_lines(res);
          n_bins = len / latt.N;
          fclose(res);

          printf("Analysing %s \n", filename);
          printf("Number of bins: %d \n", n_bins);
          printf("Number of effective bins: %d \n", n_bins / n_rebin);

          res = fopen(filename, "r");
          analyse_latt(res, &obs_eq, &latt, n_bins, n_rebin);
          fclose(res);


          filename[strlen(filename) - 1] = '\0';
          strcat(filename, "JK");
          out = fopen(filename, "w");
          write_latt_k(out, &obs_eq, &latt);
          fclose(out);

          free(obs_eq.obs_mean);
          free(obs_eq.obs_std);

          free_lattice(&latt);
          break;
        case 't':
          getline(&line, &len, info);
          fscanf(info, "x: %d\n", &x);
          fscanf(info, "y: %d\n", &y);
          fscanf(info, "n_max: %d\n", &n_max);
          fscanf(info, "beta: %lf\n", &beta);

          obs_transp.obs_mean = (double _Complex*) malloc(n_max * sizeof(double _Complex));
          obs_transp.obs_std = (double _Complex*) malloc(n_max * sizeof(double _Complex));
          for (n = 0; n < n_max; n++) {
            obs_transp.obs_mean[n] = 0.0;
            obs_transp.obs_std[n] = 0.0;
          }

          strcpy(filename, argv[i]);
          res = fopen(filename, "r");
          len = num_lines(res);
          n_bins = len / n_max;
          fclose(res);

          printf("Analysing %s \n", filename);
          printf("Number of bins: %d \n", n_bins);
          printf("Number of effective bins: %d \n", n_bins / n_rebin);

          res = fopen(filename, "r");
          analyse_transp(res, &obs_transp, n_max, n_bins, n_rebin);
          fclose(res);

          strcat(filename, "J");
          out = fopen(filename, "w");
          write_transp(out, &obs_transp, n_max, beta);
          fclose(out);

          free(obs_transp.obs_mean);
          free(obs_transp.obs_std);
          break;
      }
      
      fclose(info);
    }
  }

  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "sse/sse.h"
#include "sampling/observables.h"
#include "hamiltonian/hamiltonian.h"
#include "hamiltonian/lattice.h"
#include "io/io.h"
#include "analysis/analysis.h"

// Lattice
int Lx, Ly;
char boundary_condition;
double a_1[2], a_2[2];
Lattice latt;

// Simulation
Sim_info sim;
double beta;
double epsilon;

// Hamiltonian
double S, J_perp, J_par, h, D;
XXZ_ham ham;

// Observables
Ana_scalar obs_scal;
Ana_latt obs_eq;


int main(int argc, char** argv)
{
  int i, n;
  char filename[BUFFER];
  FILE* out;
  FILE* in;

  if (argc <= 1) {
    printf("Please provide the files to be analysed. \n");
    printf("Usage: %s file_1 file_2 ... \n ", argv[0]);
    exit(1);
  }

  read_parameters(&Lx, &Ly, &boundary_condition, &beta, &epsilon, &(sim.therm_cycles), &(sim.n_bins), &(sim.mc_sweeps), &S, &J_perp, &J_par,  &h, &D);

  // Lattice
  a_1[0] = 1.0; a_1[1] = 0.0;
  a_2[0] = 0.0; a_2[1] = 1.0;
  make_lattice(Lx, Ly, a_1, a_2, boundary_condition, &latt);

  // Hamiltonian
  init_ham(S, J_par, J_perp, h, D, epsilon, &ham);
  read_vtx_info(&(ham.vertices), &(ham.n_diagrams));
  ham.latt = &latt;
  ham.C = S*S * fabs(J_par) + 2.0 * (fabs(h) * S + fabs(D) * S*S) / latt.z + epsilon; 

  // Analysis
  for (i = 1; i < argc; i++) {
    if (strcmp(argv[i], "parameters") == 0 || strcmp(argv[i], "vertices_info") == 0) {
      continue;
    }

    if (access(argv[i], F_OK) != 0) {
      printf("%s file does not exist. \n", argv[i]);
      continue;
    }

    in = fopen(argv[i], "r");

    if (strstr(argv[i], "_scal") != NULL && strstr(argv[i], "_scalJ") == NULL) {
      obs_scal.obs_mean = 0.0;
      obs_scal.obs_std = 0.0;

      analyse_scal(in, &sim, &ham, &obs_scal);

      strcpy(filename, argv[i]);
      strcat(filename, "J");
      out = fopen(filename, "w");
      write_scal(out, &obs_scal);
      fclose(out);
    } else if (strstr(argv[i], "_eqR") != NULL) {
      obs_eq.latt = &latt;
      obs_eq.obs_mean = (double _Complex*) malloc(latt.N * sizeof(double _Complex));
      obs_eq.obs_std = (double _Complex*) malloc(latt.N * sizeof(double _Complex));
      for (n = 0; n < latt.N; n++) {
        obs_eq.obs_mean[n] = 0.0;
        obs_eq.obs_std[n] = 0.0;
      }

      analyse_latt(in, &sim, &ham, &obs_eq);

      strcpy(filename, argv[i]);
      filename[strlen(filename)-1] = '\0';
      strcat(filename, "JR");
      out = fopen(filename, "w");
      write_latt_r(out, &obs_eq);
      fclose(out);

      free(obs_eq.obs_mean);
      free(obs_eq.obs_std);
    } else if (strstr(argv[i], "_eqK") != NULL) {
      obs_eq.latt = &latt;
      obs_eq.obs_mean = (double _Complex*) malloc(latt.N * sizeof(double _Complex));
      obs_eq.obs_std = (double _Complex*) malloc(latt.N * sizeof(double _Complex));
      for (n = 0; n < latt.N; n++) {
        obs_eq.obs_mean[n] = 0.0;
        obs_eq.obs_std[n] = 0.0;
      }

      analyse_latt(in, &sim, &ham, &obs_eq);

      strcpy(filename, argv[i]);
      filename[strlen(filename)-1] = '\0';
      strcat(filename, "JK");
      out = fopen(filename, "w");
      write_latt_k(out, &obs_eq);
      fclose(out);

      free(obs_eq.obs_mean);
      free(obs_eq.obs_std);
    }

    fclose(in);
  }

  return 0;
}

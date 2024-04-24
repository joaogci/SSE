#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <unistd.h>
#include <stdbool.h>

#include "sse/sse.h"
#include "rng/pcg_basic.h"
#include "sampling/sampling.h"
#include "sampling/observables.h"
#include "sampling/transport.h"
#include "hamiltonian/hamiltonian.h"
#include "hamiltonian/lattice.h"
#include "io/io.h"

#ifndef TESTING
#define SEED (u_int64_t) time(NULL)
#else
#define SEED (u_int64_t) 2
#endif

// Lattice
int Lx, Ly;
char boundary_condition;
double a_1[2], a_2[2];
Lattice latt;

// Simulation
SSE_config state;
Sim_info sim;
double beta;
double epsilon;

// Hamiltonian
double S, J_perp, J_par, h, D;
XXZ_ham ham;

// Observables
Obs_scalar* obs_scal;
Obs_latt* obs_eq;
Obs_transport* obs_transp;
int n_scal = 6;
int n_eq = 1;
int n_transp = 3;


int main(int argc, char **argv)
{
  clock_t start_clock, end_clock;
  time_t time_;

  // Read and process inputs
  if (argc < 3) {
    printf("Please provide the number of threads and the path to the main directory. \n");
    printf("Usage: %s n_threads SSE_DIR \n", argv[0]);
    exit(1);
  }

  read_parameters(&Lx, &Ly, &boundary_condition, &beta, &epsilon, &(sim.therm_cycles), &(sim.n_bins), &(sim.mc_sweeps), &S, &J_perp, &J_par,  &h, &D);
  sim.n_threads = atoi(argv[1]);
  omp_set_num_threads(sim.n_threads);
  
  if (sim.n_bins < sim.n_threads) {
    printf("The number of bins (%d) should be equal or larger than the number of threads (%d).\n", sim.n_bins, sim.n_threads);
    exit(1);
  }

  // Lattice
  a_1[0] = 1.0; a_1[1] = 0.0;
  a_2[0] = 0.0; a_2[1] = 1.0;
  make_lattice(Lx, Ly, a_1, a_2, boundary_condition, &latt);

  // Hamiltonian
  init_ham(S, J_par, J_perp, h, D, epsilon, &ham);
  read_vtx_info(&(ham.vertices), &(ham.n_diagrams), argv[2]);
  ham.latt = &latt;
  ham.C = S*S * fabs(J_par) + 2.0 * (fabs(h) * S + fabs(D) * S*S) / latt.z + epsilon; 

  time_ = time(NULL);
  printf("Starting SSE simulation\n");
  printf("Simulation started at: %s \n", ctime(&time_));
  fflush(stdout);

  start_clock = clock();

  #pragma omp parallel private(state, obs_scal, obs_eq, obs_transp)
  {
    int thread_id, start_bin, end_bin;
    int n;    
    long t;
    char filename[BUFFER_SIZE];
    pcg32_random_t rng;

    thread_id = omp_get_thread_num();
    start_bin = (thread_id * sim.n_bins) / sim.n_threads;
    end_bin = ((thread_id + 1) * sim.n_bins) / sim.n_threads; 

    init_sse_config(beta, ham.latt->N, &state);
    reset_sse_config(ham.latt->N, ham.Sz[0], &state);
    pcg32_srandom_r(&rng, (SEED * (thread_id + 1)) ^ (intptr_t)&rng, (SEED * (thread_id + 1)));

    obs_scal = (Obs_scalar*) malloc(n_scal * sizeof(Obs_scalar));
    obs_eq = (Obs_latt*) malloc(n_eq * sizeof(Obs_latt));
    set_observables(obs_scal, n_scal, obs_eq, n_eq, &latt);
    obs_transp = (Obs_transport*) malloc(n_transp * sizeof(Obs_transport));
    set_observables_transport(obs_transp, n_transp, beta, &ham);

    // Thermalization
    sprintf(filename, "confin_%d", thread_id);
    if (access(filename, F_OK) == 0) {
      read_configuration(thread_id, latt.N, &(state));
    } else {
      for (t = 0; t < sim.therm_cycles; t++) {
        diag_update(&ham, &state, &rng);

        create_vtx_list(&ham, &state);
        loop_update(&ham, &state, &rng);

        ajust_cutoff(&state, t);
      }
    }
    

    // Measurement sweeps
    for (n = start_bin; n < end_bin; n++) {
      reset_observables(obs_scal, n_scal, obs_eq, n_eq);
      reset_observables_transport(obs_transp, n_transp);

      for (t = 0; t < sim.mc_sweeps; t++) {
        diag_update(&ham, &state, &rng);

        state.loop_size = 0;
        create_vtx_list(&ham, &state);
        loop_update(&ham, &state, &rng);

        if (n_scal > 0 || n_eq > 0) {
          sample(obs_scal, n_scal, obs_eq, n_eq, &ham, &state);
        }
        if (n_transp > 0) {
          sample_transport(obs_transp, n_transp, &ham, &state, &rng);
        }
      }

      // Write bin to file
      if (n_scal > 0 || n_eq > 0) {
        #pragma omp critical
        write_observables(obs_scal, n_scal, obs_eq, n_eq);
      }
      if (n_transp > 0) {
        #pragma omp critical
        write_transport_obeservables(obs_transp, n_transp);
      }
      #pragma omp critical
      write_configuration(thread_id, latt.N, &state);
    }

    #pragma omp critical
    write_histogram(thread_id, latt.N, &state);

    free_sse_config(&state);
    free(obs_eq);
    free(obs_scal);
    free(obs_transp);
  }

  end_clock = clock();
  sim.wall_time = ((double) (end_clock - start_clock)) / (CLOCKS_PER_SEC * sim.n_threads);

  write_sim_info(sim);

  printf("Simulation finished in %.5lfs \n", sim.wall_time);
  fflush(stdout);

  // FREE THE VARIABLES
  free_ham(&ham);
  free_lattice(&latt);

  exit(0);
}

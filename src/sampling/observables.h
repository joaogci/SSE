#ifndef OBSERVABLES_H 
#define OBSERVABLES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <stdint.h>

#include "../hamiltonian/lattice.h"

#define BUFFER 32

/*
 * Struct for scalar observables
 */
typedef struct Obs_scalar
{
  u_int64_t N;
  char filename[BUFFER];

  double _Complex obs_vec;
} Obs_scalar;

/*
 * Struct for equal time correlation functions
 * or static susceptibilities 
 */
typedef struct Obs_latt
{
  u_int64_t N;
  char filename[BUFFER];

  double _Complex** obs_latt;
  double _Complex* obs_latt0;

  double _Complex* obs_k;
  double _Complex* obs_i;

  Lattice* latt;
} Obs_latt;

/* 
 * Struct for transport coefficients
 */
typedef struct Obs_transport
{
  u_int64_t N;
  char filename[BUFFER];
  
  int x;
  int y;
  int n_max;
  double beta;
  double* omega_n;

  double _Complex* obs_transport;
} Obs_transport;


void init_obs_scalar(char* filename, Obs_scalar* obs);
void init_obs_latt(char* filename, Lattice* latt, Obs_latt* obs);
void init_obs_transport(char* filename, int x, int y, double beta, int n_max, Obs_transport* obs);

void reset_obs_scalar(Obs_scalar* obs);
void reset_obs_latt(Obs_latt* obs);
void reset_obs_transport(Obs_transport* obs);

void write_obs_scalar(FILE* out, Obs_scalar* obs);
void write_obs_latt(FILE* out_i, FILE* out_k, Obs_latt* obs);
void write_obs_transport(FILE* out, Obs_transport* obs);

void write_obs_scalar_info(FILE* info, Obs_scalar* obs);
void write_obs_eq_info(FILE* info, Obs_latt* obs);
void write_obs_transport_info(FILE* out, Obs_transport* obs);

void free_obs_latt(Obs_latt* obs);
void free_obs_transport(Obs_transport* obs);

void fourier_trans(Obs_latt* obs);
void inv_fourier_trans(Obs_latt* obs);

#endif // OBSRVABLES_H

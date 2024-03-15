#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "../sse/sse.h"
#include "../sampling/observables.h"
#include "../hamiltonian/hamiltonian.h"
#include "../hamiltonian/lattice.h"
#include "../hamiltonian/lattice_hyperbolic.h"

#define BUFFER 32

/*
 * Struct for scalar observables
 */
typedef struct Ana_scalar
{
  double _Complex obs_mean;
  double _Complex obs_std;
} Ana_scalar;

/*
 * Struct for lattice observables
 */
typedef struct Ana_latt
{
  double _Complex* obs_mean;
  double _Complex* obs_std;
} Ana_latt;

/*
 * Struct for transport observables
 */
typedef struct Ana_transp
{
  double _Complex* obs_mean;
  double _Complex* obs_std;
} Ana_transp;


/*
 * Analyses scalar observable
 */
void analyse_scal(FILE* fp, Ana_scalar* obs, int n_bins);

/*
 * Analyses lattice observable
 */
void analyse_latt(FILE* fp, Ana_latt* obs, Lattice_Hyperbolic* latt, int n_bins);

/*
 * Analyses transport observable
 */
void analyse_transp(FILE* fp, Ana_transp* obs, int n_max, int n_bins);

/*
 * Writes the results for scalar observable
 */
void write_scal(FILE* fp, Ana_scalar* obs);

/*
 * Writes the results for lattice observable
 */
void write_latt_r(FILE* fp, Ana_latt* obs, Lattice_Hyperbolic* latt);
void write_latt_k(FILE* fp, Ana_latt* obs, Lattice* latt);

/*
 * Wrties the results for transport observable
 */
void write_transp(FILE* fp, Ana_transp* obs, int n_max, double beta);

#endif // ANALYSIS_H

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

  Lattice* latt;
} Ana_latt;


/*
 * Analyses scalar observable
 */
void analyse_scal(FILE* fp, Sim_info* sim, XXZ_ham* ham, Ana_scalar* obs);

/*
 * Analyses lattice observable
 */
void analyse_latt(FILE* fp, Sim_info* sim, XXZ_ham* ham, Ana_latt* obs);

/*
 * Writes the results for scalar observable
 */
void write_scal(FILE* fp, Ana_scalar* obs);

/*
 * Writes the results for lattice observable
 */
void write_latt_r(FILE* fp, Ana_latt* obs);
void write_latt_k(FILE* fp, Ana_latt* obs);


#endif // ANALYSIS_H

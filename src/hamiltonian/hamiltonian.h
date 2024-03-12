#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <stdlib.h>

#include "lattice.h"
#include "lattice_hyperbolic.h"

#define N_LEGS      4
#define N_TYPES     3
#define N_UPDATES   2

// #define C_ ( S * S * fabs(delta) + 2.0 * S * hb_ )
// #define hb_ ( h / (2.0 * d) )

/* 
 * information about each vertex type
 */
typedef struct Vertices
{
  int indx;
  int type;
  int spin[N_LEGS];
  int new_vtx_type[N_UPDATES][N_LEGS][N_LEGS];
  double prob_exit[N_UPDATES][N_LEGS][N_LEGS];
  double weight;
} Vertices;

/* 
 *  information about the simulated system
 */
typedef struct XXZ_ham 
{
  double S;
  int n_proj;
  int* Sz;
  
  double J_par;
  double J_perp;
  double h;
  double D;
  
  double C;
  double epsilon;

  Lattice_Hyperbolic* latt;
  Vertices* vertices;
  int n_diagrams;
  int n_type;
} XXZ_ham;

/* 
 *  initializes the heisenberg_system struct
 */
void init_ham(double S, double J_par, double J_perp, double h, double D, double epsilon, XXZ_ham* ham);

/* 
 *  frees allocated memory in the xxz_ham
 */
void free_ham(XXZ_ham* ham);

#endif // HAMILTONIAN_H

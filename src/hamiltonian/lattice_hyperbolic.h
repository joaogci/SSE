#ifndef LATTICE_HYPERBOLIC_H
#define LATTICE_HYPERBOLIC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

/*
 * Hyperbolic lattice struct
 */
typedef struct Lattice_Hyperbolic
{
  int p;
  int q;
  int nl;

  double L;
  int N;
  int N_bulk;
  int Nb;

  double** r;

  int** bond_list;  

  int* bulk;
  int* sublattice;
} Lattice_Hyperbolic;

/*
 * makes hyperbolic lattice
 */
void make_lattice_hyperbolic(int p, int q, int nl, int*** adj_mat, Lattice_Hyperbolic* latt);

/*
 * frees hyperbolic lattice
 */
void free_lattice_hyperbolic(Lattice_Hyperbolic* latt);

#endif // LATTICE_HYPERBOLIC_H

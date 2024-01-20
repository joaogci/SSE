#ifndef LATTICE_H 
#define LATTICE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

/*
 * Square lattice struct
 */
typedef struct Lattice
{
  char bc;
  int L1;
  int L2;
  int N;
  int Nb;
  int d;
  int z;
  
  int** nn_list;
  int** bond_list;

  double a_1[2];
  double a_2[2];

  int** r;
  int** inv_r;
  int*** r_ij;

  double b_1[2];
  double b_2[2];

  int** k;
  int** inv_k;
  int*** k_ij;
} Lattice;

/* 
 * Initializes Lattice 
 */
void make_lattice(int L1, int L2, double* a_1, double* a_2, char bc, Lattice* lattice);

/* 
 * Frees the lattice 
 */
void free_lattice(Lattice* lattice);

#endif // LATTICE_H

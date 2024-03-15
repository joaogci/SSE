#include "lattice_hyperbolic.h"

void make_lattice_hyperbolic(int p, int q, int nl, int*** adj_mat, Lattice_Hyperbolic* latt)
{
  int i, j, b;

  latt->p = p;
  latt->q = q;
  latt->nl = nl;
  latt->L = sqrt(latt->N);
  
  latt->Nb = 0;
  for (i = 0; i < latt->N; i++) {
    for (j = 0; j <= i; j++) {
      if ((*adj_mat)[i][j] != 0) {
        latt->Nb++;
      }
    }
  }

  /* bond list */
  latt->bond_list = (int**) malloc(latt->Nb * sizeof(int*));
  for (i = 0; i < latt->Nb; i++) {
    latt->bond_list[i] = (int*) malloc(2 * sizeof(int));
  }

  b = 0;
  for (i = 0; i < latt->N; i++) {
    for (j = 0; j <= i; j++) {
      if ((*adj_mat)[i][j] != 0) {
        latt->bond_list[b][0] = i;
        latt->bond_list[b][1] = j;
        b++;
      }
    }
  }

  for (i = 0; i < latt->N; i++) {
    free((*adj_mat)[i]);
  }
  free((*adj_mat));
}

void free_lattice_hyperbolic(Lattice_Hyperbolic* latt)
{
  int i;

  for (i = 0; i < latt->Nb; i++) {
    free(latt->bond_list[i]);
  }
  free(latt->bond_list);
}



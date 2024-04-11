#include "lattice.h"

void make_lattice(int L1, int L2, double* a_1, double* a_2, char bc, Lattice* lattice) 
{
  int i, i1, j, j1, n, nn, b;
  int NN[4][2] = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
  double C;

  if (L2 > L1) {
    printf("L2 must be smaller than L1. \n");
    exit(1);
  }
  if (L2 <= 0 || L1 <= 0) {
    printf("For 1D lattices set L1 to L and L2 to 1. \n");
    exit(1);
  }
  if (bc != 'O' && bc != 'P') {
    printf("Boundary conditions must be 'O' for open or 'P' for periodic.\n");
    exit(1);
  }

  lattice->bc = bc;
  lattice->L1 = L1;
  lattice->L2 = L2;
  lattice->N = L1 * L2;
  lattice->d = 2;
  lattice->z = 4;
  if (L2 == 1) 
  { 
    lattice->d = 1;
    lattice->z = 2;
  }
  lattice->Nb = lattice->z * lattice->N / 2;

  if (bc == 'O' && lattice->d == 1) {
    lattice->Nb = lattice->Nb - 1;
  } else if (bc == 'O' && lattice->d == 2) {
    lattice->Nb = lattice->Nb - (L1 + L2);
  }

  for (i = 0; i < 2; i++) {
    lattice->a_1[i] = a_1[i];
    lattice->a_2[i] = a_2[i];
  }

  /* Real-space lattice */
  lattice->r = (int**) malloc(lattice->N * sizeof(int*));
  for (i = 0; i < lattice->N; i++) lattice->r[i] = (int*) malloc(2 * sizeof(int));
  lattice->inv_r = (int**) malloc(lattice->L1 * sizeof(int*));
  for (i = 0; i < lattice->L1; i++) lattice->inv_r[i] = (int*) malloc(lattice->L2 * sizeof(int));
  lattice->r_ij = (int***) malloc(lattice->N * sizeof(int**));
  for (i = 0; i < lattice->N; i++) { 
    lattice->r_ij[i] = (int**) malloc(lattice->N * sizeof(int*));
    for (j = 0; j < lattice->N; j++) {
      lattice->r_ij[i][j] = (int*) malloc(2 * sizeof(int));
    }
  }

  n = 0;
  for (j = 0; j < L2; j++) {
    for (i = 0; i < L1; i++) {
      lattice->r[n][0] = (i + 1) - lattice->L1/2;
      lattice->r[n][1] = (j + 1) - lattice->L2/2;
      lattice->inv_r[i][j] = n;
      n++;
    }
  }

  for (i = 0; i < lattice->N; i++) {
    for (j = 0; j < lattice->N; j++) {
      lattice->r_ij[i][j][0] = lattice->r[i][0] - lattice->r[j][0];
      lattice->r_ij[i][j][1] = lattice->r[i][1] - lattice->r[j][1];
    }
  }

  /* Reciprocal lattice from periodicity in L */
  C = 2 * M_PI / (a_1[0] * a_2[1] - a_1[1] * a_2[0]);
  lattice->b_1[0] = C * (a_2[1]) / ((double) L1);
  lattice->b_1[1] = C * (- a_2[0]) / ((double) L1);
  if (L2 != 1) {
    lattice->b_2[0] = C * (- a_1[1]) / ((double) L2);
    lattice->b_2[1] = C * (a_1[0]) / ((double) L2);
  } else {
    lattice->b_2[0] = 0.0;
    lattice->b_2[1] = 0.0;
  }

  lattice->k = (int**) malloc(lattice->N * sizeof(int*));
  for (i = 0; i < lattice->N; i++) lattice->k[i] = (int*) malloc(2 * sizeof(int));
  lattice->inv_k = (int**) malloc(lattice->L1 * sizeof(int*));
  for (i = 0; i < lattice->L1; i++) lattice->inv_k[i] = (int*) malloc(lattice->L2 * sizeof(int));
  lattice->k_ij = (int***) malloc(lattice->N * sizeof(int**));
  for (i = 0; i < lattice->N; i++) { 
    lattice->k_ij[i] = (int**) malloc(lattice->N * sizeof(int*));
    for (j = 0; j < lattice->N; j++) {
      lattice->k_ij[i][j] = (int*) malloc(2 * sizeof(int));
    }
  }

  n = 0;
  for (j = 0; j < L2; j++) {
    for (i = 0; i < L1; i++) {
      lattice->k[n][0] = (i + 1) - lattice->L1/2;
      lattice->k[n][1] = (j + 1) - lattice->L2/2;
      lattice->inv_k[i][j] = n;
      n++;
    }
  }

  for (i = 0; i < lattice->N; i++) {
    for (j = 0; j < lattice->N; j++) {
      lattice->k_ij[i][j][0] = lattice->k[i][0] - lattice->k[j][0];
      lattice->k_ij[i][j][1] = lattice->k[i][1] - lattice->k[j][1];
    }
  }

  /* NN list */
  lattice->nn_list = (int**) malloc(lattice->N * sizeof(int*));
  for (i = 0; i < lattice->N; i++) {
    lattice->nn_list[i] = (int*) malloc(lattice->z * sizeof(int));
  }

  switch (bc)
  {
  case 'P':
    for (i = 0; i < lattice->L1; i++) {
      for (j = 0; j < lattice->L2; j++) {
        for (nn = 0; nn < lattice->z; nn++) {
          i1 = (i + NN[nn][0])%lattice->L1;
          if (i1 < 0) i1 = lattice->L1 - 1;

          j1 = (j + NN[nn][1])%lattice->L2;
          if (j1 < 0) j1 = lattice->L2 - 1;

          lattice->nn_list[lattice->inv_r[i][j]][nn] = lattice->inv_r[i1][j1];
        }
      }
    }        
    break;
  case 'O':
    for (i = 0; i < lattice->L1; i++) {
      for (j = 0; j < lattice->L2; j++) {
        for (nn = 0; nn < lattice->z; nn++) {
          i1 = i + NN[nn][0];
          j1 = j + NN[nn][1];

          if (i1 < lattice->L1 && i1 >= 0 && j1 < lattice->L2 && j1 >= 0)
            lattice->nn_list[lattice->inv_r[i][j]][nn] = lattice->inv_r[i1][j1];
          else 
            lattice->nn_list[lattice->inv_r[i][j]][nn] = -1;
        }
      }
    }        
    break;
  }

  /* bonds */
  lattice->bond_list = (int**) malloc(lattice->Nb * sizeof(int*));
  for (i = 0; i < lattice->Nb; i++) {
    lattice->bond_list[i] = (int*) malloc(2 * sizeof(int));
  }

  b = 0;
  for (n = 0; n < lattice->N; n++) {
    for (i = 0; i < lattice->z / 2; i++) {
      if (lattice->nn_list[n][2*i+1] != -1) {
        lattice->bond_list[b][0] = n;
        lattice->bond_list[b][1] = lattice->nn_list[n][2*i+1];
        b++;
      }
    }
  }
}

void free_lattice(Lattice* lattice) 
{
  int i, j;

  for (i = 0; i < lattice->N; i++) {
    for (j = 0; j < lattice->N; j++) {
      free(lattice->r_ij[i][j]);
    }
    free(lattice->r_ij[i]);
    free(lattice->r[i]);
    free(lattice->k[i]);
    free(lattice->nn_list[i]);
  }   
  free(lattice->r);
  free(lattice->r_ij);
  free(lattice->k);
  free(lattice->nn_list);

  for (i = 0; i < lattice->L1; i++) {
    free(lattice->inv_r[i]);
    free(lattice->inv_k[i]);
  }
  free(lattice->inv_r);
  free(lattice->inv_k);

  for (i = 0; i < lattice->Nb; i++) {
    free(lattice->bond_list[i]);
  }
  free(lattice->bond_list);
}

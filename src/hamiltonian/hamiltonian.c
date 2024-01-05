#include "hamiltonian.h"

void init_ham(double S, double J_par, double J_perp, double h, double D, double epsilon, XXZ_ham* ham)
{
  ham->S = S;
  ham->n_proj = 2.0 * S + 1.0;
  ham->Sz = (int *) malloc(ham->n_proj * sizeof(int));
  for (int i = 0; i < ham->n_proj; i++) {
    ham->Sz[i] = - 2 * (S - i);
  }

  ham->J_par = J_par;
  ham->J_perp = J_perp;
  ham->h = h;
  ham->D = D;
  ham->epsilon = epsilon;
}

void free_ham(XXZ_ham* ham)
{
  free(ham->Sz);
  free(ham->vertices);
}


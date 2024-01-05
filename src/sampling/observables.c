#include "observables.h"

void init_obs_scalar(char* filename, Obs_scalar* obs)
{
  strcpy(obs->filename, filename);
}

void init_obs_latt(char* filename, Lattice* latt, Obs_latt* obs)
{
  strcpy(obs->filename, filename);
  obs->latt = latt;

  obs->obs_latt = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
  obs->obs_latt0 = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
}

void init_obs_transport(char* filename, int x, int y, double beta, int n_max, Obs_transport* obs)
{
  int n;

  strcpy(obs->filename, filename);
  obs->beta = beta;
  obs->n_max = n_max;
  obs->x = x;
  obs->y = y;

  obs->omega_n = (double*) malloc(n_max * sizeof(double));
  for (n = 1; n <= n_max; n++) {
    obs->omega_n[n] = 2.0 * M_PI * n / beta;
  }

  obs->obs_transport = (double _Complex*) malloc(n_max * sizeof(double _Complex));
}

void reset_obs_scalar(Obs_scalar* obs)
{
  obs->N = 0;
  obs->obs_vec = 0.0;
}

void reset_obs_latt(Obs_latt* obs)
{
  int i;

  obs->N = 0;
  for (i = 0; i < obs->latt->N; i++) {
    obs->obs_latt[i] = 0.0;
    obs->obs_latt0[i] = 0.0;
  }
}

void reset_obs_transport(Obs_transport* obs)
{
  int n;

  obs->N = 0;
  for (n = 1; n <= obs->n_max; n++) 
    obs->obs_transport[n] = 0.0;
}

void write_obs_scalar(FILE* out, Obs_scalar* obs)
{ 
  obs->obs_vec = obs->obs_vec / obs->N;
  fprintf(out, "(%lf, %lf) \n", creal(obs->obs_vec), cimag(obs->obs_vec));
}

void write_obs_latt(FILE* out, Obs_latt* obs)
{
  int i;
  
  for (i = 0; i < obs->latt->N; i++) {
    obs->obs_latt[i] = obs->obs_latt[i] / obs->N;
    obs->obs_latt0[i] = obs->obs_latt0[i] / obs->N;
  }
}

void write_obs_transport(FILE* out, Obs_transport* obs)
{

}

void free_obs_latt(Obs_latt* obs)
{
  free(obs->obs_latt);
  free(obs->obs_latt0);
}

void free_obs_transport(Obs_transport* obs)
{
  free(obs->omega_n);
  free(obs->obs_transport);
}

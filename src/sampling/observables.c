#include "observables.h"

void init_obs_scalar(char* filename, Obs_scalar* obs)
{
  strcpy(obs->filename, filename);
}

void init_obs_latt(char* filename, Lattice* latt, Obs_latt* obs)
{
  int i;

  strcpy(obs->filename, filename);
  obs->latt = latt;

  obs->obs_latt = (double _Complex**) malloc(latt->N * sizeof(double _Complex*));
  for (i = 0; i < latt->N; i++) {
    obs->obs_latt[i] = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
  }
  obs->obs_latt0 = (double _Complex*) malloc(latt->N * sizeof(double _Complex));

  obs->obs_i = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
  obs->obs_k = (double _Complex*) malloc(latt->N * sizeof(double _Complex));
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
    obs->omega_n[n - 1] = 2.0 * M_PI * n / beta;
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
  int i, j;

  obs->N = 0;
  for (i = 0; i < obs->latt->N; i++) {
    for (j = 0; j < obs->latt->N; j++) {
      obs->obs_latt[i][j] = 0.0;
    }
    obs->obs_latt0[i] = 0.0;
    obs->obs_i[i] = 0.0;
    obs->obs_k[i] = 0.0;
  }
}

void reset_obs_transport(Obs_transport* obs)
{
  int n;

  obs->N = 0;
  for (n = 0; n < obs->n_max; n++) 
    obs->obs_transport[n] = 0.0;
}

void write_obs_scalar(FILE* out, Obs_scalar* obs)
{ 
  obs->obs_vec = obs->obs_vec / obs->N;
  fprintf(out, "(%lf, %lf) \n", creal(obs->obs_vec), cimag(obs->obs_vec));
}

void write_obs_scalar_info(FILE* info, Obs_scalar* obs)
{
  fprintf(info, "-- Analysis Mode --\n");
  fprintf(info, "identity\n");
}

void write_obs_latt(FILE* out_i, FILE* out_k, Obs_latt* obs)
{
  int i, j;

  for (i = 0; i < obs->latt->N; i++) {
    for (j = 0; j < obs->latt->N; j++) {
      obs->obs_latt[i][j] = obs->obs_latt[i][j] / obs->N;
    }
    obs->obs_latt0[i] = obs->obs_latt0[i] / obs->N;
  }

  fourier_trans(obs);
  inv_fourier_trans(obs);

  for (i = 0; i < obs->latt->N; i++) {
    fprintf(out_i, "(%lf, %lf) (%lf, %lf) \n", (double) obs->latt->r[i][0], 
                                               (double) obs->latt->r[i][1], 
                                               creal(obs->obs_i[i]), 
                                               cimag(obs->obs_i[i]));
  }

  for (i = 0; i < obs->latt->N; i++) {
    fprintf(out_k, "(%lf, %lf) (%lf, %lf) \n", (double) (obs->latt->k[i][0] * obs->latt->b_1[0] + obs->latt->k[i][1] * obs->latt->b_2[0]), 
                                               (double) (obs->latt->k[i][0] * obs->latt->b_1[1] + obs->latt->k[i][1] * obs->latt->b_2[1]), 
                                               creal(obs->obs_k[i]), 
                                               cimag(obs->obs_k[i]));
  }
}

void write_obs_eq_info(FILE* info, Obs_latt* obs)
{
  fprintf(info, "-- Analysis Mode --\n");
  fprintf(info, "equal time\n");
  fprintf(info, "-- Lattice --\n");
  fprintf(info, "L1: %d\n", obs->latt->L1);
  fprintf(info, "L2: %d\n", obs->latt->L2);
  fprintf(info, "a1: (%lf, %lf)\n", obs->latt->a_1[0], obs->latt->a_1[1]);
  fprintf(info, "a2: (%lf, %lf)\n", obs->latt->a_2[0], obs->latt->a_2[1]);
  fprintf(info, "boundary condition: %c\n", obs->latt->bc);
}

void write_obs_transport(FILE* out, Obs_transport* obs)
{
  int k; 

  for (k = 0; k < obs->n_max; k++) {
    obs->obs_transport[k] = obs->obs_transport[k] / obs->N;
  }

  for (k = 0; k < obs->n_max; k++) {
    fprintf(out, "%lf (%lf, %lf) \n", obs->omega_n[k], creal(obs->obs_transport[k]), cimag(obs->obs_transport[k]));
  }
}

void write_obs_transport_info(FILE* info, Obs_transport* obs)
{
  fprintf(info, "-- Analysis Mode --\n");
  fprintf(info, "transport\n");
  fprintf(info, "-- Transport --\n");
  fprintf(info, "x: %d\n", obs->x);
  fprintf(info, "y: %d\n", obs->y);
  fprintf(info, "n_max: %d\n", obs->n_max);
  fprintf(info, "beta: %lf\n", obs->beta);
}

void free_obs_latt(Obs_latt* obs)
{
  int i;
  
  for (i = 0; i < obs->latt->N; i++) {
    free(obs->obs_latt[i]);
  }
  free(obs->obs_latt);
  free(obs->obs_latt0);
  free(obs->obs_i);
  free(obs->obs_k);
}

void free_obs_transport(Obs_transport* obs)
{
  free(obs->omega_n);
  free(obs->obs_transport);
}

void fourier_trans(Obs_latt* obs) 
{
  int i, j, n;
  double a[2], b[2];

  for (n = 0; n < obs->latt->N; n++) {
    for (i = 0; i < obs->latt->N; i++) {
      for (j = 0; j < obs->latt->N; j++) {
        a[0] = obs->latt->k[n][0] * obs->latt->b_1[0] + obs->latt->k[n][1] * obs->latt->b_2[0];
        a[1] = obs->latt->k[n][0] * obs->latt->b_1[1] + obs->latt->k[n][1] * obs->latt->b_2[1];
        b[0] = obs->latt->r_ij[i][j][0];
        b[1] = obs->latt->r_ij[i][j][1];

        obs->obs_k[n] += cexp(- I * (a[0] * b[0] + a[1] * b[1])) * (obs->obs_latt[i][j] - obs->obs_latt0[i] * obs->obs_latt0[j]) / (obs->latt->N);
      }
    }
  }
}

void inv_fourier_trans(Obs_latt* obs)
{
  int i, n;
  double a[2], b[2];

  for (i = 0; i < obs->latt->N; i++) {
    for (n = 0; n < obs->latt->N; n++) {
      a[0] = obs->latt->k[n][0] * obs->latt->b_1[0] + obs->latt->k[n][1] * obs->latt->b_2[0];
      a[1] = obs->latt->k[n][0] * obs->latt->b_1[1] + obs->latt->k[n][1] * obs->latt->b_2[1];
      b[0] = obs->latt->r[i][0];
      b[1] = obs->latt->r[i][1];

      obs->obs_i[i] += cexp(I * (a[0] * b[0] + a[1] * b[1])) * obs->obs_k[n] / (obs->latt->N);
    }
  }
}


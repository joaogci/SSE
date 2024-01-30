#include "analysis.h"

void analyse_scal(FILE* fp, Ana_scalar* obs, int n_bins)
{
  int i;
  double real, imag;
  double _Complex* measurements;

  measurements = (double _Complex*) malloc(n_bins * sizeof(double _Complex));
  for (i = 0; i < n_bins; i++) {
    fscanf(fp, "(%lf, %lf) ", &real, &imag);
    measurements[i] = real + I * imag;
  }

  for (i = 0; i < n_bins; i++) {
    obs->obs_mean += measurements[i];
  }  
  obs->obs_mean = obs->obs_mean / n_bins;

  for (i = 0; i < n_bins; i++) {
    obs->obs_std += cpow(measurements[i] - obs->obs_mean, 2.0);
  }
  obs->obs_std = csqrt(obs->obs_std / n_bins);

  free(measurements);
}

void analyse_latt(FILE* fp, Ana_latt* obs, Lattice *latt, int n_bins)
{
  int i, n;
  double real, imag, x, y;
  double _Complex** measurements;

  measurements = (double _Complex**) malloc(n_bins * sizeof(double _Complex*));
  for (i = 0; i < n_bins; i++) {
    measurements[i] = (double _Complex*) malloc(latt->N * sizeof(double _Complex));

    for (n = 0; n < latt->N; n++) {
      fscanf(fp, "(%lf, %lf) (%lf, %lf) \n", &x, &y, &real, &imag);
      measurements[i][n] = real + I * imag;
    }
  }

  for (n = 0; n < latt->N; n++) {
    for (i = 0; i < n_bins; i++) {
      obs->obs_mean[n] += measurements[i][n];
    }
    obs->obs_mean[n] = obs->obs_mean[n] / n_bins;
  }  

  for (n = 0; n < latt->N; n++) {
    for (i = 0; i < n_bins; i++) {
      obs->obs_std[n] += cpow(measurements[i][n] - obs->obs_mean[n], 2.0);
    }
    obs->obs_std[n] = csqrt(obs->obs_std[n] / n_bins);
  }

  for (i = 0; i < n_bins; i++) {
    free(measurements[i]);
  }
  free(measurements);
}

void analyse_transp(FILE* fp, Ana_transp* obs, int n_max, int n_bins)
{
  int i, n;
  double real, imag, omega_n;
  double _Complex** measurements;

  measurements = (double _Complex**) malloc(n_bins * sizeof(double _Complex*));
  for (i = 0; i < n_bins; i++) {
    measurements[i] = (double _Complex*) malloc(n_max * sizeof(double _Complex));

    for (n = 0; n < n_max; n++) {
      fscanf(fp, "%lf (%lf, %lf) \n", &omega_n, &real, &imag);
      measurements[i][n] = real + I * imag;
    }
  }

  for (n = 0; n < n_max; n++) {
    for (i = 0; i < n_bins; i++) {
      obs->obs_mean[n] += measurements[i][n];
    }
    obs->obs_mean[n] = obs->obs_mean[n] / n_bins;
  }  

  for (n = 0; n < n_max; n++) {
    for (i = 0; i < n_bins; i++) {
      obs->obs_std[n] += cpow(measurements[i][n] - obs->obs_mean[n], 2.0);
    }
    obs->obs_std[n] = csqrt(obs->obs_std[n] / n_bins);
  }

  for (i = 0; i < n_bins; i++) {
    free(measurements[i]);
  }
  free(measurements);
}

void write_scal(FILE* fp, Ana_scalar* obs)
{
  fprintf(fp, "%lf %lf %lf %lf \n", creal(obs->obs_mean), creal(obs->obs_std), cimag(obs->obs_mean), cimag(obs->obs_std));
}

void write_latt_r(FILE* fp, Ana_latt* obs, Lattice* latt)
{
  int i;

  for (i = 0; i < latt->N; i++) {
    fprintf(fp, "%lf %lf \n", (double) latt->r[i][0], (double) latt->r[i][1]);
    fprintf(fp, "%lf, %lf %lf, %lf \n", creal(obs->obs_mean[i]), creal(obs->obs_std[i]), cimag(obs->obs_mean[i]), cimag(obs->obs_std[i]));
  }
}

void write_latt_k(FILE* fp, Ana_latt* obs, Lattice* latt)
{
  int i;

  for (i = 0; i < latt->N; i++) {
    fprintf(fp, "%lf %lf \n", (double) (latt->k[i][0] * latt->b_1[0] + latt->k[i][1] * latt->b_2[0]), (double) (latt->k[i][0] * latt->b_1[1] + latt->k[i][1] * latt->b_2[1]));
    fprintf(fp, "%lf, %lf %lf, %lf \n", creal(obs->obs_mean[i]), creal(obs->obs_std[i]), cimag(obs->obs_mean[i]), cimag(obs->obs_std[i]));
  }
}

void write_transp(FILE* fp, Ana_transp* obs, int n_max, double beta)
{
  int i;

  for (i = 0; i < n_max; i++) {
    fprintf(fp, "%lf \n", 2.0 * M_PI * (i+1) / beta);
    fprintf(fp, "%lf, %lf %lf, %lf \n", creal(obs->obs_mean[i]), creal(obs->obs_std[i]), cimag(obs->obs_mean[i]), cimag(obs->obs_std[i]));
  }
}

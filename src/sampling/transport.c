#include "transport.h"

void set_observables_transport(Obs_transport* obs_transp, int n_transp, double beta, XXZ_ham* ham)
{
  int n;

  for (n = 0; n < n_transp; n++) {
    switch (n) 
    {
    case 0:
      init_obs_transport("L_SS", ham->latt->L1/2 - 1    , ham->latt->L1/2 - 1, beta, 30, &(obs_transp[n]));
      break;
    case 1:
      init_obs_transport("L_SS", ham->latt->L1/2 - 1 + 1, ham->latt->L1/2 - 1, beta, 30, &(obs_transp[n]));
      break;
    case 2:
      init_obs_transport("L_SS", ham->latt->L1/2 - 1 + 2, ham->latt->L1/2 - 1, beta, 30, &(obs_transp[n]));
      break;
    default:
      printf("Observable not found. \n");
      exit(1);
      break;
    }
  }
}

void sample_transport(Obs_transport* obs, int n_transp, XXZ_ham* ham, SSE_config* state, pcg32_random_t* rng)
{
  double _Complex** Ka,** Kb,** Ga,** Gb;
  int* Cab,* x_vals,* y_vals;
  double* omega_n;
  int n_obs, n_max, a, b, n, p, k;
  int Sa, Sb, bond, type;

  // Assign random times to the operators in the operator string
  assign_times(state, rng);

  n_obs = 3;
  n_max = obs[0].n_max;
  omega_n = (double*) malloc(n_max * sizeof(double));
  memcpy(omega_n, obs[0].omega_n, n_max * sizeof(double));

  x_vals = (int*) malloc(n_obs * sizeof(int));
  y_vals = (int*) malloc(n_obs * sizeof(int));
  for (n = 0; n < n_obs; n++) {
    x_vals[n] = obs[n].x;
    y_vals[n] = obs[n].y;
  }

  Ka = (double _Complex**) malloc(n_obs * sizeof(double _Complex*));
  Kb = (double _Complex**) malloc(n_obs * sizeof(double _Complex*));
  Ga = (double _Complex**) malloc(n_obs * sizeof(double _Complex*));
  Gb = (double _Complex**) malloc(n_obs * sizeof(double _Complex*));
  Cab = (int*) malloc(n_obs * sizeof(int));
  for (n = 0; n < n_obs; n++) {
    Ka[n] = (double _Complex*) malloc(n_max * sizeof(double _Complex));
    Kb[n] = (double _Complex*) malloc(n_max * sizeof(double _Complex));
    Ga[n] = (double _Complex*) malloc(n_max * sizeof(double _Complex));
    Gb[n] = (double _Complex*) malloc(n_max * sizeof(double _Complex));
  }

  for (n = 0; n < n_obs; n++) {
    Sa = 0.0;
    Sb = 0.0;

    for (a = x_vals[n] + 1; a < ham->latt->L1; a++) {
      Sa += state->spin_config[a];
    }
    for (b = y_vals[n] + 1; b < ham->latt->L1; b++) {
      Sb += state->spin_config[b];
    }

    Cab[n] = 0.0;
    for (k = 0; k < n_max; k++) {
      Ka[n][k] = 0.0;
      Kb[n][k] = 0.0;
      Ga[n][k] = 0.0;
      Gb[n][k] = 0.0;

      for (p = 0; p <= state->n; p++) {
        Ka[n][k] += Sa * 0.5 * (cexp(I * omega_n[k] * state->op_tau[p+1]) - cexp(I * omega_n[k] * state->op_tau[p])) / (I * omega_n[k]);
        Kb[n][k] += Sb * 0.5 * (cexp(I * omega_n[k] * state->op_tau[p+1]) - cexp(I * omega_n[k] * state->op_tau[p])) / (I * omega_n[k]);

        bond = (state->reduced_op_string[p] / N_TYPES) - 1;
        type = state->reduced_op_string[p] % N_TYPES;

        if (bond == x_vals[n] && type != 0) {
          if (type == 1) {
            Sa = Sa - 2;
          } else if (type == 2) {
            Sa = Sa + 2;
          }
        } 
        if (bond == y_vals[n] && type != 0) {
          if (type == 1) {
            Sb = Sb - 2;
          } else if (type == 2) {
            Sb = Sb + 2;
          }
        }

        if (p < state->n) {
          Ga[n][k] = (bond > x_vals[n]) ? cexp(I * omega_n[k] * state->op_tau[p+1]) : 0.0;
          Gb[n][k] = (bond > y_vals[n]) ? cexp(I * omega_n[k] * state->op_tau[p+1]) : 0.0;
          
          Cab[n] += (bond > y_vals[n] && bond > x_vals[n]) ? 1 : 0;
        }
      }
    }
  }

  for (n = 0; n < n_transp; n++) {
    obs[n].N++;
  }

  for (k = 0; k < n_max; k++) {
    obs[0].obs_transport[k] += obs[0].omega_n[k] * (Ka[0][k] * conj(Kb[0][k])) / obs[0].beta;
    obs[1].obs_transport[k] += obs[1].omega_n[k] * (Ka[1][k] * conj(Kb[1][k])) / obs[1].beta;
    obs[2].obs_transport[k] += obs[2].omega_n[k] * (Ka[2][k] * conj(Kb[2][k])) / obs[2].beta;
  }

  for (n = 0; n < n_obs; n++) {
    free(Ka[n]);
    free(Kb[n]);
    free(Ga[n]);
    free(Gb[n]);
  }
  free(Ka);
  free(Kb);
  free(Ga);
  free(Gb);
  free(Cab);
  free(omega_n);
  free(x_vals);
  free(y_vals);
}

void reset_observables_transport(Obs_transport* obs, int n_transp)
{
  int n;

  for (n = 0; n < n_transp; n++) {
    reset_obs_transport(&(obs[n]));
  }
}

void free_observables_transport(Obs_transport* obs, int n_transp)
{
  int n;

  for (n = 0; n < n_transp; n++) {
    free_obs_transport(&(obs[n]));
  }
}

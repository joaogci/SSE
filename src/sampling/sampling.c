#include "sampling.h"

void set_observables(Obs_scalar* obs_scalar, int n_scal, Obs_latt* obs_eq, int n_eq, Lattice_Hyperbolic* latt)
{
  int n;

  for (n = 0; n < n_scal; n++) {
    switch (n)
    {
    case 0:
      init_obs_scalar("E", &(obs_scalar[n]));
      break;
    case 1:
      init_obs_scalar("Sz", &(obs_scalar[n]));
      break;
    case 2:
      init_obs_scalar("Sz2", &(obs_scalar[n]));
      break;
    case 3:
      init_obs_scalar("n", &(obs_scalar[n]));
      break;
    case 4:
      init_obs_scalar("Sz_stag", &(obs_scalar[n]));
      break;
    case 5:
      init_obs_scalar("Sz2_stag", &(obs_scalar[n]));
      break;
    case 6:
      init_obs_scalar("Sz_suscep", &(obs_scalar[n]));
      break;
    case 7:
      init_obs_scalar("Sz_stag_bulk", &(obs_scalar[n]));
      break;
    case 8:
      init_obs_scalar("Sz2_stag_bulk", &(obs_scalar[n]));
      break;
    case 9:
      init_obs_scalar("Sz_stag_edge", &(obs_scalar[n]));
      break;
    case 10:
      init_obs_scalar("Sz2_stag_edge", &(obs_scalar[n]));
      break;
    // case 3:
    //   init_obs_scalar("n2", &(obs_scalar[n]));
    //   break;
    // case 4:
    //   init_obs_scalar("C", &(obs_scalar[n]));
    //   break;

    default:
      printf("Observable not found. \n");
      exit(1);
      break;
    }
  }

  for (n = 0; n < n_eq; n++) {
    switch (n)
    {
    case 0:
      init_obs_latt("Sz", latt, &(obs_eq[n]));
      break;
    default:
      printf("Observable not found. \n");
      exit(1);
      break;
    }
  }
}

void reset_observables(Obs_scalar* obs_scalar, int n_scal, Obs_latt* obs_eq, int n_eq)
{
  int n;

  for (n = 0; n < n_scal; n++) {
    reset_obs_scalar(&(obs_scalar[n]));
  }

  for (n = 0; n < n_eq; n++) {
    reset_obs_latt(&(obs_eq[n]));
  }
}

void sample(Obs_scalar* obs_scal, int n_scal, Obs_latt* obs_eq, int n_eq, XXZ_ham* ham, SSE_config* state)
{
  if (n_scal > 0) {
    sample_obs_scalar(obs_scal, n_scal, ham, state);
  }

  if (n_eq > 0) {
    sample_obs_eq(obs_eq, n_eq, ham, state);
  }
}

void sample_obs_scalar(Obs_scalar* obs, int n_scal, XXZ_ham* ham, SSE_config* state)
{
  int i;
  double Sz, Sz_stag, Sz_stag_bulk, Sz_stag_edge;

  for (i = 0; i < n_scal; i++) {
    obs[i].N++;
  }

  obs[0].obs_vec += - state->n / (state->beta * ham->latt->N) + ham->latt->Nb * ham->C / ham->latt->N;
  
  Sz = 0.0;
  Sz_stag = 0.0;
  Sz_stag_bulk = 0.0;
  Sz_stag_edge = 0.0;
  for (i = 0; i < ham->latt->N; i++) {
    Sz += state->spin_config[i] * 0.5;
    Sz_stag += ham->latt->sublattice[i] * state->spin_config[i] * 0.5;
    Sz_stag_bulk += ham->latt->bulk[i] * ham->latt->sublattice[i] * state->spin_config[i] * 0.5;
    Sz_stag_bulk += (1 - ham->latt->bulk[i]) * ham->latt->sublattice[i] * state->spin_config[i] * 0.5;
  }

  obs[1].obs_vec += Sz / ham->latt->N;
  obs[2].obs_vec += (Sz / ham->latt->N) * (Sz / ham->latt->N);

  obs[3].obs_vec += state->n;

  obs[4].obs_vec += Sz_stag / ham->latt->N;
  obs[5].obs_vec += (Sz_stag / ham->latt->N) * (Sz_stag / ham->latt->N);

  obs[6].obs_vec += state->beta * Sz * Sz / ham->latt->N;

  obs[7].obs_vec += Sz_stag_bulk / ham->latt->N;
  obs[8].obs_vec += (Sz_stag_bulk / ham->latt->N) * (Sz_stag_bulk / ham->latt->N);
  obs[9].obs_vec += Sz_stag_edge / ham->latt->N;
  obs[10].obs_vec += (Sz_stag_edge / ham->latt->N) * (Sz_stag_edge / ham->latt->N);
}

void sample_obs_eq(Obs_latt* obs, int n_eq, XXZ_ham* ham, SSE_config* state)
{
  int i, j, a, b, p;

  for (i = 0; i < n_eq; i++) {
    obs[i].N++;
  }

  for (i = 0; i < ham->latt->N; i++) {
    for (j = 0; j < ham->latt->N; j++) {
      obs[0].obs_latt[i][j] += state->spin_config[i] * state->spin_config[j] * 0.25;
    }
    obs[0].obs_latt0[i] += state->spin_config[i] * 0.5;
  }

  for (p = 0; p < state->n; p++) {
    a = state->reduced_op_string[p] % N_TYPES;
    b = (state->reduced_op_string[p] / N_TYPES) - 1;
    
    if (a == 1) {
      state->spin_config[ham->latt->bond_list[b][0]] += 2;
      state->spin_config[ham->latt->bond_list[b][1]] += -2;
    } else if (a == 2) {
      state->spin_config[ham->latt->bond_list[b][0]] += -2;
      state->spin_config[ham->latt->bond_list[b][1]] += 2;
    }

    obs[0].N++;
    for (i = 0; i < ham->latt->N; i++) {
      for (j = 0; j < ham->latt->N; j++) {
        obs[0].obs_latt[i][j] += state->spin_config[i] * state->spin_config[j] * 0.25;
      }
      obs[0].obs_latt0[i] += state->spin_config[i] * 0.5;
    }
  }
}

void free_observables(Obs_scalar* obs_scalar, int n_scal, Obs_latt* obs_eq, int n_eq)
{
  int n;

  for (n = 0; n < n_eq; n++) {
    free_obs_latt(&(obs_eq[n]));
  }
}

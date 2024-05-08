#include "sse.h"

void init_sse_config(double beta, int N, SSE_config *state)
{
  state->n = 0;
  state->M = N;
  state->n_loops = N;
  state->loop_size = 0;

  state->op_string = (int*) malloc(state->M * sizeof(int));
  state->op_tau = (double*) malloc((state->n + 2) * sizeof(double));
  state->reduced_op_string = (int*) malloc(state->n * sizeof(int));
  state->trans_op_string = (int*) malloc(state->n * sizeof(int));

  state->first = (int*) malloc(N * sizeof(int));
  state->vertex_list = (int*) malloc(state->n * sizeof(int));
  state->link = (int*) malloc(4 * state->n * sizeof(int)); 

  state->spin_config = (int*) malloc(N * sizeof(int));
  state->beta = beta;

  state->loop_size_2 = 0;
  state->n_loops_2 = 0;
  state->loop_histogram = (long long*) malloc(N * sizeof(long long));
  for (int i = 0; i < N; i++) {
    state->loop_histogram[i] = 0;
  }
}

void reset_sse_config(int N, double Sz, SSE_config *state) 
{
  memset(state->spin_config, Sz, N * sizeof(int));

  state->n = 0;
  state->M = N;
  state->loop_size = 0;
  state->n_loops = N;

  free(state->op_string);
  state->op_string = (int*) malloc(state->M * sizeof(int));
  memset(state->op_string, 0, state->M * sizeof(int));
}

void free_sse_config(SSE_config *state)
{
  free(state->spin_config);
  free(state->first);
  free(state->op_string);
  free(state->vertex_list);
  free(state->link);
  free(state->reduced_op_string);
  free(state->trans_op_string);

  if (state->op_tau != NULL) {
    free(state->op_tau);
  }
}

void diag_update(XXZ_ham* ham, SSE_config* state, pcg32_random_t* rng) 
{
  int p, a, b;
  int state1, state2;

  for (p = 0; p < state->M; p++) {
    if (state->op_string[p] == 0) {
      b = pcg32_boundedrand_r(rng, ham->latt->Nb);
      state1 = state->spin_config[ham->latt->bond_list[b][0]];
      state2 = state->spin_config[ham->latt->bond_list[b][1]];

      if (pcg32_double_r(rng) <= (ham->latt->Nb * state->beta * prob(state1, state2, ham)) / (state->M - state->n)) {
        state->op_string[p] = N_TYPES * (b + 1);
        state->n++;
      }
  } else if (state->op_string[p] % N_TYPES == 0) {
      b = (state->op_string[p] / N_TYPES) - 1;
      state1 = state->spin_config[ham->latt->bond_list[b][0]];
      state2 = state->spin_config[ham->latt->bond_list[b][1]];

      if (pcg32_double_r(rng) <= (state->M - state->n + 1) / (state->beta * ham->latt->Nb * prob(state1, state2, ham))) {
        state->op_string[p] = 0;
        state->n--;
      }
  } else {
      b = (state->op_string[p] / N_TYPES) - 1;
      a = state->op_string[p] % N_TYPES;

      if (a == 1) {
        state->spin_config[ham->latt->bond_list[b][0]] += 2;
        state->spin_config[ham->latt->bond_list[b][1]] += -2;
      } else if (a == 2) {
        state->spin_config[ham->latt->bond_list[b][0]] += -2;
        state->spin_config[ham->latt->bond_list[b][1]] += 2;
      }
    }
  }
}

void loop_update(XXZ_ham* ham, SSE_config* state, pcg32_random_t* rng, bool measure) 
{
  int loop, i, j, j0, p, li, le, l;
  int update, update_idx, old_state, new_state;
  long loop_size;
  double r;

  if (state->n == 0) {
    for (loop = 0; loop < state->n_loops; loop++) {
      for (i = 0; i < ham->latt->N; i++) {
        state->spin_config[i] = ham->Sz[pcg32_boundedrand_r(rng, ham->n_proj)];
      }
    }

    return;
  }

  for (loop = 0; loop < state->n_loops; loop++) {
    j0 = pcg32_boundedrand_r(rng, N_LEGS * state->n);
    j = j0;
    p = j / N_LEGS;
    li = j % N_LEGS;

/* SELECT A RANDOM UPDATE CHANGE THIS */
    update = 2;
    update_idx = 1;
    if (pcg32_double_r(rng) <= 0.5) { 
      update = -2; update_idx = 0; 
    }

    if (abs(ham->vertices[state->vertex_list[p]].spin[li] + update) > 2 * ham->S) {
      continue; 
    }

    loop_size = 0;
    while (true) {
      p = j / N_LEGS;
      li = j % N_LEGS;

      r = pcg32_double_r(rng);
      for (le = 0; le < N_LEGS; le++) {
        if (ham->vertices[state->vertex_list[p]].prob_exit[update_idx][li][le] == 1.0 || 
            r < ham->vertices[state->vertex_list[p]].prob_exit[update_idx][li][le]) {
          break;
        }
      } 

      old_state = ham->vertices[state->vertex_list[p]].spin[le];
      state->vertex_list[p] = ham->vertices[state->vertex_list[p]].new_vtx_type[update_idx][li][le];
      new_state = ham->vertices[state->vertex_list[p]].spin[le];
      if (new_state - old_state == -2) {
        update_idx = 0;
      } else if (new_state - old_state == 2) {
        update_idx = 1;
      } else if (new_state - old_state == 0) {
        if (update_idx == 1) { 
          update_idx = 0; 
        } else if (update_idx == 0) { 
          update_idx = 1; 
        }
      }

      if (li != le) { 
        loop_size++;
      }
      if (measure) {
        state->loop_histogram[ham->latt->bond_list[(state->reduced_op_string[p] / N_TYPES) - 1][0]]++;
        state->loop_histogram[ham->latt->bond_list[(state->reduced_op_string[p] / N_TYPES) - 1][1]]++;
      }

      j = N_LEGS * p + le;
      if (j == j0) { 
        break; 
      }
      
      j = state->link[j];
      if (j == j0) { 
        break; 
      }

      if (loop_size >= MAX_LOOP_SIZE) {
        printf("aborted loop update \n");
        fflush(stdout);
        return;
      }
    }

    if (measure) {
      state->loop_size_2 += loop_size;
      if (loop_size != 0) {
        state->n_loops_2++;
      }
    }
    state->loop_size += loop_size;
  }

  for (p = 0; p < state->n; p++) {
    state->reduced_op_string[p] = N_TYPES * (state->reduced_op_string[p] / N_TYPES) + ham->vertices[state->vertex_list[p]].type;
    state->op_string[state->trans_op_string[p]] = state->reduced_op_string[p];
  }

  for (i = 0; i < ham->latt->N; i++) {
    if (state->first[i] != -1) {
      p = state->first[i] / 4;
      l = state->first[i] % 4;
      state->spin_config[i] = ham->vertices[state->vertex_list[p]].spin[l];
    }
    else {
      state->spin_config[i] = ham->Sz[pcg32_boundedrand_r(rng, ham->n_proj)];
    }
  }
}

void ajust_cutoff(SSE_config *state, long t) 
{
  int* opstring_cpy;
  long M_new = state->n * 1.33;

  if (M_new > state->M) {
    opstring_cpy = (int*) malloc(state->M * sizeof(int));
    memcpy(opstring_cpy, state->op_string, state->M * sizeof(int));

    free(state->op_string); 
    state->op_string = (int*) malloc(M_new * sizeof(int));

    memset(state->op_string, 0, M_new * sizeof(int));
    memcpy(state->op_string, opstring_cpy, state->M * sizeof(int));
    free(opstring_cpy);

    state->M = M_new;
  }

  if ((t + 1) % 10 == 0) {
    if (state->loop_size > 10) {
      state->n_loops = 20 * state->M * state->n_loops / state->loop_size;   
    }
    if (state->n_loops < 4) {
      state->n_loops = 4;
    }

    state->loop_size = 0;
  }
}

void create_vtx_list(XXZ_ham* ham, SSE_config* state) 
{
  int p, p_red, b, a, v0, v1, v2, i, i1, i2;
  int last[ham->latt->N];
  int l[N_LEGS];
  
  free(state->vertex_list);
  free(state->link);
  free(state->reduced_op_string);
  free(state->trans_op_string);
  state->vertex_list = (int*) malloc(state->n * sizeof(int));
  state->link = (int*) malloc(4 * state->n * sizeof(int)); 
  state->reduced_op_string = (int*) malloc(state->n * sizeof(int));
  state->trans_op_string = (int*) malloc(state->n * sizeof(int));
  
  memset(last, -1, ham->latt->N * sizeof(int));
  memset(state->first, -1, ham->latt->N * sizeof(int));
  memset(state->link, -1, 4 * state->n * sizeof(int));
  memset(state->vertex_list, -1, state->n * sizeof(int));
  memset(state->reduced_op_string, 0, state->n * sizeof(int));
  memset(state->trans_op_string, 0, state->n * sizeof(int));
  
  p_red = 0;
  for (p = 0; p < state->M; p++) {
    if (state->op_string[p] != 0) {
      state->reduced_op_string[p_red] = state->op_string[p];
      state->trans_op_string[p_red] = p;
      p_red++;
    }
  }

  for (p = 0; p < state->n; p++) {
    b = (state->reduced_op_string[p] / 3) - 1;
    a = state->reduced_op_string[p] % 3;
    v0 = 4 * p;

    i1 = ham->latt->bond_list[b][0];
    i2 = ham->latt->bond_list[b][1];

    l[0] = state->spin_config[i1];
    l[1] = state->spin_config[i2];
    if (a == 1) {
      state->spin_config[i1] += 2;
      state->spin_config[i2] += -2;
    } else if (a == 2) {
      state->spin_config[i1] += -2;
      state->spin_config[i2] += 2;
    }
    l[2] = state->spin_config[i1];
    l[3] = state->spin_config[i2];

    for (i = 0; i < ham->n_diagrams; i++) {
      if (l[0] == ham->vertices[i].spin[0] && l[1] == ham->vertices[i].spin[1] && 
          l[2] == ham->vertices[i].spin[2] && l[3] == ham->vertices[i].spin[3]) { 
        state->vertex_list[p] = ham->vertices[i].indx; 
      }
    }
    
    v1 = last[i1];
    v2 = last[i2];

    if (v1 != -1) {
      state->link[v1] = v0;
      state->link[v0] = v1;
    } else {
      state->first[i1] = v0;
    }

    if (v2 != -1) {
      state->link[v2] = v0 + 1;
      state->link[v0 + 1] = v2;
    } else {
      state->first[i2] = v0 + 1;
    }

    last[i1] = v0 + 2;
    last[i2] = v0 + 3;
  }

  for (i = 0; i < ham->latt->N; i++) {
    if (state->first[i] != -1) {
      state->link[state->first[i]] = last[i];
      state->link[last[i]] = state->first[i];
    }
  }
}

double prob(int state1, int state2, XXZ_ham* ham) 
{
  for (int i = 0; i < ham->n_diagrams; i++) {
    if (ham->vertices[i].type == 0 && ham->vertices[i].spin[0] == state1 && ham->vertices[i].spin[1] == state2) {
      return ham->vertices[i].weight;
    }
  }

  return 0.0;
}

void assign_times(SSE_config* state, pcg32_random_t* rng)
{
  int i;

  free(state->op_tau);
  state->op_tau = (double*) malloc((state->n + 2) * sizeof(double));

  state->op_tau[0] = 0.0;
  for (i = 1; i <= state->n; i++) {
    state->op_tau[i] = pcg32_double_r(rng) * state->beta;    
  }
  state->op_tau[state->n + 1] = state->beta;
  qsort(state->op_tau, state->n + 2, sizeof(double), compare);
}

int compare( const void* num1, const void* num2)
{
    double a = *(double*) num1;  
    double b = *(double*) num2;  

    if(a > b)
    {  
        return 1;  
    }  
    else if(a < b)  
    {  
        return -1;  
    }  
    return 0;  
}

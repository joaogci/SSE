#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../sampling/sampling.h"
#include "../vtx/vtx_type.h"
#include "../tinyexpr/tinyexpr.h"

void read_inputs(char *file_name, int *d, int *L, double *S, double *delta, double *h, double *epsilon, long *therm_cycles, long *mc_cycles, int *n_bins, double **beta_vals, int *len_beta);
void read_vtx_info(char *file_name, vtx_element **vtx, int *n_diagrams, int *n_updates, int* n_legs);
void write_outputs(char *file_name, sampled_quantities *samples);

#endif // IO_H

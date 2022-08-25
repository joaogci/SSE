#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../sampling/sampling.h"

void read_inputs(char *file_name, int *d, int *L, double *J, double *delta, double *h, double *epsilon, long *therm_cycles, long *mc_cycles, int *n_bins, double **beta_vals, int *len_beta);
void write_outputs(char *file_name, sampled_quantities *samples);

#endif // IO_H
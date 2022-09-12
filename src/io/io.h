#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../sampling/sampling.h"
#include "../vtx/vtx_type.h"

void read_inputs(char *file_name, int *d, int *L, double *J, double *delta, double *h, double *epsilon, long *therm_cycles, long *mc_cycles, int *n_bins, double **beta_vals, int *len_beta);
void read_vtx_info(char *file_name, vtx_element **vtx, int *n_diagrams);
void write_outputs(char *file_name, sampled_quantities *samples);

#ifdef M_AUTOCORRELATION
    # include "../sampling/autocorrelation.h"
    void write_autocorrelation(char *file_name, autocorrelation *corr_series);
#endif

#endif // IO_H

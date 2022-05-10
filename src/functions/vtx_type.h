#ifndef VTX_TYPE_H
#define VTX_TYPE_H

#include <stdlib.h>

#define N_LEGS      4
#define N_DIAGRAMS  6

struct vtx_type
{
    int indx;
    int type;
    int spin[4];
    int new_vtx_type[4][4];
    double prob_exit[4][4];
    double H;
};

struct vtx_type *create_vtx_type_list(double J, double delta, double h, double C);
void spin_leg(int spin[N_LEGS], int vtx_type);
void new_vtx(int vtx[N_LEGS][N_LEGS], int vtx_type);

#endif // VTX_TYPE_H

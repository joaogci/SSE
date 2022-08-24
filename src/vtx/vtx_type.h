#ifndef VTX_TYPE_H
#define VTX_TYPE_H

#include <stdlib.h>

#define N_LEGS      4
#define N_DIAGRAMS  6

/* HEAT_BATH
    If equals 1 use heat-bath algorithm, else use the directed loop algorithm. */
#define HEAT_BATH   1

/* vtx_type
    Information about each vertex type (1 - 6).
    (int) indx - index of the vertex (1 - 6)
    (int) type - vertex type (0 - diagonal; 1 - off-diagonal)
    (int *) spin - spin value of each leg of the vertex (length 4)
    (int **) new_vtx_tpye - vertex type by flipping spins in l1 and l2 (lenght 4x4)
    (double **) prob_exit - probability of entering in leg l1 and exiting in l2 (length 4x4)
    (double) H - Hamiltonian matrix element for the vertex */
typedef struct vtx_type
{
    int indx;
    int type;
    int spin[N_LEGS];
    int new_vtx_type[N_LEGS][N_LEGS];
    double prob_exit[N_LEGS][N_LEGS];
    double H;
} vtx_element;

vtx_element *create_vtx_type_list(double J, double delta, double h, double epsilon);

void spin_leg(int spin[N_LEGS], int vtx_type);
void new_vtx(int vtx[N_LEGS][N_LEGS], int vtx_type);

#endif // VTX_TYPE_H

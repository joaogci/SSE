#ifndef VTX_TYPE_H
#define VTX_TYPE_H

#include <stdlib.h>

#define N_LEGS      4
#define N_UPDATES   2

/* 
 * struct: vtx_element
 *  information about each vertex type
 *  vertex information is generated by the python script
 *
 *  parameters:
 *      (int) indx: index of the vertex
 *      (int) type: vertex type (0 - SzSz; 1 - S+S-; 2 - S-S+)
 *      (int *) spin: spin value of each leg of the vertex (length N_LEGS)
 *      (int ***) new_vtx_tpye: vertex type by flipping spins in l1
 *          and l2 (lenght n_updates x N_LEGS x N_LEGS)
 *      (double ***) prob_exit: probability of entering in leg l1 
 *          and exiting in l2 (length n_updates x N_LEGS x N_LEGS)
 *      (double) H: Hamiltonian matrix element for the vertex 
 */
typedef struct vtx_type
{
    int indx;
    int type;
    int spin[N_LEGS];
    int new_vtx_type[N_UPDATES][N_LEGS][N_LEGS];
    double prob_exit[N_UPDATES][N_LEGS][N_LEGS];
    double H;
} vtx_element;

#endif // VTX_TYPE_H

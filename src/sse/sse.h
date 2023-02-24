#ifndef SSE_H
#define SSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>

#include "../vtx/vtx_type.h"
#include "../rng/pcg_basic.h"

#define MAX_(a, b) ((a) > (b) ? (a) : (b))
#define MIN_(a, b) ((a) < (b) ? (a) : (b))

#define C_ ( S * S * fabs(delta) + 2.0 * S * hb_ )
#define hb_ ( h / (2.0 * d) )

/* 
 * struct: heisenberg_system
 *  information about the simulated system
 *
 *  parameters
 *      (int) d: dimension
 *      (int) L: number of unit cells
 *      (int) N: number of particles (pow(L, d))
 *      (int) Nb: number of bonds 
 *      (int) boundary_cond: boundary_cond of the system
 *      (double) S: spin quantum number
 *      (int) n_proj: number of spin projections in the z-direction
 *      (int *) Sz: spin values in the z-direction
 *      (double) delta: z-axis anisotropy strength
 *      (double) h: z-axis magnetic field strength
 *      (double) epsilon: constant added to the Hamiltonian
 *      (int *) spin: spin state (length N)
 *      (int **) bond: lattice information (length Nb x 2) 
 */
typedef struct heisenberg_system 
{
    int d;
    int L;
    int N;
    int Nb;
    int boundary_cond;

    double S;
    int n_proj;
    int *Sz;
    
    double delta;
    double h;
    double epsilon;

    int *spin;
    int **bond;
} heisenberg_system;

/* 
 * struct: sse_state
 *  information about the state of the SSE simulation
 * 
 *  parameters:
 *      (int *) op_string: operator string (length M) 
 *      (int) M: maximum expansion of the partiton function
 *      (int) n: number of operators in the operator string
 *      (int) n_loops: number of loops to construct each MCS
 *      (int) loop_size: size of the loops
 *      (int *) link: linked list for the loop construction (length 4*n)
 *      (int *) first: first imaginary time which the spins appear (length N)
 *      (int *) vtx: vertex type at each imaginery time (length n)
 *      (vtx_element *) vtx_type: information about each vertex type (length n_diagrams) 
 *      (int) n_diagrams: number of available vtx_elements
 *      (int *) red_op_string: reduced operator string (length n)
 *      (int *) trans_op_string: translation between the reduced and normal operator string (length n) 
 */
typedef struct sse_state 
{
    int *op_string;
    int M;
    int n;

    int n_loops;
    int loop_size;
    int *link;
    int *first;
    int *vtx;
    vtx_element *vtx_type;
    int n_diagrams;

    int* red_op_string;
    int* trans_op_string;
} sse_state;

/* 
 * function: init_heisenberg_system 
 *  initializes the heisenberg_system struct
 * 
 *  parameters:
 *      (int) d: dimension
 *      (int) L: number of unit cells
 *      (int) boundary_cond: boundary condition of the lattice
 *      (double) S: spin quantum number
 *      (double) delta: z-axis anisotropy strength
 *      (double) h: z-axis magnetic field strength
 *      (double) epsilon: constant added to the Hamiltonian
 *      (heisenberg_system *) system: system to be initialized
 */
void init_heisenberg_system(int d, int L, int boundary_cond, double S, double delta, double h, 
    double epsilon, heisenberg_system *system);

/* 
 * function: init_sse_state 
 *  initializes the sse_state struct
 * 
 *  parameters:
 *      (uint64_t) seed: seed for the RNG
 *      (heisenberg_system *) system: simulated system
 *      (sse_state *) state: sse_state to be initialized
 */
void init_sse_state(uint64_t seed, heisenberg_system *system, sse_state *state);

/* 
 * function: reset_sse_state
 *  resets the SSE state
 *  deletes all of the operator string and spin states
 * 
 *  parameters:
 *      (heisenberg_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void reset_sse_state(heisenberg_system *system, sse_state *state);

/* 
 * function: diag_update
 *  diagonal update for the SSE MCS
 *  inserts or removes a diagonal operator from the operator string
 *   with some probablity
 * 
 *  paramters:
 *      (double) beta: temperature
 *      (heisenberg_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void diag_update(double beta, heisenberg_system *system, sse_state *state, pcg32_random_t* rng);

/* 
 * function: loop_update
 *  loop update for the SSE MCS
 *  creates a loop within the vertex list and updates the configuration
 *  off-diagonal update
 * 
 *  parameters:
 *      (heisenberg_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void loop_update(heisenberg_system *system, sse_state *state, pcg32_random_t* rng);

/* 
 * function: ajust_cutoff
 *  dinamically adjusts the expansion cutoff during the thermalization
 *   part of the simulation
 * 
 *  parameters:
 *      (sse_state *) state: SSE state
 *      (bool) adjust_loop
 */
void ajust_cutoff(sse_state *state, bool adjust_loop);

/* 
 * function: create_vtx_list
 *  tranlates the operator string in to a doubly linked list of vertecies
 * 
 *  parameters:
 *      (heisenberg_system *) system: simulated system
 *      (sse_state *) state: SSE state
 */
void create_vtx_list(heisenberg_system *system, sse_state *state);

/* 
 * function: prob
 *  computes probability for diagonal updates
 *  
 *  parameters:
 *      (int) b: bond for the insertion/removal of operator
 *      (heisenberg_system *) system: system
 *      (sse_state* ) state: SSE state
 *  
 *  returns:
 *      (double) prob: probability for the update
 */
double prob(int b, heisenberg_system *system, sse_state *state);

/* 
 * function: free_memory
 *  frees allocated memory in the heisenberg_system and sse_state structs
 * 
 *  parameters:
 *      (heisenberg_system *) system: heisenberg_system struct
 *      (sse_state *) state: sse_state struct
 */
void free_memory(heisenberg_system *system, sse_state *state);

#endif // SSE_H

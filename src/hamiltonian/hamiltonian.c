#include "hamiltonian.h"

H_mat_term **create_hamiltonian(double J, double delta, double h, double C) 
{
    H_mat_term **H = (H_mat_term **) malloc(4 * sizeof(H_mat_term *));

    for (int i = 0; i < 4; i++) {
        H[i] = (H_mat_term *) malloc(4 * sizeof(H_mat_term));
    }

    H[0][0].value = C - 0.25 * delta + 0.5 * h / J;
    H[0][0].si_p = 1;
    H[0][0].sj_p = 1;
    H[0][0].si_p_1 = 1;
    H[0][0].sj_p_1 = 1;

    H[1][1].value = C + 0.25 * delta;
    H[1][1].si_p = 1;
    H[1][1].sj_p = -1;
    H[1][1].si_p_1 = 1;
    H[1][1].sj_p_1 = -1;

    H[2][2].value = C + 0.25 * delta;
    H[2][2].si_p = -1;
    H[2][2].sj_p = 1;
    H[2][2].si_p_1 = -1;
    H[2][2].sj_p_1 = 1;

    H[3][3].value = C - 0.25 * delta - 0.5 * h / J;
    H[3][3].si_p = -1;
    H[3][3].sj_p = -1;
    H[3][3].si_p_1 = -1;
    H[3][3].sj_p_1 = -1;

    H[1][2].value = 0.5;
    H[1][2].si_p = 1;
    H[1][2].sj_p = -1;
    H[1][2].si_p_1 = -1;
    H[1][2].sj_p_1 = 1;

    H[2][1].value = 0.5;
    H[2][1].si_p = -1;
    H[2][1].sj_p = 1;
    H[2][1].si_p_1 = 1;
    H[2][1].sj_p_1 = -1;

    return H;
}

#include "vtx_type.h"

struct vtx_type *create_vtx_type_list(double J, double delta, double h, double C) {
    int i;
    struct vtx_type *vtx = (struct vtx_type *) malloc(6 * sizeof(struct vtx_type));

    for (i = 0; i < N_DIAGRAMS; i++) {
        vtx[i].indx = i + 1;
        vtx[i].type = 0;
        if (i == 3 || i == 4) {
            vtx[i].type = 1;
        }

        if (vtx[i].type == 0) {
            if (vtx[i].indx == 6) {
                vtx[i].H = C - 0.25 * delta + 0.5 * h / J;
            }
            else if (vtx[i].indx == 1) {
                vtx[i].H = C - 0.25 * delta - 0.5 * h / J;
            }
            else {
                vtx[i].H = C + 0.25 * delta;
            }
        }
        else {
            vtx[i].H = 0.5;
        }

        new_vtx(vtx[i].new_vtx_type, i + 1);
        spin_leg(vtx[i].spin, i + 1);
    }

    double denom;
    int new_vtx, new_vtx2;
    int li, le, le2;
    for (i = 0; i < N_DIAGRAMS; i++) {
        for (li = 0; li < N_LEGS; li++) {
            for (le = 0; le < N_LEGS; le++) {
                new_vtx = vtx[i].new_vtx_type[li][le] - 1;
                if (new_vtx >= 0) {
                    vtx[i].prob_exit[li][le] = vtx[new_vtx].H;

                    denom = 0.0;
                    for (le2 = 0; le2 < N_LEGS; le2++) {
                        new_vtx2 = vtx[i].new_vtx_type[li][le2] - 1;
                        if (new_vtx2 >= 0) {
                            denom += vtx[new_vtx2].H;
                        }
                    }

                    vtx[i].prob_exit[li][le] /= denom;
                }
                else 
                {
                    vtx[i].prob_exit[li][le] = 0.0;
                }
            }
        }
    }

    return vtx;
}

void spin_leg(int spin[N_LEGS], int vtx_type) {
    for (int l = 1; l <= N_LEGS; l++) {
        if (vtx_type == 1) {
            spin[l - 1] = -1;
        }
        else if (vtx_type == 2) {
            if (l % 2 == 0)
                spin[l - 1] = 1;
            else 
                spin[l - 1] = -1;
        }
        else if (vtx_type == 3) {
            if (l % 2 == 0)
                spin[l - 1] = -1;
            else 
                spin[l - 1]= 1;
        }
        else if (vtx_type == 4) {
            if (l == 1 || l == 4)
                spin[l - 1] = -1;
            else 
                spin[l - 1] = 1;
        }
        else if (vtx_type == 5) {
            if (l == 1 || l == 4)
                spin[l - 1] = 1;
            else 
                spin[l - 1] = -1;
        }
        else if (vtx_type == 6) {
            spin[l - 1] = 1;
        }
    }
}

void new_vtx(int vtx[N_LEGS][N_LEGS], int vtx_type) {
    for (int li = 1; li <= N_LEGS; li++) {
        for (int le = 1; le <= N_LEGS; le++) {
            vtx[li - 1][le - 1] = 0;
            if (li == le) {
                vtx[li - 1][le - 1] = vtx_type;
            }
            else {
                if (vtx_type == 1) {
                    if (li == 1) {
                        if (le == 3)
                            vtx[li - 1][le - 1] = 3;
                        else if (le == 4)
                            vtx[li - 1][le - 1] = 5;
                    }
                    else if (li == 2) {
                        if (le == 4) 
                            vtx[li - 1][le - 1] = 2;
                        else if (le == 3)
                            vtx[li - 1][le - 1] = 4;
                    }
                    else if (li == 3) {
                        if (le == 1)
                            vtx[li - 1][le - 1] = 3;
                        else if (le == 2)
                            vtx[li - 1][le - 1] = 4;
                    }
                    else if (li == 4) {
                        if (le == 2)
                            vtx[li - 1][le - 1] = 2;
                        else if (le == 1)
                            vtx[li - 1][le - 1] = 5;
                    }
                }
                else if (vtx_type == 2) {
                    if (li == 1) {
                        if (le == 3)
                            vtx[li - 1][le - 1] = 6;
                        else if (le == 2)
                            vtx[li - 1][le - 1] = 5;
                    }
                    else if (li == 2) {
                        if (le == 4)
                            vtx[li - 1][le - 1] = 1;
                        else if (le == 1)
                            vtx[li - 1][le - 1] = 5;
                    }
                    else if (li == 3) {
                        if (le == 1)
                            vtx[li - 1][le - 1] = 6;
                        else if (le == 4)
                            vtx[li - 1][le - 1] = 4;
                    }
                    else if (li == 4) {
                        if (le == 2)
                            vtx[li - 1][le - 1] = 1;
                        else if (le == 3)
                            vtx[li - 1][le - 1] = 4;
                    }
                } 
                else if (vtx_type == 3) {
                    if (li == 1) {
                        if (le == 3)
                            vtx[li - 1][le - 1] = 1;
                        else if (le == 2)
                            vtx[li - 1][le - 1] = 4;
                    }
                    else if (li == 2) {
                        if (le == 4)
                            vtx[li - 1][le - 1] = 6;
                        else if (le == 1)
                            vtx[li - 1][le - 1] = 4;
                    }
                    else if (li == 3) {
                        if (le == 1)
                            vtx[li - 1][le - 1] = 1;
                        else if (le == 4)
                            vtx[li - 1][le - 1] = 5;
                    }
                    else if (li == 4) {
                        if (le == 2)
                            vtx[li - 1][le - 1] = 6;
                        else if (le == 3)
                            vtx[li - 1][le - 1] = 5;
                    }
                }
                else if (vtx_type == 4) {
                    if (li == 1) {
                        if (le == 4)
                            vtx[li - 1][le - 1] = 6;
                        else if (le == 2)
                            vtx[li - 1][le - 1] = 3;
                    }
                    else if (li == 2) {
                        if (le == 3)
                            vtx[li - 1][le - 1] = 1;
                        else if (le == 1)
                            vtx[li - 1][le - 1] = 3;
                    }
                    else if (li == 3) {
                        if (le == 2)
                            vtx[li - 1][le - 1] = 1;
                        else if (le == 4)
                            vtx[li - 1][le - 1] = 2;
                    }
                    else if (li == 4) {
                        if (le == 1)
                            vtx[li - 1][le - 1] = 6;
                        else if (le == 3)
                            vtx[li - 1][le - 1] = 2;
                    }
                }
                else if (vtx_type == 5) {
                    if (li == 1) {
                        if (le == 4)
                            vtx[li - 1][le - 1] = 1;
                        else if (le == 2)
                            vtx[li - 1][le - 1] = 2;
                    }
                    else if (li == 2) {
                        if (le == 3)
                            vtx[li - 1][le - 1] = 6;
                        else if (le == 1)
                            vtx[li - 1][le - 1] = 2;
                    }
                    else if (li == 3) {
                        if (le == 2)
                            vtx[li - 1][le - 1] = 6;
                        else if (le == 4)
                            vtx[li - 1][le - 1] = 3;
                    }
                    else if (li == 4) {
                        if (le == 1)
                            vtx[li - 1][le - 1] = 1;
                        else if (le == 3)
                            vtx[li - 1][le - 1] = 3;
                    }
                }
                else if (vtx_type == 6) {
                    if (li == 1) {
                        if (le == 3)
                            vtx[li - 1][le - 1] = 2;
                        else if (le == 4)
                            vtx[li - 1][le - 1] = 4;
                    }
                    else if (li == 2) {
                        if (le == 4)
                            vtx[li - 1][le - 1] = 3;
                        else if (le == 3)
                            vtx[li - 1][le - 1] = 5;
                    }
                    else if (li == 3) {
                        if (le == 1)
                            vtx[li - 1][le - 1] = 2;
                        else if (le == 2)
                            vtx[li - 1][le - 1] = 5;
                    }
                    else if (li == 4) {
                        if (le == 2)
                            vtx[li - 1][le - 1] = 3;
                        else if (le == 1)
                            vtx[li - 1][le - 1] = 4;
                    }
                }
            }
        }
    }
}

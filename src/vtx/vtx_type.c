#include "vtx_type.h"

void free_vtx_type(vtx_element* vtx_type, int n_legs, int n_updates)
{
    free(vtx_type->spin);

    for (int j = 0; j < n_updates; j++) {
        for (int l = 0; l < n_legs; l++) {
            free(vtx_type->prob_exit_[j][l]);
            free(vtx_type->new_vtx_type_[j][l]);
        }
        free(vtx_type->prob_exit_[j]);
        free(vtx_type->new_vtx_type_[j]);
    }
    free(vtx_type->prob_exit_);
    free(vtx_type->new_vtx_type_);
}

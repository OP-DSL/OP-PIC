/* 
BSD 3-Clause License

Copyright (c) 2022, OP-DSL

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#pragma once

//*********************************************
// USER WRITTEN CODE
//*********************************************

#include "opp_lib.h"

#define KERNEL_ONE_OVER_SIX     (1.0 / 6.0)
#define KERNEL_N_PER_C          4
#define KERNEL_DET_FIELDS       4
#define KERNEL_NEIGHB_C         4
#define KERNEL_DIM              3

//*************************************************************************************************
inline void init_boundary_pot_kernel(
    const int *node_type, 
    double *n_bnd_pot
) {
    switch (*node_type) {
        case 2: // INLET: 
            *n_bnd_pot = 0; break;     
        case 3: // FIXED: 
            *n_bnd_pot = -1 * CONST_wall_potential[0]; break;
        default: // NORMAL or OPEN
            *n_bnd_pot = 0; /*default*/
    }
}

//*************************************************************************************************
inline void inject_ions_kernel(
    double *part_pos,
    double *part_vel,
    int *part_cell_connectivity,
    const int *cell_id, 
    const double *cell_ef,
    const double *iface_u,
    const double *iface_v,
    const double *iface_normal,
    const double *node_pos,
    const double* dummy_part_random
) {
    double a = dummy_part_random[0];
    double b = dummy_part_random[1];
    if ((a + b) > 1) {  // TODO : Change the random dat to avoid this
        a = (1 - a);
        b = (1 - b);
    }

    for (int i = 0; i < KERNEL_DIM; i++) {
        part_pos[i] = a * iface_u[i] + b * iface_v[i] + node_pos[i];
        
        part_vel[i] = (iface_normal[i] * CONST_ion_velocity[0]);
        part_vel[i] -= CONST_charge[0] / CONST_mass[0] * cell_ef[i] * (0.5 * CONST_dt[0]);
    }

    (*part_cell_connectivity) = (*cell_id);
}

//*************************************************************************************************
inline void calculate_new_pos_vel_kernel(
    const double *cell_ef,
    double *part_pos,
    double *part_vel 
) {
    const double coefficient1 = CONST_charge[0] / CONST_mass[0] * (CONST_dt[0]);
    for (int i = 0; i < KERNEL_DIM; i++) {
        part_vel[i] += (coefficient1 * cell_ef[i]);   
        part_pos[i] += part_vel[i] * (CONST_dt[0]);                
    }
}

//*************************************************************************************************
inline void move_kernel(
    const double *point_pos, double* point_lc, 
    const double *cell_volume, const double *cell_det
) {
    const double coefficient2 = KERNEL_ONE_OVER_SIX / (*cell_volume);

    for (int i=0; i<KERNEL_N_PER_C; i++) { /*loop over vertices*/
    
        point_lc[i] = coefficient2 * (
            cell_det[i * KERNEL_DET_FIELDS + 0] - 
            cell_det[i * KERNEL_DET_FIELDS + 1] * point_pos[0] + 
            cell_det[i * KERNEL_DET_FIELDS + 2] * point_pos[1] - 
            cell_det[i * KERNEL_DET_FIELDS + 3] * point_pos[2]);
    }  

    if (!(point_lc[0] < 0.0 || point_lc[0] > 1.0 ||
          point_lc[1] < 0.0 || point_lc[1] > 1.0 ||
          point_lc[2] < 0.0 || point_lc[2] > 1.0 ||
          point_lc[3] < 0.0 || point_lc[3] > 1.0)) { 
            
        OPP_PARTICLE_MOVE_DONE;
        return;
    } 

    // outside the last known cell, find most negative weight and 
    // use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = point_lc[0];
    
    for (int i=1; i<KERNEL_NEIGHB_C; i++) {
        if (point_lc[i] < min_lc) {
            min_lc = point_lc[i];
            min_i = i;
        }
    }

    if (opp_c2c[min_i] >= 0) { // is there a neighbor in this direction?
        (*opp_p2c) = opp_c2c[min_i];
        OPP_PARTICLE_NEED_MOVE;
    }
    else {
        (*opp_p2c) = MAX_CELL_INDEX;
        OPP_PARTICLE_NEED_REMOVE;
    }
}

//*************************************************************************************************
inline void deposit_charge_on_nodes_kernel(
    const double *part_lc,
    double *node_charge_den0,
    double *node_charge_den1,
    double *node_charge_den2,
    double *node_charge_den3
) {
    node_charge_den0[0] += part_lc[0];
    node_charge_den1[0] += part_lc[1];
    node_charge_den2[0] += part_lc[2];
    node_charge_den3[0] += part_lc[3];
}

//*************************************************************************************************
inline void compute_node_charge_density_kernel(
    double *node_charge_den,
    const double *node_volume
) {
    (*node_charge_den) *= (CONST_spwt[0] / node_volume[0]);
}

//*************************************************************************************************
inline void compute_electric_field_kernel(
    double *cell_electric_field,             
    const double *cell_shape_deriv,
    const double *node_potential0,
    const double *node_potential1,
    const double *node_potential2,
    const double *node_potential3
) {
    for (int dim = 0; dim < KERNEL_DIM; dim++) { 
        const double c1 = (cell_shape_deriv[0 * KERNEL_DIM + dim] * (*node_potential0));
        const double c2 = (cell_shape_deriv[1 * KERNEL_DIM + dim] * (*node_potential1));
        const double c3 = (cell_shape_deriv[2 * KERNEL_DIM + dim] * (*node_potential2));
        const double c4 = (cell_shape_deriv[3 * KERNEL_DIM + dim] * (*node_potential3));

        cell_electric_field[dim] -= (c1 + c2 + c3 + c4);
    }    
}

//*************************************************************************************************
inline void get_final_max_values_kernel(
    const OPP_REAL* n_charge_den,
    OPP_REAL* max_n_charge_den,
    const OPP_REAL* n_pot,
    OPP_REAL* max_n_pot
) {
    *max_n_charge_den = MAX(abs(*n_charge_den), *max_n_charge_den);
    *max_n_pot = MAX(*n_pot, *max_n_pot);
}

//*************************************************************************************************
inline void get_max_cef_kernel(
    const OPP_REAL* val,
    OPP_REAL* max_val
) {
    for (int dim = 0; dim < KERNEL_DIM; ++dim) {
        *max_val = MAX(val[dim], *max_val);
    }
}

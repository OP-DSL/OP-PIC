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

// AUTO GENERATED CODE

#include "oppic_seq.h"
#include "../fempic.h"

double OPP_VAR_spwt;                // todo : these should be opp_consts
double OPP_VAR_ion_velocity;
double OPP_VAR_dt;
double OPP_VAR_plasma_den;
double OPP_VAR_mass;
double OPP_VAR_charge;
move_status OPP_move_status = OPP_MOVE_DONE; // this should be in lib level
bool OPP_inside_cell = true;
bool OPP_iteration_one = true;

#include "../kernels.h"


//*************************************************************************************************
void oppic_seq_loop_inject__Increase_particle_count
(
    oppic_set particles_set,    // particles_set
    oppic_set set,              // inlect_face_set
    oppic_arg arg0,             // injected total global,
    oppic_arg arg1,             // iface_area,
    oppic_arg arg2,             // iface_inj_part_dist,
    oppic_arg arg3              // remainder global,
)
{ TRACE_ME;
    
    OPP_VAR_plasma_den = opp_params->get<REAL>("plasma_den"); // todo : these should be opp_consts and should be populated by runtime
    OPP_VAR_dt = opp_params->get<REAL>("dt");
    OPP_VAR_ion_velocity = opp_params->get<REAL>("ion_velocity");
    OPP_VAR_spwt = 2e2;
    OPP_VAR_mass = 2 * AMU;
    OPP_VAR_charge = 1 * QE;

    if (OP_DEBUG) printf("FEMPIC - oppic_seq_loop_inject__Increase_particle_count num_particles %d diff %d\n", set->size, set->diff);

    for (int i = 0; i < set->size; i++)
    {   
        calculate_injection_distribution(
            ((int *)arg0.data),
            &((double *)arg1.data)[i],
            &((int *)arg2.data)[i],
            ((double *)arg3.data) 
        );
    }

    oppic_increase_particle_count(particles_set, *((int *)arg0.data));

    int* part_mesh_connectivity = (int *)particles_set->cell_index_dat->data;
    int* distribution           = (int *)arg2.data;

    int start = (particles_set->size - particles_set->diff);
    int j = 0;

    for (int i = 0; i < particles_set->diff; i++)
    {
        if (i >= distribution[j]) j++; // check whether it is j or j-1    
        part_mesh_connectivity[start + i] = j;
    }   
}

//*************************************************************************************************
void oppic_par_loop_inject__InjectIons(
    oppic_set set,      // particles_set
    oppic_arg arg0,     // part_position,
    oppic_arg arg1,     // part_velocity,
    oppic_arg arg2,     // part_cell_connectivity,
    oppic_arg arg3,     // iface to cell map
    oppic_arg arg4,     // cell_ef,
    oppic_arg arg5,     // iface_u,
    oppic_arg arg6,     // iface_v,
    oppic_arg arg7,     // iface_normal,
    oppic_arg arg8      // iface_node_pos
)
{ TRACE_ME;

    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_inject__InjectIons num_particles %d diff %d\n", set->size, set->diff);

    int start = (set->size - set->diff);

    for (int i = 0; i < set->diff; i++)
    {    
        int map0idx    = ((int *)set->cell_index_dat->data)[(start + i) * set->cell_index_dat->dim]; // iface index

        const int map1idx = arg4.map_data[map0idx]; // cell index

        inject_ions__kernel(
            &((double *)arg0.data)[(start + i) * arg0.dim],    // part_position,
            &((double *)arg1.data)[(start + i) * arg1.dim],    // part_velocity,
            &((int *)arg2.data)[(start + i) * arg2.dim],       // part_cell_connectivity,
            &((int *)arg3.data)[map0idx * arg3.dim],           // iface to cell map
            &((double*)arg4.data)[map1idx * arg4.dim],         // cell_ef,
            &((double*)arg5.data)[map0idx * arg5.dim],         // iface_u,
            &((double*)arg6.data)[map0idx * arg6.dim],         // iface_v,
            &((double*)arg7.data)[map0idx * arg7.dim],         // iface_normal,
            &((double*)arg8.data)[map0idx * arg8.dim]          // iface_node_pos
        );
    }
}

//*************************************************************************************************
void oppic_par_loop_all__ResetIonDensity(
    oppic_set set,     // nodes_set
    oppic_arg arg0     // node_charge_density
)
{ TRACE_ME;
    
    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__ResetIonDensity num_nodes %d\n", set->size);

    for (int i=0; i<set->size; i++) 
    {
        reset_ion_density__kernel(
            &((double*)arg0.data)[i * arg0.dim]
        );
    }
}

//*************************************************************************************************
void oppic_par_loop_particle_all__MoveToCells(
    oppic_set set,      // particles_set
    oppic_arg arg0,     // cell_ef,
    oppic_arg arg1,     // part_pos,
    oppic_arg arg2,     // part_vel,
    oppic_arg arg3,     // part_lc,
    oppic_arg arg4,     // current_cell_index,
    oppic_arg arg5,     // current_cell_volume,
    oppic_arg arg6,     // current_cell_det,
    oppic_arg arg7,     // cell_connectivity,
    oppic_arg arg8,     // node_charge_den0,
    oppic_arg arg9,     // node_charge_den1,
    oppic_arg arg10,    // node_charge_den2,
    oppic_arg arg11     // node_charge_den3,
)
{ TRACE_ME;

    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_particle_all__MoveToCells num_particles %d diff %d\n", set->size, set->diff);

    oppic_init_particle_move(set);

    for (int i = 0; i < set->size; i++)
    {        
        OPP_move_status = OPP_NEED_MOVE;
        OPP_iteration_one = true;

        do
        { 
            OPP_inside_cell = true;

            int& map0idx      = ((int *)set->cell_index_dat->data)[i * set->cell_index_dat->dim];    // this is the cell_index

            const int map1idx = arg8.map_data[map0idx * arg8.map->dim + 0];
            const int map2idx = arg8.map_data[map0idx * arg8.map->dim + 1];
            const int map3idx = arg8.map_data[map0idx * arg8.map->dim + 2];
            const int map4idx = arg8.map_data[map0idx * arg8.map->dim + 3];

            move_all_particles_to_cell__kernel(
                &((double *)arg0.data)[map0idx * arg0.dim], // const double *cell_ef,
                &((double *)arg1.data)[i * arg1.dim],       // double *part_pos,
                &((double *)arg2.data)[i * arg2.dim],       // double *part_vel,
                &((double *)arg3.data)[i * arg3.dim],       // double *part_lc,
                &((int *)arg4.data)[i * arg4.dim],          // int* current_cell_index,
                &((double*)arg5.data)[map0idx * arg5.dim],  // const double *current_cell_volume,
                &((double*)arg6.data)[map0idx * arg6.dim],  // const double *current_cell_det,
                &((int*)arg7.data)[map0idx * arg7.dim],     // const int *cell_connectivity,
                &((double*)arg8.data)[map1idx],             // double *node_charge_den0,
                &((double*)arg8.data)[map2idx],             // double *node_charge_den1,
                &((double*)arg8.data)[map3idx],             // double *node_charge_den2,
                &((double*)arg8.data)[map4idx]              // double *node_charge_den3,
            );             
            
            OPP_iteration_one = false;

        } while (OPP_move_status == (int)OPP_NEED_MOVE);

        oppic_mark_particle_to_move(set, i, (int)OPP_move_status); // send the enum directly
    }

    oppic_finalize_particle_move(set);
}

//*************************************************************************************************
void oppic_par_loop_all__ComputeNodeChargeDensity(
    oppic_set set,     // nodes_set
    oppic_arg arg0,    // node_charge_density
    oppic_arg arg1     // node_volume
)
{ TRACE_ME;
    
    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__ComputeNodeChargeDensity num_nodes %d\n", set->size);

    for (int i=0; i<set->size; i++) 
    {
        compute_node_charge_density__kernel(
            &((double*)arg0.data)[i * arg0.dim],
            &((double*)arg1.data)[i * arg1.dim]
        );
    }
}

//*************************************************************************************************
void oppic_par_loop_all__ComputeElectricField(
    oppic_set set,      //cells_set,                                                        // cells_set
    oppic_arg arg0,     //oppic_arg_dat(cell_electric_field,                  OP_INC),      // cell_electric_field,
    oppic_arg arg1,     //oppic_arg_dat(cell_shape_deriv,                     OP_READ),     // cell_shape_deriv,
    oppic_arg arg2,     //oppic_arg_dat(node_potential, 0, cell_to_nodes_map, OP_READ),     // node_potential0,
    oppic_arg arg3,     //oppic_arg_dat(node_potential, 1, cell_to_nodes_map, OP_READ),     // node_potential1,
    oppic_arg arg4,     //oppic_arg_dat(node_potential, 2, cell_to_nodes_map, OP_READ),     // node_potential2,
    oppic_arg arg5      //oppic_arg_dat(node_potential, 3, cell_to_nodes_map, OP_READ)      // node_potential3,
)
{
    if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__ComputeElectricField size %d\n", set->size);

    for (int i = 0; i < set->size; i++)
    {
        const int map1idx = arg2.map_data[i * arg2.map->dim + 0];
        const int map2idx = arg2.map_data[i * arg2.map->dim + 1];
        const int map3idx = arg2.map_data[i * arg2.map->dim + 2];
        const int map4idx = arg2.map_data[i * arg2.map->dim + 3];

        compute_electric_field__kernel(
            &((double*)arg0.data)[i * arg0.dim],    // cell_electric_field
            &((double*)arg1.data)[i * arg1.dim],    // cell_shape_deriv
            &((double*)arg2.data)[map1idx],         // node_potential0
            &((double*)arg2.data)[map2idx],         // node_potential1
            &((double*)arg2.data)[map3idx],         // node_potential2
            &((double*)arg2.data)[map4idx]          // node_potential3
        );
    }   
}

// //*************************************************************************************************
// oppic_par_loop_all__MoveParticles(
//     particles_set,
//     oppic_arg_dat(part_position,       OP_INC),
//     oppic_arg_dat(part_velocity,       OP_INC),
//     oppic_arg_dat(cell_electric_field, OP_READ, true)
// );

// void oppic_par_loop_all__MoveParticles(
//     oppic_set set,     // particles_set
//     oppic_arg arg0,    // part_position,
//     oppic_arg arg1,    // part_velocity,
//     oppic_arg arg2     // cell_electric_field,
//     )
// { TRACE_ME;
    
//     if (OP_DEBUG) printf("FEMPIC - oppic_par_loop_all__MoveParticles num_particles %d\n", set->size);

//     for (int i = 0; i < set->size; i++)
//     {
//         int map0idx    = ((int *)set->cell_index_dat->data)[i * set->cell_index_dat->dim];

//         move_particles__kernel(
//             &((double *)arg0.data)[i * arg0.dim],          // part_position,
//             &((double *)arg1.data)[i * arg1.dim],          // part_velocity,
//             &((double *)arg2.data)[map0idx * arg2.dim]     // cell_electric_field
//         );
//     }
// }
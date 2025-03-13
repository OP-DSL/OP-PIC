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

//*********************************************
// AUTO GENERATED CODE
//*********************************************

#include "opp_seq.h"

//****************************************
// double CONST_spwt = 0;
void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    // These will be created when opp_decl_const<>() is used in simpic.cpp

    // if (!strcmp(name,"CONST_spwt"))              CONST_spwt = *((double*)data);
    // else std::cerr << "error: unknown const name" << std::endl;
}
//****************************************

#include "../kernels.h"

//*************************************************************************************************
void opp_par_loop_all__WeightFieldsToParticles(
    opp_set set,     // particles
    opp_arg arg0,    // lhs node0_field_E
    opp_arg arg1,    // rhs node1_field_E
    opp_arg arg2,    // particle0_position_x
    opp_arg arg3     // particle0_field_E
    )
{ 

    if (OPP_DBG) printf("SIMPIC - opp_par_loop_all__WeightFieldsToParticles [%d]\n", set->size);

    const int set_size = set->size;
    for (int n = 0; n < set_size; n++)
    {
        const int map0idx    = ((int *)set->mesh_relation_dat->data)[n]; // this is the cell_index

        const int map1idx = arg0.map_data[map0idx * arg0.map->dim + 0]; //LHS node
        const int map2idx = arg1.map_data[map0idx * arg1.map->dim + 1]; //RHS node

        weight_fields_to_particles__kernel(
            &((double*)arg0.data)[map1idx],
            &((double*)arg1.data)[map2idx],
            &((double*)arg2.data)[n],
            &((double*)arg3.data)[n]);
    }
}

//*************************************************************************************************
void opp_par_loop_particle_all__PushParticles( 
    opp_set set,     // particles
    opp_arg arg0,    // particle0_field_E
    opp_arg arg1,    // particle0_velocity_x
    opp_arg arg2,    // particle0_position_x
    opp_arg arg3     // particle0_cell_index
    )
{ 

    if (OPP_DBG) printf("SIMPIC - opp_par_loop_particle_all__PushParticles [%d]\n", set->size);

    opp_init_particle_move(set, 0, nullptr);

    OPP_INT* cellIdx = nullptr;

    const int set_size = set->size;
    for (int i = 0; i < set_size; i++)
    {              
        OPP_MOVE_RESET_FLAGS;

        do
        { 
            cellIdx = &(OPP_mesh_relation_data[i]);

            particle_push__kernel(
                &((double*)arg0.data)[i],
                &((double*)arg1.data)[i],
                &((double*)arg2.data)[i],
                &((int*)arg3.data)[i]
            );

        } while (opp_check_part_move_status(*cellIdx, i, set->particle_remove_count));
    }

    opp_finalize_particle_move(set);
}

//*************************************************************************************************
void opp_par_loop_all__WeightParticlesToFields(
    opp_set set,     // particles
    opp_arg arg0,    // lhs node0_field_J
    opp_arg arg1,    // rhs node1_field_J
    opp_arg arg2     // particle0_position_x
    )
{ 

    if (OPP_DBG) printf("SIMPIC - opp_par_loop_all__WeightParticlesToFields [%d]\n", set->size);

    const int set_size = set->size;
    for (int n = 0; n < set_size; n++)
    {
        const int map0idx    = ((int *)set->mesh_relation_dat->data)[n]; // this is the cell_index

        const int map1idx = arg0.map_data[map0idx * arg0.map->dim + 0]; //LHS node
        const int map2idx = arg1.map_data[map0idx * arg1.map->dim + 1]; //RHS node

        weight_particles_to_fields__kernel(
            &((double*)arg0.data)[map1idx],
            &((double*)arg1.data)[map2idx],
            &((double*)arg2.data)[n]);
    }
}

//*************************************************************************************************
void opp_par_loop_all__FieldSolveSumLaplace(
    opp_set set,   // nodes
    opp_arg arg0,  // node0_xlocal
    opp_arg arg1)  // node0_field_P
{ 

    if (OPP_DBG) printf("SIMPIC - opp_par_loop__FieldSolveSumLaplace [%d]\n", set->size);

    const int set_size = set->size;
    for(int n = 0; n < set_size; n++)
    {
        field_solve_sum_laplace__kernel(
            &((double*)arg0.data)[n],
            &((double*)arg1.data)[n]);
    }
}

//*************************************************************************************************
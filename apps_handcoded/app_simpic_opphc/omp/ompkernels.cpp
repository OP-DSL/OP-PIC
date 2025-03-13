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


#include "opp_omp.h"

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

}

//*************************************************************************************************
void opp_par_loop_all__FieldSolveSumLaplace(
    opp_set set,   // nodes
    opp_arg arg0,  // node0_xlocal
    opp_arg arg1)  // node0_field_P
{ 

    if (OPP_DBG) printf("SIMPIC - opp_par_loop__FieldSolveSumLaplace [%d]\n", set->size);

}

//*************************************************************************************************
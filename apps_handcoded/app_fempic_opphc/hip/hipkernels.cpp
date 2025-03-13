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


#include "../fempic_defs.h"
#include "opp_hip.h"

#include "hip/hip_runtime.h"
#include <thrust/execution_policy.h>

__constant__ int OPP_cells_set_size_d;
int OPP_cells_set_size;

__constant__ int OPP_comm_iteration_d;

// TODO : This should be removed - only for testing
double CONST_spwt = 0, CONST_ion_velocity = 0, CONST_dt = 0, CONST_plasma_den = 0, CONST_mass = 0, CONST_charge = 0, CONST_wall_potential = 0;

//****************************************
__constant__ double CONST_spwt_device = 0.0;
__constant__ double CONST_ion_velocity_device = 0.0;
__constant__ double CONST_dt_device = 0.0;
__constant__ double CONST_plasma_den_device = 0.0;
__constant__ double CONST_mass_device = 0.0;
__constant__ double CONST_charge_device = 0.0;
__constant__ double CONST_wall_potential_device = 0.0;

void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    if (!strcmp(name,"CONST_spwt"))              
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_spwt_device), data, dim*size));
    else if (!strcmp(name,"CONST_ion_velocity")) 
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_ion_velocity_device), data, dim*size));
    else if (!strcmp(name,"CONST_dt"))           
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_dt_device), data, dim*size));
    else if (!strcmp(name,"CONST_plasma_den"))   
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_plasma_den_device), data, dim*size));
    else if (!strcmp(name,"CONST_mass"))         
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_mass_device), data, dim*size));
    else if (!strcmp(name,"CONST_charge"))       
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_charge_device), data, dim*size));
    else if (!strcmp(name,"CONST_wall_potential"))
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_wall_potential_device), data, dim*size));
    else 
        std::cerr << "error: unknown const name" << std::endl;

    // TODO : This block should be removed - only for testing
    {   
        if (!strcmp(name,"CONST_spwt"))             
            CONST_spwt = *((double*)data);
        else if (!strcmp(name,"CONST_ion_velocity")) 
            CONST_ion_velocity = *((double*)data);
        else if (!strcmp(name,"CONST_dt"))           
            CONST_dt = *((double*)data);
        else if (!strcmp(name,"CONST_plasma_den"))   
            CONST_plasma_den = *((double*)data);
        else if (!strcmp(name,"CONST_mass"))         
            CONST_mass = *((double*)data);
        else if (!strcmp(name,"CONST_charge"))       
            CONST_charge = *((double*)data);
        else if (!strcmp(name,"CONST_wall_potential")) 
            CONST_wall_potential = *((double*)data);
        else std::cerr << "error: unknown const name" << std::endl; 
       
    } // TODO : This block should be removed
}
//****************************************


//*************************************************************************************************
#include "opp_loop_inject__inject_ions.cpp"

//*************************************************************************************************
#include "opp_loop_all__calc_pos_vel.cpp"

//*************************************************************************************************
#include "opp_particle_move__move.cpp"

//*************************************************************************************************
#include "opp_loop_all__deposit_charge_on_nodes.cpp"

//*************************************************************************************************
#include "opp_loop_all__compute_node_charge_density.cpp"

//*************************************************************************************************
#include "opp_loop_all__compute_electric_field.cpp"

//*************************************************************************************************
#include "opp_loop_all__init_boundary_pot.cpp"

//*************************************************************************************************
#include "opp_loop_all__get_max_values.cpp"
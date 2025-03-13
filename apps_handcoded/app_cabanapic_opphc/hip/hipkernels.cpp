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


#include "../cabana_defs.h"
#include "opp_hip.h"
#include "hip/hip_runtime.h"

__constant__ int OPP_cells_set_size_d;
int OPP_cells_set_size;

__constant__ int OPP_comm_iteration_d;

//****************************************
__constant__ OPP_INT CONST_DEV_c_per_dim[DIM];
__constant__ OPP_REAL CONST_DEV_dt;
__constant__ OPP_REAL CONST_DEV_qsp;
__constant__ OPP_REAL CONST_DEV_cdt_d[DIM];
__constant__ OPP_REAL CONST_DEV_p[DIM];
__constant__ OPP_REAL CONST_DEV_qdt_2mc;
__constant__ OPP_REAL CONST_DEV_dt_eps0;
__constant__ OPP_REAL CONST_DEV_acc_coef[DIM];

void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    if (!strcmp(name,"CONST_c_per_dim"))              
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEV_c_per_dim), data, dim*size));
    else if (!strcmp(name,"CONST_dt")) 
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEV_dt), data, dim*size));
    else if (!strcmp(name,"CONST_qsp"))           
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEV_qsp), data, dim*size));
    else if (!strcmp(name,"CONST_cdt_d"))   
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEV_cdt_d), data, dim*size));
    else if (!strcmp(name,"CONST_p"))         
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEV_p), data, dim*size));
    else if (!strcmp(name,"CONST_qdt_2mc"))       
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEV_qdt_2mc), data, dim*size));
    else if (!strcmp(name,"CONST_dt_eps0"))
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEV_dt_eps0), data, dim*size));
    else if (!strcmp(name,"CONST_acc_coef"))
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEV_acc_coef), data, dim*size));
    else 
        std::cerr << "error: unknown const name" << std::endl;
}
//****************************************

__constant__ OPP_REAL CONST_DEV_fourth         = (1.0 / 4.0);
__constant__ OPP_REAL CONST_DEV_half           = (1.0 / 2.0);
__constant__ OPP_REAL CONST_DEV_one            = 1.0;
__constant__ OPP_REAL CONST_DEV_one_third      = (1.0 / 3.0);
__constant__ OPP_REAL CONST_DEV_two_fifteenths = (2.0 / 15.0);

//*************************************************************************************************
#include "opp_loop_all__interpolate_mesh_fields.cpp"

//*************************************************************************************************
#include "opp_particle_move__move.cpp"

//*************************************************************************************************
#include "opp_loop_all__accumulate_current_to_cells.cpp"

// //*************************************************************************************************
#include "opp_loop_all__half_advance_b.cpp"

//*************************************************************************************************
#include "opp_loop_all__advance_e.cpp"

//*************************************************************************************************
#include "opp_loop_all__get_max_values.cpp"

//*************************************************************************************************
#include "opp_loop_all__compute_energy.cpp"

//*************************************************************************************************
#include "opp_loop_all__update_ghosts_B.cpp"

//*************************************************************************************************
#include "opp_loop_all__update_ghosts.cpp"

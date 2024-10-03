
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

#include "opp_cuda.h"


__constant__ OPP_REAL CONST_spwt_d[1];
__constant__ OPP_REAL CONST_ion_velocity_d[1];
__constant__ OPP_REAL CONST_dt_d[1];
__constant__ OPP_REAL CONST_plasma_den_d[1];
__constant__ OPP_REAL CONST_mass_d[1];
__constant__ OPP_REAL CONST_charge_d[1];
__constant__ OPP_REAL CONST_wall_potential_d[1];
    
__constant__ int OPP_cells_set_size_d;
int OPP_cells_set_size;

__constant__ int OPP_comm_iteration_d;

void opp_decl_const_impl(int dim, int size, char* data, const char* name) {
    
    if (OPP_DBG)
        opp_printf("opp_decl_const_impl", "Registering %s", name);

    if (!strcmp(name, "CONST_spwt")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_spwt_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_ion_velocity")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_ion_velocity_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_dt")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_dt_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_plasma_den")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_plasma_den_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_mass")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_mass_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_charge")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_charge_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_wall_potential")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_wall_potential_d, data, dim * size));
        return;
    }

    opp_abort(std::string("Error: unknown const name") + std::string(name));
}

#include "init_boundary_pot_kernel_loop.hpp"

#include "inject_ions_kernel_loop.hpp"

#include "calculate_new_pos_vel_kernel_loop.hpp"

#include "move_kernel_loop.hpp"

#include "deposit_charge_on_nodes_kernel_loop.hpp"

#include "compute_node_charge_density_kernel_loop.hpp"

#include "compute_electric_field_kernel_loop.hpp"

#include "get_max_cef_kernel_loop.hpp"

#include "get_final_max_values_kernel_loop.hpp"


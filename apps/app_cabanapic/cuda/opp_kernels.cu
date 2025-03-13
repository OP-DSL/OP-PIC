
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
#include "device_kernels/cuda_inline_kernels.h"

OPP_REAL CONST_dt[1];
OPP_REAL CONST_qsp[1];
OPP_REAL CONST_cdt_d[3];
OPP_REAL CONST_p[3];
OPP_REAL CONST_qdt_2mc[1];
OPP_REAL CONST_dt_eps0[1];
OPP_REAL CONST_acc_coef[3];

__constant__ OPP_REAL CONST_dt_d[1];
__constant__ OPP_REAL CONST_qsp_d[1];
__constant__ OPP_REAL CONST_cdt_d_d[3];
__constant__ OPP_REAL CONST_p_d[3];
__constant__ OPP_REAL CONST_qdt_2mc_d[1];
__constant__ OPP_REAL CONST_dt_eps0_d[1];
__constant__ OPP_REAL CONST_acc_coef_d[3];

void opp_decl_const_impl(int dim, int size, char* data, const char* name) {
    
    if (OPP_DBG)
        opp_printf("opp_decl_const_impl", "Registering %s", name);

    if (!strcmp(name, "CONST_dt")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_dt_d, data, dim * size));
        std::memcpy(&CONST_dt, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_qsp")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_qsp_d, data, dim * size));
        std::memcpy(&CONST_qsp, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_cdt_d")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_cdt_d_d, data, dim * size));
        std::memcpy(&CONST_cdt_d, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_p")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_p_d, data, dim * size));
        std::memcpy(&CONST_p, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_qdt_2mc")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_qdt_2mc_d, data, dim * size));
        std::memcpy(&CONST_qdt_2mc, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_dt_eps0")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_dt_eps0_d, data, dim * size));
        std::memcpy(&CONST_dt_eps0, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_acc_coef")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_acc_coef_d, data, dim * size));
        std::memcpy(&CONST_acc_coef, data, (size*dim));
        return;
    }

    opp_abort(std::string("Error: unknown const name") + std::string(name));
}

#include "interpolate_mesh_fields_kernel_loop.hpp"

#include "move_deposit_kernel_loop.hpp"

#include "accumulate_current_to_cells_kernel_loop.hpp"

#include "half_advance_b_kernel_loop.hpp"

#include "update_ghosts_B_kernel_loop.hpp"

#include "update_ghosts_kernel_loop.hpp"

#include "advance_e_kernel_loop.hpp"

#include "compute_energy_kernel_loop.hpp"

#include "get_max_x_values_kernel_loop.hpp"



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


__constant__ OPP_REAL CONST_lhs_voltage_d[1];
__constant__ OPP_REAL CONST_L_d[1];
__constant__ OPP_REAL CONST_xl_d[1];
__constant__ OPP_REAL CONST_xr_d[1];
__constant__ OPP_REAL CONST_dx_d[1];
__constant__ OPP_REAL CONST_qm_d[1];
__constant__ OPP_REAL CONST_qscale_d[1];
__constant__ OPP_INT CONST_neighbour_cell_d[2];
__constant__ OPP_INT CONST_rank_d[1];
__constant__ OPP_INT CONST_comm_size_d[1];
    
__constant__ int OPP_cells_set_size_d;
int OPP_cells_set_size;

__constant__ int OPP_comm_iteration_d;

void opp_decl_const_impl(int dim, int size, char* data, const char* name) {
    
    if (!strcmp(name, "CONST_lhs_voltage")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_lhs_voltage_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_L")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_L_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_xl")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_xl_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_xr")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_xr_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_dx")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_dx_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_qm")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_qm_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_qscale")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_qscale_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_neighbour_cell")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_neighbour_cell_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_rank")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_rank_d, data, dim * size));
        return;
    }
    if (!strcmp(name, "CONST_comm_size")) {
        cutilSafeCall(cudaMemcpyToSymbol(CONST_comm_size_d, data, dim * size));
        return;
    }

    opp_printf("Error: unknown const name %s", name);
    opp_abort("Error: unknown const name");
}

#include "weight_f2p_kernel_loop.hpp"

#include "move_kernel_loop.hpp"

#include "weight_p2f_kernel_loop.hpp"

#include "sum_laplace_kernel_loop.hpp"


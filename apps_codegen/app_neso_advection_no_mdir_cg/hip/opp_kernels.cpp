
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

#include "opp_hip.h"
#include "device_kernels/hip_inline_kernels.h"

OPP_REAL CONST_extents[2];
OPP_REAL CONST_dt[1];
OPP_REAL CONST_cell_width[1];
OPP_INT CONST_ndimcells[2];

__constant__ OPP_REAL CONST_extents_d[2];
__constant__ OPP_REAL CONST_dt_d[1];
__constant__ OPP_REAL CONST_cell_width_d[1];
__constant__ OPP_INT CONST_ndimcells_d[2];

void opp_decl_const_impl(int dim, int size, char* data, const char* name) {
    
    if (OPP_DBG)
        opp_printf("opp_decl_const_impl", "Registering %s", name);

    if (!strcmp(name, "CONST_extents")) {
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_extents_d), data, dim * size));
        std::memcpy(&CONST_extents, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_dt")) {
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_dt_d), data, dim * size));
        std::memcpy(&CONST_dt, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_cell_width")) {
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_cell_width_d), data, dim * size));
        std::memcpy(&CONST_cell_width, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_ndimcells")) {
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_ndimcells_d), data, dim * size));
        std::memcpy(&CONST_ndimcells, data, (size*dim));
        return;
    }

    opp_abort(std::string("Error: unknown const name") + std::string(name));
}

#include "update_pos_kernel_loop.hpp"

#include "move_kernel_loop.hpp"

#include "verify_kernel_loop.hpp"


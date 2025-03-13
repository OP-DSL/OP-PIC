
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

OPP_REAL CONST_lhs_voltage[1];
OPP_REAL CONST_L[1];
OPP_REAL CONST_xl[1];
OPP_REAL CONST_xr[1];
OPP_REAL CONST_dx[1];
OPP_REAL CONST_qm[1];
OPP_REAL CONST_dt[1];
OPP_REAL CONST_qscale[1];
OPP_INT CONST_neighbour_cell[2];
OPP_INT CONST_rank[1];
OPP_INT CONST_comm_size[1];

void opp_decl_const_impl(int dim, int size, char* data, const char* name) {
    
    if (OPP_DBG)
        opp_printf("opp_decl_const_impl", "Registering %s", name);

    if (!strcmp(name, "CONST_lhs_voltage")) {
        std::memcpy(&CONST_lhs_voltage, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_L")) {
        std::memcpy(&CONST_L, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_xl")) {
        std::memcpy(&CONST_xl, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_xr")) {
        std::memcpy(&CONST_xr, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_dx")) {
        std::memcpy(&CONST_dx, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_qm")) {
        std::memcpy(&CONST_qm, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_dt")) {
        std::memcpy(&CONST_dt, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_qscale")) {
        std::memcpy(&CONST_qscale, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_neighbour_cell")) {
        std::memcpy(&CONST_neighbour_cell, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_rank")) {
        std::memcpy(&CONST_rank, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_comm_size")) {
        std::memcpy(&CONST_comm_size, data, (size*dim));
        return;
    }

    opp_abort(std::string("Error: unknown const name") + std::string(name));
}

#include "weight_f2p_kernel_loop.hpp"

#include "move_kernel_loop.hpp"

#include "weight_p2f_kernel_loop.hpp"

#include "sum_laplace_kernel_loop.hpp"



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

#include "opp_sycl.h"

OPP_REAL* CONST_lhs_voltage_s = nullptr;
OPP_REAL* CONST_L_s = nullptr;
OPP_REAL* CONST_xl_s = nullptr;
OPP_REAL* CONST_xr_s = nullptr;
OPP_REAL* CONST_dx_s = nullptr;
OPP_REAL* CONST_qm_s = nullptr;
OPP_REAL* CONST_dt_s = nullptr;
OPP_REAL* CONST_qscale_s = nullptr;
OPP_INT* CONST_neighbour_cell_s = nullptr;
OPP_INT* CONST_rank_s = nullptr;
OPP_INT* CONST_comm_size_s = nullptr;
    
OPP_INT* cells_set_size_s = nullptr;
OPP_INT cells_set_size = -1;

OPP_INT* comm_iteration_s = nullptr;
OPP_INT comm_iteration = -1;

void opp_decl_const_impl(int dim, int size, char* data, const char* name) {
    
    if (OPP_DBG)
        opp_printf("opp_decl_const_impl", "Registering %s", name);

    if (!strcmp(name, "CONST_lhs_voltage")) {
        opp_register_const<OPP_REAL>(CONST_lhs_voltage_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_lhs_voltage_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_L")) {
        opp_register_const<OPP_REAL>(CONST_L_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_L_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_xl")) {
        opp_register_const<OPP_REAL>(CONST_xl_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_xl_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_xr")) {
        opp_register_const<OPP_REAL>(CONST_xr_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_xr_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_dx")) {
        opp_register_const<OPP_REAL>(CONST_dx_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_dx_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_qm")) {
        opp_register_const<OPP_REAL>(CONST_qm_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_qm_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_dt")) {
        opp_register_const<OPP_REAL>(CONST_dt_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_dt_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_qscale")) {
        opp_register_const<OPP_REAL>(CONST_qscale_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_qscale_s, (OPP_REAL*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_neighbour_cell")) {
        opp_register_const<OPP_INT>(CONST_neighbour_cell_s, dim);
        opp_mem::copy_host_to_dev<OPP_INT>(CONST_neighbour_cell_s, (OPP_INT*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_rank")) {
        opp_register_const<OPP_INT>(CONST_rank_s, dim);
        opp_mem::copy_host_to_dev<OPP_INT>(CONST_rank_s, (OPP_INT*)data, dim);
        return;
    }
    if (!strcmp(name, "CONST_comm_size")) {
        opp_register_const<OPP_INT>(CONST_comm_size_s, dim);
        opp_mem::copy_host_to_dev<OPP_INT>(CONST_comm_size_s, (OPP_INT*)data, dim);
        return;
    }

    opp_abort(std::string("Error: unknown const name") + std::string(name));
}

#include "weight_f2p_kernel_loop.hpp"

#include "move_kernel_loop.hpp"

#include "weight_p2f_kernel_loop.hpp"

#include "sum_laplace_kernel_loop.hpp"


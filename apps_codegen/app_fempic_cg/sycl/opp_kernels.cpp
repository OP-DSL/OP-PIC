
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
#include "device_kernels/sycl_inline_kernels.h"

OPP_REAL* CONST_spwt_s = nullptr;
OPP_REAL* CONST_ion_velocity_s = nullptr;
OPP_REAL* CONST_dt_s = nullptr;
OPP_REAL* CONST_plasma_den_s = nullptr;
OPP_REAL* CONST_mass_s = nullptr;
OPP_REAL* CONST_charge_s = nullptr;
OPP_REAL* CONST_wall_potential_s = nullptr;

OPP_REAL CONST_spwt[1];
OPP_REAL CONST_ion_velocity[1];
OPP_REAL CONST_dt[1];
OPP_REAL CONST_plasma_den[1];
OPP_REAL CONST_mass[1];
OPP_REAL CONST_charge[1];
OPP_REAL CONST_wall_potential[1];

void opp_decl_const_impl(int dim, int size, char* data, const char* name) {
    
    if (OPP_DBG)
        opp_printf("opp_decl_const_impl", "Registering %s", name);

    if (!strcmp(name, "CONST_spwt")) {
        opp_register_const<OPP_REAL>(CONST_spwt_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_spwt_s, (OPP_REAL*)data, dim);
        std::memcpy(CONST_spwt, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_ion_velocity")) {
        opp_register_const<OPP_REAL>(CONST_ion_velocity_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_ion_velocity_s, (OPP_REAL*)data, dim);
        std::memcpy(CONST_ion_velocity, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_dt")) {
        opp_register_const<OPP_REAL>(CONST_dt_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_dt_s, (OPP_REAL*)data, dim);
        std::memcpy(CONST_dt, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_plasma_den")) {
        opp_register_const<OPP_REAL>(CONST_plasma_den_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_plasma_den_s, (OPP_REAL*)data, dim);
        std::memcpy(CONST_plasma_den, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_mass")) {
        opp_register_const<OPP_REAL>(CONST_mass_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_mass_s, (OPP_REAL*)data, dim);
        std::memcpy(CONST_mass, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_charge")) {
        opp_register_const<OPP_REAL>(CONST_charge_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_charge_s, (OPP_REAL*)data, dim);
        std::memcpy(CONST_charge, data, (size*dim));
        return;
    }
    if (!strcmp(name, "CONST_wall_potential")) {
        opp_register_const<OPP_REAL>(CONST_wall_potential_s, dim);
        opp_mem::copy_host_to_dev<OPP_REAL>(CONST_wall_potential_s, (OPP_REAL*)data, dim);
        std::memcpy(CONST_wall_potential, data, (size*dim));
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


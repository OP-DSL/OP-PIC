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


#include "../advec_defs.h"
#include "opp_hip.h"

__constant__ int OPP_cells_set_size_d;
int OPP_cells_set_size;

__constant__ int OPP_comm_iteration_d;

//****************************************
// TODO : This should be removed - only for testing
OPP_REAL CONST_extents[2];
OPP_REAL CONST_dt = 0.0;
OPP_REAL CONST_cell_width = 0.0;
OPP_INT CONST_ndimcells[2];

//****************************************
__constant__ OPP_REAL CONST_DEVICE_extents[2];
__constant__ OPP_REAL CONST_DEVICE_dt;
__constant__ OPP_REAL CONST_DEVICE_cell_width;
__constant__ OPP_INT CONST_DEVICE_ndimcells[2];

void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    if (!strcmp(name,"CONST_extents"))              
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEVICE_extents), data, dim*size));
    else if (!strcmp(name,"CONST_dt")) 
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEVICE_dt), data, dim*size));
    else if (!strcmp(name,"CONST_cell_width"))           
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEVICE_cell_width), data, dim*size));
    else if (!strcmp(name,"CONST_ndimcells"))           
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(CONST_DEVICE_ndimcells), data, dim*size));
    else 
        std::cerr << "error: unknown const name" << std::endl;

    // TODO : This block should be removed - only for testing
    {
        if (!strcmp(name,"CONST_extents"))         std::memcpy(CONST_extents, data, (size*dim));
        else if (!strcmp(name,"CONST_dt"))         std::memcpy(&CONST_dt, data, (size*dim));
        else if (!strcmp(name,"CONST_cell_width")) std::memcpy(&CONST_cell_width, data, (size*dim));
        else if (!strcmp(name,"CONST_ndimcells"))  std::memcpy(CONST_ndimcells, data, (size*dim));
        else std::cerr << "error: unknown const name" << std::endl;
    } // TODO : This block should be removed
}
//****************************************

//*************************************************************************************************
#include "opp_loop_all_part__Move.cpp"

//*************************************************************************************************
#include "opp_loop_all__verify.cpp"
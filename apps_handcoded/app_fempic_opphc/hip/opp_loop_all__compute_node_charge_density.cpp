
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

#include "hip/hip_runtime.h"

//user function
//*************************************************************************************************
__device__ void compute_node_charge_density__kernel_gpu(
    double *node_charge_density,
    const double *node_volume
)
{
    (*node_charge_density) *= (CONST_spwt_device / (*node_volume));
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void opp_device_ComputeNodeChargeDensity(
    double *__restrict dir_arg0,
    const double *__restrict dir_arg1,
    int start,
    int end
) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        //user-supplied kernel call
        compute_node_charge_density__kernel_gpu(
            (dir_arg0 + n),
            (dir_arg1 + n)
        );
    }
}

//*************************************************************************************************
void opp_loop_all__compute_node_charge_density(
    opp_set set,     // nodes_set
    opp_arg arg0,    // node_charge_density
    opp_arg arg1     // node_volume
)
{ 
    
    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__compute_node_charge_density set_size %d", set->size);

    opp_profiler->start("ComputeNodeChargeDensity");

    int nargs = 2;
    opp_arg args[nargs];

    args[0] = std::move(arg0);
    args[1] = std::move(arg1);

    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    if (set_size > 0) 
    {
        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = OPP_gpu_threads_per_block;
            int nblocks  = (end - start - 1) / nthread + 1;

            opp_device_ComputeNodeChargeDensity <<<nblocks, nthread>>> (
                (double *)  args[0].data_d,
                (double *)  args[1].data_d,
                start, 
                end);
        } 

        opp_set_dirtybit_grouped(nargs, args, Device_GPU);
        cutilSafeCall(hipDeviceSynchronize());       
    }

    opp_profiler->end("ComputeNodeChargeDensity");
}
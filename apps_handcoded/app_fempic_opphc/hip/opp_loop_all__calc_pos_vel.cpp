
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

int calc_pos_vel_stride_OPP_HOST_0 = -1;
int calc_pos_vel_stride_OPP_HOST_1 = -1;
int calc_pos_vel_stride_OPP_HOST_2 = -1;

__constant__ int calc_pos_vel_stride_OPP_DEVICE_0;
__constant__ int calc_pos_vel_stride_OPP_DEVICE_1;
__constant__ int calc_pos_vel_stride_OPP_DEVICE_2;

//user function
//*************************************************************************************************
__device__ void compute_pos_vel__kernel_gpu(
    const double *cell_ef,
    double *part_pos,
    double *part_vel
)
{
    double coefficient1 = CONST_charge_device / CONST_mass_device * (CONST_dt_device);

    for (int i = 0; i < DIM; i++)
        part_vel[i * calc_pos_vel_stride_OPP_DEVICE_2] += 
            (coefficient1 * cell_ef[i * calc_pos_vel_stride_OPP_DEVICE_0]);           

    for (int i = 0; i < DIM; i++)
        part_pos[i * calc_pos_vel_stride_OPP_DEVICE_1] += 
            part_vel[i * calc_pos_vel_stride_OPP_DEVICE_2] * (CONST_dt_device); // v = u + at
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void opp_device_all_CalcPosVel(
    int *__restrict d_cell_index,
    double *__restrict ind_arg0,
    double *__restrict dir_arg1,
    double *__restrict dir_arg2,
    int start,
    int end) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        int map0idx = d_cell_index[n];

        //user-supplied kernel call
        compute_pos_vel__kernel_gpu(
            (ind_arg0 + map0idx),
            (dir_arg1 + n),
            (dir_arg2 + n)
        );
    }
}

void opp_loop_all__calc_pos_vel(
    opp_set set,      // particles_set
    opp_arg arg0,     // cell_ef,
    opp_arg arg1,     // part_pos,
    opp_arg arg2      // part_vel,
)
{ 
    
    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__calc_pos_vel set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("CalcPosVel");

    int nargs = 3;
    opp_arg args[nargs];

    args[0]  = std::move(arg0);
    args[1]  = std::move(arg1);
    args[2]  = std::move(arg2);

    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    if (set_size > 0) 
    {
        calc_pos_vel_stride_OPP_HOST_0 = args[0].dat->set->set_capacity;
        calc_pos_vel_stride_OPP_HOST_1 = args[1].dat->set->set_capacity;
        calc_pos_vel_stride_OPP_HOST_2 = args[2].dat->set->set_capacity;

        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(calc_pos_vel_stride_OPP_DEVICE_0), 
                                                    &calc_pos_vel_stride_OPP_HOST_0, sizeof(int)));
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(calc_pos_vel_stride_OPP_DEVICE_1), 
                                                    &calc_pos_vel_stride_OPP_HOST_1, sizeof(int)));
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(calc_pos_vel_stride_OPP_DEVICE_2), 
                                                    &calc_pos_vel_stride_OPP_HOST_2, sizeof(int)));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = OPP_gpu_threads_per_block;
            int nblocks = (end - start - 1) / nthread + 1;

            opp_device_all_CalcPosVel <<<nblocks, nthread>>> (
                (int *)     set->mesh_relation_dat->data_d,
                (double *)  args[0].data_d,
                (double *)  args[1].data_d,
                (double *)  args[2].data_d,
                start, 
                end);
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(hipDeviceSynchronize());

    opp_profiler->end("CalcPosVel");
}

//*************************************************************************************************
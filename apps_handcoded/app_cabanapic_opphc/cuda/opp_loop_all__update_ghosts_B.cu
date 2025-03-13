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

int ugb_OPP_HOST_0 = -1;
int ugb_OPP_HOST_1 = -1;
int ugb_OPP_HOST_2_MAP = -1;

__constant__ int ugb_OPP_DEVICE_0;
__constant__ int ugb_OPP_DEVICE_1;
__constant__ int ugb_OPP_DEVICE_2_MAP;

//user function
//*************************************************************************************************
__device__ void update_ghosts_B_kernel_gpu(
    const OPP_INT* c_mask_ugb,
    const OPP_REAL* from_cell,
    OPP_REAL* to_cell,
    const OPP_INT* m_idx)
{
    if (c_mask_ugb[*m_idx * ugb_OPP_DEVICE_0] == 1)
    {
        to_cell[Dim::x * ugb_OPP_DEVICE_1] = from_cell[Dim::x * ugb_OPP_DEVICE_1];
        to_cell[Dim::y * ugb_OPP_DEVICE_1] = from_cell[Dim::y * ugb_OPP_DEVICE_1];
        to_cell[Dim::z * ugb_OPP_DEVICE_1] = from_cell[Dim::z * ugb_OPP_DEVICE_1];
    }  
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void opp_dev_update_ghosts_B(
    const OPP_INT *__restrict dir_arg0,
    const OPP_REAL *__restrict dir_arg1,
    OPP_REAL *__restrict ind_arg2,
    const OPP_INT *__restrict arg2_map,
    const OPP_INT *__restrict dir_arg3,
    const int start,
    const int end
) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        const OPP_INT map0idx = arg2_map[n + ugb_OPP_DEVICE_2_MAP * 0];

        //user-supplied kernel call
        update_ghosts_B_kernel_gpu(
            (dir_arg0 + n),
            (dir_arg1 + n),
            (ind_arg2 + map0idx),
            dir_arg3
        );
    }
}

//*************************************************************************************************
void opp_loop_all__update_ghosts_B(
    opp_set set,     // cells set
    opp_arg arg0,    // cell_mask_ugb,       OPP_READ
    opp_arg arg1,    // cell,                OPP_READ
    opp_arg arg2,    // cell, 0, c2cugb_map, OPP_WRITE
    opp_arg arg3     // mask_idx global
)
{
    
    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__update_ghosts_B set_size %d", set->size);

    opp_profiler->start("UpGhostB");

    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    if (set_size > 0) 
    {
        const int start = 0;
        const int end   = set_size;

        ugb_OPP_HOST_0 = args[0].dat->set->set_capacity;
        ugb_OPP_HOST_1 = args[1].dat->set->set_capacity;
        ugb_OPP_HOST_2_MAP = args[2].size;

        cutilSafeCall(cudaMemcpyToSymbol(ugb_OPP_DEVICE_0, &ugb_OPP_HOST_0, sizeof(int)));
        cutilSafeCall(cudaMemcpyToSymbol(ugb_OPP_DEVICE_1, &ugb_OPP_HOST_1, sizeof(int)));
        cutilSafeCall(cudaMemcpyToSymbol(ugb_OPP_DEVICE_2_MAP, &ugb_OPP_HOST_2_MAP, sizeof(int)));

        cutilSafeCall(cudaMalloc((void**)&(args[3].data_d), sizeof(int)));
        cutilSafeCall(cudaMemcpy((int*)(args[3].data_d), (int*)(args[3].data), sizeof(int), cudaMemcpyHostToDevice));

        int nthread = OPP_gpu_threads_per_block;
        int nblocks = (end - start - 1) / nthread + 1;

        opp_dev_update_ghosts_B<<<nblocks,nthread>>> (
            (OPP_INT *) args[0].data_d,
            (OPP_REAL*) args[1].data_d,
            (OPP_REAL*) args[2].data_d,
            (OPP_INT *) args[2].map_data_d,
            (OPP_INT *) args[3].data_d,
            start, 
            end
        );

        cutilSafeCall(cudaFree(args[3].data_d));

        opp_set_dirtybit_grouped(nargs, args, Device_GPU);
        cutilSafeCall(cudaDeviceSynchronize());       
    }

    opp_profiler->end("UpGhostB");
}
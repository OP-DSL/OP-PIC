
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

__constant__ int hab_OPP_DEVICE_0;
__constant__ int hab_OPP_DEVICE_4;
__constant__ int hab_OPP_DEVICE_0_MAP;

int hab_OPP_HOST_0 = -1;
int hab_OPP_HOST_4 = -1;
int hab_OPP_HOST_0_MAP = -1;

//user function
//*************************************************************************************************
__device__ void dev_half_advance_b__kernel (
    const OPP_REAL* cell_x_e, 
    const OPP_REAL* cell_y_e, 
    const OPP_REAL* cell_z_e, 
    const OPP_REAL* cell0_e, 
    OPP_REAL* cell0_b,
    const OPP_INT* cell0_ghost)
{
    if (cell0_ghost[0] == 0) 
    {
        // No need of atomics here, since we are directly incrementing core elements 
        
        cell0_b[hab_OPP_DEVICE_4 * Dim::x] -= ( 0.5 * CONST_DEV_p[Dim::y] * 
                        ( cell_y_e[hab_OPP_DEVICE_0 * Dim::z] - cell0_e[hab_OPP_DEVICE_0 * Dim::z] ) 
                            - 0.5 * CONST_DEV_p[Dim::z] * 
                        ( cell_z_e[hab_OPP_DEVICE_0 * Dim::y] - cell0_e[hab_OPP_DEVICE_0 * Dim::y] ) );

        cell0_b[hab_OPP_DEVICE_4 * Dim::y] -= ( 0.5 * CONST_DEV_p[Dim::z] * 
                        ( cell_z_e[hab_OPP_DEVICE_0 * Dim::x] - cell0_e[hab_OPP_DEVICE_0 * Dim::x] ) 
                            - 0.5 * CONST_DEV_p[Dim::x] * 
                        ( cell_x_e[hab_OPP_DEVICE_0 * Dim::z] - cell0_e[hab_OPP_DEVICE_0 * Dim::z] ) );

        cell0_b[hab_OPP_DEVICE_4 * Dim::z] -= ( 0.5 * CONST_DEV_p[Dim::x] * 
                        ( cell_x_e[hab_OPP_DEVICE_0 * Dim::y] - cell0_e[hab_OPP_DEVICE_0 * Dim::y] ) 
                            - 0.5 * CONST_DEV_p[Dim::y] * 
                        ( cell_y_e[hab_OPP_DEVICE_0 * Dim::x] - cell0_e[hab_OPP_DEVICE_0 * Dim::x] ) );
    }
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void dev_half_advance_b(
    const OPP_REAL *__restrict arg0_ind,
    const OPP_INT *__restrict arg0_map,
    const OPP_REAL *__restrict arg1_ind,
    const OPP_REAL *__restrict arg2_ind,
    const OPP_REAL *__restrict arg3_dir,
    OPP_REAL *__restrict arg4_dir,
    const OPP_INT *__restrict arg5_dir,
    int start,
    int end) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        const OPP_INT map0idx = arg0_map[n + hab_OPP_DEVICE_0_MAP * CellMap::xu_y_z];
        const OPP_INT map1idx = arg0_map[n + hab_OPP_DEVICE_0_MAP * CellMap::x_yu_z];
        const OPP_INT map2idx = arg0_map[n + hab_OPP_DEVICE_0_MAP * CellMap::x_y_zu];

        //user-supplied kernel call
        dev_half_advance_b__kernel(
            (arg0_ind + map0idx),
            (arg1_ind + map1idx),
            (arg2_ind + map2idx),
            (arg3_dir + n),
            (arg4_dir + n),
            (arg5_dir + n)
        );
    }
}

//*************************************************************************************************
void opp_loop_all__half_advance_b(
    opp_set set,     // cells set
    opp_arg arg0,    // cell_x_e        // OPP_READ
    opp_arg arg1,    // cell_y_e        // OPP_READ
    opp_arg arg2,    // cell_z_e        // OPP_READ
    opp_arg arg3,    // cell0_e         // OPP_READ
    opp_arg arg4,    // cell0_b         // OPP_INC
    opp_arg arg5     // cell0_ghost     // OPP_READ
)
{

    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__half_advance_b set_size %d", set->size);

    opp_profiler->start("HalfAdv_B");

    const int nargs = 6;
    opp_arg args[nargs];

    args[0] = std::move(arg0);
    args[1] = std::move(arg1);
    args[2] = std::move(arg2);
    args[3] = std::move(arg3);
    args[4] = std::move(arg4);
    args[5] = std::move(arg5);

    opp_profiler->start("HAdv_B_Halo");
    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    opp_profiler->end("HAdv_B_Halo");
    
    if (set_size > 0) 
    {
        hab_OPP_HOST_0 = args[0].dat->set->set_capacity;
        hab_OPP_HOST_4 = args[4].dat->set->set_capacity;
        hab_OPP_HOST_0_MAP = args[0].size;

        cutilSafeCall(cudaMemcpyToSymbol(hab_OPP_DEVICE_0, &hab_OPP_HOST_0, sizeof(int)));
        cutilSafeCall(cudaMemcpyToSymbol(hab_OPP_DEVICE_4, &hab_OPP_HOST_4, sizeof(int)));
        cutilSafeCall(cudaMemcpyToSymbol(hab_OPP_DEVICE_0_MAP, &hab_OPP_HOST_0_MAP, sizeof(int)));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = OPP_gpu_threads_per_block;
            int nblocks = (end - start - 1) / nthread + 1;

            dev_half_advance_b <<<nblocks, nthread>>> (
                (OPP_REAL*) args[0].data_d,         // cell_x_e        // OPP_WRITE
                (OPP_INT *) args[0].map_data_d,
                (OPP_REAL*) args[1].data_d,         // cell_y_e        // OPP_READ
                (OPP_REAL*) args[2].data_d,         // cell_z_e        // OPP_READ
                (OPP_REAL*) args[3].data_d,         // cell0_e         // OPP_READ
                (OPP_REAL*) args[4].data_d,         // cell0_b         // OPP_INC
                (OPP_INT*)  args[5].data_d,         // cell0_ghost     // OPP_READ
                start, 
                end);
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());

    opp_profiler->end("HalfAdv_B");   
}

//*************************************************************************************************

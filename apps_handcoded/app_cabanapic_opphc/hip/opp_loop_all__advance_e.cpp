
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

__constant__ int ae_OPP_DEVICE_0;
__constant__ int ae_OPP_DEVICE_4;
__constant__ int ae_OPP_DEVICE_5;
__constant__ int ae_OPP_DEVICE_0_MAP;

int ae_OPP_HOST_0 = -1;
int ae_OPP_HOST_4 = -1;
int ae_OPP_HOST_5 = -1;
int ae_OPP_HOST_0_MAP = -1;

//user function
//*************************************************************************************************
__device__ void dev_advance_e__kernel(
    const OPP_REAL* cell_x_b, 
    const OPP_REAL* cell_y_b, 
    const OPP_REAL* cell_z_b, 
    const OPP_REAL* cell0_b, 
    const OPP_REAL* cell0_j, 
    OPP_REAL* cell0_e,
    const OPP_INT* iter_adv_e
)
{
    if (iter_adv_e[0] == 1)
    {
        // No need of atomics here, since we are directly incrementing core elements 

        cell0_e[ae_OPP_DEVICE_5 * Dim::x] += ( - CONST_DEV_dt_eps0 * cell0_j[ae_OPP_DEVICE_4 * Dim::x] ) + 
            ( CONST_DEV_p[Dim::y] * (cell0_b[ae_OPP_DEVICE_0 * Dim::z] - cell_y_b[ae_OPP_DEVICE_0 * Dim::z]) - 
            CONST_DEV_p[Dim::z] * (cell0_b[ae_OPP_DEVICE_0 * Dim::y] - cell_z_b[ae_OPP_DEVICE_0 * Dim::y]) );

        cell0_e[ae_OPP_DEVICE_5 * Dim::y] += ( - CONST_DEV_dt_eps0 * cell0_j[ae_OPP_DEVICE_4 * Dim::y] ) +            
            ( CONST_DEV_p[Dim::z] * (cell0_b[ae_OPP_DEVICE_0 * Dim::x] - cell_z_b[ae_OPP_DEVICE_0 * Dim::x]) - 
            CONST_DEV_p[Dim::x] * (cell0_b[ae_OPP_DEVICE_0 * Dim::z] - cell_x_b[ae_OPP_DEVICE_0 * Dim::z]) );

        cell0_e[ae_OPP_DEVICE_5 * Dim::z] += ( - CONST_DEV_dt_eps0 * cell0_j[ae_OPP_DEVICE_4 * Dim::z] ) +           
            ( CONST_DEV_p[Dim::x] * (cell0_b[ae_OPP_DEVICE_0 * Dim::y] - cell_x_b[ae_OPP_DEVICE_0 * Dim::y]) - 
            CONST_DEV_p[Dim::y] * (cell0_b[ae_OPP_DEVICE_0 * Dim::x] - cell_y_b[ae_OPP_DEVICE_0 * Dim::x]) );  
    }
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void dev_advance_e(
    const OPP_REAL *__restrict arg0_ind,
    const OPP_INT *__restrict arg0_map,
    const OPP_REAL *__restrict arg1_ind,
    const OPP_REAL *__restrict arg2_ind,
    const OPP_REAL *__restrict arg3_dir,
    const OPP_REAL *__restrict arg4_dir,
    OPP_REAL *__restrict arg5_dir,
    const OPP_INT *__restrict arg6_dir,
    int start,
    int end) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        const OPP_INT map0idx = arg0_map[n + ae_OPP_DEVICE_0_MAP * CellMap::xd_y_z];
        const OPP_INT map1idx = arg0_map[n + ae_OPP_DEVICE_0_MAP * CellMap::x_yd_z];
        const OPP_INT map2idx = arg0_map[n + ae_OPP_DEVICE_0_MAP * CellMap::x_y_zd];

        //user-supplied kernel call
        dev_advance_e__kernel(
            (arg0_ind + map0idx),
            (arg1_ind + map1idx),
            (arg2_ind + map2idx),
            (arg3_dir + n),
            (arg4_dir + n),
            (arg5_dir + n),
            (arg6_dir + n)
        );
    }
}

//*************************************************************************************************
void opp_loop_all__advance_e(
    opp_set set,     // cells set
    opp_arg arg0,    // cell_x_b        // OPP_READ
    opp_arg arg1,    // cell_y_b        // OPP_READ
    opp_arg arg2,    // cell_z_b        // OPP_READ
    opp_arg arg3,    // cell0_b         // OPP_READ
    opp_arg arg4,    // cell0_j         // OPP_READ
    opp_arg arg5,    // cell0_e         // OPP_INC
    opp_arg arg6     // iter_adv_e      // OPP_READ
)
{

    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__advance_e set_size %d", set->size);

    opp_profiler->start("Adv_E");

    const int nargs = 7;
    opp_arg args[nargs];

    args[0] = std::move(arg0);
    args[1] = std::move(arg1);
    args[2] = std::move(arg2);
    args[3] = std::move(arg3);
    args[4] = std::move(arg4);
    args[5] = std::move(arg5);
    args[6] = std::move(arg6);

    opp_profiler->start("Adv_E_Halo");
    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    opp_profiler->end("Adv_E_Halo");

    if (set_size > 0) 
    {
        ae_OPP_HOST_0 = args[0].dat->set->set_capacity;
        ae_OPP_HOST_4 = args[4].dat->set->set_capacity;
        ae_OPP_HOST_5 = args[5].dat->set->set_capacity;
        ae_OPP_HOST_0_MAP = args[0].size;

        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(ae_OPP_DEVICE_0), 
                                                    &ae_OPP_HOST_0, sizeof(int)));
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(ae_OPP_DEVICE_4), 
                                                    &ae_OPP_HOST_4, sizeof(int)));
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(ae_OPP_DEVICE_5), 
                                                    &ae_OPP_HOST_5, sizeof(int)));
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(ae_OPP_DEVICE_0_MAP), 
                                                    &ae_OPP_HOST_0_MAP, sizeof(int)));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = OPP_gpu_threads_per_block;
            int nblocks = (end - start - 1) / nthread + 1;

            dev_advance_e <<<nblocks, nthread>>> (
                (OPP_REAL*) args[0].data_d,         // cell_x_b        // OPP_READ
                (OPP_INT *) args[0].map_data_d,
                (OPP_REAL*) args[1].data_d,         // cell_y_b        // OPP_READ
                (OPP_REAL*) args[2].data_d,         // cell_z_b        // OPP_READ
                (OPP_REAL*) args[3].data_d,         // cell0_b         // OPP_READ
                (OPP_REAL*) args[4].data_d,         // cell0_j         // OPP_READ
                (OPP_REAL*) args[5].data_d,         // cell0_e         // OPP_INC
                (OPP_INT*)  args[6].data_d,         // iter_adv_e      // OPP_READ
                start, 
                end);
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(hipDeviceSynchronize());

    opp_profiler->end("Adv_E");   
}

//*************************************************************************************************

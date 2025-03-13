
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

int acc_OPP_HOST_0 = -1;
int acc_OPP_HOST_2 = -1;
int acc_OPP_HOST_2_MAP_STRIDE = -1;

__constant__ int acc_OPP_DEV_0;
__constant__ int acc_OPP_DEV_2;
__constant__ int acc_OPP_DEV_2_MAP_STRIDE;

//user function
//*************************************************************************************************
__device__ void dev_accumulate_current_to_cells__kernel(
    OPP_REAL* cell0_j, 
    const OPP_REAL* cell0_acc,
    const OPP_REAL* cell_xd_acc, 
    const OPP_REAL* cell_yd_acc, 
    const OPP_REAL* cell_zd_acc, 
    const OPP_REAL* cell_xyd_acc, 
    const OPP_REAL* cell_yzd_acc, 
    const OPP_REAL* cell_xzd_acc,
    const OPP_INT* iter_acc
)
{
    if (iter_acc[0] == 1)
    {
        cell0_j[acc_OPP_DEV_0 * Dim::x] = CONST_DEV_acc_coef[Dim::x] * 
                                        (cell0_acc[acc_OPP_DEV_2 * (CellAcc::jfx + 0)] +
                                        cell_yd_acc[acc_OPP_DEV_2 * (CellAcc::jfx + 1)] +
                                        cell_zd_acc[acc_OPP_DEV_2 * (CellAcc::jfx + 2)] +
                                        cell_yzd_acc[acc_OPP_DEV_2 * (CellAcc::jfx + 3)]);

        cell0_j[acc_OPP_DEV_0 * Dim::y] = CONST_DEV_acc_coef[Dim::y] * 
                                        (cell0_acc[acc_OPP_DEV_2 * (CellAcc::jfy + 0)] +
                                        cell_zd_acc[acc_OPP_DEV_2 * (CellAcc::jfy + 1)] +
                                        cell_xd_acc[acc_OPP_DEV_2 * (CellAcc::jfy + 2)] +
                                        cell_xzd_acc[acc_OPP_DEV_2 * (CellAcc::jfy + 3)]);

        cell0_j[acc_OPP_DEV_0 * Dim::z] = CONST_DEV_acc_coef[Dim::z] * 
                                        (cell0_acc[acc_OPP_DEV_2 * (CellAcc::jfz + 0)] +
                                        cell_xd_acc[acc_OPP_DEV_2 * (CellAcc::jfz + 1)] +
                                        cell_yd_acc[acc_OPP_DEV_2 * (CellAcc::jfz + 2)] +
                                        cell_xyd_acc[acc_OPP_DEV_2 * (CellAcc::jfz + 3)]);
    }
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void dev_accumulate_current_to_cells(
    OPP_REAL *__restrict__ arg0_dir,
    const OPP_REAL *__restrict__ arg1_dir,
    const OPP_REAL *__restrict__ arg2_ind,
    const OPP_INT  *__restrict__ arg2_map,     
    const OPP_REAL *__restrict__ arg3_ind,
    const OPP_REAL *__restrict__ arg4_ind,
    const OPP_REAL *__restrict__ arg5_ind,
    const OPP_REAL *__restrict__ arg6_ind,
    const OPP_REAL *__restrict__ arg7_ind,
    const OPP_INT *__restrict__ arg8_dir,
    int start,
    int end
) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        const OPP_INT map_2idx  = arg2_map[n + acc_OPP_DEV_2_MAP_STRIDE * CellMap::xd_y_z];
        const OPP_INT map_3idx  = arg2_map[n + acc_OPP_DEV_2_MAP_STRIDE * CellMap::x_yd_z];
        const OPP_INT map_4idx  = arg2_map[n + acc_OPP_DEV_2_MAP_STRIDE * CellMap::x_y_zd];
        const OPP_INT map_5idx  = arg2_map[n + acc_OPP_DEV_2_MAP_STRIDE * CellMap::xd_yd_z];
        const OPP_INT map_6idx  = arg2_map[n + acc_OPP_DEV_2_MAP_STRIDE * CellMap::x_yd_zd];
        const OPP_INT map_7idx  = arg2_map[n + acc_OPP_DEV_2_MAP_STRIDE * CellMap::xd_y_zd];

        //user-supplied kernel call
        dev_accumulate_current_to_cells__kernel(
            (arg0_dir + n),
            (arg1_dir + n),
            (arg2_ind + map_2idx),
            (arg3_ind + map_3idx),
            (arg4_ind + map_4idx),
            (arg5_ind + map_5idx),
            (arg6_ind + map_6idx),
            (arg7_ind + map_7idx),
            (arg8_dir + n)
        );
    }
}

//*************************************************************************************************
void opp_loop_all__accumulate_current_to_cells(
    opp_set set,     // cells set
    opp_arg arg0,    // cell0_j         // OPP_WRITE
    opp_arg arg1,    // cell0_acc       // OPP_READ
    opp_arg arg2,    // cell_xd_acc     // OPP_READ
    opp_arg arg3,    // cell_yd_acc     // OPP_READ
    opp_arg arg4,    // cell_zd_acc     // OPP_READ
    opp_arg arg5,    // cell_xyd_acc    // OPP_READ
    opp_arg arg6,    // cell_yzd_acc    // OPP_READ
    opp_arg arg7,    // cell_xzd_acc    // OPP_READ
    opp_arg arg8     // iter_acc        // OPP_READ
)
{ 
    
    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__accumulate_current_to_cells set_size %d", set->size);

    opp_profiler->start("Acc_Current");

    const int nargs = 9;
    opp_arg args[nargs];

    args[0] = std::move(arg0);
    args[1] = std::move(arg1);
    args[2] = std::move(arg2);
    args[3] = std::move(arg3);
    args[4] = std::move(arg4);
    args[5] = std::move(arg5);
    args[6] = std::move(arg6);
    args[7] = std::move(arg7);
    args[8] = std::move(arg8);

    opp_profiler->start("Acc_Halo");
    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    opp_profiler->end("Acc_Halo");

    if (set_size > 0) 
    {
        acc_OPP_HOST_0 = args[0].dat->set->set_capacity;
        acc_OPP_HOST_2 = args[2].dat->set->set_capacity;
        acc_OPP_HOST_2_MAP_STRIDE = args[2].size;

        cutilSafeCall(cudaMemcpyToSymbol(acc_OPP_DEV_0, 
                                                    &acc_OPP_HOST_0, sizeof(int)));
        cutilSafeCall(cudaMemcpyToSymbol(acc_OPP_DEV_2, 
                                                    &acc_OPP_HOST_2, sizeof(int)));
        cutilSafeCall(cudaMemcpyToSymbol(acc_OPP_DEV_2_MAP_STRIDE, 
                                                    &acc_OPP_HOST_2_MAP_STRIDE, sizeof(int)));
        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = OPP_gpu_threads_per_block;
            int nblocks  = (end - start - 1) / nthread + 1;

            dev_accumulate_current_to_cells <<<nblocks, nthread>>> (
                (OPP_REAL*) args[0].data_d,         // cell0_j         // OPP_WRITE
                (OPP_REAL*) args[1].data_d,         // cell0_acc       // OPP_READ
                (OPP_REAL*) args[2].data_d,         // cell_xd_acc     // OPP_READ
                (OPP_INT *) args[2].map_data_d,     
                (OPP_REAL*) args[3].data_d,         // cell_yd_acc     // OPP_READ
                (OPP_REAL*) args[4].data_d,         // cell_zd_acc     // OPP_READ
                (OPP_REAL*) args[5].data_d,         // cell_xyd_acc    // OPP_READ
                (OPP_REAL*) args[6].data_d,         // cell_yzd_acc    // OPP_READ
                (OPP_REAL*) args[7].data_d,         // cell_xzd_acc    // OPP_READ
                (OPP_INT*)  args[8].data_d,         // iter_acc        // OPP_READ
                start, 
                end);
        } 

        opp_set_dirtybit_grouped(nargs, args, Device_GPU);
        cutilSafeCall(cudaDeviceSynchronize());   
    }

    opp_profiler->end("Acc_Current");
}
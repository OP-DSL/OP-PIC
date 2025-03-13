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


int dc_stride_OPP_HOST_0 = -1;
int dc_stride_OPP_HOST_1 = -1;
int dc_stride_OPP_HOST_X = -1;

__constant__ int dc_stride_OPP_DEV0;
__constant__ int dc_stride_OPP_DEV1;
__constant__ int dc_stride_OPP_DEV_X;

// Segmented reduction functions START ------------------------------------------------------------

thrust::device_vector<int> sg_keys_dv;
thrust::device_vector<double> sg_values_dv;

//*************************************************************************************************
__global__ void create_key_value_pairs(
    const int *__restrict d_cell_index,
    int *__restrict keys,
    double *__restrict values,
    const double *__restrict dir_arg0,
    const int *__restrict oppDat1Map,
    const int start,
    const int end) 
{
    const int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        const int n = tid + start;

        const int map0idx = d_cell_index[n];
        
        keys[n + dc_stride_OPP_DEV_X * 0] = oppDat1Map[map0idx + dc_stride_OPP_DEV1 * 0];
        keys[n + dc_stride_OPP_DEV_X * 1] = oppDat1Map[map0idx + dc_stride_OPP_DEV1 * 1];
        keys[n + dc_stride_OPP_DEV_X * 2] = oppDat1Map[map0idx + dc_stride_OPP_DEV1 * 2];
        keys[n + dc_stride_OPP_DEV_X * 3] = oppDat1Map[map0idx + dc_stride_OPP_DEV1 * 3];

        values[n + dc_stride_OPP_DEV_X * 0] = dir_arg0[n + dc_stride_OPP_DEV0 * 0];
        values[n + dc_stride_OPP_DEV_X * 1] = dir_arg0[n + dc_stride_OPP_DEV0 * 1];
        values[n + dc_stride_OPP_DEV_X * 2] = dir_arg0[n + dc_stride_OPP_DEV0 * 2];
        values[n + dc_stride_OPP_DEV_X * 3] = dir_arg0[n + dc_stride_OPP_DEV0 * 3];
    }
}

//*************************************************************************************************
__global__ void assign_values(
    const int *__restrict keys,
    const double *__restrict values,
    double *__restrict arg1,
    const int start,
    const int end) 
{
    const int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        const int n = tid + start;

        const int map0idx = keys[n];
        
        arg1[map0idx] = values[n];
    }
}

// Segmented reduction functions END --------------------------------------------------------------

//user function
//*************************************************************************************************
__device__ void dep_node_charge__kernel_gpu(
    const double* part_lc,
    double *node_charge_den0,
    double *node_charge_den1,
    double *node_charge_den2,
    double *node_charge_den3
)
{
    atomicAdd(node_charge_den0, (part_lc[0 * dc_stride_OPP_DEV0]));
    atomicAdd(node_charge_den1, (part_lc[1 * dc_stride_OPP_DEV0]));
    atomicAdd(node_charge_den2, (part_lc[2 * dc_stride_OPP_DEV0]));
    atomicAdd(node_charge_den3, (part_lc[3 * dc_stride_OPP_DEV0]));
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void opp_dev_all_DepNodeCharge(
    int *__restrict d_cell_index,
    double *__restrict dir_arg0,
    double *__restrict ind_arg1,
    const int *__restrict oppDat1Map,
    double *__restrict ind_arg2,
    double *__restrict ind_arg3,
    double *__restrict ind_arg4,
    int start,
    int end) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    double* k_arg1 = ind_arg1;
    
    if (tid + start < end) 
    {
        int n = tid + start;

        int map0idx = d_cell_index[n];
        
        const int map1idx = oppDat1Map[map0idx + dc_stride_OPP_DEV1 * 0];
        const int map2idx = oppDat1Map[map0idx + dc_stride_OPP_DEV1 * 1];
        const int map3idx = oppDat1Map[map0idx + dc_stride_OPP_DEV1 * 2];
        const int map4idx = oppDat1Map[map0idx + dc_stride_OPP_DEV1 * 3];

        //user-supplied kernel call
        dep_node_charge__kernel_gpu(
            (dir_arg0 + n),
            (k_arg1 + map1idx),
            (k_arg1 + map2idx),
            (k_arg1 + map3idx),
            (k_arg1 + map4idx)
        );
    }
}

void opp_loop_all__deposit_charge_on_nodes(
    opp_set set,      // particles_set
    opp_arg arg0,     // part_lc,
    opp_arg arg1,     // node_charge_den0,
    opp_arg arg2,     // node_charge_den1,
    opp_arg arg3,     // node_charge_den2,
    opp_arg arg4      // node_charge_den3,
)
{ 
    
    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__deposit_charge_on_nodes set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("DepCharge");

    int nargs = 5;
    opp_arg args[nargs];

    args[0]  = std::move(arg0);
    args[1]  = std::move(arg1);
    args[2]  = std::move(arg2);
    args[3]  = std::move(arg3);
    args[4]  = std::move(arg4);

    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);

#ifdef USE_MPI
    opp_init_double_indirect_reductions_cuda(nargs, args);
#endif

    if (set_size > 0) 
    {
        dc_stride_OPP_HOST_0 = args[0].dat->set->set_capacity;
        dc_stride_OPP_HOST_1 = args[1].size;

        cudaMemcpyToSymbol(dc_stride_OPP_DEV0, &dc_stride_OPP_HOST_0, sizeof(int));
        cudaMemcpyToSymbol(dc_stride_OPP_DEV1, &dc_stride_OPP_HOST_1, sizeof(int));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = OPP_gpu_threads_per_block;
            int nblocks = (end - start - 1) / nthread + 1;

            if (!opp_params->get<OPP_BOOL>("use_reg_red")) // Do atomics ----------
            {
                cutilSafeCall(cudaDeviceSynchronize());
                opp_profiler->start("Dep_kernel");
                opp_dev_all_DepNodeCharge <<<nblocks, nthread>>> (
                    (int *)     set->mesh_relation_dat->data_d,
                    (double *)  args[0].data_d,
                    (double *)  args[1].data_d,
                    (int *)     args[1].map_data_d,
                    (double *)  args[2].data_d,
                    (double *)  args[3].data_d,
                    (double *)  args[4].data_d,
                    start, 
                    end);
                
                cutilSafeCall(cudaDeviceSynchronize());
                opp_profiler->end("Dep_kernel");
            }
            else // Do segmented reductions ----------
            {
                dc_stride_OPP_HOST_X = args[0].dat->set->size;
                cudaMemcpyToSymbol(dc_stride_OPP_DEV_X, &dc_stride_OPP_HOST_X, sizeof(int));

                opp_profiler->start("Dep_Resize1");
                const size_t operating_size = (size_t)(args[0].dat->dim * args[0].dat->set->size); //  no halo in args[0]
                const size_t resize_size = (size_t)(args[0].dat->dim * args[0].dat->set->set_capacity); 
                if (resize_size != sg_keys_dv.size()) // resize only if current vector is small
                {
                    sg_keys_dv.resize(resize_size);
                    sg_values_dv.resize(resize_size);
                }
                opp_profiler->end("Dep_Resize1");

                // Create key/value pairs by nodes (linked with particle->cell->node) with part_lc values
                opp_profiler->start("Dep_CrKeyVal");
                create_key_value_pairs<<<nblocks, nthread>>> (
                    (int *)     set->mesh_relation_dat->data_d,
                    (int *)     thrust::raw_pointer_cast(sg_keys_dv.data()),
                    (double *)  thrust::raw_pointer_cast(sg_values_dv.data()),
                    (double *)  args[0].data_d,
                    (int *)     args[1].map_data_d,
                    start, 
                    end);
                cutilSafeCall(cudaDeviceSynchronize());
                opp_profiler->end("Dep_CrKeyVal");

                // Sort by keys to bring the identical keys together
                opp_profiler->start("Dep_Sort");
                thrust::sort_by_key(sg_keys_dv.begin(), sg_keys_dv.begin() + operating_size, sg_values_dv.begin());
                opp_profiler->end("Dep_Sort");

                // Compute the unique keys and their corresponding values
                opp_profiler->start("Dep_Red");
                auto new_end = thrust::reduce_by_key(
                    sg_keys_dv.begin(), sg_keys_dv.begin() + operating_size,
                    sg_values_dv.begin(),
                    sg_keys_dv.begin(),
                    sg_values_dv.begin());        
                opp_profiler->end("Dep_Red");

                const size_t expected_size = (size_t)(args[1].dat->set->size + args[1].dat->set->exec_size + 
                                            args[1].dat->set->nonexec_size) * args[1].dat->dim;

                const size_t reduced_size = (new_end.first - sg_keys_dv.begin());
                if (reduced_size != expected_size) // Might not be equal if all cells dont have particles
                {
                    start = 0;
                    end = (int)reduced_size;

                    // Assign reduced values to the nodes using keys/values
                    opp_profiler->start("Dep_Assign");                
                    assign_values<<<nblocks, nthread>>> (
                        (int *)     thrust::raw_pointer_cast(sg_keys_dv.data()),
                        (double *)  thrust::raw_pointer_cast(sg_values_dv.data()),
                        (double *)  args[1].data_d,
                        start, 
                        end);
                    cutilSafeCall(cudaDeviceSynchronize());
                    opp_profiler->end("Dep_Assign");
                }
                else // all nodes are mapped during reduction
                { 
                    opp_profiler->start("Dep_AssignAll"); 
                    *(args[1].dat->thrust_real) = sg_values_dv; 
                    args[1].dat->data_d = (char*)thrust::raw_pointer_cast(args[1].dat->thrust_real->data());
                    opp_profiler->end("Dep_AssignAll");
                } 

                if (opp_params->get<OPP_INT>("max_iter") == (OPP_main_loop_iter + 1))
                {
                    cutilSafeCall(cudaDeviceSynchronize());
                    sg_keys_dv.clear(); sg_keys_dv.shrink_to_fit();
                    sg_values_dv.clear(); sg_values_dv.shrink_to_fit();
                }          
            }
        }
    }

    cutilSafeCall(cudaDeviceSynchronize());
    
    opp_profiler->start("DI_red");
#ifdef USE_MPI
    opp_exchange_double_indirect_reductions_cuda(nargs, args);
    opp_complete_double_indirect_reductions_cuda(nargs, args);
#endif
    opp_profiler->end("DI_red");

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());

    opp_profiler->end("DepCharge");
}

//*************************************************************************************************
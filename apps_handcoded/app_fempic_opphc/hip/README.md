# DEPOSIT CHARGE WITH BLOCK SEGMENTED REDUCTION BACKUP

//*********************************************
// AUTO GENERATED CODE
//*********************************************

#include "hip/hip_runtime.h"
#include <thrust/execution_policy.h>

int dep_charge_stride_OPP_HOST_0 = -1;
int dep_charge_stride_OPP_HOST_1 = -1;
int dep_charge_stride_OPP_HOST_X = -1;

__constant__ int dep_charge_stride_OPP_DEVICE_0;
__constant__ int dep_charge_stride_OPP_DEVICE_1;
__constant__ int dep_charge_stride_OPP_DEVICE_X;

#ifndef USE_REDUCE_BY_KEY
    #define USE_REDUCE_BY_KEY
#endif

#ifndef USE_REDUCE_BY_KEY
/* // Uncomment For Shared memory 
int dep_charge_stride_OPP_HOST_1_SH = -1;
__constant__ int dep_charge_stride_OPP_DEVICE_1_SH;
__constant__ bool use_shared_OPP_DEVICE = false;
__constant__ int gpu_threads_per_block_OPP_DEVICE;
*/

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
    atomicAdd(node_charge_den0, (part_lc[0 * dep_charge_stride_OPP_DEVICE_0]));
    atomicAdd(node_charge_den1, (part_lc[1 * dep_charge_stride_OPP_DEVICE_0]));
    atomicAdd(node_charge_den2, (part_lc[2 * dep_charge_stride_OPP_DEVICE_0]));
    atomicAdd(node_charge_den3, (part_lc[3 * dep_charge_stride_OPP_DEVICE_0]));

    // unsafeAtomicAdd(node_charge_den0, (part_lc[0 * dep_charge_stride_OPP_DEVICE_0]));
    // unsafeAtomicAdd(node_charge_den1, (part_lc[1 * dep_charge_stride_OPP_DEVICE_0]));
    // unsafeAtomicAdd(node_charge_den2, (part_lc[2 * dep_charge_stride_OPP_DEVICE_0]));
    // unsafeAtomicAdd(node_charge_den3, (part_lc[3 * dep_charge_stride_OPP_DEVICE_0]));

    // *node_charge_den0 += (part_lc[0 * dep_charge_stride_OPP_DEVICE_0]);
    // *node_charge_den1 += (part_lc[1 * dep_charge_stride_OPP_DEVICE_0]);
    // *node_charge_den2 += (part_lc[2 * dep_charge_stride_OPP_DEVICE_0]);
    // *node_charge_den3 += (part_lc[3 * dep_charge_stride_OPP_DEVICE_0]);
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void opp_device_all_DepNodeCharge(
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
    
    /* // Uncomment For Shared memory 
    extern __shared__ double shared_arg1[];
    if (use_shared_OPP_DEVICE)
    {
        for (int i = threadIdx.x; i < dep_charge_stride_OPP_DEVICE_1_SH; i+=gpu_threads_per_block_OPP_DEVICE) 
            shared_arg1[i] = 0;
        __syncthreads();

        k_arg1 = shared_arg1;
    } */

    if (tid + start < end) 
    {
        int n = tid + start;

        int map0idx = d_cell_index[n];
        
        const int map1idx = oppDat1Map[map0idx + dep_charge_stride_OPP_DEVICE_1 * 0];
        const int map2idx = oppDat1Map[map0idx + dep_charge_stride_OPP_DEVICE_1 * 1];
        const int map3idx = oppDat1Map[map0idx + dep_charge_stride_OPP_DEVICE_1 * 2];
        const int map4idx = oppDat1Map[map0idx + dep_charge_stride_OPP_DEVICE_1 * 3];

        //user-supplied kernel call
        dep_node_charge__kernel_gpu(
            (dir_arg0 + n),
            (k_arg1 + map1idx),
            (k_arg1 + map2idx),
            (k_arg1 + map3idx),
            (k_arg1 + map4idx)
        );
    }

    /* // Uncomment For Shared memory 
    __syncthreads();

    if (use_shared_OPP_DEVICE)
    {
        for (int i = threadIdx.x; i < dep_charge_stride_OPP_DEVICE_1_SH; i+=gpu_threads_per_block_OPP_DEVICE) 
        {
            atomicAdd(&(ind_arg1[i]), shared_arg1[i]);
        }
        __syncthreads();
    } */
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
    opp_init_double_indirect_reductions_hip(nargs, args);
#endif

    if (set_size > 0) 
    {
        dep_charge_stride_OPP_HOST_0 = args[0].dat->set->set_capacity;
        dep_charge_stride_OPP_HOST_1 = args[1].size;

        hipMemcpyToSymbol(HIP_SYMBOL(dep_charge_stride_OPP_DEVICE_0), &dep_charge_stride_OPP_HOST_0, sizeof(int));
        hipMemcpyToSymbol(HIP_SYMBOL(dep_charge_stride_OPP_DEVICE_1), &dep_charge_stride_OPP_HOST_1, sizeof(int));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = OPP_gpu_threads_per_block;
            int nblocks = (end - start - 1) / nthread + 1;
            size_t shared_mem = 0;

            /* // Uncomment For Shared memory 
            dep_charge_stride_OPP_HOST_1_SH = args[1].dat->set->size + args[1].dat->set->exec_size + 
                                                args[1].dat->set->nonexec_size;
            hipMemcpyToSymbol(HIP_SYMBOL(dep_charge_stride_OPP_DEVICE_1_SH), &dep_charge_stride_OPP_HOST_1_SH, sizeof(int));
            hipMemcpyToSymbol(HIP_SYMBOL(gpu_threads_per_block_OPP_DEVICE), &OPP_gpu_threads_per_block, sizeof(int));

            shared_mem += args[1].dat->size * dep_charge_stride_OPP_HOST_1_SH;

            if (shared_mem > OPP_gpu_shared_mem_per_block) 
            {
                bool use_shared = false;
                hipMemcpyToSymbol(HIP_SYMBOL(use_shared_OPP_DEVICE), &use_shared, sizeof(bool));
                shared_mem = 1;
                
                if (FP_DEBUG) 
                    opp_printf("FEMPIC", "MoveToCells Not using shared mem [avail: %zu, need: %zu]", 
                        OPP_gpu_shared_mem_per_block, shared_mem);
            }
            else
            {
                bool use_shared = true;
                hipMemcpyToSymbol(HIP_SYMBOL(use_shared_OPP_DEVICE), &use_shared, sizeof(bool));

                if (FP_DEBUG) 
                    opp_printf("FEMPIC", "MoveToCells Using shared mem [avail: %zu, need: %zu]", 
                        OPP_gpu_shared_mem_per_block, shared_mem);
            }
            */

            opp_device_all_DepNodeCharge <<<nblocks, nthread, shared_mem>>> (
                (int *)     set->mesh_relation_dat->data_d,
                (double *)  args[0].data_d,
                (double *)  args[1].data_d,
                (int *)     args[1].map_data_d,
                (double *)  args[2].data_d,
                (double *)  args[3].data_d,
                (double *)  args[4].data_d,
                start, 
                end);
        }
    }

#ifdef USE_MPI
    opp_exchange_double_indirect_reductions_hip(nargs, args);
    opp_complete_double_indirect_reductions_hip(nargs, args);
#endif

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(hipDeviceSynchronize());

    opp_profiler->end("DepCharge");
}

#else // #ifdef USE_REDUCE_BY_KEY

thrust::device_vector<int> d_keys;
thrust::device_vector<double> d_values;
bool d_vectors_initialized = false;

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
        
        keys[tid + dep_charge_stride_OPP_DEVICE_X * 0] = oppDat1Map[map0idx + dep_charge_stride_OPP_DEVICE_1 * 0];
        keys[tid + dep_charge_stride_OPP_DEVICE_X * 1] = oppDat1Map[map0idx + dep_charge_stride_OPP_DEVICE_1 * 1];
        keys[tid + dep_charge_stride_OPP_DEVICE_X * 2] = oppDat1Map[map0idx + dep_charge_stride_OPP_DEVICE_1 * 2];
        keys[tid + dep_charge_stride_OPP_DEVICE_X * 3] = oppDat1Map[map0idx + dep_charge_stride_OPP_DEVICE_1 * 3];

        values[tid + dep_charge_stride_OPP_DEVICE_X * 0] = dir_arg0[n + dep_charge_stride_OPP_DEVICE_0 * 0];
        values[tid + dep_charge_stride_OPP_DEVICE_X * 1] = dir_arg0[n + dep_charge_stride_OPP_DEVICE_0 * 1];
        values[tid + dep_charge_stride_OPP_DEVICE_X * 2] = dir_arg0[n + dep_charge_stride_OPP_DEVICE_0 * 2];
        values[tid + dep_charge_stride_OPP_DEVICE_X * 3] = dir_arg0[n + dep_charge_stride_OPP_DEVICE_0 * 3];
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
        
        arg1[map0idx] += values[n];
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
    
    if (FP_DEBUG) 
        opp_printf("FEMPIC", "opp_loop_all__deposit_charge_on_nodes set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("DepCharge");

    int nargs = 5;
    opp_arg args[nargs];

    args[0]  = std::move(arg0);
    args[1]  = std::move(arg1);
    args[2]  = std::move(arg2);
    args[3]  = std::move(arg3);
    args[4]  = std::move(arg4);

    opp_profiler->start("Dep_preSync");
    cutilSafeCall(hipDeviceSynchronize());
    opp_profiler->end("Dep_preSync");

    opp_profiler->start("Dep_ExHalo");
    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_profiler->end("Dep_ExHalo");

    opp_profiler->start("Dep_ExHaloWait");
    opp_mpi_halo_wait_all(nargs, args);
    opp_profiler->end("Dep_ExHaloWait");

#ifdef USE_MPI
    opp_init_double_indirect_reductions_hip(nargs, args);
#endif

    if (set_size > 0) 
    {
        dep_charge_stride_OPP_HOST_0 = args[0].dat->set->set_capacity;
        dep_charge_stride_OPP_HOST_1 = args[1].size;

        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(dep_charge_stride_OPP_DEVICE_0), 
                                                    &dep_charge_stride_OPP_HOST_0, sizeof(int)));
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(dep_charge_stride_OPP_DEVICE_1), 
                                                    &dep_charge_stride_OPP_HOST_1, sizeof(int)));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int block_size = opp_params->get<OPP_INT>("sed_red_block");
            if (!d_vectors_initialized) // resize once
            {
                opp_profiler->start("Dep_Resize1");
                if (block_size > args[0].dat->set->set_capacity) 
                    block_size = args[0].dat->set->set_capacity;
                const size_t resize_size = (size_t)(block_size * args[0].dat->dim); 
                d_keys.resize(resize_size);
                d_values.resize(resize_size);
                d_vectors_initialized = true;
                opp_profiler->end("Dep_Resize1");
            }

            // Do blockwise reduction since sort_by_key takes exp time when size is large
            for (int start_internal = start; start_internal < end; start_internal += block_size)
            {
                int end_internal = (start_internal + block_size);
                if (end_internal > end) end_internal = end;

                dep_charge_stride_OPP_HOST_X = (end_internal - start_internal);
                cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(dep_charge_stride_OPP_DEVICE_X), 
                                                    &dep_charge_stride_OPP_HOST_X, sizeof(int)));

                const size_t operating_size = (size_t)(dep_charge_stride_OPP_HOST_X * args[0].dat->dim); //  no halo in args[0]

                const int nthread = OPP_gpu_threads_per_block;
                const int nblocks = (dep_charge_stride_OPP_HOST_X - 1) / nthread + 1;
                const size_t shared_mem = 0;

                // Create key/value pairs by nodes (linked with particle->cell->node) with part_lc values
                opp_profiler->start("Dep_CrKeyVal");
                create_key_value_pairs<<<nblocks, nthread, shared_mem>>> (
                    (int *)     set->mesh_relation_dat->data_d,
                    (int *)     thrust::raw_pointer_cast(d_keys.data()),
                    (double *)  thrust::raw_pointer_cast(d_values.data()),
                    (double *)  args[0].data_d,
                    (int *)     args[1].map_data_d,
                    start_internal, 
                    end_internal);
                cutilSafeCall(hipDeviceSynchronize());
                opp_profiler->end("Dep_CrKeyVal");

                // Sort by keys to bring the identical keys together
                opp_profiler->start("Dep_Sort");
                thrust::sort_by_key(thrust::device, d_keys.begin(), d_keys.begin() + operating_size, d_values.begin());
                opp_profiler->end("Dep_Sort");

                // Compute the unique keys and their corresponding values
                opp_profiler->start("Dep_Red");
                auto new_end = thrust::reduce_by_key(thrust::device,
                    d_keys.begin(), d_keys.begin() + operating_size,
                    d_values.begin(),
                    d_keys.begin(),
                    d_values.begin()
                );        
                opp_profiler->end("Dep_Red");

                // Assign reduced values to the nodes using keys/values
                opp_profiler->start("Dep_Assign");                
                assign_values<<<nblocks, nthread, shared_mem>>> (
                    (int *)     thrust::raw_pointer_cast(d_keys.data()),
                    (double *)  thrust::raw_pointer_cast(d_values.data()),
                    (double *)  args[1].data_d,
                    0, 
                    (int)(new_end.first - d_keys.begin()));
                cutilSafeCall(hipDeviceSynchronize());
                opp_profiler->end("Dep_Assign");
            }
        }
    }

#ifdef USE_MPI
    opp_exchange_double_indirect_reductions_hip(nargs, args);
    opp_complete_double_indirect_reductions_hip(nargs, args);
#endif

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(hipDeviceSynchronize());

    opp_profiler->end("DepCharge");
}


#endif
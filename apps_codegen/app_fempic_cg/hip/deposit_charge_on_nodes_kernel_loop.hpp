
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k5_dat0_stride = -1;
OPP_INT opp_k5_dat1_stride = -1;
OPP_INT opp_k5_map0_stride = -1;

__constant__ OPP_INT opp_k5_dat0_stride_d;
__constant__ OPP_INT opp_k5_dat1_stride_d;
__constant__ OPP_INT opp_k5_map0_stride_d;

OPP_INT opp_k5_sr_set_stride = -1;
__constant__ OPP_INT opp_k5_sr_set_stride_d;

thrust::device_vector<OPP_INT> sr_dat1_keys_dv;
thrust::device_vector<OPP_REAL> sr_dat1_values_dv;

thrust::device_vector<OPP_INT> sr_dat1_keys_dv2;
thrust::device_vector<OPP_REAL> sr_dat1_values_dv2;

namespace opp_k5 {
__device__ inline void deposit_charge_on_nodes_kernel(
    const double *part_lc,
    double *node_charge_den0,
    double *node_charge_den1,
    double *node_charge_den2,
    double *node_charge_den3
) {
    node_charge_den0[0] += part_lc[(0) * opp_k5_dat0_stride_d];
    node_charge_den1[0] += part_lc[(1) * opp_k5_dat0_stride_d];
    node_charge_den2[0] += part_lc[(2) * opp_k5_dat0_stride_d];
    node_charge_den3[0] += part_lc[(3) * opp_k5_dat0_stride_d];
}
}

//--------------------------------------------------------------
__global__ void opp_dev_deposit_charge_on_nodes_kernel(
    const OPP_REAL *__restrict__ dat0,  // p_lc
    OPP_REAL **__restrict__ dat1, const OPP_INT dat1_arr_count, // n_charge_den
    const OPP_INT *__restrict__ p2c_map,
    const OPP_INT *__restrict__ map0,  // c2n_map
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        const OPP_INT* opp_p2c = p2c_map + n;
        
        OPP_REAL arg1_0_local[1];
        for (int d = 0; d < 1; ++d)
            arg1_0_local[d] = OPP_REAL_ZERO;

        OPP_REAL arg2_1_local[1];
        for (int d = 0; d < 1; ++d)
            arg2_1_local[d] = OPP_REAL_ZERO;

        OPP_REAL arg3_2_local[1];
        for (int d = 0; d < 1; ++d)
            arg3_2_local[d] = OPP_REAL_ZERO;

        OPP_REAL arg4_3_local[1];
        for (int d = 0; d < 1; ++d)
            arg4_3_local[d] = OPP_REAL_ZERO;

        opp_k5::deposit_charge_on_nodes_kernel(
            dat0 + n, // p_lc 
            arg1_0_local, // n_charge_den 
            arg2_1_local, // n_charge_den 
            arg3_2_local, // n_charge_den 
            arg4_3_local // n_charge_den 
        );

        OPP_REAL* tmp1 = dat1[threadIdx.x % dat1_arr_count];

        for (int d = 0; d < 1; ++d)
            atomicAdd(tmp1 + map0[opp_k5_map0_stride_d * 0 + opp_p2c[0]] + (d * opp_k5_dat1_stride_d), arg1_0_local[d]);

        for (int d = 0; d < 1; ++d)
            atomicAdd(tmp1 + map0[opp_k5_map0_stride_d * 1 + opp_p2c[0]] + (d * opp_k5_dat1_stride_d), arg2_1_local[d]);

        for (int d = 0; d < 1; ++d)
            atomicAdd(tmp1 + map0[opp_k5_map0_stride_d * 2 + opp_p2c[0]] + (d * opp_k5_dat1_stride_d), arg3_2_local[d]);

        for (int d = 0; d < 1; ++d)
            atomicAdd(tmp1 + map0[opp_k5_map0_stride_d * 3 + opp_p2c[0]] + (d * opp_k5_dat1_stride_d), arg4_3_local[d]);
    }
    
}

//--------------------------------------------------------------
__global__ void opp_dev_sr_deposit_charge_on_nodes_kernel( // Used for Segmented Reductions
    const OPP_REAL *__restrict__ dat0,  // p_lc
    OPP_INT *__restrict__ sr_dat1_keys,     // sr keys for n_charge_den
    OPP_REAL *__restrict__ sr_dat1_values,     // sr values for n_charge_den
    const OPP_INT *__restrict__ p2c_map,
    const OPP_INT *__restrict__ map0,  // c2n_map
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        const OPP_INT* opp_p2c = p2c_map + n;

        OPP_REAL arg1_0_local[1];
        for (int d = 0; d < 1; ++d)
            arg1_0_local[d] = OPP_REAL_ZERO;

        OPP_REAL arg2_1_local[1];
        for (int d = 0; d < 1; ++d)
            arg2_1_local[d] = OPP_REAL_ZERO;

        OPP_REAL arg3_2_local[1];
        for (int d = 0; d < 1; ++d)
            arg3_2_local[d] = OPP_REAL_ZERO;

        OPP_REAL arg4_3_local[1];
        for (int d = 0; d < 1; ++d)
            arg4_3_local[d] = OPP_REAL_ZERO;

        opp_k5::deposit_charge_on_nodes_kernel(
            dat0 + n, // p_lc 
            arg1_0_local, // n_charge_den 
            arg2_1_local, // n_charge_den 
            arg3_2_local, // n_charge_den 
            arg4_3_local // n_charge_den 
        );

        for (int d = 0; d < 1; ++d) {
            sr_dat1_values[4 * n + 0 + opp_k5_sr_set_stride_d * d] = arg1_0_local[d]; 
        }    
        sr_dat1_keys[4 * n + 0] = map0[opp_k5_map0_stride_d * 0 + opp_p2c[0]];
        for (int d = 0; d < 1; ++d) {
            sr_dat1_values[4 * n + 1 + opp_k5_sr_set_stride_d * d] = arg2_1_local[d]; 
        }    
        sr_dat1_keys[4 * n + 1] = map0[opp_k5_map0_stride_d * 1 + opp_p2c[0]];
        for (int d = 0; d < 1; ++d) {
            sr_dat1_values[4 * n + 2 + opp_k5_sr_set_stride_d * d] = arg3_2_local[d]; 
        }    
        sr_dat1_keys[4 * n + 2] = map0[opp_k5_map0_stride_d * 2 + opp_p2c[0]];
        for (int d = 0; d < 1; ++d) {
            sr_dat1_values[4 * n + 3 + opp_k5_sr_set_stride_d * d] = arg4_3_local[d]; 
        }    
        sr_dat1_keys[4 * n + 3] = map0[opp_k5_map0_stride_d * 3 + opp_p2c[0]];
    }
    
}

void opp_par_loop_all__deposit_charge_on_nodes_kernel(opp_set set,
    opp_arg arg0, // p_lc | OPP_READ
    opp_arg arg1, // n_charge_den | OPP_INC
    opp_arg arg2, // n_charge_den | OPP_INC
    opp_arg arg3, // n_charge_den | OPP_INC
    opp_arg arg4 // n_charge_den | OPP_INC
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const int nargs = 5;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;

    opp_profiler->start("deposit_charge_on_nodes_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__deposit_charge_on_nodes_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);

#ifdef USE_MPI
    opp_init_double_indirect_reductions_device(nargs, args);
#endif
 
 
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k5_dat0_stride_d, &opp_k5_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k5_dat1_stride_d, &opp_k5_dat1_stride, &(args[1].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k5_map0_stride_d, &opp_k5_map0_stride, &(args[1].size), 1);

#ifdef OPP_BLOCK_SIZE_5
    const int block_size = OPP_BLOCK_SIZE_5;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    int num_blocks = 200;
    if (iter_size > 0) 
    {
        const OPP_INT start = 0;
        const OPP_INT end = iter_size;
        num_blocks = (end - start - 1) / block_size + 1;

        if (!opp_use_segmented_reductions) // Do atomics ----------       
        {
            const int array_count = opp_params->get<OPP_INT>("gpu_reduction_arrays");
            OPP_REAL** arg1_dat_thread_data_d = opp_create_thread_level_data<OPP_REAL>(args[1]);
            opp_dev_deposit_charge_on_nodes_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // p_lc
                arg1_dat_thread_data_d, array_count,     // n_charge_den
                (OPP_INT *)set->mesh_relation_dat->data_d,
                args[1].map_data_d,     // c2n_map
                start,
                end
            );
            OPP_DEVICE_SYNCHRONIZE(); 
            opp_reduce_thread_level_data<OPP_REAL>(args[1]);
        }
     
        else // Do segmented reductions ----------       
        {
            size_t operating_size_dat1 = 0, resize_size_dat1 = 0;

            operating_size_dat1 += (size_t)(args[1].dat->dim);
            resize_size_dat1 += (size_t)(args[1].dat->dim);
            operating_size_dat1 += (size_t)(args[2].dat->dim);
            resize_size_dat1 += (size_t)(args[2].dat->dim);
            operating_size_dat1 += (size_t)(args[3].dat->dim);
            resize_size_dat1 += (size_t)(args[3].dat->dim);
            operating_size_dat1 += (size_t)(args[4].dat->dim);
            resize_size_dat1 += (size_t)(args[4].dat->dim);

            operating_size_dat1 *= (size_t)(iter_size);
            resize_size_dat1 *= (size_t)(set->set_capacity);

            int k5_stride = 0;
            k5_stride += iter_size * 4;
            opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k5_sr_set_stride_d, &opp_k5_sr_set_stride, &(k5_stride), 1);

            opp_sr::init_arrays<OPP_REAL>(args[1].dat->dim, operating_size_dat1, resize_size_dat1,
                        sr_dat1_keys_dv, sr_dat1_values_dv, sr_dat1_keys_dv2, sr_dat1_values_dv2);
        
            // Create key/value pairs
            opp_profiler->start("SR_CrKeyVal");
            opp_dev_sr_deposit_charge_on_nodes_kernel<<<num_blocks, block_size>>>( 
                (OPP_REAL *)args[0].data_d,     // p_lc
                opp_get_dev_raw_ptr<OPP_INT>(sr_dat1_keys_dv),     // sr keys for n_charge_den
                opp_get_dev_raw_ptr<OPP_REAL>(sr_dat1_values_dv),     // sr values for n_charge_den
                (OPP_INT *)set->mesh_relation_dat->data_d,
                args[1].map_data_d,     // c2n_map
                start,
                end
            );
            OPP_DEVICE_SYNCHRONIZE();
            opp_profiler->end("SR_CrKeyVal");

            opp_sr::do_segmented_reductions<OPP_REAL>(args[1], k5_stride,
                        sr_dat1_keys_dv, sr_dat1_values_dv, sr_dat1_keys_dv2, sr_dat1_values_dv2);
        }
    }
    opp_sr::clear_arrays<OPP_REAL>(sr_dat1_keys_dv, sr_dat1_values_dv, sr_dat1_keys_dv2, sr_dat1_values_dv2);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   

#ifdef USE_MPI    
    opp_exchange_double_indirect_reductions_device(nargs, args);
    opp_complete_double_indirect_reductions_device(nargs, args);
#endif
 
    opp_profiler->end("deposit_charge_on_nodes_kernel");
}


//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k6_dat0_stride = -1;
OPP_INT opp_k6_dat1_stride = -1;
OPP_INT opp_k6_map0_stride = -1;

__constant__ OPP_INT opp_k6_dat0_stride_d;
__constant__ OPP_INT opp_k6_dat1_stride_d;
__constant__ OPP_INT opp_k6_map0_stride_d;

OPP_INT opp_k6_sr_set_stride = -1;
__constant__ OPP_INT opp_k6_sr_set_stride_d;

thrust::device_vector<OPP_INT> sr_dat1_keys_dv;
thrust::device_vector<OPP_REAL> sr_dat1_values_dv;

thrust::device_vector<OPP_INT> sr_dat1_keys_dv2;
thrust::device_vector<OPP_REAL> sr_dat1_values_dv2;

namespace opp_k6 {
__device__ inline void update_ghosts_kernel(
    const int* c_mask_ug,
    const double* from_cell,
    double* to_cell,
    const int* m_idx,
    const int* dim)
{
    if (c_mask_ug[(*m_idx) * opp_k6_dat0_stride_d] == 1)
        to_cell[*dim] += from_cell[(*dim) * opp_k6_dat1_stride_d];
}
}

//--------------------------------------------------------------
__global__ void opp_dev_update_ghosts_kernel(
    const OPP_INT *__restrict__ dat0,  // c_mask_ug
    OPP_REAL *__restrict__ dat1,  // c_j
    const OPP_INT *__restrict__ map0,  // c2cug0_map
    OPP_INT *gbl3,
    OPP_INT *gbl4,
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        OPP_REAL arg2_0_local[3];
        for (int d = 0; d < 3; ++d)
            arg2_0_local[d] = OPP_REAL_ZERO;

        opp_k6::update_ghosts_kernel(
            dat0 + n, // c_mask_ug 
            dat1 + n, // c_j 
            arg2_0_local, // c_j 
            gbl3, // 
            gbl4 // 
        );

        for (int d = 0; d < 3; ++d)
            atomicAdd(dat1 + map0[opp_k6_map0_stride_d * 0 + n] + (d * opp_k6_dat1_stride_d), arg2_0_local[d]);
    }
    
}

//--------------------------------------------------------------
__global__ void opp_dev_sr_update_ghosts_kernel( // Used for Segmented Reductions
    const OPP_INT *__restrict__ dat0,  // c_mask_ug
    OPP_INT *__restrict__ sr_dat1_keys,     // sr keys for c_j
    OPP_REAL *__restrict__ sr_dat1_values,     // sr values for c_j
    OPP_REAL *__restrict__ dat1,  // c_j
    const OPP_INT *__restrict__ map0,  // c2cug0_map
    OPP_INT *gbl3,
    OPP_INT *gbl4,
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;

        OPP_REAL arg2_0_local[3];
        for (int d = 0; d < 3; ++d)
            arg2_0_local[d] = OPP_REAL_ZERO;

        opp_k6::update_ghosts_kernel(
            dat0 + n, // c_mask_ug 
            dat1 + n, // c_j 
            arg2_0_local, // c_j 
            gbl3, // 
            gbl4 // 
        );

        for (int d = 0; d < 3; ++d) {
            sr_dat1_values[1 * n + 0 + opp_k6_sr_set_stride_d * d] = arg2_0_local[d]; 
        }    
        sr_dat1_keys[1 * n + 0] = map0[opp_k6_map0_stride_d * 0 + n];
    }
    
}

void opp_par_loop_all__update_ghosts_kernel(opp_set set,
    opp_arg arg0, // c_mask_ug | OPP_READ
    opp_arg arg1, // c_j | OPP_READ
    opp_arg arg2, // c_j | OPP_INC
    opp_arg arg3, // | OPP_READ
    opp_arg arg4 // | OPP_READ
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const int nargs = 5;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;

    opp_profiler->start("update_ghosts_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__update_ghosts_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    OPP_INT *arg3_host_data = (OPP_INT *)args[3].data;
    OPP_INT *arg4_host_data = (OPP_INT *)args[4].data;

    int const_bytes = 0;

    const_bytes += ROUND_UP(1 * sizeof(OPP_INT));
    const_bytes += ROUND_UP(1 * sizeof(OPP_INT));

    opp_reallocConstArrays(const_bytes);
    const_bytes = 0;

    args[3].data   = OPP_consts_h + const_bytes;
    args[3].data_d = OPP_consts_d + const_bytes;

    for (int d = 0; d < 1; ++d)
        ((OPP_INT *)args[3].data)[d] = arg3_host_data[d];

    const_bytes += ROUND_UP(1 * sizeof(OPP_INT));

    args[4].data   = OPP_consts_h + const_bytes;
    args[4].data_d = OPP_consts_d + const_bytes;

    for (int d = 0; d < 1; ++d)
        ((OPP_INT *)args[4].data)[d] = arg4_host_data[d];

    const_bytes += ROUND_UP(1 * sizeof(OPP_INT));

    opp_mvConstArraysToDevice(const_bytes);

    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k6_dat0_stride_d, &opp_k6_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k6_dat1_stride_d, &opp_k6_dat1_stride, &(args[1].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k6_map0_stride_d, &opp_k6_map0_stride, &(args[2].size), 1);

#ifdef OPP_BLOCK_SIZE_6
    const int block_size = OPP_BLOCK_SIZE_6;
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
            opp_dev_update_ghosts_kernel<<<num_blocks, block_size>>>(
                (OPP_INT *)args[0].data_d,     // c_mask_ug
                (OPP_REAL *)args[1].data_d,     // c_j
                args[2].map_data_d,     // c2cug0_map
                (OPP_INT *)args[3].data_d,
                (OPP_INT *)args[4].data_d,
                start,
                end
            );
        }
     
        else // Do segmented reductions ----------       
        {
            size_t operating_size_dat1 = 0, resize_size_dat1 = 0;

            operating_size_dat1 += (size_t)(args[2].dat->dim);
            resize_size_dat1 += (size_t)(args[2].dat->dim);

            operating_size_dat1 *= (size_t)(iter_size);
            resize_size_dat1 *= (size_t)(set->set_capacity);

            int k6_stride = 0;
            k6_stride += iter_size * 1;
            opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k6_sr_set_stride_d, &opp_k6_sr_set_stride, &(k6_stride), 1);

            opp_sr::init_arrays<OPP_REAL>(args[2].dat->dim, operating_size_dat1, resize_size_dat1,
                        sr_dat1_keys_dv, sr_dat1_values_dv, sr_dat1_keys_dv2, sr_dat1_values_dv2);
        
            // Create key/value pairs
            opp_profiler->start("SR_CrKeyVal");
            opp_dev_sr_update_ghosts_kernel<<<num_blocks, block_size>>>( 
                (OPP_INT *)args[0].data_d,     // c_mask_ug
                opp_get_dev_raw_ptr<OPP_INT>(sr_dat1_keys_dv),     // sr keys for c_j
                opp_get_dev_raw_ptr<OPP_REAL>(sr_dat1_values_dv),     // sr values for c_j
                (OPP_REAL *)args[1].data_d,     // c_j
                args[2].map_data_d,     // c2cug0_map
                (OPP_INT *)args[3].data_d,
                (OPP_INT *)args[4].data_d,
                start,
                end
            );
            OPP_DEVICE_SYNCHRONIZE();
            opp_profiler->end("SR_CrKeyVal");

            opp_sr::do_segmented_reductions<OPP_REAL>(args[2], k6_stride,
                        sr_dat1_keys_dv, sr_dat1_values_dv, sr_dat1_keys_dv2, sr_dat1_values_dv2);
        }
    }
    args[3].data = (char *)arg3_host_data;
    args[4].data = (char *)arg4_host_data;

    opp_sr::clear_arrays<OPP_REAL>(sr_dat1_keys_dv, sr_dat1_values_dv, sr_dat1_keys_dv2, sr_dat1_values_dv2);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("update_ghosts_kernel");
}

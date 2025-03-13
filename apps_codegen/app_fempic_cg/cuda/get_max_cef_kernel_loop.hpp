
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k8_dat0_stride = -1;

__constant__ OPP_INT opp_k8_dat0_stride_d;



namespace opp_k8 {
__device__ inline void get_max_cef_kernel(
    const double* val,
    double* max_val
) {
    for (int dim = 0; dim < 3; ++dim) {
        *max_val = ((val[(dim) * opp_k8_dat0_stride_d] > *max_val) ? (val[(dim) * opp_k8_dat0_stride_d]) : (*max_val));
    }
}

}

//--------------------------------------------------------------
__global__ void opp_dev_get_max_cef_kernel(
    const OPP_REAL *__restrict__ dat0,  // c_ef
    OPP_REAL *gbl1,
    const OPP_INT start,
    const OPP_INT end
) 
{
    OPP_REAL gbl1_local[1];
    for (int d = 0; d < 1; ++d)
        gbl1_local[d] = gbl1[blockIdx.x * 1 + d];

    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    for (int n = thread_id; n < (end - start); n += blockDim.x * gridDim.x) {

        opp_k8::get_max_cef_kernel(
            dat0 + n, // c_ef 
            gbl1_local // 
        );
    }

    for (int d = 0; d < 1; ++d)
        opp_reduction<OPP_MAX>(gbl1 + blockIdx.x * 1 + d, gbl1_local[d]);
    
}

//--------------------------------------------------------------
void opp_par_loop_all__get_max_cef_kernel(opp_set set,
    opp_arg arg0, // c_ef | OPP_READ
    opp_arg arg1 // | OPP_MAX
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("get_max_cef_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__get_max_cef_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    OPP_REAL *arg1_host_data = (OPP_REAL *)args[1].data;

    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k8_dat0_stride_d, &opp_k8_dat0_stride, &(args[0].dat->set->set_capacity), 1);

#ifdef OPP_BLOCK_SIZE_8
    const int block_size = OPP_BLOCK_SIZE_8;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    int num_blocks = 200;

    int max_blocks = num_blocks;

    int reduction_bytes = 0;
    int reduction_size = 0;

    reduction_bytes += ROUND_UP(max_blocks * 1 * sizeof(OPP_REAL));
    reduction_size   = MAX(reduction_size, sizeof(OPP_REAL));

    opp_reallocReductArrays(reduction_bytes);
    reduction_bytes = 0;

    args[1].data   = OPP_reduct_h + reduction_bytes;
    args[1].data_d = OPP_reduct_d + reduction_bytes;

    for (int b = 0; b < max_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            ((OPP_REAL *)args[1].data)[b * 1 + d] = arg1_host_data[d];
    }

    reduction_bytes += ROUND_UP(max_blocks * 1 * sizeof(OPP_REAL));

    opp_mvReductArraysToDevice(reduction_bytes);
    
    if (iter_size > 0) 
    {
        const OPP_INT start = 0;
        const OPP_INT end = iter_size;
        {
            opp_dev_get_max_cef_kernel<<<num_blocks, block_size, (reduction_size * block_size)>>>(
                (OPP_REAL *)args[0].data_d,     // c_ef
                (OPP_REAL *)args[1].data_d,
                start,
                end
            );
        }
    }

    opp_mvReductArraysToHost(reduction_bytes);

    for (int b = 0; b < max_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            arg1_host_data[d] = MAX(arg1_host_data[d], ((OPP_REAL *)args[1].data)[b * 1 + d]);
    }

    args[1].data = (char *)arg1_host_data;
    opp_mpi_reduce(&args[1], arg1_host_data);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("get_max_cef_kernel");
}


//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k8_dat0_stride = -1;
OPP_INT opp_k8_dat1_stride = -1;

__constant__ OPP_INT opp_k8_dat0_stride_d;
__constant__ OPP_INT opp_k8_dat1_stride_d;



namespace opp_k8 {
enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

__device__ inline void compute_energy_kernel(
    const int* cell0_ghost,
    const double* cell_field,
    double* energy)
{
    if (cell0_ghost[(0) * opp_k8_dat0_stride_d] == 0)
    {
        energy[0] += cell_field[(Dim::x) * opp_k8_dat1_stride_d] * cell_field[(Dim::x) * opp_k8_dat1_stride_d] +
                    cell_field[(Dim::y) * opp_k8_dat1_stride_d] * cell_field[(Dim::y) * opp_k8_dat1_stride_d] +
                    cell_field[(Dim::z) * opp_k8_dat1_stride_d] * cell_field[(Dim::z) * opp_k8_dat1_stride_d];
    }
}

}

//--------------------------------------------------------------
__global__ void opp_dev_compute_energy_kernel(
    const OPP_INT *__restrict__ dat0,  // c_mask_ghost
    const OPP_REAL *__restrict__ dat1,  // c_e
    OPP_REAL *gbl2,
    const OPP_INT start,
    const OPP_INT end
) 
{
    OPP_REAL gbl2_local[1];
    for (int d = 0; d < 1; ++d)
        gbl2_local[d] = OPP_REAL_ZERO;

    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    for (int n = thread_id; n < (end - start); n += blockDim.x * gridDim.x) {

        opp_k8::compute_energy_kernel(
            dat0 + n, // c_mask_ghost 
            dat1 + n, // c_e 
            gbl2_local // 
        );
    }

    for (int d = 0; d < 1; ++d)
        opp_reduction<OPP_INC>(gbl2 + blockIdx.x * 1 + d, gbl2_local[d]);
    
}

//--------------------------------------------------------------
void opp_par_loop_all__compute_energy_kernel(opp_set set,
    opp_arg arg0, // c_mask_ghost | OPP_READ
    opp_arg arg1, // c_e | OPP_READ
    opp_arg arg2 // | OPP_INC
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("compute_energy_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__compute_energy_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    OPP_REAL *arg2_host_data = (OPP_REAL *)args[2].data;

    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k8_dat0_stride_d, &opp_k8_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k8_dat1_stride_d, &opp_k8_dat1_stride, &(args[1].dat->set->set_capacity), 1);

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

    args[2].data   = OPP_reduct_h + reduction_bytes;
    args[2].data_d = OPP_reduct_d + reduction_bytes;

    for (int b = 0; b < max_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            ((OPP_REAL *)args[2].data)[b * 1 + d] = OPP_REAL_ZERO;
    }

    reduction_bytes += ROUND_UP(max_blocks * 1 * sizeof(OPP_REAL));

    opp_mvReductArraysToDevice(reduction_bytes);
    
    if (iter_size > 0) 
    {
        const OPP_INT start = 0;
        const OPP_INT end = iter_size;
        {
            opp_dev_compute_energy_kernel<<<num_blocks, block_size, (reduction_size * block_size)>>>(
                (OPP_INT *)args[0].data_d,     // c_mask_ghost
                (OPP_REAL *)args[1].data_d,     // c_e
                (OPP_REAL *)args[2].data_d,
                start,
                end
            );
        }
    }

    opp_mvReductArraysToHost(reduction_bytes);

    for (int b = 0; b < max_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            arg2_host_data[d] += ((OPP_REAL *)args[2].data)[b * 1 + d];
    }

    args[2].data = (char *)arg2_host_data;
    opp_mpi_reduce(&args[2], arg2_host_data);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("compute_energy_kernel");
}

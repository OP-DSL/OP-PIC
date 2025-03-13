
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k5_dat0_stride = -1;
OPP_INT opp_k5_dat1_stride = -1;
OPP_INT opp_k5_map0_stride = -1;

__constant__ OPP_INT opp_k5_dat0_stride_d;
__constant__ OPP_INT opp_k5_dat1_stride_d;
__constant__ OPP_INT opp_k5_map0_stride_d;



namespace opp_k5 {
enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

__device__ inline void update_ghosts_B_kernel(
    const int* c_mask_ugb,
    const double* from_cell,
    double* to_cell,
    const int* m_idx)
{
    if (c_mask_ugb[(*m_idx) * opp_k5_dat0_stride_d] == 1)
    {
        to_cell[(Dim::x) * opp_k5_dat1_stride_d] = from_cell[(Dim::x) * opp_k5_dat1_stride_d];
        to_cell[(Dim::y) * opp_k5_dat1_stride_d] = from_cell[(Dim::y) * opp_k5_dat1_stride_d];
        to_cell[(Dim::z) * opp_k5_dat1_stride_d] = from_cell[(Dim::z) * opp_k5_dat1_stride_d];
    }
}

}

//--------------------------------------------------------------
__global__ void opp_dev_update_ghosts_B_kernel(
    const OPP_INT *__restrict__ dat0,  // c_mask_ugb
    OPP_REAL *__restrict__ dat1,  // c_b
    const OPP_INT *__restrict__ map0,  // c2cugb0_map
    OPP_INT *gbl3,
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        opp_k5::update_ghosts_B_kernel(
            dat0 + n, // c_mask_ugb 
            dat1 + n, // c_b 
            dat1 + map0[n + opp_k5_map0_stride_d * 0], // c_b 
            gbl3 // 
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__update_ghosts_B_kernel(opp_set set,
    opp_arg arg0, // c_mask_ugb | OPP_READ
    opp_arg arg1, // c_b | OPP_READ
    opp_arg arg2, // c_b | OPP_WRITE
    opp_arg arg3 // | OPP_READ
) 
{
    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    opp_profiler->start("update_ghosts_B_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__update_ghosts_B_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    OPP_INT *arg3_host_data = (OPP_INT *)args[3].data;

    int const_bytes = 0;

    const_bytes += ROUND_UP(1 * sizeof(OPP_INT));

    opp_reallocConstArrays(const_bytes);
    const_bytes = 0;

    args[3].data   = OPP_consts_h + const_bytes;
    args[3].data_d = OPP_consts_d + const_bytes;

    for (int d = 0; d < 1; ++d)
        ((OPP_INT *)args[3].data)[d] = arg3_host_data[d];

    const_bytes += ROUND_UP(1 * sizeof(OPP_INT));

    opp_mvConstArraysToDevice(const_bytes);

    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k5_dat0_stride_d, &opp_k5_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k5_dat1_stride_d, &opp_k5_dat1_stride, &(args[1].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k5_map0_stride_d, &opp_k5_map0_stride, &(args[2].size), 1);

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

        {
            opp_dev_update_ghosts_B_kernel<<<num_blocks, block_size>>>(
                (OPP_INT *)args[0].data_d,     // c_mask_ugb
                (OPP_REAL *)args[1].data_d,     // c_b
                args[2].map_data_d,     // c2cugb0_map
                (OPP_INT *)args[3].data_d,
                start,
                end
            );
        }
    }
    args[3].data = (char *)arg3_host_data;


    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("update_ghosts_B_kernel");
}

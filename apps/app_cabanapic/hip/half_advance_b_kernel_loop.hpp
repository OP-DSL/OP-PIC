
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k4_dat0_stride = -1;
OPP_INT opp_k4_dat1_stride = -1;
OPP_INT opp_k4_dat2_stride = -1;
OPP_INT opp_k4_map0_stride = -1;

__constant__ OPP_INT opp_k4_dat0_stride_d;
__constant__ OPP_INT opp_k4_dat1_stride_d;
__constant__ OPP_INT opp_k4_dat2_stride_d;
__constant__ OPP_INT opp_k4_map0_stride_d;



namespace opp_k4 {
enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

__device__ inline void half_advance_b_kernel (
    const double* cell_x_e,
    const double* cell_y_e,
    const double* cell_z_e,
    const double* cell0_e,
    double* cell0_b,
    const int* cell0_ghost)
{
    if (cell0_ghost[(0) * opp_k4_dat2_stride_d] == 0)
    {
        cell0_b[(Dim::x) * opp_k4_dat1_stride_d] -= (0.5 * CONST_p_d[Dim::y] * (cell_y_e[(Dim::z) * opp_k4_dat0_stride_d] - cell0_e[(Dim::z) * opp_k4_dat0_stride_d])
                            - 0.5 * CONST_p_d[Dim::z] * (cell_z_e[(Dim::y) * opp_k4_dat0_stride_d] - cell0_e[(Dim::y) * opp_k4_dat0_stride_d]));

        cell0_b[(Dim::y) * opp_k4_dat1_stride_d] -= (0.5 * CONST_p_d[Dim::z] * (cell_z_e[(Dim::x) * opp_k4_dat0_stride_d] - cell0_e[(Dim::x) * opp_k4_dat0_stride_d])
                            - 0.5 * CONST_p_d[Dim::x] * (cell_x_e[(Dim::z) * opp_k4_dat0_stride_d] - cell0_e[(Dim::z) * opp_k4_dat0_stride_d]));

        cell0_b[(Dim::z) * opp_k4_dat1_stride_d] -= (0.5 * CONST_p_d[Dim::x] * (cell_x_e[(Dim::y) * opp_k4_dat0_stride_d] - cell0_e[(Dim::y) * opp_k4_dat0_stride_d])
                            - 0.5 * CONST_p_d[Dim::y] * (cell_y_e[(Dim::x) * opp_k4_dat0_stride_d] - cell0_e[(Dim::x) * opp_k4_dat0_stride_d]));
    }
}

}

//--------------------------------------------------------------
__global__ void opp_dev_half_advance_b_kernel(
    const OPP_REAL *__restrict__ dat0,  // c_e
    OPP_REAL *__restrict__ dat1,  // c_b
    const OPP_INT *__restrict__ dat2,  // c_mask_ghost
    const OPP_INT *__restrict__ map0,  // c2c_map
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        opp_k4::half_advance_b_kernel(
            dat0 + map0[n + opp_k4_map0_stride_d * 9], // c_e 
            dat0 + map0[n + opp_k4_map0_stride_d * 7], // c_e 
            dat0 + map0[n + opp_k4_map0_stride_d * 6], // c_e 
            dat0 + n, // c_e 
            dat1 + n, // c_b 
            dat2 + n // c_mask_ghost 
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__half_advance_b_kernel(opp_set set,
    opp_arg arg0, // c_e | OPP_READ
    opp_arg arg1, // c_e | OPP_READ
    opp_arg arg2, // c_e | OPP_READ
    opp_arg arg3, // c_e | OPP_READ
    opp_arg arg4, // c_b | OPP_INC
    opp_arg arg5 // c_mask_ghost | OPP_READ
) 
{
    const int nargs = 6;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;

    opp_profiler->start("half_advance_b_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__half_advance_b_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat0_stride_d, &opp_k4_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat1_stride_d, &opp_k4_dat1_stride, &(args[4].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat2_stride_d, &opp_k4_dat2_stride, &(args[5].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_map0_stride_d, &opp_k4_map0_stride, &(args[0].size), 1);

#ifdef OPP_BLOCK_SIZE_4
    const int block_size = OPP_BLOCK_SIZE_4;
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
            opp_dev_half_advance_b_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // c_e
                (OPP_REAL *)args[4].data_d,     // c_b
                (OPP_INT *)args[5].data_d,     // c_mask_ghost
                args[0].map_data_d,     // c2c_map
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("half_advance_b_kernel");
}

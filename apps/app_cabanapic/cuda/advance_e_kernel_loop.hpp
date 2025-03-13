
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k7_dat0_stride = -1;
OPP_INT opp_k7_dat1_stride = -1;
OPP_INT opp_k7_dat2_stride = -1;
OPP_INT opp_k7_dat3_stride = -1;
OPP_INT opp_k7_map0_stride = -1;

__constant__ OPP_INT opp_k7_dat0_stride_d;
__constant__ OPP_INT opp_k7_dat1_stride_d;
__constant__ OPP_INT opp_k7_dat2_stride_d;
__constant__ OPP_INT opp_k7_dat3_stride_d;
__constant__ OPP_INT opp_k7_map0_stride_d;



namespace opp_k7 {
enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

__device__ inline void advance_e_kernel (
    const double* cell_x_b,
    const double* cell_y_b,
    const double* cell_z_b,
    const double* cell0_b,
    const double* cell0_j,
    double* cell0_e,
    const int* iter_adv_e)
{
    if (iter_adv_e[(0) * opp_k7_dat3_stride_d] == 1)
    {
        cell0_e[(Dim::x) * opp_k7_dat2_stride_d] += ( - CONST_dt_eps0_d[0] * cell0_j[(Dim::x) * opp_k7_dat1_stride_d] ) +
            ( CONST_p_d[Dim::y] * (cell0_b[(Dim::z) * opp_k7_dat0_stride_d] - cell_y_b[(Dim::z) * opp_k7_dat0_stride_d]) -
            CONST_p_d[Dim::z] * (cell0_b[(Dim::y) * opp_k7_dat0_stride_d] - cell_z_b[(Dim::y) * opp_k7_dat0_stride_d]) );

        cell0_e[(Dim::y) * opp_k7_dat2_stride_d] += ( - CONST_dt_eps0_d[0] * cell0_j[(Dim::y) * opp_k7_dat1_stride_d] ) +
            ( CONST_p_d[Dim::z] * (cell0_b[(Dim::x) * opp_k7_dat0_stride_d] - cell_z_b[(Dim::x) * opp_k7_dat0_stride_d]) -
            CONST_p_d[Dim::x] * (cell0_b[(Dim::z) * opp_k7_dat0_stride_d] - cell_x_b[(Dim::z) * opp_k7_dat0_stride_d]) );

        cell0_e[(Dim::z) * opp_k7_dat2_stride_d] += ( - CONST_dt_eps0_d[0] * cell0_j[(Dim::z) * opp_k7_dat1_stride_d] ) +
            ( CONST_p_d[Dim::x] * (cell0_b[(Dim::y) * opp_k7_dat0_stride_d] - cell_x_b[(Dim::y) * opp_k7_dat0_stride_d]) -
            CONST_p_d[Dim::y] * (cell0_b[(Dim::x) * opp_k7_dat0_stride_d] - cell_y_b[(Dim::x) * opp_k7_dat0_stride_d]) );
    }
}

}

//--------------------------------------------------------------
__global__ void opp_dev_advance_e_kernel(
    const OPP_REAL *__restrict__ dat0,  // c_b
    const OPP_REAL *__restrict__ dat1,  // c_j
    OPP_REAL *__restrict__ dat2,  // c_e
    const OPP_INT *__restrict__ dat3,  // c_mask_right
    const OPP_INT *__restrict__ map0,  // c2c_map
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        opp_k7::advance_e_kernel(
            dat0 + map0[n + opp_k7_map0_stride_d * 2], // c_b 
            dat0 + map0[n + opp_k7_map0_stride_d * 4], // c_b 
            dat0 + map0[n + opp_k7_map0_stride_d * 5], // c_b 
            dat0 + n, // c_b 
            dat1 + n, // c_j 
            dat2 + n, // c_e 
            dat3 + n // c_mask_right 
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__advance_e_kernel(opp_set set,
    opp_arg arg0, // c_b | OPP_READ
    opp_arg arg1, // c_b | OPP_READ
    opp_arg arg2, // c_b | OPP_READ
    opp_arg arg3, // c_b | OPP_READ
    opp_arg arg4, // c_j | OPP_READ
    opp_arg arg5, // c_e | OPP_INC
    opp_arg arg6 // c_mask_right | OPP_READ
) 
{
    const int nargs = 7;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;
    args[6] = arg6;

    opp_profiler->start("advance_e_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__advance_e_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k7_dat0_stride_d, &opp_k7_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k7_dat1_stride_d, &opp_k7_dat1_stride, &(args[4].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k7_dat2_stride_d, &opp_k7_dat2_stride, &(args[5].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k7_dat3_stride_d, &opp_k7_dat3_stride, &(args[6].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k7_map0_stride_d, &opp_k7_map0_stride, &(args[0].size), 1);

#ifdef OPP_BLOCK_SIZE_7
    const int block_size = OPP_BLOCK_SIZE_7;
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
            opp_dev_advance_e_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // c_b
                (OPP_REAL *)args[4].data_d,     // c_j
                (OPP_REAL *)args[5].data_d,     // c_e
                (OPP_INT *)args[6].data_d,     // c_mask_right
                args[0].map_data_d,     // c2c_map
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("advance_e_kernel");
}

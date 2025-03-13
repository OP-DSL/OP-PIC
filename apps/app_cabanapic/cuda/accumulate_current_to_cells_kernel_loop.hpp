
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k3_dat0_stride = -1;
OPP_INT opp_k3_dat1_stride = -1;
OPP_INT opp_k3_dat2_stride = -1;
OPP_INT opp_k3_map0_stride = -1;

__constant__ OPP_INT opp_k3_dat0_stride_d;
__constant__ OPP_INT opp_k3_dat1_stride_d;
__constant__ OPP_INT opp_k3_dat2_stride_d;
__constant__ OPP_INT opp_k3_map0_stride_d;



namespace opp_k3 {
enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

enum CellAcc {
    jfx = 0 * 4,
    jfy = 1 * 4,
    jfz = 2 * 4,
};

__device__ inline void accumulate_current_to_cells_kernel(
        double* cell0_j,
        const double* cell0_acc,
        const double* cell_xd_acc,
        const double* cell_yd_acc,
        const double* cell_zd_acc,
        const double* cell_xyd_acc,
        const double* cell_yzd_acc,
        const double* cell_xzd_acc,
        const int* iter_acc)
{
    if (iter_acc[(0) * opp_k3_dat2_stride_d] == 1)
    {
        cell0_j[(Dim::x) * opp_k3_dat0_stride_d] = CONST_acc_coef_d[Dim::x] * (cell0_acc[(CellAcc::jfx + 0) * opp_k3_dat1_stride_d] +
                                                    cell_yd_acc[(CellAcc::jfx + 1) * opp_k3_dat1_stride_d] +
                                                    cell_zd_acc[(CellAcc::jfx + 2) * opp_k3_dat1_stride_d] +
                                                    cell_yzd_acc[(CellAcc::jfx + 3) * opp_k3_dat1_stride_d]);

        cell0_j[(Dim::y) * opp_k3_dat0_stride_d] = CONST_acc_coef_d[Dim::y] * (cell0_acc[(CellAcc::jfy + 0) * opp_k3_dat1_stride_d] +
                                                    cell_zd_acc[(CellAcc::jfy + 1) * opp_k3_dat1_stride_d] +
                                                    cell_xd_acc[(CellAcc::jfy + 2) * opp_k3_dat1_stride_d] +
                                                    cell_xzd_acc[(CellAcc::jfy + 3) * opp_k3_dat1_stride_d]);

        cell0_j[(Dim::z) * opp_k3_dat0_stride_d] = CONST_acc_coef_d[Dim::z] * (cell0_acc[(CellAcc::jfz + 0) * opp_k3_dat1_stride_d] +
                                                    cell_xd_acc[(CellAcc::jfz + 1) * opp_k3_dat1_stride_d] +
                                                    cell_yd_acc[(CellAcc::jfz + 2) * opp_k3_dat1_stride_d] +
                                                    cell_xyd_acc[(CellAcc::jfz + 3) * opp_k3_dat1_stride_d]);
    }
}

}

//--------------------------------------------------------------
__global__ void opp_dev_accumulate_current_to_cells_kernel(
    OPP_REAL *__restrict__ dat0,  // c_j
    const OPP_REAL *__restrict__ dat1,  // c_acc
    const OPP_INT *__restrict__ dat2,  // c_mask_right
    const OPP_INT *__restrict__ map0,  // c2c_map
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        opp_k3::accumulate_current_to_cells_kernel(
            dat0 + n, // c_j 
            dat1 + n, // c_acc 
            dat1 + map0[n + opp_k3_map0_stride_d * 2], // c_acc 
            dat1 + map0[n + opp_k3_map0_stride_d * 4], // c_acc 
            dat1 + map0[n + opp_k3_map0_stride_d * 5], // c_acc 
            dat1 + map0[n + opp_k3_map0_stride_d * 0], // c_acc 
            dat1 + map0[n + opp_k3_map0_stride_d * 3], // c_acc 
            dat1 + map0[n + opp_k3_map0_stride_d * 1], // c_acc 
            dat2 + n // c_mask_right 
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__accumulate_current_to_cells_kernel(opp_set set,
    opp_arg arg0, // c_j | OPP_WRITE
    opp_arg arg1, // c_acc | OPP_READ
    opp_arg arg2, // c_acc | OPP_READ
    opp_arg arg3, // c_acc | OPP_READ
    opp_arg arg4, // c_acc | OPP_READ
    opp_arg arg5, // c_acc | OPP_READ
    opp_arg arg6, // c_acc | OPP_READ
    opp_arg arg7, // c_acc | OPP_READ
    opp_arg arg8 // c_mask_right | OPP_READ
) 
{
    const int nargs = 9;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;
    args[6] = arg6;
    args[7] = arg7;
    args[8] = arg8;

    opp_profiler->start("accumulate_current_to_cells_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__accumulate_current_to_cells_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k3_dat0_stride_d, &opp_k3_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k3_dat1_stride_d, &opp_k3_dat1_stride, &(args[1].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k3_dat2_stride_d, &opp_k3_dat2_stride, &(args[8].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k3_map0_stride_d, &opp_k3_map0_stride, &(args[2].size), 1);

#ifdef OPP_BLOCK_SIZE_3
    const int block_size = OPP_BLOCK_SIZE_3;
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
            opp_dev_accumulate_current_to_cells_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // c_j
                (OPP_REAL *)args[1].data_d,     // c_acc
                (OPP_INT *)args[8].data_d,     // c_mask_right
                args[2].map_data_d,     // c2c_map
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("accumulate_current_to_cells_kernel");
}


//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k7_dat0_stride = -1;
OPP_INT opp_k7_dat1_stride = -1;
OPP_INT opp_k7_dat2_stride = -1;
OPP_INT opp_k7_map0_stride = -1;

__constant__ OPP_INT opp_k7_dat0_stride_d;
__constant__ OPP_INT opp_k7_dat1_stride_d;
__constant__ OPP_INT opp_k7_dat2_stride_d;
__constant__ OPP_INT opp_k7_map0_stride_d;



namespace opp_k7 {
__device__ inline void compute_electric_field_kernel(
    double *cell_electric_field,
    const double *cell_shape_deriv,
    const double *node_potential0,
    const double *node_potential1,
    const double *node_potential2,
    const double *node_potential3
) {
    for (int dim = 0; dim < 3; dim++) {
        const double c1 = (cell_shape_deriv[(0 * 3 + dim) * opp_k7_dat1_stride_d] * (*node_potential0));
        const double c2 = (cell_shape_deriv[(1 * 3 + dim) * opp_k7_dat1_stride_d] * (*node_potential1));
        const double c3 = (cell_shape_deriv[(2 * 3 + dim) * opp_k7_dat1_stride_d] * (*node_potential2));
        const double c4 = (cell_shape_deriv[(3 * 3 + dim) * opp_k7_dat1_stride_d] * (*node_potential3));

        cell_electric_field[(dim) * opp_k7_dat0_stride_d] -= (c1 + c2 + c3 + c4);
    }
}

}

//--------------------------------------------------------------
__global__ void opp_dev_compute_electric_field_kernel(
    OPP_REAL *__restrict__ dat0,  // c_ef
    const OPP_REAL *__restrict__ dat1,  // c_sd
    const OPP_REAL *__restrict__ dat2,  // n_potential
    const OPP_INT *__restrict__ map0,  // c2n_map
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        opp_k7::compute_electric_field_kernel(
            dat0 + n, // c_ef 
            dat1 + n, // c_sd 
            dat2 + map0[n + opp_k7_map0_stride_d * 0], // n_potential 
            dat2 + map0[n + opp_k7_map0_stride_d * 1], // n_potential 
            dat2 + map0[n + opp_k7_map0_stride_d * 2], // n_potential 
            dat2 + map0[n + opp_k7_map0_stride_d * 3] // n_potential 
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__compute_electric_field_kernel(opp_set set,
    opp_arg arg0, // c_ef | OPP_INC
    opp_arg arg1, // c_sd | OPP_READ
    opp_arg arg2, // n_potential | OPP_READ
    opp_arg arg3, // n_potential | OPP_READ
    opp_arg arg4, // n_potential | OPP_READ
    opp_arg arg5 // n_potential | OPP_READ
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const int nargs = 6;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;

    opp_profiler->start("compute_electric_field_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__compute_electric_field_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k7_dat0_stride_d, &opp_k7_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k7_dat1_stride_d, &opp_k7_dat1_stride, &(args[1].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k7_dat2_stride_d, &opp_k7_dat2_stride, &(args[2].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k7_map0_stride_d, &opp_k7_map0_stride, &(args[2].size), 1);

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
            opp_dev_compute_electric_field_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // c_ef
                (OPP_REAL *)args[1].data_d,     // c_sd
                (OPP_REAL *)args[2].data_d,     // n_potential
                args[2].map_data_d,     // c2n_map
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("compute_electric_field_kernel");
}

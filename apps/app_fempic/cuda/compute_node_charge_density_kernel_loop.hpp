
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k6_dat0_stride = -1;
OPP_INT opp_k6_dat1_stride = -1;

__constant__ OPP_INT opp_k6_dat0_stride_d;
__constant__ OPP_INT opp_k6_dat1_stride_d;



namespace opp_k6 {
__device__ inline void compute_node_charge_density_kernel(
    double *node_charge_den,
    const double *node_volume
) {
    (*node_charge_den) *= (CONST_spwt_d[0] / node_volume[(0) * opp_k6_dat1_stride_d]);
}

}

//--------------------------------------------------------------
__global__ void opp_dev_compute_node_charge_density_kernel(
    OPP_REAL *__restrict__ dat0,  // n_charge_den
    const OPP_REAL *__restrict__ dat1,  // n_volume
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    for (int n = thread_id; n < (end - start); n += blockDim.x * gridDim.x) {

        opp_k6::compute_node_charge_density_kernel(
            dat0 + n, // n_charge_den 
            dat1 + n // n_volume 
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__compute_node_charge_density_kernel(opp_set set,
    opp_arg arg0, // n_charge_den | OPP_RW
    opp_arg arg1 // n_volume | OPP_READ
) 
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("compute_node_charge_density_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__compute_node_charge_density_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k6_dat0_stride_d, &opp_k6_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k6_dat1_stride_d, &opp_k6_dat1_stride, &(args[1].dat->set->set_capacity), 1);

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

        {
            opp_dev_compute_node_charge_density_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // n_charge_den
                (OPP_REAL *)args[1].data_d,     // n_volume
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("compute_node_charge_density_kernel");
}

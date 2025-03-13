
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k4_dat0_stride = -1;
OPP_INT opp_k4_dat1_stride = -1;

__constant__ OPP_INT opp_k4_dat0_stride_d;
__constant__ OPP_INT opp_k4_dat1_stride_d;



namespace opp_k4 {
__device__ void sum_laplace_kernel(
        const double* node0_xlocal,
        double* node0_field_P
    )
{
    double rv = 0.0;
    double lv = CONST_lhs_voltage_d[0];

    double frac = ((*node0_xlocal) / CONST_L_d[0]);
    (*node0_field_P) += (frac * rv + (1. - frac) * lv);
}

}

//--------------------------------------------------------------
__global__ void opp_dev_sum_laplace_kernel(
    const OPP_REAL *__restrict__ dat0,  // n_xlocal
    OPP_REAL *__restrict__ dat1,  // n_field_p
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    for (int n = thread_id; n < (end - start); n += blockDim.x * gridDim.x) {

        opp_k4::sum_laplace_kernel(
            dat0 + n, // n_xlocal 
            dat1 + n // n_field_p 
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__sum_laplace_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // n_xlocal | OPP_READ
    opp_arg arg1 // n_field_p | OPP_RW
) 
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("sum_laplace_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__sum_laplace_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    if (opp_k4_dat0_stride != args[0].dat->set->set_capacity) {
        opp_k4_dat0_stride = args[0].dat->set->set_capacity;
        cutilSafeCall(cudaMemcpyToSymbol(opp_k4_dat0_stride_d, &opp_k4_dat0_stride, sizeof(OPP_INT)));
    }
    if (opp_k4_dat1_stride != args[1].dat->set->set_capacity) {
        opp_k4_dat1_stride = args[1].dat->set->set_capacity;
        cutilSafeCall(cudaMemcpyToSymbol(opp_k4_dat1_stride_d, &opp_k4_dat1_stride, sizeof(OPP_INT)));
    }

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
            opp_dev_sum_laplace_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // n_xlocal
                (OPP_REAL *)args[1].data_d,     // n_field_p
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());   
 
    opp_profiler->end("sum_laplace_kernel");
}


//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k1_dat0_stride = -1;
OPP_INT opp_k1_dat1_stride = -1;
OPP_INT opp_k1_dat2_stride = -1;

__constant__ OPP_INT opp_k1_dat0_stride_d;
__constant__ OPP_INT opp_k1_dat1_stride_d;
__constant__ OPP_INT opp_k1_dat2_stride_d;



namespace opp_k1 {
__device__ inline void update_pos_kernel(const double* p_vel, double* p_pos, int* p_mdir)
{
    for (int dm = 0; dm < 2; dm++) {

        const double offset = p_vel[(dm) * opp_k1_dat0_stride_d] * CONST_dt_d[0];
        p_pos[(dm) * opp_k1_dat1_stride_d] += offset; // s1 = s0 + ut

        // correct for periodic boundary conditions
        const int n_extent_offset_int = std::abs(p_pos[(dm) * opp_k1_dat1_stride_d]) + 2.0;
        const double temp_pos = p_pos[(dm) * opp_k1_dat1_stride_d] + n_extent_offset_int * CONST_extents_d[dm];
        p_pos[(dm) * opp_k1_dat1_stride_d] = std::fmod(temp_pos, CONST_extents_d[dm]);

        p_mdir[(dm) * opp_k1_dat2_stride_d] = (offset > 0) ? 1 : -1;
    }
}

}

//--------------------------------------------------------------
__global__ void opp_dev_update_pos_kernel(
    const OPP_REAL *__restrict__ dat0,  // p_vel
    OPP_REAL *__restrict__ dat1,  // p_pos
    OPP_INT *__restrict__ dat2,  // p_mdir
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    for (int n = thread_id; n < (end - start); n += blockDim.x * gridDim.x) {

        opp_k1::update_pos_kernel(
            dat0 + n, // p_vel 
            dat1 + n, // p_pos 
            dat2 + n // p_mdir 
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__update_pos_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // p_vel | OPP_READ
    opp_arg arg1, // p_pos | OPP_RW
    opp_arg arg2 // p_mdir | OPP_WRITE
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("update_pos_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__update_pos_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    if (opp_k1_dat0_stride != args[0].dat->set->set_capacity) {
        opp_k1_dat0_stride = args[0].dat->set->set_capacity;
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k1_dat0_stride_d), &opp_k1_dat0_stride, sizeof(OPP_INT)));
    }
    if (opp_k1_dat1_stride != args[1].dat->set->set_capacity) {
        opp_k1_dat1_stride = args[1].dat->set->set_capacity;
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k1_dat1_stride_d), &opp_k1_dat1_stride, sizeof(OPP_INT)));
    }
    if (opp_k1_dat2_stride != args[2].dat->set->set_capacity) {
        opp_k1_dat2_stride = args[2].dat->set->set_capacity;
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k1_dat2_stride_d), &opp_k1_dat2_stride, sizeof(OPP_INT)));
    }

#ifdef OPP_BLOCK_SIZE_1
    const int block_size = OPP_BLOCK_SIZE_1;
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
            opp_dev_update_pos_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // p_vel
                (OPP_REAL *)args[1].data_d,     // p_pos
                (OPP_INT *)args[2].data_d,     // p_mdir
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("update_pos_kernel");
}

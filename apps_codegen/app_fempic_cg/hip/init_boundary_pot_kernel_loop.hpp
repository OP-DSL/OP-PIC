
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k1_dat0_stride = -1;
OPP_INT opp_k1_dat1_stride = -1;

__constant__ OPP_INT opp_k1_dat0_stride_d;
__constant__ OPP_INT opp_k1_dat1_stride_d;



namespace opp_k1 {
__device__ inline void init_boundary_pot_kernel(
    const int *node_type,
    double *n_bnd_pot
) {
    switch (*node_type) {
        case 2: // INLET:
            *n_bnd_pot = 0; break;
        case 3: // FIXED:
            *n_bnd_pot = -1 * CONST_wall_potential_d[0]; break;
        default: // NORMAL or OPEN
            *n_bnd_pot = 0; /*default*/
    }
}
}

//--------------------------------------------------------------
__global__ void opp_dev_init_boundary_pot_kernel(
    const OPP_INT *__restrict__ dat0,  // n_type
    OPP_REAL *__restrict__ dat1,  // n_bnd_pot
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    for (int n = thread_id; n < (end - start); n += blockDim.x * gridDim.x) {

        opp_k1::init_boundary_pot_kernel(
            dat0 + n, // n_type 
            dat1 + n // n_bnd_pot 
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__init_boundary_pot_kernel(opp_set set,
    opp_arg arg0, // n_type | OPP_READ
    opp_arg arg1 // n_bnd_pot | OPP_WRITE
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("init_boundary_pot_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__init_boundary_pot_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k1_dat0_stride_d, &opp_k1_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k1_dat1_stride_d, &opp_k1_dat1_stride, &(args[1].dat->set->set_capacity), 1);

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
            opp_dev_init_boundary_pot_kernel<<<num_blocks, block_size>>>(
                (OPP_INT *)args[0].data_d,     // n_type
                (OPP_REAL *)args[1].data_d,     // n_bnd_pot
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("init_boundary_pot_kernel");
}

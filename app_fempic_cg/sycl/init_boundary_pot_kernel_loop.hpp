
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k1_dat0_stride = -1;
OPP_INT opp_k1_dat1_stride = -1;

OPP_INT* opp_k1_dat0_stride_s = nullptr;
OPP_INT* opp_k1_dat1_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__init_boundary_pot_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // n_type | OPP_READ
    opp_arg arg1 // n_bnd_pot | OPP_WRITE
) 
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("init_boundary_pot_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__init_boundary_pot_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_set_stride(opp_k1_dat0_stride_s, opp_k1_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k1_dat1_stride_s, opp_k1_dat1_stride, args[1].dat->set->set_capacity);

#ifdef OPP_BLOCK_SIZE_1
    const int block_size = OPP_BLOCK_SIZE_1;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    if (iter_size > 0) {

        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k1_dat0_stride_sycl = opp_k1_dat0_stride_s;
            const OPP_INT* opp_k1_dat1_stride_sycl = opp_k1_dat1_stride_s;
    
            const OPP_REAL* CONST_wall_potential_sycl = CONST_wall_potential_s;

            OPP_INT* dat0_sycl = (OPP_INT*)args[0].data_d;     // n_type
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // n_bnd_pot

            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            auto  init_boundary_pot_kernel_sycl = [=](
                const int *node_type,
                double *n_bnd_pot
            )
            {
                switch (*node_type)
                {
                    case 2: // INLET:
                        *n_bnd_pot = 0;
                        break;
                    case 3: // FIXED:
                        *n_bnd_pot = -1 * CONST_wall_potential_sycl[0];
                        break;
                    default: // NORMAL or OPEN
                        *n_bnd_pot = 0; /*default*/
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                for (int n = tid; n < iter_size; n += item.get_global_range()[0]) {

                    init_boundary_pot_kernel_sycl(
                        dat0_sycl + n, // n_type 
                        dat1_sycl + n // n_bnd_pot 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class init_boundary_pot_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    opp_queue->wait();   
 
    opp_profiler->end("init_boundary_pot_kernel");
}

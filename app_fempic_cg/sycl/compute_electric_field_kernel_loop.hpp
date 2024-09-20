
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k7_dat0_stride = -1;
OPP_INT opp_k7_dat1_stride = -1;
OPP_INT opp_k7_dat2_stride = -1;
OPP_INT opp_k7_map0_stride = -1;

OPP_INT* opp_k7_dat0_stride_s = nullptr;
OPP_INT* opp_k7_dat1_stride_s = nullptr;
OPP_INT* opp_k7_dat2_stride_s = nullptr;
OPP_INT* opp_k7_map0_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__compute_electric_field_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // c_ef | OPP_INC
    opp_arg arg1, // c_sd | OPP_READ
    opp_arg arg2, // n_potential | OPP_READ
    opp_arg arg3, // n_potential | OPP_READ
    opp_arg arg4, // n_potential | OPP_READ
    opp_arg arg5 // n_potential | OPP_READ
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

    opp_profiler->start("compute_electric_field_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__compute_electric_field_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_set_stride(opp_k7_dat0_stride_s, opp_k7_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k7_dat1_stride_s, opp_k7_dat1_stride, args[1].dat->set->set_capacity);
    opp_set_stride(opp_k7_dat2_stride_s, opp_k7_dat2_stride, args[2].dat->set->set_capacity);
    opp_set_stride(opp_k7_map0_stride_s, opp_k7_map0_stride, args[2].size);

#ifdef OPP_BLOCK_SIZE_7
    const int block_size = OPP_BLOCK_SIZE_7;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    if (iter_size > 0) {

        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k7_dat0_stride_sycl = opp_k7_dat0_stride_s;
            const OPP_INT* opp_k7_dat1_stride_sycl = opp_k7_dat1_stride_s;
            const OPP_INT* opp_k7_dat2_stride_sycl = opp_k7_dat2_stride_s;
            const OPP_INT* opp_k7_map0_stride_sycl = opp_k7_map0_stride_s;
    

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // c_ef
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // c_sd
            OPP_REAL* dat2_sycl = (OPP_REAL*)args[2].data_d;     // n_potential
            const OPP_INT* map0_sycl = args[2].map_data_d;     // c2n_map

            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            auto  compute_electric_field_kernel_sycl = [=](
                double *cell_electric_field,
                const double *cell_shape_deriv,
                const double *node_potential0,
                const double *node_potential1,
                const double *node_potential2,
                const double *node_potential3
            )
            {
                for (int dim = 0; dim < 3; dim++)
                {
                    const double c1 = (cell_shape_deriv[(0 * 3 + dim) * opp_k7_dat1_stride_sycl[0]] * (*node_potential0));
                    const double c2 = (cell_shape_deriv[(1 * 3 + dim) * opp_k7_dat1_stride_sycl[0]] * (*node_potential1));
                    const double c3 = (cell_shape_deriv[(2 * 3 + dim) * opp_k7_dat1_stride_sycl[0]] * (*node_potential2));
                    const double c4 = (cell_shape_deriv[(3 * 3 + dim) * opp_k7_dat1_stride_sycl[0]] * (*node_potential3));

                    cell_electric_field[(dim) * opp_k7_dat0_stride_sycl[0]] -= (c1 + c2 + c3 + c4);
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                const int n = (tid + start);
                if (n < end) { 

                    compute_electric_field_kernel_sycl(
                        dat0_sycl + n, // c_ef 
                        dat1_sycl + n, // c_sd 
                        dat2_sycl + map0_sycl[n + opp_k7_map0_stride_sycl[0] * 0], // n_potential 
                        dat2_sycl + map0_sycl[n + opp_k7_map0_stride_sycl[0] * 1], // n_potential 
                        dat2_sycl + map0_sycl[n + opp_k7_map0_stride_sycl[0] * 2], // n_potential 
                        dat2_sycl + map0_sycl[n + opp_k7_map0_stride_sycl[0] * 3] // n_potential 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class compute_electric_field_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    opp_queue->wait();   
 
    opp_profiler->end("compute_electric_field_kernel");
}


//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k3_dat0_stride = -1;
OPP_INT opp_k3_dat1_stride = -1;
OPP_INT opp_k3_map0_stride = -1;

OPP_INT* opp_k3_dat0_stride_s = nullptr;
OPP_INT* opp_k3_dat1_stride_s = nullptr;
OPP_INT* opp_k3_map0_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__weight_p2f_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // n_field_j | OPP_INC
    opp_arg arg1, // n_field_j | OPP_INC
    opp_arg arg2 // p_pos_x | OPP_READ
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("weight_p2f_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__weight_p2f_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);

#ifdef USE_MPI
    opp_init_double_indirect_reductions_device(nargs, args);
#endif
 
 
    opp_set_stride(opp_k3_dat0_stride_s, opp_k3_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k3_dat1_stride_s, opp_k3_dat1_stride, args[2].dat->set->set_capacity);
    opp_set_stride(opp_k3_map0_stride_s, opp_k3_map0_stride, args[0].size);

#ifdef OPP_BLOCK_SIZE_3
    const int block_size = OPP_BLOCK_SIZE_3;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    if (iter_size > 0) {

        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k3_dat0_stride_sycl = opp_k3_dat0_stride_s;
            const OPP_INT* opp_k3_dat1_stride_sycl = opp_k3_dat1_stride_s;
            const OPP_INT* opp_k3_map0_stride_sycl = opp_k3_map0_stride_s;
    
            const OPP_REAL* CONST_qscale_sycl = CONST_qscale_s;
            const OPP_REAL* CONST_dx_sycl = CONST_dx_s;
            const OPP_REAL* CONST_xl_sycl = CONST_xl_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // n_field_j
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[2].data_d;     // p_pos_x
            const OPP_INT* map0_sycl = args[0].map_data_d;     // c2n_map

            const OPP_INT *p2c_map = (OPP_INT *)set->mesh_relation_dat->data_d;
            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
             weight_p2f_kernel_sycl = [=](
                    double* node0_field_J,
                    double* node1_field_J,
                    const double* particle0_position_x
                )
            {
                double xx = ((particle0_position_x[(0) * opp_k3_dat1_stride_sycl[0]] - CONST_xl_sycl[0]) / CONST_dx_sycl[0]); // Makes Global position to local position comapared to the cell
                int n = int(xx);
                double frac = (xx - n);

                (*node0_field_J) += (CONST_qscale_sycl[0] * (1.0 - frac));  // Can change qscale to be calculated from particle data
                (*node1_field_J) += (CONST_qscale_sycl[0] * frac);
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                const int n = (tid + start);
                if (n < end) { 
                    const OPP_INT* opp_p2c = p2c_map + n;
                    OPP_REAL arg0_0_local[1];
                    for (int d = 0; d < 1; ++d)
                        arg0_0_local[d] = OPP_REAL_ZERO;
                    OPP_REAL arg1_1_local[1];
                    for (int d = 0; d < 1; ++d)
                        arg1_1_local[d] = OPP_REAL_ZERO;

                    weight_p2f_kernel_sycl(
                        arg0_0_local, // n_field_j 
                        arg1_1_local, // n_field_j 
                        dat1_sycl + n // p_pos_x 
                    );

                    for (int d = 0; d < 1; ++d)
                        opp_atomic_fetch_add(dat0_sycl + map0_sycl[opp_k3_map0_stride_sycl[0] * 0 + opp_p2c[0]] + (d * opp_k3_dat0_stride_sycl[0]), arg0_0_local[d]);

                    for (int d = 0; d < 1; ++d)
                        opp_atomic_fetch_add(dat0_sycl + map0_sycl[opp_k3_map0_stride_sycl[0] * 1 + opp_p2c[0]] + (d * opp_k3_dat0_stride_sycl[0]), arg1_1_local[d]);
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class weight_p2f_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   

#ifdef USE_MPI    
    opp_exchange_double_indirect_reductions_device(nargs, args);
    opp_complete_double_indirect_reductions_device(nargs, args);
#endif
 
    opp_profiler->end("weight_p2f_kernel");
}

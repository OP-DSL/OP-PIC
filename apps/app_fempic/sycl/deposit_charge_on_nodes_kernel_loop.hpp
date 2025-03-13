
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k5_dat0_stride = -1;
OPP_INT opp_k5_dat1_stride = -1;
OPP_INT opp_k5_map0_stride = -1;

OPP_INT* opp_k5_dat0_stride_s = nullptr;
OPP_INT* opp_k5_dat1_stride_s = nullptr;
OPP_INT* opp_k5_map0_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__deposit_charge_on_nodes_kernel(opp_set set,
    opp_arg arg0, // p_lc | OPP_READ
    opp_arg arg1, // n_charge_den | OPP_INC
    opp_arg arg2, // n_charge_den | OPP_INC
    opp_arg arg3, // n_charge_den | OPP_INC
    opp_arg arg4 // n_charge_den | OPP_INC
) 
{
    const int nargs = 5;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;

    opp_profiler->start("deposit_charge_on_nodes_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__deposit_charge_on_nodes_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);

#ifdef USE_MPI
    opp_init_double_indirect_reductions_device(nargs, args);
#endif
 
 
    opp_set_stride(opp_k5_dat0_stride_s, opp_k5_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k5_dat1_stride_s, opp_k5_dat1_stride, args[1].dat->set->set_capacity);
    opp_set_stride(opp_k5_map0_stride_s, opp_k5_map0_stride, args[1].size);

#ifdef OPP_BLOCK_SIZE_5
    const int block_size = OPP_BLOCK_SIZE_5;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    if (iter_size > 0) {

        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k5_dat0_stride_sycl = opp_k5_dat0_stride_s;
            const OPP_INT* opp_k5_dat1_stride_sycl = opp_k5_dat1_stride_s;
            const OPP_INT* opp_k5_map0_stride_sycl = opp_k5_map0_stride_s;
    

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // p_lc
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // n_charge_den
            const OPP_INT* map0_sycl = args[1].map_data_d;     // c2n_map

            const OPP_INT *p2c_map = (OPP_INT *)set->mesh_relation_dat->data_d;
            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            auto  deposit_charge_on_nodes_kernel_sycl = [=](
                const double *part_lc,
                double *node_charge_den0,
                double *node_charge_den1,
                double *node_charge_den2,
                double *node_charge_den3
            ) {
                node_charge_den0[0] += part_lc[(0) * opp_k5_dat0_stride_sycl[0]];
                node_charge_den1[0] += part_lc[(1) * opp_k5_dat0_stride_sycl[0]];
                node_charge_den2[0] += part_lc[(2) * opp_k5_dat0_stride_sycl[0]];
                node_charge_den3[0] += part_lc[(3) * opp_k5_dat0_stride_sycl[0]];
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                const int n = (tid + start);
                if (n < end) { 
                    const OPP_INT* opp_p2c = p2c_map + n;
                    OPP_REAL arg1_0_local[1];
                    for (int d = 0; d < 1; ++d)
                        arg1_0_local[d] = OPP_REAL_ZERO;
                    OPP_REAL arg2_1_local[1];
                    for (int d = 0; d < 1; ++d)
                        arg2_1_local[d] = OPP_REAL_ZERO;
                    OPP_REAL arg3_2_local[1];
                    for (int d = 0; d < 1; ++d)
                        arg3_2_local[d] = OPP_REAL_ZERO;
                    OPP_REAL arg4_3_local[1];
                    for (int d = 0; d < 1; ++d)
                        arg4_3_local[d] = OPP_REAL_ZERO;

                    deposit_charge_on_nodes_kernel_sycl(
                        dat0_sycl + n, // p_lc 
                        arg1_0_local, // n_charge_den 
                        arg2_1_local, // n_charge_den 
                        arg3_2_local, // n_charge_den 
                        arg4_3_local // n_charge_den 
                    );

                    for (int d = 0; d < 1; ++d)
                        opp_atomic_fetch_add(dat1_sycl + map0_sycl[opp_k5_map0_stride_sycl[0] * 0 + opp_p2c[0]] + (d * opp_k5_dat1_stride_sycl[0]), arg1_0_local[d]);

                    for (int d = 0; d < 1; ++d)
                        opp_atomic_fetch_add(dat1_sycl + map0_sycl[opp_k5_map0_stride_sycl[0] * 1 + opp_p2c[0]] + (d * opp_k5_dat1_stride_sycl[0]), arg2_1_local[d]);

                    for (int d = 0; d < 1; ++d)
                        opp_atomic_fetch_add(dat1_sycl + map0_sycl[opp_k5_map0_stride_sycl[0] * 2 + opp_p2c[0]] + (d * opp_k5_dat1_stride_sycl[0]), arg3_2_local[d]);

                    for (int d = 0; d < 1; ++d)
                        opp_atomic_fetch_add(dat1_sycl + map0_sycl[opp_k5_map0_stride_sycl[0] * 3 + opp_p2c[0]] + (d * opp_k5_dat1_stride_sycl[0]), arg4_3_local[d]);
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class deposit_charge_on_nodes_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   

#ifdef USE_MPI    
    opp_exchange_double_indirect_reductions_device(nargs, args);
    opp_complete_double_indirect_reductions_device(nargs, args);
#endif
 
    opp_profiler->end("deposit_charge_on_nodes_kernel");
}

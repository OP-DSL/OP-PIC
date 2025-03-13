
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k9_dat0_stride = -1;
OPP_INT opp_k9_dat1_stride = -1;
OPP_INT opp_k9_dat2_stride = -1;

OPP_INT* opp_k9_dat0_stride_s = nullptr;
OPP_INT* opp_k9_dat1_stride_s = nullptr;
OPP_INT* opp_k9_dat2_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__get_max_x_values_kernel(opp_set set,
    opp_arg arg0, // c_j | OPP_READ
    opp_arg arg1, // | OPP_MAX
    opp_arg arg2, // c_e | OPP_READ
    opp_arg arg3, // | OPP_MAX
    opp_arg arg4, // c_b | OPP_READ
    opp_arg arg5 // | OPP_MAX
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

    opp_profiler->start("get_max_x_values_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__get_max_x_values_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    OPP_REAL *arg1_host_data = (OPP_REAL *)args[1].data;
    OPP_REAL *arg3_host_data = (OPP_REAL *)args[3].data;
    OPP_REAL *arg5_host_data = (OPP_REAL *)args[5].data;

    opp_set_stride(opp_k9_dat0_stride_s, opp_k9_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k9_dat1_stride_s, opp_k9_dat1_stride, args[2].dat->set->set_capacity);
    opp_set_stride(opp_k9_dat2_stride_s, opp_k9_dat2_stride, args[4].dat->set->set_capacity);

#ifdef OPP_BLOCK_SIZE_9
    const int block_size = OPP_BLOCK_SIZE_9;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    const int num_blocks = 200;

    int reduction_bytes = 0;
    int reduction_size = 0;

    reduction_bytes += ROUND_UP(num_blocks * 1 * sizeof(OPP_REAL));
    reduction_size   = MAX(reduction_size, sizeof(OPP_REAL));
    reduction_bytes += ROUND_UP(num_blocks * 1 * sizeof(OPP_REAL));
    reduction_size   = MAX(reduction_size, sizeof(OPP_REAL));
    reduction_bytes += ROUND_UP(num_blocks * 1 * sizeof(OPP_REAL));
    reduction_size   = MAX(reduction_size, sizeof(OPP_REAL));

    opp_reallocReductArrays(reduction_bytes);
    reduction_bytes = 0;

    args[1].data   = OPP_reduct_h + reduction_bytes;
    args[1].data_d = OPP_reduct_d + reduction_bytes;

    for (int b = 0; b < num_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            ((OPP_REAL *)args[1].data)[b * 1 + d] = arg1_host_data[d];
    }

    reduction_bytes += ROUND_UP(num_blocks * 1 * sizeof(OPP_REAL));
    args[3].data   = OPP_reduct_h + reduction_bytes;
    args[3].data_d = OPP_reduct_d + reduction_bytes;

    for (int b = 0; b < num_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            ((OPP_REAL *)args[3].data)[b * 1 + d] = arg3_host_data[d];
    }

    reduction_bytes += ROUND_UP(num_blocks * 1 * sizeof(OPP_REAL));
    args[5].data   = OPP_reduct_h + reduction_bytes;
    args[5].data_d = OPP_reduct_d + reduction_bytes;

    for (int b = 0; b < num_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            ((OPP_REAL *)args[5].data)[b * 1 + d] = arg5_host_data[d];
    }

    reduction_bytes += ROUND_UP(num_blocks * 1 * sizeof(OPP_REAL));

    opp_mvReductArraysToDevice(reduction_bytes);
    
    if (iter_size > 0) {
        
        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k9_dat0_stride_sycl = opp_k9_dat0_stride_s;
            const OPP_INT* opp_k9_dat1_stride_sycl = opp_k9_dat1_stride_s;
            const OPP_INT* opp_k9_dat2_stride_sycl = opp_k9_dat2_stride_s;
    

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // c_j
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[2].data_d;     // c_e
            OPP_REAL* dat2_sycl = (OPP_REAL*)args[4].data_d;     // c_b
            OPP_REAL* gbl1_sycl = (OPP_REAL*)args[1].data_d;
            sycl::accessor<OPP_REAL, 1, sycl::access::mode::read_write, sycl::access::target::local>
                                        red_OPP_REAL_1(block_size, cgh); // temp var for reduction
            OPP_REAL* gbl3_sycl = (OPP_REAL*)args[3].data_d;
            sycl::accessor<OPP_REAL, 1, sycl::access::mode::read_write, sycl::access::target::local>
                                        red_OPP_REAL_3(block_size, cgh); // temp var for reduction
            OPP_REAL* gbl5_sycl = (OPP_REAL*)args[5].data_d;
            sycl::accessor<OPP_REAL, 1, sycl::access::mode::read_write, sycl::access::target::local>
                                        red_OPP_REAL_5(block_size, cgh); // temp var for reduction

            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            auto  get_max_x_values_kernel_sycl = [=](
                const double* cell_j,
                double* max_j,
                const double* cell_e,
                double* max_e,
                const double* cell_b,
                double* max_b)
            {
                *max_j = ((*cell_j > *max_j) ? (*cell_j) : (*max_j));

                *max_e = ((*cell_e > *max_e) ? (*cell_e) : (*max_e));

                *max_b = ((*cell_b > *max_b) ? (*cell_b) : (*max_b));
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                OPP_REAL gbl1_local[1];
                for (int d = 0; d < 1; ++d)
                    gbl1_local[d] = gbl1_sycl[item.get_group(0) * 1 + d];

                OPP_REAL gbl3_local[1];
                for (int d = 0; d < 1; ++d)
                    gbl3_local[d] = gbl3_sycl[item.get_group(0) * 1 + d];

                OPP_REAL gbl5_local[1];
                for (int d = 0; d < 1; ++d)
                    gbl5_local[d] = gbl5_sycl[item.get_group(0) * 1 + d];

                const int tid = item.get_global_linear_id();
                for (int n = tid; n < iter_size; n += item.get_global_range()[0]) {

                    get_max_x_values_kernel_sycl(
                        dat0_sycl + n, // c_j 
                        gbl1_local, // 
                        dat1_sycl + n, // c_e 
                        gbl3_local, // 
                        dat2_sycl + n, // c_b 
                        gbl5_local // 
                    );
                }

                for (int d = 0; d < 1; ++d) //arg1_offset + 
                    opp_reduction<OPP_MAX, 0>(gbl1_sycl, (d + item.get_group_linear_id() * 1), 
                                    gbl1_local[d], red_OPP_REAL_1, item);

                for (int d = 0; d < 1; ++d) //arg3_offset + 
                    opp_reduction<OPP_MAX, 0>(gbl3_sycl, (d + item.get_group_linear_id() * 1), 
                                    gbl3_local[d], red_OPP_REAL_3, item);

                for (int d = 0; d < 1; ++d) //arg5_offset + 
                    opp_reduction<OPP_MAX, 0>(gbl5_sycl, (d + item.get_group_linear_id() * 1), 
                                    gbl5_local[d], red_OPP_REAL_5, item);
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class get_max_x_values_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_mvReductArraysToHost(reduction_bytes);

    for (int b = 0; b < num_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            arg1_host_data[d] = MAX(arg1_host_data[d], ((OPP_REAL *)args[1].data)[b * 1 + d]);
    }
    for (int b = 0; b < num_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            arg3_host_data[d] = MAX(arg3_host_data[d], ((OPP_REAL *)args[3].data)[b * 1 + d]);
    }
    for (int b = 0; b < num_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            arg5_host_data[d] = MAX(arg5_host_data[d], ((OPP_REAL *)args[5].data)[b * 1 + d]);
    }

    args[1].data = (char *)arg1_host_data;
    opp_mpi_reduce(&args[1], arg1_host_data);

    args[3].data = (char *)arg3_host_data;
    opp_mpi_reduce(&args[3], arg3_host_data);

    args[5].data = (char *)arg5_host_data;
    opp_mpi_reduce(&args[5], arg5_host_data);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("get_max_x_values_kernel");
}

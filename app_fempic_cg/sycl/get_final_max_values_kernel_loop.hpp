
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k9_dat0_stride = -1;
OPP_INT opp_k9_dat1_stride = -1;

OPP_INT* opp_k9_dat0_stride_s = nullptr;
OPP_INT* opp_k9_dat1_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__get_final_max_values_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // n_charge_den | OPP_READ
    opp_arg arg1, // | OPP_MAX
    opp_arg arg2, // n_potential | OPP_READ
    opp_arg arg3 // | OPP_MAX
) 
{
    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    opp_profiler->start("get_final_max_values_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__get_final_max_values_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    OPP_REAL *arg1_host_data = (OPP_REAL *)args[1].data;
    OPP_REAL *arg3_host_data = (OPP_REAL *)args[3].data;

    opp_set_stride(opp_k9_dat0_stride_s, opp_k9_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k9_dat1_stride_s, opp_k9_dat1_stride, args[2].dat->set->set_capacity);

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

    opp_mvReductArraysToDevice(reduction_bytes);
    
    if (iter_size > 0) {

        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k9_dat0_stride_sycl = opp_k9_dat0_stride_s;
            const OPP_INT* opp_k9_dat1_stride_sycl = opp_k9_dat1_stride_s;
    

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // n_charge_den
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[2].data_d;     // n_potential
            OPP_REAL* gbl1_sycl = (OPP_REAL*)args[1].data_d;
            sycl::accessor<OPP_REAL, 1, sycl::access::mode::read_write, sycl::access::target::local>
                                        red_OPP_REAL_1(block_size, cgh); // temp var for reduction
            OPP_REAL* gbl3_sycl = (OPP_REAL*)args[3].data_d;
            sycl::accessor<OPP_REAL, 1, sycl::access::mode::read_write, sycl::access::target::local>
                                        red_OPP_REAL_3(block_size, cgh); // temp var for reduction

            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            auto  get_final_max_values_kernel_sycl = [=](
                const double* n_charge_den,
                double* max_n_charge_den,
                const double* n_pot,
                double* max_n_pot)
            {
                *max_n_charge_den = ((abs(*n_charge_den) > *max_n_charge_den) ? (abs(*n_charge_den)) : (*max_n_charge_den));

                *max_n_pot = ((*n_pot > *max_n_pot) ? (*n_pot) : (*max_n_pot));
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                OPP_REAL gbl1_local[1];
                for (int d = 0; d < 1; ++d)
                    gbl1_local[d] = gbl1_sycl[item.get_group(0) * 1 + d];

                OPP_REAL gbl3_local[1];
                for (int d = 0; d < 1; ++d)
                    gbl3_local[d] = gbl3_sycl[item.get_group(0) * 1 + d];

                const int tid = item.get_global_linear_id();
                for (int n = tid; n < iter_size; n += item.get_global_range()[0]) {

                    get_final_max_values_kernel_sycl(
                        dat0_sycl + n, // n_charge_den 
                        gbl1_local, // 
                        dat1_sycl + n, // n_potential 
                        gbl3_local // 
                    );
                }

                for (int d = 0; d < 1; ++d) //arg1_offset + 
                    opp_reduction<OPP_MAX, 0>(gbl1_sycl, (d + item.get_group_linear_id() * 1), 
                                    gbl1_local[d], red_OPP_REAL_1, item);

                for (int d = 0; d < 1; ++d) //arg3_offset + 
                    opp_reduction<OPP_MAX, 0>(gbl3_sycl, (d + item.get_group_linear_id() * 1), 
                                    gbl3_local[d], red_OPP_REAL_3, item);
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class get_final_max_values_kernel>(
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

    args[1].data = (char *)arg1_host_data;
    opp_mpi_reduce(&args[1], arg1_host_data);

    args[3].data = (char *)arg3_host_data;
    opp_mpi_reduce(&args[3], arg3_host_data);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    opp_queue->wait();   
 
    opp_profiler->end("get_final_max_values_kernel");
}


//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k8_dat0_stride = -1;
OPP_INT opp_k8_dat1_stride = -1;

OPP_INT* opp_k8_dat0_stride_s = nullptr;
OPP_INT* opp_k8_dat1_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__compute_energy_kernel(opp_set set,
    opp_arg arg0, // c_mask_ghost | OPP_READ
    opp_arg arg1, // c_e | OPP_READ
    opp_arg arg2 // | OPP_INC
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("compute_energy_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__compute_energy_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    OPP_REAL *arg2_host_data = (OPP_REAL *)args[2].data;

    opp_set_stride(opp_k8_dat0_stride_s, opp_k8_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k8_dat1_stride_s, opp_k8_dat1_stride, args[1].dat->set->set_capacity);

#ifdef OPP_BLOCK_SIZE_8
    const int block_size = OPP_BLOCK_SIZE_8;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    const int num_blocks = 200;

    int reduction_bytes = 0;
    int reduction_size = 0;

    reduction_bytes += ROUND_UP(num_blocks * 1 * sizeof(OPP_REAL));
    reduction_size   = MAX(reduction_size, sizeof(OPP_REAL));

    opp_reallocReductArrays(reduction_bytes);
    reduction_bytes = 0;

    args[2].data   = OPP_reduct_h + reduction_bytes;
    args[2].data_d = OPP_reduct_d + reduction_bytes;

    for (int b = 0; b < num_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            ((OPP_REAL *)args[2].data)[b * 1 + d] = OPP_REAL_ZERO;
    }

    reduction_bytes += ROUND_UP(num_blocks * 1 * sizeof(OPP_REAL));

    opp_mvReductArraysToDevice(reduction_bytes);
    
    if (iter_size > 0) {
        
        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k8_dat0_stride_sycl = opp_k8_dat0_stride_s;
            const OPP_INT* opp_k8_dat1_stride_sycl = opp_k8_dat1_stride_s;
    

            OPP_INT* dat0_sycl = (OPP_INT*)args[0].data_d;     // c_mask_ghost
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // c_e
            OPP_REAL* gbl2_sycl = (OPP_REAL*)args[2].data_d;
            sycl::accessor<OPP_REAL, 1, sycl::access::mode::read_write, sycl::access::target::local>
                                        red_OPP_REAL_2(block_size, cgh); // temp var for reduction

            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            enum Dim {
                x = 0,
                y = 1,
                z = 2,
            };

            auto  compute_energy_kernel_sycl = [=](
                const int* cell0_ghost,
                const double* cell_field,
                double* energy)
            {
                if (cell0_ghost[(0) * opp_k8_dat0_stride_sycl[0]] == 0)
                {
                    energy[0] += cell_field[(Dim::x) * opp_k8_dat1_stride_sycl[0]] * cell_field[(Dim::x) * opp_k8_dat1_stride_sycl[0]] +
                                cell_field[(Dim::y) * opp_k8_dat1_stride_sycl[0]] * cell_field[(Dim::y) * opp_k8_dat1_stride_sycl[0]] +
                                cell_field[(Dim::z) * opp_k8_dat1_stride_sycl[0]] * cell_field[(Dim::z) * opp_k8_dat1_stride_sycl[0]];
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                OPP_REAL gbl2_local[1];
                for (int d = 0; d < 1; ++d)
                    gbl2_local[d] = OPP_REAL_ZERO;

                const int tid = item.get_global_linear_id();
                for (int n = tid; n < iter_size; n += item.get_global_range()[0]) {

                    compute_energy_kernel_sycl(
                        dat0_sycl + n, // c_mask_ghost 
                        dat1_sycl + n, // c_e 
                        gbl2_local // 
                    );
                }

                for (int d = 0; d < 1; ++d) //arg2_offset + 
                    opp_reduction<OPP_INC, 0>(gbl2_sycl, (d + item.get_group_linear_id() * 1), 
                                    gbl2_local[d], red_OPP_REAL_2, item);
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class compute_energy_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_mvReductArraysToHost(reduction_bytes);

    for (int b = 0; b < num_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            arg2_host_data[d] += ((OPP_REAL *)args[2].data)[b * 1 + d];
    }

    args[2].data = (char *)arg2_host_data;
    opp_mpi_reduce(&args[2], arg2_host_data);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("compute_energy_kernel");
}

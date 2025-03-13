
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
void opp_par_loop_all__update_ghosts_B_kernel(opp_set set,
    opp_arg arg0, // c_mask_ugb | OPP_READ
    opp_arg arg1, // c_b | OPP_READ
    opp_arg arg2, // c_b | OPP_WRITE
    opp_arg arg3 // | OPP_READ
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    opp_profiler->start("update_ghosts_B_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__update_ghosts_B_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    OPP_INT *arg3_host_data = (OPP_INT *)args[3].data;

    int const_bytes = 0;

    const_bytes += ROUND_UP(1 * sizeof(OPP_INT));

    opp_reallocConstArrays(const_bytes);
    const_bytes = 0;

    args[3].data   = OPP_consts_h + const_bytes;
    args[3].data_d = OPP_consts_d + const_bytes;

    for (int d = 0; d < 1; ++d)
        ((OPP_INT *)args[3].data)[d] = arg3_host_data[d];

    const_bytes += ROUND_UP(1 * sizeof(OPP_INT));

    opp_mvConstArraysToDevice(const_bytes);
    opp_set_stride(opp_k5_dat0_stride_s, opp_k5_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k5_dat1_stride_s, opp_k5_dat1_stride, args[1].dat->set->set_capacity);
    opp_set_stride(opp_k5_map0_stride_s, opp_k5_map0_stride, args[2].size);

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
    

            OPP_INT* dat0_sycl = (OPP_INT*)args[0].data_d;     // c_mask_ugb
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // c_b
            const OPP_INT* map0_sycl = args[2].map_data_d;     // c2cugb0_map
            OPP_INT* gbl3_sycl = (OPP_INT*)args[3].data_d;

            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            enum Dim {
                x = 0,
                y = 1,
                z = 2,
            };

            auto  update_ghosts_B_kernel_sycl = [=](
                const int* c_mask_ugb,
                const double* from_cell,
                double* to_cell,
                const int* m_idx)
            {
                if (c_mask_ugb[(*m_idx) * opp_k5_dat0_stride_sycl[0]] == 1)
                {
                    to_cell[(Dim::x) * opp_k5_dat1_stride_sycl[0]] = from_cell[(Dim::x) * opp_k5_dat1_stride_sycl[0]];
                    to_cell[(Dim::y) * opp_k5_dat1_stride_sycl[0]] = from_cell[(Dim::y) * opp_k5_dat1_stride_sycl[0]];
                    to_cell[(Dim::z) * opp_k5_dat1_stride_sycl[0]] = from_cell[(Dim::z) * opp_k5_dat1_stride_sycl[0]];
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                const int n = (tid + start);
                if (n < end) { 

                    update_ghosts_B_kernel_sycl(
                        dat0_sycl + n, // c_mask_ugb 
                        dat1_sycl + n, // c_b 
                        dat1_sycl + map0_sycl[n + opp_k5_map0_stride_sycl[0] * 0], // c_b 
                        gbl3_sycl // 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class update_ghosts_B_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }
    args[3].data = (char *)arg3_host_data;


    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("update_ghosts_B_kernel");
}

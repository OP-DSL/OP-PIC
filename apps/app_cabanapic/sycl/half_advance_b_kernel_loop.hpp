
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k4_dat0_stride = -1;
OPP_INT opp_k4_dat1_stride = -1;
OPP_INT opp_k4_dat2_stride = -1;
OPP_INT opp_k4_map0_stride = -1;

OPP_INT* opp_k4_dat0_stride_s = nullptr;
OPP_INT* opp_k4_dat1_stride_s = nullptr;
OPP_INT* opp_k4_dat2_stride_s = nullptr;
OPP_INT* opp_k4_map0_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__half_advance_b_kernel(opp_set set,
    opp_arg arg0, // c_e | OPP_READ
    opp_arg arg1, // c_e | OPP_READ
    opp_arg arg2, // c_e | OPP_READ
    opp_arg arg3, // c_e | OPP_READ
    opp_arg arg4, // c_b | OPP_INC
    opp_arg arg5 // c_mask_ghost | OPP_READ
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

    opp_profiler->start("half_advance_b_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__half_advance_b_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_set_stride(opp_k4_dat0_stride_s, opp_k4_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k4_dat1_stride_s, opp_k4_dat1_stride, args[4].dat->set->set_capacity);
    opp_set_stride(opp_k4_dat2_stride_s, opp_k4_dat2_stride, args[5].dat->set->set_capacity);
    opp_set_stride(opp_k4_map0_stride_s, opp_k4_map0_stride, args[0].size);

#ifdef OPP_BLOCK_SIZE_4
    const int block_size = OPP_BLOCK_SIZE_4;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    if (iter_size > 0) {

        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k4_dat0_stride_sycl = opp_k4_dat0_stride_s;
            const OPP_INT* opp_k4_dat1_stride_sycl = opp_k4_dat1_stride_s;
            const OPP_INT* opp_k4_dat2_stride_sycl = opp_k4_dat2_stride_s;
            const OPP_INT* opp_k4_map0_stride_sycl = opp_k4_map0_stride_s;
    
            const OPP_REAL* CONST_p_sycl = CONST_p_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // c_e
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[4].data_d;     // c_b
            OPP_INT* dat2_sycl = (OPP_INT*)args[5].data_d;     // c_mask_ghost
            const OPP_INT* map0_sycl = args[0].map_data_d;     // c2c_map

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

            auto  half_advance_b_kernel_sycl = [=] (
                const double* cell_x_e,
                const double* cell_y_e,
                const double* cell_z_e,
                const double* cell0_e,
                double* cell0_b,
                const int* cell0_ghost)
            {
                if (cell0_ghost[(0) * opp_k4_dat2_stride_sycl[0]] == 0)
                {
                    cell0_b[(Dim::x) * opp_k4_dat1_stride_sycl[0]] -= (0.5 * CONST_p_sycl[Dim::y] * (cell_y_e[(Dim::z) * opp_k4_dat0_stride_sycl[0]] - cell0_e[(Dim::z) * opp_k4_dat0_stride_sycl[0]])
                                        - 0.5 * CONST_p_sycl[Dim::z] * (cell_z_e[(Dim::y) * opp_k4_dat0_stride_sycl[0]] - cell0_e[(Dim::y) * opp_k4_dat0_stride_sycl[0]]));

                    cell0_b[(Dim::y) * opp_k4_dat1_stride_sycl[0]] -= (0.5 * CONST_p_sycl[Dim::z] * (cell_z_e[(Dim::x) * opp_k4_dat0_stride_sycl[0]] - cell0_e[(Dim::x) * opp_k4_dat0_stride_sycl[0]])
                                        - 0.5 * CONST_p_sycl[Dim::x] * (cell_x_e[(Dim::z) * opp_k4_dat0_stride_sycl[0]] - cell0_e[(Dim::z) * opp_k4_dat0_stride_sycl[0]]));

                    cell0_b[(Dim::z) * opp_k4_dat1_stride_sycl[0]] -= (0.5 * CONST_p_sycl[Dim::x] * (cell_x_e[(Dim::y) * opp_k4_dat0_stride_sycl[0]] - cell0_e[(Dim::y) * opp_k4_dat0_stride_sycl[0]])
                                        - 0.5 * CONST_p_sycl[Dim::y] * (cell_y_e[(Dim::x) * opp_k4_dat0_stride_sycl[0]] - cell0_e[(Dim::x) * opp_k4_dat0_stride_sycl[0]]));
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                const int n = (tid + start);
                if (n < end) { 

                    half_advance_b_kernel_sycl(
                        dat0_sycl + map0_sycl[n + opp_k4_map0_stride_sycl[0] * 9], // c_e 
                        dat0_sycl + map0_sycl[n + opp_k4_map0_stride_sycl[0] * 7], // c_e 
                        dat0_sycl + map0_sycl[n + opp_k4_map0_stride_sycl[0] * 6], // c_e 
                        dat0_sycl + n, // c_e 
                        dat1_sycl + n, // c_b 
                        dat2_sycl + n // c_mask_ghost 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class half_advance_b_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("half_advance_b_kernel");
}

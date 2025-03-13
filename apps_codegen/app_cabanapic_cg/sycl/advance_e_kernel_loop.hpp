
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k7_dat0_stride = -1;
OPP_INT opp_k7_dat1_stride = -1;
OPP_INT opp_k7_dat2_stride = -1;
OPP_INT opp_k7_dat3_stride = -1;
OPP_INT opp_k7_map0_stride = -1;

OPP_INT* opp_k7_dat0_stride_s = nullptr;
OPP_INT* opp_k7_dat1_stride_s = nullptr;
OPP_INT* opp_k7_dat2_stride_s = nullptr;
OPP_INT* opp_k7_dat3_stride_s = nullptr;
OPP_INT* opp_k7_map0_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__advance_e_kernel(opp_set set,
    opp_arg arg0, // c_b | OPP_READ
    opp_arg arg1, // c_b | OPP_READ
    opp_arg arg2, // c_b | OPP_READ
    opp_arg arg3, // c_b | OPP_READ
    opp_arg arg4, // c_j | OPP_READ
    opp_arg arg5, // c_e | OPP_INC
    opp_arg arg6 // c_mask_right | OPP_READ
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const int nargs = 7;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;
    args[6] = arg6;

    opp_profiler->start("advance_e_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__advance_e_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_set_stride(opp_k7_dat0_stride_s, opp_k7_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k7_dat1_stride_s, opp_k7_dat1_stride, args[4].dat->set->set_capacity);
    opp_set_stride(opp_k7_dat2_stride_s, opp_k7_dat2_stride, args[5].dat->set->set_capacity);
    opp_set_stride(opp_k7_dat3_stride_s, opp_k7_dat3_stride, args[6].dat->set->set_capacity);
    opp_set_stride(opp_k7_map0_stride_s, opp_k7_map0_stride, args[0].size);

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
            const OPP_INT* opp_k7_dat3_stride_sycl = opp_k7_dat3_stride_s;
            const OPP_INT* opp_k7_map0_stride_sycl = opp_k7_map0_stride_s;
    
            const OPP_REAL* CONST_dt_eps0_sycl = CONST_dt_eps0_s;
            const OPP_REAL* CONST_p_sycl = CONST_p_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // c_b
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[4].data_d;     // c_j
            OPP_REAL* dat2_sycl = (OPP_REAL*)args[5].data_d;     // c_e
            OPP_INT* dat3_sycl = (OPP_INT*)args[6].data_d;     // c_mask_right
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

            auto  advance_e_kernel_sycl = [=] (
                const double* cell_x_b,
                const double* cell_y_b,
                const double* cell_z_b,
                const double* cell0_b,
                const double* cell0_j,
                double* cell0_e,
                const int* iter_adv_e)
            {
                if (iter_adv_e[(0) * opp_k7_dat3_stride_sycl[0]] == 1)
                {
                    cell0_e[(Dim::x) * opp_k7_dat2_stride_sycl[0]] += ( - CONST_dt_eps0_sycl[0] * cell0_j[(Dim::x) * opp_k7_dat1_stride_sycl[0]] ) +
                        ( CONST_p_sycl[Dim::y] * (cell0_b[(Dim::z) * opp_k7_dat0_stride_sycl[0]] - cell_y_b[(Dim::z) * opp_k7_dat0_stride_sycl[0]]) -
                        CONST_p_sycl[Dim::z] * (cell0_b[(Dim::y) * opp_k7_dat0_stride_sycl[0]] - cell_z_b[(Dim::y) * opp_k7_dat0_stride_sycl[0]]) );

                    cell0_e[(Dim::y) * opp_k7_dat2_stride_sycl[0]] += ( - CONST_dt_eps0_sycl[0] * cell0_j[(Dim::y) * opp_k7_dat1_stride_sycl[0]] ) +
                        ( CONST_p_sycl[Dim::z] * (cell0_b[(Dim::x) * opp_k7_dat0_stride_sycl[0]] - cell_z_b[(Dim::x) * opp_k7_dat0_stride_sycl[0]]) -
                        CONST_p_sycl[Dim::x] * (cell0_b[(Dim::z) * opp_k7_dat0_stride_sycl[0]] - cell_x_b[(Dim::z) * opp_k7_dat0_stride_sycl[0]]) );

                    cell0_e[(Dim::z) * opp_k7_dat2_stride_sycl[0]] += ( - CONST_dt_eps0_sycl[0] * cell0_j[(Dim::z) * opp_k7_dat1_stride_sycl[0]] ) +
                        ( CONST_p_sycl[Dim::x] * (cell0_b[(Dim::y) * opp_k7_dat0_stride_sycl[0]] - cell_x_b[(Dim::y) * opp_k7_dat0_stride_sycl[0]]) -
                        CONST_p_sycl[Dim::y] * (cell0_b[(Dim::x) * opp_k7_dat0_stride_sycl[0]] - cell_y_b[(Dim::x) * opp_k7_dat0_stride_sycl[0]]) );
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                const int n = (tid + start);
                if (n < end) { 

                    advance_e_kernel_sycl(
                        dat0_sycl + map0_sycl[n + opp_k7_map0_stride_sycl[0] * 2], // c_b 
                        dat0_sycl + map0_sycl[n + opp_k7_map0_stride_sycl[0] * 4], // c_b 
                        dat0_sycl + map0_sycl[n + opp_k7_map0_stride_sycl[0] * 5], // c_b 
                        dat0_sycl + n, // c_b 
                        dat1_sycl + n, // c_j 
                        dat2_sycl + n, // c_e 
                        dat3_sycl + n // c_mask_right 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class advance_e_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("advance_e_kernel");
}

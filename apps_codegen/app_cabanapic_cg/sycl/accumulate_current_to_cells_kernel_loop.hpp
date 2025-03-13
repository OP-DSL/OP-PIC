
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k3_dat0_stride = -1;
OPP_INT opp_k3_dat1_stride = -1;
OPP_INT opp_k3_dat2_stride = -1;
OPP_INT opp_k3_map0_stride = -1;

OPP_INT* opp_k3_dat0_stride_s = nullptr;
OPP_INT* opp_k3_dat1_stride_s = nullptr;
OPP_INT* opp_k3_dat2_stride_s = nullptr;
OPP_INT* opp_k3_map0_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__accumulate_current_to_cells_kernel(opp_set set,
    opp_arg arg0, // c_j | OPP_WRITE
    opp_arg arg1, // c_acc | OPP_READ
    opp_arg arg2, // c_acc | OPP_READ
    opp_arg arg3, // c_acc | OPP_READ
    opp_arg arg4, // c_acc | OPP_READ
    opp_arg arg5, // c_acc | OPP_READ
    opp_arg arg6, // c_acc | OPP_READ
    opp_arg arg7, // c_acc | OPP_READ
    opp_arg arg8 // c_mask_right | OPP_READ
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const int nargs = 9;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;
    args[6] = arg6;
    args[7] = arg7;
    args[8] = arg8;

    opp_profiler->start("accumulate_current_to_cells_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__accumulate_current_to_cells_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_set_stride(opp_k3_dat0_stride_s, opp_k3_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k3_dat1_stride_s, opp_k3_dat1_stride, args[1].dat->set->set_capacity);
    opp_set_stride(opp_k3_dat2_stride_s, opp_k3_dat2_stride, args[8].dat->set->set_capacity);
    opp_set_stride(opp_k3_map0_stride_s, opp_k3_map0_stride, args[2].size);

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
            const OPP_INT* opp_k3_dat2_stride_sycl = opp_k3_dat2_stride_s;
            const OPP_INT* opp_k3_map0_stride_sycl = opp_k3_map0_stride_s;
    
            const OPP_REAL* CONST_acc_coef_sycl = CONST_acc_coef_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // c_j
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // c_acc
            OPP_INT* dat2_sycl = (OPP_INT*)args[8].data_d;     // c_mask_right
            const OPP_INT* map0_sycl = args[2].map_data_d;     // c2c_map

            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            enum CellAcc {
                jfx = 0 * 4,
                jfy = 1 * 4,
                jfz = 2 * 4,
            };

            enum Dim {
                x = 0,
                y = 1,
                z = 2,
            };

            auto  accumulate_current_to_cells_kernel_sycl = [=](
                    double* cell0_j,
                    const double* cell0_acc,
                    const double* cell_xd_acc,
                    const double* cell_yd_acc,
                    const double* cell_zd_acc,
                    const double* cell_xyd_acc,
                    const double* cell_yzd_acc,
                    const double* cell_xzd_acc,
                    const int* iter_acc)
            {
                if (iter_acc[(0) * opp_k3_dat2_stride_sycl[0]] == 1)
                {
                    cell0_j[(Dim::x) * opp_k3_dat0_stride_sycl[0]] = CONST_acc_coef_sycl[Dim::x] * (cell0_acc[(CellAcc::jfx + 0) * opp_k3_dat1_stride_sycl[0]] +
                                                                cell_yd_acc[(CellAcc::jfx + 1) * opp_k3_dat1_stride_sycl[0]] +
                                                                cell_zd_acc[(CellAcc::jfx + 2) * opp_k3_dat1_stride_sycl[0]] +
                                                                cell_yzd_acc[(CellAcc::jfx + 3) * opp_k3_dat1_stride_sycl[0]]);

                    cell0_j[(Dim::y) * opp_k3_dat0_stride_sycl[0]] = CONST_acc_coef_sycl[Dim::y] * (cell0_acc[(CellAcc::jfy + 0) * opp_k3_dat1_stride_sycl[0]] +
                                                                cell_zd_acc[(CellAcc::jfy + 1) * opp_k3_dat1_stride_sycl[0]] +
                                                                cell_xd_acc[(CellAcc::jfy + 2) * opp_k3_dat1_stride_sycl[0]] +
                                                                cell_xzd_acc[(CellAcc::jfy + 3) * opp_k3_dat1_stride_sycl[0]]);

                    cell0_j[(Dim::z) * opp_k3_dat0_stride_sycl[0]] = CONST_acc_coef_sycl[Dim::z] * (cell0_acc[(CellAcc::jfz + 0) * opp_k3_dat1_stride_sycl[0]] +
                                                                cell_xd_acc[(CellAcc::jfz + 1) * opp_k3_dat1_stride_sycl[0]] +
                                                                cell_yd_acc[(CellAcc::jfz + 2) * opp_k3_dat1_stride_sycl[0]] +
                                                                cell_xyd_acc[(CellAcc::jfz + 3) * opp_k3_dat1_stride_sycl[0]]);
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                const int n = (tid + start);
                if (n < end) { 

                    accumulate_current_to_cells_kernel_sycl(
                        dat0_sycl + n, // c_j 
                        dat1_sycl + n, // c_acc 
                        dat1_sycl + map0_sycl[n + opp_k3_map0_stride_sycl[0] * 2], // c_acc 
                        dat1_sycl + map0_sycl[n + opp_k3_map0_stride_sycl[0] * 4], // c_acc 
                        dat1_sycl + map0_sycl[n + opp_k3_map0_stride_sycl[0] * 5], // c_acc 
                        dat1_sycl + map0_sycl[n + opp_k3_map0_stride_sycl[0] * 0], // c_acc 
                        dat1_sycl + map0_sycl[n + opp_k3_map0_stride_sycl[0] * 3], // c_acc 
                        dat1_sycl + map0_sycl[n + opp_k3_map0_stride_sycl[0] * 1], // c_acc 
                        dat2_sycl + n // c_mask_right 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class accumulate_current_to_cells_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("accumulate_current_to_cells_kernel");
}

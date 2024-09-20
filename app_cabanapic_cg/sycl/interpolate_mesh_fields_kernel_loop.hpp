
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k1_dat0_stride = -1;
OPP_INT opp_k1_dat1_stride = -1;
OPP_INT opp_k1_dat2_stride = -1;
OPP_INT opp_k1_dat3_stride = -1;
OPP_INT opp_k1_map0_stride = -1;

OPP_INT* opp_k1_dat0_stride_s = nullptr;
OPP_INT* opp_k1_dat1_stride_s = nullptr;
OPP_INT* opp_k1_dat2_stride_s = nullptr;
OPP_INT* opp_k1_dat3_stride_s = nullptr;
OPP_INT* opp_k1_map0_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__interpolate_mesh_fields_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // c_e | OPP_READ
    opp_arg arg1, // c_b | OPP_READ
    opp_arg arg2, // c_e | OPP_READ
    opp_arg arg3, // c_e | OPP_READ
    opp_arg arg4, // c_e | OPP_READ
    opp_arg arg5, // c_e | OPP_READ
    opp_arg arg6, // c_e | OPP_READ
    opp_arg arg7, // c_e | OPP_READ
    opp_arg arg8, // c_b | OPP_READ
    opp_arg arg9, // c_b | OPP_READ
    opp_arg arg10, // c_b | OPP_READ
    opp_arg arg11, // c_interp | OPP_WRITE
    opp_arg arg12 // c_mask_ghost | OPP_READ
) 
{
    const int nargs = 13;
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
    args[9] = arg9;
    args[10] = arg10;
    args[11] = arg11;
    args[12] = arg12;

    opp_profiler->start("interpolate_mesh_fields_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__interpolate_mesh_fields_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_set_stride(opp_k1_dat0_stride_s, opp_k1_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k1_dat1_stride_s, opp_k1_dat1_stride, args[1].dat->set->set_capacity);
    opp_set_stride(opp_k1_dat2_stride_s, opp_k1_dat2_stride, args[11].dat->set->set_capacity);
    opp_set_stride(opp_k1_dat3_stride_s, opp_k1_dat3_stride, args[12].dat->set->set_capacity);
    opp_set_stride(opp_k1_map0_stride_s, opp_k1_map0_stride, args[2].size);

#ifdef OPP_BLOCK_SIZE_1
    const int block_size = OPP_BLOCK_SIZE_1;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    if (iter_size > 0) {

        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k1_dat0_stride_sycl = opp_k1_dat0_stride_s;
            const OPP_INT* opp_k1_dat1_stride_sycl = opp_k1_dat1_stride_s;
            const OPP_INT* opp_k1_dat2_stride_sycl = opp_k1_dat2_stride_s;
            const OPP_INT* opp_k1_dat3_stride_sycl = opp_k1_dat3_stride_s;
            const OPP_INT* opp_k1_map0_stride_sycl = opp_k1_map0_stride_s;
    

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // c_e
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // c_b
            OPP_REAL* dat2_sycl = (OPP_REAL*)args[11].data_d;     // c_interp
            OPP_INT* dat3_sycl = (OPP_INT*)args[12].data_d;     // c_mask_ghost
            const OPP_INT* map0_sycl = args[2].map_data_d;     // c2c_map

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

            enum CellInterp {
                ex = 0,
                dexdy,
                dexdz,
                d2exdydz,
                ey,
                deydz,
                deydx,
                d2eydzdx,
                ez,
                dezdx,
                dezdy,
                d2ezdxdy,
                cbx,
                dcbxdx,
                cby,
                dcbydy,
                cbz,
                dcbzdz,
            };

            auto  interpolate_mesh_fields_kernel_sycl = [=](
                const double* cell0_e,
                const double* cell0_b,
                const double* cell_x_e,
                const double* cell_y_e,
                const double* cell_z_e,
                const double* cell_xy_e,
                const double* cell_yz_e,
                const double* cell_xz_e,
                const double* cell_x_b,
                const double* cell_y_b,
                const double* cell_z_b,
                double* cell0_interp,
                const int* cell0_ghost)
            {
                if (cell0_ghost[(0) * opp_k1_dat3_stride_sycl[0]] == 0)
                {
                    double w0 = 0.0, w1 = 0.0, w2 = 0.0, w3 = 0.0;

                    // ex interpolation coefficients
                    w0 = cell0_e[(Dim::x) * opp_k1_dat0_stride_sycl[0]];                       // pf0->ex;
                    w1 = cell_y_e[(Dim::x) * opp_k1_dat0_stride_sycl[0]];                      // pfy->ex;
                    w2 = cell_z_e[(Dim::x) * opp_k1_dat0_stride_sycl[0]];                      // pfz->ex;
                    w3 = cell_yz_e[(Dim::x) * opp_k1_dat0_stride_sycl[0]];                     // pfyz->ex;

                    cell0_interp[(CellInterp::ex) * opp_k1_dat2_stride_sycl[0]]       = (1.0 / 4.0)*( (w3 + w0) + (w1 + w2) );
                    cell0_interp[(CellInterp::dexdy) * opp_k1_dat2_stride_sycl[0]]    = (1.0 / 4.0)*( (w3 - w0) + (w1 - w2) );
                    cell0_interp[(CellInterp::dexdz) * opp_k1_dat2_stride_sycl[0]]    = (1.0 / 4.0)*( (w3 - w0) - (w1 - w2) );
                    cell0_interp[(CellInterp::d2exdydz) * opp_k1_dat2_stride_sycl[0]] = (1.0 / 4.0)*( (w3 + w0) - (w1 + w2) );

                    // ey interpolation coefficients
                    w0 = cell0_e[(Dim::y) * opp_k1_dat0_stride_sycl[0]];
                    w1 = cell_z_e[(Dim::y) * opp_k1_dat0_stride_sycl[0]];                       // pfz->ey;
                    w2 = cell_x_e[(Dim::y) * opp_k1_dat0_stride_sycl[0]];                       // pfx->ey;
                    w3 = cell_xz_e[(Dim::y) * opp_k1_dat0_stride_sycl[0]];                      // pfzx->ey;

                    cell0_interp[(CellInterp::ey) * opp_k1_dat2_stride_sycl[0]]       = (1.0 / 4.0)*( (w3 + w0) + (w1 + w2) );
                    cell0_interp[(CellInterp::deydz) * opp_k1_dat2_stride_sycl[0]]    = (1.0 / 4.0)*( (w3 - w0) + (w1 - w2) );
                    cell0_interp[(CellInterp::deydx) * opp_k1_dat2_stride_sycl[0]]    = (1.0 / 4.0)*( (w3 - w0) - (w1 - w2) );
                    cell0_interp[(CellInterp::d2eydzdx) * opp_k1_dat2_stride_sycl[0]] = (1.0 / 4.0)*( (w3 + w0) - (w1 + w2) );

                    // ez interpolation coefficients
                    w0 = cell0_e[(Dim::z) * opp_k1_dat0_stride_sycl[0]];                       // pf0->ez;
                    w1 = cell_x_e[(Dim::z) * opp_k1_dat0_stride_sycl[0]];                      // pfx->ez;
                    w2 = cell_y_e[(Dim::z) * opp_k1_dat0_stride_sycl[0]];                      // pfy->ez;
                    w3 = cell_xy_e[(Dim::z) * opp_k1_dat0_stride_sycl[0]];                     // pfxy->ez;

                    cell0_interp[(CellInterp::ez) * opp_k1_dat2_stride_sycl[0]]       = (1.0 / 4.0)*( (w3 + w0) + (w1 + w2) );
                    cell0_interp[(CellInterp::dezdx) * opp_k1_dat2_stride_sycl[0]]    = (1.0 / 4.0)*( (w3 - w0) + (w1 - w2) );
                    cell0_interp[(CellInterp::dezdy) * opp_k1_dat2_stride_sycl[0]]    = (1.0 / 4.0)*( (w3 - w0) - (w1 - w2) );
                    cell0_interp[(CellInterp::d2ezdxdy) * opp_k1_dat2_stride_sycl[0]] = (1.0 / 4.0)*( (w3 + w0) - (w1 + w2) );

                    // bx interpolation coefficients
                    w0 = cell0_b[(Dim::x) * opp_k1_dat1_stride_sycl[0]];                      // pf0->cbx;
                    w1 = cell_x_b[(Dim::x) * opp_k1_dat1_stride_sycl[0]];                     // pfx->cbx;
                    cell0_interp[(CellInterp::cbx) * opp_k1_dat2_stride_sycl[0]]    = (1.0 / 2.0)*( w1 + w0 );
                    cell0_interp[(CellInterp::dcbxdx) * opp_k1_dat2_stride_sycl[0]] = (1.0 / 2.0)*( w1 - w0 );

                    // by interpolation coefficients
                    w0 = cell0_b[(Dim::y) * opp_k1_dat1_stride_sycl[0]];                      // pf0->cby;
                    w1 = cell_y_b[(Dim::y) * opp_k1_dat1_stride_sycl[0]];                     // pfy->cby;
                    cell0_interp[(CellInterp::cby) * opp_k1_dat2_stride_sycl[0]]    = (1.0 / 2.0)*( w1 + w0 );
                    cell0_interp[(CellInterp::dcbydy) * opp_k1_dat2_stride_sycl[0]] = (1.0 / 2.0)*( w1 - w0 );

                    // bz interpolation coefficients
                    w0 = cell0_b[(Dim::z) * opp_k1_dat1_stride_sycl[0]];                      // pf0->cbz;
                    w1 = cell_z_b[(Dim::z) * opp_k1_dat1_stride_sycl[0]];                     // pfz->cbz;
                    cell0_interp[(CellInterp::cbz) * opp_k1_dat2_stride_sycl[0]]    = (1.0 / 2.0)*( w1 + w0 );
                    cell0_interp[(CellInterp::dcbzdz) * opp_k1_dat2_stride_sycl[0]] = (1.0 / 2.0)*( w1 - w0 );
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                const int n = (tid + start);
                if (n < end) { 

                    interpolate_mesh_fields_kernel_sycl(
                        dat0_sycl + n, // c_e 
                        dat1_sycl + n, // c_b 
                        dat0_sycl + map0_sycl[n + opp_k1_map0_stride_sycl[0] * 9], // c_e 
                        dat0_sycl + map0_sycl[n + opp_k1_map0_stride_sycl[0] * 7], // c_e 
                        dat0_sycl + map0_sycl[n + opp_k1_map0_stride_sycl[0] * 6], // c_e 
                        dat0_sycl + map0_sycl[n + opp_k1_map0_stride_sycl[0] * 11], // c_e 
                        dat0_sycl + map0_sycl[n + opp_k1_map0_stride_sycl[0] * 8], // c_e 
                        dat0_sycl + map0_sycl[n + opp_k1_map0_stride_sycl[0] * 10], // c_e 
                        dat1_sycl + map0_sycl[n + opp_k1_map0_stride_sycl[0] * 9], // c_b 
                        dat1_sycl + map0_sycl[n + opp_k1_map0_stride_sycl[0] * 7], // c_b 
                        dat1_sycl + map0_sycl[n + opp_k1_map0_stride_sycl[0] * 6], // c_b 
                        dat2_sycl + n, // c_interp 
                        dat3_sycl + n // c_mask_ghost 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class interpolate_mesh_fields_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    opp_queue->wait();   
 
    opp_profiler->end("interpolate_mesh_fields_kernel");
}


//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k1_dat0_stride = -1;
OPP_INT opp_k1_dat1_stride = -1;
OPP_INT opp_k1_dat2_stride = -1;
OPP_INT opp_k1_map0_stride = -1;

OPP_INT* opp_k1_dat0_stride_s = nullptr;
OPP_INT* opp_k1_dat1_stride_s = nullptr;
OPP_INT* opp_k1_dat2_stride_s = nullptr;
OPP_INT* opp_k1_map0_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__weight_f2p_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // n_field_e | OPP_READ
    opp_arg arg1, // n_field_e | OPP_READ
    opp_arg arg2, // p_pos_x | OPP_READ
    opp_arg arg3 // p_field_e | OPP_WRITE
) 
{
    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    opp_profiler->start("weight_f2p_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__weight_f2p_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_set_stride(opp_k1_dat0_stride_s, opp_k1_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k1_dat1_stride_s, opp_k1_dat1_stride, args[2].dat->set->set_capacity);
    opp_set_stride(opp_k1_dat2_stride_s, opp_k1_dat2_stride, args[3].dat->set->set_capacity);
    opp_set_stride(opp_k1_map0_stride_s, opp_k1_map0_stride, args[0].size);

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
            const OPP_INT* opp_k1_map0_stride_sycl = opp_k1_map0_stride_s;
    
            const OPP_REAL* CONST_xl_sycl = CONST_xl_s;
            const OPP_REAL* CONST_dx_sycl = CONST_dx_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // n_field_e
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[2].data_d;     // p_pos_x
            OPP_REAL* dat2_sycl = (OPP_REAL*)args[3].data_d;     // p_field_e
            const OPP_INT* map0_sycl = args[0].map_data_d;     // c2n_map

            const OPP_INT *p2c_map = (OPP_INT *)set->mesh_relation_dat->data_d;
            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
             weight_f2p_kernel_sycl = [=](
                    const double* node0_field_E,  //LHS
                    const double* node1_field_E,  //RHS
                    const double* particle0_position_x,
                    double* particle0_field_E
                )
            {
                double xx = ((particle0_position_x[(0) * opp_k1_dat1_stride_sycl[0]] - CONST_xl_sycl[0]) / CONST_dx_sycl[0]); // Makes Global position to local position comapared to the cell
                int n = int(xx);
                double frac = (xx - n);

                particle0_field_E[(0) * opp_k1_dat2_stride_sycl[0]] = ((frac * node1_field_E[(0) * opp_k1_dat0_stride_sycl[0]]) + ((1.0 - frac) * node0_field_E[(0) * opp_k1_dat0_stride_sycl[0]]));
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                const int n = (tid + start);
                if (n < end) { 
                    const OPP_INT* opp_p2c = p2c_map + n;

                    weight_f2p_kernel_sycl(
                        dat0_sycl + map0_sycl[opp_p2c[0] + opp_k1_map0_stride_sycl[0] * 0], // n_field_e 
                        dat0_sycl + map0_sycl[opp_p2c[0] + opp_k1_map0_stride_sycl[0] * 1], // n_field_e 
                        dat1_sycl + n, // p_pos_x 
                        dat2_sycl + n // p_field_e 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class weight_f2p_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("weight_f2p_kernel");
}

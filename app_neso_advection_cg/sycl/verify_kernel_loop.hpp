
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k3_dat0_stride = -1;
OPP_INT opp_k3_dat1_stride = -1;

OPP_INT* opp_k3_dat0_stride_s = nullptr;
OPP_INT* opp_k3_dat1_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__verify_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1, // c_idx | OPP_READ
    opp_arg arg2 // | OPP_INC
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("verify_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__verify_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    OPP_INT *arg2_host_data = (OPP_INT *)args[2].data;

    opp_set_stride(opp_k3_dat0_stride_s, opp_k3_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k3_dat1_stride_s, opp_k3_dat1_stride, args[1].dat->set->set_capacity);

#ifdef OPP_BLOCK_SIZE_3
    const int block_size = OPP_BLOCK_SIZE_3;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    const int num_blocks = 200;

    int reduction_bytes = 0;
    int reduction_size = 0;

    reduction_bytes += ROUND_UP(num_blocks * 1 * sizeof(OPP_INT));
    reduction_size   = MAX(reduction_size, sizeof(OPP_INT));

    opp_reallocReductArrays(reduction_bytes);
    reduction_bytes = 0;

    args[2].data   = OPP_reduct_h + reduction_bytes;
    args[2].data_d = OPP_reduct_d + reduction_bytes;

    for (int b = 0; b < num_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            ((OPP_INT *)args[2].data)[b * 1 + d] = OPP_INT_ZERO;
    }

    reduction_bytes += ROUND_UP(num_blocks * 1 * sizeof(OPP_INT));

    opp_mvReductArraysToDevice(reduction_bytes);
    
    if (iter_size > 0) {

        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k3_dat0_stride_sycl = opp_k3_dat0_stride_s;
            const OPP_INT* opp_k3_dat1_stride_sycl = opp_k3_dat1_stride_s;
    
            const OPP_INT* CONST_ndimcells_sycl = CONST_ndimcells_s;
            const OPP_REAL* CONST_cell_width_sycl = CONST_cell_width_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // p_pos
            OPP_INT* dat1_sycl = (OPP_INT*)args[1].data_d;     // c_idx
            OPP_INT* gbl2_sycl = (OPP_INT*)args[2].data_d;
            sycl::accessor<OPP_INT, 1, sycl::access::mode::read_write, sycl::access::target::local>
                                        red_OPP_INT_2(block_size, cgh); // temp var for reduction

            const OPP_INT *p2c_map = (OPP_INT *)set->mesh_relation_dat->data_d;
            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            enum Dim {
                x = 0,
                y = 1,
            };

            auto  verify_kernel_sycl = [=](
                    const double* part_pos,
                    const int* cell_global_idx,
                    int* incorrect_part_count)
            {
                // get the cell boundaries for the current cell_index - using global cell index
                int ix = -1, iy = -1;
                int _ix, _iy; _ix  = ((*cell_global_idx)); _iy  = _ix/int(CONST_ndimcells_sycl[Dim::x]); _ix -= _iy*int(CONST_ndimcells_sycl[Dim::x]); (ix) = _ix; (iy) = _iy;;

                if (ix < 0 || iy < 0)
                {
                    // opp_printf("VERIFY", "Incorrect ix[%d] iy[%d] for global cell[%d] nx[%d]",
                    //     ix, iy, (*cell_global_idx), CONST_ndimcells[Dim::x]);
                    (*incorrect_part_count)++;
                    return;
                }

                // get the boundaries of that cell
                const double boundary_ll[2] = { (ix * CONST_cell_width_sycl[0]), (iy * CONST_cell_width_sycl[0]) };

                // check whether the current particle is within those boundaries or not!
                const double part_pos_x = part_pos[(Dim::x) * opp_k3_dat0_stride_sycl[0]];
                if (part_pos_x < boundary_ll[Dim::x] ||
                        part_pos_x > (boundary_ll[Dim::x] + CONST_cell_width_sycl[0])) {

                    (*incorrect_part_count)++;
                    return;
                }

                const double part_pos_y = part_pos[(Dim::y) * opp_k3_dat0_stride_sycl[0]];
                if (part_pos_y < boundary_ll[Dim::y] ||
                        part_pos_y > (boundary_ll[Dim::y] + CONST_cell_width_sycl[0])) {

                    (*incorrect_part_count)++;
                    return;
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                OPP_INT gbl2_local[1];
                for (int d = 0; d < 1; ++d)
                    gbl2_local[d] = OPP_INT_ZERO;

                const int tid = item.get_global_linear_id();
                for (int n = tid; n < iter_size; n += item.get_global_range()[0]) {

                    const OPP_INT* opp_p2c = p2c_map + n;

                    verify_kernel_sycl(
                        dat0_sycl + n, // p_pos 
                        dat1_sycl + opp_p2c[0], // c_idx 
                        gbl2_local // 
                    );
                }

                for (int d = 0; d < 1; ++d) //arg2_offset + 
                    opp_reduction<OPP_INC, 0>(gbl2_sycl, (d + item.get_group_linear_id() * 1), 
                                    gbl2_local[d], red_OPP_INT_2, item);
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class verify_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_mvReductArraysToHost(reduction_bytes);

    for (int b = 0; b < num_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            arg2_host_data[d] += ((OPP_INT *)args[2].data)[b * 1 + d];
    }

    args[2].data = (char *)arg2_host_data;
    opp_mpi_reduce(&args[2], arg2_host_data);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    opp_queue->wait();   
 
    opp_profiler->end("verify_kernel");
}

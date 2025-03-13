
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k1_dat0_stride = -1;
OPP_INT opp_k1_dat1_stride = -1;
OPP_INT opp_k1_dat2_stride = -1;

OPP_INT* opp_k1_dat0_stride_s = nullptr;
OPP_INT* opp_k1_dat1_stride_s = nullptr;
OPP_INT* opp_k1_dat2_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__update_pos_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // p_vel | OPP_READ
    opp_arg arg1, // p_pos | OPP_RW
    opp_arg arg2 // p_mdir | OPP_WRITE
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("update_pos_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__update_pos_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_set_stride(opp_k1_dat0_stride_s, opp_k1_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k1_dat1_stride_s, opp_k1_dat1_stride, args[1].dat->set->set_capacity);
    opp_set_stride(opp_k1_dat2_stride_s, opp_k1_dat2_stride, args[2].dat->set->set_capacity);

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
    
            const OPP_REAL* CONST_dt_sycl = CONST_dt_s;
            const OPP_REAL* CONST_extents_sycl = CONST_extents_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // p_vel
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // p_pos
            OPP_INT* dat2_sycl = (OPP_INT*)args[2].data_d;     // p_mdir

            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            auto  update_pos_kernel_sycl = [=](const double* p_vel, double* p_pos, int* p_mdir)
            {
                for (int dm = 0; dm < 2; dm++) {

                    const double offset = p_vel[(dm) * opp_k1_dat0_stride_sycl[0]] * CONST_dt_sycl[0];
                    p_pos[(dm) * opp_k1_dat1_stride_sycl[0]] += offset; // s1 = s0 + ut

                    // correct for periodic boundary conditions
                    const int n_extent_offset_int = std::abs(p_pos[(dm) * opp_k1_dat1_stride_sycl[0]]) + 2.0;
                    const double temp_pos = p_pos[(dm) * opp_k1_dat1_stride_sycl[0]] + n_extent_offset_int * CONST_extents_sycl[dm];
                    p_pos[(dm) * opp_k1_dat1_stride_sycl[0]] = std::fmod(temp_pos, CONST_extents_sycl[dm]);

                    p_mdir[(dm) * opp_k1_dat2_stride_sycl[0]] = (offset > 0) ? 1 : -1;
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                for (int n = tid; n < iter_size; n += item.get_global_range()[0]) {

                    update_pos_kernel_sycl(
                        dat0_sycl + n, // p_vel 
                        dat1_sycl + n, // p_pos 
                        dat2_sycl + n // p_mdir 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class update_pos_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("update_pos_kernel");
}

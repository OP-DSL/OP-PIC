
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k4_dat0_stride = -1;
OPP_INT opp_k4_dat1_stride = -1;

OPP_INT* opp_k4_dat0_stride_s = nullptr;
OPP_INT* opp_k4_dat1_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__sum_laplace_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // n_xlocal | OPP_READ
    opp_arg arg1 // n_field_p | OPP_RW
) 
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("sum_laplace_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__sum_laplace_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_set_stride(opp_k4_dat0_stride_s, opp_k4_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k4_dat1_stride_s, opp_k4_dat1_stride, args[1].dat->set->set_capacity);

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
    
            const OPP_REAL* CONST_L_sycl = CONST_L_s;
            const OPP_REAL* CONST_lhs_voltage_sycl = CONST_lhs_voltage_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // n_xlocal
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // n_field_p

            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
             sum_laplace_kernel_sycl = [=](
                    const double* node0_xlocal,
                    double* node0_field_P
                )
            {
                double rv = 0.0;
                double lv = CONST_lhs_voltage_sycl[0];

                double frac = ((*node0_xlocal) / CONST_L_sycl[0]);
                (*node0_field_P) += (frac * rv + (1. - frac) * lv);
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                for (int n = tid; n < iter_size; n += item.get_global_range()[0]) {

                    sum_laplace_kernel_sycl(
                        dat0_sycl + n, // n_xlocal 
                        dat1_sycl + n // n_field_p 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class sum_laplace_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("sum_laplace_kernel");
}

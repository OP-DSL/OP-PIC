
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k6_dat0_stride = -1;
OPP_INT opp_k6_dat1_stride = -1;

OPP_INT* opp_k6_dat0_stride_s = nullptr;
OPP_INT* opp_k6_dat1_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__compute_node_charge_density_kernel(opp_set set,
    opp_arg arg0, // n_charge_den | OPP_RW
    opp_arg arg1 // n_volume | OPP_READ
) 
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("compute_node_charge_density_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__compute_node_charge_density_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_set_stride(opp_k6_dat0_stride_s, opp_k6_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k6_dat1_stride_s, opp_k6_dat1_stride, args[1].dat->set->set_capacity);

#ifdef OPP_BLOCK_SIZE_6
    const int block_size = OPP_BLOCK_SIZE_6;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    if (iter_size > 0) {

        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k6_dat0_stride_sycl = opp_k6_dat0_stride_s;
            const OPP_INT* opp_k6_dat1_stride_sycl = opp_k6_dat1_stride_s;
    
            const OPP_REAL* CONST_spwt_sycl = CONST_spwt_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // n_charge_den
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // n_volume

            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            auto  compute_node_charge_density_kernel_sycl = [=](
                double *node_charge_den,
                const double *node_volume
            ) {
                (*node_charge_den) *= (CONST_spwt_sycl[0] / node_volume[(0) * opp_k6_dat1_stride_sycl[0]]);
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                for (int n = tid; n < iter_size; n += item.get_global_range()[0]) {

                    compute_node_charge_density_kernel_sycl(
                        dat0_sycl + n, // n_charge_den 
                        dat1_sycl + n // n_volume 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class compute_node_charge_density_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("compute_node_charge_density_kernel");
}

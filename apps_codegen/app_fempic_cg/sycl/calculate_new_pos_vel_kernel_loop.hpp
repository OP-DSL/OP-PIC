
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k3_dat0_stride = -1;
OPP_INT opp_k3_dat1_stride = -1;
OPP_INT opp_k3_dat2_stride = -1;

OPP_INT* opp_k3_dat0_stride_s = nullptr;
OPP_INT* opp_k3_dat1_stride_s = nullptr;
OPP_INT* opp_k3_dat2_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_all__calculate_new_pos_vel_kernel(opp_set set,
    opp_arg arg0, // c_ef | OPP_READ
    opp_arg arg1, // p_pos | OPP_WRITE
    opp_arg arg2 // p_vel | OPP_WRITE
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("calculate_new_pos_vel_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__calculate_new_pos_vel_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_set_stride(opp_k3_dat0_stride_s, opp_k3_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k3_dat1_stride_s, opp_k3_dat1_stride, args[1].dat->set->set_capacity);
    opp_set_stride(opp_k3_dat2_stride_s, opp_k3_dat2_stride, args[2].dat->set->set_capacity);

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
    
            const OPP_REAL* CONST_charge_sycl = CONST_charge_s;
            const OPP_REAL* CONST_dt_sycl = CONST_dt_s;
            const OPP_REAL* CONST_mass_sycl = CONST_mass_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // c_ef
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // p_pos
            OPP_REAL* dat2_sycl = (OPP_REAL*)args[2].data_d;     // p_vel

            const OPP_INT *p2c_map = (OPP_INT *)set->mesh_relation_dat->data_d;
            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            auto  calculate_new_pos_vel_kernel_sycl = [=](
                const double *cell_ef,
                double *part_pos,
                double *part_vel
            ) {
                const double coefficient1 = CONST_charge_sycl[0] / CONST_mass_sycl[0] * (CONST_dt_sycl[0]);
                for (int i = 0; i < 3; i++) {
                    part_vel[(i) * opp_k3_dat2_stride_sycl[0]] += (coefficient1 * cell_ef[(i) * opp_k3_dat0_stride_sycl[0]]);
                    part_pos[(i) * opp_k3_dat1_stride_sycl[0]] += part_vel[(i) * opp_k3_dat2_stride_sycl[0]] * (CONST_dt_sycl[0]);
                }
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                const int n = (tid + start);
                if (n < end) { 
                    const OPP_INT* opp_p2c = p2c_map + n;

                    calculate_new_pos_vel_kernel_sycl(
                        dat0_sycl + opp_p2c[0], // c_ef 
                        dat1_sycl + n, // p_pos 
                        dat2_sycl + n // p_vel 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class calculate_new_pos_vel_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("calculate_new_pos_vel_kernel");
}

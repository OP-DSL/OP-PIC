
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k2_dat0_stride = -1;
OPP_INT opp_k2_dat1_stride = -1;
OPP_INT opp_k2_dat2_stride = -1;
OPP_INT opp_k2_dat3_stride = -1;
OPP_INT opp_k2_dat4_stride = -1;
OPP_INT opp_k2_dat5_stride = -1;
OPP_INT opp_k2_dat6_stride = -1;
OPP_INT opp_k2_dat7_stride = -1;
OPP_INT opp_k2_dat8_stride = -1;
OPP_INT opp_k2_dat9_stride = -1;
OPP_INT opp_k2_map0_stride = -1;

OPP_INT* opp_k2_dat0_stride_s = nullptr;
OPP_INT* opp_k2_dat1_stride_s = nullptr;
OPP_INT* opp_k2_dat2_stride_s = nullptr;
OPP_INT* opp_k2_dat3_stride_s = nullptr;
OPP_INT* opp_k2_dat4_stride_s = nullptr;
OPP_INT* opp_k2_dat5_stride_s = nullptr;
OPP_INT* opp_k2_dat6_stride_s = nullptr;
OPP_INT* opp_k2_dat7_stride_s = nullptr;
OPP_INT* opp_k2_dat8_stride_s = nullptr;
OPP_INT* opp_k2_dat9_stride_s = nullptr;
OPP_INT* opp_k2_map0_stride_s = nullptr;

//--------------------------------------------------------------
void opp_par_loop_injected__inject_ions_kernel(opp_set set,
    opp_arg arg0, // p_pos | OPP_WRITE
    opp_arg arg1, // p_vel | OPP_WRITE
    opp_arg arg2, // p2c_map | OPP_RW
    opp_arg arg3, // if2c_map | OPP_READ
    opp_arg arg4, // c_ef | OPP_READ
    opp_arg arg5, // if_u_norm | OPP_READ
    opp_arg arg6, // if_v_norm | OPP_READ
    opp_arg arg7, // if_norm | OPP_READ
    opp_arg arg8, // if_n_pos | OPP_READ
    opp_arg arg9 // dp_rand | OPP_READ
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const int nargs = 10;
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

    opp_profiler->start("inject_ions_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_injected__inject_ions_kernel set_size %d", set->size);

    opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    const int iter_size = set->diff; 
    const int inj_start = (set->size - set->diff);  
 
    opp_set_stride(opp_k2_dat0_stride_s, opp_k2_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat1_stride_s, opp_k2_dat1_stride, args[1].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat2_stride_s, opp_k2_dat2_stride, args[2].size);
    opp_set_stride(opp_k2_dat3_stride_s, opp_k2_dat3_stride, args[3].size);
    opp_set_stride(opp_k2_dat4_stride_s, opp_k2_dat4_stride, args[4].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat5_stride_s, opp_k2_dat5_stride, args[5].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat6_stride_s, opp_k2_dat6_stride, args[6].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat7_stride_s, opp_k2_dat7_stride, args[7].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat8_stride_s, opp_k2_dat8_stride, args[8].dat->set->set_capacity);
    opp_set_stride(opp_k2_dat9_stride_s, opp_k2_dat9_stride, args[9].dat->set->set_capacity);
    opp_set_stride(opp_k2_map0_stride_s, opp_k2_map0_stride, args[4].size);

#ifdef OPP_BLOCK_SIZE_2
    const int block_size = OPP_BLOCK_SIZE_2;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    if (iter_size > 0) {
        
        opp_queue->submit([&](sycl::handler &cgh) {

            const OPP_INT* opp_k2_dat0_stride_sycl = opp_k2_dat0_stride_s;
            const OPP_INT* opp_k2_dat1_stride_sycl = opp_k2_dat1_stride_s;
            const OPP_INT* opp_k2_dat2_stride_sycl = opp_k2_dat2_stride_s;
            const OPP_INT* opp_k2_dat3_stride_sycl = opp_k2_dat3_stride_s;
            const OPP_INT* opp_k2_dat4_stride_sycl = opp_k2_dat4_stride_s;
            const OPP_INT* opp_k2_dat5_stride_sycl = opp_k2_dat5_stride_s;
            const OPP_INT* opp_k2_dat6_stride_sycl = opp_k2_dat6_stride_s;
            const OPP_INT* opp_k2_dat7_stride_sycl = opp_k2_dat7_stride_s;
            const OPP_INT* opp_k2_dat8_stride_sycl = opp_k2_dat8_stride_s;
            const OPP_INT* opp_k2_dat9_stride_sycl = opp_k2_dat9_stride_s;
            const OPP_INT* opp_k2_map0_stride_sycl = opp_k2_map0_stride_s;
    
            const OPP_REAL* CONST_charge_sycl = CONST_charge_s;
            const OPP_REAL* CONST_dt_sycl = CONST_dt_s;
            const OPP_REAL* CONST_ion_velocity_sycl = CONST_ion_velocity_s;
            const OPP_REAL* CONST_mass_sycl = CONST_mass_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // p_pos
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // p_vel
            OPP_INT* dat2_sycl = (OPP_INT*)args[2].data_d;     // p2c_map
            OPP_INT* dat3_sycl = (OPP_INT*)args[3].data_d;     // if2c_map
            OPP_REAL* dat4_sycl = (OPP_REAL*)args[4].data_d;     // c_ef
            OPP_REAL* dat5_sycl = (OPP_REAL*)args[5].data_d;     // if_u_norm
            OPP_REAL* dat6_sycl = (OPP_REAL*)args[6].data_d;     // if_v_norm
            OPP_REAL* dat7_sycl = (OPP_REAL*)args[7].data_d;     // if_norm
            OPP_REAL* dat8_sycl = (OPP_REAL*)args[8].data_d;     // if_n_pos
            OPP_REAL* dat9_sycl = (OPP_REAL*)args[9].data_d;     // dp_rand
            const OPP_INT* map0_sycl = args[4].map_data_d;     // if2c_map

            const OPP_INT *p2c_map = (OPP_INT *)set->mesh_relation_dat->data_d;
            const OPP_INT start = 0;
            const OPP_INT end = iter_size;

            const int num_blocks = (end - start - 1) / block_size + 1;

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            auto  inject_ions_kernel_sycl = [=](
                double *part_pos,
                double *part_vel,
                int *part_cell_connectivity,
                const int *cell_id,
                const double *cell_ef,
                const double *iface_u,
                const double *iface_v,
                const double *iface_normal,
                const double *node_pos,
                const double* dummy_part_random
            ) {
                double a = dummy_part_random[(0) * opp_k2_dat9_stride_sycl[0]];
                double b = dummy_part_random[(1) * opp_k2_dat9_stride_sycl[0]];
                if ((a + b) > 1) {  // TODO : Change the random dat to avoid this
                    a = (1 - a);
                    b = (1 - b);
                }

                for (int i = 0; i < 3; i++) {
                    part_pos[(i) * opp_k2_dat0_stride_sycl[0]] = a * iface_u[(i) * opp_k2_dat5_stride_sycl[0]] + b * iface_v[(i) * opp_k2_dat6_stride_sycl[0]] + node_pos[(i) * opp_k2_dat8_stride_sycl[0]];

                    part_vel[(i) * opp_k2_dat1_stride_sycl[0]] = (iface_normal[(i) * opp_k2_dat7_stride_sycl[0]] * CONST_ion_velocity_sycl[0]);
                    part_vel[(i) * opp_k2_dat1_stride_sycl[0]] -= CONST_charge_sycl[0] / CONST_mass_sycl[0] * cell_ef[(i) * opp_k2_dat4_stride_sycl[0]] * (0.5 * CONST_dt_sycl[0]);
                }

                (*part_cell_connectivity) = (*cell_id);
            };

            // -----------------------------------------------------------------------------------------
            auto kernel = [=](sycl::nd_item<1> item) {

                const int tid = item.get_global_linear_id();
                const int n = (tid + start);
                if (n < end) { 
                    const OPP_INT* opp_p2c = p2c_map + (inj_start + n);

                    inject_ions_kernel_sycl(
                        dat0_sycl + (inj_start + n), // p_pos 
                        dat1_sycl + (inj_start + n), // p_vel 
                        dat2_sycl + (inj_start + n), // p2c_map 
                        dat3_sycl + opp_p2c[0], // if2c_map 
                        dat4_sycl + map0_sycl[opp_p2c[0] + opp_k2_map0_stride_sycl[0] * 0], // c_ef 
                        dat5_sycl + opp_p2c[0], // if_u_norm 
                        dat6_sycl + opp_p2c[0], // if_v_norm 
                        dat7_sycl + opp_p2c[0], // if_norm 
                        dat8_sycl + opp_p2c[0], // if_n_pos 
                        dat9_sycl + n // dp_rand 
                    );
                }
    
            };
 
            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class inject_ions_kernel>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), kernel);
        });
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("inject_ions_kernel");
}

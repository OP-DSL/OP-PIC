
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

__constant__ OPP_INT opp_k2_dat0_stride_d;
__constant__ OPP_INT opp_k2_dat1_stride_d;
__constant__ OPP_INT opp_k2_dat2_stride_d;
__constant__ OPP_INT opp_k2_dat3_stride_d;
__constant__ OPP_INT opp_k2_dat4_stride_d;
__constant__ OPP_INT opp_k2_dat5_stride_d;
__constant__ OPP_INT opp_k2_dat6_stride_d;
__constant__ OPP_INT opp_k2_dat7_stride_d;
__constant__ OPP_INT opp_k2_dat8_stride_d;
__constant__ OPP_INT opp_k2_dat9_stride_d;
__constant__ OPP_INT opp_k2_map0_stride_d;



namespace opp_k2 {
__device__ inline void inject_ions_kernel(
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
    double a = dummy_part_random[(0) * opp_k2_dat9_stride_d];
    double b = dummy_part_random[(1) * opp_k2_dat9_stride_d];
    if ((a + b) > 1) {  // TODO : Change the random dat to avoid this
        a = (1 - a);
        b = (1 - b);
    }

    for (int i = 0; i < 3; i++) {
        part_pos[(i) * opp_k2_dat0_stride_d] = a * iface_u[(i) * opp_k2_dat5_stride_d] + b * iface_v[(i) * opp_k2_dat6_stride_d] + node_pos[(i) * opp_k2_dat8_stride_d];

        part_vel[(i) * opp_k2_dat1_stride_d] = (iface_normal[(i) * opp_k2_dat7_stride_d] * CONST_ion_velocity_d[0]);
        part_vel[(i) * opp_k2_dat1_stride_d] -= CONST_charge_d[0] / CONST_mass_d[0] * cell_ef[(i) * opp_k2_dat4_stride_d] * (0.5 * CONST_dt_d[0]);
    }

    (*part_cell_connectivity) = (*cell_id);
}

}

//--------------------------------------------------------------
__global__ void opp_dev_inject_ions_kernel(
    OPP_REAL *__restrict__ dat0,  // p_pos
    OPP_REAL *__restrict__ dat1,  // p_vel
    OPP_INT *__restrict__ dat2,  // p2c_map
    const OPP_INT *__restrict__ dat3,  // if2c_map
    const OPP_REAL *__restrict__ dat4,  // c_ef
    const OPP_REAL *__restrict__ dat5,  // if_u_norm
    const OPP_REAL *__restrict__ dat6,  // if_v_norm
    const OPP_REAL *__restrict__ dat7,  // if_norm
    const OPP_REAL *__restrict__ dat8,  // if_n_pos
    const OPP_REAL *__restrict__ dat9,  // dp_rand
    const OPP_INT *__restrict__ p2c_map,
    const OPP_INT *__restrict__ map0,  // if2c_map
    const OPP_INT inj_start,
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        const OPP_INT* opp_p2c = p2c_map + (inj_start + n);
        
        opp_k2::inject_ions_kernel(
            dat0 + (inj_start + n), // p_pos 
            dat1 + (inj_start + n), // p_vel 
            dat2 + (inj_start + n), // p2c_map 
            dat3 + opp_p2c[0], // if2c_map 
            dat4 + map0[opp_p2c[0] + opp_k2_map0_stride_d * 0], // c_ef 
            dat5 + opp_p2c[0], // if_u_norm 
            dat6 + opp_p2c[0], // if_v_norm 
            dat7 + opp_p2c[0], // if_norm 
            dat8 + opp_p2c[0], // if_n_pos 
            dat9 + n // dp_rand 
        );
    }
    
}

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
 
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat0_stride_d, &opp_k2_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat1_stride_d, &opp_k2_dat1_stride, &(args[1].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat2_stride_d, &opp_k2_dat2_stride, &(args[2].size), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat3_stride_d, &opp_k2_dat3_stride, &(args[3].size), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat4_stride_d, &opp_k2_dat4_stride, &(args[4].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat5_stride_d, &opp_k2_dat5_stride, &(args[5].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat6_stride_d, &opp_k2_dat6_stride, &(args[6].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat7_stride_d, &opp_k2_dat7_stride, &(args[7].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat8_stride_d, &opp_k2_dat8_stride, &(args[8].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat9_stride_d, &opp_k2_dat9_stride, &(args[9].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_map0_stride_d, &opp_k2_map0_stride, &(args[4].size), 1);

#ifdef OPP_BLOCK_SIZE_2
    const int block_size = OPP_BLOCK_SIZE_2;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    int num_blocks = 200;
    if (iter_size > 0) 
    {
        const OPP_INT start = 0;
        const OPP_INT end = iter_size;
        num_blocks = (end - start - 1) / block_size + 1;

        {
            opp_dev_inject_ions_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // p_pos
                (OPP_REAL *)args[1].data_d,     // p_vel
                (OPP_INT *)args[2].data_d,     // p2c_map
                (OPP_INT *)args[3].data_d,     // if2c_map
                (OPP_REAL *)args[4].data_d,     // c_ef
                (OPP_REAL *)args[5].data_d,     // if_u_norm
                (OPP_REAL *)args[6].data_d,     // if_v_norm
                (OPP_REAL *)args[7].data_d,     // if_norm
                (OPP_REAL *)args[8].data_d,     // if_n_pos
                (OPP_REAL *)args[9].data_d,     // dp_rand
                (OPP_INT *)set->mesh_relation_dat->data_d,
                args[4].map_data_d,     // if2c_map
                inj_start,
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("inject_ions_kernel");
}


//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k2 {
inline void inject_ions_kernel(
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
    double a = dummy_part_random[0];
    double b = dummy_part_random[1];
    if ((a + b) > 1) {  // TODO : Change the random dat to avoid this
        a = (1 - a);
        b = (1 - b);
    }

    for (int i = 0; i < 3; i++) {
        part_pos[i] = a * iface_u[i] + b * iface_v[i] + node_pos[i];

        part_vel[i] = (iface_normal[i] * CONST_ion_velocity[0]);
        part_vel[i] -= CONST_charge[0] / CONST_mass[0] * cell_ef[i] * (0.5 * CONST_dt[0]);
    }

    (*part_cell_connectivity) = (*cell_id);
}
}

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
{
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

    const int nthreads = omp_get_max_threads(); 

    opp_mpi_halo_exchanges(set, nargs, args);
 
    opp_mpi_halo_wait_all(nargs, args);

    const int iter_size = set->diff; 
    const int inj_start = (set->size - set->diff);  
    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 
  
    
  
 
    #pragma omp parallel for 
    for (int thr = 0; thr < nthreads; thr++)
    {
        const size_t start  = ((size_t)iter_size * thr) / nthreads;
        const size_t finish = ((size_t)iter_size * (thr+1)) / nthreads;
      
        for (size_t n = start; n < finish; n++)
        { 
            OPP_INT* opp_p2c = OPP_mesh_relation_data + inj_start + n;
            const OPP_INT *map0 = args[4].map_data + (opp_p2c[0] * 1);
   
            opp_k2::inject_ions_kernel(
                (OPP_REAL *)args[0].data + ((inj_start + n) * 3),
                (OPP_REAL *)args[1].data + ((inj_start + n) * 3),
                (OPP_INT *)args[2].data + ((inj_start + n) * 1),
                (const OPP_INT *)args[3].data + (opp_p2c[0] * 1),
                (const OPP_REAL *)args[4].data + (map0[0] * 3),
                (const OPP_REAL *)args[5].data + (opp_p2c[0] * 3),
                (const OPP_REAL *)args[6].data + (opp_p2c[0] * 3),
                (const OPP_REAL *)args[7].data + (opp_p2c[0] * 3),
                (const OPP_REAL *)args[8].data + (opp_p2c[0] * 9),
                (const OPP_REAL *)args[9].data + (n * 2)
            );
        }
    }
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("inject_ions_kernel");
}

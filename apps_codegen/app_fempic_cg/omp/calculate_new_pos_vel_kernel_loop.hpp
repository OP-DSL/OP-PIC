
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k3 {
inline void calculate_new_pos_vel_kernel(
    const double *cell_ef,
    double *part_pos,
    double *part_vel
) {
    const double coefficient1 = CONST_charge[0] / CONST_mass[0] * (CONST_dt[0]);
    for (int i = 0; i < 3; i++) {
        part_vel[i] += (coefficient1 * cell_ef[i]);
        part_pos[i] += part_vel[i] * (CONST_dt[0]);
    }
}
}

void opp_par_loop_all__calculate_new_pos_vel_kernel(opp_set set,
    opp_arg arg0, // c_ef | OPP_READ
    opp_arg arg1, // p_pos | OPP_WRITE
    opp_arg arg2 // p_vel | OPP_WRITE
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("calculate_new_pos_vel_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__calculate_new_pos_vel_kernel set_size %d", set->size);

    const int nthreads = omp_get_max_threads(); 

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
 
    opp_mpi_halo_wait_all(nargs, args);    
    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 
  
    
  
 
    #pragma omp parallel for 
    for (int thr = 0; thr < nthreads; thr++)
    {
        const size_t start  = ((size_t)iter_size * thr) / nthreads;
        const size_t finish = ((size_t)iter_size * (thr+1)) / nthreads;
      
        for (size_t n = start; n < finish; n++)
        { 
            OPP_INT* opp_p2c = OPP_mesh_relation_data + n;
   
            opp_k3::calculate_new_pos_vel_kernel(
                (const OPP_REAL *)args[0].data + (opp_p2c[0] * 3),
                (OPP_REAL *)args[1].data + (n * 3),
                (OPP_REAL *)args[2].data + (n * 3)
            );
        }
    }
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("calculate_new_pos_vel_kernel");
}

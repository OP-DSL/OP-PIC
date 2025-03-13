
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k9 {
inline void get_final_max_values_kernel(
    const double* n_charge_den,
    double* max_n_charge_den,
    const double* n_pot,
    double* max_n_pot
) {
    *max_n_charge_den = ((abs(*n_charge_den) > *max_n_charge_den) ? (abs(*n_charge_den)) : (*max_n_charge_den));
    *max_n_pot = ((*n_pot > *max_n_pot) ? (*n_pot) : (*max_n_pot));
}
}

void opp_par_loop_all__get_final_max_values_kernel(opp_set set,
    opp_arg arg0, // n_charge_den | OPP_READ
    opp_arg arg1, // | OPP_MAX
    opp_arg arg2, // n_potential | OPP_READ
    opp_arg arg3 // | OPP_MAX
) 
{
    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    opp_profiler->start("get_final_max_values_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__get_final_max_values_kernel set_size %d", set->size);

    const int nthreads = omp_get_max_threads(); 

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
 
    opp_mpi_halo_wait_all(nargs, args);    

    OPP_REAL arg1_l[nthreads * 1];
    for (int thr = 0; thr < nthreads; thr++)
        for (int d = 0; d < 1; d++)
            arg1_l[1 * thr + d] = OPP_REAL_ZERO;


    OPP_REAL arg3_l[nthreads * 1];
    for (int thr = 0; thr < nthreads; thr++)
        for (int d = 0; d < 1; d++)
            arg3_l[1 * thr + d] = OPP_REAL_ZERO;

  
    
  
 
    #pragma omp parallel for 
    for (int thr = 0; thr < nthreads; thr++)
    {
        const size_t start  = ((size_t)iter_size * thr) / nthreads;
        const size_t finish = ((size_t)iter_size * (thr+1)) / nthreads;
      
        for (size_t n = start; n < finish; n++)
        { 
            opp_k9::get_final_max_values_kernel(
                (const OPP_REAL *)args[0].data + (n * 1),
                arg1_l + (1 * thr),
                (const OPP_REAL *)args[2].data + (n * 1),
                arg3_l + (1 * thr)
            );
        }
    }
    
    for (int thr = 0; thr < nthreads; thr++) 
        for (int d = 0; d < 1; d++)
            ((OPP_REAL*)args[1].data)[d] = MAX(((OPP_REAL *)args[1].data)[d], arg1_l[1 * thr + d]);
#ifdef USE_MPI 
    opp_mpi_reduce(&args[1], (OPP_REAL *)args[1].data);
#endif
    
    for (int thr = 0; thr < nthreads; thr++) 
        for (int d = 0; d < 1; d++)
            ((OPP_REAL*)args[3].data)[d] = MAX(((OPP_REAL *)args[3].data)[d], arg3_l[1 * thr + d]);
#ifdef USE_MPI 
    opp_mpi_reduce(&args[3], (OPP_REAL *)args[3].data);
#endif
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("get_final_max_values_kernel");
}

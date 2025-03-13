
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k9 {
inline void get_max_x_values_kernel(
    const double* cell_j,
    double* max_j,
    const double* cell_e,
    double* max_e,
    const double* cell_b,
    double* max_b)
{
    *max_j = ((*cell_j > *max_j) ? (*cell_j) : (*max_j));

    *max_e = ((*cell_e > *max_e) ? (*cell_e) : (*max_e));

    *max_b = ((*cell_b > *max_b) ? (*cell_b) : (*max_b));
}
}

void opp_par_loop_all__get_max_x_values_kernel(opp_set set,
    opp_arg arg0, // c_j | OPP_READ
    opp_arg arg1, // | OPP_MAX
    opp_arg arg2, // c_e | OPP_READ
    opp_arg arg3, // | OPP_MAX
    opp_arg arg4, // c_b | OPP_READ
    opp_arg arg5 // | OPP_MAX
) 
{
    const int nargs = 6;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;

    opp_profiler->start("get_max_x_values_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__get_max_x_values_kernel set_size %d", set->size);

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


    OPP_REAL arg5_l[nthreads * 1];
    for (int thr = 0; thr < nthreads; thr++)
        for (int d = 0; d < 1; d++)
            arg5_l[1 * thr + d] = OPP_REAL_ZERO;

  
    
  
 
    #pragma omp parallel for 
    for (int thr = 0; thr < nthreads; thr++)
    {
        const size_t start  = ((size_t)iter_size * thr) / nthreads;
        const size_t finish = ((size_t)iter_size * (thr+1)) / nthreads;
      
        for (size_t n = start; n < finish; n++)
        { 
            opp_k9::get_max_x_values_kernel(
                (const OPP_REAL *)args[0].data + (n * 3),
                arg1_l + (1 * thr),
                (const OPP_REAL *)args[2].data + (n * 3),
                arg3_l + (1 * thr),
                (const OPP_REAL *)args[4].data + (n * 3),
                arg5_l + (1 * thr)
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
    
    for (int thr = 0; thr < nthreads; thr++) 
        for (int d = 0; d < 1; d++)
            ((OPP_REAL*)args[5].data)[d] = MAX(((OPP_REAL *)args[5].data)[d], arg5_l[1 * thr + d]);
#ifdef USE_MPI 
    opp_mpi_reduce(&args[5], (OPP_REAL *)args[5].data);
#endif
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("get_max_x_values_kernel");
}

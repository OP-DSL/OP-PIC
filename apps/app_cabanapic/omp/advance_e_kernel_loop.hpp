
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k7 {
enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

inline void advance_e_kernel (
    const double* cell_x_b,
    const double* cell_y_b,
    const double* cell_z_b,
    const double* cell0_b,
    const double* cell0_j,
    double* cell0_e,
    const int* iter_adv_e)
{
    if (iter_adv_e[0] == 1)
    {
        cell0_e[Dim::x] += ( - CONST_dt_eps0[0] * cell0_j[Dim::x] ) +
            ( CONST_p[Dim::y] * (cell0_b[Dim::z] - cell_y_b[Dim::z]) -
            CONST_p[Dim::z] * (cell0_b[Dim::y] - cell_z_b[Dim::y]) );

        cell0_e[Dim::y] += ( - CONST_dt_eps0[0] * cell0_j[Dim::y] ) +
            ( CONST_p[Dim::z] * (cell0_b[Dim::x] - cell_z_b[Dim::x]) -
            CONST_p[Dim::x] * (cell0_b[Dim::z] - cell_x_b[Dim::z]) );

        cell0_e[Dim::z] += ( - CONST_dt_eps0[0] * cell0_j[Dim::z] ) +
            ( CONST_p[Dim::x] * (cell0_b[Dim::y] - cell_x_b[Dim::y]) -
            CONST_p[Dim::y] * (cell0_b[Dim::x] - cell_y_b[Dim::x]) );
    }
}
}

void opp_par_loop_all__advance_e_kernel(opp_set set,
    opp_arg arg0, // c_b | OPP_READ
    opp_arg arg1, // c_b | OPP_READ
    opp_arg arg2, // c_b | OPP_READ
    opp_arg arg3, // c_b | OPP_READ
    opp_arg arg4, // c_j | OPP_READ
    opp_arg arg5, // c_e | OPP_INC
    opp_arg arg6 // c_mask_right | OPP_READ
) 
{
    const int nargs = 7;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;
    args[6] = arg6;

    opp_profiler->start("advance_e_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__advance_e_kernel set_size %d", set->size);

    const int nthreads = omp_get_max_threads(); 

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
 
    opp_mpi_halo_wait_all(nargs, args);    
  
    
  
 
    #pragma omp parallel for 
    for (int thr = 0; thr < nthreads; thr++)
    {
        const size_t start  = ((size_t)iter_size * thr) / nthreads;
        const size_t finish = ((size_t)iter_size * (thr+1)) / nthreads;
      
        for (size_t n = start; n < finish; n++)
        { 
            const OPP_INT *map0 = args[0].map_data + (n * 12);
   
            opp_k7::advance_e_kernel(
                (const OPP_REAL *)args[0].data + (map0[2] * 3),
                (const OPP_REAL *)args[1].data + (map0[4] * 3),
                (const OPP_REAL *)args[2].data + (map0[5] * 3),
                (const OPP_REAL *)args[3].data + (n * 3),
                (const OPP_REAL *)args[4].data + (n * 3),
                (OPP_REAL *)args[5].data + (n * 3),
                (const OPP_INT *)args[6].data + (n * 1)
            );
        }
    }
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("advance_e_kernel");
}

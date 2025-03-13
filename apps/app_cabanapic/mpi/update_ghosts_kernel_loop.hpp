
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k6 {
inline void update_ghosts_kernel(
    const int* c_mask_ug,
    const double* from_cell,
    double* to_cell,
    const int* m_idx,
    const int* dim)
{
    if (c_mask_ug[*m_idx] == 1)
        to_cell[*dim] += from_cell[*dim];
}
}

void opp_par_loop_all__update_ghosts_kernel(opp_set set,
    opp_arg arg0, // c_mask_ug | OPP_READ
    opp_arg arg1, // c_j | OPP_READ
    opp_arg arg2, // c_j | OPP_INC
    opp_arg arg3, // | OPP_READ
    opp_arg arg4 // | OPP_READ
) 
{
    const int nargs = 5;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;

    opp_profiler->start("update_ghosts_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__update_ghosts_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
 
    for (int n = 0; n < iter_size; ++n) 
    {
        if (n == set->core_size)
            opp_mpi_halo_wait_all(nargs, args);

        const OPP_INT *map0 = args[2].map_data + (n * 1);
   
        opp_k6::update_ghosts_kernel(
            (const OPP_INT *)args[0].data + (n * 6),
            (const OPP_REAL *)args[1].data + (n * 3),
            (OPP_REAL *)args[2].data + (map0[0] * 3),
            (OPP_INT *)args[3].data,
            (OPP_INT *)args[4].data
        );
    }
    if (iter_size == 0 || iter_size == set->core_size)
        opp_mpi_halo_wait_all(nargs, args);
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("update_ghosts_kernel");
}

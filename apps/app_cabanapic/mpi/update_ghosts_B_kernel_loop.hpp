
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k5 {
enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

inline void update_ghosts_B_kernel(
    const int* c_mask_ugb,
    const double* from_cell,
    double* to_cell,
    const int* m_idx)
{
    if (c_mask_ugb[*m_idx] == 1)
    {
        to_cell[Dim::x] = from_cell[Dim::x];
        to_cell[Dim::y] = from_cell[Dim::y];
        to_cell[Dim::z] = from_cell[Dim::z];
    }
}
}

void opp_par_loop_all__update_ghosts_B_kernel(opp_set set,
    opp_arg arg0, // c_mask_ugb | OPP_READ
    opp_arg arg1, // c_b | OPP_READ
    opp_arg arg2, // c_b | OPP_WRITE
    opp_arg arg3 // | OPP_READ
) 
{
    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    opp_profiler->start("update_ghosts_B_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__update_ghosts_B_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
 
    for (int n = 0; n < iter_size; ++n) 
    {
        if (n == set->core_size)
            opp_mpi_halo_wait_all(nargs, args);

        const OPP_INT *map0 = args[2].map_data + (n * 1);
   
        opp_k5::update_ghosts_B_kernel(
            (const OPP_INT *)args[0].data + (n * 6),
            (const OPP_REAL *)args[1].data + (n * 3),
            (OPP_REAL *)args[2].data + (map0[0] * 3),
            (OPP_INT *)args[3].data
        );
    }
    if (iter_size == 0 || iter_size == set->core_size)
        opp_mpi_halo_wait_all(nargs, args);
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("update_ghosts_B_kernel");
}

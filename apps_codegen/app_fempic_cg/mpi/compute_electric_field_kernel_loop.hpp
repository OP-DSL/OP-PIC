
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k7 {
inline void compute_electric_field_kernel(
    double *cell_electric_field,
    const double *cell_shape_deriv,
    const double *node_potential0,
    const double *node_potential1,
    const double *node_potential2,
    const double *node_potential3
) {
    for (int dim = 0; dim < 3; dim++) {
        const double c1 = (cell_shape_deriv[0 * 3 + dim] * (*node_potential0));
        const double c2 = (cell_shape_deriv[1 * 3 + dim] * (*node_potential1));
        const double c3 = (cell_shape_deriv[2 * 3 + dim] * (*node_potential2));
        const double c4 = (cell_shape_deriv[3 * 3 + dim] * (*node_potential3));

        cell_electric_field[dim] -= (c1 + c2 + c3 + c4);
    }
}
}

void opp_par_loop_all__compute_electric_field_kernel(opp_set set,
    opp_arg arg0, // c_ef | OPP_INC
    opp_arg arg1, // c_sd | OPP_READ
    opp_arg arg2, // n_potential | OPP_READ
    opp_arg arg3, // n_potential | OPP_READ
    opp_arg arg4, // n_potential | OPP_READ
    opp_arg arg5 // n_potential | OPP_READ
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

    opp_profiler->start("compute_electric_field_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__compute_electric_field_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
 
    for (int n = 0; n < iter_size; ++n) 
    {
        if (n == set->core_size)
            opp_mpi_halo_wait_all(nargs, args);

        const OPP_INT *map0 = args[2].map_data + (n * 4);
   
        opp_k7::compute_electric_field_kernel(
            (OPP_REAL *)args[0].data + (n * 3),
            (const OPP_REAL *)args[1].data + (n * 12),
            (const OPP_REAL *)args[2].data + (map0[0] * 1),
            (const OPP_REAL *)args[3].data + (map0[1] * 1),
            (const OPP_REAL *)args[4].data + (map0[2] * 1),
            (const OPP_REAL *)args[5].data + (map0[3] * 1)
        );
    }
    if (iter_size == 0 || iter_size == set->core_size)
        opp_mpi_halo_wait_all(nargs, args);
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("compute_electric_field_kernel");
}

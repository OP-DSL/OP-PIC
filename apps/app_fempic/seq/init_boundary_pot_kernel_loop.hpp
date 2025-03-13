
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k1 {
inline void init_boundary_pot_kernel(
    const int *node_type,
    double *n_bnd_pot
) {
    switch (*node_type) {
        case 2: // INLET:
            *n_bnd_pot = 0; break;
        case 3: // FIXED:
            *n_bnd_pot = -1 * CONST_wall_potential[0]; break;
        default: // NORMAL or OPEN
            *n_bnd_pot = 0; /*default*/
    }
}
}

void opp_par_loop_all__init_boundary_pot_kernel(opp_set set,
    opp_arg arg0, // n_type | OPP_READ
    opp_arg arg1 // n_bnd_pot | OPP_WRITE
) 
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("init_boundary_pot_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__init_boundary_pot_kernel set_size %d", set->size);

    const int iter_size = set->size;

    for (int n = 0; n < iter_size; ++n) 
    {
        opp_k1::init_boundary_pot_kernel(
            (const OPP_INT *)args[0].data + (n * 1),
            (OPP_REAL *)args[1].data + (n * 1)
        );
    }

    opp_profiler->end("init_boundary_pot_kernel");
}

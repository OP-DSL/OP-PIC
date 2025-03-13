
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k4 {
enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

inline void half_advance_b_kernel (
    const double* cell_x_e,
    const double* cell_y_e,
    const double* cell_z_e,
    const double* cell0_e,
    double* cell0_b,
    const int* cell0_ghost)
{
    if (cell0_ghost[0] == 0)
    {
        cell0_b[Dim::x] -= (0.5 * CONST_p[Dim::y] * (cell_y_e[Dim::z] - cell0_e[Dim::z])
                            - 0.5 * CONST_p[Dim::z] * (cell_z_e[Dim::y] - cell0_e[Dim::y]));

        cell0_b[Dim::y] -= (0.5 * CONST_p[Dim::z] * (cell_z_e[Dim::x] - cell0_e[Dim::x])
                            - 0.5 * CONST_p[Dim::x] * (cell_x_e[Dim::z] - cell0_e[Dim::z]));

        cell0_b[Dim::z] -= (0.5 * CONST_p[Dim::x] * (cell_x_e[Dim::y] - cell0_e[Dim::y])
                            - 0.5 * CONST_p[Dim::y] * (cell_y_e[Dim::x] - cell0_e[Dim::x]));
    }
}
}

void opp_par_loop_all__half_advance_b_kernel(opp_set set,
    opp_arg arg0, // c_e | OPP_READ
    opp_arg arg1, // c_e | OPP_READ
    opp_arg arg2, // c_e | OPP_READ
    opp_arg arg3, // c_e | OPP_READ
    opp_arg arg4, // c_b | OPP_INC
    opp_arg arg5 // c_mask_ghost | OPP_READ
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

    opp_profiler->start("half_advance_b_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__half_advance_b_kernel set_size %d", set->size);

    const int iter_size = set->size;

    for (int n = 0; n < iter_size; ++n) 
    {
        const OPP_INT *map0 = args[0].map_data + (n * 12);
   
        opp_k4::half_advance_b_kernel(
            (const OPP_REAL *)args[0].data + (map0[9] * 3),
            (const OPP_REAL *)args[1].data + (map0[7] * 3),
            (const OPP_REAL *)args[2].data + (map0[6] * 3),
            (const OPP_REAL *)args[3].data + (n * 3),
            (OPP_REAL *)args[4].data + (n * 3),
            (const OPP_INT *)args[5].data + (n * 1)
        );
    }

    opp_profiler->end("half_advance_b_kernel");
}

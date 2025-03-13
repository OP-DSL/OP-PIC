
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k3 {
enum CellAcc {
    jfx = 0 * 4,
    jfy = 1 * 4,
    jfz = 2 * 4,
};

enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

inline void accumulate_current_to_cells_kernel(
        double* cell0_j,
        const double* cell0_acc,
        const double* cell_xd_acc,
        const double* cell_yd_acc,
        const double* cell_zd_acc,
        const double* cell_xyd_acc,
        const double* cell_yzd_acc,
        const double* cell_xzd_acc,
        const int* iter_acc)
{
    if (iter_acc[0] == 1)
    {
        cell0_j[Dim::x] = CONST_acc_coef[Dim::x] * (cell0_acc[CellAcc::jfx + 0] +
                                                    cell_yd_acc[CellAcc::jfx + 1] +
                                                    cell_zd_acc[CellAcc::jfx + 2] +
                                                    cell_yzd_acc[CellAcc::jfx + 3]);

        cell0_j[Dim::y] = CONST_acc_coef[Dim::y] * (cell0_acc[CellAcc::jfy + 0] +
                                                    cell_zd_acc[CellAcc::jfy + 1] +
                                                    cell_xd_acc[CellAcc::jfy + 2] +
                                                    cell_xzd_acc[CellAcc::jfy + 3]);

        cell0_j[Dim::z] = CONST_acc_coef[Dim::z] * (cell0_acc[CellAcc::jfz + 0] +
                                                    cell_xd_acc[CellAcc::jfz + 1] +
                                                    cell_yd_acc[CellAcc::jfz + 2] +
                                                    cell_xyd_acc[CellAcc::jfz + 3]);
    }
}
}

void opp_par_loop_all__accumulate_current_to_cells_kernel(opp_set set,
    opp_arg arg0, // c_j | OPP_WRITE
    opp_arg arg1, // c_acc | OPP_READ
    opp_arg arg2, // c_acc | OPP_READ
    opp_arg arg3, // c_acc | OPP_READ
    opp_arg arg4, // c_acc | OPP_READ
    opp_arg arg5, // c_acc | OPP_READ
    opp_arg arg6, // c_acc | OPP_READ
    opp_arg arg7, // c_acc | OPP_READ
    opp_arg arg8 // c_mask_right | OPP_READ
) 
{
    const int nargs = 9;
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

    opp_profiler->start("accumulate_current_to_cells_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__accumulate_current_to_cells_kernel set_size %d", set->size);

    const int iter_size = set->size;

    for (int n = 0; n < iter_size; ++n) 
    {
        const OPP_INT *map0 = args[2].map_data + (n * 12);
   
        opp_k3::accumulate_current_to_cells_kernel(
            (OPP_REAL *)args[0].data + (n * 3),
            (const OPP_REAL *)args[1].data + (n * 12),
            (const OPP_REAL *)args[2].data + (map0[2] * 12),
            (const OPP_REAL *)args[3].data + (map0[4] * 12),
            (const OPP_REAL *)args[4].data + (map0[5] * 12),
            (const OPP_REAL *)args[5].data + (map0[0] * 12),
            (const OPP_REAL *)args[6].data + (map0[3] * 12),
            (const OPP_REAL *)args[7].data + (map0[1] * 12),
            (const OPP_INT *)args[8].data + (n * 1)
        );
    }

    opp_profiler->end("accumulate_current_to_cells_kernel");
}

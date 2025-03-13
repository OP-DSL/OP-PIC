
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k1 {
enum CellInterp {
    ex = 0,
    dexdy,
    dexdz,
    d2exdydz,
    ey,
    deydz,
    deydx,
    d2eydzdx,
    ez,
    dezdx,
    dezdy,
    d2ezdxdy,
    cbx,
    dcbxdx,
    cby,
    dcbydy,
    cbz,
    dcbzdz,
};

enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

inline void interpolate_mesh_fields_kernel(
    const double* cell0_e,
    const double* cell0_b,
    const double* cell_x_e,
    const double* cell_y_e,
    const double* cell_z_e,
    const double* cell_xy_e,
    const double* cell_yz_e,
    const double* cell_xz_e,
    const double* cell_x_b,
    const double* cell_y_b,
    const double* cell_z_b,
    double* cell0_interp,
    const int* cell0_ghost)
{
    if (cell0_ghost[0] == 0)
    {
        double w0 = 0.0, w1 = 0.0, w2 = 0.0, w3 = 0.0;

        // ex interpolation coefficients
        w0 = cell0_e[Dim::x];                       // pf0->ex;
        w1 = cell_y_e[Dim::x];                      // pfy->ex;
        w2 = cell_z_e[Dim::x];                      // pfz->ex;
        w3 = cell_yz_e[Dim::x];                     // pfyz->ex;

        cell0_interp[CellInterp::ex]       = (1.0 / 4.0)*( (w3 + w0) + (w1 + w2) );
        cell0_interp[CellInterp::dexdy]    = (1.0 / 4.0)*( (w3 - w0) + (w1 - w2) );
        cell0_interp[CellInterp::dexdz]    = (1.0 / 4.0)*( (w3 - w0) - (w1 - w2) );
        cell0_interp[CellInterp::d2exdydz] = (1.0 / 4.0)*( (w3 + w0) - (w1 + w2) );

        // ey interpolation coefficients
        w0 = cell0_e[Dim::y];
        w1 = cell_z_e[Dim::y];                       // pfz->ey;
        w2 = cell_x_e[Dim::y];                       // pfx->ey;
        w3 = cell_xz_e[Dim::y];                      // pfzx->ey;

        cell0_interp[CellInterp::ey]       = (1.0 / 4.0)*( (w3 + w0) + (w1 + w2) );
        cell0_interp[CellInterp::deydz]    = (1.0 / 4.0)*( (w3 - w0) + (w1 - w2) );
        cell0_interp[CellInterp::deydx]    = (1.0 / 4.0)*( (w3 - w0) - (w1 - w2) );
        cell0_interp[CellInterp::d2eydzdx] = (1.0 / 4.0)*( (w3 + w0) - (w1 + w2) );

        // ez interpolation coefficients
        w0 = cell0_e[Dim::z];                       // pf0->ez;
        w1 = cell_x_e[Dim::z];                      // pfx->ez;
        w2 = cell_y_e[Dim::z];                      // pfy->ez;
        w3 = cell_xy_e[Dim::z];                     // pfxy->ez;

        cell0_interp[CellInterp::ez]       = (1.0 / 4.0)*( (w3 + w0) + (w1 + w2) );
        cell0_interp[CellInterp::dezdx]    = (1.0 / 4.0)*( (w3 - w0) + (w1 - w2) );
        cell0_interp[CellInterp::dezdy]    = (1.0 / 4.0)*( (w3 - w0) - (w1 - w2) );
        cell0_interp[CellInterp::d2ezdxdy] = (1.0 / 4.0)*( (w3 + w0) - (w1 + w2) );

        // bx interpolation coefficients
        w0 = cell0_b[Dim::x];                      // pf0->cbx;
        w1 = cell_x_b[Dim::x];                     // pfx->cbx;
        cell0_interp[CellInterp::cbx]    = (1.0 / 2.0)*( w1 + w0 );
        cell0_interp[CellInterp::dcbxdx] = (1.0 / 2.0)*( w1 - w0 );

        // by interpolation coefficients
        w0 = cell0_b[Dim::y];                      // pf0->cby;
        w1 = cell_y_b[Dim::y];                     // pfy->cby;
        cell0_interp[CellInterp::cby]    = (1.0 / 2.0)*( w1 + w0 );
        cell0_interp[CellInterp::dcbydy] = (1.0 / 2.0)*( w1 - w0 );

        // bz interpolation coefficients
        w0 = cell0_b[Dim::z];                      // pf0->cbz;
        w1 = cell_z_b[Dim::z];                     // pfz->cbz;
        cell0_interp[CellInterp::cbz]    = (1.0 / 2.0)*( w1 + w0 );
        cell0_interp[CellInterp::dcbzdz] = (1.0 / 2.0)*( w1 - w0 );
    }
}
}

void opp_par_loop_all__interpolate_mesh_fields_kernel(opp_set set,
    opp_arg arg0, // c_e | OPP_READ
    opp_arg arg1, // c_b | OPP_READ
    opp_arg arg2, // c_e | OPP_READ
    opp_arg arg3, // c_e | OPP_READ
    opp_arg arg4, // c_e | OPP_READ
    opp_arg arg5, // c_e | OPP_READ
    opp_arg arg6, // c_e | OPP_READ
    opp_arg arg7, // c_e | OPP_READ
    opp_arg arg8, // c_b | OPP_READ
    opp_arg arg9, // c_b | OPP_READ
    opp_arg arg10, // c_b | OPP_READ
    opp_arg arg11, // c_interp | OPP_WRITE
    opp_arg arg12 // c_mask_ghost | OPP_READ
) 
{
    const int nargs = 13;
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
    args[9] = arg9;
    args[10] = arg10;
    args[11] = arg11;
    args[12] = arg12;

    opp_profiler->start("interpolate_mesh_fields_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__interpolate_mesh_fields_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
 
    for (int n = 0; n < iter_size; ++n) 
    {
        if (n == set->core_size)
            opp_mpi_halo_wait_all(nargs, args);

        const OPP_INT *map0 = args[2].map_data + (n * 12);
   
        opp_k1::interpolate_mesh_fields_kernel(
            (const OPP_REAL *)args[0].data + (n * 3),
            (const OPP_REAL *)args[1].data + (n * 3),
            (const OPP_REAL *)args[2].data + (map0[9] * 3),
            (const OPP_REAL *)args[3].data + (map0[7] * 3),
            (const OPP_REAL *)args[4].data + (map0[6] * 3),
            (const OPP_REAL *)args[5].data + (map0[11] * 3),
            (const OPP_REAL *)args[6].data + (map0[8] * 3),
            (const OPP_REAL *)args[7].data + (map0[10] * 3),
            (const OPP_REAL *)args[8].data + (map0[9] * 3),
            (const OPP_REAL *)args[9].data + (map0[7] * 3),
            (const OPP_REAL *)args[10].data + (map0[6] * 3),
            (OPP_REAL *)args[11].data + (n * 18),
            (const OPP_INT *)args[12].data + (n * 1)
        );
    }
    if (iter_size == 0 || iter_size == set->core_size)
        opp_mpi_halo_wait_all(nargs, args);
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("interpolate_mesh_fields_kernel");
}

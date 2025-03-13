
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k1_dat0_stride = -1;
OPP_INT opp_k1_dat1_stride = -1;
OPP_INT opp_k1_dat2_stride = -1;
OPP_INT opp_k1_dat3_stride = -1;
OPP_INT opp_k1_map0_stride = -1;

__constant__ OPP_INT opp_k1_dat0_stride_d;
__constant__ OPP_INT opp_k1_dat1_stride_d;
__constant__ OPP_INT opp_k1_dat2_stride_d;
__constant__ OPP_INT opp_k1_dat3_stride_d;
__constant__ OPP_INT opp_k1_map0_stride_d;



namespace opp_k1 {
enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

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

__device__ inline void interpolate_mesh_fields_kernel(
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
    if (cell0_ghost[(0) * opp_k1_dat3_stride_d] == 0)
    {
        double w0 = 0.0, w1 = 0.0, w2 = 0.0, w3 = 0.0;

        // ex interpolation coefficients
        w0 = cell0_e[(Dim::x) * opp_k1_dat0_stride_d];                       // pf0->ex;
        w1 = cell_y_e[(Dim::x) * opp_k1_dat0_stride_d];                      // pfy->ex;
        w2 = cell_z_e[(Dim::x) * opp_k1_dat0_stride_d];                      // pfz->ex;
        w3 = cell_yz_e[(Dim::x) * opp_k1_dat0_stride_d];                     // pfyz->ex;

        cell0_interp[(CellInterp::ex) * opp_k1_dat2_stride_d]       = (1.0 / 4.0)*( (w3 + w0) + (w1 + w2) );
        cell0_interp[(CellInterp::dexdy) * opp_k1_dat2_stride_d]    = (1.0 / 4.0)*( (w3 - w0) + (w1 - w2) );
        cell0_interp[(CellInterp::dexdz) * opp_k1_dat2_stride_d]    = (1.0 / 4.0)*( (w3 - w0) - (w1 - w2) );
        cell0_interp[(CellInterp::d2exdydz) * opp_k1_dat2_stride_d] = (1.0 / 4.0)*( (w3 + w0) - (w1 + w2) );

        // ey interpolation coefficients
        w0 = cell0_e[(Dim::y) * opp_k1_dat0_stride_d];
        w1 = cell_z_e[(Dim::y) * opp_k1_dat0_stride_d];                       // pfz->ey;
        w2 = cell_x_e[(Dim::y) * opp_k1_dat0_stride_d];                       // pfx->ey;
        w3 = cell_xz_e[(Dim::y) * opp_k1_dat0_stride_d];                      // pfzx->ey;

        cell0_interp[(CellInterp::ey) * opp_k1_dat2_stride_d]       = (1.0 / 4.0)*( (w3 + w0) + (w1 + w2) );
        cell0_interp[(CellInterp::deydz) * opp_k1_dat2_stride_d]    = (1.0 / 4.0)*( (w3 - w0) + (w1 - w2) );
        cell0_interp[(CellInterp::deydx) * opp_k1_dat2_stride_d]    = (1.0 / 4.0)*( (w3 - w0) - (w1 - w2) );
        cell0_interp[(CellInterp::d2eydzdx) * opp_k1_dat2_stride_d] = (1.0 / 4.0)*( (w3 + w0) - (w1 + w2) );

        // ez interpolation coefficients
        w0 = cell0_e[(Dim::z) * opp_k1_dat0_stride_d];                       // pf0->ez;
        w1 = cell_x_e[(Dim::z) * opp_k1_dat0_stride_d];                      // pfx->ez;
        w2 = cell_y_e[(Dim::z) * opp_k1_dat0_stride_d];                      // pfy->ez;
        w3 = cell_xy_e[(Dim::z) * opp_k1_dat0_stride_d];                     // pfxy->ez;

        cell0_interp[(CellInterp::ez) * opp_k1_dat2_stride_d]       = (1.0 / 4.0)*( (w3 + w0) + (w1 + w2) );
        cell0_interp[(CellInterp::dezdx) * opp_k1_dat2_stride_d]    = (1.0 / 4.0)*( (w3 - w0) + (w1 - w2) );
        cell0_interp[(CellInterp::dezdy) * opp_k1_dat2_stride_d]    = (1.0 / 4.0)*( (w3 - w0) - (w1 - w2) );
        cell0_interp[(CellInterp::d2ezdxdy) * opp_k1_dat2_stride_d] = (1.0 / 4.0)*( (w3 + w0) - (w1 + w2) );

        // bx interpolation coefficients
        w0 = cell0_b[(Dim::x) * opp_k1_dat1_stride_d];                      // pf0->cbx;
        w1 = cell_x_b[(Dim::x) * opp_k1_dat1_stride_d];                     // pfx->cbx;
        cell0_interp[(CellInterp::cbx) * opp_k1_dat2_stride_d]    = (1.0 / 2.0)*( w1 + w0 );
        cell0_interp[(CellInterp::dcbxdx) * opp_k1_dat2_stride_d] = (1.0 / 2.0)*( w1 - w0 );

        // by interpolation coefficients
        w0 = cell0_b[(Dim::y) * opp_k1_dat1_stride_d];                      // pf0->cby;
        w1 = cell_y_b[(Dim::y) * opp_k1_dat1_stride_d];                     // pfy->cby;
        cell0_interp[(CellInterp::cby) * opp_k1_dat2_stride_d]    = (1.0 / 2.0)*( w1 + w0 );
        cell0_interp[(CellInterp::dcbydy) * opp_k1_dat2_stride_d] = (1.0 / 2.0)*( w1 - w0 );

        // bz interpolation coefficients
        w0 = cell0_b[(Dim::z) * opp_k1_dat1_stride_d];                      // pf0->cbz;
        w1 = cell_z_b[(Dim::z) * opp_k1_dat1_stride_d];                     // pfz->cbz;
        cell0_interp[(CellInterp::cbz) * opp_k1_dat2_stride_d]    = (1.0 / 2.0)*( w1 + w0 );
        cell0_interp[(CellInterp::dcbzdz) * opp_k1_dat2_stride_d] = (1.0 / 2.0)*( w1 - w0 );
    }
}

}

//--------------------------------------------------------------
__global__ void opp_dev_interpolate_mesh_fields_kernel(
    const OPP_REAL *__restrict__ dat0,  // c_e
    const OPP_REAL *__restrict__ dat1,  // c_b
    OPP_REAL *__restrict__ dat2,  // c_interp
    const OPP_INT *__restrict__ dat3,  // c_mask_ghost
    const OPP_INT *__restrict__ map0,  // c2c_map
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        opp_k1::interpolate_mesh_fields_kernel(
            dat0 + n, // c_e 
            dat1 + n, // c_b 
            dat0 + map0[n + opp_k1_map0_stride_d * 9], // c_e 
            dat0 + map0[n + opp_k1_map0_stride_d * 7], // c_e 
            dat0 + map0[n + opp_k1_map0_stride_d * 6], // c_e 
            dat0 + map0[n + opp_k1_map0_stride_d * 11], // c_e 
            dat0 + map0[n + opp_k1_map0_stride_d * 8], // c_e 
            dat0 + map0[n + opp_k1_map0_stride_d * 10], // c_e 
            dat1 + map0[n + opp_k1_map0_stride_d * 9], // c_b 
            dat1 + map0[n + opp_k1_map0_stride_d * 7], // c_b 
            dat1 + map0[n + opp_k1_map0_stride_d * 6], // c_b 
            dat2 + n, // c_interp 
            dat3 + n // c_mask_ghost 
        );
    }
    
}

//--------------------------------------------------------------
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

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k1_dat0_stride_d, &opp_k1_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k1_dat1_stride_d, &opp_k1_dat1_stride, &(args[1].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k1_dat2_stride_d, &opp_k1_dat2_stride, &(args[11].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k1_dat3_stride_d, &opp_k1_dat3_stride, &(args[12].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k1_map0_stride_d, &opp_k1_map0_stride, &(args[2].size), 1);

#ifdef OPP_BLOCK_SIZE_1
    const int block_size = OPP_BLOCK_SIZE_1;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    int num_blocks = 200;
    if (iter_size > 0) 
    {
        const OPP_INT start = 0;
        const OPP_INT end = iter_size;
        num_blocks = (end - start - 1) / block_size + 1;

        {
            opp_dev_interpolate_mesh_fields_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // c_e
                (OPP_REAL *)args[1].data_d,     // c_b
                (OPP_REAL *)args[11].data_d,     // c_interp
                (OPP_INT *)args[12].data_d,     // c_mask_ghost
                args[2].map_data_d,     // c2c_map
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("interpolate_mesh_fields_kernel");
}

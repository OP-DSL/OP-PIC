
/* 
BSD 3-Clause License

Copyright (c) 2022, OP-DSL

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//*********************************************
// AUTO GENERATED CODE
//*********************************************

int imf_OPP_HOST_0 = -1;
int imf_OPP_HOST_1 = -1;
int imf_OPP_HOST_2_MAP_STRIDE = -1;
int imf_OPP_HOST_11 = -1;

__constant__ int imf_OPP_DEV_0;
__constant__ int imf_OPP_DEV_1;
__constant__ int imf_OPP_DEV_2_MAP_STRIDE;
__constant__ int imf_OPP_DEV_11;

//user function
//*************************************************************************************************
__device__ void dev_interpolate_mesh_fields__kernel(
    const OPP_REAL* cell0_e,  
    const OPP_REAL* cell0_b,  
    const OPP_REAL* cell_x_e, 
    const OPP_REAL* cell_y_e, 
    const OPP_REAL* cell_z_e, 
    const OPP_REAL* cell_xy_e,
    const OPP_REAL* cell_yz_e,
    const OPP_REAL* cell_xz_e,
    const OPP_REAL* cell_x_b, 
    const OPP_REAL* cell_y_b, 
    const OPP_REAL* cell_z_b, 
    OPP_REAL* cell0_interp,
    const OPP_INT* cell0_ghost
)
{
    if (cell0_ghost[0] == 0)
    {

        OPP_REAL w0 = 0.0, w1 = 0.0, w2 = 0.0, w3 = 0.0;

        // ex interpolation coefficients
        w0 = cell0_e[imf_OPP_DEV_0 * Dim::x];                       // pf0->ex;
        w1 = cell_y_e[imf_OPP_DEV_0 * Dim::x];                      // pfy->ex;
        w2 = cell_z_e[imf_OPP_DEV_0 * Dim::x];                      // pfz->ex;
        w3 = cell_yz_e[imf_OPP_DEV_0 * Dim::x];                     // pfyz->ex;

        cell0_interp[imf_OPP_DEV_11 * CellInterp::ex]       = CONST_DEV_fourth*( (w3 + w0) + (w1 + w2) );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::dexdy]    = CONST_DEV_fourth*( (w3 - w0) + (w1 - w2) );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::dexdz]    = CONST_DEV_fourth*( (w3 - w0) - (w1 - w2) );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::d2exdydz] = CONST_DEV_fourth*( (w3 + w0) - (w1 + w2) );

        // ey interpolation coefficients
        w0 = cell0_e[imf_OPP_DEV_0 * Dim::y];
        w1 = cell_z_e[imf_OPP_DEV_0 * Dim::y];                       // pfz->ey;
        w2 = cell_x_e[imf_OPP_DEV_0 * Dim::y];                       // pfx->ey;
        w3 = cell_xz_e[imf_OPP_DEV_0 * Dim::y];                      // pfzx->ey;

        cell0_interp[imf_OPP_DEV_11 * CellInterp::ey]       = CONST_DEV_fourth*( (w3 + w0) + (w1 + w2) );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::deydz]    = CONST_DEV_fourth*( (w3 - w0) + (w1 - w2) );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::deydx]    = CONST_DEV_fourth*( (w3 - w0) - (w1 - w2) );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::d2eydzdx] = CONST_DEV_fourth*( (w3 + w0) - (w1 + w2) );

        // ez interpolation coefficients
        w0 = cell0_e[imf_OPP_DEV_0 * Dim::z];                       // pf0->ez;
        w1 = cell_x_e[imf_OPP_DEV_0 * Dim::z];                      // pfx->ez;
        w2 = cell_y_e[imf_OPP_DEV_0 * Dim::z];                      // pfy->ez;
        w3 = cell_xy_e[imf_OPP_DEV_0 * Dim::z];                     // pfxy->ez;

        cell0_interp[imf_OPP_DEV_11 * CellInterp::ez]       = CONST_DEV_fourth*( (w3 + w0) + (w1 + w2) );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::dezdx]    = CONST_DEV_fourth*( (w3 - w0) + (w1 - w2) );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::dezdy]    = CONST_DEV_fourth*( (w3 - w0) - (w1 - w2) );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::d2ezdxdy] = CONST_DEV_fourth*( (w3 + w0) - (w1 + w2) );

        // bx interpolation coefficients
        w0 = cell0_b[imf_OPP_DEV_1 * Dim::x];                      // pf0->cbx;
        w1 = cell_x_b[imf_OPP_DEV_1 * Dim::x];                     // pfx->cbx;
        cell0_interp[imf_OPP_DEV_11 * CellInterp::cbx]    = CONST_DEV_half*( w1 + w0 );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::dcbxdx] = CONST_DEV_half*( w1 - w0 );

        // by interpolation coefficients
        w0 = cell0_b[imf_OPP_DEV_1 * Dim::y];                      // pf0->cby;
        w1 = cell_y_b[imf_OPP_DEV_1 * Dim::y];                     // pfy->cby;
        cell0_interp[imf_OPP_DEV_11 * CellInterp::cby]    = CONST_DEV_half*( w1 + w0 );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::dcbydy] = CONST_DEV_half*( w1 - w0 );

        // bz interpolation coefficients
        w0 = cell0_b[imf_OPP_DEV_1 * Dim::z];                      // pf0->cbz;
        w1 = cell_z_b[imf_OPP_DEV_1 * Dim::z];                     // pfz->cbz;
        cell0_interp[imf_OPP_DEV_11 * CellInterp::cbz]    = CONST_DEV_half*( w1 + w0 );
        cell0_interp[imf_OPP_DEV_11 * CellInterp::dcbzdz] = CONST_DEV_half*( w1 - w0 );
    }
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void dev_interpolate_mesh_fields(
    const OPP_REAL *__restrict__ dir_arg0,
    const OPP_REAL *__restrict__ dir_arg1,
    const OPP_REAL *__restrict__ ind_arg2,
    const OPP_INT *__restrict__ arg2_map,
    const OPP_REAL *__restrict__ ind_arg3,
    const OPP_REAL *__restrict__ ind_arg4,
    const OPP_REAL *__restrict__ ind_arg5,
    const OPP_REAL *__restrict__ ind_arg6,
    const OPP_REAL *__restrict__ ind_arg7,
    const OPP_REAL *__restrict__ ind_arg8,
    const OPP_REAL *__restrict__ ind_arg9,
    const OPP_REAL *__restrict__ ind_arg10,
    OPP_REAL *__restrict__ dir_arg11,
    const OPP_INT *__restrict__ dir_arg12,
    const int start,
    const int end) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        const int map_2idx  = arg2_map[n + imf_OPP_DEV_2_MAP_STRIDE * CellMap::xu_y_z];
        const int map_3idx  = arg2_map[n + imf_OPP_DEV_2_MAP_STRIDE * CellMap::x_yu_z];
        const int map_4idx  = arg2_map[n + imf_OPP_DEV_2_MAP_STRIDE * CellMap::x_y_zu];
        const int map_5idx  = arg2_map[n + imf_OPP_DEV_2_MAP_STRIDE * CellMap::xu_yu_z];
        const int map_6idx  = arg2_map[n + imf_OPP_DEV_2_MAP_STRIDE * CellMap::x_yu_zu];
        const int map_7idx  = arg2_map[n + imf_OPP_DEV_2_MAP_STRIDE * CellMap::xu_y_zu];
        const int map_8idx  = arg2_map[n + imf_OPP_DEV_2_MAP_STRIDE * CellMap::xu_y_z];
        const int map_9idx  = arg2_map[n + imf_OPP_DEV_2_MAP_STRIDE * CellMap::x_yu_z];
        const int map_10idx = arg2_map[n + imf_OPP_DEV_2_MAP_STRIDE * CellMap::x_y_zu];

        //user-supplied kernel call
        dev_interpolate_mesh_fields__kernel(
            (dir_arg0 + n),
            (dir_arg1 + n),
            (ind_arg2 + map_2idx),
            (ind_arg3 + map_3idx),
            (ind_arg4 + map_4idx),
            (ind_arg5 + map_5idx),
            (ind_arg6 + map_6idx),
            (ind_arg7 + map_7idx),
            (ind_arg8 + map_8idx),
            (ind_arg9 + map_9idx),
            (ind_arg10 + map_10idx),
            (dir_arg11 + n),
            (dir_arg12 + n)
        );
    }
}

//*************************************************************************************************
void opp_loop_all__interpolate_mesh_fields(
    opp_set set,        // cells_set
    opp_arg arg0,       // cell0_e,        // OPP_READ
    opp_arg arg1,       // cell0_b,        // OPP_READ
    opp_arg arg2,       // cell_x_e,       // OPP_READ
    opp_arg arg3,       // cell_y_e,       // OPP_READ
    opp_arg arg4,       // cell_z_e,       // OPP_READ
    opp_arg arg5,       // cell_yz_e,      // OPP_READ 
    opp_arg arg6,       // cell_xz_e,      // OPP_READ
    opp_arg arg7,       // cell_xy_e,      // OPP_READ
    opp_arg arg8,       // cell_x_b,       // OPP_READ
    opp_arg arg9,       // cell_y_b,       // OPP_READ
    opp_arg arg10,      // cell_z_b        // OPP_READ
    opp_arg arg11,      // cell0_interp    // OPP_WRITE
    opp_arg arg12       // cell0_ghost     // OPP_READ
)
{ 
    
    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__interpolate_mesh_fields set_size %d", set->size);

    opp_profiler->start("Interpolate");

    const int nargs = 13;
    opp_arg args[nargs];

    args[0] = std::move(arg0);
    args[1] = std::move(arg1);
    args[2] = std::move(arg2);
    args[3] = std::move(arg3);
    args[4] = std::move(arg4);
    args[5] = std::move(arg5);
    args[6] = std::move(arg6);
    args[7] = std::move(arg7);
    args[8] = std::move(arg8);
    args[9] = std::move(arg9);
    args[10] = std::move(arg10);
    args[11] = std::move(arg11);
    args[12] = std::move(arg12);

    opp_profiler->start("Interp_Halo");
    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    opp_profiler->end("Interp_Halo");
    
    if (set_size > 0) 
    {
        imf_OPP_HOST_0 = args[0].dat->set->set_capacity;
        imf_OPP_HOST_1 = args[1].dat->set->set_capacity;
        imf_OPP_HOST_2_MAP_STRIDE = args[2].size;
        imf_OPP_HOST_11 = args[11].dat->set->set_capacity;

        cutilSafeCall(cudaMemcpyToSymbol(imf_OPP_DEV_0, 
                                                    &imf_OPP_HOST_0, sizeof(int)));
        cutilSafeCall(cudaMemcpyToSymbol(imf_OPP_DEV_1, 
                                                    &imf_OPP_HOST_1, sizeof(int)));
        cutilSafeCall(cudaMemcpyToSymbol(imf_OPP_DEV_2_MAP_STRIDE, 
                                                    &imf_OPP_HOST_2_MAP_STRIDE, sizeof(int)));
        cutilSafeCall(cudaMemcpyToSymbol(imf_OPP_DEV_11, 
                                                    &imf_OPP_HOST_11, sizeof(int)));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = OPP_gpu_threads_per_block;
            int nblocks = (end - start - 1) / nthread + 1;
            
            dev_interpolate_mesh_fields<<<nblocks, nthread>>>(
                (OPP_REAL*) args[0].data_d,     // cell0_e,        // OPP_READ
                (OPP_REAL*) args[1].data_d,     // cell0_b,        // OPP_READ
                (OPP_REAL*) args[2].data_d,     // cell_x_e,       // OPP_READ
                (OPP_INT *) args[2].map_data_d,
                (OPP_REAL*) args[3].data_d,     // cell_y_e,       // OPP_READ
                (OPP_REAL*) args[4].data_d,     // cell_z_e,       // OPP_READ
                (OPP_REAL*) args[5].data_d,     // cell_yz_e,      // OPP_READ
                (OPP_REAL*) args[6].data_d,     // cell_xz_e,      // OPP_READ
                (OPP_REAL*) args[7].data_d,     // cell_xy_e,      // OPP_READ
                (OPP_REAL*) args[8].data_d,     // cell_x_b,       // OPP_READ
                (OPP_REAL*) args[9].data_d,     // cell_y_b,       // OPP_READ
                (OPP_REAL*) args[10].data_d,    // cell_z_b        // OPP_READ
                (OPP_REAL*) args[11].data_d,    // cell0_interp    // OPP_WRITE
                (OPP_INT*) args[12].data_d,     // cell0_ghost,    // OPP_READ
                start, 
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());

    opp_profiler->end("Interpolate");
}

//*************************************************************************************************
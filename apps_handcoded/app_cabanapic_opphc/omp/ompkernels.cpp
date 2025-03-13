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


#include "opp_omp.h"
#include "../cabana_defs.h"

OPP_INT CONST_c_per_dim[DIM];
OPP_REAL CONST_dt;
OPP_REAL CONST_qsp;
OPP_REAL CONST_cdt_d[DIM];
OPP_REAL CONST_p[DIM];
OPP_REAL CONST_qdt_2mc;
OPP_REAL CONST_dt_eps0;
OPP_REAL CONST_acc_coef[DIM];

//****************************************
void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    if (!strcmp(name,"CONST_c_per_dim"))     std::memcpy(CONST_c_per_dim, data, (size*dim));
    else if (!strcmp(name,"CONST_dt"))       std::memcpy(&CONST_dt, data, (size*dim));
    else if (!strcmp(name,"CONST_qsp"))      std::memcpy(&CONST_qsp, data, (size*dim));
    else if (!strcmp(name,"CONST_cdt_d"))    std::memcpy(&CONST_cdt_d, data, (size*dim));
    else if (!strcmp(name,"CONST_p"))        std::memcpy(&CONST_p, data, (size*dim));
    else if (!strcmp(name,"CONST_qdt_2mc"))  std::memcpy(&CONST_qdt_2mc, data, (size*dim));
    else if (!strcmp(name,"CONST_dt_eps0"))  std::memcpy(&CONST_dt_eps0, data, (size*dim));
    else if (!strcmp(name,"CONST_acc_coef")) std::memcpy(&CONST_acc_coef, data, (size*dim));
    else std::cerr << "error: unknown const name" << std::endl;
}
//****************************************

#include "../kernels.h"

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

    #pragma omp parallel for
    for (int n = 0; n < set->size; n++)
    {
        const int map_2idx  = arg2.map_data[n * arg2.map->dim + CellMap::xu_y_z];
        const int map_3idx  = arg3.map_data[n * arg3.map->dim + CellMap::x_yu_z];
        const int map_4idx  = arg4.map_data[n * arg4.map->dim + CellMap::x_y_zu];
        const int map_5idx  = arg5.map_data[n * arg5.map->dim + CellMap::xu_yu_z];
        const int map_6idx  = arg6.map_data[n * arg6.map->dim + CellMap::x_yu_zu];
        const int map_7idx  = arg7.map_data[n * arg7.map->dim + CellMap::xu_y_zu];
        const int map_8idx  = arg8.map_data[n * arg8.map->dim + CellMap::xu_y_z];
        const int map_9idx  = arg9.map_data[n * arg9.map->dim + CellMap::x_yu_z];
        const int map_10idx = arg10.map_data[n * arg10.map->dim + CellMap::x_y_zu];

        interpolate_mesh_fields_kernel(
            &((double*) arg0.data)[n * arg0.dim],            // cell0_e,        // OPP_READ
            &((double*) arg1.data)[n * arg1.dim],            // cell0_b,        // OPP_READ
            &((double*) arg2.data)[map_2idx * arg2.dim],     // cell_x_e,       // OPP_READ
            &((double*) arg3.data)[map_3idx * arg3.dim],     // cell_y_e,       // OPP_READ
            &((double*) arg4.data)[map_4idx * arg4.dim],     // cell_z_e,       // OPP_READ
            &((double*) arg5.data)[map_5idx * arg5.dim],     // cell_yz_e,      // OPP_READ
            &((double*) arg6.data)[map_6idx * arg6.dim],     // cell_xz_e,      // OPP_READ
            &((double*) arg7.data)[map_7idx * arg7.dim],     // cell_xy_e,      // OPP_READ
            &((double*) arg8.data)[map_8idx * arg8.dim],     // cell_x_b,       // OPP_READ
            &((double*) arg9.data)[map_9idx * arg9.dim],     // cell_y_b,       // OPP_READ
            &((double*) arg10.data)[map_10idx * arg10.dim],  // cell_z_b        // OPP_READ
            &((double*) arg11.data)[n * arg11.dim],          // cell0_interp    // OPP_WRITE
            &((int*)    arg12.data)[n * arg12.dim]           // cell0_ghost,    // OPP_READ
        );
    }

    opp_profiler->end("Interpolate");
}

#include "ompkernelgen.h"

//*************************************************************************************************
void opp_particle_move__move_deposit(
    opp_set set, opp_map c2c_map, opp_dat p2c_map,
    opp_arg arg0,       // part_vel           // OPP_RW
    opp_arg arg1,       // part_pos           // OPP_RW
    opp_arg arg2,       // part_streak_mid    // OPP_RW
    opp_arg arg3,       // part_weight        // OPP_READ
    opp_arg arg4,       // cell_inter         // OPP_READ
    opp_arg arg5        // cell_acc           // OPP_INC
)
{

    if (OPP_DBG) 
        opp_printf("CABANA", "opp_particle_move__move_deposit set_size %d diff %d", set->size, set->diff);

    opp_profiler->start("Move");

    opp_init_particle_move(set, 0, nullptr);

    opp_create_thread_level_data<OPP_REAL>(arg5, 0.0);

    const int nthreads = omp_get_max_threads();
    const int c2c_dim   = c2c_map->dim;

    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        const size_t start  = ((size_t)set->size * thr) / nthreads;
        const size_t finish = ((size_t)set->size * (thr+1)) / nthreads;

        char* arg5_dat_thread_data = (*(arg5.dat->thread_data))[thr]; // checked, all values are zero       
        char move_flag = OPPX_MOVE_DONE;
        bool iter_one_flag = true;

        for (size_t n = start; n < finish; n++)
        { 
            iter_one_flag = true;
            OPP_INT *opp_thr_c2c = nullptr, *opp_thr_p2c = nullptr;

            do
            {
                move_flag = OPPX_MOVE_DONE;
                opp_thr_p2c = &((OPP_INT*) p2c_map->data)[n];   // TODO : remove OPP_INT* after making this into a map
                opp_thr_c2c = &(c2c_map->map)[*opp_thr_p2c * c2c_dim];

                move_deposit_kernel_omp(move_flag, iter_one_flag, opp_thr_c2c, opp_thr_p2c,
                    &((double*)arg0.data)[n * arg0.dim],                    // part_vel 
                    &((double*)arg1.data)[n * arg1.dim],                    // part_pos 
                    &((double*)arg2.data)[n * arg2.dim],                    // part_streak_mid 
                    &((double*)arg3.data)[n * arg3.dim],                    // part_weight 
                    &((double*)arg4.data)[*opp_thr_p2c * arg4.dim],             // cell_interp 
                    &((double*)arg5_dat_thread_data)[*opp_thr_p2c * arg5.dim]   // cell_acc
                );

            } while (opp_check_part_move_status(move_flag, iter_one_flag, *opp_thr_p2c, n, thr)); 
        }
    }

    opp_reduce_thread_level_data<OPP_REAL>(arg5);

    opp_finalize_particle_move(set);

    opp_profiler->end("Move");
}

//*************************************************************************************************
void opp_loop_all__accumulate_current_to_cells(
    opp_set set,     // cells set
    opp_arg arg0,    // cell0_j         // OPP_WRITE
    opp_arg arg1,    // cell0_acc       // OPP_READ
    opp_arg arg2,    // cell_xd_acc     // OPP_READ
    opp_arg arg3,    // cell_yd_acc     // OPP_READ
    opp_arg arg4,    // cell_zd_acc     // OPP_READ
    opp_arg arg5,    // cell_xyd_acc    // OPP_READ
    opp_arg arg6,    // cell_yzd_acc    // OPP_READ
    opp_arg arg7,    // cell_xzd_acc    // OPP_READ
    opp_arg arg8     // iter_acc        // OPP_READ
)
{

    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__accumulate_current_to_cells set_size %d", set->size);

    opp_profiler->start("Acc_Current");

    #pragma omp parallel for
    for (int n = 0; n < set->size; n++)
    {
        const int map_2idx  = arg2.map_data[n * arg2.map->dim + CellMap::xd_y_z];
        const int map_3idx  = arg3.map_data[n * arg3.map->dim + CellMap::x_yd_z];
        const int map_4idx  = arg4.map_data[n * arg4.map->dim + CellMap::x_y_zd];
        const int map_5idx  = arg5.map_data[n * arg5.map->dim + CellMap::xd_yd_z];
        const int map_6idx  = arg6.map_data[n * arg6.map->dim + CellMap::x_yd_zd];
        const int map_7idx  = arg7.map_data[n * arg7.map->dim + CellMap::xd_y_zd];

        accumulate_current_to_cells_kernel(   
            &((double*) arg0.data)[n * arg0.dim],            // cell0_j     
            &((double*) arg1.data)[n * arg1.dim],            // cell0_acc   
            &((double*) arg2.data)[map_2idx * arg2.dim],     // cell_xd_acc 
            &((double*) arg3.data)[map_3idx * arg3.dim],     // cell_yd_acc 
            &((double*) arg4.data)[map_4idx * arg4.dim],     // cell_zd_acc
            &((double*) arg5.data)[map_5idx * arg5.dim],     // cell_xyd_acc 
            &((double*) arg6.data)[map_6idx * arg6.dim],     // cell_yzd_acc 
            &((double*) arg7.data)[map_7idx * arg7.dim],     // cell_xzd_acc
            &((int*)    arg8.data)[n * arg8.dim]             // iter_acc 
        );
    }

    opp_profiler->end("Acc_Current");
}

//*************************************************************************************************
void opp_loop_all__half_advance_b(
    opp_set set,     // cells set
    opp_arg arg0,    // cell_x_e        // OPP_READ
    opp_arg arg1,    // cell_y_e        // OPP_READ
    opp_arg arg2,    // cell_z_e        // OPP_READ
    opp_arg arg3,    // cell0_e         // OPP_READ
    opp_arg arg4,    // cell0_b         // OPP_INC
    opp_arg arg5     // cell0_ghost     // OPP_READ
)
{

    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__half_advance_b set_size %d", set->size);

    opp_profiler->start("HalfAdv_B");

    #pragma omp parallel for
    for (int n = 0; n < set->size; n++)
    {
        const int map_0idx  = arg0.map_data[n * arg0.map->dim + CellMap::xu_y_z];
        const int map_1idx  = arg1.map_data[n * arg1.map->dim + CellMap::x_yu_z];
        const int map_2idx  = arg2.map_data[n * arg2.map->dim + CellMap::x_y_zu];

        half_advance_b_kernel (
            &((double*) arg0.data)[map_0idx * arg0.dim],     // cell_x_e, 
            &((double*) arg1.data)[map_1idx * arg1.dim],     // cell_y_e, 
            &((double*) arg2.data)[map_2idx * arg2.dim],     // cell_z_e, 
            &((double*) arg3.data)[n * arg3.dim],            // cell0_e, 
            &((double*) arg4.data)[n * arg4.dim],            // cell0_b
            &((int*)    arg5.data)[n * arg5.dim]             // cell0_ghost
        );
    }

    opp_profiler->end("HalfAdv_B");
}

//*************************************************************************************************
void opp_loop_all__advance_e(
    opp_set set,     // cells set
    opp_arg arg0,    // cell_x_b        // OPP_READ
    opp_arg arg1,    // cell_y_b        // OPP_READ
    opp_arg arg2,    // cell_z_b        // OPP_READ
    opp_arg arg3,    // cell0_b         // OPP_READ
    opp_arg arg4,    // cell0_j         // OPP_READ
    opp_arg arg5,    // cell0_e         // OPP_INC
    opp_arg arg6     // iter_adv_e      // OPP_READ
)
{

    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__advance_e set_size %d", set->size);

    opp_profiler->start("Adv_E");

    #pragma omp parallel for
    for (int n = 0; n < set->size; n++)
    {
        const int map_0idx  = arg0.map_data[n * arg0.map->dim + CellMap::xd_y_z];
        const int map_1idx  = arg1.map_data[n * arg1.map->dim + CellMap::x_yd_z];
        const int map_2idx  = arg2.map_data[n * arg2.map->dim + CellMap::x_y_zd];

        advance_e_kernel (
            &((double*) arg0.data)[map_0idx * arg0.dim],     // cell_x_b  
            &((double*) arg1.data)[map_1idx * arg1.dim],     // cell_y_b  
            &((double*) arg2.data)[map_2idx * arg2.dim],     // cell_z_b  
            &((double*) arg3.data)[n * arg3.dim],            // cell0_b   
            &((double*) arg4.data)[n * arg4.dim],            // cell0_j   
            &((double*) arg5.data)[n * arg5.dim],            // cell0_e
            &((int*)    arg6.data)[n * arg6.dim]             // iter_adv_e   
        );
    }

    opp_profiler->end("Adv_E");
}

//*************************************************************************************************
void opp_loop_all__get_max_values(
    opp_set set,     // cells set
    opp_arg arg0,    // cell_j       // OPP_READ
    opp_arg arg1,    // max_j        // OPP_MAX
    opp_arg arg2,    // cell_e       // OPP_READ
    opp_arg arg3,    // max_e        // OPP_MAX
    opp_arg arg4,    // cell_b       // OPP_READ
    opp_arg arg5     // max_b        // OPP_MAX
)
{
    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__get_max set_size %d", set->size);

    opp_profiler->start("GetMax");

    const int nargs = 6;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;

    const int set_size = set->size;
    const int nthreads = omp_get_max_threads();

    OPP_REAL arg1_l[nthreads*1];
    for (int thr = 0; thr < nthreads; thr++)
    {
        for (int d = 0; d < 1; d++) // can have multiple dimension defined for global_arg
            arg1_l[1 * thr + d] = ZERO_double;
    }
    OPP_REAL arg3_l[nthreads*1];
    for (int thr = 0; thr < nthreads; thr++)
    {
        for (int d = 0; d < 1; d++) // can have multiple dimension defined for global_arg
            arg3_l[1 * thr + d] = ZERO_double;
    }   
    OPP_REAL arg5_l[nthreads*1];
    for (int thr = 0; thr < nthreads; thr++)
    {
        for (int d = 0; d < 1; d++) // can have multiple dimension defined for global_arg
            arg5_l[1 * thr + d] = ZERO_double;
    }  

    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        size_t start  = ((size_t)set_size * thr) / nthreads;
        size_t finish = ((size_t)set_size * (thr+1)) / nthreads;

        for (size_t n = start; n < finish; n++)
        { 
            get_max_x_values_kernel(
                &((double*) args[0].data)[n * args[0].dim],     // cell_j  
                &((double*) arg1_l)[1 * thr],
                &((double*) args[2].data)[n * args[2].dim],     // cell_e  
                &((double*) arg3_l)[1 * thr],
                &((double*) args[4].data)[n * args[4].dim],     // cell_b  
                &((double*) arg5_l)[1 * thr]
            );
        }
    }

    // combine reduction data -- TODO : Create a MAX API call!
    for (int thr = 0; thr < nthreads; thr++) 
    {
        for (int d = 0; d < 1; d++) // can have multiple dimension defined for global_arg
            ((double*)args[1].data)[d] = (((double*)args[1].data)[d] > arg1_l[1 * thr + d]) ? 
                                                    ((double*)args[1].data)[d] : arg1_l[1 * thr + d];
    }
    for (int thr = 0; thr < nthreads; thr++) 
    {
        for (int d = 0; d < 1; d++) // can have multiple dimension defined for global_arg
            ((double*)args[3].data)[d] = (((double*)args[3].data)[d] > arg3_l[1 * thr + d]) ? 
                                                    ((double*)args[3].data)[d] : arg3_l[1 * thr + d];
    }
    for (int thr = 0; thr < nthreads; thr++) 
    {
        for (int d = 0; d < 1; d++) // can have multiple dimension defined for global_arg
            ((double*)args[5].data)[d] = (((double*)args[5].data)[d] > arg5_l[1 * thr + d]) ? 
                                                    ((double*)args[5].data)[d] : arg5_l[1 * thr + d];
    }

    opp_profiler->end("GetMax");
}

//*************************************************************************************************
void opp_loop_all__update_ghosts_B(
    opp_set set,     // cells set
    opp_arg arg0,    // cell_mask_ugb,       OPP_READ
    opp_arg arg1,    // cell,                OPP_READ
    opp_arg arg2,    // cell, 0, c2cugb_map, OPP_WRITE
    opp_arg arg3     // mask_idx global
)
{
    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__update_ghosts_B set_size %d", set->size);

    opp_profiler->start("UpGhostB");

    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    const int set_size = set->size;
    const int nthreads = omp_get_max_threads();

    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        size_t start  = ((size_t)set_size * thr) / nthreads;
        size_t finish = ((size_t)set_size * (thr+1)) / nthreads;

        for (size_t n = start; n < finish; n++)
        { 
            const int map_0idx  = args[2].map_data[n * args[2].map->dim + 0];

            update_ghosts_B_kernel(
                &((OPP_INT*) args[0].data)[n * args[0].dim],        
                &((OPP_REAL*) args[1].data)[n * args[1].dim],       
                &((OPP_REAL*) args[2].data)[map_0idx * args[2].dim],
                (OPP_INT*) args[3].data
            );
        }
    }

    opp_profiler->end("UpGhostB");   
}

//*************************************************************************************************
void opp_loop_all__update_ghosts(
    opp_set set,     // cells set
    opp_arg arg0,    // cell_mask_ug,       OPP_READ
    opp_arg arg1,    // cell,               OPP_READ
    opp_arg arg2,    // cell, 0, c2cug_map, OPP_INC
    opp_arg arg3,    // mask_idx global
    opp_arg arg4     // dim_idx
)
{
    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__update_ghosts set_size %d", set->size);

    opp_profiler->start("UpGhost");

    const int nargs = 5;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;

    const int set_size = set->size;
    const int nthreads = omp_get_max_threads();

    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        size_t start  = ((size_t)set_size * thr) / nthreads;
        size_t finish = ((size_t)set_size * (thr+1)) / nthreads;

        for (size_t n = start; n < finish; n++)
        { 
            const int map_0idx  = args[2].map_data[n * args[2].map->dim + 0];

            update_ghosts_kernel(
                &((OPP_INT*) args[0].data)[n * args[0].dim],        
                &((OPP_REAL*) args[1].data)[n * args[1].dim],       
                &((OPP_REAL*) args[2].data)[map_0idx * args[2].dim],
                (OPP_INT*) args[3].data,
                (OPP_INT*) args[4].data
            );
        }
    }

    opp_profiler->end("UpGhost");   
}

//*************************************************************************************************
void opp_loop_all__compute_energy(
    opp_set set,     // cells set
    opp_arg arg0,    // cell0_ghost, OPP_READ
    opp_arg arg1,    // cell_field,  OPP_READ
    opp_arg arg2     // energy,      OPP_INC
)
{
    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__compute_energy set_size %d", set->size);

    opp_profiler->start("Energy");

    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    const int set_size = set->size;
    const int nthreads = omp_get_max_threads();

    OPP_REAL arg2_l[nthreads*1];
    for (int thr = 0; thr < nthreads; thr++)
    {
        for (int d = 0; d < 1; d++) // can have multiple dimension defined for global_arg
            arg2_l[1 * thr + d] = ZERO_double;
    }  

    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        size_t start  = ((size_t)set_size * thr) / nthreads;
        size_t finish = ((size_t)set_size * (thr+1)) / nthreads;

        for (size_t n = start; n < finish; n++)
        { 
            compute_energy_kernel(
                &((OPP_INT*)  args[0].data)[n * args[0].dim],        
                &((OPP_REAL*) args[1].data)[n * args[1].dim],       
                (OPP_REAL*) &(arg2_l[1 * thr])
            );
        }
    }

    for (int thr = 0; thr < nthreads; thr++) 
    {
        for (int d = 0; d < 1; d++) // can have multiple dimension defined for global_arg
            ((OPP_REAL*)args[2].data)[d] += arg2_l[1 * thr + d];
    }

    opp_profiler->end("Energy");   
}

//*************************************************************************************************
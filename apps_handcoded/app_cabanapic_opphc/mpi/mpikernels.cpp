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


#include "opp_mpi.h"
#include "../cabana_defs.h"

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
    if (!strcmp(name,"CONST_dt"))            std::memcpy(&CONST_dt, data, (size*dim));
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

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);

    opp_profiler->start("Interpolate_Kernel");
    // #pragma omp simd
    for (int n = 0; n < set_size; n++)
    {
        // if (n == set->core_size) {
        //     opp_mpi_halo_wait_all(nargs, args);
        // }

        const OPP_INT map_2idx  = args[2].map_data[n * args[2].map->dim + CellMap::xu_y_z];
        const OPP_INT map_3idx  = args[3].map_data[n * args[3].map->dim + CellMap::x_yu_z];
        const OPP_INT map_4idx  = args[4].map_data[n * args[4].map->dim + CellMap::x_y_zu];
        const OPP_INT map_5idx  = args[5].map_data[n * args[5].map->dim + CellMap::xu_yu_z];
        const OPP_INT map_6idx  = args[6].map_data[n * args[6].map->dim + CellMap::x_yu_zu];
        const OPP_INT map_7idx  = args[7].map_data[n * args[7].map->dim + CellMap::xu_y_zu];
        const OPP_INT map_8idx  = args[8].map_data[n * args[8].map->dim + CellMap::xu_y_z];
        const OPP_INT map_9idx  = args[9].map_data[n * args[9].map->dim + CellMap::x_yu_z];
        const OPP_INT map_10idx = args[10].map_data[n * args[10].map->dim + CellMap::x_y_zu];

        interpolate_mesh_fields_kernel(
            &((const OPP_REAL*) args[0].data)[n * args[0].dim],            // cell0_e,        // OPP_READ
            &((const OPP_REAL*) args[1].data)[n * args[1].dim],            // cell0_b,        // OPP_READ
            &((const OPP_REAL*) args[2].data)[map_2idx * args[2].dim],     // cell_x_e,       // OPP_READ
            &((const OPP_REAL*) args[3].data)[map_3idx * args[3].dim],     // cell_y_e,       // OPP_READ
            &((const OPP_REAL*) args[4].data)[map_4idx * args[4].dim],     // cell_z_e,       // OPP_READ
            &((const OPP_REAL*) args[5].data)[map_5idx * args[5].dim],     // cell_yz_e,      // OPP_READ
            &((const OPP_REAL*) args[6].data)[map_6idx * args[6].dim],     // cell_xz_e,      // OPP_READ
            &((const OPP_REAL*) args[7].data)[map_7idx * args[7].dim],     // cell_xy_e,      // OPP_READ
            &((const OPP_REAL*) args[8].data)[map_8idx * args[8].dim],     // cell_x_b,       // OPP_READ
            &((const OPP_REAL*) args[9].data)[map_9idx * args[9].dim],     // cell_y_b,       // OPP_READ
            &((const OPP_REAL*) args[10].data)[map_10idx * args[10].dim],  // cell_z_b        // OPP_READ
            &((OPP_REAL*)       args[11].data)[n * args[11].dim],          // cell0_interp    // OPP_WRITE
            &((const OPP_INT*)  args[12].data)[n * args[12].dim]           // cell0_ghost,    // OPP_READ
        );
    }
    opp_profiler->end("Interpolate_Kernel");

    // if (set->size == 0 || set->size == set->core_size) {
    //     opp_mpi_halo_wait_all(nargs, args);
    // }

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("Interpolate");
}

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

    const int nargs = 7;
    opp_arg args[nargs];

    args[0]  = arg0;
    args[1]  = arg1;
    args[2]  = arg2;
    args[3]  = arg3;
    args[4]  = arg4;
    args[5]  = arg5;
    args[6]  = opp_arg_dat(p2c_map, OPP_RW); // required to make dirty or should manually make it dirty

    opp_mpi_halo_exchanges(set, nargs, args);
    // unable to overlap much of computation and communication
    opp_mpi_halo_wait_all(nargs, args); 

    do // iterate until all mpi ranks say, I am done
    {
        opp_init_particle_move(set, nargs, args);
        
        if (OPP_DBG || OPP_comm_iteration > 3)
            opp_printf("MOVE", "iter %d start %d end %d - COUNT=%d", OPP_comm_iteration, 
                        OPP_iter_start, OPP_iter_end, (OPP_iter_end - OPP_iter_start));

        const int c2c_dim   = c2c_map->dim;

        opp_profiler->start("MoveDeposit_Kernel");
        const int start = OPP_iter_start;
        const int end = OPP_iter_end;
        for (int n = start; n < end; n++)
        {
            OPP_MOVE_RESET_FLAGS;

            do
            {
                opp_p2c = &((OPP_INT*) p2c_map->data)[n];   // TODO : remove OPP_INT* after making this into a map
                opp_c2c = &(c2c_map->map)[*opp_p2c * c2c_dim];

                move_deposit_kernel(
                    &((OPP_REAL*)       args[0].data)[n * args[0].dim],         // part_vel 
                    &((OPP_REAL*)       args[1].data)[n * args[1].dim],         // part_pos 
                    &((OPP_REAL*)       args[2].data)[n * args[2].dim],         // part_streak_mid 
                    &((const OPP_REAL*) args[3].data)[n * args[3].dim],         // part_weight 
                    &((const OPP_REAL*) args[4].data)[*opp_p2c * args[4].dim],  // cell_interp 
                    &((OPP_REAL*)       args[5].data)[*opp_p2c * args[5].dim]   // cell_acc
                );

            } while (opp_check_part_move_status(*opp_p2c, n, set->particle_remove_count)); 
        }
        opp_profiler->end("MoveDeposit_Kernel");

    } while (opp_finalize_particle_move(set)); // iterate until all mpi ranks say, I am done

    opp_set_dirtybit(nargs, args);

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

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);

    opp_profiler->start("Acc_Current_Kernel");
    for (int n = 0; n < set_size; n++)
    {
        // if (n == set->core_size) {
        //     opp_mpi_halo_wait_all(nargs, args);
        // }

        const OPP_INT map_2idx  = args[2].map_data[n * args[2].map->dim + CellMap::xd_y_z];
        const OPP_INT map_3idx  = args[3].map_data[n * args[3].map->dim + CellMap::x_yd_z];
        const OPP_INT map_4idx  = args[4].map_data[n * args[4].map->dim + CellMap::x_y_zd];
        const OPP_INT map_5idx  = args[5].map_data[n * args[5].map->dim + CellMap::xd_yd_z];
        const OPP_INT map_6idx  = args[6].map_data[n * args[6].map->dim + CellMap::x_yd_zd];
        const OPP_INT map_7idx  = args[7].map_data[n * args[7].map->dim + CellMap::xd_y_zd];

        accumulate_current_to_cells_kernel(   
            &((OPP_REAL*)       args[0].data)[n * args[0].dim],            // cell0_j     
            &((const OPP_REAL*) args[1].data)[n * args[1].dim],            // cell0_acc   
            &((const OPP_REAL*) args[2].data)[map_2idx * args[2].dim],     // cell_xd_acc 
            &((const OPP_REAL*) args[3].data)[map_3idx * args[3].dim],     // cell_yd_acc 
            &((const OPP_REAL*) args[4].data)[map_4idx * args[4].dim],     // cell_zd_acc
            &((const OPP_REAL*) args[5].data)[map_5idx * args[5].dim],     // cell_xyd_acc 
            &((const OPP_REAL*) args[6].data)[map_6idx * args[6].dim],     // cell_yzd_acc 
            &((const OPP_REAL*) args[7].data)[map_7idx * args[7].dim],     // cell_xzd_acc
            &((const OPP_INT*)  args[8].data)[n * args[8].dim]             // iter_acc 
        );
    }
    opp_profiler->end("Acc_Current_Kernel");

    // if (set->size == 0 || set->size == set->core_size) {
    //     opp_mpi_halo_wait_all(nargs, args);
    // }

    opp_set_dirtybit(nargs, args);

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

    const int nargs = 6;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);

    opp_profiler->start("HalfAdv_B_Kernel");
    for (int n = 0; n < set_size; n++)
    {
        // if (n == set->core_size) {
        //     opp_mpi_halo_wait_all(nargs, args);
        // }

        const OPP_INT map_0idx  = args[0].map_data[n * args[0].map->dim + CellMap::xu_y_z];
        const OPP_INT map_1idx  = args[1].map_data[n * args[1].map->dim + CellMap::x_yu_z];
        const OPP_INT map_2idx  = args[2].map_data[n * args[2].map->dim + CellMap::x_y_zu];

        half_advance_b_kernel (
            &((const OPP_REAL*) args[0].data)[map_0idx * args[0].dim],     // cell_x_e, 
            &((const OPP_REAL*) args[1].data)[map_1idx * args[1].dim],     // cell_y_e, 
            &((const OPP_REAL*) args[2].data)[map_2idx * args[2].dim],     // cell_z_e, 
            &((const OPP_REAL*) args[3].data)[n * args[3].dim],            // cell0_e, 
            &((OPP_REAL*)       args[4].data)[n * args[4].dim],            // cell0_b
            &((const OPP_INT*)  args[5].data)[n * args[5].dim]             // cell0_ghost
        );
    }
    opp_profiler->end("HalfAdv_B_Kernel");

    // if (set->size == 0 || set->size == set->core_size) {
    //     opp_mpi_halo_wait_all(nargs, args);
    // }

    opp_set_dirtybit(nargs, args);

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

    const int nargs = 7;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;
    args[5] = arg5;
    args[6] = arg6;

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);
    
    opp_profiler->start("Adv_E_Kernel");
    for (int n = 0; n < set_size; n++)
    {
        // if (n == set->core_size) {
        //     opp_mpi_halo_wait_all(nargs, args);
        // }
                
        const OPP_INT map_0idx  = args[0].map_data[n * args[0].map->dim + CellMap::xd_y_z];
        const OPP_INT map_1idx  = args[1].map_data[n * args[1].map->dim + CellMap::x_yd_z];
        const OPP_INT map_2idx  = args[2].map_data[n * args[2].map->dim + CellMap::x_y_zd];

        advance_e_kernel (
            &((const OPP_REAL*) args[0].data)[map_0idx * args[0].dim],     // cell_x_b  
            &((const OPP_REAL*) args[1].data)[map_1idx * args[1].dim],     // cell_y_b  
            &((const OPP_REAL*) args[2].data)[map_2idx * args[2].dim],     // cell_z_b  
            &((const OPP_REAL*) args[3].data)[n * args[3].dim],            // cell0_b   
            &((const OPP_REAL*) args[4].data)[n * args[4].dim],            // cell0_j   
            &((OPP_REAL*)       args[5].data)[n * args[5].dim],            // cell0_e
            &((const OPP_INT*)  args[6].data)[n * args[6].dim]             // iter_adv_e 
        );
    }
    opp_profiler->end("Adv_E_Kernel");

    // if (set->size == 0 || set->size == set->core_size) {
    //     opp_mpi_halo_wait_all(nargs, args);
    // }

    opp_set_dirtybit(nargs, args);

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

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);

    opp_profiler->start("GetMax_Kernel");
    for (int n = 0; n < set_size; n++)
    {
        // if (n == set->core_size) {
        //     opp_mpi_halo_wait_all(nargs, args);
        // }

        get_max_x_values_kernel(
            &((const OPP_REAL*) args[0].data)[n * args[0].dim],     // cell_j  
            (OPP_REAL*)         args[1].data,
            &((const OPP_REAL*) args[2].data)[n * args[2].dim],     // cell_e  
            (OPP_REAL*)         args[3].data,
            &((const OPP_REAL*) args[4].data)[n * args[4].dim],     // cell_b  
            (OPP_REAL*)         args[5].data
        );
    }
    opp_profiler->end("GetMax_Kernel");

    // if (set->size == 0 || set->size == set->core_size) {
    //     opp_mpi_halo_wait_all(nargs, args);
    // }

    opp_mpi_reduce_double(&args[1], (double*)args[1].data);
    opp_mpi_reduce_double(&args[3], (double*)args[3].data);
    opp_mpi_reduce_double(&args[5], (double*)args[5].data);

    opp_set_dirtybit(nargs, args);

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

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);

    opp_profiler->start("UpGhostB_Kernel");
    for (int n = 0; n < set_size; n++)
    {
        // if (n == set->core_size) {
        //     opp_mpi_halo_wait_all(nargs, args);
        // }

        const OPP_INT map_0idx  = args[2].map_data[n * args[2].map->dim + 0];

        update_ghosts_B_kernel(
            &((const OPP_INT*)  args[0].data)[n * args[0].dim],        
            &((const OPP_REAL*) args[1].data)[n * args[1].dim],       
            &((OPP_REAL*)       args[2].data)[map_0idx * args[2].dim],
            (const OPP_INT*)    args[3].data
        );
    }
    opp_profiler->end("UpGhostB_Kernel");

    // if (set->size == 0 || set->size == set->core_size) {
    //     opp_mpi_halo_wait_all(nargs, args);
    // }

    opp_set_dirtybit(nargs, args);

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

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);
    
    opp_profiler->start("UpGhost_Kernel");
    for (int n = 0; n < set_size; n++)
    {
        // if (n == set->core_size) {
        //     opp_mpi_halo_wait_all(nargs, args);
        // }

        const int map_0idx  = args[2].map_data[n * args[2].map->dim + 0];

        update_ghosts_kernel(
            &((const OPP_INT*)  args[0].data)[n * args[0].dim],        
            &((const OPP_REAL*) args[1].data)[n * args[1].dim],       
            &((OPP_REAL*)       args[2].data)[map_0idx * args[2].dim],
            (const OPP_INT*)    args[3].data,
            (const OPP_INT*)    args[4].data
        );
    }
    opp_profiler->end("UpGhost_Kernel");

    // if (set->size == 0 || set->size == set->core_size) {
    //     opp_mpi_halo_wait_all(nargs, args);
    // }

    opp_set_dirtybit(nargs, args);

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

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);
    
    opp_profiler->start("Energy_Kernel");
    for (int n = 0; n < set_size; n++)
    {
        compute_energy_kernel(
            &((const OPP_INT*)  args[0].data)[n * args[0].dim],        
            &((const OPP_REAL*) args[1].data)[n * args[1].dim],       
            (OPP_REAL*)         args[2].data
        );
    }
    opp_profiler->end("Energy_Kernel");

    // if (set->size == 0 || set->size == set->core_size) {
    //     opp_mpi_halo_wait_all(nargs, args);
    // }

    opp_mpi_reduce_double(&args[2], (OPP_REAL*)args[2].data);

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("Energy");   
}

//*************************************************************************************************
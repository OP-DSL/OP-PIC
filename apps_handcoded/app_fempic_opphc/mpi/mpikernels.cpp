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
#include "../fempic_defs.h"

using namespace opp;

//****************************************
OPP_REAL CONST_spwt = 0, CONST_ion_velocity = 0, CONST_dt = 0, CONST_plasma_den = 0, CONST_mass = 0, CONST_charge = 0, CONST_wall_potential = 0;

void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    if (!strcmp(name,"CONST_spwt"))                CONST_spwt = *((OPP_REAL*)data);
    else if (!strcmp(name,"CONST_ion_velocity"))   CONST_ion_velocity = *((OPP_REAL*)data);
    else if (!strcmp(name,"CONST_dt"))             CONST_dt = *((OPP_REAL*)data);
    else if (!strcmp(name,"CONST_plasma_den"))     CONST_plasma_den = *((OPP_REAL*)data);
    else if (!strcmp(name,"CONST_mass"))           CONST_mass = *((OPP_REAL*)data);
    else if (!strcmp(name,"CONST_charge"))         CONST_charge = *((OPP_REAL*)data);
    else if (!strcmp(name,"CONST_wall_potential")) CONST_wall_potential = *((OPP_REAL*)data);
    else std::cerr << "error: unknown const name" << std::endl;
}
//****************************************

#include "../kernels.h"

//*************************************************************************************************
void opp_loop_all__init_boundary_pot(
    opp_set set,      // nodes_set
    opp_arg arg0,     // node_type,
    opp_arg arg1      // node_bnd_pot,
)
{ 
    if (FP_DEBUG) 
        opp_printf("FEMPIC", "opp_loop_all__init_boundary_pot set_size %d diff %d", 
            set->size, set->diff);
    
    opp_profiler->start("InitBndPot");

    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);  

    if (set_size > 0) 
    {
        for (int i = 0; i < set_size; i++)
        {    
            init_boundary_potential(
                &((const OPP_INT*) args[0].data)[i * args[0].dim],     // node_type,
                &((OPP_REAL*)      args[1].data)[i * args[1].dim]      // node_bnd_pot,
            );
        }
    }

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("InitBndPot");
}

//*************************************************************************************************
void opp_loop_inject__inject_ions(
    opp_set set,      // particles_set
    opp_arg arg0,     // part_position,
    opp_arg arg1,     // part_velocity,
    opp_arg arg2,     // part_cell_connectivity,
    opp_arg arg3,     // part_id
    opp_arg arg4,     // iface to cell map
    opp_arg arg5,     // cell_ef,
    opp_arg arg6,     // iface_u,
    opp_arg arg7,     // iface_v,
    opp_arg arg8,     // iface_normal,
    opp_arg arg9,     // iface_node_pos
    opp_arg arg10     // dummy_part_random
)
{ 
    if (FP_DEBUG) 
        opp_printf("FEMPIC", "opp_loop_inject__inject_ions set_size %d diff %d", 
            set->size, set->diff);
    
    opp_profiler->start("InjectIons");

    const int inj_start = (set->size - set->diff);

    const int nargs = 11;
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

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args); // 

    if (set_size > 0) 
    {
        opp_profiler->start("InjectIons_Kernel");
        OPP_INT map0idx = -1, map1idx = 0;
        const int iter_size = set->diff;   
        for (int i = 0; i < iter_size; i++)
        {    
            map0idx = ((OPP_INT *)set->mesh_relation_dat->data)[inj_start + i]; // iface index
            map1idx = args[5].map_data[map0idx]; // cell index

            inject_ions__kernel(
                &((OPP_REAL*)       args[0].data)[(inj_start + i) * args[0].dim],     // part_position,
                &((OPP_REAL*)       args[1].data)[(inj_start + i) * args[1].dim],     // part_velocity,
                &((OPP_INT*)        args[2].data)[(inj_start + i) * args[2].dim],     // part_cell_connectivity,
                &((OPP_INT*)        args[3].data)[(inj_start + i) * args[3].dim],     // part_id
                &((const OPP_INT*)  args[4].data)[map0idx * args[4].dim],             // iface to cell map
                &((const OPP_REAL*) args[5].data)[map1idx * args[5].dim],             // cell_ef,
                &((const OPP_REAL*) args[6].data)[map0idx * args[6].dim],             // iface_u,
                &((const OPP_REAL*) args[7].data)[map0idx * args[7].dim],             // iface_v,
                &((const OPP_REAL*) args[8].data)[map0idx * args[8].dim],             // iface_normal,
                &((const OPP_REAL*) args[9].data)[map0idx * args[9].dim],             // iface_node_pos
                &((const OPP_REAL*) args[10].data)[i * args[10].dim]                  // dummy_part_random
            );
        }
        opp_profiler->end("InjectIons_Kernel");
    }

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("InjectIons");
}

//*************************************************************************************************
void opp_loop_all__calc_pos_vel(
    opp_set set,      // particles_set
    opp_arg arg0,     // cell_ef,
    opp_arg arg1,     // part_pos,
    opp_arg arg2      // part_vel,
) 
{
    if (FP_DEBUG) 
        opp_printf("FEMPIC", "opp_loop_all__calc_pos_vel set_size %d diff %d", 
            set->size, set->diff);

    opp_profiler->start("CalcPosVel");  

    int nargs = 3;
    opp_arg args[nargs];

    args[0]  = std::move(arg0);
    args[1]  = std::move(arg1);
    args[2]  = std::move(arg2);

    opp_mpi_halo_exchanges(set, nargs, args);

    opp_mpi_halo_wait_all(nargs, args); 

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 

    opp_profiler->start("CalcPosVel_Kernel");
    const int iter_size = set->size;
    for (int i = 0; i < iter_size; i++)
    { 
        const OPP_INT map0idx = OPP_mesh_relation_data[i];

        calculate_new_pos_vel__kernel(
            &((const OPP_REAL*) arg0.data)[map0idx * arg0.dim],  // cell_ef,
            &((OPP_REAL*)       arg1.data)[i * arg1.dim],        // part_pos,
            &((OPP_REAL*)       arg2.data)[i * arg2.dim]         // part_vel,
        );
    }
    opp_profiler->end("CalcPosVel_Kernel");

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("CalcPosVel");    
}

//*************************************************************************************************
void opp_loop_all__deposit_charge_on_nodes(
    opp_set set,      // particles_set
    opp_arg arg0,     // part_lc,
    opp_arg arg1,     // node_charge_den0,
    opp_arg arg2,     // node_charge_den1,
    opp_arg arg3,     // node_charge_den2,
    opp_arg arg4      // node_charge_den3,
)
{
    if (FP_DEBUG) 
        opp_printf("FEMPIC", "opp_loop_all__deposit_charge_on_nodes set_size %d diff %d", 
            set->size, set->diff);

    opp_profiler->start("DepCharge");  

    const int nargs = 5;
    opp_arg args[nargs];

    args[0]  = arg0;
    args[1]  = arg1;
    args[2]  = arg2;
    args[3]  = arg3;
    args[4]  = arg4;

    opp_mpi_halo_exchanges(set, nargs, args);

    opp_init_double_indirect_reductions(nargs, args);

    // unable to overlap much of computation and communication
    opp_mpi_halo_wait_all(nargs, args); 

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 

    opp_profiler->start("DepCharge_Kernel");
    const int iter_size = set->size;
    for (int i = 0; i < iter_size; i++)
    { 
        const OPP_INT map0idx = OPP_mesh_relation_data[i];

        const OPP_INT map1idx = arg1.map_data[map0idx * arg1.map->dim + 0];
        const OPP_INT map2idx = arg1.map_data[map0idx * arg1.map->dim + 1];
        const OPP_INT map3idx = arg1.map_data[map0idx * arg1.map->dim + 2];
        const OPP_INT map4idx = arg1.map_data[map0idx * arg1.map->dim + 3];

        deposit_charge_on_nodes__kernel(
            &((const OPP_REAL*) arg0.data)[i * arg0.dim],      // part_lc,
            &((OPP_REAL*)       arg1.data)[map1idx],           // node_charge_den0,
            &((OPP_REAL*)       arg1.data)[map2idx],           // node_charge_den1,
            &((OPP_REAL*)       arg1.data)[map3idx],           // node_charge_den2,
            &((OPP_REAL*)       arg1.data)[map4idx]            // node_charge_den3,
        );
    }
    opp_profiler->end("DepCharge_Kernel");

    opp_exchange_double_indirect_reductions(nargs, args);

    opp_complete_double_indirect_reductions(nargs, args);

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("DepCharge");   
}

// void opp_particle_move__move_fused(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,
//     opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
// // #define DEBUG_INTERNAL
// //*************************************************************************************************
// void opp_particle_move__move_fused(
//     opp_set set,      // particles_set
//     opp_arg arg0,     // cell_ef,
//     opp_arg arg1,     // part_pos,
//     opp_arg arg2,     // part_vel,
//     opp_arg arg3,     // part_lc,
//     opp_arg arg4,     // current_cell_index,
//     opp_arg arg5,     // current_cell_volume,
//     opp_arg arg6,     // current_cell_det,
//     opp_arg arg7,     // cell_connectivity,
//     opp_arg arg8,     // node_charge_den0,
//     opp_arg arg9,     // node_charge_den1,
//     opp_arg arg10,    // node_charge_den2,
//     opp_arg arg11     // node_charge_den3,
// )
// { 

//     if (FP_DEBUG) 
//         opp_printf("FEMPIC", "opp_particle_move__move_fused set_size %d diff %d", 
//             set->size, set->diff);

//     opp_profiler->start("MoveToCells");

// #ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
//     double kernel_time = 0.0;
//     int total_particles = 0;
//     int comm_iteration = 0;
//     int max_internal_hops = 0, internal_hops = 0;
//     std::vector<int> particle_loops_per_comm_iter(10, 0);
//     auto total_start = std::chrono::system_clock::now();
//     auto start = std::chrono::system_clock::now();
// #endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

//     int nargs = 12;
//     opp_arg args[nargs];

//     args[0]  = std::move(arg0);
//     args[1]  = std::move(arg1);
//     args[2]  = std::move(arg2);
//     args[3]  = std::move(arg3);
//     args[4]  = std::move(arg4);
//     args[5]  = std::move(arg5);
//     args[6]  = std::move(arg6);
//     args[7]  = std::move(arg7);
//     args[8]  = std::move(arg8);
//     args[9]  = std::move(arg9);
//     args[10] = std::move(arg10);
//     args[11] = std::move(arg11);

//     OPP_INT *map0idx = nullptr;

//     opp_mpi_halo_exchanges(set, nargs, args);

//     // if there is access to a dat with OPP_Map_from_Mesh_Rel and a mapping, 
//     // then we should reduce the contributions to the element containing rank
//     // Here we should make the values of that dat to zero prior loop, 
//     // execute the loop and communicate the outcome to the residing rank, like in a halo exchange, 
//     // but when received, that rank should do the reduction
//     opp_init_double_indirect_reductions(nargs, args);

//     // unable to overlap much of computation and communication
//     opp_mpi_halo_wait_all(nargs, args); 

//     do // iterate until all mpi ranks say, I am done
//     {
//         opp_init_particle_move(set, nargs, args);
        
//         if (FP_DEBUG) 
//             opp_printf("FEMPIC", "opp_loop_particle_all__MoveToCells Starting iteration %d, start[%d] end[%d]", 
//                 OPP_comm_iteration, OPP_iter_start, OPP_iter_end);

// #ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
//         start = std::chrono::system_clock::now();
// #endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

//         for (int i = OPP_iter_start; i < OPP_iter_end; i++)
//         {
// #ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
//             internal_hops = 0;
// #endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

//             OPP_MOVE_RESET_FLAGS;
            
//             do
//             { 
//                 map0idx = &(OPP_mesh_relation_data[i]);

//                 const OPP_INT map1idx = args[8].map_data[*map0idx * args[8].map->dim + 0];
//                 const OPP_INT map2idx = args[8].map_data[*map0idx * args[8].map->dim + 1];
//                 const OPP_INT map3idx = args[8].map_data[*map0idx * args[8].map->dim + 2];
//                 const OPP_INT map4idx = args[8].map_data[*map0idx * args[8].map->dim + 3];

//                 move_all_particles_to_cell__kernel(
//                     &((const OPP_REAL*) args[0].data)[*map0idx * args[0].dim],  // cell_ef,
//                     &((OPP_REAL*)       args[1].data)[i * args[1].dim],         // part_pos,
//                     &((OPP_REAL*)       args[2].data)[i * args[2].dim],         // part_vel,
//                     &((OPP_REAL*)       args[3].data)[i * args[3].dim],         // part_lc,
//                     &((OPP_INT *)       args[4].data)[i * args[4].dim],         // current_cell_index,
//                     &((const OPP_REAL*) args[5].data)[*map0idx * args[5].dim],  // current_cell_volume,
//                     &((const OPP_REAL*) args[6].data)[*map0idx * args[6].dim],  // current_cell_det,
//                     &((const OPP_INT*)  args[7].data)[*map0idx * args[7].dim],  // cell_connectivity,
//                     &((OPP_REAL*)       args[8].data)[map1idx * args[8].dim],   // node_charge_den0,
//                     &((OPP_REAL*)       args[8].data)[map2idx * args[8].dim],   // node_charge_den1,
//                     &((OPP_REAL*)       args[8].data)[map3idx * args[8].dim],   // node_charge_den2,
//                     &((OPP_REAL*)       args[8].data)[map4idx * args[8].dim]    // node_charge_den3,
//                 );

// #ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
//                 internal_hops++; // can remove
// #endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

//             // opp_part_check_status checks whether map0idx is in halo list, 
//             // if yes, pack the particle into MPI buffer and set status to NEED_REMOVE
//             } while (opp_check_part_move_status(*map0idx, i, set->particle_remove_count));

// #ifdef DEBUG_INTERNAL
//             if (max_internal_hops < internal_hops) 
//                 max_internal_hops = internal_hops; // can remove
// #endif // DEBUG_INTERNAL ----------------------------------------------------------------------------
//         }

// #ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
//         std::chrono::duration<double> diff   = std::chrono::system_clock::now() - start;
//         kernel_time += (double)diff.count();
//         total_particles += (OPP_iter_end - OPP_iter_start);
//         particle_loops_per_comm_iter[OPP_comm_iteration] = (OPP_iter_end - OPP_iter_start);
//         comm_iteration++;
// #endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

//     } while (opp_finalize_particle_move(set)); // iterate until all mpi ranks say, I am done

//     opp_exchange_double_indirect_reductions(nargs, args);

//     // if auto_sort is set, then could sort here

//     // TODO : can this be added to opp_mpi_halo_exchanges, to complete before the next usage?
//     opp_complete_double_indirect_reductions(nargs, args);

//     opp_set_dirtybit(nargs, args);

// #ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
//     std::chrono::duration<double> total_diff   = std::chrono::system_clock::now() - total_start;
//     opp_printf("MoveToCells", "TotalTime: %2.15lE KernelTime: %2.15lE | total_particles: %d | \
//         particle_loops_per_comm_iter [%d %d %d %d] | comm_iteration: %d max_internal_hops: %d", 
//         (double)total_diff.count(), kernel_time, total_particles, 
//         particle_loops_per_comm_iter[0], particle_loops_per_comm_iter[1], particle_loops_per_comm_iter[2], 
//         particle_loops_per_comm_iter[3], comm_iteration, max_internal_hops);
// #endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

//     opp_profiler->end("MoveToCells");
// }

//*************************************************************************************************
void opp_loop_all__compute_node_charge_density(
    opp_set set,     // nodes_set
    opp_arg arg0,    // node_charge_density
    opp_arg arg1     // node_volume
)
{ 

    if (FP_DEBUG) 
        opp_printf("FEMPIC", "opp_loop_all__compute_node_charge_density set_size %d", 
            set->size);

    opp_profiler->start("ComputeNodeChargeDensity");

    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args); 
    
    if (set_size > 0) 
    {
        opp_profiler->start("ComputeNChargeDen_Kernel");
        for (int i=0; i < set_size; i++) 
        {
            // if (i == set->core_size) 
            // {
            //     opp_mpi_halo_wait_all(nargs, args);
            // }

            compute_node_charge_density__kernel(
                &((OPP_REAL*)       args[0].data)[i * args[0].dim],
                &((const OPP_REAL*) args[1].data)[i * args[1].dim]
            );
        }
        opp_profiler->end("ComputeNChargeDen_Kernel");
    }  

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("ComputeNodeChargeDensity");
}

//*************************************************************************************************
void opp_loop_all__compute_electric_field(
    opp_set set,      // cells_set
    opp_arg arg0,     // cell_electric_field,
    opp_arg arg1,     // cell_shape_deriv,
    opp_arg arg2,     // node_potential0,
    opp_arg arg3,     // node_potential1,
    opp_arg arg4,     // node_potential2,
    opp_arg arg5      // node_potential3,
)
{ 

    if (FP_DEBUG) 
        opp_printf("FEMPIC", "opp_loop_all__compute_electric_field set_size %d", 
            set->size);

    opp_profiler->start("ComputeElectricField");

    int nargs = 6;
    opp_arg args[nargs];

    args[0] = std::move(arg0);
    args[1] = std::move(arg1);
    args[2] = std::move(arg2);
    args[3] = std::move(arg3);
    args[4] = std::move(arg4);
    args[5] = std::move(arg5);

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);
    
    if (set_size > 0) 
    {
        opp_profiler->start("ComputeEF_Kernel");
        for (int i = 0; i < set_size; i++)
        {
            // if (i == set->core_size) 
            // {
            //     opp_mpi_halo_wait_all(nargs, args);
            // }

            const OPP_INT map1idx = args[2].map_data[i * args[2].map->dim + 0];
            const OPP_INT map2idx = args[2].map_data[i * args[2].map->dim + 1];
            const OPP_INT map3idx = args[2].map_data[i * args[2].map->dim + 2];
            const OPP_INT map4idx = args[2].map_data[i * args[2].map->dim + 3];

            compute_electric_field__kernel(
                &((OPP_REAL*)args[0].data)[i * args[0].dim],               // cell_electric_field
                &((const OPP_REAL*)args[1].data)[i * args[1].dim],         // cell_shape_deriv
                &((const OPP_REAL*)args[2].data)[map1idx * args[2].dim],   // node_potential0
                &((const OPP_REAL*)args[2].data)[map2idx * args[2].dim],   // node_potential1
                &((const OPP_REAL*)args[2].data)[map3idx * args[2].dim],   // node_potential2
                &((const OPP_REAL*)args[2].data)[map4idx * args[2].dim]    // node_potential3
            );
        } 
        opp_profiler->end("ComputeEF_Kernel");
    }  

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("ComputeElectricField");
}





//*************************************************************************************************
inline void generate_dh_mappings(opp_set cells_set, const opp_dat global_cell_id_dat, 
        const opp_dat cellVolume_dat, const opp_dat cellDet_dat) { 

    // if (OPP_DBG) 
    if (OPP_rank == 0)            
        opp_printf("FUNC", "generate_dh_mappings cells [%s] start global grid dimensions %d %d %d",
            cells_set->name, cellMapper->globalGridDims.x, cellMapper->globalGridDims.y, cellMapper->globalGridDims.z);

    const int cellSetSizeIncHalo = cells_set->size + cells_set->exec_size + cells_set->nonexec_size;
    if (cellSetSizeIncHalo <= 0) {
        opp_printf("FUNC", "Error... cellSetSizeIncHalo <= 0 for set %s, Terminating...",
            cells_set->name);
        opp_abort("Error... FUNC cellSetSizeIncHalo <= 0");
    }

    double x = 0.0, y = 0.0, z = 0.0;
    double lc[4];
    std::map<size_t, opp_point> removedCoordinates;
    const opp_point& minGlbCoordinate = boundingBox->getGlobalMin();
    const opp_point& maxCoordinate = boundingBox->getLocalMax(); // required for GET_VERT

    auto all_cell_checker = [&](const opp_point& point, int& cellIndex) { 
        bool isInside;           
        for (int ci = 0; ci < cells_set->size; ci++) {
            isPointInCellKernel( 
                            isInside,
                            (const OPP_REAL*)&point, 
                            lc,
                            &((OPP_REAL*)cellVolume_dat->data)[ci * cellVolume_dat->dim], 
                            &((OPP_REAL*)cellDet_dat->data)[ci * cellDet_dat->dim]);
            if (isInside) {       
                cellIndex = ci;
                break;
            }
        }
    };

    // auto neighbour_cell_checker = [&](const opp_point& point, int& cellIndex) {
    //     opp_move_status m;
    //     cellIndex = cells_set->size / 2;
    //     do {
    //         m = getCellIndexKernel(
    //             (const double*)&point, 
    //             &cellIndex,
    //             lc,
    //             &((double*)cellVolume_dat->data)[cellIndex], 
    //             &((double*)cellDet_dat->data)[cellIndex * cellDet_dat->dim], 
    //             &((int*)cellConnectivity_map->map)[cellIndex * cellConnectivity_map->dim]);
    //     } while (m == OPP_NEED_MOVE && cellIndex < cellSetSizeIncHalo); 
    // };

    cellMapper->createStructMeshMappingArrays();

    // Step 1 : Get the centroids of the structured mesh cells and try to relate them to unstructured mesh cell indices
    for (int dz = cellMapper->localGridStart.z; dz < cellMapper->localGridEnd.z; dz++) {
        
        z = minGlbCoordinate.z + dz * cellMapper->gridSpacing;
        
        for (int dy = cellMapper->localGridStart.y; dy < cellMapper->localGridEnd.y; dy++) {
            
            y = minGlbCoordinate.y + dy * cellMapper->gridSpacing;
            
            for (int dx = cellMapper->localGridStart.x; dx < cellMapper->localGridEnd.x; dx++) {
                
                x = minGlbCoordinate.x + dx * cellMapper->gridSpacing;
                
                size_t index = (size_t)(dx + dy * cellMapper->globalGridDims.x + 
                    dz * cellMapper->globalGridDims.x * cellMapper->globalGridDims.y);
                
                const opp_point centroid = cellMapper->getCentroidOfBox(opp_point(x, y ,z));

                int cellIndex = MAX_CELL_INDEX;

                all_cell_checker(centroid, cellIndex); // Find in which cell this centroid lies

                if (cellIndex == MAX_CELL_INDEX) {
                    removedCoordinates.insert(std::make_pair(index, opp_point(x, y ,z)));
                    continue;
                }

                if (cellIndex < cells_set->size) { // write only if the structured cell belong to the current MPI rank
                    
                    cellMapper->enrichStructuredMesh(index, ((int*)global_cell_id_dat->data)[cellIndex], OPP_rank);
                } 
            }
        }
    }

#ifdef USE_MPI
    // Step 2 : For MPI, get the inter-node values reduced to the structured mesh

    cellMapper->reduceInterNodeMappings(1);

    // The marked structured cells from this rank might be filled by another rank, so if already filled, no need to recalculate from current rank
    for (auto it = removedCoordinates.begin(); it != removedCoordinates.end(); ) {

        size_t removedIndex = it->first;
        
        if (cellMapper->structMeshToRankMapping[removedIndex] != MAX_CELL_INDEX) {

            it = removedCoordinates.erase(it); // This structured index is already written by another rank
            // opp_printf("FUNC", "Already written %d to struct index %zu", this->structMeshToRankMapping[removedIndex], removedIndex);
        } else {
            ++it; 
        }
    }

    cellMapper->waitBarrier();
#endif

    // Step 3 : Iterate over all the NEED_REMOVE points, try to check whether atleast one vertex of the structured mesh can be within 
    //          an unstructured mesh cell. If multiple are found, get the minimum cell index to match with MPI
    for (auto& p : removedCoordinates) {

        size_t index = p.first;
        double& x = p.second.x;
        double& y = p.second.y;
        double& z = p.second.z;
            
        // still no one has claimed that this cell belongs to it

        const double gs = cellMapper->gridSpacing;
        int mostSuitableCellIndex = MAX_CELL_INDEX, mostSuitableGblCellIndex = MAX_CELL_INDEX;

        std::array<opp_point,8> vertices = {
            opp_point(GET_VERT(x,x),    GET_VERT(y,y),    GET_VERT(z,z)),
            opp_point(GET_VERT(x,x),    GET_VERT(y,y+gs), GET_VERT(z,z)),
            opp_point(GET_VERT(x,x),    GET_VERT(y,y+gs), GET_VERT(z,z+gs)),
            opp_point(GET_VERT(x,x),    GET_VERT(y,y),    GET_VERT(z,z+gs)),
            opp_point(GET_VERT(x,x+gs), GET_VERT(y,y),    GET_VERT(z,z)),
            opp_point(GET_VERT(x,x+gs), GET_VERT(y,y+gs), GET_VERT(z,z)),
            opp_point(GET_VERT(x,x+gs), GET_VERT(y,y),    GET_VERT(z,z+gs)),
            opp_point(GET_VERT(x,x+gs), GET_VERT(y,y+gs), GET_VERT(z,z+gs)),
        };

        for (const auto& point : vertices) {

            int cellIndex = MAX_CELL_INDEX;

            all_cell_checker(point, cellIndex);

            if (cellIndex != MAX_CELL_INDEX && (cellIndex < cellSetSizeIncHalo)) { 

                int gblCellIndex = ((int*)global_cell_id_dat->data)[cellIndex];

                if (mostSuitableGblCellIndex > gblCellIndex) {
                    mostSuitableGblCellIndex = gblCellIndex;
                    mostSuitableCellIndex = cellIndex;
                }
            }
        }    

        cellMapper->lockWindows();

        int alreadyAvailGblCellIndex = cellMapper->structMeshToCellMapping[index];

        // Allow neighbours to write on-behalf of the current rank, to reduce issues
        if (mostSuitableGblCellIndex != MAX_CELL_INDEX && mostSuitableGblCellIndex < alreadyAvailGblCellIndex && 
            (mostSuitableCellIndex < cells_set->size)) {
            
            cellMapper->enrichStructuredMesh(index, mostSuitableGblCellIndex, OPP_rank);       
        }

        cellMapper->unlockWindows();
    }

#ifdef USE_MPI
    // Step 4 : For MPI, get the inter-node values reduced to the structured mesh
    cellMapper->reduceInterNodeMappings(2);
    
    // Step 5 : For MPI, convert the global cell coordinates to rank local coordinates for increased performance
    cellMapper->convertToLocalMappings(global_cell_id_dat);
#endif

    // if (OPP_DBG) 
    if (OPP_rank == 0)
        opp_printf("FUNC", "generate_dh_mappings end");
}




//*******************************************************************************
void init_particle_mover(const double gridSpacing, const int dim, const opp_dat node_pos_dat, 
    const opp_dat cellVolume_dat, const opp_dat cellDet_dat, const opp_dat global_cell_id_dat) {
    
    opp_profiler->start("SetupMover");

    useGlobalMove = opp_params->get<OPP_BOOL>("opp_global_move");

    if (useGlobalMove) {
        
#ifdef USE_MPI
        comm = std::make_shared<Comm>(MPI_COMM_WORLD);
        globalMover = std::make_unique<GlobalParticleMover>(comm->comm_parent);
#endif
        boundingBox = std::make_shared<BoundingBox>(node_pos_dat, dim);

        cellMapper = std::make_shared<CellMapper>(boundingBox, gridSpacing, comm);

        generate_dh_mappings(cellVolume_dat->set, global_cell_id_dat, 
            cellVolume_dat, cellDet_dat);
    }

    opp_profiler->reg("GlbToLocal");
    opp_profiler->reg("GblMv_Move");
    opp_profiler->reg("GblMv_AllMv");
    for (int i = 0; i < 5; i++) {
        std::string profName = std::string("Mv_AllMv") + std::to_string(i);
        opp_profiler->reg(profName);
    }

    opp_profiler->end("SetupMover");
}


//*******************************************************************************
void opp_particle_move__move(
    opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // part_position,   
    opp_arg arg1, // part_lc,       
    opp_arg arg2, // cell_volume, 
    opp_arg arg3  // cell_det,   
) 
{
    if (FP_DEBUG) opp_printf("FEMPIC", "move set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("Move");

    const int nargs = 5;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

    OPP_mesh_relation_data = (OPP_INT*)p2c_map->p2c_dat->data;

    // lambda function for multi hop particle movement
    auto multihop_mover = [&](const int i) {

        opp_p2c = &((OPP_mesh_relation_data)[i]);   // TODO : remove OPP_INT* after making this into a map

        if (*opp_p2c == MAX_CELL_INDEX) {
            return;
        }

        OPP_MOVE_RESET_FLAGS;

        do {
            opp_c2c = &((c2c_map->map)[*opp_p2c * 4]);

            getCellIndexKernel(
                &((const OPP_REAL*) args[0].data)[i * 3], 
                &((OPP_REAL*)       args[1].data)[i * 4],
                &((const OPP_REAL*) args[2].data)[*opp_p2c * 1], 
                &((const OPP_REAL*) args[3].data)[*opp_p2c * 16]);

        } while (opp_check_part_move_status(*opp_p2c, i, set->particle_remove_count));
    };

    // ----------------------------------------------------------------------------
    opp_init_particle_move(set, 0, nullptr);
    
    if (useGlobalMove) {
        
        globalMover->initGlobalMove();

        opp_profiler->start("GblMv_Move");

        // check whether particles needs to be moved over global move routine
        const int start = OPP_iter_start;
        const int end = OPP_iter_end;
        for (int i = start; i < end; i++) {   
            
            opp_p2c = &(OPP_mesh_relation_data)[i];   // TODO : remove OPP_INT* after making this into a map
            const opp_point* point = (const opp_point*)&(((OPP_REAL*)args[0].data)[i * args[0].dim]);

            // check for global move, and if satisfy global move criteria, then remove the particle from current rank
            if (opp_part_checkForGlobalMove(set, *point, i, *opp_p2c)) {
                
                set->particle_remove_count++;
                continue;  
            }
        }

        opp_profiler->end("GblMv_Move");

        globalMover->communicate(set);
    }

    opp_profiler->start("Mv_AllMv0");

    // ----------------------------------------------------------------------------
    // check whether all particles not marked for global comm is within cell, 
    // and if not mark to move between cells within the MPI rank, mark for neighbour comm
    opp_profiler->start("Move_Kernel");
    const int start1 = OPP_iter_start;
    const int end1 = OPP_iter_end;
    for (int i = start1; i < end1; i++) { 
        
        multihop_mover(i);
    }
    opp_profiler->end("Move_Kernel");

    opp_profiler->end("Mv_AllMv0");

    // ----------------------------------------------------------------------------
    // finalize the global move routine and iterate over newly added particles and check whether they need neighbour comm
    if (useGlobalMove && globalMover->finalize(set) > 0) {
        
        opp_profiler->start("GblMv_AllMv");

        // need to change arg data since particle resize in globalMover::finalize could change the pointer in dat->data 
        for (int i = 0; i < nargs; i++)
            if (args[i].argtype == OPP_ARG_DAT && args[i].dat->set->is_particle)
                args[i].data = args[i].dat->data;

        // check whether the new particle is within cell, and if not move between cells within the MPI rank, 
        // mark for neighbour comm. Do only for the globally moved particles 
        opp_profiler->start("Move_Kernel");
        const int start2 = (set->size - set->diff);
        const int end2 = set->size;
        for (int i = start2; i < end2; i++) { 
                
            multihop_mover(i);                 
        }
        opp_profiler->end("Move_Kernel");

        opp_profiler->end("GblMv_AllMv");
    }

    // ----------------------------------------------------------------------------
    // Do neighbour communication and if atleast one particle is received by the currect rank, 
    // then iterate over the newly added particles
    while (opp_finalize_particle_move(set)) {
        
        const std::string profName = std::string("Mv_AllMv") + std::to_string(OPP_comm_iteration);
        opp_profiler->start(profName);
        
        opp_init_particle_move(set, 0, nullptr);

        // check whether particle is within cell, and if not move between cells within the MPI rank, mark for neighbour comm
        opp_profiler->start("Move_Kernel");
        const int start3 = OPP_iter_start;
        const int end3 = OPP_iter_end;
        for (int i = start3; i < end3; i++) {      

            multihop_mover(i);
        }
        opp_profiler->end("Move_Kernel");

        opp_profiler->end(profName);
    }

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("Move");
}

//*************************************************************************************************
void opp_loop_all__get_max_values(
    opp_set set,     // cells set
    opp_arg arg0,    // n_charge_den            // OPP_READ
    opp_arg arg1,    // global_max_n_charge_den // OPP_MAX
    opp_arg arg2,    // n_pot                   // OPP_READ
    opp_arg arg3     // global_max_n_pot        // OPP_MAX
)
{
    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__get_max_values set_size %d", set->size);

    opp_profiler->start("GetMax");

    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);  

    if (set_size > 0)
    {
        opp_profiler->start("GetMax__Kernel");
        for (int n = 0; n < set_size; n++)
        {
            get_final_max_values_kernel(
                &((const OPP_REAL*) args[0].data)[n * args[0].dim],     // n_charge_den  
                (OPP_REAL*)         args[1].data,
                &((const OPP_REAL*) args[2].data)[n * args[2].dim],     // n_pot  
                (OPP_REAL*)         args[3].data
            );
        }
        opp_profiler->end("GetMax__Kernel");
    }

    opp_mpi_reduce_double(&args[1], (OPP_REAL*)args[1].data);
    opp_mpi_reduce_double(&args[3], (OPP_REAL*)args[3].data);

    opp_set_dirtybit(nargs, args);
    
    opp_profiler->end("GetMax");
}

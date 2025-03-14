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


#include "opp_seq.h"
#include "../fempic_defs.h"

using namespace opp;

//****************************************
double CONST_spwt = 0, CONST_ion_velocity = 0, CONST_dt = 0, CONST_plasma_den = 0, CONST_mass = 0, CONST_charge = 0, CONST_wall_potential = 0;
void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    if (!strcmp(name,"CONST_spwt"))              CONST_spwt = *((double*)data);
    else if (!strcmp(name,"CONST_ion_velocity")) CONST_ion_velocity = *((double*)data);
    else if (!strcmp(name,"CONST_dt"))           CONST_dt = *((double*)data);
    else if (!strcmp(name,"CONST_plasma_den"))   CONST_plasma_den = *((double*)data);
    else if (!strcmp(name,"CONST_mass"))         CONST_mass = *((double*)data);
    else if (!strcmp(name,"CONST_charge"))       CONST_charge = *((double*)data);
    else if (!strcmp(name,"CONST_wall_potential")) CONST_wall_potential = *((double*)data);
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

    const int set_size = set->size;

    for (int i = 0; i < set_size; i++)
    {    
        init_boundary_potential(
            &((OPP_INT*) arg0.data)[i * arg0.dim],     // node_type,
            &((OPP_REAL*) arg1.data)[i * arg1.dim]   // node_bnd_pot,
        );
    }

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

    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_inject__inject_ions set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("InjectIons");

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

    const int inj_start = (set->size - set->diff);
    OPP_INT map0idx = -1, map1idx = 0;
    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 
    const int iter_size = set->diff;

    for (int i = 0; i < iter_size; i++)
    {    
        map0idx = OPP_mesh_relation_data[inj_start + i]; // iface index
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

    opp_profiler->end("InjectIons");
}

void opp_loop_all__calc_pos_vel(
    opp_set set,      // particles_set
    opp_arg arg0,     // cell_ef,
    opp_arg arg1,     // part_pos,
    opp_arg arg2      // part_vel,
) {

    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__calc_pos_vel set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("CalcPosVel");  

    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 
    const int iter_size = set->size;

    for (int i = 0; i < iter_size; i++)
    { 
        const int map0idx = OPP_mesh_relation_data[i];

        calculate_new_pos_vel__kernel(
            &((const OPP_REAL*) arg0.data)[map0idx * arg0.dim],  // cell_ef,
            &((OPP_REAL*)       arg1.data)[i * arg1.dim],        // part_pos,
            &((OPP_REAL*)       arg2.data)[i * arg2.dim]         // part_vel,
        );
    }

    opp_profiler->end("CalcPosVel");    
}

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

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 
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

    opp_profiler->end("DepCharge");   
}

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

//     if (FP_DEBUG) opp_printf("FEMPIC", "opp_particle_move__move_fused set_size %d diff %d", 
//         set->size, set->diff);

//     opp_profiler->start("MoveToCells");

//     int *map0idx = nullptr;
// int over_one = 0;
//     opp_init_particle_move(set, 0, nullptr);
// int map1idx, map2idx, map3idx, map4idx;

//     for (int i = OPP_iter_start; i < OPP_iter_end; i++)
//     {        
//         opp_move_var m = opp_get_move_var();
// int in_count = 0;
//         do
//         { 
//             map0idx = &(OPP_mesh_relation_data[i]);

//             map1idx = arg8.map_data[*map0idx * arg8.map->dim + 0];
//             map2idx = arg8.map_data[*map0idx * arg8.map->dim + 1];
//             map3idx = arg8.map_data[*map0idx * arg8.map->dim + 2];
//             map4idx = arg8.map_data[*map0idx * arg8.map->dim + 3];

//             move_all_particles_to_cell__kernel(
//                 (m),
//                 &((double*) arg0.data)[*map0idx * arg0.dim],  // cell_ef,
//                 &((double*) arg1.data)[i * arg1.dim],         // part_pos,
//                 &((double*) arg2.data)[i * arg2.dim],         // part_vel,
//                 &((double*) arg3.data)[i * arg3.dim],         // part_lc,
//                 &((int *)   arg4.data)[i * arg4.dim],         // current_cell_index,
//                 &((double*) arg5.data)[*map0idx * arg5.dim],  // current_cell_volume,
//                 &((double*) arg6.data)[*map0idx * arg6.dim],  // current_cell_det,
//                 &((int*)    arg7.data)[*map0idx * arg7.dim],  // cell_connectivity,
//                 &((double*) arg8.data)[map1idx],                 // node_charge_den0,
//                 &((double*) arg8.data)[map2idx],                 // node_charge_den1,
//                 &((double*) arg8.data)[map3idx],                 // node_charge_den2,
//                 &((double*) arg8.data)[map4idx]                  // node_charge_den3,
//             );         
// in_count++;
//         } while (opp_part_check_status(m, *map0idx, set, i, set->particle_remove_count));  
// if (in_count > 1) {
//     over_one++;
// }
//     }
// // opp_printf("FEMPIC", "opp_particle_move__move_fused Iterated=%d, parts with inner over one iter=%d", 
// //         (OPP_iter_end - OPP_iter_start), over_one);
//     opp_finalize_particle_move(set);

//     opp_profiler->end("MoveToCells");
// }

//*************************************************************************************************
void opp_loop_all__compute_node_charge_density(
    opp_set set,     // nodes_set
    opp_arg arg0,    // node_charge_density
    opp_arg arg1     // node_volume
)
{
    
    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__compute_node_charge_density set_size %d", set->size);

    opp_profiler->start("ComputeNodeChargeDensity");

    const int iter_size = set->size;

    for (int i = 0; i < iter_size; i++) 
    {
        compute_node_charge_density__kernel(
            &((OPP_REAL*)       arg0.data)[i * arg0.dim],
            &((const OPP_REAL*) arg1.data)[i * arg1.dim]
        );
    }

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

    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__compute_electric_field set_size %d", set->size);

    opp_profiler->start("ComputeElectricField");

    const int iter_size = set->size;

    for (int i = 0; i < iter_size; i++)
    {
        const OPP_INT map1idx = arg2.map_data[i * arg2.map->dim + 0];
        const OPP_INT map2idx = arg2.map_data[i * arg2.map->dim + 1];
        const OPP_INT map3idx = arg2.map_data[i * arg2.map->dim + 2];
        const OPP_INT map4idx = arg2.map_data[i * arg2.map->dim + 3];

        compute_electric_field__kernel(
            &((OPP_REAL*)arg0.data)[i * arg0.dim],          // cell_electric_field
            &((const OPP_REAL*)arg1.data)[i * arg1.dim],    // cell_shape_deriv
            &((const OPP_REAL*)arg2.data)[map1idx],         // node_potential0
            &((const OPP_REAL*)arg2.data)[map2idx],         // node_potential1
            &((const OPP_REAL*)arg2.data)[map3idx],         // node_potential2
            &((const OPP_REAL*)arg2.data)[map4idx]          // node_potential3
        );
    }

    opp_profiler->end("ComputeElectricField");   
}




//*************************************************************************************************
inline void generateStructMeshToGlobalCellMappings(opp_set cells_set, const opp_dat global_cell_id_dat, 
        const opp_dat cellVolume_dat, const opp_dat cellDet_dat) { 

    // if (OPP_DBG) 
    if (OPP_rank == 0)            
        opp_printf("FUNC", "generateStructMeshToGlobalCellMappings cells [%s] start global grid dimensions %d %d %d",
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

    // int p = 0, err = 0;
    auto all_cell_checker = [&](const opp_point& point, int& cellIndex) { 
        bool isInside;
        int ci = 0;
        for (ci = 0; ci < cells_set->size; ci++) {
            isPointInCellKernel( 
                            isInside,
                            (const double*)&point, 
                            lc,
                            &((double*)cellVolume_dat->data)[ci * cellVolume_dat->dim], 
                            &((double*)cellDet_dat->data)[ci * cellDet_dat->dim]);
            if (isInside) {       
                cellIndex = ci;
                break;
            }
        }
        // p++;
        // if (ci >= cells_set->size) {
        //     err++;
        //     opp_printf("FUNC", "Error... %d point %d", err, p);
        // }
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
    
    // if (OPP_DBG) 
    if (OPP_rank == 0)
        opp_printf("FUNC", "generateStructMeshToGlobalCellMappings end");
}

//*******************************************************************************
void init_particle_mover(const double gridSpacing, int dim, const opp_dat node_pos_dat, 
    const opp_dat cellVolume_dat, const opp_dat cellDet_dat, const opp_dat global_cell_id_dat) {
    
    opp_profiler->start("SetupMover");

    useGlobalMove = opp_params->get<OPP_BOOL>("opp_global_move");

    if (useGlobalMove) {
        
        boundingBox = std::make_shared<BoundingBox>(node_pos_dat, dim);

        cellMapper = std::make_shared<CellMapper>(boundingBox, gridSpacing, comm);

        generateStructMeshToGlobalCellMappings(cellVolume_dat->set, global_cell_id_dat, 
            cellVolume_dat, cellDet_dat);
    }

    opp_profiler->reg("GlbToLocal");
    opp_profiler->reg("GblMv_Move");
    opp_profiler->reg("GblMv_AllMv");
    for (int i = 0; i < 10; i++) {
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
    args[4] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty

    const int args0_dim = args[0].dim;
    const int args1_dim = args[1].dim;
    const int args2_dim = args[2].dim;
    const int args3_dim = args[3].dim;
    const int c2c_dim   = c2c_map->dim;

    OPP_mesh_relation_data = (OPP_INT*)p2c_map->p2c_dat->data;

    // lambda function for multi hop particle movement
    auto multihop_mover = [&](const int i) {

        opp_p2c = &(OPP_mesh_relation_data)[i];   // TODO : remove OPP_INT* after making this into a map

        if (*opp_p2c == MAX_CELL_INDEX) {
            return;
        }

        OPP_MOVE_RESET_FLAGS;

        do {        
            opp_c2c = &(c2c_map->map)[*opp_p2c * c2c_dim];
            
            getCellIndexKernel(
                &((const OPP_REAL*) args[0].data)[i * args0_dim], 
                &((OPP_REAL*)       args[1].data)[i * args1_dim],
                &((const OPP_REAL*) args[2].data)[*opp_p2c * args2_dim], 
                &((const OPP_REAL*) args[3].data)[*opp_p2c * args3_dim]);

        } while (opp_check_part_move_status(*opp_p2c, i, set->particle_remove_count));
    };

    // ----------------------------------------------------------------------------
    opp_init_particle_move(set, 0, nullptr);

    if (useGlobalMove) {

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
            }
        }

        opp_profiler->end("GblMv_Move");
    }

    // ----------------------------------------------------------------------------
    // check whether all particles not marked for global comm is within cell, 
    // and if not mark to move between cells within the MPI rank, mark for neighbour comm
    opp_profiler->start("Mv_AllMv0");

    const int start1 = OPP_iter_start;
    const int end1 = OPP_iter_end;
    for (int i = start1; i < end1; i++) { 
        
        multihop_mover(i);
    }

    opp_profiler->end("Mv_AllMv0");

    opp_finalize_particle_move(set);

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
    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__get_max set_size %d", set->size);

    opp_profiler->start("GetMax");

    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    const int iter_size = set->size;

    for (int n = 0; n < iter_size; n++)
    {
        get_final_max_values_kernel(
            &((const OPP_REAL*) args[0].data)[n * args[0].dim],     // n_charge_den  
            (OPP_REAL*) args[1].data,
            &((const OPP_REAL*) args[2].data)[n * args[2].dim],     // n_pot  
            (OPP_REAL*) args[3].data
        );
    }

    opp_profiler->end("GetMax");
}
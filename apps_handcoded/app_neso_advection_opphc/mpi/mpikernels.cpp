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
#include "../advec_defs.h"

using namespace opp;

OPP_REAL CONST_extents[2];
OPP_REAL CONST_dt = 0.0;
OPP_REAL CONST_cell_width = 0.0;
OPP_INT CONST_ndimcells[2];

//****************************************
void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    if (!strcmp(name,"CONST_extents"))         std::memcpy(CONST_extents, data, (size*dim));
    else if (!strcmp(name,"CONST_dt"))         std::memcpy(&CONST_dt, data, (size*dim));
    else if (!strcmp(name,"CONST_cell_width")) std::memcpy(&CONST_cell_width, data, (size*dim));
    else if (!strcmp(name,"CONST_ndimcells"))  std::memcpy(CONST_ndimcells, data, (size*dim));
    else std::cerr << "error: unknown const name" << std::endl;
}
//****************************************

#include "../kernels.h"


//*******************************************************************************
void opp_loop_all__verify(
    opp_set set,        // particles_set
    opp_arg arg0,       // part_mesh_rel,        OPP_RW
    opp_arg arg1,       // part_pos,             OPP_READ
    opp_arg arg2,       // cell_global_index,    OPP_READ
    opp_arg arg3        // incorrect_part_count, OPP_INC
)
{
    if (OPP_DBG) 
        opp_printf("ADVEC", "opp_loop_all__verify set_size %d diff %d", set->size, set->diff);

    opp_profiler->start("Verify");

    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);  
    
    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 

    for (int n = 0; n < set_size; n++)
    { 
        const int map0idx = OPP_mesh_relation_data[n];

        verify_kernel( 
            &((OPP_INT*)  args[0].data)[n * args[0].dim],       // part_mesh_rel,      
            &((OPP_REAL*) args[1].data)[n * args[1].dim],       // part_pos,           
            &((OPP_INT*)  args[2].data)[map0idx * args[2].dim], // cell_global_index,  
            (OPP_INT*)    args[3].data                          // incorrect_part_count
        );
    }

    opp_mpi_reduce(&args[3], (int*)args[3].data);

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("Verify");
}

#ifdef FUSE_KERNELS // ############################################################################

//*************************************************************************************************
void opp_particle_move__update_pos_move(
    opp_set set,        // particles_set
    opp_arg arg0,       // part_mesh_rel, OPP_RW
    opp_arg arg1,       // part_vel,      OPP_RW
    opp_arg arg2,       // part_pos,      OPP_RW
    opp_arg arg3,       // cell_pos_ll,   OPP_READ
    opp_arg arg4        // cell_cell_map, OPP_READ
)
{

    if (OPP_DBG) 
        opp_printf("ADVEC", "opp_particle_move__update_pos_move set_size %d diff %d", set->size, set->diff);

    opp_profiler->start("Move");

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
    double kernel_time = 0.0;
    int total_particles = 0;
    int comm_iteration = 0;
    int max_internal_hops = 0, internal_hops = 0;
    std::vector<int> particle_loops_per_comm_iter(10, 0);
    auto total_start = std::chrono::system_clock::now();
    auto start = std::chrono::system_clock::now();
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

    const int nargs = 5;
    opp_arg args[nargs];

    args[0]  = arg0;
    args[1]  = arg1;
    args[2]  = arg2;
    args[3]  = arg3;
    args[4]  = arg4;

    opp_mpi_halo_exchanges(set, nargs, args);
    
    opp_mpi_halo_wait_all(nargs, args);  // unable to overlap computation and communication

    OPP_INT* cellIdx = nullptr;

    do // iterate until all mpi ranks say, I am done
    {
        opp_init_particle_move(set, nargs, args);
        
        if (OPP_DBG) 
            opp_printf("ADVEC", "opp_particle_move__update_pos_move Starting iteration %d, start[%d] end[%d]", 
                OPP_comm_iteration, OPP_iter_start, OPP_iter_end);

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
        start = std::chrono::system_clock::now();
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

        for (int i = OPP_iter_start; i < OPP_iter_end; i++)
        {
            OPP_MOVE_RESET_FLAGS;

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
            internal_hops = 0;
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

            do
            { 
                cellIdx = &(OPP_mesh_relation_data[i]);

                push_particles_kernel( 
                    &((OPP_INT*)  args[0].data)[i * args[0].dim],        // part_cid 
                    &((OPP_REAL*) args[1].data)[i * args[1].dim],        // part_vel 
                    &((OPP_REAL*) args[2].data)[i * args[2].dim],        // part_pos 
                    &((OPP_REAL*) args[3].data)[*cellIdx * args[3].dim], // cell_pos_ll 
                    &((OPP_INT*)  args[4].data)[*cellIdx * args[4].dim]  // cell_cell_map 
                );

            } while (opp_check_part_move_status(*cellIdx, i, set->particle_remove_count));

#ifdef DEBUG_INTERNAL
            if (max_internal_hops < internal_hops) 
                max_internal_hops = internal_hops; // can remove
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------
        }

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
        std::chrono::duration<double> diff   = std::chrono::system_clock::now() - start;
        kernel_time += (double)diff.count();
        total_particles += (OPP_iter_end - OPP_iter_start);
        particle_loops_per_comm_iter[OPP_comm_iteration] = (OPP_iter_end - OPP_iter_start);
        comm_iteration++;
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

    } while (opp_finalize_particle_move(set)); // iterate until all mpi ranks say, I am done

    opp_set_dirtybit(nargs, args);

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
    std::chrono::duration<double> total_diff   = std::chrono::system_clock::now() - total_start;
    opp_printf("Move", "TotalTime: %2.15lE KernelTime: %2.15lE | total_particles: %d | \
        particle_loops_per_comm_iter [%d %d %d %d] | comm_iteration: %d max_internal_hops: %d", 
        (double)total_diff.count(), kernel_time, total_particles, 
        particle_loops_per_comm_iter[0], particle_loops_per_comm_iter[1], particle_loops_per_comm_iter[2], 
        particle_loops_per_comm_iter[3], comm_iteration, max_internal_hops);
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

    opp_profiler->end("Move");
}

#else // ################################################################################################

//*******************************************************************************
void opp_loop_all__UpdatePos(
    opp_set set,        // particles_set
    opp_arg arg0,       // part_vel,      OPP_READ
    opp_arg arg1        // part_pos,      OPP_RW      
)
{
    if (OPP_DBG) 
        opp_printf("ADVEC", "opp_loop_all__UpdatePos set_size %d diff %d", set->size, set->diff);

    opp_profiler->start("UpdatePos");

    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    const int set_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_mpi_halo_wait_all(nargs, args);  // unable to overlap computation and communication

    for (int n = 0; n < set_size; n++)
    { 
        update_pos_kernel( 
            &((OPP_REAL*) args[0].data)[n * args[0].dim],       // part_vel 
            &((OPP_REAL*) args[1].data)[n * args[1].dim]        // part_pos 
        );
    }

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("UpdatePos");
}

// #define DEBUG_INTERNAL

//*******************************************************************************
void opp_particle_mover__Move(
    opp_set set,        // particles_set
    opp_arg arg0,       // part_mesh_rel, OPP_RW
    opp_arg arg1,       // part_pos,      OPP_READ
    opp_arg arg2,       // cell_centroid, OPP_READ
    opp_arg arg3        // cell_cell_map, OPP_READ
) 
{
    if (OPP_DBG) opp_printf("CABANA", "move set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("Move");

    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    const int args0_dim = args[0].dim;
    const int args1_dim = args[1].dim;
    const int args2_dim = args[2].dim;
    const int args3_dim = args[3].dim;

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
    int global_max_comms = 0, internal_comms = 0, global_max_multi = 0, max_multi_hops = 0, max_internal_hops = 0;
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

    // lambda function for multi hop particle movement
    auto multihop_mover = [&](const int i) {

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
        max_internal_hops = 0;
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

        int& cellIdx = ((int*)args[0].data)[i];

        if (cellIdx == MAX_CELL_INDEX) {
            return;
        }

        OPP_MOVE_RESET_FLAGS;

        do {
            is_point_in_current_cell_kernel( 
                &((OPP_INT*)  args[0].data)[i * args0_dim],        // part_cid 
                &((OPP_REAL*) args[1].data)[i * args1_dim],        // part_pos 
                &((OPP_REAL*) args[2].data)[cellIdx * args2_dim], // cell_pos_ll 
                &((OPP_INT*)  args[3].data)[cellIdx * args3_dim]  // cell_cell_map  
            );

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
            max_internal_hops++;
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

        } while (opp_check_part_move_status(cellIdx, i, set->particle_remove_count));

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
        if (max_internal_hops > max_multi_hops) max_multi_hops = max_internal_hops;
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------
    };

    // ----------------------------------------------------------------------------
    opp_init_particle_move(set, 0, nullptr);
    
    if (useGlobalMove) {
        
        globalMover->initGlobalMove();

        opp_profiler->start("GblMv_Move");
        
        opp_point point;
        point.z = 0;

        // check whether particles needs to be moved over global move routine
        for (int i = OPP_iter_start; i < OPP_iter_end; i++) {   
            
            int* cellIdx = &((int*)args[0].data)[i];
            point.x = ((double*)args[1].data)[i * args1_dim + 0];
            point.y = ((double*)args[1].data)[i * args1_dim + 1];

            // const opp_point* point = (const opp_point*)&(((double*)args[1].data)[i * args1_dim]);

            // check for global move, and if satisfy global move criteria, then remove the particle from current rank
            if (opp_part_checkForGlobalMove(set, point, i, *cellIdx)) {
                
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
    const int start1 = OPP_iter_start;
    const int end1 = OPP_iter_end;
    for (int i = start1; i < end1; i++) { 
        
        multihop_mover(i);
    }

    opp_profiler->end("Mv_AllMv0");

    // ----------------------------------------------------------------------------
    // finalize the global move routine and iterate over newly added particles and check whether they need neighbour comm
    if (useGlobalMove && globalMover->finalize(set) > 0) {
        
        opp_profiler->start("GblMv_AllMv");

        // check whether the new particle is within cell, and if not move between cells within the MPI rank, 
        // mark for neighbour comm. Do only for the globally moved particles 
        const int start2 = (set->size - set->diff);
        const int end2 = set->size;
        for (int i = start2; i < end2; i++) { 
                
            multihop_mover(i);                 
        }

        opp_profiler->end("GblMv_AllMv");
    }

    // ----------------------------------------------------------------------------
    // Do neighbour communication and if atleast one particle is received by the currect rank, 
    // then iterate over the newly added particles
    while (opp_finalize_particle_move(set)) {

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
        internal_comms++;
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------

        std::string profName = std::string("Mv_AllMv") + std::to_string(OPP_comm_iteration);
        opp_profiler->start(profName);
        
        opp_init_particle_move(set, 0, nullptr);

        // check whether particle is within cell, and if not move between cells within the MPI rank, mark for neighbour comm
        const int start3 = OPP_iter_start;
        const int end3 = OPP_iter_end;
        for (int i = start3; i < end3; i++) {  
            
            multihop_mover(i);
        }

        opp_profiler->end(profName);
    }

#ifdef DEBUG_INTERNAL // ----------------------------------------------------------------------------
    MPI_Reduce(&internal_comms, &global_max_comms, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&max_multi_hops, &global_max_multi, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    if (OPP_rank == 0) opp_printf("Move", "Max comms %d Max multi %d", global_max_comms, global_max_multi);
#endif // DEBUG_INTERNAL ----------------------------------------------------------------------------
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("Move");
}

//*************************************************************************************************
inline void generate_struct_mesh_to_gbl_cell_maps(const opp_dat cell_pos_ll_dat, const opp_dat global_cid_dat,
                                                    const opp_map cell_cell_map) {
    
    const opp_set cells_set = cell_pos_ll_dat->set;

    // if (OPP_DBG) 
    if (OPP_rank == 0)
        opp_printf("FUNC", "generate_struct_mesh_to_gbl_cell_maps cells [%s] start global grid dimensions %d %d %d",
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

    auto neighbour_cell_checker = [&](const opp_point& point, int& cellIndex) { 
        opp_move_status m;
        cellIndex = 0;
        do {

            m = is_point_in_current_cell_kernel( 
                            &cellIndex,
                            (const OPP_REAL*)&point, 
                            &((OPP_REAL*)cell_pos_ll_dat->data)[cellIndex * cell_pos_ll_dat->dim], 
                            &((OPP_INT*)cell_cell_map->map)[cellIndex * cell_cell_map->dim]);

        } while (m == OPP_NEED_MOVE && cellIndex < cells_set->size); 
    };

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

                neighbour_cell_checker(centroid, cellIndex); // Find in which cell this centroid lies

                if (cellIndex == MAX_CELL_INDEX) {
                    removedCoordinates.insert(std::make_pair(index, opp_point(x, y ,z)));
                    continue;
                }

                if (cellIndex < cells_set->size) { // write only if the structured cell belong to the current MPI rank
                    
                    cellMapper->enrichStructuredMesh(index, ((int*)global_cid_dat->data)[cellIndex], OPP_rank);
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

            neighbour_cell_checker(point, cellIndex);

            if (cellIndex != MAX_CELL_INDEX && (cellIndex < cellSetSizeIncHalo)) { 

                int gblCellIndex = ((int*)global_cid_dat->data)[cellIndex];

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
    cellMapper->convertToLocalMappings(global_cid_dat);
#endif

    // if (OPP_DBG) 
    if (OPP_rank == 0)
        opp_printf("FUNC", "generate_struct_mesh_to_gbl_cell_maps end");
}

//*******************************************************************************
void initialize_particle_mover(const double grid_spacing, int dim, 
    const opp_dat cell_pos_ll_dat, const opp_dat global_cid_dat, const opp_map cell_cell_map) {
// void initializeParticleMover(const double gridSpacing, int dim, const opp_dat node_pos_dat, 
//     const opp_dat cellVolume_dat, const opp_dat cellDet_dat, const opp_dat global_cell_id_dat) {
    
    opp_profiler->start("SetupMover");

    useGlobalMove = opp_params->get<OPP_BOOL>("opp_global_move");

    if (useGlobalMove) {
        
#ifdef USE_MPI
        comm = std::make_shared<Comm>(MPI_COMM_WORLD);
        globalMover = std::make_unique<GlobalParticleMover>(comm->comm_parent);
#endif

        boundingBox = std::make_shared<BoundingBox>(cell_pos_ll_dat, dim, comm);

        cellMapper = std::make_shared<CellMapper>(boundingBox, grid_spacing, comm);

        // generateStructMeshToGlobalCellMappings(cellVolume_dat->set, global_cell_id_dat, 
        //     cellVolume_dat, cellDet_dat);
        generate_struct_mesh_to_gbl_cell_maps(cell_pos_ll_dat, global_cid_dat, cell_cell_map);
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

#endif
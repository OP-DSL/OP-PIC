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

#pragma once

#include <opp_sycl.h>

OPP_INT* OPP_cells_set_size_d = nullptr;
OPP_INT OPP_cells_set_size = -1;

OPP_INT* OPP_comm_iteration_d = nullptr;
OPP_INT OPP_comm_iteration_h = -1;

OPP_INT* cellMapper_pos_stride_d = nullptr;
OPP_INT cellMapper_pos_stride = -1;

OPP_INT* OPP_rank_d = nullptr;
OPP_INT OPP_rank_h = -1;

OPP_INT* OPP_comm_size_d = nullptr;
OPP_INT OPP_comm_size_h = -1;

size_t* opp_maxSavedDHGrid_d = nullptr;
size_t* opp_minSavedDHGrid_d = nullptr;

//*******************************************************************************
// Returns true only if another hop is required by the current rank
inline bool opp_part_check_status_device(char& move_flag, bool& iter_one_flag, 
        int* cell_id, const int particle_index, int* remove_count, int *remove_part_indices, 
        int *move_part_indices, int *move_cell_indices, int *move_count, const OPP_INT* cell_set_size_d)
{
    iter_one_flag = false;

    if (move_flag == OPP_MOVE_DONE) {
        return false;
    }
    else if (move_flag == OPP_NEED_REMOVE) {
        cell_id[0] = MAX_CELL_INDEX;
        const int removeArrayIndex = opp_atomic_fetch_add(remove_count, 1);
        remove_part_indices[removeArrayIndex] = particle_index;

        return false;
    }
    else if (cell_id[0] >= cell_set_size_d[0])
    {
        // cell_id cell is not owned by the current mpi rank, need to communicate
        const int moveArrayIndex = opp_atomic_fetch_add(move_count, 1);
        move_part_indices[moveArrayIndex] = particle_index;
        move_cell_indices[moveArrayIndex] = *cell_id;

        // To be removed from the current rank, packing will be done prior exchange & removal
        move_flag = OPP_NEED_REMOVE; // not required
        const int removeArrayIndex = opp_atomic_fetch_add(remove_count, 1);
        remove_part_indices[removeArrayIndex] = particle_index;

        return false;
    }

    // cell_id is an own cell and move_flag == OPP_NEED_MOVE
    return true;
}

//*******************************************************************************
template <opp_access reduction, int intel, class T, class out_acc, class local_acc>
void opp_reduction(out_acc dat_g, int offset, T dat_l, local_acc temp, sycl::nd_item<1> &item) {
    T dat_t;

    /* important to finish all previous activity */
    item.barrier(sycl::access::fence_space::local_space); 

    size_t tid = item.get_local_id(0);
    temp[tid] = dat_l;

    for (size_t d = item.get_local_range(0) / 2; d > 0; d >>= 1) {
        item.barrier(sycl::access::fence_space::local_space);
        if (tid < d) {
        dat_t = temp[tid + d];

        switch (reduction) {
        case OPP_INC:
            dat_l = dat_l + dat_t;
            break;
        case OPP_MIN:
            if (dat_t < dat_l)
                dat_l = dat_t;
            break;
        case OPP_MAX:
            if (dat_t > dat_l)
                dat_l = dat_t;
            break;
        }
        temp[tid] = dat_l;
        }
    }

    if (tid == 0) {
        switch (reduction) {
        case OPP_INC:
            dat_g[offset] = dat_g[offset] + dat_l;
            break;
        case OPP_MIN:
            if (dat_l < dat_g[offset])
            dat_g[offset] = dat_l;
            break;
        case OPP_MAX:
            if (dat_l > dat_g[offset])
            dat_g[offset] = dat_l;
            break;
        }
    }
}

//*******************************************************************************
inline size_t opp_dev_findStructuredCellIndex2D(const OPP_REAL* pos, const OPP_REAL* oneOverGridSpacing, 
            const OPP_INT* pos_stride, const OPP_REAL* minGlbCoordinate, const size_t* globalGridDims, 
            const size_t* globalGridSize, const size_t* minSavedDHGrid, const size_t* maxSavedDHGrid)
{
    // Round to the nearest integer to minimize rounding errors
    const size_t xIndex = (size_t)((pos[0 * pos_stride[0]] - minGlbCoordinate[0]) * 
                                            oneOverGridSpacing[0]);
    const size_t yIndex = (size_t)((pos[1 * pos_stride[0]] - minGlbCoordinate[1]) * 
                                            oneOverGridSpacing[0]);

    const bool isOutOfCuboid = ((xIndex >= maxSavedDHGrid[0] || xIndex < minSavedDHGrid[0]) ||
                                (yIndex >= maxSavedDHGrid[1] || yIndex < minSavedDHGrid[1]));

    // Calculate the cell index mapping index
    const size_t index = xIndex + (yIndex * globalGridDims[0]);

    return (isOutOfCuboid) ? MAX_CELL_INDEX : index;
}

//*******************************************************************************
inline size_t opp_dev_findStructuredCellIndex3D(const OPP_REAL* pos, const OPP_REAL* oneOverGridSpacing, 
            const OPP_INT* pos_stride, const OPP_REAL* minGlbCoordinate, const size_t* globalGridDims, 
            const size_t* globalGridSize, const size_t* minSavedDHGrid, const size_t* maxSavedDHGrid)
{
    // Round to the nearest integer to minimize rounding errors
    const size_t xIndex = (size_t)((pos[0 * pos_stride[0]] - minGlbCoordinate[0]) * 
                                            oneOverGridSpacing[0]);
    const size_t yIndex = (size_t)((pos[1 * pos_stride[0]] - minGlbCoordinate[1]) * 
                                            oneOverGridSpacing[0]);
    const size_t zIndex = (size_t)((pos[2 * pos_stride[0]] - minGlbCoordinate[2]) * 
                                            oneOverGridSpacing[0]);

    const bool isOutOfCuboid = ((xIndex >= maxSavedDHGrid[0] || xIndex < minSavedDHGrid[0]) ||
                                (yIndex >= maxSavedDHGrid[1] || yIndex < minSavedDHGrid[1]) ||
                                (zIndex >= maxSavedDHGrid[2] || zIndex < minSavedDHGrid[2]));

    // Calculate the cell index mapping index
    const size_t index = xIndex + (yIndex * globalGridDims[0]) + (zIndex * globalGridDims[3]);
    
    // printf("%lf %lf %lf - %lld %lld %lld -- MAX %lld %lld %lld -- MIN %lld %lld %lld -- out %d -- %lld\n", 
    //     pos[0 * pos_stride[0]], pos[1 * pos_stride[0]], pos[2 * pos_stride[0]],
    //     xIndex, yIndex, zIndex, 
    //     maxSavedDHGrid[0],maxSavedDHGrid[1],maxSavedDHGrid[2],
    //     minSavedDHGrid[0],minSavedDHGrid[1],minSavedDHGrid[2],
    //     isOutOfCuboid ? 1 : 0, index);

    return (isOutOfCuboid) ? MAX_CELL_INDEX : index;
}

//*******************************************************************************
inline void opp_dev_remove_dh_particle(OPP_INT& p2c, const OPP_INT part_idx,
            OPP_INT *remove_count, OPP_INT *remove_part_indices) 
{
    p2c = MAX_CELL_INDEX;
    const int removeArrayIndex = opp_atomic_fetch_add(remove_count, 1);
    remove_part_indices[removeArrayIndex] = part_idx;
}

//*******************************************************************************
inline void opp_dev_move_dh_particle(const OPP_INT part_idx, const OPP_INT rank, 
            const OPP_INT cid, OPP_INT *move_count, OPP_INT *arr_part_idx, 
            OPP_INT *arr_cid, OPP_INT *arr_rank) 
{
    const int moveArrayIndex = opp_atomic_fetch_add(move_count, 1);
    arr_part_idx[moveArrayIndex] = part_idx;
    arr_rank[moveArrayIndex] = rank;
    arr_cid[moveArrayIndex] = cid;
}

//*******************************************************************************
inline void opp_dev_checkForGlobalMove(
    const size_t struct_cell_id, const int part_idx, OPP_INT* p2c_map, const int* rank,
    const OPP_INT* struct_mesh_to_cell_map, const OPP_INT* struct_mesh_to_rank_map,
    OPP_INT* remove_count, OPP_INT*  remove_part_indices,
    OPP_INT* move_part_indices, OPP_INT* move_cell_indices, OPP_INT* move_rank_indices, OPP_INT* move_count) 
{
    if (struct_cell_id == MAX_CELL_INDEX) { // This happens when point is out of the unstructured mesh      
        opp_dev_remove_dh_particle(p2c_map[part_idx], part_idx, remove_count, remove_part_indices);
        return;
    }

#ifdef USE_MPI
    const OPP_INT struct_cell_rank = struct_mesh_to_rank_map[struct_cell_id];

    // Check whether the paticles need global moving, if yes start global moving process, 
    // if no, move to the closest local cell
    if (struct_cell_rank != rank[0]) {

        if (struct_cell_rank == MAX_CELL_INDEX) { 
            opp_dev_remove_dh_particle(p2c_map[part_idx], part_idx, remove_count, remove_part_indices);
            return;
        }

        // Due to renumbering local cell indices will be different to global, hence do global comm with global indices
        const OPP_INT unstruct_cell_id = struct_mesh_to_cell_map[struct_cell_id];

        if (unstruct_cell_id == MAX_CELL_INDEX) {
            opp_dev_remove_dh_particle(p2c_map[part_idx], part_idx, remove_count, remove_part_indices);
            return;
        }

        // if the new rank is not the current rank, mark the particle to be sent via global comm
        opp_dev_move_dh_particle(part_idx, struct_cell_rank, unstruct_cell_id, move_count, 
                                    move_part_indices, move_cell_indices, move_rank_indices);
        
        // remove the dh moved particle from the current rank
        opp_dev_remove_dh_particle(p2c_map[part_idx], part_idx, remove_count, remove_part_indices);
        return;
    }
    else
#endif
    {
        int x = struct_mesh_to_cell_map[struct_cell_id];  
        p2c_map[part_idx] = struct_mesh_to_cell_map[struct_cell_id];           
        if (p2c_map[part_idx] == MAX_CELL_INDEX) { // Particle is outside the mesh, need to remove
            // printf("x");
            opp_dev_remove_dh_particle(p2c_map[part_idx], part_idx, remove_count, remove_part_indices);
            return;
        }
    }
}

//*******************************************************************************
inline void opp_dev_checkForGlobalMove_sycl(opp_set set, const OPP_REAL *p_pos, OPP_INT *p2c_map) 
{
    opp_set_stride(cellMapper_pos_stride_d, cellMapper_pos_stride, set->set_capacity);
    opp_set_stride(OPP_rank_d, OPP_rank_h, OPP_rank);

    opp_register_const<size_t>(opp_minSavedDHGrid_d, 3);
    opp_register_const<size_t>(opp_maxSavedDHGrid_d, 3);

    opp_mem::copy_host_to_dev<size_t>(opp_minSavedDHGrid_d, opp_minSavedDHGrid, 3);
    opp_mem::copy_host_to_dev<size_t>(opp_maxSavedDHGrid_d, opp_maxSavedDHGrid, 3);

#ifdef OPP_BLOCK_SIZE_4
    const int block_size = OPP_BLOCK_SIZE_4;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif
    
    const int num_blocks = (OPP_iter_end - OPP_iter_start - 1) / block_size + 1;

    opp_queue->submit([&](sycl::handler &cgh) {

        const OPP_INT* struct_mesh_to_cell_map = cellMapper->structMeshToCellMapping_d;
        const OPP_INT* struct_mesh_to_rank_map = cellMapper->structMeshToRankMapping_d;
        const OPP_REAL* oneOverGridSpacing = cellMapper->oneOverGridSpacing_d;
        const OPP_REAL* minGlbCoordinate = cellMapper->minGlbCoordinate_d;
        const size_t* globalGridDims = cellMapper->globalGridDims_d;
        const size_t* globalGridSize = cellMapper->globalGridSize_d;
        OPP_INT* remove_count = set->particle_remove_count_d;
        OPP_INT* remove_part_indices = OPP_remove_particle_indices_d;
        OPP_INT* move_part_indices = dh_indices_d.part_indices;
        OPP_INT* move_cell_indices = dh_indices_d.cell_indices;
        OPP_INT* move_rank_indices = dh_indices_d.rank_indices;
        OPP_INT* move_count = dh_indices_d.move_count;
        const size_t* minSavedDHGrid = opp_minSavedDHGrid_d;
        const size_t* maxSavedDHGrid = opp_maxSavedDHGrid_d;
        const int iter_start = OPP_iter_start;
        const int iter_end = OPP_iter_end;
        const OPP_INT* pos_stride = cellMapper_pos_stride_d;
        const int* rank = OPP_rank_d;

        const OPP_REAL* p_pos_sycl = p_pos;
        OPP_INT* p2c_map_sycl = p2c_map;

        // -----------------------------------------------------------------------------------------
        auto opp_dh_kernel2D = [=](sycl::nd_item<1> item) {
            
            const int tid = item.get_global_linear_id();
            const int part_idx = tid + iter_start;

            if (part_idx < iter_end) {

                const size_t struct_cell_id = opp_dev_findStructuredCellIndex2D(p_pos_sycl + part_idx, 
                                                oneOverGridSpacing, pos_stride, minGlbCoordinate, 
                                                globalGridDims, globalGridSize, minSavedDHGrid, maxSavedDHGrid);
                
                opp_dev_checkForGlobalMove(struct_cell_id, part_idx, p2c_map_sycl, rank,
                                    struct_mesh_to_cell_map, struct_mesh_to_rank_map,
                                    remove_count, remove_part_indices,
                                    move_part_indices, move_cell_indices, move_rank_indices, move_count);
            }        
        };

        // -----------------------------------------------------------------------------------------
        auto opp_dh_kernel3D = [=](sycl::nd_item<1> item) {
            
            const int tid = item.get_global_linear_id();
            const int part_idx = tid + iter_start;

            if (part_idx < iter_end) {

                const size_t struct_cell_id = opp_dev_findStructuredCellIndex3D(p_pos_sycl + part_idx, 
                                                oneOverGridSpacing, pos_stride, minGlbCoordinate, 
                                                globalGridDims, globalGridSize, minSavedDHGrid, maxSavedDHGrid);
                
                opp_dev_checkForGlobalMove(struct_cell_id, part_idx, p2c_map_sycl, rank,
                                    struct_mesh_to_cell_map, struct_mesh_to_rank_map,
                                    remove_count, remove_part_indices,
                                    move_part_indices, move_cell_indices, move_rank_indices, move_count);
            }        
        };

        // -----------------------------------------------------------------------------------------
        switch (cellMapper->boundingBox->getDim())
        {
        case 2:
            cgh.parallel_for<class opp_dh_kernel2D_sycl>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), opp_dh_kernel2D);
            break;
        case 3:
            cgh.parallel_for<class opp_dh_kernel3D_sycl>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), opp_dh_kernel3D);
            break;        
        default:
            opp_abort("DH is not implemented for the current dimension");
            break;
        }
    });

    OPP_DEVICE_SYNCHRONIZE();
}
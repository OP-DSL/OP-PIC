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

#include <opp_cuda.h>

__constant__ OPP_INT OPP_cells_set_size_d;
OPP_INT OPP_cells_set_size;

__constant__ int OPP_comm_iteration_d;

__constant__ OPP_INT cellMapper_pos_stride_d;
OPP_INT cellMapper_pos_stride = -1;

__constant__ OPP_INT OPP_rank_d;
__constant__ OPP_INT OPP_comm_size_d;

__constant__ size_t opp_maxSavedDHGrid_d[3];
__constant__ size_t opp_minSavedDHGrid_d[3];

//*******************************************************************************
// Returns true only if another hop is required by the current rank
__inline__ __device__ bool opp_part_check_status_device(char& move_flag, bool& iter_one_flag, 
        int* cell_id, int particle_index, int& remove_count, int *remove_part_indices, 
        int *move_part_indices, int *move_cell_indices, int *move_count) 
{
    iter_one_flag = false;

    if (move_flag == OPP_MOVE_DONE)
    {
        return false;
    }
    else if (move_flag == OPP_NEED_REMOVE)
    {
        *cell_id = MAX_CELL_INDEX;
        const int removeArrayIndex = atomicAdd(&remove_count, 1);
        remove_part_indices[removeArrayIndex] = particle_index;

        return false;
    }
    else if (*cell_id >= OPP_cells_set_size_d)
    {
        // cell_id cell is not owned by the current mpi rank, need to communicate
        const int moveArrayIndex = atomicAdd(move_count, 1);
        move_part_indices[moveArrayIndex] = particle_index;
        move_cell_indices[moveArrayIndex] = *cell_id;

        // Needs to be removed from the current rank, 
        // particle packing will be done just prior exchange and removal
        move_flag = OPP_NEED_REMOVE; 
        const int removeArrayIndex = atomicAdd(&remove_count, 1);
        remove_part_indices[removeArrayIndex] = particle_index;

        return false;
    }

    // cell_id is an own cell and move_flag == OPP_NEED_MOVE
    return true;
}

//*******************************************************************************
template <opp_access reduction, class T>
__inline__ __device__ void opp_reduction(volatile T *dat_g, T dat_l) 
{
    extern __shared__ volatile double temp2[];
    __shared__ volatile T *temp;
    temp = (T *)temp2;
    T dat_t;

    __syncthreads(); /* important to finish all previous activity */

    int tid = threadIdx.x;
    temp[tid] = dat_l;

    // first, cope with blockDim.x perhaps not being a power of 2
    __syncthreads();

    int d = 1 << (31 - __clz(((int)blockDim.x - 1)));
    // d = blockDim.x/2 rounded up to nearest power of 2
    if (tid + d < blockDim.x) {
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

    // second, do reductions involving more than one warp
    for (d >>= 1; d > warpSize; d >>= 1) {
        __syncthreads();

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

    // third, do reductions involving just one warp
    __syncthreads();

    if (tid < warpSize) {
        for (; d > 0; d >>= 1) {
            __syncwarp();
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

        // finally, update global reduction variable
        if (tid == 0) {
            switch (reduction) {
            case OPP_INC:
                *dat_g = *dat_g + dat_l;
                break;
            case OPP_MIN:
                if (dat_l < *dat_g)
                *dat_g = dat_l;
                break;
            case OPP_MAX:
                if (dat_l > *dat_g)
                *dat_g = dat_l;
                break;
            }
        }
    }
}

//*******************************************************************************
__inline__ __device__ size_t opp_dev_findStructuredCellIndex2D(const OPP_REAL* pos, 
            const OPP_REAL* oneOverGridSpace, const OPP_REAL* minGlbCoordinate, 
            const size_t* globalGridDims, const size_t* globalGridSize)
{
    // Round to the nearest integer to minimize rounding errors
    const size_t xIndex = (size_t)((pos[0 * cellMapper_pos_stride_d] - minGlbCoordinate[0]) * 
                                            oneOverGridSpace[0]);
    const size_t yIndex = (size_t)((pos[1 * cellMapper_pos_stride_d] - minGlbCoordinate[1]) * 
                                            oneOverGridSpace[0]);

    const bool isOutOfCuboid = ((xIndex >= opp_maxSavedDHGrid_d[0] || xIndex < opp_minSavedDHGrid_d[0]) ||
                                (yIndex >= opp_maxSavedDHGrid_d[1] || yIndex < opp_minSavedDHGrid_d[1]));

    // Calculate the cell index mapping index
    const size_t index = xIndex + (yIndex * globalGridDims[0]);

    return (isOutOfCuboid) ? MAX_CELL_INDEX : index;
}

//*******************************************************************************
__inline__ __device__ size_t opp_dev_findStructuredCellIndex3D(const OPP_REAL* pos, 
            const OPP_REAL* oneOverGridSpace, const OPP_REAL* minGlbCoordinate, 
            const size_t* globalGridDims, const size_t* globalGridSize)
{
    // Round to the nearest integer to minimize rounding errors
    const size_t xIndex = (size_t)((pos[0 * cellMapper_pos_stride_d] - minGlbCoordinate[0]) * 
                                            oneOverGridSpace[0]);
    const size_t yIndex = (size_t)((pos[1 * cellMapper_pos_stride_d] - minGlbCoordinate[1]) * 
                                            oneOverGridSpace[0]);
    const size_t zIndex = (size_t)((pos[2 * cellMapper_pos_stride_d] - minGlbCoordinate[2]) * 
                                            oneOverGridSpace[0]);

    const bool isOutOfCuboid = ((xIndex >= opp_maxSavedDHGrid_d[0] || xIndex < opp_minSavedDHGrid_d[0]) ||
                                (yIndex >= opp_maxSavedDHGrid_d[1] || yIndex < opp_minSavedDHGrid_d[1]) ||
                                (zIndex >= opp_maxSavedDHGrid_d[2] || zIndex < opp_minSavedDHGrid_d[2]));

    // Calculate the cell index mapping index
    const size_t index = xIndex + (yIndex * globalGridDims[0]) + (zIndex * globalGridDims[3]);
    
    // printf("%lf %lf %lf - %lld %lld %lld -- MAX %lld %lld %lld -- MIN %lld %lld %lld -- out %d -- %lld\n", 
    //     pos[0 * cellMapper_pos_stride_d], pos[1 * cellMapper_pos_stride_d], pos[2 * cellMapper_pos_stride_d],
    //     xIndex, yIndex, zIndex, 
    //     opp_maxSavedDHGrid_d[0],opp_maxSavedDHGrid_d[1],opp_maxSavedDHGrid_d[2],
    //     opp_minSavedDHGrid_d[0],opp_minSavedDHGrid_d[1],opp_minSavedDHGrid_d[2],
    //     isOutOfCuboid ? 1 : 0, index);

    return (isOutOfCuboid) ? MAX_CELL_INDEX : index;
}

//*******************************************************************************
__inline__ __device__ void opp_dev_remove_dh_particle(OPP_INT& p2c, const OPP_INT part_idx,
            OPP_INT *__restrict__ remove_count, OPP_INT *__restrict__ remove_part_indices) 
{
    p2c = MAX_CELL_INDEX;
    const int removeArrayIndex = atomicAdd(remove_count, 1);
    remove_part_indices[removeArrayIndex] = part_idx;
}

//*******************************************************************************
__inline__ __device__ void opp_dev_move_dh_particle(const OPP_INT part_idx, const OPP_INT rank, 
            const OPP_INT cid, OPP_INT *__restrict__ move_count, OPP_INT*__restrict__ arr_part_idx, 
            OPP_INT*__restrict__ arr_cid, OPP_INT*__restrict__ arr_rank) 
{
    const int moveArrayIndex = atomicAdd(move_count, 1);
    arr_part_idx[moveArrayIndex] = part_idx;
    arr_rank[moveArrayIndex] = rank;
    arr_cid[moveArrayIndex] = cid;
}

//*******************************************************************************
__inline__ __device__ void opp_dev_checkForGlobalMove(
    const size_t struct_cell_id,
    const int part_idx, 
    OPP_INT* __restrict__ p2c_map,
    const OPP_INT* __restrict__ struct_mesh_to_cell_map, 
    const OPP_INT* __restrict__ struct_mesh_to_rank_map,
    OPP_INT* __restrict__ remove_count, 
    OPP_INT* __restrict__ remove_part_indices,
    OPP_INT* __restrict__ move_part_indices, 
    OPP_INT* __restrict__ move_cell_indices, 
    OPP_INT* __restrict__ move_rank_indices, 
    OPP_INT* __restrict__ move_count) 
{
    if (struct_cell_id == MAX_CELL_INDEX) { // This happens when point is out of the unstructured mesh      
        opp_dev_remove_dh_particle(p2c_map[part_idx], part_idx, remove_count, remove_part_indices);
        return;
    }

#ifdef USE_MPI
    const OPP_INT struct_cell_rank = struct_mesh_to_rank_map[struct_cell_id];

    // Check whether the paticles need global moving, if yes start global moving process, 
    // if no, move to the closest local cell
    if (struct_cell_rank != OPP_rank_d) {

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
        p2c_map[part_idx] = struct_mesh_to_cell_map[struct_cell_id];           
        if (p2c_map[part_idx] == MAX_CELL_INDEX) { // Particle is outside the mesh, need to remove
            opp_dev_remove_dh_particle(p2c_map[part_idx], part_idx, remove_count, remove_part_indices);
            return;
        }
    }
}

//*******************************************************************************
__global__ void opp_dev_checkForGlobalMove2D_kernel(
    const OPP_REAL* __restrict__ p_pos, 
    OPP_INT* __restrict__ p2c_map,
    const OPP_INT* __restrict__ struct_mesh_to_cell_map, 
    const OPP_INT* __restrict__ struct_mesh_to_rank_map,
    const OPP_REAL* __restrict__ oneOverGridSpacing, 
    const OPP_REAL* __restrict__ minGlbCoordinate, 
    const size_t* __restrict__ globalGridDims, 
    const size_t* __restrict__ globalGridSize,
    OPP_INT* __restrict__ remove_count, 
    OPP_INT* __restrict__ remove_part_indices,
    OPP_INT* __restrict__ move_part_indices, 
    OPP_INT* __restrict__ move_cell_indices, 
    OPP_INT* __restrict__ move_rank_indices, 
    OPP_INT* __restrict__ move_count, 
    const OPP_INT start, 
    const OPP_INT end) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    const int part_idx = thread_id + start;  

    if (part_idx < end) {
        const size_t struct_cell_id = opp_dev_findStructuredCellIndex2D(p_pos + part_idx,  
                                                    oneOverGridSpacing, minGlbCoordinate, 
                                                    globalGridDims, globalGridSize);
        
        opp_dev_checkForGlobalMove(struct_cell_id, part_idx, p2c_map, struct_mesh_to_cell_map, 
            struct_mesh_to_rank_map, remove_count, remove_part_indices, move_part_indices, 
            move_cell_indices, move_rank_indices, move_count);
    }
}

//*******************************************************************************
__global__ void opp_dev_checkForGlobalMove3D_kernel(
    const OPP_REAL* __restrict__ p_pos, 
    OPP_INT* __restrict__ p2c_map,
    const OPP_INT* __restrict__ struct_mesh_to_cell_map, 
    const OPP_INT* __restrict__ struct_mesh_to_rank_map,
    const OPP_REAL* __restrict__ oneOverGridSpacing, 
    const OPP_REAL* __restrict__ minGlbCoordinate, 
    const size_t* __restrict__ globalGridDims, 
    const size_t* __restrict__ globalGridSize,
    OPP_INT* __restrict__ remove_count, 
    OPP_INT* __restrict__ remove_part_indices,
    OPP_INT* __restrict__ move_part_indices, 
    OPP_INT* __restrict__ move_cell_indices, 
    OPP_INT* __restrict__ move_rank_indices, 
    OPP_INT* __restrict__ move_count, 
    const OPP_INT start, 
    const OPP_INT end
) {
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    const int part_idx = thread_id + start;  

    if (part_idx < end) {
        const size_t struct_cell_id = opp_dev_findStructuredCellIndex3D(p_pos + part_idx,  
                                                    oneOverGridSpacing, minGlbCoordinate, 
                                                    globalGridDims, globalGridSize);
        
        opp_dev_checkForGlobalMove(struct_cell_id, part_idx, p2c_map, struct_mesh_to_cell_map, 
            struct_mesh_to_rank_map, remove_count, remove_part_indices, move_part_indices, 
            move_cell_indices, move_rank_indices, move_count);
    }
}

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

#include "hip/hip_runtime.h"

//*******************************************************************************
void init_particle_mover(const double gridSpacing, int dim, const opp_dat node_pos_dat, 
    const opp_dat cellVolume_dat, const opp_dat cellDet_dat, const opp_dat global_cell_id_dat) {
    
    opp_profiler->start("SetupMover");

    useGlobalMove = opp_params->get<OPP_BOOL>("opp_global_move");
    useGlobalMove = false; // Comment this if direct move is implemented on GPU

    // if (useGlobalMove) {
        
    //     boundingBox = std::make_shared<BoundingBox>(node_pos_dat, dim, comm);

    //     cellMapper = std::make_shared<CellMapper>(boundingBox, gridSpacing, comm);

    //     generateStructMeshToGlobalCellMappings(cellVolume_dat->set, global_cell_id_dat, 
    //         cellVolume_dat, cellDet_dat);
    // }

    // opp_profiler->reg("GlbToLocal");
    // opp_profiler->reg("GblMv_Move");
    // opp_profiler->reg("GblMv_AllMv");
    // for (int i = 0; i < 10; i++) {
    //     std::string profName = std::string("Mv_AllMv") + std::to_string(i);
    //     opp_profiler->reg(profName);
    // }

    opp_profiler->end("SetupMover");
}

int move_stride_OPP_HOST_0 = -1;
int move_stride_OPP_HOST_2 = -1;
int move_stride_OPP_HOST_4 = -1;
int move_stride_OPP_HOST_5 = -1;

__constant__ int move_stride_OPP_DEVICE_0;
__constant__ int move_stride_OPP_DEVICE_2;
__constant__ int move_stride_OPP_DEVICE_4;
__constant__ int move_stride_OPP_DEVICE_5;

__constant__ double C_ONE_OVER_SIX = (1.0 / 6.0);

//user function
//*************************************************************************************************
__device__ void move_all_particles_to_cell__kernel(
    opp_move_var& m,
    const double *part_pos,
    int* cell_index,
    double *part_lc,   
    const double *cell_volume,
    const double *cell_det,
    const int *cell_connectivity
)
{
    bool inside = true;
    double coefficient2 = C_ONE_OVER_SIX / (*cell_volume) ;

    part_lc[0 * move_stride_OPP_DEVICE_2] = coefficient2 * (
        cell_det[(0 * DET_FIELDS + 0) * move_stride_OPP_DEVICE_4] - 
        cell_det[(0 * DET_FIELDS + 1) * move_stride_OPP_DEVICE_4] * part_pos[0 * move_stride_OPP_DEVICE_0] + 
        cell_det[(0 * DET_FIELDS + 2) * move_stride_OPP_DEVICE_4] * part_pos[1 * move_stride_OPP_DEVICE_0] - 
        cell_det[(0 * DET_FIELDS + 3) * move_stride_OPP_DEVICE_4] * part_pos[2 * move_stride_OPP_DEVICE_0]);

    part_lc[1 * move_stride_OPP_DEVICE_2] = coefficient2 * (
        cell_det[(1 * DET_FIELDS + 0) * move_stride_OPP_DEVICE_4] - 
        cell_det[(1 * DET_FIELDS + 1) * move_stride_OPP_DEVICE_4] * part_pos[0 * move_stride_OPP_DEVICE_0] + 
        cell_det[(1 * DET_FIELDS + 2) * move_stride_OPP_DEVICE_4] * part_pos[1 * move_stride_OPP_DEVICE_0] - 
        cell_det[(1 * DET_FIELDS + 3) * move_stride_OPP_DEVICE_4] * part_pos[2 * move_stride_OPP_DEVICE_0]);
    
    part_lc[2 * move_stride_OPP_DEVICE_2] = coefficient2 * (
        cell_det[(2 * DET_FIELDS + 0) * move_stride_OPP_DEVICE_4] - 
        cell_det[(2 * DET_FIELDS + 1) * move_stride_OPP_DEVICE_4] * part_pos[0 * move_stride_OPP_DEVICE_0] + 
        cell_det[(2 * DET_FIELDS + 2) * move_stride_OPP_DEVICE_4] * part_pos[1 * move_stride_OPP_DEVICE_0] - 
        cell_det[(2 * DET_FIELDS + 3) * move_stride_OPP_DEVICE_4] * part_pos[2 * move_stride_OPP_DEVICE_0]);
    
    part_lc[3 * move_stride_OPP_DEVICE_2] = coefficient2 * (
        cell_det[(3 * DET_FIELDS + 0) * move_stride_OPP_DEVICE_4] - 
        cell_det[(3 * DET_FIELDS + 1) * move_stride_OPP_DEVICE_4] * part_pos[0 * move_stride_OPP_DEVICE_0] + 
        cell_det[(3 * DET_FIELDS + 2) * move_stride_OPP_DEVICE_4] * part_pos[1 * move_stride_OPP_DEVICE_0] - 
        cell_det[(3 * DET_FIELDS + 3) * move_stride_OPP_DEVICE_4] * part_pos[2 * move_stride_OPP_DEVICE_0]);

    if (part_lc[0 * move_stride_OPP_DEVICE_2] < 0.0 || 
        part_lc[0 * move_stride_OPP_DEVICE_2] > 1.0 ||
        part_lc[1 * move_stride_OPP_DEVICE_2] < 0.0 || 
        part_lc[1 * move_stride_OPP_DEVICE_2] > 1.0 ||
        part_lc[2 * move_stride_OPP_DEVICE_2] < 0.0 || 
        part_lc[2 * move_stride_OPP_DEVICE_2] > 1.0 ||
        part_lc[3 * move_stride_OPP_DEVICE_2] < 0.0 || 
        part_lc[3 * move_stride_OPP_DEVICE_2] > 1.0) 
            inside = false;   

    if (inside) {
        m.move_status = OPP_MOVE_DONE;
        return;
    }

    // outside the last known cell, find most negative weight and 
    // use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = part_lc[0 * move_stride_OPP_DEVICE_2];
    
    for (int i=1; i<NEIGHB_C; i++)
    {
        if (part_lc[i * move_stride_OPP_DEVICE_2] < min_lc) 
        {
            min_lc = part_lc[i * move_stride_OPP_DEVICE_2];
            min_i = i;
        }
    }

    if (cell_connectivity[min_i * move_stride_OPP_DEVICE_5] >= 0) // is there a neighbor in this direction?
    {
        (*cell_index) = cell_connectivity[min_i * move_stride_OPP_DEVICE_5];
        m.move_status = OPP_NEED_MOVE;
    }
    else
    {
        (*cell_index) = MAX_CELL_INDEX;
        m.move_status = OPP_NEED_REMOVE;
    }
}

//*******************************************************************************
// Returns true only if another hop is required by the current rank
__device__ bool opp_part_check_status_device(opp_move_var& m, int* map0idx, int particle_index, 
        int& remove_count, int *remove_particle_indices, int *move_particle_indices, int *move_cell_indices, int *move_count) 
{
    m.iteration_one = false;

    if (m.move_status == OPP_MOVE_DONE)
    {
        return false;
    }
    else if (m.move_status == OPP_NEED_REMOVE)
    {
        *map0idx = MAX_CELL_INDEX;
        int removeArrayIndex = atomicAdd(&remove_count, 1);
        remove_particle_indices[removeArrayIndex] = particle_index;

        return false;
    }
    else if (*map0idx >= OPP_cells_set_size_d)
    {
        // map0idx cell is not owned by the current mpi rank (it is in the import exec halo region), need to communicate
        int moveArrayIndex = atomicAdd(move_count, 1);
        move_particle_indices[moveArrayIndex] = particle_index;
        move_cell_indices[moveArrayIndex] = *map0idx;

        // Needs to be removed from the current rank, bdw particle packing will be done just prior exchange and removal
        m.move_status = OPP_NEED_REMOVE; 
        int removeArrayIndex = atomicAdd(&remove_count, 1);
        remove_particle_indices[removeArrayIndex] = particle_index;

        return false;
    }

    // map0idx is an own cell and m.move_status == OPP_NEED_MOVE
    return true;
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void opp_device_all_MoveToCells(
    int *__restrict d_cell_index,
    const double *__restrict dir_arg0,      // part_pos,
    int *__restrict dir_arg1,               // cell_index,
    double *__restrict dir_arg2,            // part_lc,
    const double *__restrict ind_arg3,      // cell_volume,
    const double *__restrict ind_arg4,      // cell_det,
    const int *__restrict ind_arg5,         // cell_connectivity,
    int *__restrict particle_remove_count,
    int *__restrict particle_remove_indices,
    int *__restrict move_particle_indices,
    int *__restrict move_cell_indices,
    int *__restrict move_count,
    int start,
    int end) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        opp_move_var m;
        m.iteration_one = (OPP_comm_iteration_d > 0) ? false : true;
        int* map0idx = nullptr; //MAX_CELL_INDEX;

        do
        {
            map0idx = &(dir_arg1[n]);
            // int map0idx = d_cell_index[n]; // TODO : I dont know why this isn't working ??? 
            // dir_arg2 and d_cell_index has same pointer values, but this get stuck!

            //user-supplied kernel call
            move_all_particles_to_cell__kernel(
                (m),
                (dir_arg0 + n),          // part_pos,
                (dir_arg1 + n),          // cell_index,
                (dir_arg2 + n),          // part_lc,
                (ind_arg3 + *map0idx),   // cell_volume,
                (ind_arg4 + *map0idx),   // cell_det,
                (ind_arg5 + *map0idx)    // cell_connectivity,
            );                

        } while (opp_part_check_status_device(m, map0idx, n, 
            *particle_remove_count, particle_remove_indices, move_particle_indices, move_cell_indices, move_count));
    }
}

//*******************************************************************************
void opp_particle_move__move(
    opp_set set,        // particle_set,
    opp_arg arg0,       // part_position,   
    opp_arg arg1,       // part_mesh_rel,   
    opp_arg arg2,       // part_lc,         
    opp_arg arg3,       // cell_volume,     
    opp_arg arg4,       // cell_det,        
    opp_arg arg5        // cell_v_cell_map
) 
{ 

    if (FP_DEBUG) opp_printf("FEMPIC", "opp_particle_move__move set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("MoveToCells");

    int nargs = 6;
    opp_arg args[nargs];

    args[0]  = std::move(arg0);
    args[1]  = std::move(arg1);
    args[2]  = std::move(arg2);
    args[3]  = std::move(arg3);
    args[4]  = std::move(arg4);
    args[5]  = std::move(arg5);

    opp_profiler->start("Mv_halo_exchanges");    
    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU); 
    opp_mpi_halo_wait_all(nargs, args);
    opp_profiler->end("Mv_halo_exchanges");

    if (set_size > 0) 
    {
        do {

            move_stride_OPP_HOST_0 = args[0].dat->set->set_capacity;
            move_stride_OPP_HOST_2 = args[2].dat->set->set_capacity;
            move_stride_OPP_HOST_4 = args[4].dat->set->set_capacity; 
            move_stride_OPP_HOST_5 = args[5].size;
            OPP_cells_set_size = set->cells_set->size; 
         
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(OPP_cells_set_size_d), 
                                                        &OPP_cells_set_size, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(move_stride_OPP_DEVICE_0), 
                                                        &move_stride_OPP_HOST_0, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(move_stride_OPP_DEVICE_2), 
                                                        &move_stride_OPP_HOST_2, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(move_stride_OPP_DEVICE_4), 
                                                        &move_stride_OPP_HOST_4, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(move_stride_OPP_DEVICE_5), 
                                                        &move_stride_OPP_HOST_5, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(OPP_comm_iteration_d), 
                                                        &OPP_comm_iteration, sizeof(int)));

            opp_profiler->start("Mv_init_part");
            opp_init_particle_move(set, nargs, args);
            opp_profiler->end("Mv_init_part");

            if (OPP_iter_end - OPP_iter_start > 0) 
            {

                if (OPP_DBG) 
                    opp_printf("MOVE", "iter %d start %d end %d", OPP_comm_iteration, OPP_iter_start, OPP_iter_end);

                int nthread = OPP_gpu_threads_per_block;
                int nblocks = (OPP_iter_end - OPP_iter_start - 1) / nthread + 1;

                cutilSafeCall(hipDeviceSynchronize());
                opp_profiler->start("Mv_OnlyMoveKernel");
                
                opp_device_all_MoveToCells<<<nblocks, nthread>>>(
                    (int *)           set->mesh_relation_dat->data_d,
                    (const double *)  args[0].data_d,                   // part_position,   
                    (int *)           args[1].data_d,                   // part_mesh_rel,   
                    (double *)        args[2].data_d,                   // part_lc,         
                    (const double *)  args[3].data_d,                   // cell_volume,     
                    (const double *)  args[4].data_d,                   // cell_det,        
                    (const int *)     args[5].data_d,                   // cell_v_cell_map
                    (int *)           set->particle_remove_count_d,
                    (int *)           OPP_remove_particle_indices_d,
                    (int*)            OPP_move_particle_indices_d,
                    (int*)            OPP_move_cell_indices_d,
                    (int*)            OPP_move_count_d,
                    OPP_iter_start, 
                    OPP_iter_end);

                cutilSafeCall(hipDeviceSynchronize());
                opp_profiler->end("Mv_OnlyMoveKernel");

            }

        } while (opp_finalize_particle_move(set)); 

    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);

    opp_profiler->end("MoveToCells");
}

//*************************************************************************************************
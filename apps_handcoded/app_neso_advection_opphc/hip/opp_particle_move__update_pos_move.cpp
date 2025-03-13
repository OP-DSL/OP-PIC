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

int move_stride_OPP_HOST_1 = -1;
int move_stride_OPP_HOST_2 = -1;
int move_stride_OPP_HOST_3 = -1;
int move_stride_OPP_HOST_4 = -1;

__constant__ int move_stride_OPP_DEVICE_1;
__constant__ int move_stride_OPP_DEVICE_2;
__constant__ int move_stride_OPP_DEVICE_3;
__constant__ int move_stride_OPP_DEVICE_4;

//user function
//*************************************************************************************************
__device__ void move_all_particles_to_cell__kernel(opp_move_var& m, 
    OPP_INT* part_cid, 
    OPP_REAL* part_vel, 
    OPP_REAL* part_pos, 
    const OPP_REAL* cell_pos_ll, 
    const OPP_INT* cell_cell_map
)
{
    if (m.iteration_one) {
        
        for (int dm = 0; dm < DIM; dm++) {
            
            part_pos[dm * move_stride_OPP_DEVICE_2] += part_vel[dm * move_stride_OPP_DEVICE_1] * CONST_DEVICE_dt; // s1 = s0 + ut
            
            // correct for periodic boundary conditions
            const OPP_INT n_extent_offset_int = std::abs(part_pos[dm * move_stride_OPP_DEVICE_2]) + 2.0;
            const OPP_REAL temp_pos = part_pos[dm * move_stride_OPP_DEVICE_2] + n_extent_offset_int * CONST_DEVICE_extents[dm];
            part_pos[dm * move_stride_OPP_DEVICE_2] = std::fmod(temp_pos, CONST_DEVICE_extents[dm]);
        }
    }

    // check for x direction movement
    const OPP_REAL part_pos_x = part_pos[Dim::x * move_stride_OPP_DEVICE_2];
    if (part_pos_x < cell_pos_ll[Dim::x * move_stride_OPP_DEVICE_3]) {
        part_cid[0] = cell_cell_map[CellMap::xd_y * move_stride_OPP_DEVICE_4];
        m.move_status = OPP_NEED_MOVE;
        return;
    }
    if (part_pos_x > (cell_pos_ll[Dim::x * move_stride_OPP_DEVICE_3] + CONST_DEVICE_cell_width)) {
        part_cid[0] = cell_cell_map[CellMap::xu_y * move_stride_OPP_DEVICE_4];
        m.move_status = OPP_NEED_MOVE;
        return;
    }

    // check for y direction movement
    const OPP_REAL part_pos_y = part_pos[Dim::y * move_stride_OPP_DEVICE_2];
    if (part_pos_y < cell_pos_ll[Dim::y * move_stride_OPP_DEVICE_3]) {
        part_cid[0] = cell_cell_map[CellMap::x_yd * move_stride_OPP_DEVICE_4];
        m.move_status = OPP_NEED_MOVE;
        return;
    }
    if (part_pos_y > (cell_pos_ll[Dim::y * move_stride_OPP_DEVICE_3] + CONST_DEVICE_cell_width)) {
        part_cid[0] = cell_cell_map[CellMap::x_yu * move_stride_OPP_DEVICE_4];
        m.move_status = OPP_NEED_MOVE;
        return;
    }

    m.move_status = OPP_MOVE_DONE;
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

// CUDA kernel function
//*************************************************************************************************
__global__ void opp_device_all_Move(
    OPP_INT *__restrict dir_arg0,           // part_cid 
    OPP_REAL *__restrict dir_arg1,          // part_vel 
    OPP_REAL *__restrict dir_arg2,          // part_pos 
    const OPP_REAL *__restrict ind_arg3,    // cell_pos_ll 
    const OPP_INT *__restrict ind_arg4,     // cell_cell_map
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

        OPP_INT* cellIdx = nullptr; //MAX_CELL_INDEX;

        do
        {
            cellIdx = &(dir_arg0[n]);

            //user-supplied kernel call
            move_all_particles_to_cell__kernel(
                (m),
                (dir_arg0 + n),          // part_cid 
                (dir_arg1 + n),          // part_vel 
                (dir_arg2 + n),          // part_pos 
                (ind_arg3 + *cellIdx),   // cell_pos_ll 
                (ind_arg4 + *cellIdx)    // cell_cell_map
            );                

        } while (opp_part_check_status_device(m, cellIdx, n, 
            *particle_remove_count, particle_remove_indices, move_particle_indices, move_cell_indices, move_count));
    }
}

//*******************************************************************************
void opp_particle_move__update_pos_move(
    opp_set set,        // particles_set
    opp_arg arg0,       // part_mesh_rel, OPP_RW
    opp_arg arg1,       // part_vel,      OPP_RW
    opp_arg arg2,       // part_pos,      OPP_RW
    opp_arg arg3,       // cell_pos_ll,   OPP_READ
    opp_arg arg4        // cell_cell_map, OPP_READ
) 
{ 

    if (OPP_DBG) opp_printf("ADVEC", "opp_particle_move__update_pos_move set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("Move");

    int nargs = 5;
    opp_arg args[nargs];

    args[0]  = std::move(arg0);
    args[1]  = std::move(arg1);
    args[2]  = std::move(arg2);
    args[3]  = std::move(arg3);
    args[4]  = std::move(arg4);

    opp_profiler->start("FMv_halo_exchanges");    
    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU); 
    
    opp_mpi_halo_wait_all(nargs, args);  // unable to overlap computation and communication
    opp_profiler->end("FMv_halo_exchanges");

    if (set_size > 0) 
    {
        do {

            move_stride_OPP_HOST_1 = args[1].dat->set->set_capacity;
            move_stride_OPP_HOST_2 = args[2].dat->set->set_capacity;
            move_stride_OPP_HOST_3 = args[3].dat->set->set_capacity; 
            move_stride_OPP_HOST_4 = args[4].size; 
            OPP_cells_set_size = set->cells_set->size; 
         
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(OPP_cells_set_size_d), 
                                                        &OPP_cells_set_size, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(move_stride_OPP_DEVICE_1), 
                                                        &move_stride_OPP_HOST_1, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(move_stride_OPP_DEVICE_2), 
                                                        &move_stride_OPP_HOST_2, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(move_stride_OPP_DEVICE_3), 
                                                        &move_stride_OPP_HOST_3, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(move_stride_OPP_DEVICE_4), 
                                                        &move_stride_OPP_HOST_4, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(OPP_comm_iteration_d), 
                                                        &OPP_comm_iteration, sizeof(int)));

            opp_profiler->start("FMv_init_part");
            opp_init_particle_move(set, nargs, args);
            opp_profiler->end("FMv_init_part");

            if (OPP_iter_end - OPP_iter_start > 0) 
            {

                if (OPP_DBG || OPP_comm_iteration > 3)
                    opp_printf("MOVE", "iter %d start %d end %d - COUNT=%d", OPP_comm_iteration, 
                                    OPP_iter_start, OPP_iter_end, (OPP_iter_end - OPP_iter_start));

                int nthread = OPP_gpu_threads_per_block;
                int nblocks = (OPP_iter_end - OPP_iter_start - 1) / nthread + 1;

                cutilSafeCall(hipDeviceSynchronize());
                opp_profiler->start("FMv_OnlyMoveKernel");
                
                opp_device_all_Move<<<nblocks, nthread>>>(
                    (OPP_INT*)        args[0].data_d,         // part_cid 
                    (OPP_REAL*)       args[1].data_d,         // part_vel 
                    (OPP_REAL*)       args[2].data_d,         // part_pos 
                    (const OPP_REAL*) args[3].data_d,         // cell_pos_ll 
                    (const OPP_INT*)  args[4].data_d,         // cell_cell_map
                    (int *)           set->particle_remove_count_d,
                    (int *)           OPP_remove_particle_indices_d,
                    (int*)            OPP_move_particle_indices_d,
                    (int*)            OPP_move_cell_indices_d,
                    (int*)            OPP_move_count_d,
                    OPP_iter_start, 
                    OPP_iter_end);

                cutilSafeCall(hipDeviceSynchronize());
                opp_profiler->end("FMv_OnlyMoveKernel");

            }

        } while (opp_finalize_particle_move(set)); 
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);

    opp_profiler->end("Move");
}

//*************************************************************************************************
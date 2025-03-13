// /* 
// BSD 3-Clause License

// Copyright (c) 2022, OP-DSL

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// */

// //*********************************************
// // AUTO GENERATED CODE
// //*********************************************


// int move_stride_OPP_HOST_0 = -1;
// int move_stride_OPP_HOST_1 = -1;
// int move_stride_OPP_HOST_2 = -1;
// int move_stride_OPP_HOST_3 = -1;
// int move_stride_OPP_HOST_6 = -1;
// int move_stride_OPP_HOST_7 = -1;
// int move_stride_OPP_HOST_8 = -1;
// int move_size_OPP_HOST_8   = -1;

// __constant__ int move_stride_OPP_CUDA_0;
// __constant__ int move_stride_OPP_CUDA_1;
// __constant__ int move_stride_OPP_CUDA_2;
// __constant__ int move_stride_OPP_CUDA_3;
// __constant__ int move_stride_OPP_CUDA_6;
// __constant__ int move_stride_OPP_CUDA_7;
// __constant__ int move_stride_OPP_CUDA_8;
// __constant__ int move_size_OPP_CUDA_8;
// __constant__ double ONE_OVER_SIX = (1.0 / 6.0);

// __constant__ bool use_shared_OPP_CUDA = true;
// __constant__ int gpu_threads_per_block_OPP_CUDA;

// //user function
// //*************************************************************************************************
// __device__ void move_all_particles_to_cell_iteration1__kernel(
//     const double *cell_ef,
//     double *part_pos,
//     double *part_vel
// )
// {
//     double coefficient1 = CONST_charge_cuda / CONST_mass_cuda * (CONST_dt_cuda);

//     part_vel[0 * move_stride_OPP_CUDA_2] += 
//         (coefficient1 * cell_ef[0 * move_stride_OPP_CUDA_0]);           
//     part_vel[1 * move_stride_OPP_CUDA_2] += 
//         (coefficient1 * cell_ef[1 * move_stride_OPP_CUDA_0]);   
//     part_vel[2 * move_stride_OPP_CUDA_2] += 
//         (coefficient1 * cell_ef[2 * move_stride_OPP_CUDA_0]);   

//     part_pos[0 * move_stride_OPP_CUDA_1] += 
//         part_vel[0 * move_stride_OPP_CUDA_2] * (CONST_dt_cuda); // v = u + at
//     part_pos[1 * move_stride_OPP_CUDA_1] += 
//         part_vel[1 * move_stride_OPP_CUDA_2] * (CONST_dt_cuda); // v = u + at
//     part_pos[2 * move_stride_OPP_CUDA_1] += 
//         part_vel[2 * move_stride_OPP_CUDA_2] * (CONST_dt_cuda); // v = u + at
// }

// __device__ void move_all_particles_to_cell__kernel(
//     int& move_status,
//     double *part_pos,
//     double *part_lc,
//     int* cell_index,
//     const double *cell_volume,
//     const double *cell_det,
//     const int *cell_connectivity,
//     double *node_charge_den0,
//     double *node_charge_den1,
//     double *node_charge_den2,
//     double *node_charge_den3
// )
// {
//     bool inside = true;
//     double coefficient2 = ONE_OVER_SIX / (*cell_volume) ;
//     for (int i=0; i<N_PER_C; i++) /*loop over vertices*/
//     {
//         part_lc[i * move_stride_OPP_CUDA_3] = coefficient2 * (
//             cell_det[(i * DET_FIELDS + 0) * move_stride_OPP_CUDA_6] - 
//             cell_det[(i * DET_FIELDS + 1) * move_stride_OPP_CUDA_6] * part_pos[0 * move_stride_OPP_CUDA_1] + 
//             cell_det[(i * DET_FIELDS + 2) * move_stride_OPP_CUDA_6] * part_pos[1 * move_stride_OPP_CUDA_1] - 
//             cell_det[(i * DET_FIELDS + 3) * move_stride_OPP_CUDA_6] * part_pos[2 * move_stride_OPP_CUDA_1]
//                 );
        
//         if (part_lc[i * move_stride_OPP_CUDA_3]<0 || 
//             part_lc[i * move_stride_OPP_CUDA_3]>1.0) 
//                 inside = false;
//     }    

//     if (inside)
//     {
//         move_status = OPP_MOVE_DONE;

//         atomicAdd(node_charge_den0, (part_lc[0 * move_stride_OPP_CUDA_3]));
//         atomicAdd(node_charge_den1, (part_lc[1 * move_stride_OPP_CUDA_3]));
//         atomicAdd(node_charge_den2, (part_lc[2 * move_stride_OPP_CUDA_3]));
//         atomicAdd(node_charge_den3, (part_lc[3 * move_stride_OPP_CUDA_3]));

//         return;
//     }

//     // outside the last known cell, find most negative weight and 
//     // use that cell_index to reduce computations
//     int min_i = 0;
//     double min_lc = part_lc[0 * move_stride_OPP_CUDA_3];
    
//     for (int i=1; i<NEIGHB_C; i++)
//     {
//         if (part_lc[i * move_stride_OPP_CUDA_3] < min_lc) 
//         {
//             min_lc = part_lc[i * move_stride_OPP_CUDA_3];
//             min_i = i;
//         }
//     }

//     if (cell_connectivity[min_i * move_stride_OPP_CUDA_7] >= 0) // is there a neighbor in this direction?
//     {
//         (*cell_index) = cell_connectivity[min_i * move_stride_OPP_CUDA_7];
//         move_status = OPP_NEED_MOVE;
//     }
//     else
//     {
//         move_status = OPP_NEED_REMOVE;
//     }
// }

// // CUDA kernel function
// //*************************************************************************************************
// __global__ void opp_cuda_all_MoveToCells(
//     int *__restrict d_cell_index,
//     double *__restrict ind_arg0,
//     double *__restrict dir_arg1,
//     double *__restrict dir_arg2,
//     double *__restrict dir_arg3,
//     int *__restrict dir_arg4,
//     const double *__restrict ind_arg5,
//     const double *__restrict ind_arg6,
//     const int *__restrict ind_arg7,
//     double *__restrict ind_arg8,
//     const int *__restrict opDat8Map,
//     double *__restrict ind_arg9,
//     double *__restrict ind_arg10,
//     double *__restrict ind_arg11,
//     int *__restrict particle_remove_count,
//     int start,
//     int end) 
// {
//     extern __shared__ double shared_arg8[];
//     double* k_arg8 = ind_arg8;

//     int tid = threadIdx.x + blockIdx.x * blockDim.x;

//     if (use_shared_OPP_CUDA)
//     {
//         for (int i = threadIdx.x; i < move_size_OPP_CUDA_8; i+=gpu_threads_per_block_OPP_CUDA) 
//             shared_arg8[i] = 0;
//         __syncthreads();

//         k_arg8 = shared_arg8;
//     }

//     if (tid + start < end) 
//     {
//         int n = tid + start;

//         int& map0idx = dir_arg4[n];
//         // int map0idx = d_cell_index[n]; // TODO : I dont know why this isn't working ??? 
//         // dir_arg2 and d_cell_index has same pointer values, but this get stuck!

//         move_all_particles_to_cell_iteration1__kernel(
//             (ind_arg0 + map0idx),   // cell_ef,
//             (dir_arg1 + n),         // part_pos,
//             (dir_arg2 + n)          // part_vel,
//         );

//         int move_status = OPP_MOVE_DONE;

//         do
//         {
//             const int map1idx = opDat8Map[map0idx + move_stride_OPP_CUDA_8 * 0];
//             const int map2idx = opDat8Map[map0idx + move_stride_OPP_CUDA_8 * 1];
//             const int map3idx = opDat8Map[map0idx + move_stride_OPP_CUDA_8 * 2];
//             const int map4idx = opDat8Map[map0idx + move_stride_OPP_CUDA_8 * 3];

//             move_all_particles_to_cell__kernel(
//                 (move_status),
//                 (dir_arg1 + n),         // part_pos,
//                 (dir_arg3 + n),         // part_lc,
//                 (dir_arg4 + n),         // cell_index,
//                 (ind_arg5 + map0idx),   // cell_volume,
//                 (ind_arg6 + map0idx),   // cell_det,
//                 (ind_arg7 + map0idx),   // cell_connectivity,
//                 (k_arg8 + map1idx),
//                 (k_arg8 + map2idx),
//                 (k_arg8 + map3idx),
//                 (k_arg8 + map4idx)
//             );                

//         } while (move_status == (int)OPP_NEED_MOVE);

//         if (move_status == OPP_NEED_REMOVE) 
//         {
//             atomicAdd(particle_remove_count, 1);
//             dir_arg4[n] = MAX_CELL_INDEX;
//         }
//     }

//     __syncthreads();

//     if (use_shared_OPP_CUDA)
//     {
//         for (int i = threadIdx.x; i < move_size_OPP_CUDA_8; i+=gpu_threads_per_block_OPP_CUDA) 
//         {
//             atomicAdd(&(ind_arg8[i]), shared_arg8[i]);
//         }
//         __syncthreads();
//     }
// }

// void opp_particle_move__move_fused(
//     opp_set set,      // particles_set
//     opp_arg arg0,     // cell_ef,
//     opp_arg arg1,     // part_pos,
//     opp_arg arg2,     // part_vel,
//     opp_arg arg3,     // part_lc,
//     opp_arg arg4,     // cell_index,
//     opp_arg arg5,     // cell_volume,
//     opp_arg arg6,     // cell_det,
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
    
//     opp_init_particle_move(set, nargs, args);

//     int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
//      opp_mpi_halo_wait_all(nargs, args);
//     opp_profiler->start("MoveToCells1");
//     if (set_size > 0) 
//     {
//         move_stride_OPP_HOST_0 = args[0].dat->set->set_capacity;
//         move_stride_OPP_HOST_1 = args[1].dat->set->set_capacity;
//         move_stride_OPP_HOST_2 = args[2].dat->set->set_capacity;
//         move_stride_OPP_HOST_3 = args[3].dat->set->set_capacity;
//         move_stride_OPP_HOST_6 = args[6].dat->set->set_capacity; 
//         move_stride_OPP_HOST_7 = args[7].size;
//         move_stride_OPP_HOST_8 = args[8].map->from->size;
//         move_size_OPP_HOST_8 = args[8].dat->set->size;

//         cudaMemcpyToSymbol(move_stride_OPP_CUDA_0, &move_stride_OPP_HOST_0, sizeof(int));
//         cudaMemcpyToSymbol(move_stride_OPP_CUDA_1, &move_stride_OPP_HOST_1, sizeof(int));
//         cudaMemcpyToSymbol(move_stride_OPP_CUDA_2, &move_stride_OPP_HOST_2, sizeof(int));
//         cudaMemcpyToSymbol(move_stride_OPP_CUDA_3, &move_stride_OPP_HOST_3, sizeof(int));
//         cudaMemcpyToSymbol(move_stride_OPP_CUDA_6, &move_stride_OPP_HOST_6, sizeof(int));
//         cudaMemcpyToSymbol(move_stride_OPP_CUDA_7, &move_stride_OPP_HOST_7, sizeof(int));
//         cudaMemcpyToSymbol(move_stride_OPP_CUDA_8, &move_stride_OPP_HOST_8, sizeof(int));
//         cudaMemcpyToSymbol(move_size_OPP_CUDA_8, &move_size_OPP_HOST_8, sizeof(int));
//         cudaMemcpyToSymbol(gpu_threads_per_block_OPP_CUDA, &OPP_gpu_threads_per_block, sizeof(int));

//         int start       = 0;
//         int end         = set->size;

//         if (end - start > 0) 
//         {
//             int nthread = OPP_gpu_threads_per_block;
//             int nblocks = (end - start - 1) / nthread + 1;

//             size_t shared_mem = args[8].dat->size * args[8].dat->set->size;
//             if (shared_mem > OPP_gpu_shared_mem_per_block) 
//             {
//                 bool use_shared = false;
//                 cudaMemcpyToSymbol(use_shared_OPP_CUDA, &use_shared, sizeof(bool));
//                 shared_mem = 1;
//             }

//             opp_cuda_all_MoveToCells<<<nblocks, nthread, shared_mem>>>(
//                 (int *)     set->mesh_relation_dat->data_d,
//                 (double *)  args[0].data_d,
//                 (double *)  args[1].data_d,
//                 (double *)  args[2].data_d,
//                 (double *)  args[3].data_d,
//                 (int *)     args[4].data_d,
//                 (double *)  args[5].data_d,
//                 (double *)  args[6].data_d,
//                 (int *)     args[7].data_d,
//                 (double *)  args[8].data_d,
//                 (int *)     args[8].map_data_d,
//                 (double *)  args[9].data_d,
//                 (double *)  args[10].data_d,
//                 (double *)  args[11].data_d,
//                 (int *)     set->particle_remove_count_d,
//                 start, 
//                 end);
//         }
//     }

//     cutilSafeCall(cudaDeviceSynchronize());
//     opp_profiler->end("MoveToCells1");

//     opp_profiler->start("finalize_move");
//     opp_finalize_particle_move(set);
//     opp_profiler->end("finalize_move");

//     opp_set_dirtybit_grouped(nargs, args, Device_GPU);
//     cutilSafeCall(cudaDeviceSynchronize());

//     opp_profiler->end("MoveToCells");
// }

// //*************************************************************************************************







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


int move_stride_OPP_HOST_0 = -1;
int move_stride_OPP_HOST_1 = -1;
int move_stride_OPP_HOST_2 = -1;
int move_stride_OPP_HOST_3 = -1;
int move_stride_OPP_HOST_6 = -1;
int move_stride_OPP_HOST_7 = -1;
int move_stride_OPP_HOST_8 = -1;
int move_size_OPP_HOST_8   = -1;

__constant__ int move_stride_OPP_CUDA_0;
__constant__ int move_stride_OPP_CUDA_1;
__constant__ int move_stride_OPP_CUDA_2;
__constant__ int move_stride_OPP_CUDA_3;
__constant__ int move_stride_OPP_CUDA_6;
__constant__ int move_stride_OPP_CUDA_7;
__constant__ int move_stride_OPP_CUDA_8;
__constant__ int move_size_OPP_CUDA_8;
__constant__ double ONE_OVER_SIX = (1.0 / 6.0);

__constant__ bool use_shared_OPP_CUDA = true;
__constant__ int gpu_threads_per_block_OPP_CUDA;

//user function
//*************************************************************************************************
__device__ void move_all_particles_to_cell__kernel(
    opp_move_var& m,
    const double *cell_ef,
    double *part_pos,
    double *part_vel,
    double *part_lc,
    int* cell_index,
    const double *cell_volume,
    const double *cell_det,
    const int *cell_connectivity,
    double *node_charge_den0,
    double *node_charge_den1,
    double *node_charge_den2,
    double *node_charge_den3
)
{
    if (m.iteration_one)
    {
        double coefficient1 = CONST_charge_cuda / CONST_mass_cuda * (CONST_dt_cuda);

        for (int i = 0; i < DIM; i++)
            part_vel[i * move_stride_OPP_CUDA_2] += 
                (coefficient1 * cell_ef[i * move_stride_OPP_CUDA_0]);           
   
        for (int i = 0; i < DIM; i++)
            part_pos[i * move_stride_OPP_CUDA_1] += 
                part_vel[i * move_stride_OPP_CUDA_2] * (CONST_dt_cuda); // v = u + at
    }

    bool inside = true;
    double coefficient2 = ONE_OVER_SIX / (*cell_volume) ;
    for (int i=0; i<N_PER_C; i++) /*loop over vertices*/
    {
        part_lc[i * move_stride_OPP_CUDA_3] = coefficient2 * (
            cell_det[(i * DET_FIELDS + 0) * move_stride_OPP_CUDA_6] - 
            cell_det[(i * DET_FIELDS + 1) * move_stride_OPP_CUDA_6] * part_pos[0 * move_stride_OPP_CUDA_1] + 
            cell_det[(i * DET_FIELDS + 2) * move_stride_OPP_CUDA_6] * part_pos[1 * move_stride_OPP_CUDA_1] - 
            cell_det[(i * DET_FIELDS + 3) * move_stride_OPP_CUDA_6] * part_pos[2 * move_stride_OPP_CUDA_1]
                );
        
        if (part_lc[i * move_stride_OPP_CUDA_3]<0 || 
            part_lc[i * move_stride_OPP_CUDA_3]>1.0) 
                inside = false;
    }    

    if (inside)
    {
        m.move_status = OPP_MOVE_DONE;

        atomicAdd(node_charge_den0, (part_lc[0 * move_stride_OPP_CUDA_3]));
        atomicAdd(node_charge_den1, (part_lc[1 * move_stride_OPP_CUDA_3]));
        atomicAdd(node_charge_den2, (part_lc[2 * move_stride_OPP_CUDA_3]));
        atomicAdd(node_charge_den3, (part_lc[3 * move_stride_OPP_CUDA_3]));

        return;
    }

    // outside the last known cell, find most negative weight and 
    // use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = part_lc[0 * move_stride_OPP_CUDA_3];
    
    for (int i=1; i<NEIGHB_C; i++)
    {
        if (part_lc[i * move_stride_OPP_CUDA_3] < min_lc) 
        {
            min_lc = part_lc[i * move_stride_OPP_CUDA_3];
            min_i = i;
        }
    }

    if (cell_connectivity[min_i * move_stride_OPP_CUDA_7] >= 0) // is there a neighbor in this direction?
    {
        (*cell_index) = cell_connectivity[min_i * move_stride_OPP_CUDA_7];
        m.move_status = OPP_NEED_MOVE;
    }
    else
    {
        m.move_status = OPP_NEED_REMOVE;
    }
}

// CUDA kernel function
//*************************************************************************************************
__global__ void opp_cuda_all_MoveToCells(
    int *__restrict d_cell_index,
    double *__restrict ind_arg0,
    double *__restrict dir_arg1,
    double *__restrict dir_arg2,
    double *__restrict dir_arg3,
    int *__restrict dir_arg4,
    const double *__restrict ind_arg5,
    const double *__restrict ind_arg6,
    const int *__restrict ind_arg7,
    double *__restrict ind_arg8,
    const int *__restrict opDat8Map,
    double *__restrict ind_arg9,
    double *__restrict ind_arg10,
    double *__restrict ind_arg11,
    int *__restrict particle_remove_count,
    int start,
    int end) 
{
    extern __shared__ double shared_arg8[];
    double* k_arg8 = ind_arg8;

    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (use_shared_OPP_CUDA)
    {
        for (int i = threadIdx.x; i < move_size_OPP_CUDA_8; i+=gpu_threads_per_block_OPP_CUDA) 
            shared_arg8[i] = 0;
        __syncthreads();

        k_arg8 = shared_arg8;
    }

    if (tid + start < end) 
    {
        int n = tid + start;

        opp_move_var m;

        do
        {
            int map0idx = dir_arg4[n];
            // int map0idx = d_cell_index[n]; // TODO : I dont know why this isn't working ??? 
            // dir_arg2 and d_cell_index has same pointer values, but this get stuck!

            const int map1idx = opDat8Map[map0idx + move_stride_OPP_CUDA_8 * 0];
            const int map2idx = opDat8Map[map0idx + move_stride_OPP_CUDA_8 * 1];
            const int map3idx = opDat8Map[map0idx + move_stride_OPP_CUDA_8 * 2];
            const int map4idx = opDat8Map[map0idx + move_stride_OPP_CUDA_8 * 3];

            move_all_particles_to_cell__kernel(
                (m),
                (ind_arg0 + map0idx),   // cell_ef,
                (dir_arg1 + n),         // part_pos,
                (dir_arg2 + n),         // part_vel,
                (dir_arg3 + n),         // part_lc,
                (dir_arg4 + n),         // cell_index,
                (ind_arg5 + map0idx),   // cell_volume,
                (ind_arg6 + map0idx),   // cell_det,
                (ind_arg7 + map0idx),   // cell_connectivity,
                (k_arg8 + map1idx),     // node_charge_den0,
                (k_arg8 + map2idx),     // node_charge_den1,
                (k_arg8 + map3idx),     // node_charge_den2,
                (k_arg8 + map4idx)      // node_charge_den3,  
            );                
            
            m.iteration_one = false;

        } while (m.move_status == (int)OPP_NEED_MOVE);

        if (m.move_status == OPP_NEED_REMOVE) 
        {
            atomicAdd(particle_remove_count, 1);
            dir_arg4[n] = MAX_CELL_INDEX;
        }
    }

    __syncthreads();

    if (use_shared_OPP_CUDA)
    {
        for (int i = threadIdx.x; i < move_size_OPP_CUDA_8; i+=gpu_threads_per_block_OPP_CUDA) 
        {
            atomicAdd(&(ind_arg8[i]), shared_arg8[i]);
        }
        __syncthreads();
    }
}

void opp_particle_move__move_fused(
    opp_set set,      // particles_set
    opp_arg arg0,     // cell_ef,
    opp_arg arg1,     // part_pos,
    opp_arg arg2,     // part_vel,
    opp_arg arg3,     // part_lc,
    opp_arg arg4,     // cell_index,
    opp_arg arg5,     // cell_volume,
    opp_arg arg6,     // cell_det,
    opp_arg arg7,     // cell_connectivity,
    opp_arg arg8,     // node_charge_den0,
    opp_arg arg9,     // node_charge_den1,
    opp_arg arg10,    // node_charge_den2,
    opp_arg arg11     // node_charge_den3,
)
{ 
    
    if (FP_DEBUG) opp_printf("FEMPIC", "opp_particle_move__move_fused set_size %d diff %d", 
        set->size, set->diff);

    opp_profiler->start("MoveToCells");

    int nargs = 12;
    opp_arg args[nargs];

    args[0]  = std::move(arg0);
    args[1]  = std::move(arg1);
    args[2]  = std::move(arg2);
    args[3]  = std::move(arg3);
    args[4]  = std::move(arg4);
    args[5]  = std::move(arg5);
    args[6]  = std::move(arg6);
    args[7]  = std::move(arg7);
    args[8]  = std::move(arg8);
    args[9]  = std::move(arg9);
    args[10] = std::move(arg10);
    args[11] = std::move(arg11);
    
    opp_init_particle_move(set, nargs, args);

    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    opp_profiler->start("MoveToCells1");
    if (set_size > 0) 
    {
        move_stride_OPP_HOST_0 = args[0].dat->set->set_capacity;
        move_stride_OPP_HOST_1 = args[1].dat->set->set_capacity;
        move_stride_OPP_HOST_2 = args[2].dat->set->set_capacity;
        move_stride_OPP_HOST_3 = args[3].dat->set->set_capacity;
        move_stride_OPP_HOST_6 = args[6].dat->set->set_capacity; 
        move_stride_OPP_HOST_7 = args[7].size;
        move_stride_OPP_HOST_8 = args[8].map->from->size;
        move_size_OPP_HOST_8 = args[8].dat->set->size;

        cudaMemcpyToSymbol(move_stride_OPP_CUDA_0, &move_stride_OPP_HOST_0, sizeof(int));
        cudaMemcpyToSymbol(move_stride_OPP_CUDA_1, &move_stride_OPP_HOST_1, sizeof(int));
        cudaMemcpyToSymbol(move_stride_OPP_CUDA_2, &move_stride_OPP_HOST_2, sizeof(int));
        cudaMemcpyToSymbol(move_stride_OPP_CUDA_3, &move_stride_OPP_HOST_3, sizeof(int));
        cudaMemcpyToSymbol(move_stride_OPP_CUDA_6, &move_stride_OPP_HOST_6, sizeof(int));
        cudaMemcpyToSymbol(move_stride_OPP_CUDA_7, &move_stride_OPP_HOST_7, sizeof(int));
        cudaMemcpyToSymbol(move_stride_OPP_CUDA_8, &move_stride_OPP_HOST_8, sizeof(int));
        cudaMemcpyToSymbol(move_size_OPP_CUDA_8, &move_size_OPP_HOST_8, sizeof(int));
        cudaMemcpyToSymbol(gpu_threads_per_block_OPP_CUDA, &OPP_gpu_threads_per_block, sizeof(int));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = OPP_gpu_threads_per_block;
            int nblocks = (end - start - 1) / nthread + 1;

            size_t shared_mem = 0;
            shared_mem += args[8].dat->size * args[8].dat->set->size;
            
            if (shared_mem > OPP_gpu_shared_mem_per_block) 
            {
                bool use_shared = false;
                cudaMemcpyToSymbol(use_shared_OPP_CUDA, &use_shared, sizeof(bool));
                shared_mem = 1;
                
                if (FP_DEBUG) 
                    opp_printf("FEMPIC", "MoveToCells Not using shared mem [avail: %zu, need: %zu]", 
                        OPP_gpu_shared_mem_per_block, shared_mem);
            }
            else
            {
                if (FP_DEBUG) 
                    opp_printf("FEMPIC", "MoveToCells Using shared mem [avail: %zu, need: %zu]", 
                        OPP_gpu_shared_mem_per_block, shared_mem);
            }

            opp_cuda_all_MoveToCells<<<nblocks, nthread, shared_mem>>>(
                (int *)     set->mesh_relation_dat->data_d,
                (double *)  args[0].data_d,
                (double *)  args[1].data_d,
                (double *)  args[2].data_d,
                (double *)  args[3].data_d,
                (int *)     args[4].data_d,
                (double *)  args[5].data_d,
                (double *)  args[6].data_d,
                (int *)     args[7].data_d,
                (double *)  args[8].data_d,
                (int *)     args[8].map_data_d,
                (double *)  args[9].data_d,
                (double *)  args[10].data_d,
                (double *)  args[11].data_d,
                (int *)     set->particle_remove_count_d,
                start, 
                end);
        }
    }

    cutilSafeCall(cudaDeviceSynchronize());
    opp_profiler->end("MoveToCells1");

    opp_profiler->start("finalize_move");
    opp_finalize_particle_move(set);
    opp_profiler->end("finalize_move");

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());

    opp_profiler->end("MoveToCells");
}

//*************************************************************************************************
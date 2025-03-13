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

#include "opp_mpi_defs.h"

/*******************************************************************************
* Inline Utility Functions
*******************************************************************************/
inline int opp_get_uniform_local_size(int global_size, int rank) 
{
    int local_size = global_size / OPP_comm_size;
    int remainder = (int)fmod(global_size, OPP_comm_size);

    if (rank < remainder) 
    {
        local_size = local_size + 1;
    }

    return local_size;
}

//*******************************************************************************
inline int opp_get_uniform_local_size(int global_size) 
{
    return opp_get_uniform_local_size(global_size, OPP_rank);
}

//*************************************************************************************************
inline std::vector<int> get_local_cell_count_array(const int num_cells, const int comm_size)
{
    std::vector<int> local_counts(comm_size);

    int old_value = 0;
    for (int i = 0; i < comm_size; i++) {
        local_counts[i] = opp_get_uniform_local_size(num_cells, i) + old_value;
        old_value = local_counts[i];
    }

    return local_counts;
}

template <typename T> 
inline void opp_uniform_scatter_array(T *g_array, T *l_array, int g_size, int l_size, int elem_size) 
{
    std::vector<int64_t> sendcnts(OPP_comm_size);
    std::vector<int64_t> displs(OPP_comm_size);
    int64_t disp = 0;

    for (int i = 0; i < OPP_comm_size; i++) {
        sendcnts[i] = static_cast<int64_t>(elem_size) * opp_get_uniform_local_size(g_size, i) * sizeof(T);
    }
    for (int i = 0; i < OPP_comm_size; i++) {
        displs[i] = disp;
        disp += sendcnts[i];
    }

    std::vector<MPI_Request> send_req(OPP_comm_size);  // Include all processes
    std::vector<MPI_Request> recv_req(1);  // Only one receive request for the local process

    if (OPP_rank == OPP_ROOT) {
        for (int send_rank = 0; send_rank < OPP_comm_size; send_rank++) {
            MPI_Isend((char*)g_array + displs[send_rank], sendcnts[send_rank], MPI_CHAR, send_rank, 
                        11010, OPP_MPI_WORLD, &send_req[send_rank]);
        }  
    }

    MPI_Irecv((char*)l_array, static_cast<int64_t>(l_size * elem_size * sizeof(T)), 
                    MPI_CHAR, OPP_ROOT, 11010, OPP_MPI_WORLD, &recv_req[0]);

    MPI_Waitall(recv_req.size(), recv_req.data(), MPI_STATUSES_IGNORE);
    if (OPP_rank == OPP_ROOT) 
        MPI_Waitall(send_req.size(), send_req.data(), MPI_STATUSES_IGNORE);

    MPI_Barrier(OPP_MPI_WORLD);
}

/*******************************************************************************
* Functions declared in opp_mpi_partition.cpp
*******************************************************************************/
void opp_partition(std::string lib_name, opp_set prime_set, opp_map prime_map = nullptr, opp_dat data = nullptr);

void print_dat_to_txtfile_mpi(opp_dat dat, const char *file_name);
void opp_mpi_print_dat_to_txtfile(opp_dat dat, const char *file_name);

#if __cplusplus
extern "C" {
#endif
void opp_partition_kway(opp_map primary_map);
void opp_partition_external(opp_set primary_set, opp_dat partvec);
void opp_partition_geom(opp_dat coords);
void opp_partition_destroy();

void partition_all(opp_set primary_set, int my_rank, int comm_size);
void renumber_maps(int my_rank, int comm_size);
void migrate_all(int my_rank, int comm_size);

int frequencyof(int value, int *array, int size);
int find_mode(int *array, int size);
int compare_all_sets(opp_set target_set, opp_set other_sets[], int size);
int *create_exp_list_2(opp_set set, int *temp_list, halo_list h_list, int *part_list, 
                        int size, int comm_size, int my_rank);
void create_imp_list_2(opp_set set, int *temp_list, halo_list h_list, int total_size, int *ranks, 
                        int *sizes, int ranks_size, int comm_size, int my_rank);
int partition_from_set(opp_map map, int my_rank, int comm_size, int **part_range);
int partition_to_set(opp_map map, int my_rank, int comm_size, int **part_range);
#if __cplusplus
}
#endif

/*******************************************************************************
* Functions declared in opp_mpi_particle_comm.cpp
*******************************************************************************/
void opp_part_comm_init();
void opp_part_set_comm_init(opp_set set);
void opp_part_exchange(opp_set set);
bool opp_part_check_all_done(opp_set set);
void opp_part_wait_all(opp_set set);
void opp_part_comm_destroy();

/******************************************************************************** 
  * opp_part_mark_move() will add the opp_part_move_info (particle index and local index of the receiving rank)
  * to a send rank based vector, to be directly used during particle send
*/
inline void opp_part_mark_move(opp_set set, int particle_index, opp_particle_comm_data& comm_data)
{
    // if (OPP_DBG) 
    //     opp_printf("opp_part_mark_move", "commIter[%d] part_id[%d] send_rank[%d] foreign_rank_index[%d]", 
    //         OPP_comm_iteration, particle_index, comm_data.cell_residing_rank, comm_data.local_index);

    std::vector<opp_part_move_info>& vec = opp_part_move_indices[set->index][comm_data.cell_residing_rank];
    vec.emplace_back(particle_index, comm_data.local_index);
}

// returns true, if the current particle needs to be removed from the rank
bool opp_part_checkForGlobalMove(opp_set set, const opp_point& point, const int partIndex, int& cellIdx);

/*******************************************************************************
* Functions declared in opp_mpi_halo_core.cpp
*******************************************************************************/
void __opp_halo_create();
void __opp_halo_destroy();
void __opp_mpi_host_halo_exchange(opp_arg *arg, int exec_flag);
void __opp_mpi_host_halo_wait_all(int nargs, opp_arg *args);

void decl_partition(opp_set set, int *g_index, int *partition);
void get_part_range(int **part_range, int my_rank, int comm_size, MPI_Comm Comm);
int get_partition(int global_index, int *part_range, int *local_index, int comm_size);
int get_global_index(int local_index, int partition, int *part_range, int comm_size);
int get_global_index(int local_index, int partition, int *part_range, int comm_size);
void find_neighbors_set(halo_list List, int *neighbors, int *sizes, int *ranks_size, int my_rank, int comm_size, MPI_Comm Comm);
void create_list(int *list, int *ranks, int *disps, int *sizes, int *ranks_size, int *total, int *temp_list, 
        int size, int comm_size, int my_rank);
void create_export_list(opp_set set, int *temp_list, halo_list h_list, int size, int comm_size, int my_rank);
void create_import_list(opp_set set, int *temp_list, halo_list h_list, int total_size, int *ranks, int *sizes, 
        int ranks_size, int comm_size, int my_rank);

void create_nonexec_import_list(opp_set set, int *temp_list, halo_list h_list, int size, int comm_size, int my_rank);
  
void create_nonexec_export_list(opp_set set, int *temp_list, halo_list h_list, int total_size,
        int *ranks, int *sizes, int ranks_size, int comm_size, int my_rank);
/*******************************************************************************/

void opp_partition_core(std::string lib_name, opp_set prime_set, opp_map prime_map, opp_dat data);

void opp_sanitize_all_maps();
void opp_desanitize_all_maps();

// utility functions
void opp_get_start_end(opp_set set, opp_reset reset, int& start, int& end);
opp_dat opp_mpi_get_data(opp_dat dat);

void opp_process_marked_particles(opp_set set);

//*************************************************************************************************
// ndim : can be 2 or 3
// cell_counts : cell_counts in each direction
// cell_index : cell_index dat which holds global numbering
// cell_colors : local cell_colors dat to colour with most appropriate MPI rank
void __opp_colour_cartesian_mesh(const int ndim, std::vector<int> cell_counts, opp_dat cell_index, 
                            const opp_dat cell_colors, const int cell_ghosts);

void opp_mpi_reduce_double(opp_arg *args, double *data);
void opp_mpi_reduce_int(opp_arg *args, int *data);

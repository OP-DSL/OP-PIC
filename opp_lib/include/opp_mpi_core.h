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

#include <mpi.h>
#include <opp_lib.h>

/** extern variables for halo creation and exchange**/
extern MPI_Comm OPP_MPI_WORLD;

/*******************************************************************************
* MPI halo list data type
*******************************************************************************/
typedef struct 
{
    opp_set set;  // set related to this list
    int size; // number of elements in this list 
    int *ranks; // MPI ranks to be exported to or imported from 
    int ranks_size; // number of MPI neighbors to be exported to or imported from
    int *disps; // displacements for the starting point of each rank's element list 
    int *sizes; // number of elements exported to or imported from each ranks 
    int *list; // the list of all elements
} halo_list_core;

typedef halo_list_core *halo_list;

/*******************************************************************************
* Data structures related to MPI level partitioning
*******************************************************************************/
typedef struct 
{  
    opp_set set; // set to which this partition info blongs to  
    int *g_index;  // global index of each element held in this MPI process 
    int *elem_part;  // partition to which each element belongs  
    int is_partitioned; // indicates if this set is partitioned 1 if partitioned 0 if not
} part_core;

typedef part_core *part;

/*******************************************************************************
* Data structure to hold mpi communications of an opp_dat
*******************************************************************************/
#define NAMESIZE 20
typedef struct 
{ 
    char name[NAMESIZE]; // name of this opp_dat 
    int size; // size of this opp_dat 
    int index; // index of this opp_dat
    int count; // total number of times this opp_dat was halo exported
    int bytes; // total number of bytes halo exported for this opp_dat in this kernel
} opp_dat_mpi_comm_info_core;

typedef opp_dat_mpi_comm_info_core *opp_dat_mpi_comm_info;

/*******************************************************************************
* Buffer struct used in non-blocking mpi halo sends/receives
*******************************************************************************/
typedef struct 
{ 
    char *buf_exec; // buffer holding exec halo to be exported; 
    char *buf_nonexec; // buffer holding nonexec halo to be exported; 
    MPI_Request *s_req; // pointed to hold the MPI_Reqest for sends  
    MPI_Request *r_req; // pointed to hold the MPI_Reqest for receives  
    int s_num_req; // number of send MPI_Reqests in flight at a given time for this opp_dat 
    int r_num_req; // number of receive MPI_Reqests in flight at a given time for this opp_dat
} op_mpi_buffer_core;

typedef op_mpi_buffer_core *op_mpi_buffer;

/*******************************************************************************
* Buffer struct used in non-blocking mpi particle sends/receives
*******************************************************************************/
struct opp_part_neigh_buffers 
{ 
    char *buf_import = nullptr;
    int64_t buf_import_capacity = 0;
    int64_t buf_import_index = 0;       // not used
    char *buf_export = nullptr;
    int64_t buf_export_capacity = 0;    // init to -1
    int64_t buf_export_index = 0;       // init to zero
};

struct opp_part_all_neigh_comm_data 
{ 
    int64_t total_recv = 0;
    std::map<int,opp_part_neigh_buffers> buffers; // rank -> opp_part_neigh_buffers
    std::map<int,int64_t> import_counts; // rank -> count
    std::map<int,int64_t> export_counts; // rank -> count
    std::vector<MPI_Request> send_req;
    std::vector<MPI_Request> recv_req;
    // MPI_Request *r_req;
    // int r_num_req;
    std::vector<int> neighbours;
};

struct opp_particle_comm_data
{ 
    int cell_residing_rank = MAX_CELL_INDEX;
    int local_index = MAX_CELL_INDEX;
};

struct opp_part_move_info
{
    int local_index = MAX_CELL_INDEX;
    int foreign_cell_index = MAX_CELL_INDEX;

    opp_part_move_info(int _local_idx, int _foreign_cell_idx) : 
        local_index(_local_idx), 
        foreign_cell_index(_foreign_cell_idx) {}
};

/** external variables **/
extern int OPP_part_index;
extern part *OPP_part_list;
extern int **orig_part_range;

/** export list on the device **/
extern int **export_exec_list_d;
extern int **export_exec_list_disps_d;
extern int **export_nonexec_list_d;
extern int **export_nonexec_list_disps_d;
// extern int **export_nonexec_list_partial_d;
// extern int **import_nonexec_list_partial_d;
extern int *set_import_buffer_size;
extern int **import_exec_list_disps_d;
extern int **import_nonexec_list_disps_d;

extern halo_list *OPP_export_exec_list; // EEH list
extern halo_list *OPP_import_exec_list; // IEH list

extern halo_list *OPP_import_nonexec_list; // INH list
extern halo_list *OPP_export_nonexec_list; // ENH list

extern int OPP_part_index;
extern part *OPP_part_list;
extern int **orig_part_range;

extern std::map<opp_set, std::map<int, opp_particle_comm_data>> opp_part_comm_neighbour_data; 
extern std::map<int, std::map<int, std::vector<opp_part_move_info>>> opp_part_move_indices;

/*******************************************************************************
* Data Type to hold sliding planes info
*******************************************************************************/

typedef struct 
{
    int index;
    int coupling_group_size;
    int *coupling_proclist;

    int num_ifaces;
    int *iface_list;

    int *nprocs_per_int;
    int **proclist_per_int;
    int **nodelist_send_size;
    int ***nodelist_send;

    int max_data_size;
    char ***send_buf;
    MPI_Request **requests;
    MPI_Status **statuses;

    char *OPP_global_buffer;
    int OPP_global_buffer_size;

    int gbl_num_ifaces;
    int *gbl_iface_list;
    int *nprocs_per_gint;
    int **proclist_per_gint;

    int gbl_offset;
    opp_map cellsToNodes;
    opp_dat coords;
    opp_dat mark;
} op_export_core;

typedef op_export_core *op_export_handle;

typedef struct 
{
    int index;
    int nprocs;
    int *proclist;
    int gbl_offset;
    opp_dat coords;
    opp_dat mark;
    int max_dat_size;
    int num_my_ifaces;
    int *iface_list;
    int *nprocs_per_int;
    int **proclist_per_int;
    int *node_size_per_int;
    int **nodelist_per_int;
    char ***recv_buf;
    int *recv2int;
    int *recv2proc;
    MPI_Request *requests;
    MPI_Status *statuses;
    double *interp_dist;
} op_import_core;

typedef op_import_core *op_import_handle;

namespace opp {

    //*******************************************************************************
    class Comm {

    public:    
        Comm(MPI_Comm comm_parent);;
        ~Comm();

    public:
        /// Parent (i.e. global for the simulation) MPI communicator.
        MPI_Comm comm_parent;
        /// Communicator between one rank on each shared memory region.
        MPI_Comm comm_inter;
        /// Communicator between the ranks in a shared memory region.
        MPI_Comm comm_intra;
        /// MPI rank in the parent communicator.
        int rank_parent = -1;
        /// MPI rank in the inter shared memory region communicator.
        int rank_inter = -1;
        /// MPI rank within the shared memory region communicator.
        int rank_intra = -1;
        /// Size of the parent communicator.
        int size_parent = -1;
        /// Size of the inter shared memory communicator.
        int size_inter = -1;
        /// Size of the intra shared memory communicator.
        int size_intra = -1;
    };

    //*******************************************************************************
    class dh_particle_packer {

    public:
        dh_particle_packer(std::map<int, std::map<int, std::vector<opp_part_move_info>>>* part_move_data);
        ~dh_particle_packer();
        void pack(opp_set set);
        char* get_buffer(const opp_set set, const int send_rank);
        void unpack(opp_set set, const std::map<int, std::vector<char>>& particleRecvBuffers,
                            int64_t totalParticlesToRecv, const std::vector<int64_t>& recvRankPartCounts);

    private: 
        std::map<int, std::map<int, std::vector<opp_part_move_info>>>* part_move_data = nullptr;
        std::map<int, std::map<int, std::vector<char>>> buffers;
    };


    //*******************************************************************************
    class GlobalParticleMover {

    public:
        GlobalParticleMover(MPI_Comm comm);
        ~GlobalParticleMover();
        
        void markParticleToMove(opp_set set, int partIndex, int rankToBeMoved, int finalGlobalCellIndex);
        void initGlobalMove();
        void communicateParticleSendRecvRankCounts();
        void communicate(opp_set set);
        int64_t finalize(opp_set set);
    
    private:
        // this translate to std::map<particle_set_index, std::map<send_rank, std::vector<opp_part_move_info>>>
        std::map<int, std::map<int, std::vector<opp_part_move_info>>> dh_part_move_data;
        MPI_Comm comm;

        int *recv_win_data;
        MPI_Win recv_win;
        MPI_Request mpi_request;

        std::vector<MPI_Request> h_send_requests;
        std::vector<MPI_Request> h_recv_requests;
        std::vector<MPI_Request> h_send_data_requests;
        std::vector<MPI_Request> h_recv_data_requests;
        std::vector<MPI_Status> h_recv_status;

        std::vector<int> h_send_ranks;
        std::vector<int64_t> h_send_rank_npart;
        std::vector<int> h_recv_ranks;
        std::vector<int64_t> h_recv_rank_npart;

        std::map<int, std::vector<char>> particleRecvBuffers;

        int64_t totalParticlesToSend = 0;
        int64_t totalParticlesToRecv = 0;
        int numRemoteRecvRanks = 0;
        int numRemoteSendRanks = 0;

        std::unique_ptr<dh_particle_packer> packer;
    };
};



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

// //*******************************************************************************
// template <typename T> 
// inline void opp_uniform_scatter_array(T *g_array, T *l_array, int g_size, int l_size, int elem_size) 
// {
//     int64_t *sendcnts = (int64_t *)opp_host_malloc(OPP_comm_size * sizeof(int64_t));
//     int64_t *displs = (int64_t *)opp_host_malloc(OPP_comm_size * sizeof(int64_t));
//     int64_t disp = 0;

//     for (int i = 0; i < OPP_comm_size; i++) 
//     {
//         sendcnts[i] = (int64_t)elem_size * opp_get_uniform_local_size(g_size, i) * sizeof(T);
//         //printf("RANK %d %d\t| sendcount %d | %d\n", OPP_rank, g_size, sendcnts[i], opp_get_uniform_local_size(g_size, i));
//     }
//     for (int i = 0; i < OPP_comm_size; i++) 
//     {
//         displs[i] = disp;
//         disp = disp + sendcnts[i];
//         //printf("RANK %d %d\t| displs %d\n", OPP_rank, g_size, displs[i]);
//     }

//     MPI_Scatterv((char*)g_array, sendcnts, displs, MPI_CHAR, 
//         (char*)l_array, (int64_t)(l_size * elem_size * sizeof(T)), MPI_CHAR, OPP_ROOT, MPI_COMM_WORLD);

//     opp_host_free(sendcnts);
//     opp_host_free(displs);
// }

template <typename T> 
inline void opp_uniform_scatter_array(T *g_array, T *l_array, int g_size, int l_size, int elem_size) 
{
    std::vector<int64_t> sendcnts(OPP_comm_size);
    std::vector<int64_t> displs(OPP_comm_size);
    int64_t disp = 0;

    for (int i = 0; i < OPP_comm_size; i++) 
    {
        sendcnts[i] = static_cast<int64_t>(elem_size) * opp_get_uniform_local_size(g_size, i) * sizeof(T);
    }
    for (int i = 0; i < OPP_comm_size; i++) 
    {
        displs[i] = disp;
        disp += sendcnts[i];
    }

    std::vector<MPI_Request> send_req(OPP_comm_size);  // Include all processes
    std::vector<MPI_Request> recv_req(1);  // Only one receive request for the local process

    if (OPP_rank == OPP_ROOT) 
    {
        for (int send_rank = 0; send_rank < OPP_comm_size; send_rank++) 
        {
            MPI_Isend((char*)g_array + displs[send_rank], sendcnts[send_rank], MPI_CHAR, send_rank, 
                        11010, MPI_COMM_WORLD, &send_req[send_rank]);
        }  
    }

    MPI_Irecv((char*)l_array, static_cast<int64_t>(l_size * elem_size * sizeof(T)), 
                    MPI_CHAR, OPP_ROOT, 11010, MPI_COMM_WORLD, &recv_req[0]);
    MPI_Waitall(recv_req.size(), recv_req.data(), MPI_STATUSES_IGNORE);
    if (OPP_rank == OPP_ROOT) 
        MPI_Waitall(send_req.size(), send_req.data(), MPI_STATUSES_IGNORE);

    MPI_Barrier(MPI_COMM_WORLD);
}

//*******************************************************************************

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
int *create_exp_list_2(opp_set set, int *temp_list, halo_list h_list, int *part_list, int size, int comm_size, int my_rank);
void create_imp_list_2(opp_set set, int *temp_list, halo_list h_list, int total_size, int *ranks, int *sizes, int ranks_size, int comm_size, int my_rank);
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

extern std::vector<int> opp_move_part_indices;
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
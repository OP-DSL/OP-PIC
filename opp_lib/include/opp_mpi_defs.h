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

extern std::vector<int> opp_move_part_indices;

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
        dh_particle_packer(std::map<int, std::vector<OPP_INT>>& local_part_indices, 
                            std::map<int, std::vector<OPP_INT>>& foreign_cell_indices);
        virtual ~dh_particle_packer();

        char* get_buffer(const opp_set set, const int send_rank);
        
        virtual void pack(opp_set set) = 0;        
        virtual void unpack(opp_set set, const std::map<int, std::vector<char>>& particleRecvBuffers,
                            int64_t totalParticlesToRecv) = 0;
    
    protected:
        std::map<int, std::map<int, std::vector<char>>> buffers;
        std::map<int, std::vector<OPP_INT>>& local_part_ids;
        std::map<int, std::vector<OPP_INT>>& foreign_cell_ids;
    };

    class dh_particle_packer_cpu : public dh_particle_packer { //********************

    public:
        dh_particle_packer_cpu(std::map<int, std::vector<OPP_INT>>& local_part_indices, 
                                std::map<int, std::vector<OPP_INT>>& foreign_cell_indices);
        ~dh_particle_packer_cpu();

        void pack(opp_set set) override;
        void unpack(opp_set set, const std::map<int, std::vector<char>>& particleRecvBuffers,
                            int64_t totalParticlesToRecv) override;
    };

    class dh_particle_packer_gpu : public dh_particle_packer { //********************

    public:
        dh_particle_packer_gpu(std::map<int, std::vector<OPP_INT>>& local_part_indices, 
                                std::map<int, std::vector<OPP_INT>>& foreign_cell_indices);
        ~dh_particle_packer_gpu();

        void pack(opp_set set) override;
        void unpack(opp_set set, const std::map<int, std::vector<char>>& particleRecvBuffers,
                            int64_t totalParticlesToRecv) override;
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
        // std::map<int, std::map<int, std::vector<opp_part_move_info>>> dh_part_move_data;

        std::map<int, std::vector<OPP_INT>> dh_local_part_indices;
        std::map<int, std::vector<OPP_INT>> dh_foreign_cell_indices;

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

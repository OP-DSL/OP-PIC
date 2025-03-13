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

#include "opp_cuda.h"

constexpr bool verbose_profile = false;

void opp_part_pack_device(opp_set set);
void opp_part_unpack_device(opp_set set);
void particle_hole_fill_device(opp_set set);
void opp_part_pack_and_exchange_device_direct(opp_set set);
void opp_part_unpack_device_direct(opp_set set);

//*******************************************************************************

OPP_INT *OPP_move_count_d = nullptr;
OPP_INT OPP_move_count_h = 0;

thrust::device_vector<OPP_INT> OPP_move_particle_indices_dv;
OPP_INT *OPP_move_particle_indices_d = nullptr;

thrust::device_vector<OPP_INT> OPP_move_cell_indices_dv;
OPP_INT *OPP_move_cell_indices_d = nullptr;

thrust::device_vector<OPP_INT> OPP_remove_particle_indices_dv;
OPP_INT *OPP_remove_particle_indices_d = nullptr;

struct CopyMaxCellIndexFunctor {
    OPP_INT* B;
    CopyMaxCellIndexFunctor(OPP_INT* B) : B(B) {}

    __host__ __device__
    void operator()(OPP_INT index) const {
        B[index] = MAX_CELL_INDEX;
    }
};

//*******************************************************************************
void opp_init_particle_move(opp_set set, int nargs, opp_arg *args)
{ 
    opp_init_particle_move_core(set);

    opp_mem::copy_host_to_dev<OPP_INT>(set->particle_remove_count_d, &(set->particle_remove_count), 1);

    // const size_t buffer_alloc_size = (size_t)(set->set_capacity / 2);
    const size_t buffer_alloc_size = (size_t)(set->set_capacity);
    if (buffer_alloc_size > OPP_move_particle_indices_dv.size()) {   

        OPP_move_particle_indices_dv.resize(buffer_alloc_size);
        OPP_move_particle_indices_d = opp_get_dev_raw_ptr(OPP_move_particle_indices_dv);

        OPP_move_cell_indices_dv.resize(buffer_alloc_size);
        OPP_move_cell_indices_d = opp_get_dev_raw_ptr(OPP_move_cell_indices_dv);

        OPP_remove_particle_indices_dv.resize(buffer_alloc_size);
        OPP_remove_particle_indices_d = opp_get_dev_raw_ptr(OPP_remove_particle_indices_dv);
    }

    if (OPP_move_count_d == nullptr)
        OPP_move_count_d = opp_mem::dev_malloc<OPP_INT>(1);
    OPP_move_count_h = 0;
    opp_mem::copy_host_to_dev<OPP_INT>(OPP_move_count_d, &OPP_move_count_h, 1);

    if (OPP_comm_iteration == 0) {
        OPP_iter_start = 0;
        OPP_iter_end   = set->size;
        OPP_part_comm_count_per_iter = 0;    
    }
    else {
        // need to change the arg pointers since communication could change due to realloc
        for (int i = 0; i < nargs; i++) {
            if (args[i].argtype == OPP_ARG_DAT && args[i].dat->set->is_particle) {
                args[i].data = args[i].dat->data;
                args[i].data_d = args[i].dat->data_d;
            }
        }
    }

    if (OPP_DBG) opp_printf("opp_init_particle_move", "comm_iter=%d start=%d end=%d", 
                    OPP_comm_iteration, OPP_iter_start, OPP_iter_end);

    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 
    OPP_mesh_relation_data_d = ((OPP_INT *)set->mesh_relation_dat->data_d); 
}

//*******************************************************************************
bool opp_finalize_particle_move(opp_set set)
{ 
    OPP_DEVICE_SYNCHRONIZE();

    opp_profiler->start("Mv_Finalize");

    // this is the exchange particle count
    OPP_move_count_h = -1;
    opp_mem::copy_dev_to_host<int>(&OPP_move_count_h, OPP_move_count_d, 1);

    OPP_part_comm_count_per_iter += OPP_move_count_h;

    // remove count is the addition of removed particles and the exchange count
    set->particle_remove_count = -1;
    opp_mem::copy_dev_to_host<int>(&(set->particle_remove_count), set->particle_remove_count_d, 1);

    if (OPP_DBG)
        opp_printf("opp_finalize_particle_move", "set [%s][%d] remove_count [%d] move count [%d]", 
            set->name, set->size, set->particle_remove_count, OPP_move_count_h);

#ifdef USE_MPI
    // At this stage, particles of device is clean
    if (OPP_gpu_direct) {
        opp_part_pack_and_exchange_device_direct(set);
    }
    else {
        // download only the required particles to send and pack them in rank based mpi buffers
        opp_part_pack_device(set);

        opp_part_exchange(set);  // send the counts and send the particle data 
    }
#endif

    opp_profiler->start("Mv_fill");
    if (set->particle_remove_count > 0) {
    
        set->size -= set->particle_remove_count;

        if (OPP_fill_type == OPP_HoleFill_All || 
            (OPP_fill_type == OPP_Sort_Periodic || OPP_fill_type == OPP_Shuffle_Periodic) && 
                (OPP_main_loop_iter % OPP_fill_period != 0 || OPP_comm_iteration != 0)) {
        
            if (OPP_DBG) 
                opp_printf("opp_finalize_particle_move", "hole fill set [%s]", set->name);
            
            opp_profiler->start("Mv_holefill");
            particle_hole_fill_device(set);
            opp_profiler->end("Mv_holefill");
        }
        else if (OPP_fill_type == OPP_Sort_All || OPP_fill_type == OPP_Sort_Periodic) {
        
            if (OPP_DBG)
                opp_printf("opp_finalize_particle_move", "sort set [%s]", set->name);
            
            opp_profiler->start("Mv_sort");
            opp_particle_sort(set);
            opp_profiler->end("Mv_sort");
        }
        else if (OPP_fill_type == OPP_Shuffle_All || OPP_fill_type == OPP_Shuffle_Periodic) {
        
            if (OPP_DBG) 
                opp_printf("opp_finalize_particle_move", "shuffle set [%s]", set->name);
            
            opp_profiler->start("Mv_shuffle");
            particle_sort_device(set, true); // true will shuffle the particles
            opp_profiler->end("Mv_shuffle");
        }
    }
    opp_profiler->end("Mv_fill");

#ifdef USE_MPI
    if (opp_part_check_all_done(set)) {
    
        if (OPP_max_comm_iteration < OPP_comm_iteration)
            OPP_max_comm_iteration = OPP_comm_iteration;

        OPP_comm_iteration = 0; // reset for the next par loop
        
        OPP_DEVICE_SYNCHRONIZE();
        opp_profiler->end("Mv_Finalize");
        return false; // all mpi ranks do not have anything to communicate to any rank
    }

    // particles are expected to be received by atleast one MPI rank, wait till all are communicated
    opp_part_wait_all(set);

    if (OPP_DBG)
        opp_printf("opp_finalize_particle_move", "set [%s] size prior unpack %d", set->name, set->size);
    
    OPP_DEVICE_SYNCHRONIZE();

    // if current MPI rank received particles, increase the particle count if required and 
    // unpack the communicated particles to separate dats
    if (OPP_gpu_direct) {
        opp_part_unpack_device_direct(set);  
    }
    else {
        opp_part_unpack_device(set);    
    }

    OPP_iter_start = set->size - set->diff;
    OPP_iter_end   = set->size;  

    OPP_comm_iteration++;  

    opp_profiler->end("Mv_Finalize");
    return true; // need to run another communication iteration (particle move loop)
#else
    opp_profiler->end("Mv_Finalize");
    return false;
#endif
}

// Cannot use multiple packs before sending them, if opp_part_pack() is called multiple times 
// with PACK_SOA, the communication data may get currupted
//*******************************************************************************
void opp_part_pack_device(opp_set set)
{
    if (OPP_DBG) opp_printf("opp_part_pack_device", "start");

#ifdef USE_MPI
    opp_profiler->start("Mv_Pack");

    if (OPP_move_count_h <= 0) {
        opp_profiler->end("Mv_Pack");
        return;
    }

    // Uncomment below if comparing cuda_mpi with sycl_mpi or hip_mpi
    if (debugger) {
        thrust::sort_by_key(
            OPP_move_particle_indices_dv.begin(), 
            OPP_move_particle_indices_dv.begin() + OPP_move_count_h, 
            OPP_move_cell_indices_dv.begin());
    }

    if (verbose_profile) opp_profiler->start("Mv_Pack1");
    thrust::host_vector<int> send_part_cell_idx_hv(OPP_move_count_h);
    thrust::copy(OPP_move_cell_indices_dv.begin(), 
        OPP_move_cell_indices_dv.begin() + OPP_move_count_h, 
        send_part_cell_idx_hv.begin());
    if (verbose_profile) opp_profiler->end("Mv_Pack1");

    std::map<int, std::vector<char>> move_dat_data_map;

    // copy and download the particles to send (to host arrays)
    if (verbose_profile) opp_profiler->start("Mv_Pack3");
    for (auto& dat : *(set->particle_dats)) {

        const size_t bytes_to_copy = (OPP_move_count_h * dat->size);   
        auto& move_dat_data = move_dat_data_map[dat->index];
        move_dat_data.resize(bytes_to_copy);

        if (strcmp(dat->type, "double") == 0) {
            copy_according_to_index(dat->thrust_real, dat->thrust_real_sort, 
                OPP_move_particle_indices_dv, dat->set->set_capacity, 
                OPP_move_count_h, OPP_move_count_h, dat->dim);

            opp_mem::copy_dev_to_host<char>(move_dat_data.data(), 
                (char*)opp_get_dev_raw_ptr<OPP_REAL>(*(dat->thrust_real_sort)), bytes_to_copy);
        }
        else if (strcmp(dat->type, "int") == 0) {
            copy_according_to_index(dat->thrust_int, dat->thrust_int_sort, 
                OPP_move_particle_indices_dv, dat->set->set_capacity, 
                OPP_move_count_h, OPP_move_count_h, dat->dim);
            
            opp_mem::copy_dev_to_host<char>(move_dat_data.data(), 
                (char*)opp_get_dev_raw_ptr<OPP_INT>(*(dat->thrust_int_sort)), bytes_to_copy);
        }
    }      
    if (verbose_profile) opp_profiler->end("Mv_Pack3");

    if (verbose_profile) opp_profiler->start("Mv_Pack2");
    // enrich the particles to communicate with the correct external cell index and mpi rank
    std::map<int, opp_particle_comm_data>& set_part_com_data = opp_part_comm_neighbour_data[set];
    for (int index = 0; index < OPP_move_count_h; index++) {

        int map0idx = send_part_cell_idx_hv[index];

        auto it = set_part_com_data.find(map0idx);
        if (it == set_part_com_data.end()) {
            opp_printf("opp_part_pack_device", 
                "Error: cell %d cannot be found in opp_part_comm_neighbour_data map [%d/%d]", 
                map0idx, index, OPP_move_count_h);
            continue; // unlikely, need opp_abort() instead!
        }

        opp_part_mark_move(set, index, it->second); // it->second is the local cell index in foreign rank
    }
    if (verbose_profile) opp_profiler->end("Mv_Pack2");

    opp_part_all_neigh_comm_data* send_buffers = (opp_part_all_neigh_comm_data*)set->mpi_part_buffers;

    // increase the sizes of MPI buffers
    if (verbose_profile) opp_profiler->start("Mv_Pack4");
    for (auto& move_indices_per_rank : opp_part_move_indices[set->index]) {

        const int send_rank = move_indices_per_rank.first;
        std::vector<opp_part_move_info>& move_indices_vector = move_indices_per_rank.second;

        opp_part_neigh_buffers& send_rank_buf = send_buffers->buffers[send_rank];
        int64_t req_buf_size = (move_indices_vector.size() * (int64_t)set->particle_size);

        // resize the export buffer if required
        if (send_rank_buf.buf_export_index + req_buf_size >= send_rank_buf.buf_export_capacity) {
        
            if (send_rank_buf.buf_export == nullptr) {
            
                send_rank_buf.buf_export_capacity = OPP_mpi_part_alloc_mult * req_buf_size;
                send_rank_buf.buf_export_index = 0;
                send_rank_buf.buf_export = opp_mem::host_malloc<char>(send_rank_buf.buf_export_capacity);

                // opp_printf("opp_part_pack", "alloc buf_export cap %d", send_rank_buf.buf_export_capacity);
            }
            else {
                // Assume that there are some particles left already, increase capacity beyond buf_export_index
                send_rank_buf.buf_export_capacity = send_rank_buf.buf_export_index + 
                                                            OPP_mpi_part_alloc_mult * req_buf_size;
                opp_mem::host_realloc<char>(send_rank_buf.buf_export, send_rank_buf.buf_export_capacity);
                
                // opp_printf("opp_part_pack", "realloc buf_export cap %d", send_rank_buf.buf_export_capacity);
            }        
        }
    }
    if (verbose_profile) opp_profiler->end("Mv_Pack4");

    // iterate over all the ranks and pack to mpi buffers using SOA
    if (verbose_profile) opp_profiler->start("Mv_Pack5");
    for (auto& move_indices_per_rank : opp_part_move_indices[set->index]) {
    
        const int send_rank = move_indices_per_rank.first;
        std::vector<opp_part_move_info>& move_indices_vec = move_indices_per_rank.second;

        opp_part_neigh_buffers& send_rank_buf = send_buffers->buffers[send_rank];

        int64_t disp = 0;
        for (auto& dat : *(set->particle_dats)) {
            auto& move_dat_data = move_dat_data_map[dat->index];

            const int64_t dat_size = (int64_t)dat->size;
            const int64_t element_size = (int64_t)(dat->size / dat->dim);

            if (dat->is_cell_index) {
                for (const auto& move_info : move_indices_vec) {

                    // Need to copy the cell index of the foreign rank, to correctly unpack in the foreign rank
                    memcpy(&(send_rank_buf.buf_export[send_rank_buf.buf_export_index + disp]), 
                        &move_info.foreign_cell_index, dat->size);
                 
                    disp += dat_size;
                }
            }
            else {
                for (int d = 0; d < dat->dim; d++) {
                    for (const auto& move_info : move_indices_vec) {

                        // copy the multi dimensional dat value to the send buffer
                        memcpy(&(send_rank_buf.buf_export[send_rank_buf.buf_export_index + disp]), 
                            &(move_dat_data[(d * OPP_move_count_h + move_info.local_index) * element_size]), 
                            element_size);
                        
                        disp += element_size;
                    }
                }                
            }
        }

        send_rank_buf.buf_export_index = (int64_t)(set->particle_size * move_indices_vec.size()); // Unused
        (send_buffers->export_counts)[send_rank] = (int64_t)move_indices_vec.size();

        move_indices_vec.clear();
    }
    if (verbose_profile) opp_profiler->end("Mv_Pack5");

    // This particle is already packed, hence needs to be removed from the current rank
    if (verbose_profile) opp_profiler->start("Mv_Pack6");
    CopyMaxCellIndexFunctor copyMaxCellIndexFunctor((OPP_INT*)set->mesh_relation_dat->data_d);
    thrust::for_each(OPP_move_particle_indices_dv.begin(), 
        OPP_move_particle_indices_dv.begin() + OPP_move_count_h, 
        copyMaxCellIndexFunctor);
    if (verbose_profile) opp_profiler->end("Mv_Pack6");

    opp_profiler->end("Mv_Pack");
#endif

    if (OPP_DBG) opp_printf("opp_part_pack_device", "end");
}

//*******************************************************************************
void opp_part_unpack_device(opp_set set)
{
    if (OPP_DBG) opp_printf("opp_part_unpack_device", "set [%s]", set->name);

#ifdef USE_MPI
    opp_profiler->start("Mv_Unpack");

    std::vector<opp_dat>& particle_dats = *(set->particle_dats);
    int64_t new_part_count = 0;

    opp_part_all_neigh_comm_data* recv_buffers = (opp_part_all_neigh_comm_data*)set->mpi_part_buffers;
    std::vector<int>& neighbours = recv_buffers->neighbours;

    // count the number of particles to be received from all ranks
    for (size_t i = 0; i < neighbours.size(); i++) {
        int neighbour_rank = neighbours[i];
        new_part_count += (recv_buffers->import_counts)[neighbour_rank];
    }

    if (new_part_count > 0) {
        const int64_t recv_part_start_idx = (int64_t)(set->size);

        opp_increase_particle_count(set, (int)new_part_count);

        // create a continuous memory in host to copy to device
        std::vector<std::vector<char>> temp_data_vec;
        for (size_t d = 0; d < particle_dats.size(); d++) {
            temp_data_vec.push_back(std::vector<char>(particle_dats[d]->size * new_part_count));
        }

        int current_recv_count = 0;
        for (int i = 0; i < (int)neighbours.size(); i++) {
        
            const int recv_rank = neighbours[i];

            opp_part_neigh_buffers& recv_rank_buf = recv_buffers->buffers[recv_rank];

            const int64_t recv_count = recv_buffers->import_counts[recv_rank];
            int64_t displacement = 0;

            for (size_t dat_idx = 0; dat_idx < particle_dats.size(); dat_idx++) {
            
                opp_dat dat = particle_dats[dat_idx];
                std::vector<char>& temp_data = temp_data_vec[dat_idx];

                const int element_size = dat->size / dat->dim;

                for (int i = 0; i < dat->dim; i++) {
                    memcpy(&(temp_data[element_size * (i * new_part_count + current_recv_count)]), 
                        &(recv_rank_buf.buf_import[displacement]), element_size * recv_count);

                    displacement += element_size * recv_count; 
                }                
            }

            current_recv_count += recv_count;
        }

        // copy to device
        for (size_t dat_idx = 0; dat_idx < particle_dats.size(); dat_idx++) {
        
            opp_dat dat = particle_dats[dat_idx];
            const std::vector<char>& temp_data = temp_data_vec[dat_idx];

            const size_t bytes_to_copy_per_dim = new_part_count * dat->size / dat->dim;
            const int element_size = dat->size / dat->dim;

            for (int64_t d = 0; d < dat->dim; d++) {
                 
                const size_t data_d_offset = (recv_part_start_idx + d * set->set_capacity) * element_size;
                const size_t data_h_offset = d * new_part_count * element_size;

                char* data_d = dat->data_d + data_d_offset;
                opp_mem::copy_host_to_dev<char>(data_d, &(temp_data[data_h_offset]), 
                                                bytes_to_copy_per_dim);   
            }
        }
        OPP_DEVICE_SYNCHRONIZE();
    }

    opp_profiler->end("Mv_Unpack");
#endif

    if (OPP_DBG) opp_printf("opp_part_unpack_device", "END");    
}






// Below is for MPI GPU direct comms
//*******************************************************************************
#define MPI_COUNT_EXCHANGE 0
#define MPI_TAG_PART_EX 1


__global__ void copy_intX(const int* in_dat_d, int* out_dat_d, const int* indices,
                    const int in_stride, const int out_stride, const int dim, const int size) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;

    if (tid < size) 
    {
        const int idx = indices[tid];
        for (int d = 0; d < dim; d++)
        {
            out_dat_d[tid + d * out_stride] = in_dat_d[idx + d * in_stride];

            // printf("copting index=%d value %d [d=%d] to out_index=%d\n", idx, in_dat_d[idx + d * in_stride], d, tid + d * out_stride);
        }
    }
}

__global__ void copy_doubleX(const double* in_dat_d, double* out_dat_d, const int* indices, 
                    const int in_stride, const int out_stride, const int dim, const int size) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;

    if (tid < size) 
    {
        const int idx = indices[tid];
        for (int d = 0; d < dim; d++)
        {
            out_dat_d[tid + d * out_stride] = in_dat_d[idx + d * in_stride];
        }
    }
}

__global__ void copy_doubleY(const double* in_dat_d, double* out_dat_d, int in_stride, 
                                    int out_stride, int dim, int size) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;

    if (tid < size) 
    {
        for (int d = 0; d < dim; d++)
        {
            out_dat_d[tid + d * out_stride] = in_dat_d[tid + d * in_stride];
        }
    }
}
__global__ void copy_intY(const int* in_dat_d, int* out_dat_d, int in_stride, 
                                    int out_stride, int dim, int size, int x) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;

    if (tid < size) 
    {
        for (int d = 0; d < dim; d++)
        {
            // printf("unpacking %d index=%d value %d [d=%d] to out_index=%d\n", x, tid, in_dat_d[tid + d * in_stride], d, tid + d * out_stride);
            out_dat_d[tid + d * out_stride] = in_dat_d[tid + d * in_stride];
        }
    }
}

__global__ void setArrayToMaxCID(int* array, int* indices, int size) 
{
    int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;
    if (tid < size) {
        int idx = indices[tid];
        array[idx] = MAX_CELL_INDEX;
    }
}

std::map<int, thrust::host_vector<OPP_INT>> particle_indices_hv;    // particle ids to send, arrange according to rank
std::map<int, thrust::host_vector<OPP_INT>> cell_indices_hv;        // cellid in the foreign rank, arrange according to rank
std::map<int, thrust::device_vector<OPP_INT>> particle_indices_dv;
std::map<int, thrust::device_vector<char>> send_data;
std::map<int, thrust::device_vector<char>> recv_data;
const int threads = 64;
const double opp_comm_buff_resize_multiple = 1.5;

// Cannot use multiple packs before sending them, if opp_part_pack() is called multiple times with PACK_SOA, 
// the communication data may get currupted
//*******************************************************************************
void opp_part_pack_and_exchange_device_direct(opp_set set)
{
    if (OPP_DBG) opp_printf("opp_part_pack_and_exchange_device_direct", "OPP_move_count_h %d", OPP_move_count_h);

#ifdef USE_MPI
    opp_profiler->start("Mv_PackExDir");

    std::map<int, cudaStream_t> streams;

    cudaStreamCreate(&(streams[-1]));
    cudaStreamCreate(&(streams[-2]));

    thrust::host_vector<int> tmp_cell_indices_hv(OPP_move_count_h);
    cudaMemcpyAsync(thrust::raw_pointer_cast(tmp_cell_indices_hv.data()), 
            thrust::raw_pointer_cast(OPP_move_cell_indices_dv.data()),
            OPP_move_count_h * sizeof(int), cudaMemcpyDeviceToHost, streams[-1]);
    
    thrust::host_vector<int> tmp_particle_indices_hv(OPP_move_count_h);
    cudaMemcpyAsync(thrust::raw_pointer_cast(tmp_particle_indices_hv.data()), 
            thrust::raw_pointer_cast(OPP_move_particle_indices_dv.data()),
            OPP_move_count_h * sizeof(int), cudaMemcpyDeviceToHost, streams[-2]);

    opp_part_all_neigh_comm_data* mpi_buffers = (opp_part_all_neigh_comm_data*)set->mpi_part_buffers;
    const std::vector<int>& neighbours = mpi_buffers->neighbours;
    const int neighbour_count = neighbours.size();
    mpi_buffers->total_recv = 0;
    for (auto it = mpi_buffers->import_counts.begin(); it != mpi_buffers->import_counts.end(); it++)
        it->second = 0;
    
    for (auto it = particle_indices_hv.begin(); it != particle_indices_hv.end(); it++) it->second.clear();
    for (auto it = cell_indices_hv.begin(); it != cell_indices_hv.end(); it++) it->second.clear();
    for (auto it = particle_indices_dv.begin(); it != particle_indices_dv.end(); it++) it->second.clear();
    for (auto it = send_data.begin(); it != send_data.end(); it++) it->second.clear();
    for (auto it = recv_data.begin(); it != recv_data.end(); it++) it->second.clear();

    mpi_buffers->recv_req.clear();
    mpi_buffers->send_req.clear();
    std::vector<MPI_Request> send_req_count(neighbour_count);
    std::vector<MPI_Request> recv_req_count(neighbour_count);
    double total_send_size = 0.0;

    std::map<int, opp_particle_comm_data>& set_part_com_data = opp_part_comm_neighbour_data[set];

    cudaStreamSynchronize(streams[-1]);
    cudaStreamSynchronize(streams[-2]);

    // enrich and arrange the particles to communicate with the correct external cell index and mpi rank
    for (int index = 0; index < OPP_move_count_h; index++)
    {
        auto it = set_part_com_data.find(tmp_cell_indices_hv[index]);
        if (it == set_part_com_data.end()) 
        {
            opp_printf("opp_part_pack_and_exchange_device_direct", 
                "Error: cell %d cannot be found in opp_part_comm_neighbour_data map", tmp_cell_indices_hv[index]);
            continue; // unlikely, need opp_abort() instead!
        }

        const auto& comm_data = it->second;
        particle_indices_hv[comm_data.cell_residing_rank].push_back(tmp_particle_indices_hv[index]);
        cell_indices_hv[comm_data.cell_residing_rank].push_back(comm_data.local_index); // convert cid to local cid of recv rank
    }
    
    // copy particle_indices_dv to device asynchronously 
    for (const auto& x : particle_indices_hv)
    {
        const int rank = x.first;
        cudaStreamCreate(&(streams[rank])); 
        const size_t tmp_cpy_size = x.second.size();

        if (tmp_cpy_size > particle_indices_dv[rank].capacity()) 
            particle_indices_dv[rank].reserve(tmp_cpy_size * opp_comm_buff_resize_multiple);
        particle_indices_dv[rank].resize(tmp_cpy_size);

        cudaMemcpyAsync(thrust::raw_pointer_cast(particle_indices_dv[rank].data()), 
            x.second.data(), (tmp_cpy_size * sizeof(int)), cudaMemcpyHostToDevice, streams[rank]);
        
        mpi_buffers->export_counts[rank] = tmp_cpy_size;
    }

    // send/receive send_counts to/from all immediate neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        const int64_t& send_count = mpi_buffers->export_counts[neighbours[i]];
        MPI_Isend((void*)&send_count, 1, MPI_INT64_T, neighbours[i], MPI_COUNT_EXCHANGE, 
            OPP_MPI_WORLD, &(send_req_count[i]));

        const int64_t& recv_count = mpi_buffers->import_counts[neighbours[i]];
        MPI_Irecv((void*)&recv_count, 1, MPI_INT64_T, neighbours[i], MPI_COUNT_EXCHANGE, 
            OPP_MPI_WORLD, &(recv_req_count[i]));
    }

    // pack the send data to device memory arranged according to rank asynchronously
    for (auto& x : particle_indices_dv)
    {
        const int send_rank = x.first;
        const int64_t particle_count = x.second.size();
        const int64_t send_bytes = particle_count * set->particle_size;
        auto& send_data_dv = send_data[send_rank];

        if (send_bytes > send_data_dv.capacity()) 
            send_data_dv.reserve(send_bytes * opp_comm_buff_resize_multiple);
        send_data_dv.resize(send_bytes);

        char* send_buff = (char*)thrust::raw_pointer_cast(send_data_dv.data());
    
        int* particle_indices = (int*)thrust::raw_pointer_cast(particle_indices_dv[send_rank].data());
        const int nblocks = (particle_count - 1) / threads + 1;
        int64_t offset = 0;

        // thrust::host_vector<int> h_vec = particle_indices_dv[send_rank];
        // std::string log = "";
        // for (int i = 0; i < h_vec.size(); ++i) log += std::to_string(h_vec[i]) + " ";
        // opp_printf("TO_SEND", "%s", log.c_str());

        cudaStreamSynchronize(streams[send_rank]);

        for (auto& dat : *(set->particle_dats)) 
        {
            const int64_t dat_bytes_to_copy = (particle_count * dat->size);

            if (dat->is_cell_index)
            {
                // cell indices relative to the receiving rank is copied here
                cudaMemcpyAsync((send_buff + offset), 
                    (char*)cell_indices_hv[send_rank].data(), dat_bytes_to_copy, 
                    cudaMemcpyHostToDevice, streams[send_rank]);
            }
            else if (strcmp(dat->type, "double") == 0)
            {
                copy_doubleX<<<nblocks,threads,0,streams[send_rank]>>> (
                    (OPP_REAL*)dat->data_d,
                    (OPP_REAL*)(send_buff + offset),
                    particle_indices,
                    set->set_capacity, particle_count,
                    dat->dim, particle_count);
            }
            else if (strcmp(dat->type, "int") == 0)
            {
                copy_intX<<<nblocks,threads,0,streams[send_rank]>>> (
                    (OPP_INT*)dat->data_d,
                    (OPP_INT*)(send_buff + offset),
                    particle_indices,
                    set->set_capacity, particle_count,
                    dat->dim, particle_count);
            }
            else
            {
                opp_printf("", "Error: %s type unimplemented in opp_part_pack_and_exchange_device_direct", dat->type);
                opp_abort("datatype not implemented in opp_part_pack_and_exchange_device_direct");
            }

            offset += dat_bytes_to_copy;
        }
    }

    // since move particles ids are extracted already, mark cell index as MAX_CELL_ID to remove from current rank
    const int nblocks = (OPP_move_count_h - 1) / threads + 1;
    setArrayToMaxCID<<<nblocks,threads>>> (
        (OPP_INT*)set->mesh_relation_dat->data_d,
        (OPP_INT*)thrust::raw_pointer_cast(OPP_move_particle_indices_dv.data()),
        OPP_move_count_h);

    // send the particle data only to immediate neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        const int send_rank = neighbours[i];
        const int64_t send_count = mpi_buffers->export_counts[send_rank];

        if (send_count <= 0) {
            if (OPP_DBG) opp_printf("opp_part_pack_and_exchange_device_direct", "nothing to send to rank %d", send_rank);
            continue;
        }
        else {   
            char* send_buff = (char*)thrust::raw_pointer_cast(send_data[send_rank].data());
            if (OPP_DBG) 
                opp_printf("opp_part_pack_and_exchange_device_direct", "sending %lld particle/s (size: %lld) to rank %d | %p", 
                send_count, (int64_t)(send_count*set->particle_size), send_rank, send_buff);
        }

        MPI_Request req;
        const int64_t send_size = set->particle_size * send_count;

        cudaStreamSynchronize(streams[send_rank]); // wait till cuda aync copy is done

        char* send_buff = (char*)thrust::raw_pointer_cast(send_data[send_rank].data());
        MPI_Isend(send_buff, send_size, MPI_CHAR, send_rank, MPI_TAG_PART_EX, OPP_MPI_WORLD, &req);
        mpi_buffers->send_req.push_back(req);

        total_send_size += (send_size * 1.0f);
    }

    // wait for the counts to receive only from neighbours
    MPI_Waitall(neighbour_count, &recv_req_count[0], MPI_STATUSES_IGNORE);

    // create/resize data structures and receive particle data from neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        const int recv_rank = neighbours[i];
        const int64_t recv_bytes = (int64_t)set->particle_size * mpi_buffers->import_counts[recv_rank];
        mpi_buffers->total_recv += mpi_buffers->import_counts[recv_rank];

        if (recv_bytes <= 0)
        {
            if (OPP_DBG) 
                opp_printf("opp_part_pack_and_exchange_device_direct", "nothing to receive from rank %d", recv_rank);
            continue;
        }

        auto& recv_data_dv = recv_data[recv_rank];
        if (recv_bytes > recv_data_dv.capacity()) 
            recv_data_dv.reserve(recv_bytes * opp_comm_buff_resize_multiple);
        recv_data_dv.resize(recv_bytes);
        
        MPI_Request req;
        MPI_Irecv((char*)thrust::raw_pointer_cast(recv_data_dv.data()), recv_bytes, MPI_CHAR, 
            recv_rank, MPI_TAG_PART_EX, OPP_MPI_WORLD, &req);
        mpi_buffers->recv_req.push_back(req);
    }

    // reset the export counts for another iteration
    for (auto it = mpi_buffers->export_counts.begin(); it != mpi_buffers->export_counts.end(); it++)
    {
        it->second = 0; // make the export count to zero for the next iteration
        mpi_buffers->buffers[it->first].buf_export_index = 0; // make export indices to zero for next iteration
    }

    // for (const auto& x : streams) cudaStreamDestroy(x.second);
    OPP_DEVICE_SYNCHRONIZE();

    opp_profiler->end("Mv_PackExDir");
#endif

    if (OPP_DBG) opp_printf("opp_part_pack_device_direct", "end");
}


void opp_part_unpack_device_direct(opp_set set)
{
    if (OPP_DBG) opp_printf("opp_part_unpack_device_direct", "set [%s] size %d", set->name, set->size);

#ifdef USE_MPI
    opp_profiler->start("Mv_UnpackDir");

    opp_part_all_neigh_comm_data* recv_buffers = (opp_part_all_neigh_comm_data*)set->mpi_part_buffers;
    const auto& neighbours = recv_buffers->neighbours;
    int64_t num_new_particles = 0;
    std::map<int,int64_t> particle_start;

    // count the number of particles to be received from all ranks
    for (size_t i = 0; i < neighbours.size(); i++)
    {
        const int rank = neighbours[i];

        num_new_particles += (recv_buffers->import_counts)[rank];
        if (i == 0) 
            particle_start[i] = set->size;
        else 
            particle_start[i] = particle_start[i-1] + (recv_buffers->import_counts)[neighbours[i-1]];
    }

    if (num_new_particles > 0)
    {
        opp_increase_particle_count(set, (int)num_new_particles);

        for (int i = 0; i < (int)neighbours.size(); i++)
        {
            const int recv_rank = neighbours[i];
            const int recv_count = (int)recv_buffers->import_counts[recv_rank];

            if (recv_count <= 0) continue;

            char* recv_buff = (char*)thrust::raw_pointer_cast(recv_data[recv_rank].data()); 
            const int nblocks = (recv_count - 1) / threads + 1;
            int64_t offset = 0;

            // opp_printf("opp_part_unpack_device_direct", "recv count %d || to dat starting from %lld", 
            //     recv_count, particle_start[i]);

            for (auto& dat : *(set->particle_dats)) 
            {
                const int64_t dat_bytes = (recv_count * dat->size);
                const int64_t dat_per_dim_size = dat->size / dat->dim;

                if (strcmp(dat->type, "double") == 0)
                {
                    copy_doubleY<<<nblocks,threads>>> (
                        (OPP_REAL*)(recv_buff + offset),
                        (OPP_REAL*)(dat->data_d + particle_start[i] * dat_per_dim_size),
                        recv_count, set->set_capacity, dat->dim, recv_count);
                }
                else if (strcmp(dat->type, "int") == 0)
                {
                    int x = 0;

                    if (strcmp(dat->name, "p_index") == 0) x = 111;

                    copy_intY<<<nblocks,threads>>> (
                        (OPP_INT*)(recv_buff + offset),
                        (OPP_INT*)(dat->data_d + particle_start[i] * dat_per_dim_size),
                        recv_count, set->set_capacity, dat->dim, recv_count, x);
                }
                else
                {
                    opp_printf("", "Error: %s type unimplemented in opp_part_unpack_device_direct", dat->type);
                    opp_abort("datatype not implemented in opp_part_unpack_device_direct");
                }

                offset += dat_bytes;
            }         
        }

        OPP_DEVICE_SYNCHRONIZE();
    }

    opp_profiler->end("Mv_UnpackDir");
#endif

    if (OPP_DBG) opp_printf("opp_part_unpack_device_direct", "END");    
}

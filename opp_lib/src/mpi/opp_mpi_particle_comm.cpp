
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

#include <opp_mpi.h>

#define MPI_MH_COUNT_EXCHANGE 0
#define MPI_MH_TAG_PART_EX 1

#define PACK_AOS  false // PACK_AOS == false, is performing better

// this translate to std::map<opp_set, std::map<local_cell_index, opp_particle_comm_data>>
std::map<opp_set, std::map<int, opp_particle_comm_data>> opp_part_comm_neighbour_data; 

// this translate to std::map<particle_set, std::map<send_rank, std::vector<opp_part_move_info>>>
std::map<int, std::map<int, std::vector<opp_part_move_info>>> opp_part_move_indices;

// This vector will get enriched during move loop, if a particle needs to be communicated
std::vector<int> opp_move_part_indices; 

/*******************************************************************************
 * opp_part_pack() will get the opp_part_move_indices and pack the particle data into send rank based buffers
 * Used in Multi-hop
*/
void opp_part_pack(opp_set set)
{
    if (OPP_DBG) 
        opp_printf("opp_part_pack", "set [%s]", set->name);

    opp_profiler->start("Mv_Pack");

    opp_part_all_neigh_comm_data* part_send_data = (opp_part_all_neigh_comm_data*)set->mpi_part_buffers;

    // iterate over all the ranks to sent particles to
    for (auto& move_indices_per_rank : opp_part_move_indices[set->index])
    {
        const int send_rank = move_indices_per_rank.first;
        std::vector<opp_part_move_info>& move_indices_vector = move_indices_per_rank.second;

        opp_part_neigh_buffers& send_rank_buffer = part_send_data->buffers[send_rank];
        const int64_t required_buffer_size = (move_indices_vector.size() * (int64_t)set->particle_size);

        // resize the export buffer if required
        if (send_rank_buffer.buf_export_index + required_buffer_size >= send_rank_buffer.buf_export_capacity)
        {
            if (send_rank_buffer.buf_export == nullptr)
            {
                send_rank_buffer.buf_export_capacity  = required_buffer_size;
                send_rank_buffer.buf_export_index     = 0;
                send_rank_buffer.buf_export           = (char *)opp_host_malloc(send_rank_buffer.buf_export_capacity);

                // opp_printf("opp_part_pack", "alloc buf_export capacity %" PRId64, send_rank_buffer.buf_export_capacity);
            }
            else 
            {
                // Assume that there are some particles left already, increase capacity beyond buf_export_index
                send_rank_buffer.buf_export_capacity  = send_rank_buffer.buf_export_index + required_buffer_size;
                send_rank_buffer.buf_export           = (char *)opp_host_realloc(send_rank_buffer.buf_export, 
                                                                        send_rank_buffer.buf_export_capacity);
                
                // opp_printf("opp_part_pack", "realloc buf_export capacity %" PRId64, send_rank_buffer.buf_export_capacity);
            }        
        }

#if PACK_AOS

        // iterate over all the particles to be sent to this send_rank
        for (auto& part : move_indices_vector)
        {
            int displacement = 0;

            // pack the particle data from dats into the export buffer 
            for (auto& dat : *(set->particle_dats))
            {
                if (dat->is_cell_index)
                {
                    // we need to copy the cell index of the foreign rank, to correctly unpack in the foreign rank
                    memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), 
                        &part.foreign_cell_index, dat->size);
                }
                else
                {
                    // copy the dat value to the send buffer
                    memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), 
                        &(dat->data[part.local_index * dat->size]), dat->size);
                }

                displacement += dat->size;
            } 

            send_rank_buffer.buf_export_index += set->particle_size;
            (part_send_data->export_counts)[send_rank] += 1;
        }

#else // PACK_SOA

        // Cannot use multiple packs before sending them, if opp_part_pack() is called multiple times with PACK_SOA, 
        // the communication data will get currupted
        int64_t displacement = 0;
        for (auto& dat : *(set->particle_dats))
        {
            int64_t dat_size = (int64_t)dat->size;

            if (dat->is_cell_index)
            {
                for (const auto& part : move_indices_vector)
                {
                    // we need to copy the cell index of the foreign rank, to correctly unpack in the foreign rank
                    memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), 
                        &part.foreign_cell_index, dat->size);
                 
                    displacement += dat_size;
                }
            }
            else
            {
                for (const auto& part : move_indices_vector)
                {
                    // copy the dat value to the send buffer
                    memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), 
                        &(dat->data[part.local_index * dat->size]), dat->size);
                    
                    displacement += dat_size;
                }                
            }
        }

        send_rank_buffer.buf_export_index = (int64_t)(set->particle_size * move_indices_vector.size()); // Not used
        (part_send_data->export_counts)[send_rank] = (int64_t)move_indices_vector.size();

#endif

        move_indices_vector.clear();
    }

    opp_profiler->end("Mv_Pack");

    if (OPP_DBG) 
        opp_printf("opp_part_pack", "set [%s] END", set->name);
}

/*******************************************************************************
 * opp_part_unpack() will npack the received particle data and adds to the end of the set data
 * Used in Multi-hop
*/
void opp_part_unpack(opp_set set)
{
    if (OPP_DBG) opp_printf("opp_part_unpack", "set [%s]", set->name);

    opp_profiler->start("Mv_Unpack");

    int64_t num_particles = 0;

    opp_part_all_neigh_comm_data* part_recv_data = (opp_part_all_neigh_comm_data*)set->mpi_part_buffers;
    std::vector<int>& neighbours = part_recv_data->neighbours;

    // count the number of particles to be received from all ranks
    for (size_t i = 0; i < neighbours.size(); i++)
    {
        const int neighbour_rank = neighbours[i];
        num_particles += (part_recv_data->import_counts)[neighbour_rank];
    }

    if (num_particles > 0)
    {
        if (!opp_increase_particle_count_core(set, (int)num_particles)) // TODO : change this to int64_t
        {
            opp_printf("opp_part_unpack", "Error: Failed to increase particle count of particle set [%s]", set->name);
            opp_abort("opp_part_unpack");
        }

        int64_t new_part_index = (int64_t)(set->size - set->diff);

        for (int i = 0; i < (int)neighbours.size(); i++)
        {
            const int recv_rank = neighbours[i];

            opp_part_neigh_buffers& recv_rank_buffer = part_recv_data->buffers[recv_rank];
            const int64_t receive_count = part_recv_data->import_counts[recv_rank];

#if PACK_AOS
            std::vector<opp_dat>& particle_dats = *(set->particle_dats);
            int64_t particle_size = set->particle_size;

            // unpack the received buffer from rank 'recv_rank' in to particle dats
            for (int part = 0; part < receive_count; part++)
            {
                const char* part_buffer = &(recv_rank_buffer.buf_import[particle_size * part]);
                int displacement = 0;

                for (int i = 0; i < (int)particle_dats.size(); i++)
                {
                    opp_dat& dat = particle_dats[i];
                    memcpy(dat->data + new_part_index * dat->size, part_buffer + displacement, dat->size);
                    displacement += dat->size;     
                }
                new_part_index++;
            }

#else // PACK_SOA

            int64_t displacement = 0;
            for (auto& dat : *(set->particle_dats))
            {
                memcpy(dat->data + new_part_index * dat->size, recv_rank_buffer.buf_import + displacement, 
                    dat->size * receive_count);              
                displacement += dat->size * receive_count;            
            }

            new_part_index += receive_count;
#endif
        }
    }

    opp_profiler->end("Mv_Unpack");

    if (OPP_DBG) opp_printf("opp_part_unpack", "END");
}

/*******************************************************************************
 * During particle move loop, only the particle index is added (to opp_move_part_indices), 
 * opp_process_marked_particles() will get the required info for sending
 * Used in Multi-hop
*/
void opp_process_marked_particles(opp_set set) 
{
    for (int particle_index : opp_move_part_indices)
    {
        std::map<int, opp_particle_comm_data>& set_part_com_data = opp_part_comm_neighbour_data[set];
        const int map0idx = OPP_mesh_relation_data[particle_index];
        
        // Removed from the current rank
        OPP_mesh_relation_data[particle_index] = MAX_CELL_INDEX;
        set->particle_remove_count += 1;

        auto it = set_part_com_data.find(map0idx);
        if (it == set_part_com_data.end())
        {
            opp_printf("opp_part_check_status", 
                "Error: cell %d cannot be found in opp_part_comm_neighbour_data map particle_index=%d", 
                    map0idx, particle_index);
            opp_abort("opp_part_check_status Error: cell cannot be found in opp_part_comm_neighbour_data map");
        }

        opp_part_mark_move(set, particle_index, it->second);
    }
}

std::vector<MPI_Request> send_count_reqs;
std::vector<MPI_Request> recv_count_reqs;

/*******************************************************************************
 * This routine will exchange particles within neighbours, will return before receiving the particles
 * should follow with opp_part_wait_all() if particles are sent/received
 * Used in Multi-hop
*/
void opp_part_exchange(opp_set set)
{
    if (OPP_DBG) opp_printf("opp_part_exchange", "set[%s] particle size[%d]", set->name, set->particle_size);

    opp_profiler->start("Mv_Exchange");

    opp_part_all_neigh_comm_data* mpi_part_data = (opp_part_all_neigh_comm_data*)set->mpi_part_buffers;

    const std::vector<int>& neighbours = mpi_part_data->neighbours;
    const int neighbour_count = neighbours.size();
    mpi_part_data->total_recv = 0;

    for (auto it = mpi_part_data->import_counts.begin(); it != mpi_part_data->import_counts.end(); it++)
        it->second = 0;

    mpi_part_data->recv_req.clear();
    mpi_part_data->send_req.clear();

    send_count_reqs.clear();
    recv_count_reqs.clear();

    send_count_reqs.resize(neighbour_count);
    recv_count_reqs.resize(neighbour_count);

    opp_profiler->startMpiComm("", opp::OPP_Particle);

    const int count_ex_tag = OPP_main_loop_iter * 10 + OPP_comm_iteration;
    const int part_ex_tag = 100000 + OPP_main_loop_iter * 10 + OPP_comm_iteration;

    // send/receive send_counts to/from only to neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        const int neighbour_rank = neighbours[i];

        const int64_t& send_count = mpi_part_data->export_counts[neighbour_rank];
        MPI_Isend((void*)&send_count, 1, MPI_INT64_T, neighbour_rank, count_ex_tag, // MPI_MH_COUNT_EXCHANGE, 
                    OPP_MPI_WORLD, &(send_count_reqs[i]));

        const int64_t& recv_count = mpi_part_data->import_counts[neighbour_rank];
        MPI_Irecv((void*)&recv_count, 1, MPI_INT64_T, neighbour_rank, count_ex_tag, // MPI_MH_COUNT_EXCHANGE, 
                    OPP_MPI_WORLD, &(recv_count_reqs[i]));
        
        // opp_printf("Exchange", "%d->%d Send %" PRId64 " particles", OPP_rank, neighbour_rank, send_count);
    }

    double total_send_bytes = 0.0;

    // send the particle data only to neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        const int neighbour_rank = neighbours[i];
        const int64_t send_bytes = set->particle_size * mpi_part_data->export_counts[neighbour_rank];

        if (send_bytes <= 0)
        {
            if (OPP_DBG) opp_printf("opp_part_exchange", "nothing to send to rank %d", neighbour_rank);
            continue;
        }
        else
        {   
            if (OPP_DBG) opp_printf("opp_part_exchange", "sending %lld particle/s (size: %lld) to rank %d", 
                (send_bytes/ set->particle_size), send_bytes, neighbour_rank);
            total_send_bytes += (send_bytes * 1.0f);
        }

        char* send_rank_buffer = mpi_part_data->buffers[neighbour_rank].buf_export;

        MPI_Request req;
        MPI_Isend(send_rank_buffer, send_bytes, MPI_CHAR, neighbour_rank, part_ex_tag, OPP_MPI_WORLD, &req);
        mpi_part_data->send_req.push_back(req); 
    }

    opp_profiler->addTransferSize("", opp::OPP_Particle, total_send_bytes, (size_t)(total_send_bytes / set->particle_size));
    opp_profiler->end("Mv_Exchange");

    // wait for the counts to receive only from neighbours
    const std::string profName = std::string("Mv_WaitExCnt") + std::to_string(OPP_comm_iteration);
    opp_profiler->start(profName);
    std::vector<MPI_Status> send_req_statuses(neighbour_count);
    std::vector<MPI_Status> recv_req_statuses(neighbour_count);
    MPI_Waitall(neighbour_count, &send_count_reqs[0], send_req_statuses.data());
    MPI_Waitall(neighbour_count, &recv_count_reqs[0], recv_req_statuses.data()); 
    bool error = false;
    for (size_t i = 0; i < send_req_statuses.size(); ++i) {
        if (send_req_statuses[i].MPI_ERROR != MPI_SUCCESS) {
            opp_printf("opp_part_exchange Recv", "Error in send request ", i, send_req_statuses[i].MPI_ERROR);
            error = true;
        }
    }
    for (size_t i = 0; i < recv_req_statuses.size(); ++i) {
        if (recv_req_statuses[i].MPI_ERROR != MPI_SUCCESS) {
            opp_printf("opp_part_exchange Recv", "Error in recv request ", i, recv_req_statuses[i].MPI_ERROR);
            error = true;
        }
    }
    if (error)
        opp_abort("Error in opp_part_exchange");
    opp_profiler->end(profName);

    opp_profiler->start("Mv_Exchange");

    // create/resize data structures and receive particle data from neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        const int neighbour_rank = neighbours[i];
        const int64_t recv_bytes = (int64_t)set->particle_size * mpi_part_data->import_counts[neighbour_rank];
        mpi_part_data->total_recv += mpi_part_data->import_counts[neighbour_rank];

        if (recv_bytes <= 0)
        {
            if (OPP_DBG) opp_printf("opp_part_exchange", "nothing to receive from rank %d", neighbour_rank);
            continue;
        }

        opp_part_neigh_buffers& recv_buffer = mpi_part_data->buffers[neighbour_rank];

        // resize the receive buffer is required
        if (recv_bytes >= recv_buffer.buf_import_capacity)
        {
            if (recv_buffer.buf_import == nullptr)
            {
                recv_buffer.buf_import_capacity  = recv_bytes;           
                recv_buffer.buf_import = (char *)opp_host_malloc(recv_buffer.buf_import_capacity);          
            }
            else
            {
                recv_buffer.buf_import_capacity += recv_bytes;
                recv_buffer.buf_import = (char *)opp_host_realloc(recv_buffer.buf_import, recv_buffer.buf_import_capacity);
            }
        }

        // opp_printf("opp_part_exchange", "neighbour_rank: %d, recv_bytes: %d, total_recv: %d, buff_import_capacity: %d", 
        //     neighbour_rank, recv_bytes, mpi_part_data->total_recv, recv_buffer.buf_import_capacity);

        MPI_Request req;
        MPI_Irecv(recv_buffer.buf_import, recv_bytes, MPI_CHAR, neighbour_rank, part_ex_tag, OPP_MPI_WORLD, &req);
        mpi_part_data->recv_req.push_back(req);
    }

    // reset the export counts for another iteration
    for (auto it = mpi_part_data->export_counts.begin(); it != mpi_part_data->export_counts.end(); it++)
    {
        it->second = 0; // make the export count to zero for the next iteration
        mpi_part_data->buffers[it->first].buf_export_index = 0; // make export index of that rank to zero for next iteration
    }

    opp_profiler->endMpiComm("", opp::OPP_Particle);

    opp_profiler->end("Mv_Exchange");

    if (OPP_DBG) opp_printf("opp_part_exchange", "END");
}

/*******************************************************************************
 * opp_part_wait_all() will wait all the Multip-hop comminications are done, this is a blocking call around neighbours
*/
void opp_part_wait_all(opp_set set)
{
    if (OPP_DBG) opp_printf("opp_part_wait_all", "START");

    opp_profiler->start("Mv_WaitAll");

    opp_part_all_neigh_comm_data* mpi_part_data = (opp_part_all_neigh_comm_data*)set->mpi_part_buffers;

    std::vector<MPI_Request>& send_req = mpi_part_data->send_req;
    std::vector<MPI_Request>& recv_req = mpi_part_data->recv_req;

    opp_profiler->startMpiComm("", opp::OPP_Particle);

    // wait till all the particles from all the ranks are received
    std::string profName = std::string("Mv_WaitExRecv") + std::to_string(OPP_comm_iteration);
    opp_profiler->start(profName);

    std::vector<MPI_Status> send_req_statuses(send_req.size());
    std::vector<MPI_Status> recv_req_statuses(recv_req.size());

    MPI_Waitall(send_req.size(), &(send_req[0]), send_req_statuses.data());
    MPI_Waitall(recv_req.size(), &(recv_req[0]), recv_req_statuses.data());
    opp_profiler->end(profName);

    bool error = false;
    for (size_t i = 0; i < send_req_statuses.size(); ++i) {
        if (send_req_statuses[i].MPI_ERROR != MPI_SUCCESS) {
            opp_printf("opp_part_wait_all", "Error in send request ", i, send_req_statuses[i].MPI_ERROR);
            error = true;
        }
    }
    for (size_t i = 0; i < recv_req_statuses.size(); ++i) {
        if (recv_req_statuses[i].MPI_ERROR != MPI_SUCCESS) {
            opp_printf("opp_part_wait_all", "Error in recv request ", i, recv_req_statuses[i].MPI_ERROR);
            error = true;
        }
    }
    if (error)
        opp_abort("Error in opp_part_wait_all");

    opp_profiler->endMpiComm("", opp::OPP_Particle); // started at opp_part_exchange()

    send_req.clear();
    recv_req.clear();

    opp_profiler->end("Mv_WaitAll");

    if (OPP_DBG) opp_printf("opp_part_wait_all", "END");
}

/*******************************************************************************
 * opp_part_check_all_done() will wait till all the MPI ranks are done with their work, this is synchronous
 * returns true if all MPI ranks are done with their work, i.e. no particles sent or received, else false
*/
bool opp_part_check_all_done(opp_set set)
{
    if (OPP_DBG) opp_printf("opp_part_check_all_done", "START");

    const std::string profName = std::string("Mv_WaitDone") + std::to_string(OPP_comm_iteration);
    opp_profiler->start(profName);

    const int64_t total_recv_parts = ((opp_part_all_neigh_comm_data*)set->mpi_part_buffers)->total_recv;
    char imported_parts = 1;

    if (total_recv_parts == 0)
    {
        imported_parts = 0;
        set->diff = 0;
    }

    opp_profiler->startMpiComm("", opp::OPP_Particle);

    // gather from all MPI ranks to see whether atleast one rank needs to communicate to another
    std::vector<char> buffer_recv(OPP_comm_size, 0);
    MPI_Allgather(&imported_parts, 1, MPI_C_BOOL, buffer_recv.data(), 1, MPI_C_BOOL, OPP_MPI_WORLD);

    opp_profiler->endMpiComm("", opp::OPP_Particle);

    bool bool_ret = false;
    std::string log = "";
    for (int i = 0; i < OPP_comm_size; i++)
    {
        bool_ret = (bool_ret || buffer_recv[i]);

        if (OPP_DBG) 
            log += std::string(" ") + std::to_string(i) + (buffer_recv[i] ? "A" : "D");
    }

    if (OPP_DBG) 
        opp_printf("opp_part_check_all_done", "iteration %d recv %lld - %s -%s", OPP_comm_iteration,
            total_recv_parts, (bool_ret ? "ITER AGAIN" : "ALL DONE"), log.c_str());

    opp_profiler->end(profName);

    return !bool_ret;
}

/*******************************************************************************
 * Initialize the particle communication buffers and halo mappings, done at opp_init()
*/
void opp_part_comm_init()
{
    if (OPP_DBG) opp_printf("opp_part_comm_init", "START");

    for (opp_set& set : opp_sets)
    {
        if (!set->is_particle) continue;

        opp_part_set_comm_init(set);
    }

    opp_profiler->reg("Mv_Pack");
    opp_profiler->reg("Mv_Exchange");
    opp_profiler->reg("Mv_Unpack");
    opp_profiler->reg("Mv_WaitAll");

    std::string profName = "";
    for (int i = 0; i < 5; i++) {
        profName = std::string("Mv_WaitDone") + std::to_string(i);
        opp_profiler->reg(profName);
        profName = std::string("Mv_WaitExCnt") + std::to_string(i);
        opp_profiler->reg(profName);
        profName = std::string("Mv_WaitExRecv") + std::to_string(i);
        opp_profiler->reg(profName);
    }

    if (OPP_DBG) opp_printf("opp_part_comm_init", "END");
}

/*******************************************************************************
 * Initialize the particle communication buffers per particle set and halo mappings, done at opp_init()
*/
void opp_part_set_comm_init(opp_set set)
{
    // TODO : can use the same mappings for all particle sets with the same cells set, instead of communicating again
    const halo_list exp_exec_list = OPP_export_exec_list[set->cells_set->index];
    const halo_list imp_exec_list = OPP_import_exec_list[set->cells_set->index];

    const size_t exp_exec_rank_count = (size_t)(exp_exec_list->ranks_size);
    const size_t imp_exec_rank_count = (size_t)(imp_exec_list->ranks_size);

    if (OPP_DBG) opp_printf("opp_part_set_comm_init", "set: %s exp_exec_ranks:%zu imp_exec_ranks:%zu", 
                    set->name, exp_exec_rank_count, imp_exec_rank_count);

    std::vector<MPI_Request> send_reqs;
    std::vector<MPI_Request> recv_reqs;    
    std::vector<std::vector<int>> recv_buffers(imp_exec_rank_count);

    // send the local index of the export elements of the set to all neighbours
    for (size_t i = 0; i < exp_exec_rank_count; i++) 
    {  
        const int neighbour_rank = exp_exec_list->ranks[i];
        const int* send_buffer = &(exp_exec_list->list[exp_exec_list->disps[i]]);
        const int send_size = exp_exec_list->sizes[i];
    
        MPI_Request req;
        MPI_Isend(send_buffer, send_size, MPI_INT, neighbour_rank, OPP_rank, OPP_MPI_WORLD, &req);
        send_reqs.push_back(req);  

        //print the per rank send buffers
        if (OPP_DBG)
        {
            std::string log = "";
            for (int k = 0; k < send_size; k++)
                log += std::to_string(send_buffer[k]) + " ";        
            opp_printf("opp_part_set_comm_init", "%s SEND neighbour_rank %d send_size %d -> %s", 
                set->cells_set->name, neighbour_rank, send_size, log.c_str());
        } 
    }

    // receive the foreign ranks local index of the import elements of the set from all neighbours
    for (size_t i = 0; i < imp_exec_rank_count; i++) 
    {
        // opp_printf("opp_part_set_comm_init", "imp_exec_rank_count: %zu %zu", imp_exec_rank_count, i);
        const int neighbour_rank = imp_exec_list->ranks[i];
        const int recv_size = imp_exec_list->sizes[i];
        
        std::vector<int>& recv_buffer = recv_buffers[i];
        recv_buffer.resize(recv_size);
        
        MPI_Request req;
        MPI_Irecv(&recv_buffer[0], recv_size, MPI_INT, neighbour_rank, neighbour_rank, OPP_MPI_WORLD, &req);
        recv_reqs.push_back(req);  
    }

    std::vector<MPI_Status> recv_reqs_statuses(recv_reqs.size());
    MPI_Waitall(recv_reqs.size(), &recv_reqs[0], recv_reqs_statuses.data());
    bool error = false;
    for (size_t i = 0; i < recv_reqs_statuses.size(); ++i) {
        if (recv_reqs_statuses[i].MPI_ERROR != MPI_SUCCESS) {
            opp_printf("opp_part_set_comm_init", "Error in recv request ", i, recv_reqs_statuses[i].MPI_ERROR);
            error = true;
        }
    }
    if (error)
        opp_abort("Error in opp_part_set_comm_init");

    // print the per rank received buffers
    if (OPP_DBG)
    {
        for (size_t i = 0; i < recv_buffers.size(); i++) 
        {
            // what I have (mappings) in import exec buffers, mappings before renumbering from that rank
            const int* imp_buffer = &(imp_exec_list->list[imp_exec_list->disps[i]]); 

            std::string log = "";
            for (size_t k = 0; k < recv_buffers[i].size(); k++)
                log += std::to_string((recv_buffers[i])[k]) + "|" + std::to_string(imp_buffer[k]) + " ";  

            opp_printf("opp_part_set_comm_init", "%s RECEIVE neighbour_rank %d recv_size %d (new|old) -> %s", 
                set->cells_set->name, imp_exec_list->ranks[i], recv_buffers[i].size(), log.c_str());
        }
    }

    if (true) // TODO : this might break existing OP2 functionality, check for issues
    { 
        for (size_t i = 0; i < recv_buffers.size(); i++) 
        {
            int* imp_buffer = &(imp_exec_list->list[imp_exec_list->disps[i]]); 
            
            for (int k = 0; k < (int)recv_buffers[i].size(); k++)
                imp_buffer[k] = (recv_buffers[i])[k];  
        }
    }

    // create mappings of neighbour ranks cell information for easy access during particle communication
    for (size_t i = 0; i < imp_exec_rank_count; i++) 
    {
        const int neighbour_rank = imp_exec_list->ranks[i];
        const std::vector<int>& neighbour_rank_local_idxs = recv_buffers[i];

        for (int k = 0; k < imp_exec_list->sizes[i]; k++) 
        {       
            opp_particle_comm_data comm_data;
            comm_data.cell_residing_rank = neighbour_rank;
            comm_data.local_index = neighbour_rank_local_idxs[k];

            const int local_index = (set->cells_set->size + imp_exec_list->disps[i] + k);

            opp_part_comm_neighbour_data[set].insert({local_index, comm_data});

            // opp_printf("opp_part_comm_init", "set:[%s] li:[%d] nr:[%d] ni:[%d]", 
            //     set->cells_set->name, local_index, comm_data.cell_residing_rank, comm_data.local_index);
        }
    }  

    // create communication buffers in set
    opp_part_all_neigh_comm_data* mpi_buffers = new opp_part_all_neigh_comm_data();
    mpi_buffers->total_recv = 0;
    mpi_buffers->send_req.clear();
    mpi_buffers->recv_req.clear();

    for (size_t i = 0; i < imp_exec_rank_count; i++) 
    {
        const int neighbour_rank = imp_exec_list->ranks[i];

        mpi_buffers->neighbours.push_back(neighbour_rank);
        mpi_buffers->import_counts[neighbour_rank] = 0;

        opp_part_neigh_buffers& part_buffer = mpi_buffers->buffers[neighbour_rank];
        part_buffer.buf_import           = nullptr;
        part_buffer.buf_import_capacity  = -1;
        part_buffer.buf_import_index     = 0;
    }

    for (size_t i = 0; i < exp_exec_rank_count; i++) 
    {
        const int neighbour_rank = exp_exec_list->ranks[i];

        if (std::find(mpi_buffers->neighbours.begin(), mpi_buffers->neighbours.end(), neighbour_rank) 
                        == mpi_buffers->neighbours.end())
        {
            if (OPP_DBG)
                opp_printf("opp_part_set_comm_init", "UNLIKELY, %d is on export list, but not on import list of set %s", 
                    neighbour_rank, set->name); 
            mpi_buffers->neighbours.push_back(neighbour_rank);
        }

        mpi_buffers->export_counts[neighbour_rank] = 0;

        opp_part_neigh_buffers& part_buffer = mpi_buffers->buffers[neighbour_rank];
        part_buffer.buf_export           = nullptr;
        part_buffer.buf_export_capacity  = -1;
        part_buffer.buf_export_index     = 0;
    }

    set->mpi_part_buffers = (void*)mpi_buffers;

    opp_part_move_indices[set->index].clear();

    if (OPP_DBG) opp_printf("opp_part_set_comm_init", "set: %s END", set->name); 
}

//*******************************************************************************
void opp_part_comm_destroy()
{
    if (OPP_DBG) opp_printf("opp_part_comm_destroy", "START"); 

    for (int i = 0; i < (int)opp_sets.size(); i++)
    {
        opp_set set = opp_sets[i];

        if (!set->is_particle) continue;

        opp_part_all_neigh_comm_data* all_buffers = (opp_part_all_neigh_comm_data*)set->mpi_part_buffers;

        for (const int rank : all_buffers->neighbours)
        {
            opp_part_neigh_buffers& buffer = all_buffers->buffers[rank];

            if (buffer.buf_export) opp_host_free(buffer.buf_export);
            if (buffer.buf_import) opp_host_free(buffer.buf_import);
        }

        delete all_buffers;
    }

    opp_part_comm_neighbour_data.clear();
    opp_part_move_indices.clear();

    if (OPP_DBG) opp_printf("opp_part_comm_destroy", "END"); 
}

//*******************************************************************************
//*******************************************************************************
using namespace opp;


dh_particle_packer::dh_particle_packer(std::map<int, std::map<int, std::vector<opp_part_move_info>>>* part_move_data) 
    : part_move_data(part_move_data) {

}

//*******************************************************************************
dh_particle_packer::~dh_particle_packer() {

    this->part_move_data = nullptr;
    this->buffers.clear();
};

//*******************************************************************************
void dh_particle_packer::pack(opp_set set) {

    if (OPP_DBG) 
        opp_printf("dh_particle_packer", "pack set [%s]", set->name);

    if (part_move_data == nullptr) {
        opp_abort(std::string("part_move_data is NULL in dh_particle_packer"));
    }

    std::map<int, std::vector<char>>& buffers_of_set = this->buffers[set->index];
    for (auto& x : buffers_of_set) // try to keep the allocated vectors as it is, without deleting
        x.second.clear();

    if (part_move_data->at(set->index).size() == 0) {
        if (OPP_DBG) 
            opp_printf("dh_particle_packer", "Nothing to be sent for set [%s]", set->name);
        return;
    }     

    opp_profiler->start("MvDH_Pack");

    for (auto& move_idxs_per_rank : part_move_data->at(set->index)) {

        const int send_rank = move_idxs_per_rank.first;
        const std::vector<opp_part_move_info>& move_idxs_vec = move_idxs_per_rank.second;

        const size_t bytes_per_rank = (size_t)set->particle_size * move_idxs_vec.size();
        
        std::vector<char>& buffer_per_rank = buffers_of_set[send_rank];
        buffer_per_rank.resize(bytes_per_rank, 0);

        int displacement = 0;
        for (auto& dat : *(set->particle_dats)) {

            int dat_size = dat->size;

            if (dat->is_cell_index) {

                for (const auto& part : move_idxs_vec) {

                    memcpy(&(buffer_per_rank[displacement]), &part.foreign_cell_index, dat->size);
                    displacement += dat_size;
                }
            }
            else {

                for (const auto& part : move_idxs_vec) {

                    // copy the dat value to the send buffer
                    memcpy(&(buffer_per_rank[displacement]), 
                        &(dat->data[part.local_index * dat->size]), dat->size);
                    displacement += dat_size;
                }                
            }
        }

        if (OPP_DBG)
            opp_printf("dh_particle_packer", "Packed %zu parts to send to rank %d, displacement %d", 
                buffer_per_rank.size(), send_rank, displacement);
    }

    opp_profiler->end("MvDH_Pack");
}

//*******************************************************************************
char* dh_particle_packer::get_buffer(const opp_set set, const int send_rank) {

    auto set_it = this->buffers.find(set->index);

    if (set_it == this->buffers.end()) {
        opp_abort(std::string("Set not found in dh_particle_packer buffers"));
    } 

    auto set_itRank = set_it->second.find(send_rank);
    if (set_itRank == set_it->second.end()) {
        
        if (OPP_DBG) 
            opp_printf("dh_particle_packer", "get_buffer set [%s] Rank [%d] does not have a buffer created", 
                set->name, send_rank);
        return nullptr;
    } 

    if (OPP_DBG)
        opp_printf("dh_particle_packer", "get_buffer set [%s] Rank [%d] size %zu bytes", 
                set->name, send_rank, set_itRank->second.size());

    return &(set_itRank->second[0]);
}

//*******************************************************************************
void dh_particle_packer::unpack(opp_set set, const std::map<int, std::vector<char>>& particleRecvBuffers,
                    int64_t totalParticlesToRecv, const std::vector<int64_t>& recvRankPartCounts) {

    if (OPP_DBG) 
        opp_printf("dh_particle_packer", "unpack set [%s]", set->name);

    opp_profiler->start("MvDH_Unpack");

    if (totalParticlesToRecv > 0)
    {
        std::vector<opp_dat>& particle_dats = *(set->particle_dats);

        int64_t particle_size = set->particle_size;

        if (!opp_increase_particle_count_core(set, (int)totalParticlesToRecv)) // TODO : make int to int64_t
        {
            opp_printf("dh_particle_packer Unpack", "Error: Failed to increase particle count of particle set [%s]", 
                        set->name);
            opp_abort("dh_particle_packer Unpack error");
        }

        int64_t newPartIndex = (int64_t)(set->size - set->diff);
        // int rankx = 0;

        for (const auto& x : particleRecvBuffers) {

            // int recvRank = x.first;
            const std::vector<char>& buffer = x.second;

            int64_t recvCount = ((int64_t)buffer.size() / particle_size) ; // recvRankPartCounts[rankx++];
            int64_t displacement = 0;

            for (auto& dat : particle_dats)
            {
                memcpy(dat->data + newPartIndex * dat->size, &(buffer[displacement]), 
                    dat->size * recvCount);
                
                displacement += dat->size * recvCount; //(dat->size * recvCount);            
            }

            newPartIndex += recvCount;
        }
    }

    opp_profiler->end("MvDH_Unpack");

    if (OPP_DBG) 
        opp_printf("dh_particle_packer", "Unpack END");
}


//*******************************************************************************
GlobalParticleMover::GlobalParticleMover(MPI_Comm comm) 
    : comm(comm) {

    CHECK(MPI_Win_allocate(sizeof(int), sizeof(int), MPI_INFO_NULL, this->comm,
                    &this->recv_win_data, &this->recv_win));
    
    dh_part_move_data.clear();

    packer = std::make_unique<dh_particle_packer>(&dh_part_move_data);

    opp_profiler->reg("MvDH_WaitRanks");
    opp_profiler->reg("MvDH_Init");
    opp_profiler->reg("MvDH_Comm");
    opp_profiler->reg("MvDH_Finalize");
    opp_profiler->reg("MvDH_Pack");
    opp_profiler->reg("MvDH_Unpack");
    opp_profiler->reg("MvDH_WaitEx1");
    opp_profiler->reg("MvDH_WaitEx2");
    opp_profiler->reg("MvDH_WaitEx3");
    opp_profiler->reg("MvDH_WaitFin1");
    opp_profiler->reg("MvDH_WaitFin2");
}

//*******************************************************************************
GlobalParticleMover::~GlobalParticleMover() { 
    
    MPI_Barrier(this->comm);
    CHECK(MPI_Win_free(&this->recv_win));
}

//*******************************************************************************
void GlobalParticleMover::markParticleToMove(opp_set set, int partIndex, int rankToBeMoved, int finalGlobalCellIndex) {
    
    // These validations should be already done
    // if (finalGlobalCellIndex == MAX_CELL_INDEX) {
    //     opp_printf("GlobalParticleMover", "Error markParticleToMove particle %d will be moved to rank %d but global index is invalid",
    //         partIndex, rankToBeMoved);
    //     return;
    // }

    // if (rankToBeMoved == MAX_CELL_INDEX) {
    //     opp_printf("GlobalParticleMover", "Error markParticleToMove particle %d will be moved to finalGlobalCellIndex %d but rank is invalid",
    //         partIndex, rankToBeMoved);
    //     return;
    // }

    std::vector<opp_part_move_info>& vec = this->dh_part_move_data[set->index][rankToBeMoved];
    vec.emplace_back(partIndex, finalGlobalCellIndex);
}

//*******************************************************************************
void GlobalParticleMover::initGlobalMove() {
    
    opp_profiler->start("MvDH_Init");

    this->h_send_ranks.clear();
    this->h_send_rank_npart.clear();
    this->h_recv_ranks.clear();
    this->h_recv_rank_npart.clear();
    this->h_send_requests.clear();
    this->h_recv_requests.clear();
    this->h_recv_status.clear();

    this->totalParticlesToSend = 0;
    this->totalParticlesToRecv = 0;
    this->numRemoteRecvRanks = 0;
    this->numRemoteSendRanks = 0;

    for (auto& x : particleRecvBuffers)
        x.second.clear();

    this->recv_win_data[0] = 0;
    CHECK(MPI_Ibarrier(this->comm, &this->mpi_request));

    opp_profiler->end("MvDH_Init");
}

//*******************************************************************************
void GlobalParticleMover::communicateParticleSendRecvRankCounts() {
    
    opp_profiler->start("MvDH_WaitRanks");

    const int one[1] = {1};
    int recv[1];

    CHECK(MPI_Wait(&this->mpi_request, MPI_STATUS_IGNORE));

    for (auto& x : this->h_send_ranks) {
        
        const int rank = x;

        if (OPP_DBG)
            opp_printf("GlobalParticleMover", "locking rank %d", rank);

        CHECK(MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, this->recv_win));
        CHECK(MPI_Get_accumulate(one, 1, MPI_INT, recv, 1, MPI_INT, rank, 0, 1,
                                    MPI_INT, MPI_SUM, this->recv_win));
        CHECK(MPI_Win_unlock(rank, this->recv_win));
    }

    CHECK(MPI_Ibarrier(this->comm, &this->mpi_request));   

    opp_profiler->end("MvDH_WaitRanks");
}

//*******************************************************************************
void GlobalParticleMover::communicate(opp_set set) {
    
    opp_profiler->start("MvDH_Comm");

    std::map<int, std::vector<opp_part_move_info>>& rankVsPartData = dh_part_move_data[set->index]; 
    this->numRemoteSendRanks = rankVsPartData.size();
    this->h_send_requests.resize(this->numRemoteSendRanks);
    this->h_send_ranks.resize(this->numRemoteSendRanks);
    this->h_send_rank_npart.resize(this->numRemoteSendRanks);

    int rankx = 0;
    for (auto& x : rankVsPartData) {

        if (x.first >= OPP_comm_size || x.first < 0) {
            opp_printf("GlobalParticleMover", "ERROR locking rank %d [size %zu] from rank %d", x.first, OPP_rank, x.second.size());
            this->numRemoteSendRanks -= 1;
            continue;
        }

        this->h_send_ranks[rankx] = x.first;
        this->h_send_rank_npart[rankx] = (int64_t)x.second.size();
        this->totalParticlesToSend += (int64_t)x.second.size();
        rankx++;
    }

    packer->pack(set);

    // (1) Gather the number of ranks that is going to send particles to the current rank ------------

    communicateParticleSendRecvRankCounts();

    opp_profiler->start("MvDH_WaitEx1");
    CHECK(MPI_Wait(&this->mpi_request, MPI_STATUS_IGNORE));
    opp_profiler->end("MvDH_WaitEx1");

    // At this point, the current rank knows how many particles to recv
    this->numRemoteRecvRanks = this->recv_win_data[0];

    this->h_recv_rank_npart.resize(this->numRemoteRecvRanks);
    this->h_recv_requests.resize(this->numRemoteRecvRanks);
    
    // (2) Exchange particle counts ----------------------------------------------------------------  

    // non-blocking recv of particle counts
    for (int rankx = 0; rankx < this->numRemoteRecvRanks; rankx++) {

        CHECK(MPI_Irecv(&(this->h_recv_rank_npart[rankx]), 1, MPI_INT64_T, MPI_ANY_SOURCE, 
            42, this->comm, &(this->h_recv_requests[rankx])));
    }

    // non-blocking send of particle counts
    for (int rankx = 0; rankx < this->numRemoteSendRanks; rankx++) {
        CHECK(MPI_Isend(&this->h_send_rank_npart[rankx], 1, MPI_INT64_T, this->h_send_ranks[rankx], 
            42, this->comm, &this->h_send_requests[rankx]));
    }

    this->h_recv_ranks.resize(this->numRemoteRecvRanks);
    this->h_recv_status.resize(this->numRemoteRecvRanks);

    opp_profiler->start("MvDH_WaitEx2");
    CHECK(MPI_Waitall(this->numRemoteRecvRanks, &(this->h_recv_requests[0]), &(this->h_recv_status[0])));
    opp_profiler->end("MvDH_WaitEx2");

    for (int rankx = 0; rankx < this->numRemoteRecvRanks; rankx++) {

        const int remote_rank = this->h_recv_status[rankx].MPI_SOURCE;

        if (OPP_DBG) { 

            if (remote_rank < 0 && remote_rank < OPP_comm_size) 
                opp_abort(std::string("Recv rank is invalid"));
            if (this->h_recv_rank_npart[rankx] <= 0) 
                opp_abort(std::string("A remote rank is trying to send 0 (or fewer) particles to this rank"));
        }

        this->h_recv_ranks[rankx] = remote_rank;
        this->totalParticlesToRecv += this->h_recv_rank_npart[rankx];

        std::vector<char>& partRecvBufferPerRank = particleRecvBuffers[remote_rank];
        int64_t BytesToRecv = this->h_recv_rank_npart[rankx] * (int64_t)set->particle_size;
        partRecvBufferPerRank.resize(BytesToRecv);
    }

    if (OPP_DBG)
        opp_printf("GlobalMove", "communicate %lld", this->totalParticlesToRecv);

    opp_profiler->start("MvDH_WaitEx3");
    CHECK(MPI_Waitall(this->numRemoteSendRanks, &(this->h_send_requests[0]), MPI_STATUSES_IGNORE));
    opp_profiler->end("MvDH_WaitEx3");

    // (3) send and receive the particles -----------------------------------------------------------
    
    this->h_send_requests.clear();
    this->h_recv_requests.clear();
    this->h_send_requests.resize(this->numRemoteSendRanks);            
    this->h_recv_requests.resize(this->numRemoteRecvRanks);
    
    for (int rankx = 0; rankx < this->numRemoteRecvRanks; rankx++) {
        
        const int recvRank = this->h_recv_ranks[rankx];
        std::vector<char>& partRecvBufferPerRank = this->particleRecvBuffers[recvRank];
        int64_t recvSize = (this->h_recv_rank_npart[rankx] * set->particle_size);

        // opp_printf("Communicate", "Expected to recv %d bytes from rank %d, buffer size %zu", 
        //     recvSize, recvRank, partRecvBufferPerRank.size());

        CHECK(MPI_Irecv(&(partRecvBufferPerRank[0]), recvSize, MPI_CHAR, recvRank, 43, 
                this->comm, &this->h_recv_requests[rankx]));
    }

    // non-blocking send of particle data
    for (int rankx = 0; rankx < this->numRemoteSendRanks; rankx++) {
        
        const int sendRank = this->h_send_ranks[rankx];
        char* sendBuffer = packer->get_buffer(set, sendRank);
        int64_t sendSize = (this->h_send_rank_npart[rankx] * set->particle_size);

        // opp_printf("Communicate", "Trying to send %d bytes to rank %d buffer %p", 
        //     sendSize, sendRank, sendBuffer);

        CHECK(MPI_Isend(sendBuffer, sendSize, MPI_CHAR, sendRank, 43, this->comm, 
                &this->h_send_requests[rankx]));
    }   

    // (4) Once sent, map could be cleared for the set, keep the allocations if possible -----------

    for (auto& x : dh_part_move_data[set->index])
        x.second.clear();       

    opp_profiler->end("MvDH_Comm");  
}

//*******************************************************************************
int64_t GlobalParticleMover::finalize(opp_set set) {
    
    opp_profiler->start("MvDH_Finalize");

    opp_profiler->start("MvDH_WaitFin1");
    CHECK(MPI_Waitall(this->numRemoteRecvRanks, &(this->h_recv_requests[0]), &(this->h_recv_status[0])));
    opp_profiler->end("MvDH_WaitFin1");

    // if (OPP_DBG) 
    {
        // Check this rank recv'd the correct number of bytes from each remote
        for (int rankx = 0; rankx < this->numRemoteRecvRanks; rankx++) {
        
            int bytesRecvd = -1;
            CHECK(MPI_Get_count(&this->h_recv_status[rankx], MPI_CHAR, &bytesRecvd));
            int bytesExpected = (this->h_recv_rank_npart[rankx] * set->particle_size);

            if (bytesRecvd != bytesExpected) {
                
                opp_printf("GlobalParticleMover", "recv'd incorrect num of bytes Expected %d Received %d",
                    bytesExpected, bytesRecvd);
                opp_abort(std::string("recv'd incorrect number of bytes"));
            }
        }
    }

    packer->unpack(set, this->particleRecvBuffers, this->totalParticlesToRecv, this->h_recv_rank_npart);

    // Note : The mesh relations received will be global cell indices and need to be converted to local

    opp_profiler->start("MvDH_WaitFin2");
    CHECK(MPI_Waitall(this->numRemoteSendRanks, &(this->h_send_requests[0]), MPI_STATUSES_IGNORE));
    opp_profiler->end("MvDH_WaitFin2");

    // Since packer->unpack might realloc dats, need to set the correct OPP_mesh_relation_data ptr
    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 

    opp_profiler->end("MvDH_Finalize");

    return this->totalParticlesToRecv;
}
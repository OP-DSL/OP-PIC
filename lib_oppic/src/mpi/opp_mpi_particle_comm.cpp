
#include <opp_mpi.h>

#define MPI_COUNT_EXCHANGE 0
#define MPI_TAG_PART_EX 1

// below translate to std::map<set_index, std::map<local_cell_index, opp_particle_comm_data>>
std::map<int, std::map<int, opp_particle_comm_data>> opp_part_comm_neighbour_data; 

//*******************************************************************************
void opp_pack_particle(oppic_set set, int index, int send_rank)
{
    if (OP_DEBUG) opp_printf("opp_pack_particle", " set [%s] | index %d | send_rank %d", set->name, index, send_rank);

    opp_all_mpi_part_buffers* send_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;

    // check whether send_rank is a neighbour or not
    if (OP_DEBUG)
    {
        if (std::find(send_buffers->neighbours.begin(), send_buffers->neighbours.end(), send_rank) == send_buffers->neighbours.end())
        {
            opp_printf("opp_pack_particle", "Error: send_rank %d is not a neighbour, cannot send index %d of set [%s]",
                send_rank, index, set->name);
            MPI_Abort(OP_MPI_WORLD, 1);
        }
    }

    opp_mpi_part_buffer& send_rank_buffer = send_buffers->buffers[send_rank];

    // resize the export buffer if required
    if (send_rank_buffer.buf_export_index >= send_rank_buffer.buf_export_capacity)
    {
        if (send_rank_buffer.buf_export == nullptr)
        {
            send_rank_buffer.buf_export_capacity  = OPP_mpi_part_alloc_mult * set->particle_size;
            send_rank_buffer.buf_export_index     = 0;
            send_rank_buffer.buf_export           = (char *)malloc(send_rank_buffer.buf_export_capacity);
            //memset(send_rank_buffer.buf_export, 0, send_rank_buffer.buf_export_capacity); // not essential, can remove
        }
        else
        {
            send_rank_buffer.buf_export_capacity += OPP_mpi_part_alloc_mult * set->particle_size;
            send_rank_buffer.buf_export           = (char *)realloc(send_rank_buffer.buf_export, send_rank_buffer.buf_export_capacity);
            // memset(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_capacity - OPP_mpi_part_alloc_mult * set->particle_size]),
            //     0, OPP_mpi_part_alloc_mult * set->particle_size); // not essential, can remove
        }
    }

    std::vector<oppic_dat>& particle_dats = *(set->particle_dats);
    int displacement = 0;

    // pack the particle data from dats into the export buffer 
    for (int i = 0; i < particle_dats.size(); i++)
    {
        oppic_dat& dat = particle_dats[i];

        memcpy(&(send_rank_buffer.buf_export[send_rank_buffer.buf_export_index + displacement]), &(dat->data[index * dat->size]), dat->size);

        if (OP_DEBUG)
        {
            char* part_buffer = &send_rank_buffer.buf_export[send_rank_buffer.buf_export_index];
            std::string log = "";
            if (strcmp(dat->type, "double") == 0)
            {
                double* d = (double*)(part_buffer + displacement);
                for (int l = 0; l < dat->dim; l++) log += " " + std::to_string(d[l]);
            }
            else if (strcmp(dat->type, "int") == 0)
            {
                int* d = (int*)(part_buffer + displacement);
                for (int l = 0; l < dat->dim; l++) log += " " + std::to_string(d[l]);
            }

            opp_printf("opp_pack_particle", "%s from index %d -%s", dat->name, index, log.c_str());
        }

        displacement += dat->size;
    } 

    send_rank_buffer.buf_export_index += set->particle_size;
    (send_buffers->export_counts)[send_rank] += 1;

    if (OP_DEBUG) opp_printf("opp_pack_particle", "END send_rank %d exported count %d", send_rank, (send_buffers->export_counts)[send_rank]);
}

//*******************************************************************************
void opp_unpack_particles(oppic_set set)
{
    if (OP_DEBUG) opp_printf("opp_unpack_particles", "set [%s]", set->name);

    std::vector<oppic_dat>& particle_dats = *(set->particle_dats);
    int num_particles = 0;

    opp_all_mpi_part_buffers* receive_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;
    std::vector<int>& neighbours = receive_buffers->neighbours;

    // count the number of particles to be received from all ranks
    for (int i = 0; i < neighbours.size(); i++)
    {
        int neighbour_rank = neighbours[i];
        num_particles += (receive_buffers->import_counts)[neighbour_rank];
    }

    if (num_particles > 0)
    {
        int particle_size = set->particle_size;

        oppic_increase_particle_count_core(set, num_particles);
        int new_part_index = (set->size - set->diff);

        for (int i = 0; i < neighbours.size(); i++)
        {
            int neighbour_rank = neighbours[i];

            opp_mpi_part_buffer& receive_rank_buffer = receive_buffers->buffers[neighbour_rank];
            int receive_count = receive_buffers->import_counts[neighbour_rank];

            // unpack the received buffer from rank 'neighbour_rank' in to particle dats
            for (int part = 0; part < receive_count; part++)
            {
                char* part_buffer = &(receive_rank_buffer.buf_import[particle_size * part]);
                int displacement = 0;

                for (int i = 0; i < particle_dats.size(); i++)
                {
                    oppic_dat& dat = particle_dats[i];

                    memcpy(dat->data + new_part_index * dat->size, part_buffer + displacement, dat->size);

                    if (OP_DEBUG)
                    {
                        std::string log = "";
                        if (strcmp(dat->type, "double") == 0)
                        {
                            double* d = (double*)(part_buffer + displacement);
                            for (int l = 0; l < dat->dim; l++) log += " " + std::to_string(d[l]);
                        }
                        else if (strcmp(dat->type, "int") == 0)
                        {
                            int* d = (int*)(part_buffer + displacement);
                            for (int l = 0; l < dat->dim; l++) log += " " + std::to_string(d[l]);
                        }

                        opp_printf("opp_unpack_particles", "%s to index %d -%s", dat->name, new_part_index, log.c_str());
                    }

                    displacement += dat->size;     
                }

                new_part_index++;
            }
        }
    }

    if (OP_DEBUG) opp_printf("opp_unpack_particles", "END");
}

//*******************************************************************************
bool opp_check_part_need_comm(int map0idx, oppic_set set, int particle_index, opp_move_var& m) 
{

    if (m.OPP_move_status != OPP_NEED_MOVE)
    {
        return false; // The particle is already found to be either MOVE_DONE or NEED_REMOVE
    }
    else if (map0idx >= set->cells_set->size)
    {
        // map0idx cell is not owned by the current mpi rank (it is in the import exec halo region), need to communicate

        std::map<int, opp_particle_comm_data>& set_part_com_data = opp_part_comm_neighbour_data[set->index];

        auto it = set_part_com_data.find(map0idx);
        if (it == set_part_com_data.end())
        {
            opp_printf("opp_check_part_need_comm", "Error: cell %d cannot be found in opp_part_comm_neighbour_data map", map0idx);
            return false; // unlikely, need exit(-1) to abort instead!
        }

        opp_particle_comm_data& comm_data = it->second;
        
        // change the cell_index to reflect the correct neighbour ranks local cell index before packing
        ((int*)set->mesh_relation_dat->data)[particle_index] = comm_data.local_index;

        opp_pack_particle(set, particle_index, comm_data.cell_residing_rank);
        
        // This particle is already packed, hence needs to be removed from the current rank
        m.OPP_move_status = OPP_NEED_REMOVE; 

        return true;
    }

    // map0idx is an own cell
    return false;
}

//*******************************************************************************
void opp_exchange_particles(oppic_set set)
{
    if (OP_DEBUG) opp_printf("opp_exchange_particles", "set %s - particle size", set->name, set->particle_size);

    opp_all_mpi_part_buffers* mpi_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;

    std::vector<int>& neighbours = mpi_buffers->neighbours;
    int neighbour_count = neighbours.size();
    mpi_buffers->total_recv = 0;

    for (auto it = mpi_buffers->import_counts.begin(); it != mpi_buffers->import_counts.end(); it++)
        it->second = 0;

    mpi_buffers->recv_req.clear();
    mpi_buffers->send_req.clear();

    std::vector<MPI_Request> send_req_countEx(neighbour_count);
    std::vector<MPI_Request> recv_req_countEx(neighbour_count);

    // send/receive send_counts to/from only to neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        int neighbour_rank = neighbours[i];

        int& send_count = mpi_buffers->export_counts[neighbour_rank];
        MPI_Isend(&send_count, 1, MPI_INT, neighbour_rank, MPI_COUNT_EXCHANGE, OP_MPI_WORLD, &(send_req_countEx[i]));

        int& recv_count = mpi_buffers->import_counts[neighbour_rank];
        MPI_Irecv(&recv_count, 1, MPI_INT, neighbour_rank, MPI_COUNT_EXCHANGE, OP_MPI_WORLD, &(recv_req_countEx[i]));
    }

    // send the particle data only to neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        int neighbour_rank = neighbours[i];
        int send_size = set->particle_size * mpi_buffers->export_counts[neighbour_rank];

        if (send_size <= 0)
        {
            if (OP_DEBUG) opp_printf("opp_exchange_particles", "nothing to send to rank %d", neighbour_rank);
            continue;
        }
        else
        {   
            if (OP_DEBUG) opp_printf("opp_exchange_particles", "sending %d particle/s (size: %d) to rank %d", 
                (send_size/ set->particle_size), send_size, neighbour_rank);
        }

        char* send_buffer = mpi_buffers->buffers[neighbour_rank].buf_export;

        // opp_printf("opp_exchange_particles", "neighbour_rank: %d, send_size: %d", neighbour_rank, send_size);

        MPI_Request req;
        MPI_Isend(send_buffer, send_size, MPI_CHAR, neighbour_rank, MPI_TAG_PART_EX, OP_MPI_WORLD, &req);
        mpi_buffers->send_req.push_back(req);
    }

    // wait for the counts to receive only from neighbours
    MPI_Waitall(neighbour_count, &recv_req_countEx[0], MPI_STATUSES_IGNORE);

    // create/resize data structures and receive particle data from neighbours
    for (int i = 0; i < neighbour_count; i++)
    {
        int neighbour_rank = neighbours[i];
        int recv_size = set->particle_size * mpi_buffers->import_counts[neighbour_rank];
        mpi_buffers->total_recv += mpi_buffers->import_counts[neighbour_rank];

        if (recv_size <= 0)
        {
            if (OP_DEBUG) opp_printf("opp_exchange_particles", "nothing to receive from rank %d", neighbour_rank);
            continue;
        }

        // opp_printf("opp_exchange_particles", "neighbour_rank: %d, recv_size: %d, total_recv: %d", 
        //     neighbour_rank, recv_size, mpi_buffers->total_recv);

        opp_mpi_part_buffer& receive_buffer = mpi_buffers->buffers[neighbour_rank];

        if (recv_size >= receive_buffer.buf_import_capacity)
        {
            if (receive_buffer.buf_import == nullptr)
            {
                receive_buffer.buf_import_capacity  = OPP_mpi_part_alloc_mult * set->particle_size;           
                receive_buffer.buf_import           = (char *)malloc(receive_buffer.buf_import_capacity);
                // memset(receive_buffer.buf_import, 0 ,receive_buffer.buf_import_capacity); // not essential, can remove            
            }
            else
            {
                receive_buffer.buf_import_capacity += OPP_mpi_part_alloc_mult * set->particle_size;
                receive_buffer.buf_import           = (char *)realloc(receive_buffer.buf_import, receive_buffer.buf_import_capacity);
                // memset(&(receive_buffer.buf_import[receive_buffer.buf_import_capacity - OPP_mpi_part_alloc_mult * set->particle_size]),
                //     0, OPP_mpi_part_alloc_mult * set->particle_size); // not essential, can remove
            }
        }

        MPI_Request req;
        MPI_Irecv(receive_buffer.buf_import, recv_size, MPI_CHAR, neighbour_rank, MPI_TAG_PART_EX, OP_MPI_WORLD, &req);
        mpi_buffers->recv_req.push_back(req);
    }

    // reset the export counts for another iteration
    for (auto it = mpi_buffers->export_counts.begin(); it != mpi_buffers->export_counts.end(); it++)
    {
        it->second = 0; // make the export count to zero for the next iteration
        mpi_buffers->buffers[it->first].buf_export_index = 0; // make the export index of that rank to zero for the next iteration
    }

    if (OP_DEBUG) opp_printf("opp_exchange_particles", "END");
}

//*******************************************************************************
void opp_wait_all_particles(oppic_set set)
{
    if (OP_DEBUG) opp_printf("opp_wait_all_particles", "START");

    opp_all_mpi_part_buffers* mpi_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;

    std::vector<MPI_Request>& send_req = mpi_buffers->send_req;
    std::vector<MPI_Request>& recv_req = mpi_buffers->recv_req;

    // std::vector<MPI_Status> recv_status(recv_req.size());
    // std::vector<MPI_Status> send_status(send_req.size());

    // wait till all the particles from all the ranks are received
    MPI_Waitall(send_req.size(), &(send_req[0]), MPI_STATUSES_IGNORE); // &(send_status[0])); //
    MPI_Waitall(recv_req.size(), &(recv_req[0]), MPI_STATUSES_IGNORE); // &(recv_status[0])); //

    send_req.clear();
    recv_req.clear();

    // increase the particle count if required and unpack the communicated particle buffer in to separate particle dats
    opp_unpack_particles(set);

    if (OP_DEBUG) opp_printf("opp_wait_all_particles", "END");
}

//*******************************************************************************
bool opp_check_all_done(oppic_set set)
{
    if (OP_DEBUG) opp_printf("opp_check_all_done", "START");

    opp_all_mpi_part_buffers* mpi_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;
    bool imported_parts = true;
    
    if (mpi_buffers->total_recv == 0)
    {
        imported_parts = false;
        set->diff = 0;
    }

    bool bool_ret = false;
    bool* buffer_recv = (bool *)malloc(OPP_comm_size * sizeof(bool));

    // gather from all MPI ranks to see whether atleast one rank needs to communicate to another
    MPI_Allgather(&imported_parts, 1, MPI_C_BOOL, buffer_recv, 1, MPI_C_BOOL, OP_MPI_WORLD);

    std::string log = "";
    for (int i = 0; i < OPP_comm_size; i++)
    {
        bool_ret = (bool_ret || buffer_recv[i]);

        if (OP_DEBUG) 
            log += std::string(" R") + std::to_string(i) + (buffer_recv[i] ? "-T" : "-F");
    }

    if (OP_DEBUG) 
        opp_printf("opp_check_all_done", "%s -%s", (bool_ret ? "ITER AGAIN" : "ALL DONE"), log.c_str());

    return !bool_ret;
}

//*******************************************************************************
void opp_particle_comm_init()
{
    if (OP_DEBUG) opp_printf("opp_particle_comm_init", "START");

    for (oppic_set& set : oppic_sets)
    {
        if (!set->is_particle) continue;

        opp_particle_set_comm_init(set);
    }

    if (OP_DEBUG) opp_printf("opp_particle_comm_init", "END");
}

//*******************************************************************************
void opp_particle_set_comm_init(oppic_set set)
{
    // TODO : can use the same mappings for all particle sets with the same cells set, instead of communicating again

    if (OP_DEBUG) opp_printf("opp_particle_set_comm_init", "set: %s", set->name);

    halo_list exp_exec_list = OP_export_exec_list[set->cells_set->index];
    halo_list imp_exec_list = OP_import_exec_list[set->cells_set->index];

    std::vector<MPI_Request> send_req;
    std::vector<MPI_Request> recv_req;    

    // send the local index of the export elements of the set to all neighbours
    for (int i = 0; i < exp_exec_list->ranks_size; i++) 
    {  
        int neighbour_rank = exp_exec_list->ranks[i];
        int* send_buffer = &(exp_exec_list->list[exp_exec_list->disps[i]]);
        int send_size = exp_exec_list->sizes[i];
    
        MPI_Request req;
        MPI_Isend(send_buffer, send_size, MPI_INT, neighbour_rank, OPP_my_rank, OP_MPI_WORLD, &req);
        send_req.push_back(req);  

        //print the per rank send buffers
        if (OP_DEBUG)
        {
            std::string log = "";
            for (int k = 0; k < send_size; k++)
                log += std::to_string(send_buffer[k]) + " ";        
            opp_printf("opp_particle_set_comm_init", "%s SEND neighbour_rank %d send_size %d -> %s", 
                set->cells_set->name, neighbour_rank, send_size, log.c_str());
        } 
    }

    std::vector<std::vector<int>> recv_buffers(imp_exec_list->ranks_size);

    // receive the foreign ranks local index of the import elements of the set from all neighbours
    for (int i = 0; i < imp_exec_list->ranks_size; i++) 
    {
        int neighbour_rank = imp_exec_list->ranks[i];
        int recv_size = imp_exec_list->sizes[i];
        
        std::vector<int>& recv_buffer = recv_buffers[i];
        recv_buffer.resize(recv_size);
        
        MPI_Request req;
        MPI_Irecv(&recv_buffer[0], recv_size, MPI_INT, neighbour_rank, neighbour_rank, OP_MPI_WORLD, &req);
        recv_req.push_back(req);  
    }

    MPI_Waitall(recv_req.size(), &recv_req[0], MPI_STATUSES_IGNORE);

    // print the per rank received buffers
    if (OP_DEBUG)
    {
        for (int i = 0; i < recv_buffers.size(); i++) 
        {
            std::string log = "";
            for (int k = 0; k < recv_buffers[i].size(); k++)
                log += std::to_string((recv_buffers[i])[k]) + " ";  
            opp_printf("opp_particle_set_comm_init", "%s RECEIVE neighbour_rank %d recv_size %d -> %s", 
                set->cells_set->name, imp_exec_list->ranks[i], recv_buffers[i].size(), log.c_str());
        }
    }

    std::map<int, opp_particle_comm_data>& set_part_com_data = opp_part_comm_neighbour_data[set->index];

    // create mappings of neighbour ranks cell information for easy access during particle communication
    for (int i = 0; i < imp_exec_list->ranks_size; i++) 
    {
        int neighbour_rank = imp_exec_list->ranks[i];
        std::vector<int>& neighbour_rank_local_idxs = recv_buffers[i];

        for (int k = 0; k < imp_exec_list->sizes[i]; k++) 
        {       
            opp_particle_comm_data comm_data;
            comm_data.cell_residing_rank = neighbour_rank;
            comm_data.local_index = neighbour_rank_local_idxs[k];

            int local_index = (set->cells_set->size + imp_exec_list->disps[i] + k);

            set_part_com_data.insert({local_index, comm_data});

            // opp_printf("opp_particle_comm_init", OPP_my_rank, "cset:%s li:%d nr:%d ni:%d", 
            //     set->cells_set->name, local_index, comm_data.cell_residing_rank, comm_data.local_index);
        }
    }  

    // create communication buffers in set
    opp_all_mpi_part_buffers* mpi_buffers = new opp_all_mpi_part_buffers();
    
    mpi_buffers->total_recv = 0;
    mpi_buffers->send_req.clear();
    mpi_buffers->recv_req.clear();

    for (int i = 0; i < imp_exec_list->ranks_size; i++) 
    {
        mpi_buffers->neighbours.push_back(imp_exec_list->ranks[i]);

        mpi_buffers->import_counts[imp_exec_list->ranks[i]] = 0;

        opp_mpi_part_buffer& part_buffer = mpi_buffers->buffers[i];
        part_buffer.buf_import           = nullptr;
        part_buffer.buf_import_capacity  = -1;
        part_buffer.buf_import_index     = 0;
    }

    for (int i = 0; i < exp_exec_list->ranks_size; i++) 
    {
        if (std::find(mpi_buffers->neighbours.begin(), mpi_buffers->neighbours.end(), exp_exec_list->ranks[i]) == mpi_buffers->neighbours.end())
        {
            opp_printf("opp_particle_set_comm_init", "UNLIKELY, %d is on export list, but not on import list of set %s", 
                exp_exec_list->ranks[i], set->name); 
            mpi_buffers->neighbours.push_back(exp_exec_list->ranks[i]);
        }

        mpi_buffers->export_counts[exp_exec_list->ranks[i]] = 0;

        opp_mpi_part_buffer& part_buffer = mpi_buffers->buffers[i];
        part_buffer.buf_export           = nullptr;
        part_buffer.buf_export_capacity  = -1;
        part_buffer.buf_export_index     = 0;
    }

    set->mpi_part_buffers = (void*)mpi_buffers;

    if (OP_DEBUG) opp_printf("opp_particle_set_comm_init", "set: %s END", set->name); 
}

//*******************************************************************************
void opp_particle_comm_destroy()
{
    if (OP_DEBUG) opp_printf("opp_particle_comm_destroy", "START"); 

    for (int i = 0; i < oppic_sets.size(); i++)
    {
        oppic_set set = oppic_sets[i];

        if (!set->is_particle) continue;

        opp_all_mpi_part_buffers* all_buffers = (opp_all_mpi_part_buffers*)set->mpi_part_buffers;

        for (int rank : all_buffers->neighbours)
        {
            opp_mpi_part_buffer& buffer = all_buffers->buffers[rank];

            if (buffer.buf_export) 
                free(buffer.buf_export);

            if (buffer.buf_import) 
                free(buffer.buf_import);
        }

        delete all_buffers;
    }

    opp_part_comm_neighbour_data.clear();

    if (OP_DEBUG) opp_printf("opp_particle_comm_destroy", "END"); 
}
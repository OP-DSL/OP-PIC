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

#include "opp_omp.h"

void opp_part_pack(opp_set set);
void opp_part_unpack(opp_set set);
void opp_finalize_particle_move_omp(opp_set set);

//****************************************
void opp_init_particle_move(opp_set set, int nargs, opp_arg *args)
{
    part_remove_count_per_thr.resize(opp_nthreads);
    std::fill(part_remove_count_per_thr.begin(), part_remove_count_per_thr.end(), 0);

#ifdef USE_MPI
    opp_move_part_indices.clear();
    opp_move_part_indices.reserve(20000);

    move_part_indices_per_thr.resize(opp_nthreads);
    for (int t = 0; t < opp_nthreads; t++) 
    {
        move_part_indices_per_thr[t].clear();
        move_part_indices_per_thr[t].reserve(1000);
    }
    // std::fill(move_part_indices_per_thr.begin(), move_part_indices_per_thr.end(), std::vector<int>());

    gbl_move_indices_per_thr.resize(opp_nthreads);
    for (int t = 0; t < opp_nthreads; t++) 
    {
        gbl_move_indices_per_thr[t].clear();
        gbl_move_indices_per_thr[t].reserve(1000);
    }    
#endif

    if (OPP_comm_iteration == 0)
    {
        OPP_iter_start = 0;
        OPP_iter_end   = set->size;          
    }
    else
    {
        // need to change the arg data since particle communication could change the pointer in realloc dat->data
        for (int i = 0; i < nargs; i++)
        {
            if (args[i].argtype == OPP_ARG_DAT && args[i].dat->set->is_particle)
            {
                args[i].data = args[i].dat->data;
                if (OPP_DBG) 
                    opp_printf("opp_init_particle_move", "dat %s", args[i].dat->name);
            }
        }
    }

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 

    opp_init_particle_move_core(set);
}

//****************************************
bool opp_finalize_particle_move(opp_set set)
{ 
    if (OPP_DBG) 
        opp_printf("opp_finalize_particle_move", "Start particle set [%s]", set->name);

    opp_profiler->start("Mv_Finalize");

    for (int i = 0; i < opp_nthreads; i++) 
        set->particle_remove_count += part_remove_count_per_thr[i];

#ifdef USE_MPI
    // process indices of each thread now
    for (int i = 0; i < opp_nthreads; i++)
    {
        opp_move_part_indices.insert(opp_move_part_indices.end(), 
            move_part_indices_per_thr[i].begin(), move_part_indices_per_thr[i].end());
    }

    opp_process_marked_particles(set);

    opp_part_pack(set);

    // send the counts and send the particles  
    opp_part_exchange(set); 
#endif

    opp_finalize_particle_move_omp(set);

#ifdef USE_MPI
    if (opp_part_check_all_done(set))
    {
        if (OPP_max_comm_iteration < OPP_comm_iteration)
            OPP_max_comm_iteration = OPP_comm_iteration;

        OPP_comm_iteration = 0; // reset for the next par loop
        
        opp_profiler->end("Mv_Finalize");
        return false; // all mpi ranks do not have anything to communicate to any rank
    }

    opp_part_wait_all(set); // wait till all the particles are communicated

    if (OPP_DBG)
        opp_printf("opp_finalize_particle_move", "set [%s] size prior unpack %d", set->name, set->size);

    // increase the particle count if required and unpack the communicated particle buffer 
    // in to separate particle dats
    opp_part_unpack(set);    

    OPP_iter_start = set->size - set->diff;
    OPP_iter_end   = set->size;  

    OPP_comm_iteration++;  

    opp_profiler->end("Mv_Finalize");
    return true;
#else

    opp_profiler->end("Mv_Finalize");
    return false;
#endif
}

//****************************************
void opp_finalize_particle_move_omp(opp_set set)
{
    if (OPP_DBG) 
        opp_printf("opp_finalize_particle_move_omp", "set [%s] with particle_remove_count [%d]\n", 
        set->name, set->particle_remove_count);

    if (set->particle_remove_count <= 0) {
        return;
    }

    int *mesh_relation_data = (int *)set->mesh_relation_dat->data;
    std::vector<std::pair<size_t, size_t>> swap_indices;    // contain hole index and the index from back to swap

    // Idea: The last available element should be copied to the hole
    // In the below scope we try to calculate the element to be swapped with the hole
    {
        // set->particle_remove_count   // the particle count that should be removed
        int removed_count = 0;          // how many elements currently being removed
        int skip_count = 0;             // how many elements from the back is skipped ..
                                        // .. due to that element is also to be removed

        for (size_t j = 0; j < (size_t)set->size; j++) {
            
            // skip if the current index is not to be removed
            if (mesh_relation_data[j] != MAX_CELL_INDEX) 
                continue;

            // handle if the element from the back is also to be removed
            while ((set->size - removed_count - skip_count - 1 >= 0) && 
                    (mesh_relation_data[set->size - removed_count - skip_count - 1] == MAX_CELL_INDEX)) {
                skip_count++;
            }

            // check whether the holes are at the back!
            if ((set->size - removed_count - skip_count - 1 < 0) ||
                    (j >= (size_t)(set->size - removed_count - skip_count - 1))) {
                if (OPP_DBG) 
                    opp_printf("opp_finalize_particle_move_core", 
                    "Current Iteration index [%d] and replacement index %d; hence breaking [%s]", 
                    j, (set->size - removed_count - skip_count - 1), set->name);
                break;
            }

            swap_indices.push_back(std::make_pair(j, (size_t)(set->size - removed_count - skip_count - 1)));

            removed_count++;
        }
    }

    // For all the dats, fill the holes using the swap_indices
    #pragma omp parallel for
    for (size_t i = 0; i < set->particle_dats->size(); i++) {
        
        opp_dat dat = set->particle_dats->at(i);
        for (const auto& x : swap_indices) {
            char* dat_removed_ptr = (char *)(dat->data + (x.first * dat->size));
            size_t offset_byte = x.second * dat->size;
            char* dat_to_replace_ptr = (char *)(dat->data + offset_byte);
            
            memcpy(dat_removed_ptr, dat_to_replace_ptr, dat->size); 
        }
    }

    set->size -= set->particle_remove_count;
    set->particle_remove_count = 0;
}

// #endif

//****************************************
void opp_particle_sort(opp_set set)
{ 
    if (OPP_DBG) 
        opp_printf("opp_particle_sort", "set [%s]\n", set->name);
    
    OPP_INT* mesh_relation_data = (OPP_INT*)set->mesh_relation_dat->data;

    std::vector<size_t> idx_before_sort = sort_indexes(mesh_relation_data, set->size);

    #pragma omp parallel for
    for (size_t i = 0; i < set->particle_dats->size(); i++) {
            
        auto& dat = set->particle_dats->at(i);
        char *new_data = (char *)opp_host_malloc(set->set_capacity * dat->size);
        char *old_data = (char*)dat->data;
        
        for (int j = 0; j < set->size; j++) {
            memcpy(new_data + j * dat->size, old_data + idx_before_sort[j] * dat->size, dat->size);
        }

        opp_host_free(dat->data);
        dat->data = new_data;

        // Once sorted, below will gather the start and end indices of particles per cell index
        // if (dat->is_cell_index) { 
        //     int* cell_index_array = (int*)dat->data;
        //     int current_cell_index = -1, previous_cell_index = -1;
        //     std::map<int, part_index>& map = *(set->cell_index_v_part_index_map);
        //     map.clear();
        //     for (int j = 0; j < set->size; j++) {    
        //         current_cell_index = cell_index_array[j];
        //         if ((current_cell_index != previous_cell_index) && (current_cell_index >= 0)) {
        //             part_index& pi = map[current_cell_index];
        //             pi.start = j;
        //             if (previous_cell_index >= 0) 
        //                 map[previous_cell_index].end = (j - 1);
        //         }
        //         previous_cell_index = current_cell_index;
        //     }
        //     map[previous_cell_index].end = (set->size - 1);
        // }
    }
}

// //****************************************
// void opp_mark_particle_to_remove(opp_set set, int particle_index)
// {
//     opp_mark_particle_to_remove_core(set, particle_index);
// }

// //****************************************
// void opp_remove_marked_particles_from_set(opp_set set)
// {
//     opp_remove_marked_particles_from_set_core(set);
// }
// void opp_remove_marked_particles_from_set(opp_set set, std::vector<int>& idx_to_remove)
// {
//     opp_remove_marked_particles_from_set_core(set, idx_to_remove);
// }


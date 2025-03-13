
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

#include <opp_lib_core.h>

void particle_hole_fill_core(opp_set set);
void opp_particle_sort_core(opp_set set, bool shuffle);

//****************************************
void opp_init_particle_move_core(opp_set set)
{
    if (OPP_DBG) opp_printf("opp_init_particle_move_core", "set [%s]", set->name);

    OPP_part_cells_set_size = set->cells_set->size;
    
    set->particle_remove_count = 0;
}

//****************************************
void opp_finalize_particle_move_core(opp_set set)
{
    if (OPP_DBG) 
        opp_printf("opp_finalize_particle_move_core", 
            "set [%s] size[%d] with particle_remove_count [%d] diff [%d]", 
            set->name, set->size, set->particle_remove_count, set->diff);

    // return if there are no particles to be removed
    if (set->particle_remove_count <= 0) {
        return;
    }

    if (OPP_fill_type == OPP_HoleFill_All || 
        ((OPP_fill_type == OPP_Sort_Periodic || OPP_fill_type == OPP_Shuffle_Periodic) && 
            (OPP_main_loop_iter % OPP_fill_period != 0 || OPP_comm_iteration != 0))) {

        if (OPP_DBG) 
            opp_printf("opp_finalize_particle_move", "hole fill set [%s]", set->name);
        
        opp_profiler->start("Mv_holefill_core");
        particle_hole_fill_core(set);
        opp_profiler->end("Mv_holefill_core");
    }
    else if (OPP_fill_type == OPP_Sort_All || OPP_fill_type == OPP_Sort_Periodic) {
        
        if (OPP_DBG)
            opp_printf("opp_finalize_particle_move", "sort set [%s]", set->name);
        
        opp_profiler->start("Mv_sort_core");
        opp_particle_sort_core(set, false);
        opp_profiler->end("Mv_sort_core");
    }
    else if (OPP_fill_type == OPP_Shuffle_All || OPP_fill_type == OPP_Shuffle_Periodic) {
        
        if (OPP_DBG) 
            opp_printf("opp_finalize_particle_move", "shuffle set [%s]", set->name);
        
        opp_profiler->start("Mv_shuffle_core");
        opp_particle_sort_core(set, true); // true will shuffle the particles
        opp_profiler->end("Mv_shuffle_core");
    }
    else {
        opp_abort("OPP_fill_type is undefined");
    }

    set->size -= set->particle_remove_count;
}

//****************************************
void opp_particle_sort_core(opp_set set, bool shuffle)
{ 
    if (OPP_DBG) printf("\topp_particle_sort set [%s]\n", set->name);
    
    // opp_profiler->start("PartSort");

    std::vector<OPP_INT> from_indices = sort_iota_by_key<OPP_INT>(
                                            (OPP_INT*)set->mesh_relation_dat->data, set->size);

    for (opp_dat& dat : *(set->particle_dats)) { 

        char *new_data = (char *)opp_host_malloc(set->set_capacity * dat->size);
        char *old_data = dat->data;
        
        for (size_t j = 0; j < (size_t)set->size; j++) {
            memcpy(new_data + j * dat->size, old_data + from_indices[j] * dat->size, dat->size);
        }

        opp_host_free(dat->data);
        dat->data = new_data;
    }

    // opp_profiler->end("PartSort");
}

//****************************************
void particle_hole_fill_core(opp_set set)
{
    // opp_profiler->start("Core_HoleFill");

    OPP_INT* mesh_relation_data = (OPP_INT*)set->mesh_relation_dat->data;
    std::vector<std::pair<size_t, size_t>> swap_indices;  // Contains hole index and the index from back to swap
    swap_indices.reserve(set->particle_remove_count);     // Pre-allocate memory if the number of removals is known

        // // Idea: The last available element should be copied to the hole
        // // In the below scope we try to calculate the element to be swapped with the hole
        // {
        //     // set->particle_remove_count   // the particle count that should be removed
        //     int removed_count = 0;          // how many elements currently being removed
        //     int skip_count = 0;             // how many elements from the back is skipped ..
        //                                     // .. due to that element is also to be removed

        //     for (size_t j = 0; j < (size_t)set->size; j++)
        //     {
        //         // skip if the current index is not to be removed
        //         if (mesh_relation_data[j] != MAX_CELL_INDEX) 
        //             continue;

        //         // handle if the element from the back is also to be removed
        //         while ((set->size - removed_count - skip_count - 1 >= 0) && 
        //             (mesh_relation_data[set->size - removed_count - skip_count - 1] == MAX_CELL_INDEX))
        //         {
        //             skip_count++;
        //         }

        //         // check whether the holes are at the back!
        //         if ((set->size - removed_count - skip_count - 1 < 0) ||
        //             (j >= (size_t)(set->size - removed_count - skip_count - 1))) 
        //         {
        //             if (OPP_DBG) 
        //                 opp_printf("opp_finalize_particle_move_core", 
        //                 "Current Iteration index [%d] and replacement index %d; hence breaking [%s]", 
        //                 j, (set->size - removed_count - skip_count - 1), set->name);
        //             break;
        //         }

        //         swap_indices.push_back(std::make_pair(j, (size_t)(set->size - removed_count - skip_count - 1)));

        //         removed_count++;
        //     }
        // }
    // REMOVE ABOVE IF NEW VERSION IS PERFORMING WELL

    // Idea: The last available element should be copied to the hole
    {
        int removed_count = 0;  // How many elements have been removed
        int back_index = set->size - 1; // Calculate the initial index for back-swap target

        for (int j = 0; j < set->size; j++) { // TODO : Do we know only the remove indices? if so, we could reduce this operation
            
            if (mesh_relation_data[j] != MAX_CELL_INDEX) // Skip if the current index is not a hole (not to be removed)
                continue;

            // Find a valid element from the back to swap with
            while (back_index >= 0 && mesh_relation_data[back_index] == MAX_CELL_INDEX) {
                back_index--;
            }

            // Check if the hole is at or past the back index or whether all holes are filled
            if (back_index <= j || removed_count >= set->particle_remove_count) {
                if (OPP_DBG) {
                    opp_printf("opp_finalize_particle_move_core",
                            "Current Iteration index [%d] and replacement index %d; breaking [%s]",
                            j, back_index, set->name);
                }
                break;
            }

            swap_indices.emplace_back((size_t)j, (size_t)back_index); // Found a valid hole and a fill data 

            back_index--; // Move to the next available back element for the next swap
            removed_count++;
        }
    }

    // For all the dats, fill the holes using the swap_indices
    for (opp_dat& dat : *(set->particle_dats)) {
        for (const auto& x : swap_indices) {
            char* dat_removed_ptr = dat->data + (x.first * dat->size);
            char* dat_to_replace_ptr = dat->data + (x.second * dat->size);
            
            memcpy(dat_removed_ptr, dat_to_replace_ptr, dat->size); 

            // for (int l = 0; l < dat->size; l++) {
            //     dat_removed_ptr[l] = dat_to_replace_ptr[l];
            // }
        }
    }

    // opp_profiler->end("Core_HoleFill");
}

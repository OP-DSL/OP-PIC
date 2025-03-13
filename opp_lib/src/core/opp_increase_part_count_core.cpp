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

//****************************************
bool opp_increase_particle_count_core(opp_set set, const OPP_INT insert_count)
{
    if (OPP_DBG && !set->is_particle) {
        opp_abort("Cannot increase_particle_count of a non particle set @opp_increase_particle_count_core");
    }

    if (OPP_DBG) 
        opp_printf("opp_increase_particle_count_core", "set [%s] with size [%d]", set->name, insert_count);        

    if (insert_count <= 0) {
        return true;
    }

    const OPP_INT new_set_size = (set->size + insert_count);

    // if the new particle set size is less or equal to set capacity, then just set new sizes instead of resizing
    if (set->set_capacity >= new_set_size) {
        if (OPP_DBG) 
            opp_printf("opp_increase_particle_count_core", "set [%s] No need to realloc, new size[%d] capacity[%d]", 
                set->name, new_set_size, set->set_capacity);        
        
        set->size = new_set_size;
        set->diff = insert_count;   
        return true;
    }

    // if the set needs resizing, then use alloc multiple to increase capacity to reduce regular resizing
    const size_t new_set_capacity = set->size + ((size_t)insert_count * OPP_part_alloc_mult);
    bool return_flag = true;

    if (OPP_DBG) //  || set->size != 0
        opp_printf("opp_increase_particle_count_core", "new_set_capacity %zu set_size %d num_dats_in_set %d", 
            new_set_capacity, set->size, set->particle_dats->size());

    // iterate over all the particle dats of that set and resize the arrays as necessary
    for (auto& dat : *(set->particle_dats)) {
        
        if (dat->data == NULL) {
            dat->data = (char *)opp_host_malloc((size_t)(new_set_capacity * dat->size));
            // opp_printf("opp_increase_particle_count_core", "opp_host_malloc name %s %p size %d", 
            //     dat->name, dat->data, (new_set_capacity * dat->size));

            if (dat->data == nullptr) {
                opp_printf("opp_increase_particle_count_core", "Error... alloc of dat name %s failed (size %zu)", 
                    dat->name, (size_t)(new_set_capacity * dat->size));
                return_flag = false;
            }
        }
        else {
            // char* old = dat->data;
            dat->data = (char *)opp_host_realloc(dat->data, (size_t)(new_set_capacity * dat->size));
            // opp_printf("opp_increase_particle_count_core", "realloc %p name %s %p size %d", 
            //     old, dat->name, dat->data, (new_set_capacity * dat->size));

            if (dat->data == nullptr) {
                opp_printf("opp_increase_particle_count_core", "Error... realloc of dat name %s failed (size %zu)", 
                    dat->name, (size_t)(new_set_capacity * dat->size));
                return_flag = false;
            }
        }

        if (dat->is_cell_index && (dat->data != nullptr)) {
            OPP_INT* mesh_rel_array = (OPP_INT *)dat->data;
            // #pragma code_align 32
            for (size_t i = (size_t)set->size; i < new_set_capacity; i++)
                mesh_rel_array[i] = MAX_CELL_INDEX;
        }
        // Note : The remainder in array from set size to set capacity will be garbage, except in mesh_relation_dat
    }
    
    set->size         = new_set_size;
    set->set_capacity = new_set_capacity;
    set->diff         = insert_count;

    return return_flag;
}

//****************************************
// Not quite necessary - 
// If opp_reset_num_particles_to_insert isn't called, set->diff will not reset and opp_par_loopp_inject__ will 
// loop only on the injected particles; However it will reset upon opp_increase_particle_count call
void opp_reset_num_particles_to_insert_core(opp_set set)
{
    set->diff = 0;
}

//****************************************
bool opp_inc_part_count_with_distribution_core(opp_set set, OPP_INT insert_count, opp_dat part_dist)
{
    opp_profiler->start("IncPartCountWithDistribution");

    if (!opp_increase_particle_count_core(set, insert_count)) {
        opp_profiler->end("IncPartCountWithDistribution");
        return false;
    }

    OPP_INT* mesh_connectivity = (OPP_INT *)set->mesh_relation_dat->data;
    const OPP_INT* distribution = (OPP_INT *)part_dist->data;

    const OPP_INT start = (set->size - set->diff);
    OPP_INT j = 0;

    for (OPP_INT i = 0; i < set->diff; i++) {
        if (i >= distribution[j]) 
            j++; // check whether it is j or j-1    
        mesh_connectivity[start + i] = j;
    } 

    opp_profiler->end("IncPartCountWithDistribution");

    return true;
}

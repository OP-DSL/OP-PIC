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

#include "opp_mpi.h"

void opp_part_pack(opp_set set);
void opp_part_unpack(opp_set set);

//****************************************
void opp_init_particle_move(opp_set set, int nargs, opp_arg *args)
{ 
    opp_init_particle_move_core(set);

    opp_move_part_indices.clear();
    opp_move_part_indices.reserve(20000);

    if (OPP_comm_iteration == 0) {
        OPP_iter_start = 0;
        OPP_iter_end   = set->size;          
        OPP_part_comm_count_per_iter = 0; 
    }
    else {
        // need to change the arg data since particle communication could change the pointer in realloc dat->data
        for (int i = 0; i < nargs; i++) {
            if (args[i].argtype == OPP_ARG_DAT && args[i].dat->set->is_particle) {          
                if (OPP_DBG) 
                    opp_printf("opp_init_particle_move", "dat %s", args[i].dat->name);
                args[i].data = args[i].dat->data;
            }
        }
    }

    if (OPP_DBG) opp_printf("opp_init_particle_move", "comm_iter=%d start=%d end=%d", 
                OPP_comm_iteration, OPP_iter_start, OPP_iter_end);

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 
}

//****************************************
bool opp_finalize_particle_move(opp_set set)
{ 
    if (OPP_DBG) 
        opp_printf("opp_finalize_particle_move", "Start particle set [%s]", set->name);

    opp_profiler->start("Mv_Finalize");

    OPP_part_comm_count_per_iter += (int)opp_move_part_indices.size();

    opp_process_marked_particles(set); 

    opp_part_pack(set);
    
    // send the counts and send the particles  
    opp_part_exchange(set);  

    // Can fill the holes here, since the communicated particles will be added at the end
    opp_finalize_particle_move_core(set);

    if (opp_part_check_all_done(set)) {

        if (OPP_max_comm_iteration < OPP_comm_iteration) {
            OPP_max_comm_iteration = OPP_comm_iteration;
        }

        OPP_comm_iteration = 0; // reset for the next par loop
        
        opp_profiler->end("Mv_Finalize");

        return false; // all mpi ranks do not have anything to communicate to any rank
    }
        
    opp_part_wait_all(set); // wait till all the particles are communicated
    
    // increase the particle count if required and unpack the communicated particle buffer 
    // in to separate particle dats
    opp_part_unpack(set);

    OPP_iter_start = set->size - set->diff;
    OPP_iter_end   = set->size;  

    OPP_comm_iteration++;  

    opp_profiler->end("Mv_Finalize");

    return true;
}

//****************************************
void opp_particle_sort(opp_set set) // unused
{ 
    opp_particle_sort_core(set);
}

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

//****************************************
void opp_increase_particle_count(opp_set part_set, const int insert_count)
{
    if (!opp_increase_particle_count_core(part_set, insert_count)) {
        opp_printf("opp_increase_particle_count", 
            "Error: opp_increase_particle_count_core failed for particle set [%s]", part_set->name);
        opp_abort("opp_increase_particle_count");        
    }
}

//****************************************
void opp_inc_part_count_with_distribution(opp_set part_set, int insert_count, opp_dat part_dist, bool calc_new)
{
    if (OPP_DBG) opp_printf("opp_inc_part_count_with_distribution", "insert_count [%d]", insert_count);

    if (!opp_inc_part_count_with_distribution_core(part_set, insert_count, part_dist)) {
        opp_printf("opp_inc_part_count_with_distribution", 
            "Error: opp_inc_part_count_with_distribution_core failed for particle set [%s]", part_set->name);
        opp_abort("opp_inc_part_count_with_distribution_core");        
    }
}
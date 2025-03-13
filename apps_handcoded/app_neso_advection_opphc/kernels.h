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

//*********************************************
// USER WRITTEN CODE
//*********************************************

#include "opp_lib.h"
#include <cmath>

//*************************************************************************************************
inline void update_pos_kernel(
    OPP_REAL* part_vel, 
    OPP_REAL* part_pos
)
{
    for (int dm = 0; dm < DIM; dm++) {
        
        part_pos[dm] += part_vel[dm] * CONST_dt; // s1 = s0 + ut
        
        // correct for periodic boundary conditions
        const OPP_INT n_extent_offset_int = std::abs(part_pos[dm]) + 2.0;
        const OPP_REAL temp_pos = part_pos[dm] + n_extent_offset_int * CONST_extents[dm];
        part_pos[dm] = std::fmod(temp_pos, CONST_extents[dm]);
    }
}

//*************************************************************************************************
inline void is_point_in_current_cell_kernel(
    OPP_INT* part_cid, 
    const OPP_REAL* part_pos, 
    const OPP_REAL* cell_pos_ll, 
    const OPP_INT* cell_cell_map)
{
    // check for x direction movement
    const OPP_REAL part_pos_x = part_pos[Dim::x];
    if (part_pos_x < cell_pos_ll[Dim::x]) {
        part_cid[0] = cell_cell_map[CellMap::xd_y];

        OPP_PARTICLE_NEED_MOVE; return;
    }
    if (part_pos_x > (cell_pos_ll[Dim::x] + CONST_cell_width)) {
        part_cid[0] = cell_cell_map[CellMap::xu_y];
        OPP_PARTICLE_NEED_MOVE; return;
    }

    // check for y direction movement
    const OPP_REAL part_pos_y = part_pos[Dim::y];
    if (part_pos_y < cell_pos_ll[Dim::y]) {
        part_cid[0] = cell_cell_map[CellMap::x_yd];
        OPP_PARTICLE_NEED_MOVE; return;
    }
    if (part_pos_y > (cell_pos_ll[Dim::y] + CONST_cell_width)) {
        part_cid[0] = cell_cell_map[CellMap::x_yu];
        OPP_PARTICLE_NEED_MOVE; return;
    }

    OPP_PARTICLE_MOVE_DONE;
}

//*************************************************************************************************
inline void push_particles_kernel(
    OPP_INT* part_cid, 
    OPP_REAL* part_vel, 
    OPP_REAL* part_pos, 
    const OPP_REAL* cell_pos_ll, 
    const OPP_INT* cell_cell_map)
{
    if (OPP_DO_ONCE) {
        update_pos_kernel(
            part_vel, 
            part_pos);
    }

    is_point_in_current_cell_kernel(
        part_cid,
        part_pos, 
        cell_pos_ll, 
        cell_cell_map);
}

//*************************************************************************************************
inline void verify_kernel(
        const OPP_INT* part_cid, 
        const OPP_REAL* part_pos,
        const OPP_INT* cell_global_idx,
        OPP_INT* incorrect_part_count)
{
    // get the cell boundaries for the current cell_index - using global cell index 
    int ix = -1, iy = -1;
    RANK_TO_INDEX((*cell_global_idx), ix, iy, CONST_ndimcells[Dim::x]); 
    
    if (ix < 0 || iy < 0)
    {
        opp_printf("VERIFY", "Incorrect ix[%d] iy[%d] for global cell[%d] nx[%d]",
            ix, iy, (*cell_global_idx), CONST_ndimcells[Dim::x]);
        (*incorrect_part_count)++;
        return;
    }
    
    // get the boundaries of that cell
    const OPP_REAL boundary_ll[DIM] = { (ix * CONST_cell_width), (iy * CONST_cell_width) };
 
    // check whether the current particle is within those boundaries or not!
    const OPP_REAL part_pos_x = part_pos[Dim::x];
    if (part_pos_x < boundary_ll[Dim::x] ||
            part_pos_x > (boundary_ll[Dim::x] + CONST_cell_width)) {
        
        (*incorrect_part_count)++;
        return;
    }

    const OPP_REAL part_pos_y = part_pos[Dim::y];
    if (part_pos_y < boundary_ll[Dim::y] ||
            part_pos_y > (boundary_ll[Dim::y] + CONST_cell_width)) {
        
        (*incorrect_part_count)++;
        return;
    }
}

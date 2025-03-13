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

// //*************************************************************************************************
// inline void update_pos_kernel(const OPP_REAL* p_vel, OPP_REAL* p_pos, OPP_INT* p_mdir)
// {
//     for (int dm = 0; dm < DIM; dm++) {
        
//         const OPP_REAL offset = p_vel[dm] * CONST_dt[0];
//         p_pos[dm] += offset; // s1 = s0 + ut
        
//         // correct for periodic boundary conditions
//         const OPP_INT n_extent_offset_int = std::abs(p_pos[dm]) + 2.0;
//         const OPP_REAL temp_pos = p_pos[dm] + n_extent_offset_int * CONST_extents[dm];
//         p_pos[dm] = std::fmod(temp_pos, CONST_extents[dm]);

//         p_mdir[dm] = (offset > 0) ? 1 : -1;
//     }
// }

// //*************************************************************************************************
// inline void move_kernel(const OPP_REAL* p_pos, OPP_INT* p_mdir, const OPP_REAL* c_pos_ll)
// {
//     // check for x direction movement
//     const OPP_REAL p_pos_x_diff = (p_pos[Dim::x] - c_pos_ll[Dim::x]);
//     if ((p_pos_x_diff >= 0.0) && (p_pos_x_diff <= CONST_cell_width[0])) {
//         p_mdir[Dim::x] = 0; // within cell in x direction
//     }
//     else if (p_mdir[Dim::x] > 0) {
//         opp_p2c[0] = opp_c2c[CellMap::xu_y];
//         OPP_PARTICLE_NEED_MOVE; return;
//     }
//     else if (p_mdir[Dim::x] < 0) {
//         opp_p2c[0] = opp_c2c[CellMap::xd_y];
//         OPP_PARTICLE_NEED_MOVE; return;
//     }

//     // check for y direction movement
//     const OPP_REAL p_pos_y_diff = (p_pos[Dim::y] - c_pos_ll[Dim::y]);
//     if ((p_pos_y_diff >= 0.0) && (p_pos_y_diff <= CONST_cell_width[0])) { 
//         p_mdir[Dim::y] = 0; // within cell in y direction
//     }
//     else if (p_mdir[Dim::y] > 0) {
//         opp_p2c[0] = opp_c2c[CellMap::x_yu];
//         OPP_PARTICLE_NEED_MOVE; return;
//     }
//     else if (p_mdir[Dim::y] < 0) {
//         opp_p2c[0] = opp_c2c[CellMap::x_yd];
//         OPP_PARTICLE_NEED_MOVE; return;
//     }

//     OPP_PARTICLE_MOVE_DONE;
// }

//*************************************************************************************************
inline void update_pos_kernel(const OPP_REAL* part_vel, OPP_REAL* p_pos)
{
    for (int dm = 0; dm < DIM; dm++) {
        
        p_pos[dm] += part_vel[dm] * CONST_dt[0]; // s1 = s0 + ut
        
        // correct for periodic boundary conditions
        const OPP_INT n_extent_offset_int = std::abs(p_pos[dm]) + 2.0;
        const OPP_REAL temp_pos = p_pos[dm] + n_extent_offset_int * CONST_extents[dm];
        p_pos[dm] = std::fmod(temp_pos, CONST_extents[dm]);
    }
}

//*************************************************************************************************
inline void move_kernel(const OPP_REAL* p_pos, const OPP_REAL* c_pos_ll)
{
    // check for x direction movement
    const OPP_REAL p_pos_x = p_pos[Dim::x];
    if (p_pos_x < c_pos_ll[Dim::x]) {
        opp_p2c[0] = opp_c2c[CellMap::xd_y];
        OPP_PARTICLE_NEED_MOVE; return;
    }
    if (p_pos_x > (c_pos_ll[Dim::x] + CONST_cell_width[0])) {
        opp_p2c[0] = opp_c2c[CellMap::xu_y];
        OPP_PARTICLE_NEED_MOVE; return;
    }

    // check for y direction movement
    const OPP_REAL p_pos_y = p_pos[Dim::y];
    if (p_pos_y < c_pos_ll[Dim::y]) {
        opp_p2c[0] = opp_c2c[CellMap::x_yd];
        OPP_PARTICLE_NEED_MOVE; return;
    }
    if (p_pos_y > (c_pos_ll[Dim::y] + CONST_cell_width[0])) {
        opp_p2c[0] = opp_c2c[CellMap::x_yu];
        OPP_PARTICLE_NEED_MOVE; return;
    }

    OPP_PARTICLE_MOVE_DONE;
}

//*************************************************************************************************
inline void verify_kernel(
        const OPP_REAL* p_pos,
        const OPP_INT* c_gbl_idx,
        OPP_INT* incorrect_count)
{
    // get the cell boundaries for the current cell_index - using global cell index 
    int ix = -1, iy = -1;
    RANK_TO_INDEX((*c_gbl_idx), ix, iy, CONST_ndimcells[Dim::x]); 
    
    if (ix < 0 || iy < 0)
    {
        // opp_printf("VERIFY", "Incorrect ix[%d] iy[%d] for global cell[%d] nx[%d]",
        //     ix, iy, (*c_gbl_idx), CONST_ndimcells[Dim::x]);
        (*incorrect_count)++;
        return;
    }
    
    // get the boundaries of that cell
    const OPP_REAL boundary_ll[DIM] = { (ix * CONST_cell_width[0]), (iy * CONST_cell_width[0]) };
 
    // check whether the current particle is within those boundaries or not!
    const OPP_REAL p_pos_x = p_pos[Dim::x];
    if (p_pos_x < boundary_ll[Dim::x] ||
            p_pos_x > (boundary_ll[Dim::x] + CONST_cell_width[0])) {
        
        (*incorrect_count)++;
        return;
    }

    const OPP_REAL p_pos_y = p_pos[Dim::y];
    if (p_pos_y < boundary_ll[Dim::y] ||
            p_pos_y > (boundary_ll[Dim::y] + CONST_cell_width[0])) {
        
        (*incorrect_count)++;
        return;
    }
}

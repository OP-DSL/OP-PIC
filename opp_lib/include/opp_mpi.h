
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

#include <opp_mpi_core.h>

/*******************************************************************************/
void opp_halo_create();
void opp_halo_destroy();

int opp_mpi_halo_exchanges(opp_set set, int nargs, opp_arg *args);
void opp_mpi_halo_exchange(opp_arg *arg, int exec_flag);
void opp_mpi_halo_wait_all(int nargs, opp_arg *args);
/*******************************************************************************/

void opp_exchange_double_indirect_reductions(opp_dat dat, opp_reduc_comm reduc_comm);
void opp_complete_double_indirect_reductions(opp_dat dat);

inline void opp_mpi_reduce(opp_arg *args, double *data) 
{
    opp_mpi_reduce_double(args, data);
}

inline void opp_mpi_reduce(opp_arg *args, int *data) 
{
    opp_mpi_reduce_int(args, data);
}

extern char opp_move_status_flag;
extern bool opp_move_hop_iter_one_flag;
extern OPP_INT* opp_p2c;
extern OPP_INT* opp_c2c;

#define OPP_PARTICLE_MOVE_DONE { opp_move_status_flag = OPP_MOVE_DONE; }
#define OPP_PARTICLE_NEED_MOVE { opp_move_status_flag = OPP_NEED_MOVE; }
#define OPP_PARTICLE_NEED_REMOVE { opp_move_status_flag = OPP_NEED_REMOVE; }
#define OPP_DO_ONCE (opp_move_hop_iter_one_flag)
#define OPP_MOVE_RESET_FLAGS { opp_move_status_flag = OPP_MOVE_DONE; opp_move_hop_iter_one_flag = true; }

//*************************************************************************************************
inline bool opp_check_part_move_status(const OPP_INT map0idx, const OPP_INT particle_index, int& remove_count) 
{
    opp_move_hop_iter_one_flag = false;

    if (opp_move_status_flag == OPP_MOVE_DONE) {
        return false;
    }
    else if (opp_move_status_flag == OPP_NEED_REMOVE) {
        remove_count += 1;
        OPP_mesh_relation_data[particle_index] = MAX_CELL_INDEX;

        return false;
    }
    else if (map0idx >= OPP_part_cells_set_size) {
        // map0idx cell is not owned by the current mpi rank (it is in the import exec halo region), need to communicate
        opp_move_part_indices.push_back(particle_index);
        return false;
    }

    // map0idx is an own cell and opp_move_status_flag == OPP_NEED_MOVE
    return true;
}

//*******************************************************************************
// returns true, if the current particle needs to be removed from the rank
inline bool opp_part_checkForGlobalMove_util(opp_set set, const opp_point& point, const int partIndex, int& cellIdx, 
                const size_t structCellIdx) 
{             
    if (structCellIdx == MAX_CELL_INDEX || structCellIdx == OPP_OUT_OF_SAVED_DOMAIN) { // This happens when point is out of the unstructured mesh
        if (OPP_DBG)
            opp_printf("opp_part_checkForGlobalMove", 
            "Remove %d [Struct cell index invalid - strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                partIndex, structCellIdx, point.x, point.y, point.z);
        cellIdx = MAX_CELL_INDEX;
        set->particle_remove_count++;
        return true;
    }

    const int structCellRank = cellMapper->findClosestCellRank(structCellIdx);

    // Check whether the paticles need global moving, if yes start global moving process, 
    // if no, move to the closest local cell
    if (structCellRank != OPP_rank) {

        if (structCellRank == MAX_CELL_INDEX) {
            if (OPP_DBG)
                opp_printf("opp_part_checkForGlobalMove", 
                "Remove %d [Rank invalid - cellRank:%d loclCellIdx:%zu strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                    partIndex, structCellRank, cellMapper->findClosestCellIndex(structCellIdx), structCellIdx, 
                    point.x, point.y, point.z);
            set->particle_remove_count++;
            cellIdx = MAX_CELL_INDEX;
            return true;
        }

        // Due to renumbering local cell indices will be different to global, hence do global comm with global indices
        const size_t globalCellIndex = cellMapper->findClosestCellIndex(structCellIdx);

        if (globalCellIndex == MAX_CELL_INDEX) {
            if (OPP_DBG)
                opp_printf("opp_part_checkForGlobalMove", 
                "Remove %d [CellIdx invalid - cellRank:%d loclCellIdx:%zu strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                    partIndex, structCellRank, globalCellIndex, structCellIdx, point.x, point.y, point.z);
            set->particle_remove_count++;
            cellIdx = MAX_CELL_INDEX;
            return true;
        }

        // if the new rank is not the current rank, mark the particle to be sent via global comm
        globalMover->markParticleToMove(set, partIndex, structCellRank, globalCellIndex);

        // if (OPP_DBG)
        //     opp_printf("opp_part_checkForGlobalMove", "Mark part %d [Move to rank %d gblCellIdx %d]", 
        //         partIndex, structCellRank, globalCellIndex);
        set->particle_remove_count++;
        cellIdx = MAX_CELL_INDEX;
        return true;
    }
    else {
        
        // Due to renumbering local cell indices will be different to global, hence do global comm with global indices
        cellIdx = cellMapper->findClosestCellIndex(structCellIdx);

        if (OPP_DBG && (cellIdx < 0 || cellIdx >= set->cells_set->size)) {
            opp_printf("opp_part_checkForGlobalMove", 
                "Error... Particle %d assigned to current rank but invalid cell index %d [strCellIdx:%zu]", 
                    partIndex, cellIdx, structCellIdx);
            opp_abort("opp_part_checkForGlobalMove Error... Particle assigned to current rank but invalid cell index");
        }
    }
                
    return false;
}

//*******************************************************************************
// returns true, if the current particle needs to be removed from the rank
inline bool opp_part_checkForGlobalMove2D(opp_set set, const opp_point& point, const int partIndex, int& cellIdx) 
{
    const size_t structCellIdx = cellMapper->findStructuredCellIndex2D(point);
    return opp_part_checkForGlobalMove_util(set, point, partIndex, cellIdx, structCellIdx);
}

//*******************************************************************************
// returns true, if the current particle needs to be removed from the rank
inline bool opp_part_checkForGlobalMove3D(opp_set set, const opp_point& point, const int partIndex, int& cellIdx) 
{
    const size_t structCellIdx = cellMapper->findStructuredCellIndex3D(point);
    return opp_part_checkForGlobalMove_util(set, point, partIndex, cellIdx, structCellIdx);
}
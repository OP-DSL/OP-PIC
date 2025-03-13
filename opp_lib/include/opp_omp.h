
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

#include <opp_lib.h>
#include <omp.h>

#ifdef USE_MPI
    #include <opp_mpi_core.h>
#endif

struct opp_tmp_gbl_move_info
{
    opp_set set = nullptr;
    OPP_INT part_idx = -1;
    int struct_cell_rank = MAX_CELL_INDEX;
    OPP_INT gbl_cell_idx = MAX_CELL_INDEX;

    opp_tmp_gbl_move_info(opp_set _set, OPP_INT _part_idx, int _struct_cell_rank, OPP_INT _gbl_cell_idx) : 
        set(_set), part_idx(_part_idx), struct_cell_rank(_struct_cell_rank), gbl_cell_idx(_gbl_cell_idx) {}
};

extern int opp_nthreads;
extern std::vector<int> part_remove_count_per_thr;
extern std::vector<std::vector<int>> move_part_indices_per_thr;
extern std::vector<std::vector<opp_tmp_gbl_move_info>> gbl_move_indices_per_thr;

// TODO : May need to write with array capcity and realloc always if partticle dats are used
template <class T> 
void opp_create_thread_level_data(opp_arg arg, T init_value)
{
    opp_dat dat = arg.dat;
    opp_set set = dat->set;
    const int nthreads = omp_get_max_threads();

    if (dat->set->is_particle) {
        std::cerr << "Cannot create thread level data for particle dat [" << dat->name << 
            "] (dat in a dynamic set)" << std::endl;
        exit(-1);
    }

    if (OPP_DBG) printf("opp_create_thread_level_data template[%d]\n", nthreads);

    const size_t dat_set_size = (size_t)(set->size) + (size_t)(set->exec_size) + (size_t)(set->nonexec_size);

    if (dat->thread_data->size() <= 0) {
        dat->thread_data->push_back(dat->data);

        for (int thr = 1; thr < nthreads; thr++) {
            char* thr_data = (char *)opp_host_malloc((size_t)dat->size * dat_set_size * sizeof(char));;
            dat->thread_data->push_back(thr_data);
        }
    }

    if ((int)dat->thread_data->size() != nthreads) {
        std::cerr << "opp_create_thread_level_data dat [" << dat->name << 
            "] thread_data not properly created [(int)dat->thread_data.size():" << (int)dat->thread_data->size() << 
            " nthreads:" << nthreads << std::endl;
        return;
    }

    #pragma omp parallel for
    for (int thr = 1; thr < nthreads; thr++) {
        std::fill_n((T*)(dat->thread_data->at(thr)), (dat->dim * dat_set_size), init_value);
    }
}

// TODO : May need to write with array capcity and realloc always if partticle dats are used
template <class T> 
void opp_reduce_thread_level_data(opp_arg arg)
{
    opp_dat dat = arg.dat;
    opp_set set = dat->set;
    std::vector<char *>& thread_data = *(dat->thread_data);

    const int nthreads = omp_get_max_threads();

    if (OPP_DBG) printf("opp_reduce_thread_level_data dat [%s] nthreads [%d]\n", dat->name, nthreads);

    if (set->size > 0) {
        const size_t dat_set_size = (size_t)(set->size) + (size_t)(set->exec_size) + (size_t)(set->nonexec_size);

        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++) {
            const int start  = ((dat->dim * dat_set_size)* thr)/nthreads;
            const int finish = ((dat->dim * dat_set_size)*(thr+1))/nthreads;
            
            // if (OPP_DBG) printf("opp_reduce_thread_level_data THREAD [%d] %d %d\n", thr, start, finish);
            for (int n = start; n < finish; n++)
            {
                for (int array_num = 1; array_num < nthreads; array_num++)
                {
                    T* td = (T*)thread_data[array_num];

                    switch (arg.acc) {
                        case OPP_INC:
                            ((T*)dat->data)[n] += td[n];
                            break;
                        default:
                            std::cerr << "opp_reduce_thread_level_data dat [" << dat->name << "] acc [" << 
                                            (int)arg.acc << "] not implemented" << std::endl;
                    }
                }
            }
            // if (OPP_DBG) printf("opp_reduce_thread_level_data THREAD [%d] END\n", thr);
        }
    }
}

/*******************************************************************************/
void opp_halo_create();
void opp_halo_destroy();

/*******************************************************************************/
void print_dat_to_txtfile_mpi(opp_dat dat, const char *file_name);
void opp_mpi_print_dat_to_txtfile(opp_dat dat, const char *file_name);

/*******************************************************************************/
inline void opp_mpi_reduce(opp_arg *args, double *data) 
{
#ifdef USE_MPI
    opp_mpi_reduce_double(args, data);
#else
    (void)args;
    (void)data;
#endif
}

inline void opp_mpi_reduce(opp_arg *args, int *data) 
{
#ifdef USE_MPI
    opp_mpi_reduce_int(args, data);
#else
    (void)args;
    (void)data;
#endif
}

#define OPP_PARTICLE_MOVE_DONE {  }
#define OPP_PARTICLE_NEED_MOVE {  }
#define OPP_PARTICLE_NEED_REMOVE {  }
#define OPP_DO_ONCE ( true )
#define OPP_MOVE_RESET_FLAGS {  }

//****************************************
inline bool opp_check_part_move_status(const char move_flag, bool& iter_one_flag, const int map0idx,
    const int particle_index, int thread)
{
    iter_one_flag = false;

    if (move_flag == OPP_MOVE_DONE) {
        return false;
    }
    else if (move_flag == OPP_NEED_REMOVE) {
        part_remove_count_per_thr[thread] += 1;
        OPP_mesh_relation_data[particle_index] = MAX_CELL_INDEX;

        return false;
    }
#ifdef USE_MPI
    else if (map0idx >= OPP_part_cells_set_size) {
        // map0idx cell is not owned by the current mpi rank (it is in the import exec halo region), need to communicate
        move_part_indices_per_thr[thread].push_back(particle_index);
        return false;
    }
#endif

    return true;
}

//*******************************************************************************
// returns true, if the current particle needs to be removed from the rank
inline bool opp_part_checkForGlobalMove_util(opp_set set, const opp_point& point, const int partIndex, 
            int& cellIdx, const size_t structCellIdx, const int thread) {

    if (structCellIdx == MAX_CELL_INDEX || structCellIdx == OPP_OUT_OF_SAVED_DOMAIN) { // This happens when point is out of the unstructured mesh
        if (false && OPP_DBG)
            opp_printf("opp_part_checkForGlobalMove", 
            "Remove %d [Struct cell index invalid - strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                partIndex, structCellIdx, point.x, point.y, point.z);

        cellIdx = MAX_CELL_INDEX;
        return true;
    }

#ifdef USE_MPI  
    const int structCellRank = cellMapper->findClosestCellRank(structCellIdx);

    // Check whether the paticles need global moving, if yes start global moving process, 
    // if no, move to the closest local cell
    if (structCellRank != OPP_rank) {

        if (structCellRank == MAX_CELL_INDEX) {
            if (false && OPP_DBG)
                opp_printf("opp_part_checkForGlobalMove", 
                "Remove %d [Rank invalid - strCellRank:%d loclCellIdx:%zu strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                    partIndex, structCellRank, cellMapper->findClosestCellIndex(structCellIdx), structCellIdx, 
                    point.x, point.y, point.z);
            cellIdx = MAX_CELL_INDEX;
            return true;
        }

        // Due to renumbering local cell indices will be different to global, hence do global comm with global indices
        const size_t globalCellIndex = cellMapper->findClosestCellIndex(structCellIdx);

        if (globalCellIndex == MAX_CELL_INDEX) {
            if (false && OPP_DBG)
                opp_printf("opp_part_checkForGlobalMove", 
                "Remove %d [CellIdx invalid - strCellRank:%d loclCellIdx:%zu strCellIdx:%zu] [%2.16lE, %2.16lE, %2.16lE]", 
                    partIndex, structCellRank, globalCellIndex, structCellIdx, point.x, point.y, point.z);
            cellIdx = MAX_CELL_INDEX;
            return true;
        }

        // if the new rank is not the current rank, mark the particle to be sent via global comm
        // globalMover->markParticleToMove(set, partIndex, structCellRank, globalCellIndex);
        gbl_move_indices_per_thr[thread].push_back(
                opp_tmp_gbl_move_info(set, partIndex, structCellRank, globalCellIndex));
        // if (OPP_DBG)
        //     opp_printf("opp_part_checkForGlobalMove", "Mark part %d [Move to rank %d gblCellIdx %d]", 
        //         partIndex, structCellRank, globalCellIndex);
        cellIdx = MAX_CELL_INDEX;
        return true;
    }
    else
#endif 
    {  
        // Due to renumbering local cell indices will be different to global, hence do global comm with global indices
        cellIdx = cellMapper->findClosestCellIndex(structCellIdx);

        if (false && OPP_DBG && (cellIdx < 0 || cellIdx >= set->cells_set->size)) {
            opp_printf("opp_part_checkForGlobalMove", 
                "Error... Particle %d assigned to current rank but invalid cell index %d [strCellIdx:%zu]", 
                    partIndex, cellIdx, structCellIdx);
            opp_abort("opp_part_checkForGlobalMove Error... Particle assigned to current rank but invalid cell index");
        }

        if (cellIdx == MAX_CELL_INDEX) { // Particle is outside the mesh, need to remove
            return true;
        }
    }
                
    return false;
}

//*******************************************************************************
// returns true, if the current particle needs to be removed from the rank
inline bool opp_part_checkForGlobalMove2D(opp_set set, const opp_point& point, const int partIndex, 
        int& cellIdx, const int thread) 
{
    const size_t structCellIdx = cellMapper->findStructuredCellIndex2D(point);
    return opp_part_checkForGlobalMove_util(set, point, partIndex, cellIdx, structCellIdx, thread);
}

//*******************************************************************************
// returns true, if the current particle needs to be removed from the rank
inline bool opp_part_checkForGlobalMove3D(opp_set set, const opp_point& point, const int partIndex, 
        int& cellIdx, const int thread) 
{
    const size_t structCellIdx = cellMapper->findStructuredCellIndex3D(point);
    return opp_part_checkForGlobalMove_util(set, point, partIndex, cellIdx, structCellIdx, thread);
}

//*******************************************************************************
// gathers all per thread global move information into the global mover for communication
inline void opp_gather_gbl_move_indices()
{
#ifdef USE_MPI
    for (auto& per_thread_info : gbl_move_indices_per_thr) {
        for (auto& p_info : per_thread_info) {
            globalMover->markParticleToMove(p_info.set, p_info.part_idx, 
                p_info.struct_cell_rank, p_info.gbl_cell_idx);
        }
    }
#endif
}

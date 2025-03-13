
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

#include <opp_defs.h>
#include <opp_util.h>
#include <opp_params.h>
#include <opp_profiler.h>

//*************************************************************************************************
void opp_init_core(int argc, char **argv);
void opp_exit_core();
void opp_abort(std::string s = "");

void opp_set_args_core(char *argv);

opp_set opp_decl_set_core(int size, char const *name);
opp_set opp_decl_particle_set_core(int size, char const *name, opp_set cells_set);
opp_set opp_decl_particle_set_core(char const *name, opp_set cells_set);

opp_map opp_decl_map_core(opp_set from, opp_set to, int dim, int *imap, char const *name);

opp_dat opp_decl_dat_core(opp_set set, int dim, char const *type, int size, char *data, char const *name);

void opp_decl_const_impl(int dim, int size, char* data, const char* name);

opp_arg opp_arg_dat_core(opp_dat dat, int idx, opp_map map, int dim, const char *typ, opp_dat p2c_map, opp_access acc);
opp_arg opp_arg_dat_core(opp_map data_map, int idx, opp_map map, opp_dat p2c_map, opp_access acc);

template <typename T>
opp_arg opp_arg_gbl_core(T *data, int dim, char const *typ, opp_access acc)
{
    opp_arg arg;
    arg.argtype     = OPP_ARG_GBL;

    arg.dat         = NULL;
    arg.map         = NULL;
    arg.dim         = dim;
    arg.idx         = -1;
    arg.size        = dim * sizeof(T);
    arg.data        = (char*)data;
    arg.map_data    = NULL;
    arg.type        = typ;
    arg.acc         = acc;
    arg.opt         = 1;

    return arg;
}

bool opp_increase_particle_count_core(opp_set particles_set, const int num_particles_to_insert);

bool opp_inc_part_count_with_distribution_core(opp_set particles_set, int num_particles_to_insert, opp_dat part_dist);

void opp_reset_num_particles_to_insert_core(opp_set set);

void opp_init_particle_move_core(opp_set set);

void opp_finalize_particle_move_core(opp_set set);

void opp_particle_sort_core(opp_set set, bool shuffle = false);

void opp_reset_dat_impl(opp_dat dat, char* val, opp_reset reset = OPP_Reset_All);

void opp_print_dat_to_txtfile_core(opp_dat dat, const char *file_name_prefix, const char *file_name_suffix);
void opp_print_map_to_txtfile_core(opp_map map, const char *file_name_prefix, const char *file_name_suffix);
void opp_dump_dat_core(opp_dat data);

/*******************************************************************************/
void opp_set_dirtybit(int nargs, opp_arg *args);
void opp_set_dirtybit_device(int nargs, opp_arg *args);
void opp_set_dirtybit_grouped(int nargs, opp_arg *args, DeviceType device);

//*******************************************************************************
inline void* opp_host_malloc(size_t size)
{
    return malloc(size);
}
inline void* opp_host_realloc(void* ptr, size_t new_size)
{
    return realloc(ptr, new_size);
}
inline void opp_host_free(void* ptr)
{
    free(ptr);
    ptr = nullptr;
}

//*******************************************************************************
void* opp_load_from_file_core(const char* file_name, int set_size, int dim, char const *type, int size);

template <typename T> 
inline void opp_reduce_dat_element(T* out_dat, const T* in_dat, int dim, opp_reduc_comm reduc_comm)
{
    for (int d = 0; d < dim; d++) {
        switch (reduc_comm) {
            case OPP_Reduc_SUM_Comm: 
                out_dat[d] += in_dat[d];           
                break;
            case OPP_Reduc_MAX_Comm: 
                out_dat[d] = MAX(out_dat[d], in_dat[d]);
                break;
            case OPP_Reduc_MIN_Comm: 
                out_dat[d] = MIN(out_dat[d], in_dat[d]);
                break;
            default:
                opp_printf("opp_reduce_dat_element", "Unhandled reduction type");
        }  
    }   
}

//*******************************************************************************
template <typename T>
inline void opp_write_array_to_file(const T* array, size_t size, const std::string& filename) {
    
    std::stringstream modified_file_name;
    modified_file_name << "files/" << filename << "_r" << OPP_rank << "_i" << OPP_main_loop_iter;
    std::ofstream outFile(modified_file_name.str());
    if (!outFile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    if constexpr (std::is_same<T, double>::value)
        outFile << std::setprecision(25);
    outFile << size << " 1 -- 0 0\n";
    for (size_t i = 0; i < size; ++i) {
        outFile << " " << array[i] << "\n";
    }
    outFile.close();
}

//*******************************************************************************

#ifdef USE_MPI
    #include <opp_mpi_core.h>
#endif

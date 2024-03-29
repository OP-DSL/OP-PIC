
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

#include <oppic_lib.h>
#include <omp.h>

#ifdef USE_MPI
    #include <opp_mpi_core.h>
#endif

extern std::vector<int> part_remove_count_per_thr;

// TODO : May need to write with array capcity and realloc always if partticle dats are used
template <class T> 
void opp_create_thread_level_data(oppic_arg arg, T init_value)
{
    oppic_dat dat = arg.dat;
    int nthreads = omp_get_max_threads();

    if (dat->set->is_particle)
    {
        std::cerr << "Cannot create thread level data for particle dat [" << dat->name << "] (dat in a dynamic set)" << std::endl;
        exit(-1);
    }

    if (OP_DEBUG) printf("opp_create_thread_level_data template[%d]\n", nthreads);

    if (dat->thread_data->size() <= 0)
    {
        dat->thread_data->push_back(dat->data);

        for (int thr = 1; thr < nthreads; thr++)
        {
            char* thr_data = (char *)opp_host_malloc((size_t)dat->size * (size_t)(dat->set->size) * sizeof(char));;
            dat->thread_data->push_back(thr_data);
        }
    }

    if ((int)dat->thread_data->size() != nthreads)
    {
        std::cerr << "opp_create_thread_level_data dat [" << dat->name << "] thread_data not properly created [(int)dat->thread_data.size():" << (int)dat->thread_data->size() << " nthreads:" << nthreads << std::endl;
        return;
    }

    for (int thr = 1; thr < nthreads; thr++)
    {
        std::fill_n((T*)(dat->thread_data->at(thr)), (dat->dim * dat->set->size), init_value);
    }
}

// TODO : May need to write with array capcity and realloc always if partticle dats are used
template <class T> 
void opp_reduce_thread_level_data(oppic_arg arg)
{
    oppic_dat dat = arg.dat;
    oppic_set set = dat->set;
    std::vector<char *>& thread_data = *(dat->thread_data);

    int nthreads = omp_get_max_threads();

    if (OP_DEBUG) printf("opp_reduce_thread_level_data dat [%s] nthreads [%d]\n", dat->name, nthreads);

    if (set->size > 0) 
    {
        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++)
        {
            int start  = ((dat->dim * set->size)* thr)/nthreads;
            int finish = ((dat->dim * set->size)*(thr+1))/nthreads;
            
            for (int n = start; n < finish; n++)
            {
                for (int array_num = 1; array_num < nthreads; array_num++)
                {
                    T* td = (T*)thread_data[array_num];

                    switch (arg.acc)
                    {
                        case OP_INC:
                            ((T*)dat->data)[n] += td[n];
                            break;
                        default:
                            std::cerr << "opp_reduce_thread_level_data dat [" << dat->name << "] acc [" << (int)arg.acc << "] not implemented" << std::endl;
                    }
                }
            }
        }
    }
}

void oppic_finalize_particle_move_omp(oppic_set set);

bool opp_part_check_status_omp(opp_move_var& m, int map0idx, oppic_set set, 
    int particle_index, int& remove_count, int thread);
/*******************************************************************************/

void opp_halo_create();
void opp_halo_destroy();

/*******************************************************************************/

void print_dat_to_txtfile_mpi(op_dat dat, const char *file_name);
void opp_mpi_print_dat_to_txtfile(op_dat dat, const char *file_name);

/*******************************************************************************/

inline void opp_mpi_reduce(opp_arg *args, double *data) 
{
    (void)args;
    (void)data;
}

inline void opp_mpi_reduce(opp_arg *args, int *data) 
{
    (void)args;
    (void)data;
}
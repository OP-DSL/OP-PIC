
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

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>
// #include <device_launch_parameters.h>

#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <iostream>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/generate.h>
#include <thrust/random.h>

#ifdef USE_MPI
    #include <opp_mpi_core.h>
#endif

#define cutilSafeCall(err) __hipSafeCall(err, __FILE__, __LINE__)
#define cutilCheckMsg(msg) __cutilCheckMsg(msg, __FILE__, __LINE__)

#define OPP_GPU_THREADS_PER_BLOCK 32

extern int* opp_saved_mesh_relation_d;
extern size_t opp_saved_mesh_relation_size;
extern thrust::device_vector<int> cellIdx_dv;
extern thrust::device_vector<int> i_dv;
extern char *OPP_need_remove_flags_d;

extern int *OPP_move_indices_d;
extern int *OPP_move_count_d;
extern thrust::device_vector<int> OPP_thrust_move_indices_d;

//*************************************************************************************************

void __hipSafeCall(hipError_t err, const char *file, const int line);

void __cutilCheckMsg(const char *errorMessage, const char *file, const int line);


void oppic_hip_exit();

void cutilDeviceInit(int argc, char **argv);

void oppic_upload_dat(opp_dat dat);
void oppic_upload_map(opp_map map, bool create_new = false);

void oppic_download_dat(opp_dat dat);

void oppic_download_particle_set(opp_set particles_set, bool force_download = false);

void oppic_upload_particle_set(opp_set particles_set, bool realloc = false);

/*******************************************************************************/

void opp_halo_create();
void opp_halo_destroy();

/*******************************************************************************/

void opp_init_double_indirect_reductions_hip(int nargs, oppic_arg *args);
void opp_exchange_double_indirect_reductions_hip(int nargs, oppic_arg *args) ;
void opp_complete_double_indirect_reductions_hip(int nargs, oppic_arg *args);

/*******************************************************************************/


void print_dat_to_txtfile_mpi(op_dat dat, const char *file_name);
void opp_mpi_print_dat_to_txtfile(op_dat dat, const char *file_name);

void print_last_hip_error();

void oppic_cpHostToDevice(void **data_d, void **data_h, size_t copy_size, size_t alloc_size = 0, 
    bool create_new = false);

void oppic_create_device_arrays(oppic_dat dat, bool create_new = false);

void oppic_finalize_particle_move_hip(oppic_set set);

/*******************************************************************************/

void sort_dat_according_to_index_int(oppic_dat dat, const thrust::device_vector<int>& new_idx_dv, 
    int set_capacity, int size, bool hole_filling, int out_start_idx);
void sort_dat_according_to_index_double(oppic_dat dat, const thrust::device_vector<int>& new_idx_dv, 
    int set_capacity, int size, bool hole_filling, int out_start_idx);

/*
This function arranges the multi dimensional values in input array to output array according to the indices provided
    in_dat_dv - Input array with correct values
    out_dat_dv - Output array to have the sorted values
    new_idx_dv - indices of the input array to be arranged in the output array
    in_capacity - capacity of the input array (useful for multi dimensional dats)
    out_capacity - capacity of the output array (useful for multi dimensional dats)
    in_offset - start offset if the input array
    out_offset - start offset if the output array
    size - number of dat elements to arrange (usually this is the size of new_idx_dv)
    dimension - dimension of the dat
*/
template <class T> 
void copy_according_to_index(thrust::device_vector<T>* in_dat_dv, thrust::device_vector<T>* out_dat_dv, 
    const thrust::device_vector<int>& new_idx_dv, int in_capacity, int out_capacity, int in_offset, int out_offset, 
    int size, int dimension)
{
    switch (dimension)
    {
        case 1:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(in_dat_dv->begin() + in_offset)
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(out_dat_dv->begin() + out_offset)));
            break;
        case 2:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        in_dat_dv->begin() + in_offset, 
                        (in_dat_dv->begin() + in_offset + in_capacity)
                    )
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        out_dat_dv->begin() + out_offset, 
                        (out_dat_dv->begin() + out_offset + out_capacity))));
            break;
        case 3:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        in_dat_dv->begin() + in_offset, 
                        (in_dat_dv->begin() + in_offset + in_capacity), 
                        (in_dat_dv->begin() + in_offset + (2 * in_capacity))
                    )
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        out_dat_dv->begin() + out_offset, 
                        (out_dat_dv->begin() + out_offset + out_capacity), 
                        (out_dat_dv->begin() + out_offset + (2 * out_capacity)))));
            break;
        case 4:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        in_dat_dv->begin() + in_offset, 
                        (in_dat_dv->begin() + in_offset + in_capacity), 
                        (in_dat_dv->begin() + in_offset + (2 * in_capacity)), 
                        (in_dat_dv->begin() + in_offset + (3 * in_capacity))
                    )
                ), 
                new_idx_dv.begin()), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        out_dat_dv->begin() + out_offset, 
                        (out_dat_dv->begin() + out_offset + out_capacity), 
                        (out_dat_dv->begin() + out_offset + (2 * out_capacity)), 
                        (out_dat_dv->begin() + out_offset + (3 * out_capacity)))));
            break;
        default:
            std::cerr << "particle_sort_hip not implemented for dim " << dimension << std::endl;
            exit(-1);
    }
}

template <class T> 
void copy_according_to_index(thrust::device_vector<T>* in_dat_dv, thrust::device_vector<T>* out_dat_dv, 
    const thrust::device_vector<int>& new_idx_dv, int in_capacity, int out_capacity, int size, int dimension)
{
    copy_according_to_index<T>(in_dat_dv, out_dat_dv, new_idx_dv, in_capacity, out_capacity, 
        0, 0, size, dimension);
}

/*******************************************************************************/
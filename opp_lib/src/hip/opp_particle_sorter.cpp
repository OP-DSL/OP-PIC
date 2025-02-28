
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

#include "opp_hip.h"

//****************************************
thrust::device_vector<OPP_INT> ps_from_indices_dv;
thrust::device_vector<OPP_INT> hf_from_indices_dv;
thrust::device_vector<OPP_INT> hf_sequence_dv;

struct RandomFunctor
{
    mutable thrust::default_random_engine rng;
    RandomFunctor(unsigned int seed) : rng(seed) {}

    __device__ OPP_INT operator()(OPP_INT value) const {
        if (value == MAX_CELL_INDEX) {
            return MAX_CELL_INDEX;
        }
        else {
            thrust::uniform_int_distribution<OPP_INT> dist(0, 1000000000);
            return dist(rng);
        }
    }
};

//****************************************
void opp_particle_sort(opp_set set)
{ 
    particle_sort_device(set, false);
}

//****************************************
__global__ void copy_int(const int*__restrict in_dat_d, int*__restrict out_dat_d, const int*__restrict indices, 
                        int in_stride, int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;
    if (tid < size) {
        const int idx = indices[tid];
        for (int d = 0; d < dim; d++) {
            out_dat_d[out_offset + tid + d * out_stride] = in_dat_d[in_offset + idx + d * in_stride];
        }
    }
}
__global__ void copy_real(const double*__restrict in_dat_d, double*__restrict out_dat_d, const int*__restrict indices, 
                        int in_stride, int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;
    if (tid < size) {
        const int idx = indices[tid];
        for (int d = 0; d < dim; d++) {
            out_dat_d[out_offset + tid + d * out_stride] = in_dat_d[in_offset + idx + d * in_stride];
        }
    }
}
__global__ void copy_from_to_int(const int*__restrict in_dat_d, int*__restrict out_dat_d, const int*__restrict from_idx, 
            const int*__restrict to_idx, int in_stride, int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;
    if (tid < size) {
        const int f_idx = from_idx[tid];
        const int t_idx = to_idx[tid];
        for (int d = 0; d < dim; d++) {
            out_dat_d[out_offset + t_idx + d * out_stride] = in_dat_d[in_offset + f_idx + d * in_stride];
        }
    }
}
__global__ void copy_from_to_real(const double*__restrict in_dat_d, double*__restrict out_dat_d, const int*__restrict from_idx, 
            const int*__restrict to_idx, int in_stride, int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;
    if (tid < size) {
        const int f_idx = from_idx[tid];
        const int t_idx = to_idx[tid];
        for (int d = 0; d < dim; d++) {
            out_dat_d[out_offset + t_idx + d * out_stride] = in_dat_d[in_offset + f_idx + d * in_stride];
        }
    }
}

//****************************************
void sort_dat_according_to_index_int(opp_dat dat, const thrust::device_vector<int>& new_idx_dv, 
                                        int set_capacity, int size, bool shuffle, int out_start_idx)
{ 
    // if hole filling and second communicator iteration: arrange only the newly added particles
    // in to sorted_dat array and copy to the dat array 
    // else: arrange all and swap the array pointers of the dat

    // NOTE: Both commented thrust routine and device_kernel function has approx same performance
    // copy_according_to_index<int>(dat->thrust_int, dat->thrust_int_sort, new_idx_dv, 
    //         set_capacity, set_capacity, 0, out_start_idx, size, dat->dim);
    const int nblocks  = (size - 1) / opp_const_threads_per_block + 1;
    copy_int <<<nblocks, opp_const_threads_per_block>>> (
        opp_get_dev_raw_ptr(*(dat->thrust_int)),
        opp_get_dev_raw_ptr(*(dat->thrust_int_sort)),
        opp_get_dev_raw_ptr(new_idx_dv),
        set_capacity, set_capacity,
        0, out_start_idx,
        dat->dim, size);

    if (shuffle && OPP_comm_iteration != 0) {
        OPP_DEVICE_SYNCHRONIZE();
        for (int d = 0; d < dat->dim; d++) 
            thrust::copy(
                (dat->thrust_int_sort->begin() + d * set_capacity + out_start_idx), 
                (dat->thrust_int_sort->begin() + d * set_capacity + out_start_idx + size), 
                (dat->thrust_int->begin() + d * set_capacity + out_start_idx));
    }
    else {
        thrust::device_vector<OPP_INT>* tmp = dat->thrust_int;
        dat->thrust_int = dat->thrust_int_sort;
        dat->thrust_int_sort = tmp;

        dat->data_d = (char*)opp_get_dev_raw_ptr(*(dat->thrust_int));
        dat->data_swap_d = (char*)opp_get_dev_raw_ptr(*(dat->thrust_int_sort));
    }
}

//****************************************
void sort_dat_according_to_index_double(opp_dat dat, const thrust::device_vector<int>& new_idx_dv, 
                                        int set_capacity, int size, bool shuffle, int out_start_idx)
{ 
    // if hole filling and second communicator iteration: arrange only the newly added particles
    // in to sorted_dat array and copy to the dat array 
    // else: arrange all and swap the array pointers of the dat
    
    // NOTE: Both commented thrust routine and device_kernel function has approx same performance
    // copy_according_to_index<double>(dat->thrust_real, dat->thrust_real_sort, new_idx_dv, 
    //         set_capacity, set_capacity, 0, out_start_idx, size, dat->dim);
    const int nblocks  = (size - 1) / opp_const_threads_per_block + 1;
    copy_real <<<nblocks, opp_const_threads_per_block>>> (
        opp_get_dev_raw_ptr(*(dat->thrust_real)),
        opp_get_dev_raw_ptr(*(dat->thrust_real_sort)),
        opp_get_dev_raw_ptr(new_idx_dv),
        set_capacity, set_capacity,
        0, out_start_idx,
        dat->dim, size);

    if (shuffle && OPP_comm_iteration != 0) {
        OPP_DEVICE_SYNCHRONIZE();
        for (int d = 0; d < dat->dim; d++) 
            thrust::copy(
                (dat->thrust_real_sort->begin() + d * set_capacity + out_start_idx), 
                (dat->thrust_real_sort->begin() + d * set_capacity + out_start_idx + size), 
                (dat->thrust_real->begin() + d * set_capacity + out_start_idx));
    }
    else {
        thrust::device_vector<OPP_REAL>* tmp = dat->thrust_real;
        dat->thrust_real = dat->thrust_real_sort;
        dat->thrust_real_sort = tmp;

        dat->data_d = (char*)opp_get_dev_raw_ptr(*(dat->thrust_real));
        dat->data_swap_d = (char*)opp_get_dev_raw_ptr(*(dat->thrust_real_sort));
    }
}

//****************************************
// This assumes all the device data to be valid
void particle_sort_device(opp_set set, bool shuffle)
{ 
    const OPP_INT set_capacity = set->set_capacity;
    const OPP_INT set_size_plus_removed = set->size + set->particle_remove_count;
    
    OPP_INT sort_start_index = 0;
    OPP_INT sort_size = set_size_plus_removed;

    // in the second iteration of Move, sort only the newly added particles
    if (shuffle && OPP_comm_iteration != 0) {
        sort_start_index = set_size_plus_removed - set->diff;
        sort_size = set->diff;
    }

    if (OPP_DBG) 
        opp_printf("particle_sort_device", 
            "setSize[%d] setDiff[%d] setCap[%d] size+rem[%d] fill[%s] comIter[%d] startIdx[%d] sortSize[%d]", 
            set->size, set->diff, set_capacity, set_size_plus_removed, 
            (shuffle ? "TRUE" : "FALSE"), OPP_comm_iteration, sort_start_index, sort_size);

    thrust::device_vector<int>* cell_index_dv = set->mesh_relation_dat->thrust_int;
    thrust::device_vector<int>* ps_cell_index_dv = set->mesh_relation_dat->thrust_int_sort;

    if (shuffle) {
        // randomize the cell indices to minimize shared memory issues in later PIC routines
        // The below will create random numbers for each index, and MAX_CELL_INDEX for removed, 
        opp_profiler->start("PS_Shuffle");
        const unsigned int shuffle_seed = 123;
        thrust::transform(cell_index_dv->begin() + sort_start_index, cell_index_dv->begin() + set_size_plus_removed, 
                            ps_cell_index_dv->begin(), RandomFunctor(shuffle_seed));
        opp_profiler->end("PS_Shuffle");
    }
    else {
        // copy the cell index to the thrust vector for sorting
        opp_profiler->start("PS_CopyCID");
        thrust::copy(cell_index_dv->begin() + sort_start_index, cell_index_dv->begin() + set_size_plus_removed, 
                        ps_cell_index_dv->begin());
        opp_profiler->end("PS_CopyCID");
    }

    // Create a sequence of numbers starting from sort_start_index to be used as swap indices
    opp_profiler->start("PS_Sequence");
    ps_from_indices_dv.reserve(set->set_capacity);
    ps_from_indices_dv.resize(sort_size);
    thrust::sequence(ps_from_indices_dv.begin(), ps_from_indices_dv.end(), sort_start_index);
    opp_profiler->end("PS_Sequence");

    // sort ps_swap_indices_dp, with the key ps_cell_index_dp. Both keys and values will be sorted
    opp_profiler->start("PS_SortKey");
    thrust::sort_by_key(ps_cell_index_dv->begin(), ps_cell_index_dv->begin() + sort_size, 
                        ps_from_indices_dv.begin());
    opp_profiler->end("PS_SortKey");

    // Reorder the dats according to ps_swap_indices_dp (from mapping)
    opp_profiler->start("PS_Dats");
    for (opp_dat& dat : *(set->particle_dats)) {    

        if (!(strstr(dat->type, ":soa") != NULL || OPP_auto_soa || (dat->dim > 1))) {
            std::cerr << "particle_sort_device not implemented for non SOA data structures [dat " << 
                dat->name << "]" << std::endl;
        }

        if (strcmp(dat->type, "int") == 0) {

            sort_dat_according_to_index_int(dat, ps_from_indices_dv, 
                                        set_capacity, sort_size, shuffle, sort_start_index);
        }
        else if (strcmp(dat->type, "double") == 0) {

            sort_dat_according_to_index_double(dat, ps_from_indices_dv, 
                                        set_capacity, sort_size, shuffle, sort_start_index);
        }
        else {
            std::cerr << "particle_sort_device not implemented for type " << dat->type << " [dat " << 
                dat->name << "]" << std::endl;
            opp_abort();
        }
    }
    opp_profiler->end("PS_Dats");
}

//****************************************
// This assumes all the device data to be valid
void particle_hole_fill_device(opp_set set)
{ 
    const OPP_INT set_capacity = set->set_capacity;
    const OPP_INT part_remove_count = set->particle_remove_count;
    const OPP_INT set_size_plus_removed = set->size + part_remove_count;
    const OPP_INT nblocks  = (part_remove_count - 1) / opp_const_threads_per_block + 1;

    if (OPP_DBG) 
        opp_printf("particle_hole_fill_device", "remove=%d set_size+removed=%d capacity=%d", 
            part_remove_count, set_size_plus_removed, set_capacity);

    const OPP_INT sort_start_index = 0;
    const OPP_INT sort_size = set_size_plus_removed;

    // TODO : Optimization... Adjust the sort size when OPP_comm_iteration != 0
    // // in the second iteration of Move, sort only the newly added particles
    // if (OPP_comm_iteration != 0) 
    // {
    //     sort_start_index = set_size_plus_removed - set->diff;
    //     sort_size = set->diff;
    // }

    // sort OPP_remove_particle_indices_dv since it can be shuffled
    opp_profiler->start("HF_SORT");
    thrust::sort(thrust::device, OPP_remove_particle_indices_dv.begin(), 
                    OPP_remove_particle_indices_dv.begin() + part_remove_count);
    opp_profiler->end("HF_SORT");

    // resize hf_sequence_dv and hf_from_indices_dv if required
    if (hf_sequence_dv.capacity() < sort_size) {
        hf_sequence_dv.resize(set_capacity);
        thrust::sequence(hf_sequence_dv.begin(), hf_sequence_dv.end(), 0);

        hf_from_indices_dv.reserve(set_size_plus_removed);
    }
    hf_from_indices_dv.resize(set_size_plus_removed);

    // TODO : There may be a better implementation in sycl backend!!!
    // get the particle indices in reverse order whose cell index is not MAX_CELL_INDEX
    opp_profiler->start("HF_COPY_IF");
    auto end_iter1 = thrust::copy_if(thrust::device, 
        thrust::make_reverse_iterator(hf_sequence_dv.begin() + set_size_plus_removed), 
        thrust::make_reverse_iterator(hf_sequence_dv.begin() + sort_start_index), 
        thrust::make_reverse_iterator(set->mesh_relation_dat->thrust_int->begin() + set_size_plus_removed), 
        hf_from_indices_dv.begin(),
        [] __device__(int i) { return i != MAX_CELL_INDEX; });
    hf_from_indices_dv.resize(part_remove_count);
    opp_profiler->end("HF_COPY_IF");

    opp_profiler->start("HF_Dats");
    // For all the dats, fill the holes using the swap_indices
    for (opp_dat& dat : *(set->particle_dats)) {

        if (!(strstr(dat->type, ":soa") != NULL || OPP_auto_soa || (dat->dim > 1))) {
            std::cerr << "particle_hole_fill_device not implemented for non SOA data structures [dat " << 
                dat->name << "]" << std::endl;
        }

        const OPP_INT* from_indices = opp_get_dev_raw_ptr(hf_from_indices_dv);
        const OPP_INT* to_indices = opp_get_dev_raw_ptr(OPP_remove_particle_indices_dv);

        if (strcmp(dat->type, "int") == 0) {
            
            copy_from_to_int <<<nblocks, opp_const_threads_per_block>>> (
                opp_get_dev_raw_ptr(*(dat->thrust_int)),
                opp_get_dev_raw_ptr(*(dat->thrust_int)),
                from_indices, to_indices,
                set_capacity, set_capacity,
                0, 0,
                dat->dim, part_remove_count);
        }
        else if (strcmp(dat->type, "double") == 0) {
            
            copy_from_to_real <<<nblocks, opp_const_threads_per_block>>> (
                opp_get_dev_raw_ptr(*(dat->thrust_real)),
                opp_get_dev_raw_ptr(*(dat->thrust_real)),
                from_indices, to_indices,
                set_capacity, set_capacity,
                0, 0,
                dat->dim, part_remove_count);
        }
        else {
            std::cerr << "particle_hole_fill_device not implemented for type " << dat->type << " [dat " << 
                dat->name << "]" << std::endl;
        }
    }
    OPP_DEVICE_SYNCHRONIZE();
    opp_profiler->end("HF_Dats");
}

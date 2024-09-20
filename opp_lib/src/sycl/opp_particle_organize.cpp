
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

#include <opp_sycl.h>

size_t hf_from_indices_size = 0;
OPP_INT* hf_from_indices_dp = nullptr; // holefill from indices - starting in reverse order
size_t ps_swap_indices_size = 0;
OPP_INT* ps_swap_indices_dp = nullptr; // swap (from mapping) indices required for sorting/swapping

//****************************************
void opp_particle_sort(opp_set set)
{ 
    particle_sort_device(set, false);
}

//****************************************
template <typename T>
void sort_dat_according_to_index(opp_dat dat, const OPP_INT* swap_indices, 
    const int set_capacity, const int size, bool shuffle, const int out_start_idx)
{ 
    // if hole filling and second communicator iteration: arrange only the newly added particles
    // in to sorted_dat array and copy to the dat array 
    // else: arrange all and swap the array pointers of the dat
    
    // NOTE: Both commented thrust routine and device_kernel function has approx same performance
    // copy_according_to_index<T>((const T*)dat->data_d, (T*)dat->data_swap_d, 
    //         swap_indices, set_capacity, 
    //         set_capacity, 0, out_start_idx, size, dat->dim);

    T *dat_data = (T *)(dat->data_d);
    T *dat_swap_data = (T *)(dat->data_swap_d);
    const int *swap_indices_ptr = swap_indices; 
    const int dat_dim = dat->dim;

    const int nblocks  = (size - 1) / opp_const_threads_per_block + 1;

    // data will be swapped from dat->data_d to dat->data_swap_d according to swap_indices
    opp_queue->submit([&](sycl::handler &cgh) {
        cgh.parallel_for(
            sycl::nd_range<1>(opp_const_threads_per_block*nblocks,opp_const_threads_per_block),
            [=](sycl::nd_item<1> item) {
                copy_from<T>(
                    dat_data, dat_swap_data,
                    swap_indices_ptr,
                    set_capacity, set_capacity, 
                    0, out_start_idx,
                    dat_dim, size, 
                    item);
            });
    });

    if (shuffle && OPP_comm_iteration != 0) {
        OPP_DEVICE_SYNCHRONIZE();

        // Now the valid data is in dat->data_swap_d, but since it is not the 1st communicator iteration
        // we only swap the newly added indices (since newly added particle count is small compared to 
        // existing particles), simply copy the required indices from dat->data_swap_d to dat->data_d
        for (int d = 0; d < dat->dim; d++)
            std::copy(dpl::execution::make_device_policy(*opp_queue),
                (dat_swap_data + d * set_capacity + out_start_idx),
                (dat_swap_data + d * set_capacity + out_start_idx + size),
                (dat_data + d * set_capacity + out_start_idx));
    }
    else {
        // Now the valid data is in dat->data_swap_d, simply switch the data_swap_d and data_d pointers
        char *tmp = dat->data_swap_d;
        dat->data_swap_d = dat->data_d;
        dat->data_d = tmp;
    }
}

//****************************************
// This assumes all the device data to be valid
void particle_sort_device(opp_set set, bool shuffle)
{ 
    const int set_capacity = set->set_capacity;
    const int set_size_plus_removed = set->size + set->particle_remove_count;
    
    int sort_start_index = 0;
    int sort_size = set_size_plus_removed;

    // in the second iteration of Move, sort only the newly added particles
    if (shuffle && OPP_comm_iteration != 0) {
        sort_start_index = set_size_plus_removed - set->diff;
        sort_size = set->diff;
    }

    if (OPP_DBG) 
        opp_printf("particle_sort_device", 
            "size[%d] diff[%d] cap[%d] size+rem[%d] shuffle[%s] comIter[%d] startIdx[%d] sortSize[%d]", 
            set->size, set->diff, set_capacity, set_size_plus_removed, 
            (shuffle ? "TRUE" : "FALSE"), OPP_comm_iteration, sort_start_index, sort_size);
    
    OPP_INT* cellIdx_dp = (OPP_INT*)set->mesh_relation_dat->data_d; // this is cell id data of particle set
    OPP_INT* ps_cell_index_dp = (OPP_INT*)set->mesh_relation_dat->data_swap_d; // temp using swap array

    if (shuffle) {
        // randomize the cell indices to minimize shared memory issues in later PIC routines
        // The below will create random numbers for each index, and MAX_CELL_INDEX for removed, 
        opp_profiler->start("PS_Shuffle");
        const uint32_t shuffle_seed = 123;
        const int size = (set_size_plus_removed - sort_start_index);
        opp_queue->submit([&](sycl::handler& cgh) {
            cgh.parallel_for(sycl::nd_range<1>(OPP_gpu_threads_per_block * const_blocks, OPP_gpu_threads_per_block),
                [=](sycl::nd_item<1> item) {
                    dpl::uniform_int_distribution<int> dist(0, 1000);
                    for (size_t i = item.get_global_linear_id(); i < size; i += item.get_global_range()[0] ){
                        dpl::minstd_rand engine(shuffle_seed, i);
                        ps_cell_index_dp[i] = (cellIdx_dp[i + sort_start_index] == MAX_CELL_INDEX) ? 
                                                MAX_CELL_INDEX : dist(engine);
                    }
                }
            );
        });
        OPP_DEVICE_SYNCHRONIZE();
        opp_profiler->end("PS_Shuffle");
    }
    else {
        // copy the cell index to the thrust vector for sorting
        opp_profiler->start("PS_CopyCID");
        const int size = (set_size_plus_removed - sort_start_index);
        opp_queue->submit([&](sycl::handler& cgh) {
            cgh.parallel_for(sycl::nd_range<1>(OPP_gpu_threads_per_block * const_blocks, OPP_gpu_threads_per_block),
                [=](sycl::nd_item<1> item) {
                    for (int i = item.get_global_linear_id(); i < size; i += item.get_global_range()[0] ){
                        ps_cell_index_dp[i] = cellIdx_dp[i + sort_start_index];
                    }
                }
            );
        });
        OPP_DEVICE_SYNCHRONIZE();
        opp_profiler->end("PS_CopyCID");
    }
    
    // Create a sequence of numbers starting from sort_start_index to be used as swap indices
    opp_profiler->start("PS_Sequence");
    opp_mem::dev_resize<OPP_INT>(ps_swap_indices_dp, ps_swap_indices_size, (size_t)set->set_capacity);
    opp_queue->submit([&](sycl::handler& cgh) {
        OPP_INT* ps_swap_indices = ps_swap_indices_dp;
        cgh.parallel_for(sycl::nd_range<1>(OPP_gpu_threads_per_block * const_blocks, OPP_gpu_threads_per_block),
            [=](sycl::nd_item<1> item) {
                for (int i = item.get_global_linear_id(); i < sort_size; i += item.get_global_range()[0] ){
                    ps_swap_indices[i] = i;
                }
            }
        );
    });
    OPP_DEVICE_SYNCHRONIZE();
    opp_profiler->end("PS_Sequence");

    // sort ps_swap_indices_dp, with the key ps_cell_index_dp. Both keys and values will be sorted
    opp_profiler->start("PS_SortKey");
    dpl::sort_by_key(dpl::execution::make_device_policy(*opp_queue), 
        ps_cell_index_dp, ps_cell_index_dp + sort_size,
        ps_swap_indices_dp);
    opp_profiler->end("PS_SortKey");

    // Reorder the dats according to ps_swap_indices_dp (from mapping)
    opp_profiler->start("PS_Dats");
    for (opp_dat& dat : *(set->particle_dats)) {

        if (!(strstr(dat->type, ":soa") != NULL || OPP_auto_soa || (dat->dim > 1))) {
            std::cerr << "particle_sort_device not implemented for non SOA data structures [dat "
                << dat->name << "]" << std::endl;
            opp_abort();
        }

        if (strcmp(dat->type, "int") == 0) {

            sort_dat_according_to_index<OPP_INT>(dat, ps_swap_indices_dp, 
                                    set_capacity, sort_size, shuffle, sort_start_index);
        }
        else if (strcmp(dat->type, "double") == 0) {

            sort_dat_according_to_index<OPP_REAL>(dat, ps_swap_indices_dp, 
                                    set_capacity, sort_size, shuffle, sort_start_index);
        }
        else {
            std::cerr << "particle_sort_device not implemented for type " << dat->type << " [dat " 
                << dat->name << "]" << std::endl;
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

    // sort OPP_remove_particle_indices_d since it can be in shuffled state
    opp_profiler->start("HF_SORT");
    dpl::sort(dpl::execution::make_device_policy(*opp_queue),
        OPP_remove_particle_indices_d, OPP_remove_particle_indices_d + part_remove_count);
    opp_profiler->end("HF_SORT");

    // resize only if hf_from_indices_size < set_capacity
    opp_mem::dev_resize<OPP_INT>(hf_from_indices_dp, hf_from_indices_size, (size_t)set_capacity);

    OPP_INT* swap_indices_dp = (OPP_INT *)set->particle_remove_count_d; // using dev ptr temporarily
    OPP_INT tmp = set_size_plus_removed - 1;
    opp_mem::copy_host_to_dev<int>(swap_indices_dp, &tmp, 1);

    // Find the hole fill from indices (order may not be preserved)
    opp_profiler->start("HF_COPY_IF");
    opp_queue->submit([&](sycl::handler& cgh) {
        const OPP_INT* cid_dp = (OPP_INT*)set->mesh_relation_dat->data_d;
        OPP_INT* hf_from_indices = hf_from_indices_dp;
        cgh.parallel_for(
            sycl::nd_range<1>(opp_const_threads_per_block*nblocks,opp_const_threads_per_block),
            [=](sycl::nd_item<1> item) {   
                const OPP_INT idx = item.get_global_linear_id();
                if (idx < part_remove_count) {
                    while (true) {
                        const OPP_INT pos = opp_atomic_fetch_add(swap_indices_dp, -1);
                        if (cid_dp[pos] != MAX_CELL_INDEX) {
                            hf_from_indices[idx] = pos;
                            break;
                        }
                    }
                }         
            }
        );
    });
    OPP_DEVICE_SYNCHRONIZE();
    opp_profiler->end("HF_COPY_IF");

    // Sort the hole fill from indices to avoid hole-filling issues
    opp_profiler->start("HF_COPY_IF_SORT");
    dpl::sort(dpl::execution::make_device_policy(*opp_queue), 
                hf_from_indices_dp, hf_from_indices_dp + part_remove_count, 
                std::greater<int>());
    opp_profiler->end("HF_COPY_IF_SORT");

    // For all the dats, fill the holes using the swap_indices
    opp_profiler->start("HF_Dats");
    for (opp_dat& dat : *(set->particle_dats)) {

        if (!(strstr(dat->type, ":soa") != NULL || OPP_auto_soa || (dat->dim > 1))) {
            std::cerr << "particle_hole_fill_device not implemented for non SOA data structures [dat " 
                << dat->name << "]" << std::endl;
            opp_abort();
        }

        if (strcmp(dat->type, "int") == 0) {

            opp_queue->submit([&](sycl::handler &cgh) {
                OPP_INT *dat_data = (OPP_INT *)(dat->data_d);
                const int *from_indices = hf_from_indices_dp;
                const int *remove_indices = OPP_remove_particle_indices_d;
                const int dat_dim = dat->dim;
                cgh.parallel_for(
                    sycl::nd_range<1>(opp_const_threads_per_block*nblocks,opp_const_threads_per_block),
                    [=](sycl::nd_item<1> item) {
                        copy_from_to<OPP_INT>(
                            dat_data, dat_data,
                            from_indices, remove_indices,
                            set_capacity, set_capacity, 
                            0, 0,
                            dat_dim, part_remove_count, 
                            item);
                    });
            });
        }
        else if (strcmp(dat->type, "double") == 0) {

            opp_queue->submit([&](sycl::handler &cgh) {
                OPP_REAL *dat_data = (OPP_REAL *)(dat->data_d);
                const int *from_indices = hf_from_indices_dp;
                const int *remove_indices = OPP_remove_particle_indices_d;
                const int dat_dim = dat->dim;
                cgh.parallel_for(
                    sycl::nd_range<1>(opp_const_threads_per_block*nblocks,opp_const_threads_per_block),
                    [=](sycl::nd_item<1> item) { 
                        copy_from_to<OPP_REAL>(
                            dat_data, dat_data,
                            from_indices, remove_indices,
                            set_capacity, set_capacity, 
                            0, 0,
                            dat_dim, part_remove_count, 
                            item);
                    });
            });
        }
        else {
            std::cerr << "particle_hole_fill_device not implemented for type " << dat->type << " [dat " 
                << dat->name << "]" << std::endl;
            opp_abort();
        }
    }
    OPP_DEVICE_SYNCHRONIZE();
    opp_profiler->end("HF_Dats");
}
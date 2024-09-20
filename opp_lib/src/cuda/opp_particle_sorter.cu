
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

#include <opp_cuda.h>

//****************************************
unsigned int seed = 123; // Seed for random number generator
struct RandomFunctor
{
    mutable thrust::default_random_engine rng;

    RandomFunctor(unsigned int seed) : rng(seed) {}

    __device__ int operator()(int value) const
    {
        if (value == MAX_CELL_INDEX)
        {
            return value; // Keep the value of MAX_CELL_INDEX unchanged
        }
        else
        {
            thrust::uniform_int_distribution<int> dist(0, 1000000000);
            return dist(rng); // Assign random number < 100

            // thrust::normal_distribution<float> dist(50.0f, 10.0f);
            // return static_cast<int>(dist(rng)); // Convert to integer
        }
    }
};

__global__ void copy_int(const int* in_dat_d, int* out_dat_d, const int* indices, int in_stride, 
                                    int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < size) 
    {
        const int idx = indices[tid];
        for (int d = 0; d < dim; d++)
        {
            out_dat_d[out_offset + tid + d * out_stride] = in_dat_d[in_offset + idx + d * in_stride];
        }
    }
}

__global__ void copy_double(const double* in_dat_d, double* out_dat_d, const int* indices, int in_stride, 
                                    int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < size) 
    {
        const int idx = indices[tid];
        for (int d = 0; d < dim; d++)
        {
            out_dat_d[out_offset + tid + d * out_stride] = in_dat_d[in_offset + idx + d * in_stride];
        }
    }
}

//****************************************

thrust::device_vector<int> cellIdx_dv;
thrust::device_vector<int> i_dv;

//****************************************
// This assumes all the device data to be valid
void particle_sort_device(opp_set set, bool hole_filling)
{ 
    int set_capacity = set->set_capacity;
    int set_size_plus_removed = set->size + set->particle_remove_count;
    
    int sort_start_index = 0;
    int sort_size = set_size_plus_removed;

    // in the second iteration of Move, sort only the newly added particles
    if (hole_filling && OPP_comm_iteration != 0) 
    {
        sort_start_index = set_size_plus_removed - set->diff;
        sort_size = set->diff;
    }

    if (OPP_DBG) 
        opp_printf("particle_sort_device", 
            "setSize[%d] setDiff[%d] setCap[%d] size+rem[%d] fill[%s] comIter[%d] startIdx[%d] sortSize[%d]", 
            set->size, set->diff, set_capacity, set_size_plus_removed, 
            (hole_filling ? "TRUE" : "FALSE"), OPP_comm_iteration, sort_start_index, sort_size);
    
    opp_profiler->start("ZS_Resize");
    cellIdx_dv.reserve(set->set_capacity);
    cellIdx_dv.resize(sort_size);
    opp_profiler->end("ZS_Resize");

    opp_profiler->start("ZS_Copy");
    // copy the cell index to the thrust vector
    int* cellIdx_dp = (int*)set->mesh_relation_dat->data_d;
    thrust::copy(cellIdx_dp + sort_start_index, cellIdx_dp + set_size_plus_removed, cellIdx_dv.begin());
    opp_profiler->end("ZS_Copy");

    opp_profiler->start("ZS_HoleFill");
    if (hole_filling)
    {
        // in hole filling, randomize the cell indices to minimize shared memory issues
        // The below will create random numbers for each index, and MAX_CELL_INDEX for removed, 
        // ideally this should not be called cell_Idx_dv, better naming would be something like, random ordering
        thrust::transform(cellIdx_dv.begin(), cellIdx_dv.end(), 
            cellIdx_dv.begin(), RandomFunctor(seed));

        // thrust::replace_if(
        //     cellIdx_dv.begin(), cellIdx_dv.end(),
        //     thrust::placeholders::_1 != MAX_CELL_INDEX,  // Replace all values != MAX_CELL_INDEX
        //     0                                            // New value to assign (zero in this case)
        // );
    }
    opp_profiler->end("ZS_HoleFill");

    opp_profiler->start("ZS_Sequence");
    i_dv.reserve(set->set_capacity);
    i_dv.resize(sort_size);
    thrust::sequence(i_dv.begin(), i_dv.end(), sort_start_index);
    opp_profiler->end("ZS_Sequence");

    // int dis = (int)thrust::distance(i_dv.begin(), i_dv.end());
    // int dis2 = (int)thrust::distance(cellIdx_dv.begin(), cellIdx_dv.end());
    // opp_printf("SORT", "set->size=%d set_size_plus_removed=%d | size %d capacity %d i_dv=%d cellIdx_dv=%d", set->size, set_size_plus_removed, sort_size, set->set_capacity, dis, dis2);

    // opp_profiler->start("XSortKey");
    if (OPP_comm_iteration == 0) opp_profiler->start("ZS_SortKey0");
    else if (OPP_comm_iteration == 1) opp_profiler->start("ZS_SortKey1");
    else opp_profiler->start("ZS_SortKey");
    thrust::sort_by_key(cellIdx_dv.begin(), cellIdx_dv.end(), i_dv.begin());
    // opp_profiler->end("XSortKey");
    if (OPP_comm_iteration == 0) opp_profiler->end("ZS_SortKey0");
    else if (OPP_comm_iteration == 1) opp_profiler->end("ZS_SortKey1");
    else opp_profiler->end("ZS_SortKey");

    opp_profiler->start("ZS_Dats");
    for (int i = 0; i < (int)set->particle_dats->size(); i++)
    {    
        opp_dat& dat = set->particle_dats->at(i);

        if (!(strstr(dat->type, ":soa") != NULL || OPP_auto_soa || (dat->dim > 1)))
        {
            std::cerr << "particle_sort_device not implemented for non SOA data structures [dat " << 
                dat->name << "]" << std::endl;
        }

        if (strcmp(dat->type, "int") == 0)
        {
            sort_dat_according_to_index_int(dat, i_dv, set_capacity, sort_size, 
                hole_filling, sort_start_index);
        }
        else if (strcmp(dat->type, "double") == 0)
        {
            sort_dat_according_to_index_double(dat, i_dv, set_capacity, sort_size, 
                hole_filling, sort_start_index);
        }
        else
        {
            std::cerr << "particle_sort_device not implemented for type " << dat->type << " [dat " << 
                dat->name << "]" << std::endl;
        }
    }
    opp_profiler->end("ZS_Dats");
}

//****************************************
void sort_dat_according_to_index_int(opp_dat dat, const thrust::device_vector<int>& new_idx_dv, 
                                        int set_capacity, int size, bool hole_filling, int out_start_idx)
{ 

    // if hole filling and second communicator iteration: arrange only the newly added particles
    // in to sorted_dat array and copy to the dat array 
    // else: arrange all and swap the array pointers of the dat

    // NOTE: Both commented thrust routine and device_kernel function has approx same performance
    // copy_according_to_index<int>(dat->thrust_int, dat->thrust_int_sort, new_idx_dv, 
    //         set_capacity, set_capacity, 0, out_start_idx, size, dat->dim);
    const int nblocks  = (size - 1) / 192 + 1;
    copy_int <<<nblocks, 192>>> (
        (int*)thrust::raw_pointer_cast(dat->thrust_int->data()),
        (int*)thrust::raw_pointer_cast(dat->thrust_int_sort->data()),
        (int*)thrust::raw_pointer_cast(new_idx_dv.data()),
        set_capacity,
        set_capacity,
        0,
        out_start_idx,
        dat->dim,
        size);

    if (hole_filling && OPP_comm_iteration != 0) 
    {
        cutilSafeCall(cudaDeviceSynchronize());
        for (int d = 0; d < dat->dim; d++) 
            thrust::copy(
                (dat->thrust_int_sort->begin() + d * set_capacity + out_start_idx), 
                (dat->thrust_int_sort->begin() + d * set_capacity + out_start_idx + size), 
                (dat->thrust_int->begin() + d * set_capacity + out_start_idx));
    }
    else
    {
        thrust::device_vector<int>* tmp = dat->thrust_int;
        dat->thrust_int = dat->thrust_int_sort;
        dat->thrust_int_sort = tmp;

        dat->data_d = (char*)thrust::raw_pointer_cast(dat->thrust_int->data());
    }
}

//****************************************
void sort_dat_according_to_index_double(opp_dat dat, const thrust::device_vector<int>& new_idx_dv, 
                                        int set_capacity, int size, bool hole_filling, int out_start_idx)
{ 

    // if hole filling and second communicator iteration: arrange only the newly added particles
    // in to sorted_dat array and copy to the dat array 
    // else: arrange all and swap the array pointers of the dat
    
    // NOTE: Both commented thrust routine and device_kernel function has approx same performance
    // copy_according_to_index<double>(dat->thrust_real, dat->thrust_real_sort, new_idx_dv, 
    //         set_capacity, set_capacity, 0, out_start_idx, size, dat->dim);
    const int nblocks  = (size - 1) / 192 + 1;
    copy_double <<<nblocks, 192>>> (
        (double*)thrust::raw_pointer_cast(dat->thrust_real->data()),
        (double*)thrust::raw_pointer_cast(dat->thrust_real_sort->data()),
        (int*)thrust::raw_pointer_cast(new_idx_dv.data()),
        set_capacity,
        set_capacity,
        0,
        out_start_idx,
        dat->dim,
        size);

    if (hole_filling && OPP_comm_iteration != 0) 
    {
        cutilSafeCall(cudaDeviceSynchronize());
        for (int d = 0; d < dat->dim; d++) 
            thrust::copy(
                (dat->thrust_real_sort->begin() + d * set_capacity + out_start_idx), 
                (dat->thrust_real_sort->begin() + d * set_capacity + out_start_idx + size), 
                (dat->thrust_real->begin() + d * set_capacity + out_start_idx));
    }
    else
    {
        thrust::device_vector<double>* tmp = dat->thrust_real;
        dat->thrust_real = dat->thrust_real_sort;
        dat->thrust_real_sort = tmp;

        dat->data_d = (char*)thrust::raw_pointer_cast(dat->thrust_real->data());
    }
}

//****************************************
//****************************************
__global__ void copy_intFromTo(const int* in_dat_d, int* out_dat_d, const int* from_idx, const int* to_idx, 
                        int in_stride, int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < size) 
    {
        const int f_idx = from_idx[tid];
        const int t_idx = to_idx[tid];
        for (int d = 0; d < dim; d++)
        {
            out_dat_d[out_offset + t_idx + d * out_stride] = in_dat_d[in_offset + f_idx + d * in_stride];
        }
    }
}

__global__ void copy_doubleFromTo(const double* in_dat_d, double* out_dat_d, const int* from_idx, 
            const int* to_idx, int in_stride, int out_stride, int in_offset, int out_offset, int dim, int size) 
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < size) 
    {
        const int f_idx = from_idx[tid];
        const int t_idx = to_idx[tid];
        for (int d = 0; d < dim; d++)
        {
            out_dat_d[out_offset + t_idx + d * out_stride] = in_dat_d[in_offset + f_idx + d * in_stride];
        }
    }
}

__global__ void matchToFrom(const int* cids, int* from_idx, int* from, const int size, const int max_remove_idx) 
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < size) 
    {
        int tmp = 0;
        while (cids[tmp] == MAX_CELL_INDEX && tmp > max_remove_idx) {
            tmp = atomicAdd(from, -1);
        }
        from_idx[tid] = tmp;
    }
}

//****************************************


thrust::device_vector<int> ps_to_indices_dv;
thrust::device_vector<int> ps_from_indices_dv;
thrust::device_vector<int> ps_sequence_dv;

//****************************************
// This assumes all the device data to be valid
void particle_hole_fill_device(opp_set set)
{ 
    const int set_capacity = set->set_capacity;
    const int part_remove_count = set->particle_remove_count;
    const int set_size_plus_removed = set->size + part_remove_count;
    const int nblocks  = (part_remove_count - 1) / 192 + 1;

    if (OPP_DBG) opp_printf("particle_hole_fill_device", "remove_count=%d set_size+removed=%d set_capacity=%d", 
        part_remove_count, set_size_plus_removed, set_capacity);

    int sort_start_index = 0;
    int sort_size = set_size_plus_removed;

    // TODO : Optimization... Adjust the sort size when OPP_comm_iteration != 0
    // // in the second iteration of Move, sort only the newly added particles
    // if (OPP_comm_iteration != 0) 
    // {
    //     sort_start_index = set_size_plus_removed - set->diff;
    //     sort_size = set->diff;
    // }

    opp_profiler->start("ZF_SORT");
    // sort OPP_thrust_remove_particle_indices_d since it can be shuffled
    thrust::sort(thrust::device, OPP_thrust_remove_particle_indices_d.begin(), 
                    OPP_thrust_remove_particle_indices_d.begin() + part_remove_count);
    opp_profiler->end("ZF_SORT");

    // resize ps_sequence_dv and ps_from_indices_dv if required
    if (ps_sequence_dv.capacity() < sort_size) {
        ps_sequence_dv.resize(set_capacity);
        thrust::sequence(ps_sequence_dv.begin(), ps_sequence_dv.end(), 0);

        ps_from_indices_dv.reserve(set_size_plus_removed);
    }
    ps_from_indices_dv.resize(set_size_plus_removed);

    // get the particle indices in reverse order whose cell index is not MAX_CELL_INDEX
    opp_profiler->start("ZF_COPY_IF");
    auto end_iter1 = thrust::copy_if(thrust::device, 
        thrust::make_reverse_iterator(ps_sequence_dv.begin() + set_size_plus_removed), 
        thrust::make_reverse_iterator(ps_sequence_dv.begin() + sort_start_index), 
        thrust::make_reverse_iterator(set->mesh_relation_dat->thrust_int->begin() + set_size_plus_removed), 
        ps_from_indices_dv.begin(),
        [] __device__(int i) { return i != MAX_CELL_INDEX; });
    ps_from_indices_dv.resize(part_remove_count);
    opp_profiler->end("ZF_COPY_IF");

    opp_profiler->start("ZF_Dats");
    // For all the dats, fill the holes using the swap_indices
    for (opp_dat& dat : *(set->particle_dats))
    {
        if (!(strstr(dat->type, ":soa") != NULL || OPP_auto_soa || (dat->dim > 1)))
        {
            std::cerr << "particle_hole_fill_device not implemented for non SOA data structures [dat " << 
                dat->name << "]" << std::endl;
        }

        if (strcmp(dat->type, "int") == 0)
        {
            copy_intFromTo <<<nblocks, 192>>> (
                (int*)thrust::raw_pointer_cast(dat->thrust_int->data()),
                (int*)thrust::raw_pointer_cast(dat->thrust_int->data()),
                (int*)thrust::raw_pointer_cast(ps_from_indices_dv.data()),
                (int*)thrust::raw_pointer_cast(OPP_thrust_remove_particle_indices_d.data()),
                set_capacity,
                set_capacity,
                0,
                sort_start_index,
                dat->dim,
                part_remove_count);
        }
        else if (strcmp(dat->type, "double") == 0)
        {
            copy_doubleFromTo <<<nblocks, 192>>> (
                (double*)thrust::raw_pointer_cast(dat->thrust_real->data()),
                (double*)thrust::raw_pointer_cast(dat->thrust_real->data()),
                (int*)thrust::raw_pointer_cast(ps_from_indices_dv.data()),
                (int*)thrust::raw_pointer_cast(OPP_thrust_remove_particle_indices_d.data()),
                set_capacity,
                set_capacity,
                0,
                sort_start_index,
                dat->dim,
                part_remove_count);
        }
        else
        {
            std::cerr << "particle_hole_fill_device not implemented for type " << dat->type << " [dat " << 
                dat->name << "]" << std::endl;
        }
    }
    cutilSafeCall(cudaDeviceSynchronize());
    opp_profiler->end("ZF_Dats");
}
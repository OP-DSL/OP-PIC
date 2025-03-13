
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

#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include <dpct/dpl_utils.hpp>
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/random>

#include <opp_lib.h>

#ifdef USE_MPI
    #include <opp_mpi_core.h>
#endif

#define trunc sycl::trunc

struct opp_dh_indices {
    OPP_INT* move_count = nullptr;
    OPP_INT* part_indices = nullptr;
    OPP_INT* rank_indices = nullptr;
    OPP_INT* cell_indices = nullptr;
    OPP_INT capacity = 0;
};

constexpr bool debug_mem = false;
constexpr bool debugger = false;
constexpr int opp_const_threads_per_block = 192;
constexpr int const_blocks = 200;

#define OPP_GPU_THREADS_PER_BLOCK 32

#define OPP_PARTICLE_MOVE_DONE { opp_move_status_flag = OPP_MOVE_DONE; }
#define OPP_PARTICLE_NEED_MOVE { opp_move_status_flag = OPP_NEED_MOVE; }
#define OPP_PARTICLE_NEED_REMOVE { opp_move_status_flag = OPP_NEED_REMOVE; }
#define OPP_DO_ONCE (opp_move_hop_iter_one_flag)
#define OPP_MOVE_RESET_FLAGS { opp_move_status_flag = OPP_MOVE_DONE; opp_move_hop_iter_one_flag = true; }

#define OPP_DEVICE_SYNCHRONIZE()  opp_queue->wait()

#define OPP_DEVICE_GLOBAL_LINEAR_ID (item.get_global_linear_id())
#define OPP_GLOBAL_FUNCTION inline
#define OPP_DEVICE_FUNCTION 
#define ADDITIONAL_PARAMETERS , sycl::nd_item<1> item
#define OPP_ATOMIC_FETCH_ADD(address, value) opp_atomic_fetch_add(address, value)

/*******************************************************************************/
extern int* opp_saved_mesh_relation_d;

extern OPP_INT* hf_from_indices_dp;
extern OPP_INT* ps_from_indices_dp;

extern int *OPP_move_particle_indices_d;
extern int *OPP_move_cell_indices_d;
extern int *OPP_move_count_d;
extern int *OPP_remove_particle_indices_d;

extern std::map<int, std::vector<OPP_INT>> cell_indices_hv; // cellid in the foreign rank, arrange according to rank
extern std::map<int, std::vector<OPP_INT>> particle_indices_hv; // particle ids to send, arrange according to rank
extern std::map<int, dpct::device_vector<OPP_INT>> particle_indices_dv;
extern std::map<int, dpct::device_vector<char>> send_data;
extern std::map<int, dpct::device_vector<char>> recv_data;

// arrays for global constants and reductions
extern int OPP_consts_bytes;
extern int OPP_reduct_bytes;
extern char *OPP_reduct_h, *OPP_reduct_d;
extern char *OPP_consts_h, *OPP_consts_d;

extern std::vector<char*> opp_consts;

extern char opp_move_status_flag;
extern bool opp_move_hop_iter_one_flag;
extern OPP_INT* opp_p2c;
extern OPP_INT* opp_c2c;

extern opp_dh_indices dh_indices_d;
extern opp_dh_indices dh_indices_h;

//*************************************************************************************************
void opp_device_init(int argc, char **argv);
void opp_device_exit();

void opp_upload_map(opp_map map, bool create_new = false);
void opp_create_dat_device_arrays(opp_dat dat, bool create_new = false);

/*******************************************************************************/
void opp_halo_create();
void opp_halo_destroy();

/*******************************************************************************/
void opp_init_double_indirect_reductions_device(int nargs, opp_arg *args);
void opp_exchange_double_indirect_reductions_device(int nargs, opp_arg *args) ;
void opp_complete_double_indirect_reductions_device(int nargs, opp_arg *args);

/*******************************************************************************/
void print_dat_to_txtfile_mpi(opp_dat dat, const char *file_name);
void opp_mpi_print_dat_to_txtfile(opp_dat dat, const char *file_name);

/*******************************************************************************/
void particle_sort_device(opp_set set, bool hole_filling);

/*******************************************************************************/
// template<typename T>
// inline void opp_mpi_reduce(opp_arg *args, T *data) 
// {
// #ifdef USE_MPI
//     if constexpr (std::is_same<T, double>::value) {
//         opp_mpi_reduce_double(args, data);
//     } else if constexpr (std::is_same<T, int>::value) {
//         opp_mpi_reduce_int(args, data);
//     } else {
//         static_assert(std::is_same<T, double>::value || std::is_same<T, int>::value, 
//                       "Unsupported data type for opp_mpi_reduce.");
//     }
// #else
//     (void)args;
//     (void)data;
// #endif
// }

template<typename T, typename std::enable_if<std::is_same<T, OPP_REAL>::value, int>::type = 0>
void opp_mpi_reduce(opp_arg *args, T *data) 
{
#ifdef USE_MPI
    opp_mpi_reduce_double(args, data);
#else
    (void)args;
    (void)data;
#endif
}

template<typename T, typename std::enable_if<std::is_same<T, OPP_INT>::value, int>::type = 0>
void opp_mpi_reduce(opp_arg *args, T *data) 
{
#ifdef USE_MPI
    opp_mpi_reduce_int(args, data);
#else
    (void)args;
    (void)data;
#endif
}

/*******************************************************************************/
class opp_mem {

public:

    // Allocate host memory
    template <typename T>
    inline static T* host_malloc(size_t count) {
        return (T*)malloc(count * sizeof(T));
    }

    // Free host memory
    template <typename T>
    inline static void host_free(T* ptr) {
        if (ptr)
            free(ptr);
        ptr = nullptr;
    }

    // Reallocate host memory
    template <typename T>
    inline static void host_realloc(T*& ptr, size_t new_size) {
        T* tmp_ptr = (T*)realloc(ptr, new_size);
        ptr = tmp_ptr;
    }

    // Allocate device memory
    template <typename T>
    inline static T* dev_malloc(size_t count) {
        if (count <= 0) return nullptr;
        try {
            T* ptr = sycl::malloc_device<T>(count * sizeof(T), *opp_queue);
            if (debug_mem) opp_printf("dev_malloc", "[%p][%zu]", ptr, count);
            return ptr;
        }
        catch (const sycl::exception &e) {
            throw std::runtime_error(std::string("dev_malloc: ") + e.what());
        }
    }

    // Free device memory
    template <typename T>
    inline static void dev_free(T*& ptr) {
        if (ptr) {
            if (debug_mem) opp_printf("dev_free", "[%p]", ptr);
            sycl::free(ptr, *opp_queue);
            ptr = nullptr;
        }
    }

    // Copy memory from one device pointer to another
    template <typename T>
    inline static void dev_memcpy(T*& dst, const T* src, size_t cpy_count) {
        if (debug_mem) opp_printf("dev_memcpy", "[%p]->[%p] cpy_count[%zu]", src, dst, cpy_count);
        opp_queue->memcpy(dst, src, cpy_count * sizeof(T)).wait();
    }

    // Resize device memory
    template <typename T>
    inline static void dev_realloc(T*& ptr, size_t& current_size, const size_t& new_size) {
        if (new_size <= 0) 
            throw std::runtime_error("dev_realloc: New Realloc size invalid - " + std::to_string(new_size));
        T* new_ptr = opp_mem::dev_malloc<T>(new_size);
        if (debug_mem) opp_printf("dev_realloc", "created [%p] old [%p]", new_ptr, ptr);
        if (ptr) {
            const size_t copy_size = std::min(current_size, new_size);
            opp_mem::dev_memcpy<T>(new_ptr, ptr, copy_size);
            opp_mem::dev_free<T>(ptr);
            current_size = new_size;
        }
        if (debug_mem) opp_printf("dev_realloc", "[%p]->[%p] cpy_count[%zu]", ptr, new_ptr, current_size);
        ptr = new_ptr;
    }

    // Resize device memory (only increasing size)
    template <typename T>
    inline static void dev_resize(T*& ptr, size_t& current_size, const size_t& new_size) {
        if (debug_mem) opp_printf("dev_resize", "[%p] %zu -> %zu", ptr, current_size, new_size);
        if (new_size > current_size) {
            opp_mem::dev_realloc<T>(ptr, current_size, new_size);
        }
    }

    // initialize device memory with a specific value
    template <typename T>
    inline static void dev_memset(T* ptr, size_t count, T value) {
        opp_queue->fill(ptr, value, count).wait();
    }

    // Allocate and initialize device memory with a specific value
    template <typename T>
    inline static T* dev_malloc_set(size_t count, T value) {
        T* ptr = opp_mem::dev_malloc<T>(count);
        opp_mem::dev_memset<T>(ptr, count, value);
        return ptr;
    }

    // Copy data from host to device, create new device arrays if requested
    template <typename T>
    inline static void copy_host_to_dev(T*& data_d, const T *data_h, size_t copy_count,
                        bool no_wait = false, bool create_new = false, size_t alloc_count = 0) {
        try {
            if (create_new) {
                if (data_d != nullptr)  
                    opp_mem::dev_free<T>(data_d);
                data_d = opp_mem::dev_malloc<T>(alloc_count);
            }
            opp_queue->memcpy(data_d, data_h, copy_count * sizeof(T));
            if (!no_wait) 
                opp_queue->wait(); 
            if (debug_mem) opp_printf("copy_host_to_dev", "[%p]->[%p] copy_count[%zu]", data_h, data_d, copy_count);
        }
        catch (const sycl::exception &e) {
            throw std::runtime_error(std::string("opp_mem::copy_host_to_dev: ") + e.what());
        }
    }

    // Copy data from device to host, no dot create new host arrays since it can be allocated differently, 
    // like malloc, new, stack, std::vector<>, hence free mechanism is unknown
    template <typename T>
    inline static void copy_dev_to_host(T* data_h, const T *data_d, size_t copy_count, 
                        bool no_wait = false) {
        try {
            opp_queue->memcpy(data_h, data_d, copy_count * sizeof(T));
            if (!no_wait) 
                opp_queue->wait();
            if (debug_mem) opp_printf("copy_dev_to_host", "[%p]->[%p] copy_count[%zu]", data_d, data_h, copy_count);
        }
        catch (const sycl::exception &e) {
            throw std::runtime_error(std::string("opp_mem::copy_dev_to_host: ") + e.what());
        }
    }

    // Copy data from device to device, create new device arrays if requested
    template <typename T>
    inline static void copy_dev_to_dev(T*& data_d, const T *data_h, size_t copy_count,
                        bool no_wait = false, bool create_new = false, size_t alloc_count = 0) {
        try {
            return opp_mem::copy_host_to_dev<T>(data_d, data_h, copy_count, no_wait, create_new, alloc_count);
        }
        catch (const sycl::exception &e) {
            throw std::runtime_error(std::string("opp_mem::copy_dev_to_dev: ") + e.what());
        }
    }
};

/*******************************************************************************/
// routines to resize constant/reduct arrays, if necessary

void opp_reallocReductArrays(int reduct_bytes);
void opp_mvReductArraysToDevice(int reduct_bytes);
void opp_mvReductArraysToHost(int reduct_bytes);

void opp_reallocConstArrays(int consts_bytes);
void opp_mvConstArraysToDevice(int consts_bytes);
void opp_mvConstArraysToHost(int consts_bytes);

/*******************************************************************************/
template <typename T>
void copy_from(
    const T* in_dat_d, T* out_dat_d, 
    const int* from_idx_map, 
    const int in_stride, const int out_stride, 
    const int in_offset, const int out_offset, 
    const int dim, const int size,
    const sycl::nd_item<1> &item) 
{
    const int tid = item.get_global_linear_id();
    if (tid < size) {
        const int idx = from_idx_map[tid];
        for (int d = 0; d < dim; d++) {
            out_dat_d[out_offset + tid + d * out_stride] = 
                            in_dat_d[in_offset + idx + d * in_stride];
        }
    }
}

/*******************************************************************************/
template <typename T>
void copy_from_to(
    const T* in_dat_d, T* out_dat_d, 
    const int* from_idx_map, const int* to_idx_map, 
    const int in_stride, const int out_stride, 
    const int in_offset, const int out_offset, 
    const int dim, const int size,
    const sycl::nd_item<1> &item) 
{
    const int tid = item.get_global_linear_id();
    if (tid < size) {
        const int f_idx = from_idx_map[tid];
        const int t_idx = to_idx_map[tid];
        for (int d = 0; d < dim; d++) {
            out_dat_d[out_offset + t_idx + d * out_stride] = 
                            in_dat_d[in_offset + f_idx + d * in_stride];
        }
    }
}

/*******************************************************************************/
/*
This function arranges the multi dimensional values in input array to output array 
according to the indices provided
    in_dat_dv - Input array with correct values
    out_dat_dv - Output array to have the sorted values
    new_idx_dv - indices of the input array to be arranged in the output array
    in_capacity - capacity of the input array (useful for multi dimensional dats)
    out_capacity - capacity of the output array (useful for multi dimensional dats)
    in_offset - start offset if the input array
    out_offset - start offset if the output array
    size - number of dat elements to arrange (usually this is the size of new_idx_dv)
    dim - dimension of the dat
*/
template <class T>
void copy_according_to_index(const T *in_dat_d, T *out_dat_d, const OPP_INT* from_idx_map,
        const int in_stride, const int out_stride, const int in_offset, const int out_offset, 
        const int size, const int dim)
{
    const int nblocks  = (size - 1) / opp_const_threads_per_block + 1;

    opp_queue->submit([&](sycl::handler &cgh) {
        cgh.parallel_for(
            sycl::nd_range<1>(opp_const_threads_per_block * nblocks, opp_const_threads_per_block),
            [=](sycl::nd_item<1> item) {
                copy_from<T>(
                    in_dat_d, out_dat_d,
                    from_idx_map,
                    in_stride, out_stride, 
                    in_offset, out_offset,
                    dim, size, 
                    item);
            });
    });
}

template <class T>
void copy_according_to_index(const T *in_dat_dv, T *out_dat_dv, const OPP_INT* new_idx_dv,
        int in_capacity, int out_capacity, int size, int dim)
{
    copy_according_to_index(in_dat_dv, out_dat_dv, new_idx_dv, in_capacity, out_capacity, 
        0, 0, size, dim);
}

/*******************************************************************************/
template <class T>
void opp_register_const(T*& ptr, const size_t count) {
    if (ptr == nullptr) {
        ptr = opp_mem::dev_malloc<T>(count);
        opp_consts.push_back((char*)ptr);
    }
}

/*******************************************************************************/
template <typename T>
T opp_atomic_fetch_add(T* address, T value) {
    return dpct::atomic_fetch_add<sycl::access::address_space::generic_space>(address, value);
}

/*******************************************************************************/
inline void opp_set_stride(OPP_INT*& data_d, OPP_INT& data_h, OPP_INT new_data) {
    opp_register_const<OPP_INT>(data_d, 1);
    if (data_h != new_data) {
        data_h = new_data;
        opp_mem::copy_host_to_dev<OPP_INT>(data_d, &data_h, 1);
    }    
}

/*******************************************************************************/
template <typename T>
inline void write_array_to_file(const T* array, size_t size, const std::string& filename, bool is_host = true) {
    const T* internal_array = nullptr;
    std::vector<T> host_vec;
    if (!is_host) {
        host_vec.resize(size);
        opp_mem::copy_dev_to_host<T>(host_vec.data(), array, size);
        internal_array = host_vec.data();
    }
    else {
        internal_array = array;
    }
    const std::string modified_file_name = filename + "_sycl_r"+ std::to_string(OPP_rank) + 
                                                "_i" + std::to_string(OPP_main_loop_iter);
    std::ofstream outFile(modified_file_name);
    if (!outFile) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    if constexpr (std::is_same<T, double>::value)
        outFile << std::setprecision(25);
    outFile << size << " 1 -- 0 0\n";
    for (int i = 0; i < size; ++i) {
        outFile << " " << internal_array[i] << "\n";
    }
    outFile.close();
}

/*******************************************************************************/
// in opp_direct_hop_sycl.cpp
void opp_init_dh_device(opp_set set); 
void opp_gather_dh_move_indices(opp_set set); // gathers all device global move info to the global mover for communication

//*******************************************************************************
template <typename T>
inline T** opp_create_thread_level_data(opp_arg arg) 
{
    opp_dat dat = arg.dat;
    const int array_count = opp_params->get<OPP_INT>("gpu_reduction_arrays");
    const int array_size = (dat->set->size + dat->set->exec_size + dat->set->nonexec_size) * dat->dim;
    const int num_blocks = (array_size - 1) / OPP_gpu_threads_per_block + 1;
    if (OPP_main_loop_iter == 0) {
        dat->thread_data->resize(array_count);
        dat->thread_data->at(0) = dat->data_d;
        if (array_count > 1)
            dat->thread_data->at(1) = dat->data_swap_d;
        for (int i = 2; i < array_count; ++i) {
            dat->thread_data->at(i) = (char*)opp_mem::dev_malloc<T>(array_size);
        }
        dat->thread_data_d = (char**)opp_mem::dev_malloc<T*>(array_count);
        opp_queue->memcpy(dat->thread_data_d, dat->thread_data->data(), array_count * sizeof(char*)).wait();
    }
    else if (dat->set->is_particle) {
        dat->thread_data->at(0) = dat->data_d;  // Move loop can swap device arrays
        if (array_count > 1)
            dat->thread_data->at(1) = dat->data_swap_d;
        opp_queue->memcpy(dat->thread_data_d, dat->thread_data->data(), array_count * sizeof(char*)).wait();
    }

    opp_queue->submit([&](sycl::handler &cgh) {

        T** arrays_d = (T**)dat->thread_data_d;

        auto kernel = [=](sycl::nd_item<1> item) {
            const int idx = item.get_global_linear_id();
            if (idx < array_size) {
                for (int i = 1; i < array_count; i++) {
                    arrays_d[i][idx] = 0.0;
                }
            }
        };

        cgh.parallel_for<class opp_create_thread_level_data>(
            sycl::nd_range<1>(OPP_gpu_threads_per_block * num_blocks, OPP_gpu_threads_per_block), kernel);
    });

    OPP_DEVICE_SYNCHRONIZE();  

    return (T**)dat->thread_data_d;
}

//*******************************************************************************
template <typename T>
inline void opp_reduce_thread_level_data(opp_arg arg) 
{
    OPP_DEVICE_SYNCHRONIZE(); 
    
    opp_queue->submit([&](sycl::handler &cgh) {

        opp_dat dat = arg.dat;
        const int array_count = opp_params->get<OPP_INT>("gpu_reduction_arrays");
        const int array_size = (dat->set->size + dat->set->exec_size + dat->set->nonexec_size) * dat->dim;        
        const int num_blocks = (array_size - 1) / OPP_gpu_threads_per_block + 1;

        T** arrays_d = (T**)dat->thread_data_d;
        
        auto kernel = [=](sycl::nd_item<1> item) {
            const int idx = item.get_global_linear_id();
            if (idx < array_size) {
                for (int j = 1; j < array_count; ++j) {  // Start from the second array
                    arrays_d[0][idx] += arrays_d[j][idx];
                }
            }
        };

        cgh.parallel_for<class opp_reduce_thread_level_data>(
            sycl::nd_range<1>(OPP_gpu_threads_per_block * num_blocks, OPP_gpu_threads_per_block), kernel);
    });

    OPP_DEVICE_SYNCHRONIZE();  
}


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

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>

#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/copy.h>
#include <thrust/sequence.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/generate.h>
#include <thrust/random.h>

#ifdef USE_MPI
    #include <opp_mpi_core.h>
#endif

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

#define OPP_DEVICE_SYNCHRONIZE() \
    do { \
        /*cudaError_t err = cudaStreamSynchronize(*opp_stream); */ cudaError_t err = cudaDeviceSynchronize();\
        if (cudaSuccess != err) { \
            std::string log = std::string(__FILE__) + "(" + std::to_string(__LINE__) + \
                                std::string(") Error : ") + cudaGetErrorString(err); \
            opp_abort(log.c_str()); \
        } \
    } while (0)

#define OPP_DEV_CHECK(err) \
    do { \
        if (cudaSuccess != err) { \
            std::string log = std::string(__FILE__) + "(" + std::to_string(__LINE__) + \
                                std::string(") Error : ") + cudaGetErrorString(err); \
            opp_abort(log.c_str()); \
        } \
    } while (0)

#define OPP_DEVICE_GLOBAL_LINEAR_ID (blockIdx.x * blockDim.x + threadIdx.x)
#define OPP_GLOBAL_FUNCTION __global__ 
#define OPP_DEVICE_FUNCTION __device__ 
#define ADDITIONAL_PARAMETERS 
#define OPP_ATOMIC_FETCH_ADD(address, value) atomicAdd(address, value)

/*******************************************************************************/
extern int* opp_saved_mesh_relation_d;
extern size_t opp_saved_mesh_relation_size;

// device vectors in opp_particle_sorter.cu
extern thrust::device_vector<OPP_INT> ps_from_indices_dv;
extern thrust::device_vector<OPP_INT> hf_from_indices_dv;
extern thrust::device_vector<OPP_INT> hf_sequence_dv;

extern OPP_INT *OPP_move_particle_indices_d;
extern OPP_INT *OPP_move_cell_indices_d;
extern OPP_INT *OPP_move_count_d;
extern thrust::device_vector<OPP_INT> OPP_move_particle_indices_dv;
extern thrust::device_vector<OPP_INT> OPP_move_cell_indices_dv;

extern OPP_INT *OPP_remove_particle_indices_d;
extern thrust::device_vector<OPP_INT> OPP_remove_particle_indices_dv;

extern std::map<int, thrust::host_vector<OPP_INT>> cell_indices_hv;     // cellid in the foreign rank, arrange according to rank
extern std::map<int, thrust::host_vector<OPP_INT>> particle_indices_hv; // particle ids to send, arrange according to rank
extern std::map<int, thrust::device_vector<OPP_INT>> particle_indices_dv;
extern std::map<int, thrust::device_vector<char>> send_data;
extern std::map<int, thrust::device_vector<char>> recv_data;

// arrays for global constants and reductions
extern OPP_INT OPP_consts_bytes;
extern OPP_INT OPP_reduct_bytes;
extern char *OPP_reduct_h, *OPP_reduct_d;
extern char *OPP_consts_h, *OPP_consts_d;

extern char opp_move_status_flag;
extern bool opp_move_hop_iter_one_flag;
extern OPP_INT* opp_p2c;
extern OPP_INT* opp_c2c;

extern opp_dh_indices dh_indices_d;
extern opp_dh_indices dh_indices_h;

//*************************************************************************************************
void opp_device_exit();
void opp_device_init(int argc, char **argv);

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
// routines to resize constant/reduct arrays, if necessary

void opp_reallocReductArrays(int reduct_bytes);
void opp_mvReductArraysToDevice(int reduct_bytes);
void opp_mvReductArraysToHost(int reduct_bytes);

void opp_reallocConstArrays(int consts_bytes);
void opp_mvConstArraysToDevice(int consts_bytes);
void opp_mvConstArraysToHost(int consts_bytes);

/*******************************************************************************/
template <typename T>
T* opp_get_dev_raw_ptr(thrust::device_vector<T>& dv) {
    return (T*)thrust::raw_pointer_cast(dv.data());
}
template <typename T>
const T* opp_get_dev_raw_ptr(const thrust::device_vector<T>& dv) {
    return (const T*)thrust::raw_pointer_cast(dv.data());
}

/*******************************************************************************/
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
    const thrust::device_ptr<const OPP_INT> new_idx_dp, int in_capacity, int out_capacity, int in_offset, int out_offset, 
    int size, int dimension)
{
    switch (dimension) {
        case 1:
            thrust::copy_n(thrust::make_permutation_iterator(
                thrust::make_zip_iterator(
                    thrust::make_tuple(in_dat_dv->begin() + in_offset)
                ), 
                new_idx_dp), 
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
                new_idx_dp), 
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
                new_idx_dp), 
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
                new_idx_dp), 
                size, 
                thrust::make_zip_iterator(
                    thrust::make_tuple(
                        out_dat_dv->begin() + out_offset, 
                        (out_dat_dv->begin() + out_offset + out_capacity), 
                        (out_dat_dv->begin() + out_offset + (2 * out_capacity)), 
                        (out_dat_dv->begin() + out_offset + (3 * out_capacity)))));
            break;
        default:
            std::cerr << "copy_according_to_index not implemented for dim " << dimension << std::endl;
            exit(-1);
    }
}

template <class T> 
void copy_according_to_index(thrust::device_vector<T>* in_dat_dv, thrust::device_vector<T>* out_dat_dv, 
    const thrust::device_vector<int>& new_idx_dv, int in_capacity, int out_capacity, 
    int in_offset, int out_offset, int size, int dimension)
{
    copy_according_to_index<T>(in_dat_dv, out_dat_dv, new_idx_dv.data(), 
        in_capacity, out_capacity, in_offset, out_offset, size, dimension);
}

template <class T> 
void copy_according_to_index(thrust::device_vector<T>* in_dat_dv, thrust::device_vector<T>* out_dat_dv, 
    const thrust::device_vector<int>& new_idx_dv, int in_capacity, int out_capacity, int size, int dimension)
{
    copy_according_to_index<T>(in_dat_dv, out_dat_dv, new_idx_dv.data(), 
        in_capacity, out_capacity, 0, 0, size, dimension);
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
        T* ptr;
        cudaError_t err = cudaMalloc((void**)&ptr, count * sizeof(T));
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string("dev_malloc: ") + cudaGetErrorString(err));
        }
        return ptr;
    }

    // Free device memory
    template <typename T>
    inline static void dev_free(T*& ptr) {
        if (ptr) {
            cudaError_t err = cudaFree(ptr);
            if (err != cudaSuccess) {
                throw std::runtime_error(std::string("dev_free: ") + cudaGetErrorString(err));
            }
        }
        ptr = nullptr;
    }

    // Copy memory from one device pointer to another
    template <typename T>
    inline static void dev_memcpy(T* dst, const T* src, size_t copy_count) {
        cudaError_t err = cudaMemcpy(dst, src, copy_count * sizeof(T), cudaMemcpyDeviceToDevice);
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string("dev_memcpy: ") + cudaGetErrorString(err));
        }
    }

    // Resize device memory
    template <typename T>
    inline static void dev_realloc(T*& ptr, size_t& current_size, const size_t& new_size) {
        if (new_size <= 0) 
            throw std::runtime_error("dev_realloc: New Realloc size invalid - " + std::to_string(new_size));
        T* new_ptr = opp_mem::dev_malloc<T>(new_size);
        if (ptr) {
            const size_t copy_size = std::min(current_size, new_size);
            opp_mem::dev_memcpy<T>(new_ptr, ptr, copy_size);
            opp_mem::dev_free<T>(ptr);
            current_size = new_size;
        }
        ptr = new_ptr;
    }

    // Resize device memory (only increasing size)
    template <typename T>
    inline static void dev_resize(T*& ptr, size_t& current_size, const size_t& new_size) {
        if (new_size > current_size) {
            opp_mem::dev_realloc<T>(ptr, current_size, new_size);
        }
    }

    // Initialize device memory with a specific value
    template <typename T>
    inline static void dev_memset(T* ptr, size_t count, T value) {
        cudaError_t err = cudaMemset(ptr, value, count * sizeof(T));
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string("dev_memset: ") + cudaGetErrorString(err));
        }
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
        if (create_new) {
            if (data_d != nullptr)  
                opp_mem::dev_free<T>(data_d);
            data_d = opp_mem::dev_malloc<T>(alloc_count);
        }
        cudaError_t err = cudaMemcpy(data_d, data_h, copy_count * sizeof(T), cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string("copy_host_to_dev: ") + cudaGetErrorString(err));
        }
    }

    // Copy data from device to host
    template <typename T>
    inline static void copy_dev_to_host(T* data_h, const T *data_d, size_t copy_count) {
        cudaError_t err = cudaMemcpy(data_h, data_d, copy_count * sizeof(T), cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string("copy_dev_to_host: ") + cudaGetErrorString(err));
        }
    }

    // Copy data from device to device, create new device arrays if requested
    template <typename T>
    inline static void copy_dev_to_dev(T*& data_d, const T *data_h, size_t copy_count,
                        bool no_wait = false, bool create_new = false, size_t alloc_count = 0) {
        if (create_new) {
            if (data_d != nullptr)  
                opp_mem::dev_free<T>(data_d);
            data_d = opp_mem::dev_malloc<T>(alloc_count);
        }        
        cudaError_t err = cudaMemcpy(data_d, data_h, copy_count * sizeof(T), cudaMemcpyDeviceToDevice);
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string("copy_dev_to_dev: ") + cudaGetErrorString(err));
        }
    }

    // Copy host data to device symbol
    template <typename T>
    inline static void dev_copy_to_symbol(const T& symbol, T* data_h, const OPP_INT* new_data, 
                                            size_t copy_count) {    
        bool copy_symbol = false;
        for (int i = 0; i < copy_count; i++) {
            if (data_h[i] != new_data[i])
                copy_symbol = true;
        }
        if (copy_symbol) {
            memcpy(data_h, new_data, copy_count * sizeof(T));  
            cudaError_t err = cudaMemcpyToSymbol(symbol, new_data, copy_count * sizeof(T));
            if (err != cudaSuccess) {
                throw std::runtime_error(std::string("dev_copy_to_symbol: ") + cudaGetErrorString(err));
            }
        }
    }
    template <typename T>
    inline static void dev_copy_to_symbol(const T& symbol, const T* data_h, size_t copy_count) {    
        cudaError_t err = cudaMemcpyToSymbol(symbol, data_h, copy_count * sizeof(T));
        if (err != cudaSuccess) {
            throw std::runtime_error(std::string("dev_copy_to_symbol: ") + cudaGetErrorString(err));
        }
    }
};

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
    const std::string modified_file_name = filename + "_cuda_r" + std::to_string(OPP_rank) + 
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
// in opp_direct_hop_cuda.cu
void opp_init_dh_device(opp_set set); 
void opp_gather_dh_move_indices(opp_set set); // gathers all device global move info to the global mover for communication

//*******************************************************************************
template <typename T>
__global__ void zero_set_array_except_first_kernel(T** __restrict__ arrays, int num_arrays, int size) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        for (int i = 1; i < num_arrays; i++) {
            arrays[i][idx] = 0.0;
        }
    }
}
template <typename T>
inline T** opp_create_thread_level_data(opp_arg arg) 
{
    opp_dat dat = arg.dat;
    const int array_count = opp_params->get<OPP_INT>("gpu_reduction_arrays");
    const int array_size = (dat->set->size + dat->set->exec_size + dat->set->nonexec_size) * dat->dim;
    if (OPP_main_loop_iter == 0) {
        dat->thread_data->resize(array_count);
        dat->thread_data->at(0) = dat->data_d;
        if (array_count > 1)
            dat->thread_data->at(1) = dat->data_swap_d;
        for (int i = 2; i < array_count; ++i) {
            dat->thread_data->at(i) = (char*)opp_mem::dev_malloc<T>(array_size);
        }
        dat->thread_data_d = (char**)opp_mem::dev_malloc<T*>(array_count);
        cudaMemcpy(dat->thread_data_d, dat->thread_data->data(), array_count * sizeof(char*), cudaMemcpyHostToDevice);
    }
    else if (dat->set->is_particle) {
        dat->thread_data->at(0) = dat->data_d; // Why? because move loop can swap device arays
        if (array_count > 1)
            dat->thread_data->at(1) = dat->data_swap_d;
        cudaMemcpy(dat->thread_data_d, dat->thread_data->data(), array_count * sizeof(char*), cudaMemcpyHostToDevice);
    }

    const int num_blocks = (array_size - 1) / OPP_gpu_threads_per_block + 1;
    zero_set_array_except_first_kernel<<<num_blocks, OPP_gpu_threads_per_block>>>(
        (T**)dat->thread_data_d, array_count, array_size);
    OPP_DEVICE_SYNCHRONIZE();

    return (T**)dat->thread_data_d;
}

//*******************************************************************************
template <typename T>
__global__ void reduce_arrays_kernel(T** __restrict__ arrays, const int num_arrays, const int size) 
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < size) {
        for (int i = 1; i < num_arrays; ++i) {  // Start from the second array
            arrays[0][idx] += arrays[i][idx];
        }
    }
}
template <typename T>
inline void opp_reduce_thread_level_data(opp_arg arg) 
{
    opp_dat dat = arg.dat;
    const int array_count = opp_params->get<OPP_INT>("gpu_reduction_arrays");
    const int array_size = (dat->set->size + dat->set->exec_size + dat->set->nonexec_size) * dat->dim;
    const int num_blocks = (array_size - 1) / OPP_gpu_threads_per_block + 1;
    reduce_arrays_kernel<<<num_blocks, OPP_gpu_threads_per_block>>>(
        (T**)dat->thread_data_d, array_count, array_size);
    OPP_DEVICE_SYNCHRONIZE();    
}
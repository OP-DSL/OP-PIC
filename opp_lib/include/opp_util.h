
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

#ifndef MIN
#define MIN(a, b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b) ((a > b) ? (a) : (b))
#endif

#define ROUND_UP(bytes) (((bytes) + 15) & ~15)
#define ROUND_UP_64(bytes) (((bytes) + 63) & ~63)

//********************************************************************************
std::string get_time_str();

//********************************************************************************
char *copy_str(char const *src);

//********************************************************************************
std::vector<size_t> sort_indexes(const int* cell_indices, int size);

//********************************************************************************

int min(int array[], int size);

int binary_search(int a[], int value, int low, int high);

int linear_search(int a[], int value, int low, int high);

void quickSort(int arr[], int left, int right);

void quickSort_2(int arr1[], int arr2[], int left, int right);

void quickSort_dat(int arr[], char dat[], int left, int right, int elem_size);

void quickSort_map(int arr[], int map[], int left, int right, int dim);

int removeDups(int a[], int array_size);

int compare_sets(opp_set set1, opp_set set2);

void op_timers(double *cpu, double *et); 

template <typename T>
std::string str(T value, std::string format) {
    char buffer[32];
    //int numChars = snprintf(buffer, sizeof(buffer), format.c_str(), value);
    snprintf(buffer, sizeof(buffer), format.c_str(), value);
    return std::string(buffer);
}

double rnd();

void reset_seed();

double* get_dandom_distriution(int count, int dim);

int file_exist(char const *filename);

bool opp_type_equivalence(const char *a, const char *b);

void opp_compress_write(const std::string &filename, 
                        const int* data, const size_t count);
void opp_decompress_read(const std::string &filename, size_t originalSize, 
                        int* data);
   
//*************************************************************************************************
template <typename RNG>
inline std::vector<std::vector<double>>
    get_normal_distribution(const int N, const int ndim, const double mu, const double sigma, RNG &rng) {

    std::normal_distribution<> d{mu, sigma};
    std::vector<std::vector<double>> array(ndim);
    for (int dimx = 0; dimx < ndim; dimx++) {
        array[dimx] = std::vector<double>(N);
        for (int px = 0; px < N; px++) {
            array[dimx][px] = d(rng);
        }
    }

    return array;
}

//*************************************************************************************************
template <typename RNG>
inline std::vector<std::vector<double>>
    get_uniform_within_extents(const int N, const int ndim, const double *extents, RNG &rng) {

    std::uniform_real_distribution<double> uniform_rng(0.0, 1.0);
    std::vector<std::vector<double>> positions(ndim);

    for (int dimx = 0; dimx < ndim; dimx++) {
        positions[dimx] = std::vector<double>(N);
        const double ex = extents[dimx];
        for (int px = 0; px < N; px++) {
            positions[dimx][px] = ex * uniform_rng(rng);
        }
    }

    return positions;
}

//*************************************************************************************************
template <typename T>
inline void uniform_within_cartesian_cells(int ndim, const double* extents, const T* cell_pos_ll, 
    const int64_t cells_set_size, const int64_t npart_per_cell, std::vector<std::vector<double>> &positions, 
    std::vector<int> &cells, std::mt19937 rng) {

    const int64_t npart_total = npart_per_cell * cells_set_size;
    
    cells.resize(npart_total);

    positions.resize(ndim);
    for (int dx = 0; dx < ndim; dx++) {
        positions[dx] = std::vector<double>(npart_total);
    }

    for (int cx = 0; cx < cells_set_size; cx++) {

        const int index_start = cx * npart_per_cell;
        const int index_end = (cx + 1) * npart_per_cell;

        auto positions_ref_cell = get_uniform_within_extents(npart_per_cell, ndim, extents, rng);

        int index = 0;
        for (int ex = index_start; ex < index_end; ex++) {

            cells.at(ex) = cx;

            for (int dx = 0; dx < ndim; dx++) {
                positions.at(dx).at(ex) =
                    cell_pos_ll[cx * ndim + dx] + positions_ref_cell.at(dx).at(index);
            }
            index++;
        }
    }
}

//*************************************************************************************************
template <typename T>
inline std::vector<T> reverse_argsort(const std::vector<T> &array) 
{
    std::vector<T> indices(array.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
                [&array](int left, int right) -> bool {
                    return array[left] > array[right];
                });

    return indices;
}

//*************************************************************************************************
template <typename T>
inline void get_decomp_1d(const T N_compute_units, const T N_work_items,
                   const T work_unit, T *rstart, T *rend) 
{
    const auto pq = std::div(N_work_items, N_compute_units);
    const T i = work_unit;
    const T p = pq.quot;
    const T q = pq.rem;
    const T n = (i < q) ? (p + 1) : p;
    const T start = (MIN(i, q) * (p + 1)) + ((i > q) ? (i - q) * p : 0);
    const T end = start + n;

    *rstart = start;
    *rend = end;
}

//*************************************************************************************************
template <typename T>
std::vector<T> sort_iota_by_key(const T* keys, size_t size) 
{ 
    std::vector<T> idx(size);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(),
        [keys](T i1, T i2) { return keys[i1] < keys[i2]; });

    return idx;
}

//*************************************************************************************************
template <typename T>
std::vector<T> shuffle_iota_by_key(const T* keys, size_t size) 
{
    std::vector<T> idx(size);
    std::iota(idx.begin(), idx.end(), 0);

    // Partition the indices based on whether the key is MAX_CELL_INDEX
    auto partition_point = std::partition(idx.begin(), idx.end(), 
        [keys](T i) { return keys[i] != MAX_CELL_INDEX; });
    idx.resize(partition_point - idx.begin());

    // TODO : Since relative order is not preserved in std::partition, 
    // we might be able to remove below shuffle
    // std::random_device rd;
    // std::mt19937 g(rd());
    // std::shuffle(idx.begin(), idx.end(), g);

    return idx;
}

//*************************************************************************************************
inline void getDatTypeSize(opp_data_type dtype, std::string& type, int& size)
{
    if (dtype == DT_REAL) {
        type = "double";
        size = sizeof(OPP_REAL);
    }
    else if (dtype == DT_INT) {
        type = "int";
        size = sizeof(OPP_INT);       
    }
    else {
        std::cerr << "Data type in Dat not supported" << std::endl;
    }
}

//*************************************************************************************************
inline void opp_printf(const char* function, const char *format, ...)
{
    char buf[LOG_STR_LEN];
    va_list args;
    va_start(args, format);
    vsnprintf(buf, LOG_STR_LEN, format, args);
    va_end(args);

    printf("%s[%d][%d] - %s\n", function, OPP_rank, OPP_main_loop_iter, buf);
    fflush(stdout);
}

//*************************************************************************************************

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

#include <opp_util.h>
#include <opp_lib_core.h>
#include <chrono>
#include <numeric>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <random>

//********************************************************************************
char *copy_str(char const *src) 
{
    size_t src_len = strlen(src); // Calculate the actual length of src
    size_t dest_len = (src_len > 100) ? 100 : src_len; // Limit the destination length to 100 characters

    char *dest = (char *)opp_host_malloc((dest_len + 1) * sizeof(char));
    if (dest) {
        memcpy(dest, src, dest_len);
        dest[dest_len] = '\0'; // Ensure the destination string is null-terminated
    }
    return dest;
}

//********************************************************************************
std::string getTimeStr()
{
    std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::string s(30, '\0');
    std::strftime(&s[0], s.size(), "%Y_%m_%d__%H_%M_%S", std::localtime(&now));
    return s;
}

//********************************************************************************
std::vector<size_t> sort_indexes(const int* cell_indices, int size) 
{ 

    std::vector<size_t> idx(size);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(), // std::stable_sort(idx.begin(), idx.end(), 
        [cell_indices](size_t i1, size_t i2) 
        { 
            return cell_indices[i1] < cell_indices[i2];
        });

    return idx;
}

/*******************************************************************************
* Return the index of the min value in an array
*******************************************************************************/

int min(int array[], int size) 
{
    int min = 99; // initialized to 99 .. should check op_mpi_part_core and fix
    int index = -1;
    for (int i = 0; i < size; i++) 
    {
        if (array[i] < min) 
        {
            index = i;
            min = array[i];
        }
    }
    return index;
}

/*******************************************************************************
* Binary search an array for a given value
*******************************************************************************/

int binary_search(int a[], int value, int low, int high) 
{
    if (high < low)
        return -1; // not found
    else if (high == low) 
    {
        if (a[low] == value)
            return low;
        else
            return -1;
    } 
    else if (high == (low + 1)) 
    {
        if (a[low] == value)
            return low;
        else if (a[high] == value)
            return high;
        else
            return -1;
    }

    int mid = low + (high - low) / 2;
    if (a[mid] > value)
        return binary_search(a, value, low, mid - 1);
    else if (a[mid] < value)
        return binary_search(a, value, mid + 1, high);
    else
        return mid; // found
}

/*******************************************************************************
* Linear search an array for a given value
*******************************************************************************/

int linear_search(int a[], int value, int low, int high) 
{
    for (int i = low; i <= high; i++) {
        if (a[i] == value)
        return i;
    }
    return -1;
}

/*******************************************************************************
* Quicksort an array
*******************************************************************************/

void quickSort(int arr[], int left, int right) 
{
    int i = left;
    int j = right;
    int tmp;
    if (left==right) return;
    int pivot = arr[(left + right) / 2];

    // partition
    while (i <= j) 
    {
        while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;
        if (i <= j) 
        {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
            i++;
            j--;
        }
    };
    // recursion
    if (left < j)
        quickSort(arr, left, j);
    if (i < right)
        quickSort(arr, i, right);
}

/*******************************************************************************
* Quick sort arr1 and organise arr2 elements according to the sorted arr1 order
*******************************************************************************/

void quickSort_2(int arr1[], int arr2[], int left, int right) 
{
    int i = left;
    int j = right;
    int tmp1, tmp2;
    if (left==right) 
        return;
    int pivot = arr1[(left + right) / 2];

    // partition
    while (i <= j) 
    {
        while (arr1[i] < pivot)
            i++;
        while (arr1[j] > pivot)
            j--;
        if (i <= j) 
        {
            tmp1 = arr1[i];
            arr1[i] = arr1[j];
            arr1[j] = tmp1;

            tmp2 = arr2[i];
            arr2[i] = arr2[j];
            arr2[j] = tmp2;
            i++;
            j--;
        }
    };
    // recursion
    if (left < j)
        quickSort_2(arr1, arr2, left, j);
    if (i < right)
        quickSort_2(arr1, arr2, i, right);
}

/*******************************************************************************
* Quick sort arr and organise dat[] elements according to the sorted arr order
*******************************************************************************/

void quickSort_dat(int arr[], char dat[], int left, int right, int elem_size2) 
{
    if (left < 0 || right <= 0)
        return;
    if (left==right) return;
    size_t elem_size = elem_size2;
    int i = left, j = right;
    int tmp;
    char *tmp_dat = (char *)opp_host_malloc(sizeof(char) * elem_size);
    int pivot = arr[(left + right) / 2];

    // partition
    while (i <= j) 
    {
        while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;

        if (i < j) 
        {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;

            // tmp_dat = dat[i];
            memcpy(tmp_dat, (void *)&dat[i * elem_size], elem_size);
            // dat[i] = dat[j];
            memcpy(&dat[i * elem_size], (void *)&dat[j * elem_size], elem_size);
            // dat[j] = tmp_dat;
            memcpy(&dat[j * elem_size], (void *)tmp_dat, elem_size);
            i++;
            j--;
        } 
        else if (i == j) 
        {
            i++;
            j--;
        }
    };

    // recursion
    if (left < j)
        quickSort_dat(arr, dat, left, j, elem_size);
    if (i < right)
        quickSort_dat(arr, dat, i, right, elem_size);
    
    opp_host_free(tmp_dat);
}

/*******************************************************************************
* Quick sort arr and organise map[] elements according to the sorted arr order
*******************************************************************************/

void quickSort_map(int arr[], int map[], int left, int right, int dim) 
{
    if (left==right) 
        return;

    int i = left, j = right;
    int tmp;
    int *tmp_map = (int *)opp_host_malloc(sizeof(int) * dim);
    int pivot = arr[(left + right) / 2];

    // partition
    while (i <= j) 
    {
        while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;

        if (i < j) 
        {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;

            // tmp_dat = dat[i];
            memcpy(tmp_map, (void *)&map[i * dim], dim * sizeof(int));
            // dat[i] = dat[j];
            memcpy(&map[i * dim], (void *)&map[j * dim], dim * sizeof(int));
            // dat[j] = tmp_dat;
            memcpy(&map[j * dim], (void *)tmp_map, dim * sizeof(int));
            i++;
            j--;
        } else if (i == j) 
        {
            i++;
            j--;
        }
    };

    // recursion
    if (left < j)
        quickSort_map(arr, map, left, j, dim);
    if (i < right)
        quickSort_map(arr, map, i, right, dim);
    
    opp_host_free(tmp_map);
}

/*******************************************************************************
* Remove duplicates in an array
*******************************************************************************/

int removeDups(int a[], int array_size) 
{
    int i, j;
    j = 0;
    // Remove the duplicates ...
    for (i = 1; i < array_size; i++) 
    {
        if (a[i] != a[j]) 
        {
            j++;
            a[j] = a[i]; // Move it to the front
        }
    }
    // The new array size..
    array_size = (j + 1);
    return array_size;
}

//********************************************************************************
int compare_sets(opp_set set1, opp_set set2) 
{
    if (set1->size == set2->size && strcmp(set1->name, set2->name) == 0)
        return 1;
    else
        return 0;
}

void op_timers(double *cpu, double *et) 
{
    (void)cpu;
    struct timeval t;
    gettimeofday(&t, (struct timezone *)0);
    *et = t.tv_sec + t.tv_usec * 1.0e-6;
}

//********************************************************************************
double* get_dandom_distriution(int count, int dim)
{
    double *dist = new double[count * dim];

    for (int i = 0; i < count * dim; i++)
    {   
        dist[i] = rnd();
    }

    return dist;
}
//********************************************************************************
std::mt19937 mt_gen(0);        /*seed*/
std::uniform_real_distribution<double> rnd_dist(0, 1.0);

double rnd() 
{
    return rnd_dist(mt_gen);
}

//********************************************************************************
void reset_seed() 
{ 
    mt_gen.seed(0); 
}

//********************************************************************************
/*******************************************************************************
* Check if a file exists
*******************************************************************************/
int file_exist(char const *filename) {
    struct stat buffer;
    return (stat(filename, &buffer) == 0);
}

const char *doubles[] = {"double", "double:soa", "real(8)", "double precision"};
const char *floats[] = {"float", "float:soa", "real(4)", "real"};
const char *ints[] = {"int", "int:soa", "integer(4)", "integer"};

bool opp_type_equivalence(const char *a, const char *b) {

    for (int i = 0; i < 4; i++) {
        if (strcmp(a, doubles[i]) == 0) {
            for (int j = 0; j < 4; j++) {
                if (strcmp(b, doubles[j]) == 0) {
                return true;
                }
            }
        }
    }
    for (int i = 0; i < 4; i++) {
        if (strcmp(a, floats[i]) == 0) {
            for (int j = 0; j < 4; j++) {
                if (strcmp(b, floats[j]) == 0) {
                    return true;
                }
            }
        }
    }
    for (int i = 0; i < 4; i++) {
        if (strcmp(a, ints[i]) == 0) {
            for (int j = 0; j < 4; j++) {
                if (strcmp(b, ints[j]) == 0) {
                    return true;
                }
            }
        }
    }
    return false;
}

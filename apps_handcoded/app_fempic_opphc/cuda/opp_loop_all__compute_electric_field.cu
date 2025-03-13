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

//*********************************************
// AUTO GENERATED CODE
//*********************************************


__constant__ int computeEF_stride_OPP_CUDA_0;
__constant__ int computeEF_stride_OPP_CUDA_1;
__constant__ int computeEF_stride_OPP_CUDA_2_MAP;

int computeEF_stride_OPP_HOST_0 = -1;
int computeEF_stride_OPP_HOST_1 = -1;
int computeEF_stride_OPP_HOST_2_MAP = -1;

//user function
//*************************************************************************************************
__device__ void compute_electric_field__kernel_gpu(
    double *cell_ef,             
    const double *cell_shape_deriv,
    const double *node_potential0,
    const double *node_potential1,
    const double *node_potential2,
    const double *node_potential3
)
{
    for (int dim = 0; dim < DIM; dim++)
    {
        cell_ef[dim * computeEF_stride_OPP_CUDA_0] -= 
            (cell_shape_deriv[(0 * DIM + dim) * computeEF_stride_OPP_CUDA_1] * (*node_potential0));

        cell_ef[dim * computeEF_stride_OPP_CUDA_0] -= 
            (cell_shape_deriv[(1 * DIM + dim) * computeEF_stride_OPP_CUDA_1] * (*node_potential1));

        cell_ef[dim * computeEF_stride_OPP_CUDA_0] -= 
            (cell_shape_deriv[(2 * DIM + dim) * computeEF_stride_OPP_CUDA_1] * (*node_potential2));

        cell_ef[dim * computeEF_stride_OPP_CUDA_0] -= 
            (cell_shape_deriv[(3 * DIM + dim) * computeEF_stride_OPP_CUDA_1] * (*node_potential3));
    }    
}

// CUDA kernel function
//*************************************************************************************************
__global__ void opp_cuda_ComputeElectricField(
    double *__restrict dir_arg0,
    const double *__restrict dir_arg1,
    const double *__restrict ind_arg2,
    const int *__restrict opDat2Map,
    const double *__restrict ind_arg3,
    const double *__restrict ind_arg4,
    const double *__restrict ind_arg5,
    int start,
    int end) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        const int map1idx = opDat2Map[n + computeEF_stride_OPP_CUDA_2_MAP * 0];
        const int map2idx = opDat2Map[n + computeEF_stride_OPP_CUDA_2_MAP * 1];
        const int map3idx = opDat2Map[n + computeEF_stride_OPP_CUDA_2_MAP * 2];
        const int map4idx = opDat2Map[n + computeEF_stride_OPP_CUDA_2_MAP * 3];

        //user-supplied kernel call
        compute_electric_field__kernel_gpu(
            (dir_arg0 + n),
            (dir_arg1 + n),
            (ind_arg2 + map1idx),
            (ind_arg2 + map2idx),
            (ind_arg2 + map3idx),
            (ind_arg2 + map4idx)
        );
    }
}

//*************************************************************************************************
void opp_loop_all__compute_electric_field(
    opp_set set,      // cells_set
    opp_arg arg0,     // cell_ef,
    opp_arg arg1,     // cell_shape_deriv,
    opp_arg arg2,     // node_potential0,
    opp_arg arg3,     // node_potential1,
    opp_arg arg4,     // node_potential2,
    opp_arg arg5      // node_potential3,
)
{ 

    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__compute_electric_field set_size %d", set->size);

    opp_profiler->start("ComputeElectricField");

    int nargs = 6;
    opp_arg args[nargs];

    args[0] = std::move(arg0);
    args[1] = std::move(arg1);
    args[2] = std::move(arg2);
    args[3] = std::move(arg3);
    args[4] = std::move(arg4);
    args[5] = std::move(arg5);

    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    if (set_size > 0) 
    {
        computeEF_stride_OPP_HOST_0 = args[0].dat->set->set_capacity;
        computeEF_stride_OPP_HOST_1 = args[1].dat->set->set_capacity;
        computeEF_stride_OPP_HOST_2_MAP = args[2].size;

        cudaMemcpyToSymbol(computeEF_stride_OPP_CUDA_0, &computeEF_stride_OPP_HOST_0, sizeof(int));
        cudaMemcpyToSymbol(computeEF_stride_OPP_CUDA_1, &computeEF_stride_OPP_HOST_1, sizeof(int));
        cudaMemcpyToSymbol(computeEF_stride_OPP_CUDA_2_MAP, &computeEF_stride_OPP_HOST_2_MAP, sizeof(int));

        int start = 0;
        int end   = set->size;

        if (end - start > 0) 
        {
            int nthread = OPP_gpu_threads_per_block;
            int nblocks = (end - start - 1) / nthread + 1;

            opp_cuda_ComputeElectricField <<<nblocks, nthread>>> (
                (double *)  args[0].data_d,
                (double *)  args[1].data_d,       
                (double *)  args[2].data_d,
                (int *)     args[2].map_data_d,
                (double *)  args[3].data_d,
                (double *)  args[4].data_d,
                (double *)  args[5].data_d,
                start, 
                end);
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());

    opp_profiler->end("ComputeElectricField");   
}

//*************************************************************************************************

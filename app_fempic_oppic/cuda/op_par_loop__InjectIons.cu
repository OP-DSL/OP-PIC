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

// AUTO GENERATED CODE

int opDat0_InjectIons_stride_OPPIC_HOST = -1;
int opDat1_InjectIons_stride_OPPIC_HOST = -1;
int opDat2_InjectIons_stride_OPPIC_HOST = -1;
int opDat3_InjectIons_stride_OPPIC_HOST = -1;

__constant__ int opDat0_InjectIons_stride_OPPIC_CONSTANT;
__constant__ int opDat1_InjectIons_stride_OPPIC_CONSTANT;
__constant__ int opDat2_InjectIons_stride_OPPIC_CONSTANT;
__constant__ int opDat3_InjectIons_stride_OPPIC_CONSTANT;


//user function
//*************************************************************************************************
__device__ double rnd_gpu() 
{
    return 0.1; // TODO: Make this a normally distributed random number generator
}

__device__ void inject_ions__kernel_gpu(
    double *pos,
    double *vel,
    double *ef ,
    double *lc, 
    int *cell_index
)
{
    /*sample random position on the inlet face*/
    pos[0 * opDat0_InjectIons_stride_OPPIC_CONSTANT] = -0.1 + 0.2 * rnd_gpu();
    pos[1 * opDat0_InjectIons_stride_OPPIC_CONSTANT] = -0.1 + 0.2 * rnd_gpu();
    pos[2 * opDat0_InjectIons_stride_OPPIC_CONSTANT] = ZERO_double;

    /*injecting cold beam*/
    vel[0 * opDat1_InjectIons_stride_OPPIC_CONSTANT] = ZERO_double;
    vel[1 * opDat1_InjectIons_stride_OPPIC_CONSTANT] = ZERO_double;
    vel[2 * opDat1_InjectIons_stride_OPPIC_CONSTANT] = ION_VELOCITY;

    ef[0 * opDat2_InjectIons_stride_OPPIC_CONSTANT] = ZERO_double;
    ef[1 * opDat2_InjectIons_stride_OPPIC_CONSTANT] = ZERO_double;
    ef[2 * opDat2_InjectIons_stride_OPPIC_CONSTANT] = ZERO_double;

    lc[0 * opDat3_InjectIons_stride_OPPIC_CONSTANT] = ZERO_double;
    lc[1 * opDat3_InjectIons_stride_OPPIC_CONSTANT] = ZERO_double;
    lc[2 * opDat3_InjectIons_stride_OPPIC_CONSTANT] = ZERO_double;
    lc[3 * opDat3_InjectIons_stride_OPPIC_CONSTANT] = ZERO_double;

    *cell_index = 0;
}


// CUDA kernel function
//*************************************************************************************************
__global__ void op_cuda_InjectIons(
    double *__restrict dir_arg0,
    double *__restrict dir_arg1,
    double *__restrict dir_arg2,
    double *__restrict dir_arg3,
    int *__restrict dir_arg4,
    int start,
    int end,
    int set_size
    ) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        //user-supplied kernel call
        inject_ions__kernel_gpu(
            (dir_arg0 + n),
            (dir_arg1 + n),
            (dir_arg2 + n),
            (dir_arg3 + n),
            (dir_arg4 + n)
        );
    }
}


//*************************************************************************************************
void op_par_loop_inject__InjectIons(
    op_set set,      // particles_set
    op_arg arg0,     // part_position,
    op_arg arg1,     // part_velocity,
    op_arg arg2,     // part_electric_field,
    op_arg arg3,     // part_weights,
    op_arg arg4      // part_cell_index,
    )
{ TRACE_ME;

    if (OP_DEBUG) printf("FEMPIC - op_par_loop_inject__InjectIons num_particles %d diff %d\n", set->size, set->diff);

    opDat0_InjectIons_stride_OPPIC_HOST = arg0.dat->set->size;
    opDat1_InjectIons_stride_OPPIC_HOST = arg1.dat->set->size;
    opDat2_InjectIons_stride_OPPIC_HOST = arg2.dat->set->size;
    opDat3_InjectIons_stride_OPPIC_HOST = arg3.dat->set->size;

    cudaMemcpyToSymbol(opDat0_InjectIons_stride_OPPIC_CONSTANT, &opDat0_InjectIons_stride_OPPIC_HOST, sizeof(int));
    cudaMemcpyToSymbol(opDat1_InjectIons_stride_OPPIC_CONSTANT, &opDat1_InjectIons_stride_OPPIC_HOST, sizeof(int));
    cudaMemcpyToSymbol(opDat2_InjectIons_stride_OPPIC_CONSTANT, &opDat2_InjectIons_stride_OPPIC_HOST, sizeof(int));
    cudaMemcpyToSymbol(opDat3_InjectIons_stride_OPPIC_CONSTANT, &opDat3_InjectIons_stride_OPPIC_HOST, sizeof(int));

    int start = (set->size - set->diff);
    int end   = set->size;

    if (end - start > 0) 
    {
        int nthread = GPU_THREADS_PER_BLOCK;
        int nblocks = (end - start - 1) / nthread + 1;

        op_cuda_InjectIons <<<nblocks, nthread>>> (
            (double *)  arg0.data_d,
            (double *)  arg1.data_d,
            (double *)  arg2.data_d,
            (double *)  arg3.data_d,
            (int *)     arg4.data_d,
            start, 
            end, 
            set->size);
    }
}

//*************************************************************************************************
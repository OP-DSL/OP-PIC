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

__constant__ int opDat0_EnrichVelocity_stride_OPPIC_CONSTANT;
__constant__ int opDat1_EnrichVelocity_stride_OPPIC_CONSTANT;

int opDat0_EnrichVelocity_stride_OPPIC_HOST =- 1;
int opDat1_EnrichVelocity_stride_OPPIC_HOST =- 1;


//user function
//*************************************************************************************************
__device__ void enrich_velocity__kernel_gpu(
    double *vel,
    const double *cell_ef,
    const double *dt
)
{
    vel[0 * opDat0_EnrichVelocity_stride_OPPIC_CONSTANT] -= OP_CONST_CUDA_charge / OP_CONST_CUDA_mass * 
                                                        cell_ef[0 * opDat1_EnrichVelocity_stride_OPPIC_CONSTANT] * (0.5 * (*dt));
    
    vel[1 * opDat0_EnrichVelocity_stride_OPPIC_CONSTANT] -= OP_CONST_CUDA_charge / OP_CONST_CUDA_mass * 
                                                        cell_ef[1 * opDat1_EnrichVelocity_stride_OPPIC_CONSTANT] * (0.5 * (*dt));
    
    vel[2 * opDat0_EnrichVelocity_stride_OPPIC_CONSTANT] -= OP_CONST_CUDA_charge / OP_CONST_CUDA_mass * 
                                                        cell_ef[2 * opDat1_EnrichVelocity_stride_OPPIC_CONSTANT] * (0.5 * (*dt));
}


// CUDA kernel function
//*************************************************************************************************
__global__ void op_cuda_EnrichVelocity(
    const int *__restrict d_cell_index,
    double *__restrict dir_arg0,
    const double *__restrict ind_arg1,
    const double *__restrict dir_arg2,
    int start,
    int end,
    int set_size
    ) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;
        
        int map0idx = d_cell_index[n];

        //user-supplied kernel call
        enrich_velocity__kernel_gpu(
            (dir_arg0 + n),
            (ind_arg1 + map0idx),
            dir_arg2
        );
    }
}


//*************************************************************************************************
void op_par_loop_inject__EnrichVelocity(
    op_set set,     // particles_set
    op_arg arg0,    // part_velocity,
    op_arg arg1,    // cell_electric_field,
    op_arg arg2     // const dt,    
    )
{ TRACE_ME;
    
    if (OP_DEBUG) printf("FEMPIC - op_par_loop_inject__EnrichVelocity num_particles %d\n", set->size);

    opDat0_EnrichVelocity_stride_OPPIC_HOST = arg0.dat->set->size;
    opDat1_EnrichVelocity_stride_OPPIC_HOST = arg1.dat->set->size;
    
    cudaMemcpyToSymbol(opDat0_EnrichVelocity_stride_OPPIC_CONSTANT, &opDat0_EnrichVelocity_stride_OPPIC_HOST, sizeof(int));
    cudaMemcpyToSymbol(opDat1_EnrichVelocity_stride_OPPIC_CONSTANT, &opDat1_EnrichVelocity_stride_OPPIC_HOST, sizeof(int));

    int start = (set->size - set->diff);
    int end   = set->size;

    if (end - start > 0) 
    {
        int nthread = GPU_THREADS_PER_BLOCK;
        int nblocks = (end - start - 1) / nthread + 1;

        op_cuda_EnrichVelocity <<<nblocks, nthread>>> (
            (int *)     set->cell_index_dat->data_d,
            (double *)  arg0.data_d,
            (double *)  arg1.data_d,
            (double *)  arg2.data_d,
            start, 
            end, 
            set->size);
    }
}

//*************************************************************************************************
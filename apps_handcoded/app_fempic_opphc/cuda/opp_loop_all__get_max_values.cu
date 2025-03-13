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


//user function
//*************************************************************************************************
__device__ void GetFinalMaxValues__kernel_gpu(
    const OPP_REAL* n_charge_den,
    OPP_REAL* max_n_charge_den,
    const OPP_REAL* n_pot,
    OPP_REAL* max_n_pot
)
{
    *max_n_charge_den = MAX(abs(*n_charge_den), *max_n_charge_den);
    
    *max_n_pot = MAX(*n_pot, *max_n_pot);
}

// CUDA kernel function
//*************************************************************************************************
__global__ void opp_cuda_GetFinalMaxValues(
    const double *__restrict dir_arg0,
    double *__restrict dir_arg1,
    const double *__restrict dir_arg2,
    double *__restrict dir_arg3,
    const int start,
    const int end
) 
{
    double arg1_l[1];
    for ( int d=0; d<1; d++ )
    {
        arg1_l[d] = ZERO_double;
    }
    double arg3_l[1];
    for ( int d=0; d<1; d++ )
    {
        arg3_l[d] = ZERO_double;
    }

    //process set elements
    for (int n = threadIdx.x+blockIdx.x*blockDim.x; n < (end-start); n+= blockDim.x*gridDim.x)
    {
        //user-supplied kernel call
        GetFinalMaxValues__kernel_gpu(
            (dir_arg0 + n),
            arg1_l,
            (dir_arg2 + n),
            arg3_l
        );
    }

    for (int d = 0; d < 1; d++)
    {
        opp_reduction<OPP_MAX>(&dir_arg1[d + blockIdx.x * 1], arg1_l[d]);
    }
    for (int d = 0; d < 1; d++)
    {
        opp_reduction<OPP_MAX>(&dir_arg3[d + blockIdx.x * 1], arg3_l[d]);
    }
}

//*************************************************************************************************
void opp_loop_all__get_max_values(
    opp_set set,     // cells set
    opp_arg arg0,    // n_charge_den            // OPP_READ
    opp_arg arg1,    // global_max_n_charge_den // OPP_MAX
    opp_arg arg2,    // n_pot                   // OPP_READ
    opp_arg arg3     // global_max_n_pot        // OPP_MAX
)
{ 
    
    if (FP_DEBUG) opp_printf("FEMPIC", "opp_loop_all__get_max_values set_size %d", set->size);

    opp_profiler->start("GetMax");

    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    double* arg1h = (double *)args[1].data;
    double* arg3h = (double *)args[3].data;

    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    if (set_size > 0) 
    {

        //set CUDA execution parameters
        #ifdef OPP_BLOCK_SIZE_4
        int nthread = OPP_BLOCK_SIZE_4;
        #else
        int nthread = 32; // OPP_block_size; // TODO : CHECK this
        #endif

        int nblocks = 200; // why? TODO : Check

        //transfer global reduction data to GPU
        int maxblocks = nblocks;
        int reduct_bytes = 0;
        int reduct_size  = 0;
        
        reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double)); // for global arg 1
        reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double)); // for global arg 3
        reduct_size   = MAX(reduct_size,sizeof(double));
        
        opp_reallocReductArrays(reduct_bytes);
        
        reduct_bytes = 0;
        args[1].data   = OPP_reduct_h + reduct_bytes;
        args[1].data_d = OPP_reduct_d + reduct_bytes;
        for (int b = 0; b < maxblocks; b++)
        {
            for (int d = 0; d < 1; d++)
            {
                ((double *)args[1].data)[d+b*1] = ZERO_double;
            }
        }
        reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));

        args[3].data   = OPP_reduct_h + reduct_bytes;
        args[3].data_d = OPP_reduct_d + reduct_bytes;
        for (int b = 0; b < maxblocks; b++)
        {
            for (int d = 0; d < 1; d++)
            {
                ((double *)args[3].data)[d+b*1] = ZERO_double;
            }
        }
        reduct_bytes += ROUND_UP(maxblocks*1*sizeof(double));

        opp_mvReductArraysToDevice(reduct_bytes);

        int nshared = reduct_size*nthread;

        const int start = 0;
        const int end   = set->size;

        opp_cuda_GetFinalMaxValues <<<nblocks,nthread,nshared>>> (
            (double *)  args[0].data_d,
            (double *)  args[1].data_d,
            (double *)  args[2].data_d,
            (double *)  args[3].data_d,
            start, 
            end
        );

        opp_mvReductArraysToHost(reduct_bytes);

        for (int b = 0; b < maxblocks; b++)
        {
            for (int d = 0; d < 1; d++)
            {
                arg1h[d] = MAX(arg1h[d], ((double *)args[1].data)[d+b*1]);
            }
        }
        args[1].data = (char *)arg1h;
        opp_mpi_reduce(&args[1],arg1h);
        
        for (int b = 0; b < maxblocks; b++)
        {
            for (int d = 0; d < 1; d++)
            {
                arg3h[d] = MAX(arg3h[d], ((double *)args[3].data)[d+b*1]);
            }
        }
        args[3].data = (char *)arg3h;

        opp_mpi_reduce(&args[3],arg3h);

        opp_set_dirtybit_grouped(nargs, args, Device_GPU);
        cutilSafeCall(cudaDeviceSynchronize());       
    }

    opp_profiler->end("GetMax");
}
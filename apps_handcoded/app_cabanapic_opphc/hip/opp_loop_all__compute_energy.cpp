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

int ce_OPP_HOST_1 = -1;
__constant__ int ce_OPP_DEVICE_1;

//user function
//*************************************************************************************************
__device__ void compute_energy__kernel_gpu(
    const OPP_INT* cell0_ghost, 
    const OPP_REAL* cell_field, 
    OPP_REAL* energy)
{
    if (cell0_ghost[0] == 0) 
    {
        // energy ptr is thread local, hence no need of atomics here
        
        energy[0] += cell_field[Dim::x * ce_OPP_DEVICE_1] * cell_field[Dim::x * ce_OPP_DEVICE_1] +
                    cell_field[Dim::y * ce_OPP_DEVICE_1] * cell_field[Dim::y * ce_OPP_DEVICE_1] +
                    cell_field[Dim::z * ce_OPP_DEVICE_1] * cell_field[Dim::z * ce_OPP_DEVICE_1];
    }
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void opp_dev_compute_energy(
    const OPP_INT *__restrict dir_arg0,
    const OPP_REAL *__restrict dir_arg1,
    OPP_REAL *__restrict dir_arg2,
    const int start,
    const int end
) 
{
    double arg2_l[1];
    for ( int d=0; d<1; d++ )
    {
        arg2_l[d] = ZERO_double;
    }

    //process set elements
    for (int n = threadIdx.x+blockIdx.x*blockDim.x; n < (end-start); n+= blockDim.x*gridDim.x)
    {
        //user-supplied kernel call
        compute_energy__kernel_gpu(
            (dir_arg0 + n),
            (dir_arg1 + n),
            arg2_l
        );
    }

    for (int d = 0; d < 1; d++)
    {
        opp_reduction<OPP_INC>(&dir_arg2[d + blockIdx.x * 1], arg2_l[d]);
    }
}

//*************************************************************************************************
void opp_loop_all__compute_energy(
    opp_set set,     // cells set
    opp_arg arg0,    // cell0_ghost, OPP_READ
    opp_arg arg1,    // cell_field,  OPP_READ
    opp_arg arg2     // energy,      OPP_INC
)
{
    if (OPP_DBG) opp_printf("CABANA", "opp_loop_all__compute_energy set_size %d", set->size);

    opp_profiler->start("Energy");

    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    OPP_REAL* arg2h = (OPP_REAL *)args[2].data;

    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    if (set_size > 0) 
    {
        ce_OPP_HOST_1 = args[1].dat->set->set_capacity; 
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(ce_OPP_DEVICE_1), &ce_OPP_HOST_1, sizeof(int)));

        //set DEVICE execution parameters
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
        
        reduct_bytes += ROUND_UP(maxblocks*1*sizeof(OPP_REAL)); // for global arg 2
        reduct_size   = MAX(reduct_size,sizeof(OPP_REAL));
        
        opp_reallocReductArrays(reduct_bytes);
        
        reduct_bytes = 0;
        args[2].data   = OPP_reduct_h + reduct_bytes;
        args[2].data_d = OPP_reduct_d + reduct_bytes;
        for (int b = 0; b < maxblocks; b++)
        {
            for (int d = 0; d < 1; d++)
            {
                ((OPP_REAL *)args[2].data)[d+b*1] = ZERO_double;
            }
        }
        reduct_bytes += ROUND_UP(maxblocks*1*sizeof(OPP_REAL));

        opp_mvReductArraysToDevice(reduct_bytes);

        int nshared = reduct_size*nthread;

        const int start = 0;
        const int end   = set->size;

        opp_dev_compute_energy<<<nblocks,nthread,nshared>>> (
            (OPP_INT *)  args[0].data_d,
            (OPP_REAL *) args[1].data_d,
            (OPP_REAL *) args[2].data_d,
            start, 
            end
        );

        opp_mvReductArraysToHost(reduct_bytes);

        for (int b = 0; b < maxblocks; b++)
        {
            for (int d = 0; d < 1; d++)
            {
                arg2h[d] += ((OPP_REAL *)args[2].data)[d+b*1];
            }
        }
        args[2].data = (char *)arg2h;
        opp_mpi_reduce(&args[2],arg2h);

        opp_set_dirtybit_grouped(nargs, args, Device_GPU);
        cutilSafeCall(hipDeviceSynchronize());       
    }

    opp_profiler->end("Energy");
}
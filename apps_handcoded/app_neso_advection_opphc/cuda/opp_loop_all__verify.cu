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

int verify_stride_OPP_HOST_1 = -1;
__constant__ int verify_stride_OPP_DEVICE_1;

//user function
//*************************************************************************************************
__device__ void verify__kernel_gpu(
        const OPP_INT* part_cid, 
        const OPP_REAL* part_pos,
        const OPP_INT* cell_global_idx,
        OPP_INT* incorrect_part_count
)
{
    // get the cell boundaries for the current cell_index - using global cell index 
    int ix = -1, iy = -1;
    RANK_TO_INDEX((*cell_global_idx), ix, iy, CONST_DEVICE_ndimcells[Dim::x]); 
    
    if (ix < 0 || iy < 0)
    {
        // opp_printf("VERIFY", "Incorrect ix[%d] iy[%d] for global cell[%d] nx[%d]",
        //     ix, iy, (*cell_global_idx), CONST_ndimcells[Dim::x]);
        (*incorrect_part_count)++;
        return;
    }
    
    // get the boundaries of that cell
    const OPP_REAL boundary_ll[DIM] = { (ix * CONST_DEVICE_cell_width), (iy * CONST_DEVICE_cell_width) };
 
    // check whether the current particle is within those boundaries or not!
    const OPP_REAL part_pos_x = part_pos[Dim::x * verify_stride_OPP_DEVICE_1];
    if (part_pos_x < boundary_ll[Dim::x] ||
            part_pos_x > (boundary_ll[Dim::x] + CONST_DEVICE_cell_width)) {
        
        (*incorrect_part_count)++;
        return;
    }

    const OPP_REAL part_pos_y = part_pos[Dim::y * verify_stride_OPP_DEVICE_1];
    if (part_pos_y < boundary_ll[Dim::y] ||
            part_pos_y > (boundary_ll[Dim::y] + CONST_DEVICE_cell_width)) {
        
        (*incorrect_part_count)++;
        return;
    }
}

// HIP kernel function
//*************************************************************************************************
__global__ void opp_device_Verify(
    const OPP_INT *__restrict dir_arg0,
    const OPP_REAL *__restrict dir_arg1,
    const OPP_INT *__restrict ind_arg2,
    OPP_INT *__restrict dir_arg3,
    const int start,
    const int end
) 
{
    OPP_INT arg3_l[1];
    for ( int d=0; d<1; d++ )
    {
        arg3_l[d] = 0;
    }

    //process set elements
    for (int n = threadIdx.x+blockIdx.x*blockDim.x; n < (end-start); n+= blockDim.x*gridDim.x)
    {
        const int map0idx = dir_arg0[n];

        //user-supplied kernel call
        verify__kernel_gpu(
            (dir_arg0 + n),             // part_mesh_rel,      
            (dir_arg1 + n),             // part_pos,           
            (ind_arg2 + map0idx),       // cell_global_index,  
            arg3_l                      // incorrect_part_count
        );
    }

    for (int d = 0; d < 1; d++)
    {
        opp_reduction<OPP_INC>(&dir_arg3[d + blockIdx.x * 1], arg3_l[d]);
    }
}

//*************************************************************************************************
void opp_loop_all__verify(
    opp_set set,     // particles_set
    opp_arg arg0,    // part_mesh_rel,        OPP_RW
    opp_arg arg1,    // part_pos,             OPP_READ
    opp_arg arg2,    // cell_global_index,    OPP_READ
    opp_arg arg3     // incorrect_part_count, OPP_INC
)
{ 
    
    if (OPP_DBG) 
        opp_printf("ADVEC", "opp_loop_all__verify set_size %d diff %d", set->size, set->diff);

    opp_profiler->start("Verify");

    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    OPP_INT* arg3h = (OPP_INT *)args[3].data;

    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    if (set_size > 0) 
    {
        verify_stride_OPP_HOST_1 = args[1].dat->set->set_capacity;
        cutilSafeCall(cudaMemcpyToSymbol(verify_stride_OPP_DEVICE_1, 
                                                    &verify_stride_OPP_HOST_1, sizeof(int)));

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
        
        reduct_bytes += ROUND_UP(maxblocks*1*sizeof(OPP_INT)); // for global arg 3
        reduct_size   = MAX(reduct_size,sizeof(OPP_INT));
        
        opp_reallocReductArrays(reduct_bytes);
        
        reduct_bytes = 0;
        args[3].data   = OPP_reduct_h + reduct_bytes;
        args[3].data_d = OPP_reduct_d + reduct_bytes;
        for (int b = 0; b < maxblocks; b++)
        {
            for (int d = 0; d < 1; d++)
            {
                ((OPP_INT *)args[3].data)[d+b*1] = 0;
            }
        }
        reduct_bytes += ROUND_UP(maxblocks*1*sizeof(OPP_INT));

        opp_mvReductArraysToDevice(reduct_bytes);

        int nshared = reduct_size*nthread;

        const int start = 0;
        const int end   = set->size;

        opp_device_Verify <<<nblocks,nthread,nshared>>> (
            (OPP_INT *) args[0].data_d,                 // part_mesh_rel,      
            (OPP_REAL*) args[1].data_d,                 // part_pos,           
            (OPP_INT *) args[2].data_d,                 // cell_global_index,  
            (OPP_INT *) args[3].data_d,                 // incorrect_part_count
            start, 
            end
        );

        opp_mvReductArraysToHost(reduct_bytes);

        for (int b = 0; b < maxblocks; b++)
        {
            for (int d = 0; d < 1; d++)
            {
                arg3h[d] += ((OPP_INT *)args[3].data)[d+b*1];
            }
        }
        args[3].data = (char *)arg3h;

        opp_mpi_reduce(&args[3],arg3h);

        opp_set_dirtybit_grouped(nargs, args, Device_GPU);
        cutilSafeCall(cudaDeviceSynchronize());       
    }

    opp_profiler->end("Verify");
}
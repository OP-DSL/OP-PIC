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


#include "oppic_omp.h"
#include "../advec_defs.h"

OPP_REAL CONST_extents[2];
OPP_REAL CONST_dt = 0.0;
OPP_REAL CONST_cell_width = 0.0;
OPP_INT CONST_ndimcells[2];

//****************************************
void opp_decl_const_impl(int dim, int size, char* data, const char* name)
{
    if (!strcmp(name,"CONST_extents"))         std::memcpy(CONST_extents, data, (size*dim));
    else if (!strcmp(name,"CONST_dt"))         std::memcpy(&CONST_dt, data, (size*dim));
    else if (!strcmp(name,"CONST_cell_width")) std::memcpy(&CONST_cell_width, data, (size*dim));
    else if (!strcmp(name,"CONST_ndimcells"))  std::memcpy(CONST_ndimcells, data, (size*dim));
    else std::cerr << "error: unknown const name" << std::endl;
}
//****************************************

#include "../kernels.h"

//*************************************************************************************************
void opp_particle_mover__UpdatePosMove(
    opp_set set,        // particles_set
    opp_arg arg0,       // part_mesh_rel, OP_RW
    opp_arg arg1,       // part_vel,      OP_RW
    opp_arg arg2,       // part_pos,      OP_RW
    opp_arg arg3,       // cell_pos_ll,   OP_READ
    opp_arg arg4        // cell_cell_map, OP_READ
)
{

    if (OP_DEBUG) 
        opp_printf("ADVEC", "opp_particle_mover__UpdatePosMove set_size %d diff %d", set->size, set->diff);

    opp_profiler->start("Move");

    int nargs = 5;
    opp_arg args[nargs];

    args[0]  = std::move(arg0);
    args[1]  = std::move(arg1);
    args[2]  = std::move(arg2);
    args[3]  = std::move(arg3);
    args[4]  = std::move(arg4);

    const int nthreads = omp_get_max_threads();
 
    int set_size = opp_mpi_halo_exchanges(set, nargs, args); 
    
    opp_mpi_halo_wait_all(nargs, args);  // unable to overlap computation and communication

    if (set_size > 0) 
    {
        do {

            opp_init_particle_move(set, nargs, args);

            if (OP_DEBUG) 
                opp_printf("ADVEC", "opp_particle_mover__UpdatePosMove Starting iteration %d, start[%d] end[%d] nthreads %d", 
                    OPP_comm_iteration, OPP_iter_start, OPP_iter_end, nthreads);

            #pragma omp parallel for
            for (int thr = 0; thr < nthreads; thr++)
            {
                size_t start  = (size_t)OPP_iter_start + ((size_t)(OPP_iter_end - OPP_iter_start) * thr) / nthreads;
                size_t finish = (size_t)OPP_iter_start + ((size_t)(OPP_iter_end - OPP_iter_start) * (thr+1)) / nthreads;

                int *cellIdx = nullptr;

                for (size_t n = start; n < finish; n++)
                {          
                    opp_move_var m = opp_get_move_var(thr);

                    do
                    { 
                        cellIdx = &(OPP_mesh_relation_data[n]);

                        push_particles_kernel(m, 
                            &((OPP_INT*)  args[0].data)[n * args[0].dim],        // part_cid 
                            &((OPP_REAL*) args[1].data)[n * args[1].dim],        // part_vel 
                            &((OPP_REAL*) args[2].data)[n * args[2].dim],        // part_pos 
                            &((OPP_REAL*) args[3].data)[*cellIdx * args[3].dim], // cell_pos_ll 
                            &((OPP_INT*)  args[4].data)[*cellIdx * args[4].dim]  // cell_cell_map 
                        );

                    } while (opp_part_check_status_omp(m, *cellIdx, set, n, thr, thr));  
                }
            }

        } while(opp_finalize_particle_move(set));
    } 
    

    opp_profiler->end("Move");
}

//*******************************************************************************
void opp_loop_all__Verify(
    opp_set set,        // particles_set
    opp_arg arg0,       // part_mesh_rel,        OP_RW
    opp_arg arg1,       // part_pos,             OP_READ
    opp_arg arg2,       // cell_global_index,    OP_READ
    opp_arg arg3        // incorrect_part_count, OP_INC
)
{
    if (OP_DEBUG) 
        opp_printf("ADVEC", "opp_loop_all__Verify set_size %d diff %d", set->size, set->diff);

    opp_profiler->start("Verify");

    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    const int set_size = set->size;
    const int nthreads = omp_get_max_threads();

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 

    OPP_REAL arg3_l[nthreads*1];
    for (int thr = 0; thr < nthreads; thr++)
    {
        for (int d = 0; d < 1; d++) // can have multiple dimension defined for global_arg
            arg3_l[1 * thr + d] = ZERO_double;
    }  

    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        size_t start  = ((size_t)set_size * thr) / nthreads;
        size_t finish = ((size_t)set_size * (thr+1)) / nthreads;

        for (size_t n = start; n < finish; n++)
        { 
            const int map0idx = OPP_mesh_relation_data[n];

            verify_kernel( 
                &((OPP_INT*)  arg0.data)[n * arg0.dim],       // part_mesh_rel,      
                &((OPP_REAL*) arg1.data)[n * arg1.dim],       // part_pos,           
                &((OPP_INT*)  arg2.data)[map0idx * arg2.dim], // cell_global_index,  
                (int*) arg3.data                              // incorrect_part_count
            );
        }
    }

    for (int thr = 0; thr < nthreads; thr++) 
    {
        for (int d = 0; d < 1; d++) // can have multiple dimension defined for global_arg
            ((double*)args[3].data)[d] += arg3_l[1 * thr + d];
    }

    opp_profiler->end("Verify");
}
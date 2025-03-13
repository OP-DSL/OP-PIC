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
// USER WRITTEN CODE
//*********************************************

#include "opp_templates.h"

OPP_REAL CONST_lhs_voltage[1];
OPP_REAL CONST_L[1];
OPP_REAL CONST_xl[1];
OPP_REAL CONST_xr[1];
OPP_REAL CONST_dx[1];
OPP_REAL CONST_qm[1];
OPP_REAL CONST_dt[1];
OPP_REAL CONST_qscale[1];
OPP_INT CONST_neighbour_cell[2];
OPP_INT CONST_rank[1];
OPP_INT CONST_comm_size[1];

#include "simpic_defs.h"
#include "kernels.h"

//*********************************************MAIN****************************************************
int main(int argc, char* argv[]) 
{
    opp_init(argc, argv);

    {    
    opp_profiler->start("Setup");

        OPP_INT neighbour_cell[2] = { MAX_CELL_INDEX, MAX_CELL_INDEX };

        std::unique_ptr<DataPointers> m = Load();

        opp_set nodes_set = opp_decl_set(m->n_nodes, "mesh_nodes");
        opp_set cells_set = opp_decl_set(m->n_cells, "mesh_cells");
        opp_set parts_set = opp_decl_particle_set(m->n_particles, "particles", cells_set);  

        opp_map c2c_map   = opp_decl_map(cells_set, cells_set, 2, m->cell_to_cells_tmp, "c2c_map");
        opp_map c2n_map   = opp_decl_map(cells_set, nodes_set, 2, m->cell_to_nodes_tmp, "c2n_map");
        
        opp_dat n_field_e = opp_decl_dat(nodes_set, 1, DT_REAL, m->node_field_E_tmp, "n_field_e");
        opp_dat n_field_j = opp_decl_dat(nodes_set, 1, DT_REAL, m->node_field_J_tmp, "n_field_j");
        opp_dat n_field_p = opp_decl_dat(nodes_set, 1, DT_REAL, m->node_field_P_tmp, "n_field_p");
        opp_dat n_xlocal  = opp_decl_dat(nodes_set, 1, DT_REAL, m->node_xlocal_tmp,  "n_xlocal"); 

        opp_map p2c_map   = opp_decl_map(parts_set, cells_set, 1, m->part_cell_index_tmp, "p2c_map");

        opp_dat p_pos_x   = opp_decl_dat(parts_set, 1, DT_REAL, m->part_position_x_tmp, "p_pos_x");      
        opp_dat p_vel_x   = opp_decl_dat(parts_set, 1, DT_REAL, m->part_velocity_x_tmp, "p_vel_x");      
        opp_dat p_field_e = opp_decl_dat(parts_set, 1, DT_REAL, m->part_field_E_tmp,    "p_field_e");         

        opp_decl_const<OPP_REAL>(1, &lhsvoltage,   "CONST_lhs_voltage");
        opp_decl_const<OPP_REAL>(1, &L,            "CONST_L");
        opp_decl_const<OPP_REAL>(1, &xl,           "CONST_xl");
        opp_decl_const<OPP_REAL>(1, &xr,           "CONST_xr");
        opp_decl_const<OPP_REAL>(1, &dx,           "CONST_dx"); 
        opp_decl_const<OPP_REAL>(1, &qm,           "CONST_qm");
        opp_decl_const<OPP_REAL>(1, &dt,           "CONST_dt");
        opp_decl_const<OPP_REAL>(1, &qscale,       "CONST_qscale");
        opp_decl_const<OPP_INT>(2, neighbour_cell, "CONST_neighbour_cell");
        opp_decl_const<OPP_INT>(1, &OPP_rank,      "CONST_rank");
        opp_decl_const<OPP_INT>(1, &OPP_comm_size, "CONST_comm_size");

        m->DeleteValues();

        opp_printf("Setup", "Cells[%d] Nodes[%d] Particles[%d] dt[%d]", cells_set->size, 
            nodes_set->size, parts_set->size, dt);

    opp_profiler->end("Setup");

    opp_profiler->start("MainLoop");
        for (OPP_main_loop_iter = 0; OPP_main_loop_iter < ntimesteps; ++OPP_main_loop_iter)
        {
            // STEP 1 - Weight mesh field values to particle positions ************************************
            opp_par_loop(weight_f2p_kernel, "weight_f2p_kernel", parts_set, OPP_ITERATE_ALL,
                opp_arg_dat(n_field_e, 0, c2n_map, p2c_map, OPP_READ),  // lhs n_field_e
                opp_arg_dat(n_field_e, 1, c2n_map, p2c_map, OPP_READ),  // rhs n_field_e  
                opp_arg_dat(p_pos_x,                        OPP_READ),    
                opp_arg_dat(p_field_e,                      OPP_WRITE)     
            );

            // STEP 2 - Move the particles with the influence of the fields *******************************
            opp_particle_move(move_kernel, "move", parts_set, c2c_map, p2c_map,
                opp_arg_dat(p_field_e, OPP_READ),
                opp_arg_dat(p_vel_x,   OPP_RW),
                opp_arg_dat(p_pos_x,   OPP_RW)
            );

            // STEP 3 - Gather the contribution from particle movement to the mesh ************************
            opp_reset_dat(n_field_j, (char*)opp_zero_double16);

            opp_par_loop(weight_p2f_kernel, "weight_p2f_kernel", parts_set, OPP_ITERATE_ALL,
                opp_arg_dat(n_field_j, 0, c2n_map, p2c_map, OPP_INC),  // lhs n_field_j
                opp_arg_dat(n_field_j, 1, c2n_map, p2c_map, OPP_INC),  // rhs n_field_j
                opp_arg_dat(p_pos_x,                        OPP_READ)  
            );
        
            // STEP 4 - Solve field values on the mesh points ************************************
            seq_field_solve_poissons_equation(nodes_set, // TODO : This needs to be converted to kernels or to Petsc
                opp_arg_dat(n_field_j, OPP_RW),
                opp_arg_dat(n_field_p, OPP_RW)
            );
            opp_par_loop(sum_laplace_kernel, "sum_laplace_kernel", nodes_set, OPP_ITERATE_ALL,
                opp_arg_dat(n_xlocal,  OPP_READ),
                opp_arg_dat(n_field_p, OPP_RW)
            );
            seq_field_solve_get_potential_gradient(nodes_set, // TODO : This needs to be converted to kernels or to Petsc
                opp_arg_dat(n_field_e, OPP_WRITE),
                opp_arg_dat(n_field_p, OPP_READ)
            ); 

            OPP_RUN_ON_ROOT()
                opp_printf("Main", "loop %d Done | particles: %d ****", OPP_main_loop_iter, parts_set->size);
        }
    opp_profiler->end("MainLoop");

        OPP_RUN_ON_ROOT()
            opp_printf("Main", "Main loop completed after %d iterations ****", OPP_main_loop_iter);
    }

    opp_exit();

    return 0;
}

//*************************************************************************************************
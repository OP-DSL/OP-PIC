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


#include "simpic.h"

//*********************************************MAIN****************************************************
int main(int argc, char* argv[]) 
{
    opp_init(argc, argv);
    opp_params->write(std::cout);
    
    std::unique_ptr<Loader> loader = Load();

    int loop_count = 0;
    double max_time = t;

    {    
        opp_set nodes_set         = opp_decl_set(loader->n_nodes, "mesh_nodes");
        opp_set cells_set         = opp_decl_set(loader->n_cells, "mesh_cells");
        opp_set particles_set     = opp_decl_particle_set(loader->n_particles, "particles", cells_set);  
    
        opp_map cell_to_nodes_map = opp_decl_mesh_map(cells_set, nodes_set, 2, loader->cell_to_nodes_tmp, "cell_to_edges_map");
        
        opp_dat node_field_E      = opp_decl_dat(nodes_set, 1, DT_REAL, loader->node_field_E_tmp, "node_field_E");
        opp_dat node_field_J      = opp_decl_dat(nodes_set, 1, DT_REAL, loader->node_field_J_tmp, "node_field_J");
        opp_dat node_field_P      = opp_decl_dat(nodes_set, 1, DT_REAL, loader->node_field_P_tmp, "node_field_P");
        opp_dat node_xlocal       = opp_decl_dat(nodes_set, 1, DT_REAL, loader->node_xlocal_tmp,  "node_xlocal"); 

        opp_dat part_position_x   = opp_decl_dat(particles_set, 1, DT_REAL, loader->part_position_x_tmp, "part_position_x");      
        opp_dat part_velocity_x   = opp_decl_dat(particles_set, 1, DT_REAL, loader->part_velocity_x_tmp, "part_velocity_x");      
        opp_dat part_field_E      = opp_decl_dat(particles_set, 1, DT_REAL, loader->part_field_E_tmp,    "part_field_E");         
        opp_dat part_cell_index   = opp_decl_dat(particles_set, 1, DT_INT,  loader->part_cell_index_tmp, "part_cell_index", true);

        // Set constants here, using opp_decl_const<>() if constant values are used in kernels

        loader->DeleteValues();

        auto start = std::chrono::system_clock::now();

        for (double tt = 0; tt < max_time; tt += dt, ++loop_count)
        {
            // STEP 1 - Weight mesh field values to particle positions ************************************
            opp_par_loop_all__WeightFieldsToParticles(
                particles_set,
                opp_arg_dat(node_field_E,    0, cell_to_nodes_map, OPP_READ, OPP_Map_from_Mesh_Rel),  // lhs node_field_E
                opp_arg_dat(node_field_E,    1, cell_to_nodes_map, OPP_READ, OPP_Map_from_Mesh_Rel),  // rhs node_field_E  
                opp_arg_dat(part_position_x,                       OPP_READ),    
                opp_arg_dat(part_field_E,                          OPP_WRITE)     
            );

            // STEP 2 - Move the particles with the influence of the fields ************************************
            opp_par_loop_particle_all__PushParticles(
                particles_set, 
                opp_arg_dat(part_field_E,    OPP_READ),
                opp_arg_dat(part_velocity_x, OPP_RW),
                opp_arg_dat(part_position_x, OPP_RW),
                opp_arg_dat(part_cell_index, OPP_WRITE)
            );

            // STEP 3 - Gather the contribution from particle movement to the mesh ************************************
            opp_reset_dat(node_field_J, (char*)opp_zero_double16);
            opp_par_loop_all__WeightParticlesToFields(
                particles_set,
                opp_arg_dat(node_field_J, 0, cell_to_nodes_map, OPP_INC, OPP_Map_from_Mesh_Rel),  // lhs node_field_J
                opp_arg_dat(node_field_J, 1, cell_to_nodes_map, OPP_INC, OPP_Map_from_Mesh_Rel),  // rhs node_field_J
                opp_arg_dat(part_position_x,                    OPP_READ)  
            );
        
            // STEP 4 - Solve field values on the mesh points ************************************
            seq_field_solve_poissons_equation( // TODO : This needs to be converted to kernels or to Petsc
                nodes_set->size,
                (OPP_REAL*)node_field_J->data,
                (OPP_REAL*)node_field_P->data
            );
            opp_par_loop_all__FieldSolveSumLaplace(
                nodes_set,
                opp_arg_dat(node_xlocal,  OPP_READ),
                opp_arg_dat(node_field_P, OPP_WRITE)
            );
            seq_field_solve_get_potential_gradient( // TODO : This needs to be converted to kernels or to Petsc
                nodes_set->size,
                (OPP_REAL*)node_field_E->data,
                (OPP_REAL*)node_field_P->data
            ); 

            printf("Looping dt = %+2.15lE, tmax = %+2.15lE loop %d ******\n", tt, max_time, loop_count); 
        }

        // opp_print_dat_to_txtfile(part_position_x, "F", "part_position_x.dat");
        // opp_print_map_to_txtfile(cell_to_nodes_map, "F", "cell_to_nodes_map.dat");

        std::chrono::duration<double> diff = std::chrono::system_clock::now() - start;
        std::cout << "\nSIMPIC - Time to iterate " << (loop_count-1) << " <sec>: " << diff.count() << "\n\n";

        opp_exit();
    } // End Scope for oppic


    return 0;
}

//*************************************************************************************************
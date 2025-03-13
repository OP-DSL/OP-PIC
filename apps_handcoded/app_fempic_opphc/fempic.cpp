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

#include "opp_seq.h"
#include "fempic.h"
#include "kernels.h"

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << std::endl;
        exit(-1);
    }

    opp::Params params(argv[1]);
    params.write(std::cout);

    FieldPointers mesh = LoadMesh(params, argc, argv);

    double plasma_den = params.get<REAL>("plasma_den");
    double dt = params.get<REAL>("dt");
    double ion_velocity = params.get<REAL>("ion_velocity");
    double spwt = params.get<REAL>("spwt");
    double mass = 2 * AMU;
    double charge = 1 * QE;
    int max_iter = params.get<INT>("max_iter");
    double remainder = 0.0;
    int ts = 0;

    { // Start Scope for oppic
        opp_init(argc, argv, &params);
 
        opp_set nodes_set            = opp_decl_set(mesh.n_nodes, "mesh_nodes");
        opp_set cells_set            = opp_decl_set(mesh.n_cells, "mesh_cells");
        opp_set ifaces_set           = opp_decl_set(mesh.n_ifaces, "inlet_faces_cells");
        opp_set particles_set        = opp_decl_particle_set("particles", cells_set); 

        opp_map cell_to_nodes_map    = opp_decl_map(cells_set,  nodes_set, NODES_PER_CELL,  mesh.cell_to_nodes,  "cell_to_nodes_map");
        opp_map cell_to_cell_map     = opp_decl_map(cells_set,  cells_set, NEIGHBOUR_CELLS, mesh.cell_to_cell,   "cell_to_cell_map"); 
        opp_map iface_to_cell_map    = opp_decl_map(ifaces_set, cells_set, 1,               mesh.iface_to_cell,  "iface_to_cell_map"); 

        opp_dat cell_determinants    = opp_decl_dat(cells_set, (NEIGHBOUR_CELLS*DET_FIELDS), OPP_REAL, (char*)mesh.cell_det,         "cell_determinants");  
        opp_dat cell_volume          = opp_decl_dat(cells_set, 1,                            OPP_REAL, (char*)mesh.cell_volume,      "cell_volume");        
        opp_dat cell_electric_field  = opp_decl_dat(cells_set, DIMENSIONS,                   OPP_REAL, (char*)mesh.cell_ef,          "cell_electric_field");
        opp_dat cell_shape_deriv     = opp_decl_dat(cells_set, (NODES_PER_CELL*DIMENSIONS),  OPP_REAL, (char*)mesh.cell_shape_deriv, "cell_shape_deriv"); 
    
        opp_dat node_volume          = opp_decl_dat(nodes_set, 1,          OPP_REAL, (char*)mesh.node_volume,  "node_volume");        
        opp_dat node_potential       = opp_decl_dat(nodes_set, 1,          OPP_REAL, (char*)mesh.node_pot,     "node_potential");     
        opp_dat node_charge_density  = opp_decl_dat(nodes_set, 1,          OPP_REAL, (char*)mesh.node_ion_den, "node_charge_density");

        opp_dat iface_v_normal       = opp_decl_dat(ifaces_set, DIMENSIONS, OPP_REAL, (char*)mesh.iface_v_normal,      "iface_v_normal");        
        opp_dat iface_u_normal       = opp_decl_dat(ifaces_set, DIMENSIONS, OPP_REAL, (char*)mesh.iface_u_normal,      "iface_u_normal"); 
        opp_dat iface_normal         = opp_decl_dat(ifaces_set, DIMENSIONS, OPP_REAL, (char*)mesh.iface_normal,        "iface_normal");     
        opp_dat iface_area           = opp_decl_dat(ifaces_set, 1,          OPP_REAL, (char*)mesh.iface_area,          "iface_area");
        opp_dat iface_inj_part_dist  = opp_decl_dat(ifaces_set, 1,          OPP_INT,  (char*)mesh.iface_inj_part_dist, "iface_inj_part_dist");
        opp_dat iface_node_pos       = opp_decl_dat(ifaces_set, DIMENSIONS, OPP_REAL, (char*)mesh.iface_node_pos,      "iface_node_pos"); 

        opp_dat part_position        = opp_decl_particle_dat(particles_set, DIMENSIONS,     OPP_REAL, nullptr, "part_position");
        opp_dat part_velocity        = opp_decl_particle_dat(particles_set, DIMENSIONS,     OPP_REAL, nullptr, "part_velocity");    
        opp_dat part_lc              = opp_decl_particle_dat(particles_set, NODES_PER_CELL, OPP_REAL, nullptr, "part_lc");
        opp_dat part_mesh_relation   = opp_decl_particle_dat(particles_set, 1,              OPP_INT,  nullptr, "part_mesh_relation", true); // new cell index field

        opp_set dummy_part_set       = opp_decl_particle_set(mesh.n_approx_injected, "dummy particles", cells_set); 
        opp_dat dum_part_random      = opp_decl_dat(dummy_part_set, 2, OPP_REAL, (char*)mesh.dummy_part_random, "dum_part_random");

        opp_decl_const<double>(1, &spwt,         "CONST_spwt");
        opp_decl_const<double>(1, &ion_velocity, "CONST_ion_velocity");
        opp_decl_const<double>(1, &dt,           "CONST_dt");
        opp_decl_const<double>(1, &plasma_den,   "CONST_plasma_den");
        opp_decl_const<double>(1, &mass,         "CONST_mass");
        opp_decl_const<double>(1, &charge,       "CONST_charge");

        mesh.DeleteValues();

        auto start = std::chrono::system_clock::now();
        auto start_iter1 = std::chrono::system_clock::now();

        for (ts = 0; ts < max_iter; ts++)
        {
            if (ts == 1) start_iter1 = std::chrono::system_clock::now();

            int injected_count = 0;
            opp_inject__Increase_particle_count(
                particles_set,                                                                // inject to particles_set
                ifaces_set,                                                                   // inlect_face_set
                opp_arg_gbl(&(injected_count), 1, "int",   OPP_RW),                          // injected total global,
                opp_arg_dat(iface_area,                    OPP_READ),                        // iface_area,
                opp_arg_dat(iface_inj_part_dist,           OPP_WRITE, OPP_Map_to_Mesh_Rel),  // iface_inj_part_dist,
                opp_arg_gbl(&(remainder),     1, "double", OPP_RW)                           // remainder global,
            );

            int old_nparts = particles_set->size;
            opp_par_loop_inject__InjectIons(
                inject_ions__kernel, "inject_ions__kernel",
                particles_set, OPP_ITERATE_INJECTED,                                                          // particles_set
                opp_arg_dat(part_position,                             OPP_WRITE),                          // part_position,
                opp_arg_dat(part_velocity,                             OPP_WRITE),                          // part_velocity,
                opp_arg_dat(part_mesh_relation,                        OPP_RW),                             // part_cell_connectivity,
                opp_arg_dat(iface_to_cell_map,                         OPP_READ, OPP_Map_from_Mesh_Rel),    // iface to cell map
                opp_arg_dat(cell_electric_field, 0, iface_to_cell_map, OPP_READ, OPP_Map_from_Mesh_Rel),    // cell_ef,
                opp_arg_dat(iface_u_normal,                            OPP_READ, OPP_Map_from_Mesh_Rel),    // iface_u,
                opp_arg_dat(iface_v_normal,                            OPP_READ, OPP_Map_from_Mesh_Rel),    // iface_v,
                opp_arg_dat(iface_normal,                              OPP_READ, OPP_Map_from_Mesh_Rel),    // iface_normal,
                opp_arg_dat(iface_node_pos,                            OPP_READ, OPP_Map_from_Mesh_Rel),    // iface_node_pos
                opp_arg_dat(dum_part_random,                           OPP_READ, OPP_Map_from_Inj_part)     // dum_part_random
            );

            opp_reset_dat(node_charge_density, (char*)opp_zero_double16);
            opp_par_loop_particle_all__MoveToCells(
                move_all_particles_to_cell__kernel, "move_all_particles_to_cell__kernel",
                particles_set, OPP_ITERATE_ALL,                                                                // particles_set
                opp_arg_dat(cell_electric_field,                       OPP_READ, OPP_Map_from_Mesh_Rel),     // cell_ef,
                opp_arg_dat(part_position,                             OPP_RW),                              // part_pos,
                opp_arg_dat(part_velocity,                             OPP_RW),                              // part_vel,
                opp_arg_dat(part_lc,                                   OPP_RW),                              // part_lc,
                opp_arg_dat(part_mesh_relation,                        OPP_RW),                              // current_cell_index,
                opp_arg_dat(cell_volume,                               OPP_READ, OPP_Map_from_Mesh_Rel),     // current_cell_volume,
                opp_arg_dat(cell_determinants,                         OPP_READ, OPP_Map_from_Mesh_Rel),     // current_cell_det,
                opp_arg_dat(cell_to_cell_map,                          OPP_READ, OPP_Map_from_Mesh_Rel),     // cell_connectivity,
                opp_arg_dat(node_charge_density, 0, cell_to_nodes_map, OPP_INC,  OPP_Map_from_Mesh_Rel),     // node_charge_den0,
                opp_arg_dat(node_charge_density, 1, cell_to_nodes_map, OPP_INC,  OPP_Map_from_Mesh_Rel),     // node_charge_den1,
                opp_arg_dat(node_charge_density, 2, cell_to_nodes_map, OPP_INC,  OPP_Map_from_Mesh_Rel),     // node_charge_den2,
                opp_arg_dat(node_charge_density, 3, cell_to_nodes_map, OPP_INC,  OPP_Map_from_Mesh_Rel)      // node_charge_den3,
            );

            opp_par_loop_all__ComputeNodeChargeDensity(
                compute_node_charge_density__kernel, "compute_node_charge_density__kernel",
                nodes_set, OPP_ITERATE_ALL,                       // nodes_set
                opp_arg_dat(node_charge_density,  OPP_RW),      // node_charge_density
                opp_arg_dat(node_volume,          OPP_READ)     // node_volume
            );

            mesh.solver->computePhi(  // TODO: Change this to kernel calls
                mesh.fesolver_method,
                opp_arg_dat(node_potential,      OPP_WRITE),
                opp_arg_dat(node_charge_density, OPP_READ)
            );

            opp_reset_dat(cell_electric_field, (char*)opp_zero_double16); 
            opp_par_loop_all__ComputeElectricField(
                compute_electric_field__kernel, "compute_electric_field__kernel",
                cells_set, OPP_ITERATE_ALL,                                        // cells_set
                opp_arg_dat(cell_electric_field,                  OPP_INC),      // cell_electric_field,
                opp_arg_dat(cell_shape_deriv,                     OPP_READ),     // cell_shape_deriv,
                opp_arg_dat(node_potential, 0, cell_to_nodes_map, OPP_READ),     // node_potential0,
                opp_arg_dat(node_potential, 1, cell_to_nodes_map, OPP_READ),     // node_potential1,
                opp_arg_dat(node_potential, 2, cell_to_nodes_map, OPP_READ),     // node_potential2,
                opp_arg_dat(node_potential, 3, cell_to_nodes_map, OPP_READ)      // node_potential3,
            );

            {
                double max_den = 0.0, max_phi = 0.0;
                if (params.get<BOOL>("check_max_values"))
                {
                    for (int n = 0; n< mesh.n_nodes; n++) // ideally, need to copy data from device to host, but at this point host has correct data
                    {
                        if (((double*)node_charge_density->data)[n] > max_den) max_den = ((double*)node_charge_density->data)[n];
                        if (abs(((double*)node_potential->data)[n]) > max_phi) max_phi = abs(((double*)node_potential->data)[n]);   
                    }
                }

                std::cout<<"ts: " << ts << "\t np: " << particles_set->size 
                    << " (" <<  injected_count << " added, "<< old_nparts - particles_set->size << " removed)"                
                    << "\t max den: " << max_den << "\t max |phi|: " << max_phi << "\t ****" << std::endl;
            }
        }

        std::chrono::duration<double> diff   = std::chrono::system_clock::now() - start;
        std::chrono::duration<double> diff_1 = std::chrono::system_clock::now() - start_iter1;
        std::cout << "\nFEMPIC - Time to iterate ts=" << max_iter << " <sec>: " << diff.count() << " and (ts-1) <sec>:" << diff_1.count() << "\n\n";

        opp_exit();
    } // End Scope for oppic

    return 0;
}

//*************************************************************************************************
// std::string f = std::string("F_") + std::to_string(ts + 1);
// opp_print_map_to_txtfile(cell_to_nodes_map  , f.c_str(), "cell_to_nodes_map.dat");
// opp_print_dat_to_txtfile(node_charge_density, f.c_str(), "node_charge_density.dat");
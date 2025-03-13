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

OPP_REAL CONST_spwt[1];
OPP_REAL CONST_ion_velocity[1];
OPP_REAL CONST_dt[1];
OPP_REAL CONST_plasma_den[1];
OPP_REAL CONST_mass[1];
OPP_REAL CONST_charge[1];
OPP_REAL CONST_wall_potential[1];

#include "fempic_misc_mesh_loader.h"
#include "fempic_misc.h"
#include "kernels.h"
#include "field_solver.h"

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    opp_init(argc, argv);

    {
    opp_profiler->start("Setup");

        OPP_REAL plasma_den            = opp_params->get<OPP_REAL>("plasma_den");
        OPP_REAL dt                    = opp_params->get<OPP_REAL>("dt");
        OPP_REAL ion_velocity          = opp_params->get<OPP_REAL>("ion_velocity");
        OPP_REAL spwt                  = opp_params->get<OPP_REAL>("spwt");
        OPP_REAL wall_potential        = opp_params->get<OPP_REAL>("wall_potential");
        OPP_REAL grid_spacing          = opp_params->get<OPP_REAL>("grid_spacing");
        OPP_REAL mass                  = 2 * AMU;
        OPP_REAL charge                = opp_params->get<OPP_INT>("charge_multiple") * QE;
        OPP_INT max_iter               = opp_params->get<OPP_INT>("num_steps");   
        std::string log                = "";
        const OPP_BOOL print_final_log = opp_params->get<OPP_BOOL>("print_final");
        int64_t total_part_iter        = 0;
        const opp_point expansion(4*grid_spacing, 4*grid_spacing, 4*grid_spacing);

        std::shared_ptr<DataPointers> m = load_mesh();

        opp_set node_set       = opp_decl_set(m->n_nodes, "mesh_nodes");
        opp_set cell_set       = opp_decl_set(m->n_cells, "mesh_cells");
        opp_set iface_set      = opp_decl_set(m->n_ifaces, "inlet_faces_cells");
        opp_set particle_set   = opp_decl_particle_set("particles", cell_set); 
        opp_set dummy_part_set = opp_decl_particle_set("dummy particles", cell_set); 

        opp_map c2n_map  = opp_decl_map(cell_set,  node_set, N_PER_C,  m->c_to_n, "c_v_n_map");
        opp_map c2c_map  = opp_decl_map(cell_set,  cell_set, NEIGHB_C, m->c_to_c,  "c_v_c_map"); 
        opp_map if2c_map = opp_decl_map(iface_set, cell_set, ONE,      m->if_to_c, "if_v_c_map"); 
        opp_map if2n_map = opp_decl_map(iface_set, node_set, N_PER_IF, m->if_to_n, "if_v_n_map");

        opp_dat c_det       = opp_decl_dat(cell_set, ALL_DET,     DT_REAL, m->c_det,      "c_det");  
        opp_dat c_volume    = opp_decl_dat(cell_set, ONE,         DT_REAL, m->c_vol,      "c_volume");        
        opp_dat c_ef        = opp_decl_dat(cell_set, DIM,         DT_REAL, m->c_ef,       "c_ef");
        opp_dat c_sd        = opp_decl_dat(cell_set, N_PER_C*DIM, DT_REAL, m->c_sd,       "c_shape_deri"); 
        opp_dat c_gbl_id    = opp_decl_dat(cell_set, ONE,         DT_INT,  m->c_id,       "c_gbl_id"); 
        opp_dat c_colors    = opp_decl_dat(cell_set, ONE,         DT_INT,  m->c_col,      "c_colors");
        opp_dat c_centroids = opp_decl_dat(cell_set, DIM,         DT_REAL, m->c_centroid, "c_centroids");

        opp_dat n_volume     = opp_decl_dat(node_set, ONE, DT_REAL, m->n_vol,     "n_vol");        
        opp_dat n_potential  = opp_decl_dat(node_set, ONE, DT_REAL, m->n_pot,     "n_potential");     
        opp_dat n_charge_den = opp_decl_dat(node_set, ONE, DT_REAL, m->n_ion_den, "n_charge_den");
        opp_dat n_pos        = opp_decl_dat(node_set, DIM, DT_REAL, m->n_pos,     "n_pos");     
        opp_dat n_type       = opp_decl_dat(node_set, ONE, DT_INT,  m->n_type,    "n_type");
        opp_dat n_bnd_pot    = opp_decl_dat(node_set, ONE, DT_REAL, m->n_bnd_pot, "n_bnd_pot");

        opp_dat if_v_norm  = opp_decl_dat(iface_set, DIM,          DT_REAL, m->if_v_norm, "iface_v_norm");
        opp_dat if_u_norm  = opp_decl_dat(iface_set, DIM,          DT_REAL, m->if_u_norm, "iface_u_norm");
        opp_dat if_norm    = opp_decl_dat(iface_set, DIM,          DT_REAL, m->if_norm,   "iface_norm");  
        opp_dat if_area    = opp_decl_dat(iface_set, ONE,          DT_REAL, m->if_area,   "iface_area");
        opp_dat if_distrib = opp_decl_dat(iface_set, ONE,          DT_INT,  m->if_dist,   "iface_dist");
        opp_dat if_n_pos   = opp_decl_dat(iface_set, N_PER_IF*DIM, DT_REAL, m->if_n_pos,  "iface_n_pos");

        opp_dat p_pos   = opp_decl_dat(particle_set, DIM,     DT_REAL, nullptr, "p_position");
        opp_dat p_vel   = opp_decl_dat(particle_set, DIM,     DT_REAL, nullptr, "p_velocity");
        opp_dat p_lc    = opp_decl_dat(particle_set, N_PER_C, DT_REAL, nullptr, "p_lc");

        opp_map p2c_map = opp_decl_map(particle_set, cell_set, ONE, nullptr, "p2c_map");
        
        opp_dat dp_rand = opp_decl_dat(dummy_part_set, 2, DT_REAL, nullptr, "dummy_part_rand");

        opp_decl_const<OPP_REAL>(ONE, &spwt,           "CONST_spwt");
        opp_decl_const<OPP_REAL>(ONE, &ion_velocity,   "CONST_ion_velocity");
        opp_decl_const<OPP_REAL>(ONE, &dt,             "CONST_dt");
        opp_decl_const<OPP_REAL>(ONE, &plasma_den,     "CONST_plasma_den");
        opp_decl_const<OPP_REAL>(ONE, &mass,           "CONST_mass");
        opp_decl_const<OPP_REAL>(ONE, &charge,         "CONST_charge");
        opp_decl_const<OPP_REAL>(ONE, &wall_potential, "CONST_wall_potential");

        m->DeleteValues();

#ifdef USE_MPI
        fempic_color_block(c_colors, c_centroids, if_n_pos, if2n_map);

        // opp_partition(std::string("PARMETIS_KWAY"), cell_set, c2n_map);
        // opp_partition(std::string("PARMETIS_GEOM"), cell_set, nullptr, n_pos);
        opp_partition(std::string("EXTERNAL"), cell_set, nullptr, c_colors);
#endif
        
        opp_par_loop(init_boundary_pot_kernel, "init_boundary_pot", node_set, OPP_ITERATE_ALL,
            opp_arg_dat(n_type,    OPP_READ), 
            opp_arg_dat(n_bnd_pot, OPP_WRITE));

        const int inject_count = init_inject_distributions(if_distrib, if_area, dp_rand);

        opp_inc_part_count_with_distribution(particle_set, inject_count, if_distrib, false);

        // these two lines are only required if we plan to use direct_hop
        opp::BoundingBox bounding_box(n_pos, DIM, expansion);
        opp_init_direct_hop(grid_spacing, c_gbl_id, bounding_box);

        auto field_solver = std::make_unique<FESolver>(c2n_map, n_type, n_pos, n_bnd_pot);
        field_solver->enrich_cell_shape_deriv(c_sd);

    opp_profiler->end("Setup");

    opp_profiler->start("MainLoop");
        for (OPP_main_loop_iter = 0; OPP_main_loop_iter < max_iter; OPP_main_loop_iter++)
        {
            OPP_RUN_ON_ROOT(&& OPP_DBG) opp_printf("Main", "Loop %d start ***", OPP_main_loop_iter);

            if (OPP_main_loop_iter != 0)
                opp_inc_part_count_with_distribution(particle_set, inject_count, if_distrib, false);

            const int old_nparts = particle_set->size;
            opp_par_loop(inject_ions_kernel, "inject_ions", particle_set, OPP_ITERATE_INJECTED,
                opp_arg_dat(p_pos,                      OPP_WRITE),                      
                opp_arg_dat(p_vel,                      OPP_WRITE),                      
                opp_arg_dat(p2c_map,                    OPP_RW),
                opp_arg_dat(if2c_map,          p2c_map, OPP_READ),
                opp_arg_dat(c_ef, 0, if2c_map, p2c_map, OPP_READ),
                opp_arg_dat(if_u_norm,         p2c_map, OPP_READ),
                opp_arg_dat(if_v_norm,         p2c_map, OPP_READ),
                opp_arg_dat(if_norm,           p2c_map, OPP_READ),
                opp_arg_dat(if_n_pos,          p2c_map, OPP_READ),
                opp_arg_dat(dp_rand,                    OPP_READ, false)); // offset=false flag will iterate this from begining

            opp_reset_dat(n_charge_den, opp_zero_double16);

            opp_par_loop(calculate_new_pos_vel_kernel, "calculate_new_pos_vel", particle_set, OPP_ITERATE_ALL,
                opp_arg_dat(c_ef, p2c_map, OPP_READ),
                opp_arg_dat(p_pos,         OPP_WRITE),
                opp_arg_dat(p_vel,         OPP_WRITE));

            opp_particle_move(move_kernel, "move", particle_set, c2c_map, p2c_map,
                opp_arg_dat(p_pos,             OPP_READ),
                opp_arg_dat(p_lc,              OPP_WRITE),
                opp_arg_dat(c_volume, p2c_map, OPP_READ),
                opp_arg_dat(c_det,    p2c_map, OPP_READ));

            opp_par_loop(deposit_charge_on_nodes_kernel, "deposit_charge_on_nodes", particle_set, OPP_ITERATE_ALL,
                opp_arg_dat(p_lc,                              OPP_READ),
                opp_arg_dat(n_charge_den, 0, c2n_map, p2c_map, OPP_INC),
                opp_arg_dat(n_charge_den, 1, c2n_map, p2c_map, OPP_INC),
                opp_arg_dat(n_charge_den, 2, c2n_map, p2c_map, OPP_INC),
                opp_arg_dat(n_charge_den, 3, c2n_map, p2c_map, OPP_INC));

            opp_par_loop(compute_node_charge_density_kernel, "compute_node_charge_density", node_set, OPP_ITERATE_ALL,
                opp_arg_dat(n_charge_den,  OPP_RW), 
                opp_arg_dat(n_volume,      OPP_READ));

            field_solver->compute_phi(  // TODO: Change this to kernel calls
                opp_arg_dat(n_potential,  OPP_WRITE),
                opp_arg_dat(n_charge_den, OPP_READ),
                opp_arg_dat(n_bnd_pot,    OPP_READ));

            opp_reset_dat(c_ef, opp_zero_double16); 

            opp_par_loop(compute_electric_field_kernel, "compute_electric_field", cell_set, OPP_ITERATE_ALL,
                opp_arg_dat(c_ef,                    OPP_INC), 
                opp_arg_dat(c_sd,                    OPP_READ),
                opp_arg_dat(n_potential, 0, c2n_map, OPP_READ),
                opp_arg_dat(n_potential, 1, c2n_map, OPP_READ),
                opp_arg_dat(n_potential, 2, c2n_map, OPP_READ),
                opp_arg_dat(n_potential, 3, c2n_map, OPP_READ));

            if (print_final_log)
            {
                OPP_REAL max_n_chg_den = 0.0, max_n_pot = 0.0, max_c_ef = 0.0;

                opp_par_loop(get_max_cef_kernel, "get_max_cef", cell_set, OPP_ITERATE_ALL,
                    opp_arg_dat(c_ef, OPP_READ),
                    opp_arg_gbl(&max_c_ef, 1, "double", OPP_MAX));
                opp_par_loop(get_final_max_values_kernel, "get_final_max_values", node_set, OPP_ITERATE_ALL,
                    opp_arg_dat(n_charge_den, OPP_READ),
                    opp_arg_gbl(&max_n_chg_den, 1, "double", OPP_MAX),
                    opp_arg_dat(n_potential, OPP_READ),
                    opp_arg_gbl(&max_n_pot, 1, "double", OPP_MAX));

                log = get_global_level_log(max_c_ef, max_n_pot, particle_set->size, inject_count, 
                    (old_nparts - particle_set->size));
            }

            total_part_iter += particle_set->size;  

            OPP_RUN_ON_ROOT() opp_printf("Main", "ts: %d %s ****", OPP_main_loop_iter, log.c_str());
        }
    opp_profiler->end("MainLoop");

        const int64_t global_parts_iterated = opp_get_global_value(total_part_iter);
        OPP_RUN_ON_ROOT()
            opp_printf("Main", "Loop completed : %d iterations with %" PRId64 " particle iterations ****", 
                max_iter, global_parts_iterated);  
    }

    opp_exit();

    return 0;
}

//*************************************************************************************************
// std::string f = std::string("F_") + std::to_string(OPP_main_loop_iter + 1);
// opp_print_map_to_txtfile(c2n_map  , f.c_str(), "c2n_map.dat");
// opp_print_dat_to_txtfile(n_charge_den, f.c_str(), "n_charge_den.dat");
// opp_mpi_print_dat_to_txtfile(c_sd, "c_sd.dat");
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

#include "fempic_misc_mesh_loader.h"
#include "fempic_misc.h"
#include "FESolver.h"

void opp_loop_all__init_boundary_pot(opp_set,opp_arg,opp_arg);
void opp_loop_inject__inject_ions(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,
    opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__compute_node_charge_density(opp_set,opp_arg,opp_arg);
void opp_loop_all__compute_electric_field(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__calc_pos_vel(opp_set,opp_arg,opp_arg,opp_arg);
void opp_loop_all__deposit_charge_on_nodes(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_particle_move__move(opp_set,opp_map,opp_map,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__get_max_values(opp_set,opp_arg,opp_arg,opp_arg,opp_arg);

void init_particle_mover(const double gridSpacing, int dim, const opp_dat n_pos_dat, 
    const opp_dat cellVolume_dat, const opp_dat cellDet_dat, const opp_dat c_gbl_id_dat);

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    opp_init(argc, argv);

    {
        opp_profiler->start("Setup");

        double plasma_den     = opp_params->get<OPP_REAL>("plasma_den");
        double dt             = opp_params->get<OPP_REAL>("dt");
        double ion_velocity   = opp_params->get<OPP_REAL>("ion_velocity");
        double spwt           = opp_params->get<OPP_REAL>("spwt");
        double wall_potential = opp_params->get<OPP_REAL>("wall_potential");
        double grid_spacing   = opp_params->get<OPP_REAL>("grid_spacing");
        double mass           = 2 * AMU;
        double charge         = 1 * QE;
        int max_iter          = opp_params->get<OPP_INT>("max_iter");   
        std::string log       = "";
        const bool print_final_log = opp_params->get<OPP_BOOL>("print_final");
        int64_t total_part_iter = 0;

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

        opp_dat p_pos   = opp_decl_dat(particle_set, DIM,     DT_REAL, nullptr, "part_position");
        opp_dat p_vel   = opp_decl_dat(particle_set, DIM,     DT_REAL, nullptr, "part_velocity");
        opp_dat p_lc    = opp_decl_dat(particle_set, N_PER_C, DT_REAL, nullptr, "part_lc");

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
        // opp_partition(std::string("PARMETIS_GEOM"), iface_set, nullptr, if_n_pos);
        // opp_partition(std::string("EXTERNAL"), node_set, nullptr, node_colors);
        opp_partition(std::string("EXTERNAL"), cell_set, nullptr, c_colors);
#endif
        
        opp_loop_all__init_boundary_pot(node_set, 
            opp_arg_dat(n_type,    OPP_READ), 
            opp_arg_dat(n_bnd_pot, OPP_WRITE));

        const int inject_count = init_inject_distributions(if_distrib, if_area, dp_rand);

        init_particle_mover(grid_spacing, DIM, n_pos, c_volume, c_det, c_gbl_id);

        std::unique_ptr<FESolver> field_solver = std::make_unique<FESolver>(c2n_map, n_type,
                                                            n_pos, n_bnd_pot, argc, argv);
            
        field_solver->enrich_cell_shape_deriv(c_sd);

        opp_inc_part_count_with_distribution(particle_set, inject_count, if_distrib, false);

        opp_profiler->end("Setup");

#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        opp_profiler->start("MainLoop");
        for (OPP_main_loop_iter = 0; OPP_main_loop_iter < max_iter; OPP_main_loop_iter++)
        {
            OPP_RUN_ON_ROOT(&& OPP_DBG) opp_printf("Main", "Loop %d start ***", OPP_main_loop_iter);

            if (OPP_main_loop_iter != 0)
                opp_inc_part_count_with_distribution(particle_set, inject_count, if_distrib, false);

            int old_nparts = particle_set->size;
            opp_loop_inject__inject_ions(particle_set,                                                                           
                opp_arg_dat(p_pos,                      OPP_WRITE),                      
                opp_arg_dat(p_vel,                      OPP_WRITE),                      
                opp_arg_dat(p2c_map,                    OPP_RW),
                opp_arg_dat(p2c_map,                    OPP_RW), // TODO: remove
                opp_arg_dat(if2c_map,          p2c_map, OPP_READ),
                opp_arg_dat(c_ef, 0, if2c_map, p2c_map, OPP_READ),
                opp_arg_dat(if_u_norm,         p2c_map, OPP_READ),
                opp_arg_dat(if_v_norm,         p2c_map, OPP_READ),
                opp_arg_dat(if_norm,           p2c_map, OPP_READ),
                opp_arg_dat(if_n_pos,          p2c_map, OPP_READ),
                opp_arg_dat(dp_rand,                    OPP_READ));

            opp_reset_dat(n_charge_den, opp_zero_double16);

            opp_loop_all__calc_pos_vel(particle_set,                                                                           
                opp_arg_dat(c_ef, p2c_map, OPP_READ),
                opp_arg_dat(p_pos,         OPP_WRITE),                         
                opp_arg_dat(p_vel,         OPP_WRITE));

            opp_particle_move__move(particle_set, c2c_map, p2c_map,
                opp_arg_dat(p_pos,             OPP_READ),
                opp_arg_dat(p_lc,              OPP_WRITE),
                opp_arg_dat(c_volume, p2c_map, OPP_READ),
                opp_arg_dat(c_det,    p2c_map, OPP_READ));

            opp_loop_all__deposit_charge_on_nodes(particle_set, 
                opp_arg_dat(p_lc,                              OPP_READ),
                opp_arg_dat(n_charge_den, 0, c2n_map, p2c_map, OPP_INC),
                opp_arg_dat(n_charge_den, 1, c2n_map, p2c_map, OPP_INC),
                opp_arg_dat(n_charge_den, 2, c2n_map, p2c_map, OPP_INC),
                opp_arg_dat(n_charge_den, 3, c2n_map, p2c_map, OPP_INC));

            opp_loop_all__compute_node_charge_density(node_set,                            
                opp_arg_dat(n_charge_den,  OPP_RW), 
                opp_arg_dat(n_volume,      OPP_READ));

            field_solver->computePhi(  // TODO: Change this to kernel calls
                opp_arg_dat(n_potential,  OPP_WRITE),
                opp_arg_dat(n_charge_den, OPP_READ),
                opp_arg_dat(n_bnd_pot,    OPP_READ));

            opp_reset_dat(c_ef, opp_zero_double16); 

            opp_loop_all__compute_electric_field(cell_set,                                                   
                opp_arg_dat(c_ef,                    OPP_INC), 
                opp_arg_dat(c_sd,                    OPP_READ),
                opp_arg_dat(n_potential, 0, c2n_map, OPP_READ),
                opp_arg_dat(n_potential, 1, c2n_map, OPP_READ),
                opp_arg_dat(n_potential, 2, c2n_map, OPP_READ),
                opp_arg_dat(n_potential, 3, c2n_map, OPP_READ));

            if (print_final_log)
            {
                OPP_REAL max_n_chg_den = 0.0, max_n_pot = 0.0;

                opp_loop_all__get_max_values(node_set,
                    opp_arg_dat(n_charge_den, OPP_READ),
                    opp_arg_gbl(&max_n_chg_den, 1, "double", OPP_MAX),
                    opp_arg_dat(n_potential, OPP_READ),
                    opp_arg_gbl(&max_n_pot, 1, "double", OPP_MAX));

                log = get_global_level_log(max_n_chg_den, max_n_pot, particle_set->size, 
                    inject_count, (old_nparts - particle_set->size));
            }

            total_part_iter += particle_set->size;  

            OPP_RUN_ON_ROOT() opp_printf("Main", "ts: %d %s ****", OPP_main_loop_iter, log.c_str());
        }
        opp_profiler->end("MainLoop");

        const int64_t global_parts_iterated = get_global_parts_iterated(total_part_iter);
        OPP_RUN_ON_ROOT()
            opp_printf("Main", "Loop completed : %d iterations with %" PRId64 " particle iterations ****", 
                max_iter, global_parts_iterated);  
    }

    opp_exit();

    return 0;
}

//*************************************************************************************************
// std::string f = std::string("F_") + std::to_string(ts + 1);
// opp_print_map_to_txtfile(c2n_map  , f.c_str(), "c2n_map.dat");
// opp_print_dat_to_txtfile(n_charge_den, f.c_str(), "n_charge_den.dat");
// opp_mpi_print_dat_to_txtfile(c_sd, "c_sd.dat");
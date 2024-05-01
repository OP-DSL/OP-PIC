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

#include "cabana_defs.h"

void opp_loop_all__interpolate_mesh_fields(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_particle_move__move_deposit(opp_set,opp_map,opp_dat,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__accumulate_current_to_cells(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__half_advance_b(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__advance_e(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__get_max_values(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__compute_energy(opp_set,opp_arg,opp_arg,opp_arg);
void opp_loop_all__update_ghosts_B(opp_set,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_loop_all__update_ghosts(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);

using namespace opp;

#include "cabana_misc.cpp"

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << std::endl;
        exit(-1);
    }

    opp_init(argc, argv);
    opp_params->write(std::cout);

    {
        opp_profiler->start("Setup");

        Deck deck(*opp_params);
        if (OPP_rank == OPP_ROOT) deck.print();

        std::shared_ptr<DataPointers> m = load_mesh(deck);

        opp_set c_set       = opp_decl_set(m->n_cells, "mesh_cells");
        opp_set p_set        = opp_decl_particle_set("particles", c_set);  // Zero particles, inject after partitioning

        opp_map c2c_map     = opp_decl_map(c_set, c_set, NEIGHBOURS, m->c2c_map,     "c2c_map");
        opp_map c2ngc_map   = opp_decl_map(c_set, c_set, FACES,      m->c2ngc_map,   "c2c_non_ghost_map");
        opp_map c2cug0_map  = opp_decl_map(c_set, c_set, 1,          m->c2cug0_map,  "c2c_ug_0_map");
        opp_map c2cug1_map  = opp_decl_map(c_set, c_set, 1,          m->c2cug1_map,  "c2c_ug_1_map");
        opp_map c2cug2_map  = opp_decl_map(c_set, c_set, 1,          m->c2cug2_map,  "c2c_ug_2_map");
        opp_map c2cug3_map  = opp_decl_map(c_set, c_set, 1,          m->c2cug3_map,  "c2c_ug_3_map");
        opp_map c2cug4_map  = opp_decl_map(c_set, c_set, 1,          m->c2cug4_map,  "c2c_ug_4_map");
        opp_map c2cug5_map  = opp_decl_map(c_set, c_set, 1,          m->c2cug5_map,  "c2c_ug_5_map");
        opp_map c2cugb0_map = opp_decl_map(c_set, c_set, 1,          m->c2cugb0_map, "c2c_ugb_0_map");
        opp_map c2cugb1_map = opp_decl_map(c_set, c_set, 1,          m->c2cugb1_map, "c2c_ugb_1_map");
        opp_map c2cugb2_map = opp_decl_map(c_set, c_set, 1,          m->c2cugb2_map, "c2c_ugb_2_map");
        opp_map c2cugb3_map = opp_decl_map(c_set, c_set, 1,          m->c2cugb3_map, "c2c_ugb_3_map");
        opp_map c2cugb4_map = opp_decl_map(c_set, c_set, 1,          m->c2cugb4_map, "c2c_ugb_4_map");
        opp_map c2cugb5_map = opp_decl_map(c_set, c_set, 1,          m->c2cugb5_map, "c2c_ugb_5_map");

        opp_dat c_index      = opp_decl_dat(c_set, ONE,        DT_INT,  m->c_index,      "c_index");
        opp_dat c_e          = opp_decl_dat(c_set, DIM,        DT_REAL, m->c_e,          "c_e");  
        opp_dat c_b          = opp_decl_dat(c_set, DIM,        DT_REAL, m->c_b,          "c_b");        
        opp_dat c_j          = opp_decl_dat(c_set, DIM,        DT_REAL, m->c_j,          "c_j");
        opp_dat c_acc        = opp_decl_dat(c_set, ACC_LEN,    DT_REAL, m->c_acc,        "c_acc"); 
        opp_dat c_interp     = opp_decl_dat(c_set, INTERP_LEN, DT_REAL, m->c_interp,     "c_interp"); 
        opp_dat c_pos_ll     = opp_decl_dat(c_set, DIM,        DT_REAL, m->c_pos_ll,     "c_pos_ll");       
        opp_dat c_colors     = opp_decl_dat(c_set, ONE,        DT_INT,  m->c_colours,    "c_colors"); // init dummy
        opp_dat c_mask_ghost = opp_decl_dat(c_set, ONE,        DT_INT,  m->c_ghost,      "c_mask_ghost");
        opp_dat c_mask_right = opp_decl_dat(c_set, ONE,        DT_INT,  m->c_mask_right, "c_mask_right");
        opp_dat c_mask_ug    = opp_decl_dat(c_set, 6,          DT_INT,  m->c_mask_ug,    "c_mask_ug");
        opp_dat c_mask_ugb   = opp_decl_dat(c_set, 6,          DT_INT,  m->c_mask_ugb,   "c_mask_ugb");
   
        opp_dat p_index      = nullptr; // opp_decl_dat(p_set, ONE, DT_INT,  nullptr, "p_index"); // Unused
        opp_dat p_pos        = opp_decl_dat(p_set, DIM, DT_REAL, nullptr, "p_pos");
        opp_dat p_vel        = opp_decl_dat(p_set, DIM, DT_REAL, nullptr, "p_vel");    
        opp_dat p_streak_mid = opp_decl_dat(p_set, DIM, DT_REAL, nullptr, "p_streak_mid");
        opp_dat p_weight     = opp_decl_dat(p_set, ONE, DT_REAL, nullptr, "p_weight");
        opp_dat p2c_map      = opp_decl_dat(p_set, ONE, DT_INT,  nullptr, "p2c_map", true);

        OPP_REAL cdt_d[DIM]        = { deck.cdt_dx, deck.cdt_dy, deck.cdt_dz };
        OPP_INT cells_per_dim[DIM] = { deck.nx, deck.ny, deck.nz };
        OPP_REAL p[DIM]            = { deck.px, deck.py, deck.pz };
        OPP_REAL acc_coef[DIM]     = { deck.acc_coefx, deck.acc_coefy, deck.acc_coefz };

        opp_decl_const<OPP_REAL>(ONE, &(deck.dt),      "CONST_dt");
        opp_decl_const<OPP_REAL>(ONE, &(deck.qsp),     "CONST_qsp");
        opp_decl_const<OPP_REAL>(DIM, cdt_d,           "CONST_cdt_d");
        opp_decl_const<OPP_REAL>(DIM, p,               "CONST_p");
        opp_decl_const<OPP_REAL>(ONE, &(deck.qdt_2mc), "CONST_qdt_2mc");
        opp_decl_const<OPP_REAL>(ONE, &(deck.dt_eps0), "CONST_dt_eps0");
        opp_decl_const<OPP_REAL>(DIM, acc_coef,        "CONST_acc_coef");

        m->DeleteValues();

        // opp_colour_cartesian_mesh does nothing for non-mpi runs
        if (opp_params->get<OPP_STRING>("cluster") == "block")
            cabana_color_block_x(deck, c_index, c_colors);
        else if (opp_params->get<OPP_STRING>("cluster") == "pencil")
            cabana_color_pencil_x(deck, c_index, c_colors);
        else if (opp_params->get<OPP_STRING>("cluster") == "cart")
            opp_colour_cartesian_mesh(DIM, std::vector<OPP_INT>(cells_per_dim, cells_per_dim + DIM), 
                                        c_index, c_colors, 1);
        else
            opp_abort("cluster type not supported");

#ifdef USE_MPI
        opp_partition(std::string("EXTERNAL"), c_set, nullptr, c_colors);
#endif

        init_particles(deck, p_index, p_pos, p_vel, p_streak_mid, p_weight, p2c_map, 
                        c_pos_ll, c_index, c_mask_ghost);

        FILE *fptre = nullptr;
        if (OPP_rank == OPP_ROOT) fptre = fopen("energy.csv", "w");
        if (OPP_rank == OPP_ROOT) fprintf(fptre, "timestep,e_energy,b_energy\n");

        int ghosts = 0;
        for (int i = 0; i < c_set->size; i++) if (((int*)c_mask_ghost->data)[i] == 1) ghosts++;
        opp_printf("Setup", "Cells[%d][ghosts:%d] Particles[%d]", c_set->size, ghosts, p_set->size);
        
        opp_profiler->end("Setup");

#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        opp_profiler->start("MainLoop");

        for (OPP_main_loop_iter = 0; OPP_main_loop_iter < deck.num_steps; OPP_main_loop_iter++)
        {
            if (OPP_DBG && OPP_rank == OPP_ROOT) 
                opp_printf("Main", "Starting main loop iteration %d ****", OPP_main_loop_iter);

            opp_loop_all__interpolate_mesh_fields(c_set,
                opp_arg_dat(c_e,                            OPP_READ),
                opp_arg_dat(c_b,                            OPP_READ),
                opp_arg_dat(c_e, CellMap::xu_y_z,  c2c_map, OPP_READ),
                opp_arg_dat(c_e, CellMap::x_yu_z,  c2c_map, OPP_READ),
                opp_arg_dat(c_e, CellMap::x_y_zu,  c2c_map, OPP_READ),
                opp_arg_dat(c_e, CellMap::xu_yu_z, c2c_map, OPP_READ),
                opp_arg_dat(c_e, CellMap::x_yu_zu, c2c_map, OPP_READ),
                opp_arg_dat(c_e, CellMap::xu_y_zu, c2c_map, OPP_READ),
                opp_arg_dat(c_b, CellMap::xu_y_z,  c2c_map, OPP_READ),
                opp_arg_dat(c_b, CellMap::x_yu_z,  c2c_map, OPP_READ),
                opp_arg_dat(c_b, CellMap::x_y_zu,  c2c_map, OPP_READ),
                opp_arg_dat(c_interp,                       OPP_WRITE),
                opp_arg_dat(c_mask_ghost,                   OPP_READ));

            opp_reset_dat(c_acc, (char*)opp_zero_double16);

            opp_particle_move__move_deposit(p_set, c2ngc_map, p2c_map,
                opp_arg_dat(p_vel,             OPP_RW),
                opp_arg_dat(p_pos,             OPP_RW),
                opp_arg_dat(p_streak_mid,      OPP_RW),
                opp_arg_dat(p_weight,          OPP_READ),
                opp_arg_dat(c_interp, p2c_map, OPP_READ),
                opp_arg_dat(c_acc,    p2c_map, OPP_INC));

            opp_loop_all__accumulate_current_to_cells(c_set,
                opp_arg_dat(c_j,                              OPP_WRITE),
                opp_arg_dat(c_acc,                            OPP_READ),
                opp_arg_dat(c_acc, CellMap::xd_y_z , c2c_map, OPP_READ),
                opp_arg_dat(c_acc, CellMap::x_yd_z , c2c_map, OPP_READ),
                opp_arg_dat(c_acc, CellMap::x_y_zd , c2c_map, OPP_READ),
                opp_arg_dat(c_acc, CellMap::xd_yd_z, c2c_map, OPP_READ),
                opp_arg_dat(c_acc, CellMap::x_yd_zd, c2c_map, OPP_READ),
                opp_arg_dat(c_acc, CellMap::xd_y_zd, c2c_map, OPP_READ),
                opp_arg_dat(c_mask_right,                     OPP_READ));

            // Leap frog method
            opp_loop_all__half_advance_b(c_set,
                opp_arg_dat(c_e, CellMap::xu_y_z, c2c_map, OPP_READ),
                opp_arg_dat(c_e, CellMap::x_yu_z, c2c_map, OPP_READ), 
                opp_arg_dat(c_e, CellMap::x_y_zu, c2c_map, OPP_READ), 
                opp_arg_dat(c_e,                           OPP_READ),
                opp_arg_dat(c_b,                           OPP_INC),
                opp_arg_dat(c_mask_ghost,                  OPP_READ));

            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb0_map, OPP_WRITE), opp_arg_gbl(&m0, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb1_map, OPP_WRITE), opp_arg_gbl(&m1, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb2_map, OPP_WRITE), opp_arg_gbl(&m2, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb3_map, OPP_WRITE), opp_arg_gbl(&m3, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb4_map, OPP_WRITE), opp_arg_gbl(&m4, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb5_map, OPP_WRITE), opp_arg_gbl(&m5, 1, "int", OPP_READ));

            opp_loop_all__update_ghosts(c_set, opp_arg_dat(c_mask_ug, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cug0_map, OPP_INC), opp_arg_gbl(&m0, 1, "int", OPP_READ), 
                opp_arg_gbl(&m0, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts(c_set, opp_arg_dat(c_mask_ug, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cug1_map, OPP_INC), opp_arg_gbl(&m1, 1, "int", OPP_READ), 
                opp_arg_gbl(&m0, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts(c_set, opp_arg_dat(c_mask_ug, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cug2_map, OPP_INC), opp_arg_gbl(&m2, 1, "int", OPP_READ), 
                opp_arg_gbl(&m1, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts(c_set, opp_arg_dat(c_mask_ug, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cug3_map, OPP_INC), opp_arg_gbl(&m3, 1, "int", OPP_READ), 
                opp_arg_gbl(&m1, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts(c_set, opp_arg_dat(c_mask_ug, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cug4_map, OPP_INC), opp_arg_gbl(&m4, 1, "int", OPP_READ), 
                opp_arg_gbl(&m2, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts(c_set, opp_arg_dat(c_mask_ug, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cug5_map, OPP_INC), opp_arg_gbl(&m5, 1, "int", OPP_READ), 
                opp_arg_gbl(&m2, 1, "int", OPP_READ));

            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cugb0_map, OPP_WRITE), opp_arg_gbl(&m0, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cugb1_map, OPP_WRITE), opp_arg_gbl(&m1, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cugb2_map, OPP_WRITE), opp_arg_gbl(&m2, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cugb3_map, OPP_WRITE), opp_arg_gbl(&m3, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cugb4_map, OPP_WRITE), opp_arg_gbl(&m4, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_j, OPP_READ),
                opp_arg_dat(c_j, 0, c2cugb5_map, OPP_WRITE), opp_arg_gbl(&m5, 1, "int", OPP_READ));

            opp_loop_all__advance_e(c_set,
                opp_arg_dat(c_b, CellMap::xd_y_z, c2c_map, OPP_READ),
                opp_arg_dat(c_b, CellMap::x_yd_z, c2c_map, OPP_READ),
                opp_arg_dat(c_b, CellMap::x_y_zd, c2c_map, OPP_READ),
                opp_arg_dat(c_b,                           OPP_READ),
                opp_arg_dat(c_j,                           OPP_READ), 
                opp_arg_dat(c_e,                           OPP_INC),
                opp_arg_dat(c_mask_right,                  OPP_READ));

            opp_loop_all__half_advance_b(c_set,
                opp_arg_dat(c_e, CellMap::xu_y_z, c2c_map, OPP_READ),
                opp_arg_dat(c_e, CellMap::x_yu_z, c2c_map, OPP_READ), 
                opp_arg_dat(c_e, CellMap::x_y_zu, c2c_map, OPP_READ), 
                opp_arg_dat(c_e,                           OPP_READ),
                opp_arg_dat(c_b,                           OPP_INC),
                opp_arg_dat(c_mask_ghost,                  OPP_READ));

            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb0_map, OPP_WRITE), opp_arg_gbl(&m0, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb1_map, OPP_WRITE), opp_arg_gbl(&m1, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb2_map, OPP_WRITE), opp_arg_gbl(&m2, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb3_map, OPP_WRITE), opp_arg_gbl(&m3, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb4_map, OPP_WRITE), opp_arg_gbl(&m4, 1, "int", OPP_READ));
            opp_loop_all__update_ghosts_B(c_set, opp_arg_dat(c_mask_ugb, OPP_READ), opp_arg_dat(c_b, OPP_READ),
                opp_arg_dat(c_b, 0, c2cugb5_map, OPP_WRITE), opp_arg_gbl(&m5, 1, "int", OPP_READ));


            std::string log = ""; // TODO : print some unseful information to verify
            if (opp_params->get<OPP_BOOL>("print_energy"))
            {
                OPP_REAL e_energy = 0.0, b_energy = 0.0;
                opp_loop_all__compute_energy(c_set,
                    opp_arg_dat(c_mask_ghost, OPP_READ),
                    opp_arg_dat(c_e,          OPP_READ),
                    opp_arg_gbl(&e_energy, 1, "double", OPP_INC));
                opp_loop_all__compute_energy(c_set,
                    opp_arg_dat(c_mask_ghost, OPP_READ),
                    opp_arg_dat(c_b,          OPP_READ),
                    opp_arg_gbl(&b_energy, 1, "double", OPP_INC));
                log += str(e_energy*0.5, " e_energy: \t%.15f");
                log += str(b_energy*0.5, " b_energy: \t%.15f");

                if (OPP_rank == OPP_ROOT) 
                    fprintf(fptre,"%d, %.15f, %.15f\n", OPP_main_loop_iter, e_energy*0.5, b_energy*0.5);
            }
            if (opp_params->get<OPP_BOOL>("print_final"))
            {
                OPP_REAL max_j = 0.0, max_e = 0.0, max_b = 0.0;
                opp_loop_all__get_max_values(c_set, // plan is to get only the x values reduced here
                    opp_arg_dat(c_j, OPP_READ),
                    opp_arg_gbl(&max_j, 1, "double", OPP_MAX),
                    opp_arg_dat(c_e, OPP_READ),
                    opp_arg_gbl(&max_e, 1, "double", OPP_MAX),
                    opp_arg_dat(c_b, OPP_READ),
                    opp_arg_gbl(&max_b, 1, "double", OPP_MAX));
                log += str(max_j, " max_jx: %.10f"); // %2.15lE");
                log += str(max_e, " max_ex: %.10f"); // %2.15lE");
                log += str(max_b, " max_bx: %.16f"); // %2.15lE");
            }
            if (OPP_rank == OPP_ROOT) 
                opp_printf("Main", "ts: %d %s ****", OPP_main_loop_iter, log.c_str());
        }
        opp_profiler->end("MainLoop");
        
        if (OPP_rank == OPP_ROOT) fclose(fptre);

        if (OPP_rank == OPP_ROOT) 
            opp_printf("Main", "Main loop completed after %d iterations ****", deck.num_steps);
    }

    opp_exit();

    return 0;
}


// std::string f = std::string("F_") + std::to_string(OPP_main_loop_iter);    
// opp_print_map_to_txtfile(c_v_nodes_map  , f.c_str(), "c_v_nodes_map.dat");
// opp_print_dat_to_txtfile(p_pos, f.c_str(), "p_pos.dat");
// std::string f1 = midString + std::to_string(OPP_main_loop_iter) + "_c_e.dat";
// opp_mpi_print_dat_to_txtfile(c_e, f1.c_str());
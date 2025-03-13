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

#include "advec_defs.h"

using namespace opp;

void initialize_particle_mover(const double grid_spacing, int dim, 
    const opp_dat cell_pos_ll_dat, const opp_dat global_cid_dat, const opp_map cell_cell_map);

#ifdef FUSE_KERNELS   
void opp_particle_move__update_pos_move(opp_set,opp_arg,opp_arg,opp_arg,opp_arg,opp_arg);
#else
void opp_loop_all__update_pos(opp_set,opp_arg,opp_arg);
void opp_particle_move__move(opp_set,opp_arg,opp_arg,opp_arg,opp_arg);
#endif
void opp_loop_all__verify(opp_set,opp_arg,opp_arg,opp_arg,opp_arg);

#include "advec_misc.cpp"

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

        OPP_INT max_iter      = opp_params->get<OPP_INT>("max_iter");   
        OPP_REAL dt           = opp_params->get<OPP_REAL>("dt");
        OPP_REAL cell_width   = opp_params->get<OPP_REAL>("cell_width");
        OPP_REAL extents[2]   = {opp_params->get<OPP_INT>("nx")*cell_width, opp_params->get<OPP_INT>("ny")*cell_width};
        OPP_INT ndimcells[2]  = {opp_params->get<OPP_INT>("nx"), opp_params->get<OPP_INT>("ny")};
        OPP_BOOL verify_parts = opp_params->get<OPP_BOOL>("verify_particles");
        OPP_REAL grid_spacing = opp_params->get<OPP_REAL>("grid_spacing");

        std::shared_ptr<DataPointers> m = LoadData();

        opp_set cell_set      = opp_decl_set(m->n_cells, "mesh_cells");
        opp_map cell_cell_map = opp_decl_map(cell_set, cell_set, NEIGHBOURS, m->cell_cell_map, "c_c_map");
        opp_dat cell_index    = opp_decl_dat(cell_set, ONE, DT_INT,  m->c_index,  "c_index");
        opp_dat cell_pos_ll   = opp_decl_dat(cell_set, DIM, DT_REAL, m->c_pos_ll, "c_pos_ll");       
        opp_dat cell_colors   = opp_decl_dat(cell_set, ONE, DT_INT,  m->c_colors, "c_colors"); // used only with MPI

        opp_set part_set      = opp_decl_particle_set("particles", cell_set); // Zero particles, inject after partitioning
        opp_dat part_pos      = opp_decl_dat(part_set, DIM, DT_REAL, nullptr, "p_pos");
        opp_dat part_vel      = opp_decl_dat(part_set, DIM, DT_REAL, nullptr, "p_vel");    
        opp_dat part_mesh_rel = opp_decl_dat(part_set, ONE, DT_INT,  nullptr, "p_mesh_rel", true);
        opp_dat part_index    = opp_decl_dat(part_set, ONE, DT_INT,  nullptr, "p_index"); // Unused in the simulation

        opp_decl_const<OPP_REAL>(TWO, extents,     "CONST_extents");
        opp_decl_const<OPP_REAL>(ONE, &dt,         "CONST_dt");
        opp_decl_const<OPP_REAL>(ONE, &cell_width, "CONST_cell_width");
        opp_decl_const<OPP_INT>(TWO,  ndimcells,   "CONST_ndimcells");

        m->DeleteValues();

        // ideally opp_colour_cartesian_mesh is not required for non-mpi runs
        std::vector<int> counts = { opp_params->get<OPP_INT>("nx"), opp_params->get<OPP_INT>("ny") };
        opp_colour_cartesian_mesh(DIM, counts, cell_index, cell_colors);

#ifdef USE_MPI
        opp_partition(std::string("EXTERNAL"), cell_set, nullptr, cell_colors);
#endif
        
        init_particles(part_index, part_pos, part_vel, part_mesh_rel, cell_pos_ll, cell_index);

#ifndef FUSE_KERNELS 
        initialize_particle_mover(grid_spacing, DIM, cell_pos_ll, cell_index, cell_cell_map);
#endif

        opp_printf("Setup", "Cells[%d] Particles[%d] max_iter[%d]", cell_set->size, part_set->size, max_iter);

        opp_profiler->end("Setup");

        opp_profiler->start("MainLoop");
        for (OPP_main_loop_iter = 0; OPP_main_loop_iter < max_iter; OPP_main_loop_iter++) // Start Main loop
        {

#ifdef FUSE_KERNELS        
            opp_particle_move__update_pos_move(
                part_set,
                opp_arg_dat(part_mesh_rel, OPP_RW),
                opp_arg_dat(part_vel,      OPP_RW),
                opp_arg_dat(part_pos,      OPP_RW),
                opp_arg_dat(cell_pos_ll,   OPP_READ, OPP_Map_from_Mesh_Rel),
                opp_arg_dat(cell_cell_map, OPP_READ, OPP_Map_from_Mesh_Rel)
            );
#else
            opp_loop_all__update_pos(
                part_set,
                opp_arg_dat(part_vel,      OPP_READ),
                opp_arg_dat(part_pos,      OPP_RW)
            );

            opp_particle_move__move(
                part_set,
                opp_arg_dat(part_mesh_rel, OPP_RW),
                opp_arg_dat(part_pos,      OPP_READ),
                opp_arg_dat(cell_pos_ll,   OPP_READ, OPP_Map_from_Mesh_Rel),
                opp_arg_dat(cell_cell_map, OPP_READ, OPP_Map_from_Mesh_Rel)
            );
#endif
            std::string log = "";
            if (verify_parts)
            {
                int incorrect_part_count = 0;
                opp_loop_all__verify(
                    part_set,
                    opp_arg_dat(part_mesh_rel, OPP_READ),
                    opp_arg_dat(part_pos,      OPP_READ),
                    opp_arg_dat(cell_index,    OPP_READ, OPP_Map_from_Mesh_Rel),
                    opp_arg_gbl(&incorrect_part_count, 1, "int", OPP_INC)
                );
                log += str(incorrect_part_count, "%d Errors");
            }

            if (OPP_rank == OPP_ROOT) 
                opp_printf("Main", "ts: %d | %s ****", OPP_main_loop_iter, log.c_str());
        
        } // End Main loop
        opp_profiler->end("MainLoop");
        
        int64_t total_glb_particles = (int64_t)(part_set->size), local_particles = (int64_t)(part_set->size);
#ifdef USE_MPI
        MPI_Reduce(&local_particles, &total_glb_particles, 1, MPI_INT64_T, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
#endif
        if (OPP_rank == OPP_ROOT) 
            opp_printf("Main", "Main loop completed after %d iterations {particles=%d} ****", 
                max_iter, total_glb_particles);
    }

    opp_exit();

    return 0;
}

/*
NOTE: ------------------------------------------------------------------
    cell_set     	        28bytes
    part_set     	        40bytes

    Assume, 1024 x 1024 mesh 
    cell_set->size          1,048,576     29,360,128bytes
    
    Assume 16GB of GPU memory and we have two thrust vectors per opp_dat
    max particles           ~374,000,000
    max particles per cell  ~175
------------------------------------------------------------------------
*/


// opp_print_dat_to_txtfile(part_mesh_rel, "FINAL", "part_mesh_rel.dat");
// opp_print_dat_to_txtfile(part_pos, "FINAL", "part_pos.dat");
// opp_print_map_to_txtfile(cell_cell_map, "INIT_A", "cell_cell_map.dat");
// opp_print_dat_to_txtfile(part_mesh_rel, "INIT_A", "part_mesh_rel.dat");
// opp_print_dat_to_txtfile(part_pos, "INIT_A", "part_pos.dat");

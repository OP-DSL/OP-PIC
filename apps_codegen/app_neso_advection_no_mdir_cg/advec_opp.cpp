
// Auto-generated at 2025-02-03 16:37:08.638646 by opp-translator
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

#include "opp_lib.h"

void opp_par_loop_all__update_pos_kernel(opp_set,opp_arg,opp_arg);
void opp_particle_move__move_kernel(opp_set,opp_map,opp_map,opp_arg,opp_arg);
void opp_par_loop_all__verify_kernel(opp_set,opp_arg,opp_arg,opp_arg);
void opp_init_direct_hop_cg(double,const opp_dat,const opp::BoundingBox&,opp_map,opp_map,opp_arg,opp_arg);

#include "advec_misc.h"
// #include "kenels.h" // codegen commented...

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    opp_init(argc, argv);

    {
    opp_profiler->start("Setup");

        OPP_INT max_iter            = opp_params->get<OPP_INT>("num_steps");   
        OPP_REAL dt                 = opp_params->get<OPP_REAL>("dt");
        OPP_REAL cell_width         = opp_params->get<OPP_REAL>("cell_width");
        OPP_REAL extents[2]         = {opp_params->get<OPP_INT>("nx")*cell_width, opp_params->get<OPP_INT>("ny")*cell_width};
        OPP_INT ndimcells[2]        = {opp_params->get<OPP_INT>("nx"), opp_params->get<OPP_INT>("ny")};
        OPP_BOOL verify_parts       = opp_params->get<OPP_BOOL>("verify_particles");
        OPP_REAL grid_spacing       = opp_params->get<OPP_REAL>("grid_spacing");
        const opp_point expansion(cell_width*2, cell_width*2, 0.0);
        int64_t total_part_iter     = 0;
        int incorrect_part_count    = 0;

        std::shared_ptr<DataPointers> m = LoadData();

        opp_set cell_set = opp_decl_set(m->n_cells, "mesh_cells");
        opp_map c2c_map  = opp_decl_map(cell_set, cell_set, NEIGHBOURS, m->cell_cell_map, "c_c_map");
        
        opp_dat c_idx    = opp_decl_dat(cell_set, ONE, DT_INT,  m->c_index,  "c_index");
        opp_dat c_pos_ll = opp_decl_dat(cell_set, DIM, DT_REAL, m->c_pos_ll, "c_pos_ll");       
        opp_dat c_colors = opp_decl_dat(cell_set, ONE, DT_INT,  m->c_colors, "c_colors"); // used only with MPI

        opp_set part_set = opp_decl_particle_set("particles", cell_set); // Zero particles, inject after partitioning
        opp_map p2c_map  = opp_decl_map(part_set, cell_set, 1, nullptr, "p_mesh_rel");

        opp_dat p_pos    = opp_decl_dat(part_set, DIM, DT_REAL, nullptr, "p_pos");
        opp_dat p_vel    = opp_decl_dat(part_set, DIM, DT_REAL, nullptr, "p_vel");
        opp_dat p_mdir   = opp_decl_dat(part_set, DIM, DT_INT,  nullptr, "p_move_dir");
        opp_dat p_idx    = opp_decl_dat(part_set, ONE, DT_INT,  nullptr, "p_index"); // Unused in the simulation

        opp_decl_const<OPP_REAL>(TWO, extents,     "CONST_extents");
        opp_decl_const<OPP_REAL>(ONE, &dt,         "CONST_dt");
        opp_decl_const<OPP_REAL>(ONE, &cell_width, "CONST_cell_width");
        opp_decl_const<OPP_INT>(TWO,  ndimcells,   "CONST_ndimcells");

        m->DeleteValues();

#ifdef USE_MPI
        // std::vector<int> counts = { opp_params->get<OPP_INT>("nx"), opp_params->get<OPP_INT>("ny") };
        // opp_colour_cartesian_mesh(DIM, counts, c_idx, c_colors);
        advec_color_block(c_idx, c_colors);
        
        opp_partition(std::string("EXTERNAL"), cell_set, nullptr, c_colors);
#endif
        
        init_particles(p_idx, p_pos, p_vel, p_mdir, p2c_map, c_pos_ll, c_idx);
        
        // these two lines are only required if we plan to use direct_hop
        opp::BoundingBox bounding_box(c_pos_ll, DIM, expansion);
        opp_init_direct_hop_cg(grid_spacing, c_idx, bounding_box, c2c_map, p2c_map,
			opp_arg_dat(p_pos, OPP_READ),
			opp_arg_dat(c_pos_ll, p2c_map, OPP_READ));

        opp_printf("Setup Completed", "Cells[%d] Particles[%d] max_iter[%d]", cell_set->size, part_set->size, max_iter);

    opp_profiler->end("Setup");

    opp_profiler->start("MainLoop");
        for (OPP_main_loop_iter = 0; OPP_main_loop_iter < max_iter; OPP_main_loop_iter++)
        {
            opp_par_loop_all__update_pos_kernel(part_set, 
                opp_arg_dat(p_vel,  OPP_READ),
                opp_arg_dat(p_pos,  OPP_RW)
            );

            opp_particle_move__move_kernel(part_set, c2c_map, p2c_map,
                opp_arg_dat(p_pos,             OPP_READ),
                opp_arg_dat(c_pos_ll, p2c_map, OPP_READ)
            );

            if (verify_parts)
            {
                incorrect_part_count = 0;
                opp_par_loop_all__verify_kernel(part_set, 
                    opp_arg_dat(p_pos,          OPP_READ),
                    opp_arg_dat(c_idx, p2c_map, OPP_READ),
                    opp_arg_gbl(&incorrect_part_count, 1, "int", OPP_INC)
                );
            }

            std::string log = get_global_level_log(part_set, incorrect_part_count);
            OPP_RUN_ON_ROOT()
                opp_printf("Main", "ts: %d | %s", OPP_main_loop_iter, log.c_str());  
        
            total_part_iter += part_set->size; 
        }
    opp_profiler->end("MainLoop");
        
        OPP_RUN_ON_ROOT()
            opp_printf("Main", "Completed %d iterations [%" PRId64 " parts iterated in ROOT]", 
                max_iter, total_part_iter);
    }

    opp_exit();

    return 0;
}

/*
NOTE: ------------------------------------------------------------------
    cell_set     	        40 bytes
    part_set     	        48 bytes + 12 bytes (3 ints used in particle move)

    Assume, 512 x 512 mesh 
    cell_set->size          262,144     10,485,760 bytes
    
    Assume 16GB of GPU memory and we have two thrust vectors per opp_dat
    max particles           ~125,000,000
    max particles per cell  ~475

    Assume 64GB of GPU memory and we have two thrust vectors per opp_dat
    max particles           ~525,000,000
    max particles per cell  ~2,002

    Assume, 1024 x 1024 mesh 
    cell_set->size          1,048,576     41,943,040 bytes

    Assume 64GB of GPU memory and we have two thrust vectors per opp_dat
    max particles           ~646,000,000
    max particles per cell  ~500

    or reduce the mesh size to increase particles per cell
------------------------------------------------------------------------
ESTIMATED --
mesh      -> cells      | ppc  -> parts       -> bytes  | ppc  -> parts       -> bytes  | ppc  -> parts       -> bytes  |
256x256   -> 65,536	    | 6800 -> 445,644,800 -> 49.8GB | 3400 -> 222,822,400 -> 24.9GB | 1700 -> 111,411,200 -> 12.5GB |
512x256   -> 131,072	| 3400 -> 445,644,800 -> 49.8GB | 1700 -> 222,822,400 -> 24.9GB | 850  -> 111,411,200 -> 12.5GB |
512x512   -> 262,144	| 1700 -> 445,644,800 -> 49.8GB | 850  -> 222,822,400 -> 24.9GB | 425  -> 111,411,200 -> 12.5GB |
1024x1024 -> 1,048,576	| 425  -> 445,644,800 -> 49.8GB | 212  -> 222,298,112 -> 24.8GB |                               |

ACTUAL --
256x256-1700 - 111,411,200 - 14.6GB GPU - H100 5.7s
256x256-6800 - 445,644,800 - 56.6GB GPU - H100 22.7s
*/
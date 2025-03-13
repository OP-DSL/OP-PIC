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

#include "opp_hdf5.h"
#include "cabana_misc.h"

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
        opp_map c2c_map     = opp_decl_mesh_map(c_set, c_set, NEIGHBOURS, m->c2c_map,     "c2c_map");
        opp_map c2ngc_map   = opp_decl_mesh_map(c_set, c_set, FACES,      m->c2ngc_map,   "c2c_non_ghost_map");
        opp_map c2cug0_map  = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cug0_map,  "c2c_ug_0_map");
        opp_map c2cug1_map  = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cug1_map,  "c2c_ug_1_map");
        opp_map c2cug2_map  = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cug2_map,  "c2c_ug_2_map");
        opp_map c2cug3_map  = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cug3_map,  "c2c_ug_3_map");
        opp_map c2cug4_map  = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cug4_map,  "c2c_ug_4_map");
        opp_map c2cug5_map  = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cug5_map,  "c2c_ug_5_map");
        opp_map c2cugb0_map = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cugb0_map, "c2c_ugb_0_map");
        opp_map c2cugb1_map = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cugb1_map, "c2c_ugb_1_map");
        opp_map c2cugb2_map = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cugb2_map, "c2c_ugb_2_map");
        opp_map c2cugb3_map = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cugb3_map, "c2c_ugb_3_map");
        opp_map c2cugb4_map = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cugb4_map, "c2c_ugb_4_map");
        opp_map c2cugb5_map = opp_decl_mesh_map(c_set, c_set, 1,          m->c2cugb5_map, "c2c_ugb_5_map");

        opp_dat c_index      = opp_decl_dat(c_set, ONE,        DT_INT,  m->c_index,      "c_index");
        opp_dat c_pos_ll     = opp_decl_dat(c_set, DIM,        DT_REAL, m->c_pos_ll,     "c_pos_ll");       
        opp_dat c_mask_ghost = opp_decl_dat(c_set, ONE,        DT_INT,  m->c_ghost,      "c_mask_ghost");
        opp_dat c_mask_right = opp_decl_dat(c_set, ONE,        DT_INT,  m->c_mask_right, "c_mask_right");
        opp_dat c_mask_ug    = opp_decl_dat(c_set, 6,          DT_INT,  m->c_mask_ug,    "c_mask_ug");
        opp_dat c_mask_ugb   = opp_decl_dat(c_set, 6,          DT_INT,  m->c_mask_ugb,   "c_mask_ugb");

        m->DeleteValues();

        MPI_Barrier(MPI_COMM_WORLD);
        opp_printf("CABANA_CONVERT_MESH", "Init opp structures DONE, dumping to HDF5");

        std::string file_out = opp_params->get<OPP_STRING>("hdf_filename");
        opp_dump_to_hdf5(file_out.c_str());

        opp_printf("CABANA_CONVERT_MESH", "HDF5 dumping done");

        // Expect to crash since we do not partition, and opp_exit tries to delete halos etc...
    }

    opp_exit();

    return 0;
}

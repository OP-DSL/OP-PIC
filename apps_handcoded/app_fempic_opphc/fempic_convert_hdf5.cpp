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

#include "fempic_loader.h"
#include "fempic_misc.cpp"
#include "opp_hdf5.h"

using namespace opp;

//*********************************************MAIN****************************************************
int main(int argc, char **argv) 
{
    if (argc < 2) 
    {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << 
        "(Eg. /ext-home/zl/phd/OP-PIC/scripts/fempic_tests/configs/coarse_1.param)" << std::endl;
        exit(-1);
    }

    opp_init(argc, argv);
    opp_params->write(std::cout);

    if (OPP_comm_size > 1) 
        opp_printf("Main", "Warnining: Run with 1 MPI rank to avoid intermittent issues...");

    {
        opp_profiler->start("Setup");

        opp_profiler->start("LoadMesh");

        MPI_Barrier(MPI_COMM_WORLD);
        if (OPP_rank == OPP_ROOT) opp_printf("Main", "Load mesh Start");

        std::shared_ptr<FieldPointers> g_m, m; // g_m - global mesh, m - local mesh
        g_m = std::make_shared<FieldPointers>();
        
        if (OPP_rank == OPP_ROOT) // Load using the original FemPIC loaders and distribute
        {  
            g_m = LoadMesh(opp_params->get<OPP_STRING>("global_mesh"), opp_params->get<OPP_STRING>("inlet_mesh"),
                                opp_params->get<OPP_STRING>("wall_mesh"), opp_params->get<OPP_STRING>("cluster"));
            
            opp_printf("FEMPIC_CONVERT_MESH", "Global counts - Nodes[%d] Cells[%d] IFaces[%d]", 
                g_m->n_nodes, g_m->n_cells, g_m->n_ifaces);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        if (OPP_rank == OPP_ROOT) opp_printf("Main", "DistributeMeshOverRanks Start");

        DistributeMeshOverRanks(g_m, m);

        MPI_Barrier(MPI_COMM_WORLD);
        if (OPP_rank == OPP_ROOT) opp_printf("Main", "Creating Sets");

        opp_profiler->end("LoadMesh");

        opp_set node_set         = opp_decl_mesh_set(m->n_nodes, "mesh_nodes");
        opp_set cell_set         = opp_decl_mesh_set(m->n_cells, "mesh_cells");
        opp_set iface_set        = opp_decl_mesh_set(m->n_ifaces, "inlet_faces_cells");
        opp_set particle_set     = opp_decl_part_set("particles", cell_set); 

        MPI_Barrier(MPI_COMM_WORLD);
        if (OPP_rank == OPP_ROOT) opp_printf("Main", "Creating Maps");

        opp_map cell_v_nodes_map = opp_decl_mesh_map(cell_set,  node_set, N_PER_C,  m->c_to_n, "c_v_n_map");
        opp_map cell_v_cell_map  = opp_decl_mesh_map(cell_set,  cell_set, NEIGHB_C, m->c_to_c,  "c_v_c_map"); 
        opp_map iface_v_cell_map = opp_decl_mesh_map(iface_set, cell_set, ONE,      m->if_to_c, "if_v_c_map"); 
        opp_map iface_v_node_map = opp_decl_mesh_map(iface_set, node_set, N_PER_IF, m->if_to_n, "if_v_n_map");

        MPI_Barrier(MPI_COMM_WORLD);
        if (OPP_rank == OPP_ROOT) opp_printf("Main", "Creating Mesh Dats");

        opp_dat cell_det         = opp_decl_dat(cell_set, ALL_DET,     DT_REAL, m->c_det,      "c_det");  
        opp_dat cell_volume      = opp_decl_dat(cell_set, ONE,         DT_REAL, m->c_vol,      "c_volume");        
        opp_dat cell_ef          = opp_decl_dat(cell_set, DIM,         DT_REAL, m->c_ef,       "c_ef");
        opp_dat cell_shape_deriv = opp_decl_dat(cell_set, N_PER_C*DIM, DT_REAL, m->c_sd,       "c_shape_deri"); 
        opp_dat global_cell_id   = opp_decl_dat(cell_set, ONE,         DT_INT,  m->c_id,       "c_gbl_id"); 
        opp_dat cell_colors      = opp_decl_dat(cell_set, ONE,         DT_INT,  m->c_col,      "c_colors");
        opp_dat cell_centroids   = opp_decl_dat(cell_set, DIM,         DT_REAL, m->c_centroid, "c_centroids");

        opp_dat node_volume      = opp_decl_dat(node_set, ONE, DT_REAL, m->n_vol,     "n_vol");        
        opp_dat node_potential   = opp_decl_dat(node_set, ONE, DT_REAL, m->n_pot,     "n_potential");     
        opp_dat node_charge_den  = opp_decl_dat(node_set, ONE, DT_REAL, m->n_ion_den, "n_charge_den");
        opp_dat node_pos         = opp_decl_dat(node_set, DIM, DT_REAL, m->n_pos,     "n_pos");     
        opp_dat node_type        = opp_decl_dat(node_set, ONE, DT_INT,  m->n_type,    "n_type");
        opp_dat node_bnd_pot     = opp_decl_dat(node_set, ONE, DT_REAL, m->n_bnd_pot, "n_bnd_pot");
        opp_dat node_id          = opp_decl_dat(node_set, ONE, DT_INT,  m->n_id,      "n_id"); 
        opp_dat node_colors      = opp_decl_dat(node_set, ONE, DT_INT,  m->n_color,   "n_colors");

        opp_dat iface_v_norm     = opp_decl_dat(iface_set, DIM,          DT_REAL, m->if_v_norm, "iface_v_norm");        
        opp_dat iface_u_norm     = opp_decl_dat(iface_set, DIM,          DT_REAL, m->if_u_norm, "iface_u_norm"); 
        opp_dat iface_norm       = opp_decl_dat(iface_set, DIM,          DT_REAL, m->if_norm,   "iface_norm");     
        opp_dat iface_area       = opp_decl_dat(iface_set, ONE,          DT_REAL, m->if_area,   "iface_area");
        opp_dat iface_dist       = opp_decl_dat(iface_set, ONE,          DT_INT,  m->if_dist,   "iface_dist");
        opp_dat iface_n_pos      = opp_decl_dat(iface_set, N_PER_IF*DIM, DT_REAL, m->if_n_pos,  "iface_n_pos"); 
        opp_dat iface_id         = opp_decl_dat(iface_set, ONE,          DT_INT,  m->if_id,     "iface_id"); 

        MPI_Barrier(MPI_COMM_WORLD);
        if (OPP_rank == OPP_ROOT) opp_printf("Main", "Creating Particle Dats");

        opp_dat part_position    = opp_decl_dat(particle_set, DIM,     DT_REAL, nullptr, "part_position");
        opp_dat part_velocity    = opp_decl_dat(particle_set, DIM,     DT_REAL, nullptr, "part_velocity");    
        opp_dat part_lc          = opp_decl_dat(particle_set, N_PER_C, DT_REAL, nullptr, "part_lc");
        opp_dat part_mesh_rel    = opp_decl_dat(particle_set, ONE,     DT_INT,  nullptr, "part_mesh_rel", true);
        opp_dat part_id          = opp_decl_dat(particle_set, ONE,     DT_INT,  nullptr, "part_id");

        opp_set dummy_part_set   = opp_decl_part_set("dummy particles", cell_set); 
        opp_dat dummy_part_rand  = opp_decl_dat(dummy_part_set, 2, DT_REAL, nullptr, "dummy_part_rand");

        MPI_Barrier(MPI_COMM_WORLD);
        if (OPP_rank == OPP_ROOT) opp_printf("Main", "Deleting Loader");

        m->DeleteValues();

        opp_profiler->end("Setup");

        MPI_Barrier(MPI_COMM_WORLD);
        opp_printf("FEMPIC_CONVERT_MESH", "Init opp structures DONE, dumping to HDF5");

        std::string file_out = opp_params->get<OPP_STRING>("hdf_filename");
        opp_dump_to_hdf5(file_out.c_str());

        opp_printf("FEMPIC_CONVERT_MESH", "HDF5 dumping done");

#ifdef USE_MPI
        // since we are using mpi, app will crash if we do not partition, since it tries to delete halos etc...
        genColoursForBlockPartition(cell_colors, cell_centroids, iface_n_pos, iface_v_node_map);
        opp_partition(std::string("EXTERNAL"), cell_set, nullptr, cell_colors);
#endif
    }

    opp_exit();

    return 0;
}

//*********************************************MAIN****************************************************
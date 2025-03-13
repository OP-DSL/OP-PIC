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

#pragma once

//*********************************************
// USER WRITTEN CODE
//*********************************************

#include "fempic_defs.h"
#include "minifempic_funcs.h"

void init_mesh(std::shared_ptr<DataPointers> mesh);
void distribute_data_over_ranks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m);

/***************************************************************************************************
 * @brief Initialize the rank specific mesh data to a DataPointers utility class shared pointer
 * @return std::shared_ptr<DataPointers>
 */
inline std::shared_ptr<DataPointers> load_mesh()
{
    std::shared_ptr<DataPointers> g_m(new DataPointers());

    OPP_RUN_ON_ROOT() init_mesh(g_m);

    std::shared_ptr<DataPointers> m = std::make_shared<DataPointers>();;
    distribute_data_over_ranks(g_m, m);

    return m;
}

/***************************************************************************************************
 * @brief Initializes the mesh using 2D (nx,ny) and cell_width values in the config file
 *          Expect this to run only on the ROOT MPI rank
 * @param m std::shared_ptr<DataPointers> loaded with mesh data
 * @return (void)
 */
void init_mesh(std::shared_ptr<DataPointers> mesh) { OPP_RETURN_IF_INVALID_PROCESS;

    std::shared_ptr<Volume> volume(new Volume());
    if (!LoadVolumeMesh(opp_params->get<OPP_STRING>("global_mesh"), *(volume.get())) ||
        !(LoadSurfaceMesh(opp_params->get<OPP_STRING>("inlet_mesh"), *(volume.get()),INLET, false)) ||
        !(LoadSurfaceMesh(opp_params->get<OPP_STRING>("wall_mesh"), *(volume.get()), FIXED, false))) 
    {    
        return;
    }

    volume->summarize(std::cout);

    mesh->n_nodes  = volume->nodes.size();
    mesh->n_cells  = volume->elements.size();
    mesh->n_ifaces = volume->inlet_faces.size();

    opp_printf("load_mesh", "nodes %d cells %d ifaces %d", mesh->n_nodes, mesh->n_cells, mesh->n_ifaces);

    mesh->CreateMeshArrays();

    for (int n=0; n<mesh->n_nodes; n++)
    {
        mesh->n_pot[n] = 0.0f;
        mesh->n_ion_den[n] = 0.0f;

        Node &node = volume->nodes[n];
    
        for (int dim=0; dim<DIM; dim++)
            mesh->n_pos[n * DIM + dim] = node.pos[dim];
    
        mesh->n_vol[n] = node.volume;
        mesh->n_type[n]   = (int)node.type;
        mesh->n_id[n]     = n;

        if (n % 10000 == 0)
            opp_printf("load_mesh", "Mesh Nodes at %d", n);  
    }

    opp_printf("load_mesh", "Mesh Nodes Loaded");

    for (int cID=0; cID<mesh->n_cells; cID++)
    {
        Tetra &tet = volume->elements[cID];
        
        for (int nodeCon=0; nodeCon<N_PER_C; nodeCon++)
        {
            mesh->c_to_n[cID * N_PER_C + nodeCon] = tet.con[nodeCon];

            mesh->c_sd[cID * (N_PER_C*DIM) + nodeCon * DIM + 0 ] = 0.0;
            mesh->c_sd[cID * (N_PER_C*DIM) + nodeCon * DIM + 1 ] = 0.0;
            mesh->c_sd[cID * (N_PER_C*DIM) + nodeCon * DIM + 2 ] = 0.0;
        }

        for (int cConn=0; cConn<NEIGHB_C; cConn++)
        {
            mesh->c_to_c[cID * NEIGHB_C + cConn]     = tet.cell_con[cConn];

            mesh->c_det[(cID * NEIGHB_C + cConn) * DET_FIELDS + 0] = (tet.alpha[cConn]);
            mesh->c_det[(cID * NEIGHB_C + cConn) * DET_FIELDS + 1] = (tet.beta[cConn]);
            mesh->c_det[(cID * NEIGHB_C + cConn) * DET_FIELDS + 2] = (tet.gamma[cConn]);
            mesh->c_det[(cID * NEIGHB_C + cConn) * DET_FIELDS + 3] = (tet.delta[cConn]);
        }

        const opp::Point3D *node0_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[0]].pos);
        const opp::Point3D *node1_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[1]].pos);
        const opp::Point3D *node2_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[2]].pos);
        const opp::Point3D *node3_pos = reinterpret_cast<opp::Point3D*>(volume->nodes[tet.con[3]].pos);
        const opp::Point3D c_pos = opp::getTetraCentroid3D(*node0_pos, *node1_pos, *node2_pos, *node3_pos);

        mesh->c_centroid[cID * DIM + 0] = c_pos.x;
        mesh->c_centroid[cID * DIM + 1] = c_pos.y;
        mesh->c_centroid[cID * DIM + 2] = c_pos.z;

        mesh->c_vol[cID] = tet.volume;
        mesh->c_id[cID] = cID;

        mesh->c_ef[cID * DIM + 0] = 0.0;
        mesh->c_ef[cID * DIM + 1] = 0.0;
        mesh->c_ef[cID * DIM + 2] = 0.0;        

        if (cID % 10000 == 0)
            opp_printf("load_mesh", "Mesh Cells at %d", cID);   
    }

    opp_printf("load_mesh", "Mesh Cells Loaded");

    for (int faceID=0; faceID<mesh->n_ifaces; faceID++)
    {
        Face &face = volume->inlet_faces[faceID];
        
        mesh->if_to_c[faceID] = face.cell_con;
        mesh->if_area[faceID] = face.area;
        mesh->if_dist[faceID] = 0;
        mesh->if_id[faceID] = faceID;

        for (int dim=0; dim<3; dim++)
        {
            mesh->if_v_norm[faceID * 3 + dim] = face.v[dim];   
            mesh->if_u_norm[faceID * 3 + dim] = face.u[dim];   
            mesh->if_norm[faceID * 3 + dim]   = face.normal[dim]; 
        }

        for (int n=0; n< N_PER_IF; n++)
        {
            mesh->if_to_n[faceID * 3 + n] = face.con[n]; 

            Node &node = volume->nodes[face.con[n]];
            for (int dim=0; dim<3; dim++)
            {
                mesh->if_n_pos[faceID * (N_PER_IF * DIM) + DIM * n + dim] = node.pos[dim];
            }
        }
    }

    opp_printf("Main", "Global Nodes[%d] Cells[%d] IFaces[%d]", mesh->n_nodes, mesh->n_cells, mesh->n_ifaces);
}

//*************************************************************************************************
inline void distribute_data_over_ranks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m)
{ OPP_RETURN_IF_INVALID_PROCESS;
#ifdef USE_MPI
    MPI_Bcast(&(g_m->n_nodes), 1, MPI_INT, OPP_ROOT, OPP_MPI_WORLD);
    MPI_Bcast(&(g_m->n_cells), 1, MPI_INT, OPP_ROOT, OPP_MPI_WORLD);
    MPI_Bcast(&(g_m->n_ifaces), 1, MPI_INT, OPP_ROOT, OPP_MPI_WORLD);

    m->n_nodes  = opp_get_uniform_local_size(g_m->n_nodes);
    m->n_cells  = opp_get_uniform_local_size(g_m->n_cells);
    m->n_ifaces = opp_get_uniform_local_size(g_m->n_ifaces);
    
    m->CreateMeshArrays();

    opp_uniform_scatter_array(g_m->c_ef       , m->c_ef       , g_m->n_cells , m->n_cells , DIM);
    opp_uniform_scatter_array(g_m->c_to_n     , m->c_to_n     , g_m->n_cells , m->n_cells , N_PER_C); 
    opp_uniform_scatter_array(g_m->c_to_c     , m->c_to_c     , g_m->n_cells , m->n_cells , NEIGHB_C); 
    opp_uniform_scatter_array(g_m->c_det      , m->c_det      , g_m->n_cells , m->n_cells , ALL_DET); 
    opp_uniform_scatter_array(g_m->c_vol      , m->c_vol      , g_m->n_cells , m->n_cells , ONE); 
    opp_uniform_scatter_array(g_m->c_sd       , m->c_sd       , g_m->n_cells , m->n_cells , N_PER_C*DIM); 
    opp_uniform_scatter_array(g_m->c_col      , m->c_col      , g_m->n_cells , m->n_cells , ONE);
    opp_uniform_scatter_array(g_m->c_id       , m->c_id       , g_m->n_cells , m->n_cells , ONE); 
    opp_uniform_scatter_array(g_m->c_centroid , m->c_centroid , g_m->n_cells , m->n_cells , DIM); 
    opp_uniform_scatter_array(g_m->n_bnd_pot  , m->n_bnd_pot  , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->n_pot      , m->n_pot      , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->n_ion_den  , m->n_ion_den  , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->n_pos      , m->n_pos      , g_m->n_nodes , m->n_nodes , DIM); 
    opp_uniform_scatter_array(g_m->n_vol      , m->n_vol      , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_type     , m->n_type     , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_color    , m->n_color    , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_bnd_pot  , m->n_bnd_pot  , g_m->n_nodes , m->n_nodes , ONE);
    opp_uniform_scatter_array(g_m->n_id       , m->n_id       , g_m->n_nodes , m->n_nodes , ONE); 
    opp_uniform_scatter_array(g_m->if_to_c    , m->if_to_c    , g_m->n_ifaces, m->n_ifaces, ONE); 
    opp_uniform_scatter_array(g_m->if_to_n    , m->if_to_n    , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_v_norm  , m->if_v_norm  , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_u_norm  , m->if_u_norm  , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_norm    , m->if_norm    , g_m->n_ifaces, m->n_ifaces, DIM); 
    opp_uniform_scatter_array(g_m->if_area    , m->if_area    , g_m->n_ifaces, m->n_ifaces, ONE); 
    opp_uniform_scatter_array(g_m->if_dist    , m->if_dist    , g_m->n_ifaces, m->n_ifaces, ONE); 
    opp_uniform_scatter_array(g_m->if_n_pos   , m->if_n_pos   , g_m->n_ifaces, m->n_ifaces, N_PER_IF*DIM);       
    opp_uniform_scatter_array(g_m->if_id      , m->if_id      , g_m->n_ifaces, m->n_ifaces, ONE); 

    if (OPP_rank == OPP_ROOT)
        g_m->DeleteValues();
#else
    m = g_m;
#endif
}


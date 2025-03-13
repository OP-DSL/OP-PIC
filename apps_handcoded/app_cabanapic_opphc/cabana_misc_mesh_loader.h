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

#include "cabana_defs.h"

void init_mesh(const Deck& deck, std::shared_ptr<DataPointers> m);
void distribute_data_over_ranks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m);

/***************************************************************************************************
 * @brief Initialize the rank specific mesh data to a DataPointers utility class shared pointer
 * @return std::shared_ptr<DataPointers>
 */
std::shared_ptr<DataPointers> load_mesh(const Deck& deck) {

    std::shared_ptr<DataPointers> g_m(new DataPointers());

    OPP_RUN_ON_ROOT() init_mesh(deck, g_m);

    std::shared_ptr<DataPointers> m;
    distribute_data_over_ranks(g_m, m);

    return m;
}

/***************************************************************************************************
 * @brief Initializes the mesh using 2D (nx,ny) and cell_width values in the config file
 *          Expect this to run only on the ROOT MPI rank
 * @param m std::shared_ptr<DataPointers> loaded with mesh data
 * @return (void)
 */
void init_mesh(const Deck& deck, std::shared_ptr<DataPointers> m) {

    const OPP_INT nx             = deck.nx;
    const OPP_INT ny             = deck.ny;
    const OPP_INT nz             = deck.nz;
    const OPP_REAL c_widths[DIM] = { deck.dx, deck.dy, deck.dz };

    m->n_cells   = (nx + 2*NG) * (ny + 2*NG) * (nz + 2*NG);

    opp_printf("Setup", "init_mesh global n_cells=%d nx=%d ny=%d nz=%d", m->n_cells, nx, ny, nz);

    m->CreateMeshCommArrays();

    for (int n = 0; n < m->n_cells; n++) {
        
        for (int i = 0; i < 6; i++)
            m->c_mask_ugb[n * 6 + i] = 0;

        for (int i = 0; i < 6; i++)
            m->c_mask_ug[n * 6 + i] = 0;

        for (int i = 0; i < NEIGHBOURS; i++)
            m->c2c_map[n * NEIGHBOURS + i] = -1;

        for (int i = 0; i < FACES; i++)
            m->c2ngc_map[n * FACES + i] = n;
        
        for (int i = 0; i < 6; i++)
            m->c2cug_map[n * 6 + i] = n;
        
        for (int i = 0; i < 6; i++)
            m->c2cugb_map[n * 6 + i] = n;

        m->c2cug0_map[n] = n;
        m->c2cug1_map[n] = n;
        m->c2cug2_map[n] = n;
        m->c2cug3_map[n] = n;
        m->c2cug4_map[n] = n;
        m->c2cug5_map[n] = n;

        m->c2cugb0_map[n] = n;
        m->c2cugb1_map[n] = n;
        m->c2cugb2_map[n] = n;
        m->c2cugb3_map[n] = n;
        m->c2cugb4_map[n] = n;
        m->c2cugb5_map[n] = n;

        m->c_index[n]      = n;
		m->c_ghost[n]      = 1;
		m->c_mask_right[n] = 1;
    }

	for (int x = 0; x < nx+2*NG; x++) {
		for (int y = 0; y < ny+2*NG; y++) {
			for (int z = 0; z < nz+2*NG; z++) {

                const int i = VOXEL(x,y,z, nx,ny,nz);

                m->c_pos_ll[i*DIM + Dim::x] = x * c_widths[Dim::x];
                m->c_pos_ll[i*DIM + Dim::y] = y * c_widths[Dim::y];
                m->c_pos_ll[i*DIM + Dim::z] = z * c_widths[Dim::z];

                VOXEL_MAP(x-1, y-1, z  , nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::xd_yd_z]); 
                VOXEL_MAP(x-1, y  , z-1, nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::xd_y_zd]); 
                VOXEL_MAP(x-1, y  , z  , nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::xd_y_z]);  
                VOXEL_MAP(x  , y-1, z-1, nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::x_yd_zd]); 
                VOXEL_MAP(x  , y-1, z  , nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::x_yd_z]);  
                VOXEL_MAP(x  , y  , z-1, nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::x_y_zd]);  
                VOXEL_MAP(x  , y  , z+1, nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::x_y_zu]);  
                VOXEL_MAP(x  , y+1, z  , nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::x_yu_z]);  
                VOXEL_MAP(x  , y+1, z+1, nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::x_yu_zu]); 
                VOXEL_MAP(x+1, y  , z  , nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::xu_y_z]);  
                VOXEL_MAP(x+1, y  , z+1, nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::xu_y_zu]); 
                VOXEL_MAP(x+1, y+1, z  , nx, ny, nz, m->c2c_map[i * NEIGHBOURS + CellMap::xu_yu_z]); 

                if ( !(x < 1 || y < 1 || z < 1 || x >= nx+1 || y >= ny+1 || z >= nz+1) ) {
                    m->c_ghost[i] = 0;
                    VOXEL_MAP_NON_GHOST_PERIODIC(x-1, y  , z  , nx, ny, nz, m->c2ngc_map[i * FACES + Face::xd]);
                    VOXEL_MAP_NON_GHOST_PERIODIC(x  , y-1, z  , nx, ny, nz, m->c2ngc_map[i * FACES + Face::yd]);
                    VOXEL_MAP_NON_GHOST_PERIODIC(x  , y  , z-1, nx, ny, nz, m->c2ngc_map[i * FACES + Face::zd]);
                    VOXEL_MAP_NON_GHOST_PERIODIC(x+1, y  , z  , nx, ny, nz, m->c2ngc_map[i * FACES + Face::xu]);
                    VOXEL_MAP_NON_GHOST_PERIODIC(x  , y+1, z  , nx, ny, nz, m->c2ngc_map[i * FACES + Face::yu]);
                    VOXEL_MAP_NON_GHOST_PERIODIC(x  , y  , z+1, nx, ny, nz, m->c2ngc_map[i * FACES + Face::zu]);
                }

                if (x < 1 || y < 1 || z < 1) {
					m->c_mask_right[i] = 0;
				}
            }
        }
    }

    int from = -1;  
    for (int z = 1; z < nz+1; z++) // _zy_boundary
    {
        for (int y = 1; y < ny+1; y++)
        {
            from = VOXEL(1   , y, z, nx, ny, nz);
            m->c2cugb_map[from * 6 + 0] = VOXEL(nx+1, y, z, nx, ny, nz); 
            m->c2cugb0_map[from] = VOXEL(nx+1, y, z, nx, ny, nz); 
            m->c_mask_ugb[from * 6 + 0] = 1;

            from = VOXEL(nx  , y, z, nx, ny, nz);
            m->c2cugb_map[from * 6 + 1] = VOXEL(0   , y, z, nx, ny, nz);
            m->c2cugb1_map[from] = VOXEL(0   , y, z, nx, ny, nz);
            m->c_mask_ugb[from * 6 + 1] = 1;
        }
    }   
    for (int x = 0; x < nx+2; x++) // _xz_boundary
    {
        for (int z = 1; z < nz+1; z++)
        {
            from = VOXEL(x,    1, z, nx, ny, nz);
            m->c2cugb_map[from * 6 + 2] = VOXEL(x, ny+1, z, nx, ny, nz);
            m->c2cugb2_map[from] = VOXEL(x, ny+1, z, nx, ny, nz);
            m->c_mask_ugb[from * 6 + 2] = 1;

            from = VOXEL(x, ny  , z, nx, ny, nz);
            m->c2cugb_map[from * 6 + 3] = VOXEL(x, 0   , z, nx, ny, nz);
            m->c2cugb3_map[from] = VOXEL(x, 0   , z, nx, ny, nz);
            m->c_mask_ugb[from * 6 + 3] = 1;
        }
    }    
    for (int x = 0; x < nx+2; x++) // _yx_boundary
    {
        for (int y = 0; y < ny+2; y++)
        {
            from = VOXEL(x, y, 1   , nx, ny, nz);
            m->c2cugb_map[from * 6 + 4] = VOXEL(x, y, nz+1, nx, ny, nz);
            m->c2cugb4_map[from] = VOXEL(x, y, nz+1, nx, ny, nz);
            m->c_mask_ugb[from * 6 + 4] = 1;

            from = VOXEL(x, y, nz  , nx, ny, nz);
            m->c2cugb_map[from * 6 + 5] = VOXEL(x, y, 0   , nx, ny, nz);
            m->c2cugb5_map[from] = VOXEL(x, y, 0   , nx, ny, nz);
            m->c_mask_ugb[from * 6 + 5] = 1;
        }
    }

    for (int x = 0; x < nx+2; x++) // _x_boundary
    {
        for(int z = 1; z <= nz+1; z++){
            from = VOXEL(x, ny+1, z, nx, ny, nz);
            m->c2cug_map[from * 6 + 0] = VOXEL(x, 1   , z, nx, ny, nz);
            m->c2cug0_map[from] = VOXEL(x, 1   , z, nx, ny, nz);
            m->c_mask_ug[from * 6 + 0] = 1;
        }

        for(int y = 1; y <= ny+1; y++){
            from = VOXEL(x, y, nz+1, nx, ny, nz);
            m->c2cug_map[from * 6 + 1] = VOXEL(x, y, 1   , nx, ny, nz);
            m->c2cug1_map[from] = VOXEL(x, y, 1   , nx, ny, nz);
            m->c_mask_ug[from * 6 + 1] = 1;
        }
    }
    for (int y = 1; y < ny+1; y++) //_y_boundary
    {
        for (int x = 1; x <= nx+1; x++){
            from = VOXEL(x   , y, nz+1, nx, ny, nz);
            m->c2cug_map[from * 6 + 2] = VOXEL(x   , y, 1   , nx, ny, nz);
            m->c2cug2_map[from] = VOXEL(x   , y, 1   , nx, ny, nz);
            m->c_mask_ug[from * 6 + 2] = 1;
        }

        for (int z = 1; z <= nz+1; z++){
            from = VOXEL(nx+1, y, z   , nx, ny, nz);
            m->c2cug_map[from * 6 + 3] = VOXEL(1   , y, z   , nx, ny, nz);
            m->c2cug3_map[from] = VOXEL(1   , y, z   , nx, ny, nz);
            m->c_mask_ug[from * 6 + 3] = 1;
        }
    }
    for (int z = 1; z < nz+1; z++) //_z_boundary
    {
        for (int y = 1; y <= ny+1; y++){
            from = VOXEL(nx+1, y   , z, nx, ny, nz);
            m->c2cug_map[from * 6 + 4] = VOXEL(1   , y   , z, nx, ny, nz);
            m->c2cug4_map[from] = VOXEL(1   , y   , z, nx, ny, nz);
            m->c_mask_ug[from * 6 + 4] = 1;
        }

        for (int x = 1; x <= nx+1; x++){
            from = VOXEL(x   , ny+1, z, nx, ny, nz);
            m->c2cug_map[from * 6 + 5] = VOXEL(x   , 1   , z, nx, ny, nz);
            m->c2cug5_map[from] = VOXEL(x   , 1   , z, nx, ny, nz);
            m->c_mask_ug[from * 6 + 5] = 1;
        }
    }

    // TODO : Sanity Check : Check whether any mapping is still negative or not

    opp_printf("Setup", "init_mesh DONE");
}

/***************************************************************************************************
 * @brief This block distributes temporary DataPointers from ROOT rank to other ranks
 * @param g_m - Global mesh of temporary shared pointer of DataPointers, Root Rank should have data
 * @param m - rank specific block partitioned mesh of temporary shared pointer of DataPointers
 * @return (void)
 */
inline void distribute_data_over_ranks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m)
{ 
#ifdef USE_MPI
    MPI_Bcast(&(g_m->n_cells), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Bcast(&(g_m->n_particles), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);

    m = std::make_shared<DataPointers>();

    m->n_cells     = opp_get_uniform_local_size(g_m->n_cells);
    m->CreateMeshCommArrays();

    opp_uniform_scatter_array(g_m->c_index , m->c_index , g_m->n_cells, m->n_cells, ONE);
    opp_uniform_scatter_array(g_m->c_pos_ll, m->c_pos_ll, g_m->n_cells, m->n_cells, DIM);
    opp_uniform_scatter_array(g_m->c_ghost , m->c_ghost , g_m->n_cells, m->n_cells, ONE);
    opp_uniform_scatter_array(g_m->c_mask_right, m->c_mask_right , g_m->n_cells, m->n_cells, ONE);
    opp_uniform_scatter_array(g_m->c_mask_ug, m->c_mask_ug , g_m->n_cells, m->n_cells, 6);
    opp_uniform_scatter_array(g_m->c_mask_ugb, m->c_mask_ugb , g_m->n_cells, m->n_cells, 6);

    opp_uniform_scatter_array(g_m->c2c_map  , m->c2c_map  , g_m->n_cells , m->n_cells, NEIGHBOURS); 
    opp_uniform_scatter_array(g_m->c2ngc_map, m->c2ngc_map, g_m->n_cells , m->n_cells, FACES); 
    opp_uniform_scatter_array(g_m->c2cug_map  , m->c2cug_map  , g_m->n_cells , m->n_cells, 6); 
    opp_uniform_scatter_array(g_m->c2cugb_map  , m->c2cugb_map  , g_m->n_cells , m->n_cells, 6); 

    opp_uniform_scatter_array(g_m->c2cug0_map  , m->c2cug0_map  , g_m->n_cells , m->n_cells, 1); 
    opp_uniform_scatter_array(g_m->c2cug1_map  , m->c2cug1_map  , g_m->n_cells , m->n_cells, 1); 
    opp_uniform_scatter_array(g_m->c2cug2_map  , m->c2cug2_map  , g_m->n_cells , m->n_cells, 1); 
    opp_uniform_scatter_array(g_m->c2cug3_map  , m->c2cug3_map  , g_m->n_cells , m->n_cells, 1); 
    opp_uniform_scatter_array(g_m->c2cug4_map  , m->c2cug4_map  , g_m->n_cells , m->n_cells, 1); 
    opp_uniform_scatter_array(g_m->c2cug5_map  , m->c2cug5_map  , g_m->n_cells , m->n_cells, 1); 

    opp_uniform_scatter_array(g_m->c2cugb0_map  , m->c2cugb0_map  , g_m->n_cells , m->n_cells, 1); 
    opp_uniform_scatter_array(g_m->c2cugb1_map  , m->c2cugb1_map  , g_m->n_cells , m->n_cells, 1); 
    opp_uniform_scatter_array(g_m->c2cugb2_map  , m->c2cugb2_map  , g_m->n_cells , m->n_cells, 1); 
    opp_uniform_scatter_array(g_m->c2cugb3_map  , m->c2cugb3_map  , g_m->n_cells , m->n_cells, 1); 
    opp_uniform_scatter_array(g_m->c2cugb4_map  , m->c2cugb4_map  , g_m->n_cells , m->n_cells, 1); 
    opp_uniform_scatter_array(g_m->c2cugb5_map  , m->c2cugb5_map  , g_m->n_cells , m->n_cells, 1); 

    if (OPP_rank == OPP_ROOT)
        g_m->DeleteValues();
#else
    m = g_m;
#endif

    m->CreateMeshNonCommArrays();
}

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
#include "cabana_part_distribution.h"

using namespace opp;

int m0 = 0;
int m1 = 1;
int m2 = 2;
int m3 = 3;
int m4 = 4;
int m5 = 5;

void init_mesh(const Deck& deck, std::shared_ptr<DataPointers> m);
void distribute_data_over_ranks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m);

//*************************************************************************************************
/**
 * @brief Initialize the rank specific mesh data to a DataPointers utility class shared pointer
 * @return std::shared_ptr<DataPointers>
 */
std::shared_ptr<DataPointers> load_mesh(const Deck& deck) {

    std::shared_ptr<DataPointers> g_m(new DataPointers());

    if (OPP_rank == OPP_ROOT)      
        init_mesh(deck, g_m);

    std::shared_ptr<DataPointers> m;
    distribute_data_over_ranks(g_m, m);

    return m;
}

//*************************************************************************************************
/**
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

//*************************************************************************************************
/**
 * @brief Initializes the particles in to the particle dats in the arguments, using the cell_pos_ll dat
 *          Expect this to run on every MPI rank
 * @param part_index - opp_dat : Particle index relative to rank. TODO: make this global
 * @param part_pos - opp_dat : Particle 3D position (x,y,z)
 * @param part_vel - opp_dat : Particle 3D velocity (x,y,z)
 * @param part_streak_mid - opp_dat : Particle 3D temporary position (x,y,z)
 * @param part_weight - opp_dat : Particle weight
 * @param part_mesh_rel - opp_dat : Particle belonging cell index 
 * @param cell_pos_ll - opp_dat : Lower left 2D position coordicate of the cell
 * @return (void)
 */
void init_particles(const Deck& deck, opp_dat part_index, opp_dat part_pos, opp_dat part_vel, opp_dat part_streak_mid,
                    opp_dat part_weight, opp_dat part_mesh_rel, opp_dat cell_pos_ll, opp_dat cell_cgid, opp_dat cell_ghost) 
{
    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles START");

    const OPP_INT npart_per_cell = deck.nppc;

    const int all_cell_count = cell_pos_ll->set->size;
    int non_ghost_cell_count = 0;

    for (int i = 0; i < all_cell_count; i++) {
        if (((int*)cell_ghost->data)[i] == 0)
            non_ghost_cell_count++;
    }

    const int rank_npart = npart_per_cell * non_ghost_cell_count;
    int rank_part_start  = 0;

    if (rank_npart <= 0) {
        opp_printf("Setup", "Error No particles to add in rank %d", OPP_rank);
        return;
    }

#ifdef USE_MPI // canculate the starting particle index incase of MPI
    {
        std::vector<OPP_INT> temp(OPP_comm_size, 0);
        MPI_Allgather(&rank_npart, 1, MPI_INT, temp.data(), 1, MPI_INT, MPI_COMM_WORLD);
        for (int i = 0; i < OPP_rank; ++i) rank_part_start += temp[i];
    }  
#endif

    if (OPP_DBG)
        opp_printf("Setup", "%d parts to add in rank %d [part_start=%d]", rank_npart, OPP_rank, rank_part_start);

    // Host/Device space to store the particles.
    opp_increase_particle_count(part_pos->set, rank_npart);

    if (opp_params->get<OPP_STRING>("part_enrich") == "two_stream")
        enrich_particles_two_stream(deck, all_cell_count, (OPP_REAL*)part_pos->data, (OPP_REAL*)part_vel->data, 
            (OPP_REAL*)part_streak_mid->data, (OPP_INT*)part_mesh_rel->data, nullptr, // (OPP_INT*)part_index->data, 
            (OPP_REAL*)part_weight->data, (OPP_INT*)cell_cgid->data, (OPP_INT*)cell_ghost->data);
    else
        enrich_particles_random(deck, all_cell_count, (OPP_REAL*)part_pos->data, (OPP_REAL*)part_vel->data, 
            (OPP_REAL*)part_streak_mid->data, (OPP_INT*)part_mesh_rel->data, nullptr, // (OPP_INT*)part_index->data, 
            (OPP_REAL*)part_weight->data, (OPP_REAL*)cell_pos_ll->data, rank_part_start, (OPP_INT*)cell_ghost->data);

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles Uploading Start");

    opp_upload_particle_set(part_pos->set);

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles END");
}


//*************************************************************************************************
/**
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

//*************************************************************************************************
/**
 * @brief This uses MPI routines to colour the cell dats like a bundle of pencils along the direction of X
 *        this directional partitioning minimize the particle MPI communications
 * @param deck - 
 * @param cell_index - cell index dat of the current rank, this includes the global cell indices
 * @param cell_index - cell colours to be enriched for partitioning
 * @return (void)
 */
void cabana_color_pencil_x(const Deck& deck, opp_dat cell_index, const opp_dat cell_colors) 
{
#ifdef USE_MPI
    if (OPP_rank == OPP_ROOT) opp_printf("Setup", "x=%d y=%d z=%d", deck.nx, deck.ny, deck.nz);

    MPI_Comm comm_cart;
    int mpi_dims[3] = {0, 0, 0};
    int periods[3] = {1, 1, 1};
    int coords[3] = {0, 0, 0};
    int cell_starts[3] = {0, 0, 0}; // Holds the first cell this rank owns in each dimension.
    int cell_ends[3] = {1, 1, 1}; // Holds the last cell+1 this ranks owns in each dimension.
    
    const OPP_INT ndim = 2;
    std::vector<int> cell_counts = { deck.ny, deck.nz }; // x is the dir of velocity

    MPI_Dims_create(OPP_comm_size, ndim, mpi_dims);

    std::vector<int> cell_count_ordering = reverse_argsort(cell_counts); // direction with most cells first to match mpi_dims order

    std::vector<int> mpi_dims_reordered(ndim);
    for (int dimx = 0; dimx < ndim; dimx++)  // reorder the mpi_dims to match the actual domain
        mpi_dims_reordered[cell_count_ordering[dimx]] = mpi_dims[dimx];

    MPI_Cart_create(MPI_COMM_WORLD, ndim, mpi_dims_reordered.data(), periods, 1, &comm_cart);
    MPI_Cart_get(comm_cart, ndim, mpi_dims, periods, coords);

    for (int dimx = 0; dimx < ndim; dimx++) 
        get_decomp_1d(mpi_dims[dimx], cell_counts[dimx], coords[dimx], &cell_starts[dimx], &cell_ends[dimx]);

    for (int dimx = 0; dimx < ndim; dimx++) {
        cell_starts[dimx] += NG;
        cell_ends[dimx] += NG;
    }

    std::vector<int> all_cell_starts(OPP_comm_size * ndim);
    std::vector<int> all_cell_ends(OPP_comm_size * ndim);

    MPI_Allgather(cell_starts, ndim, MPI_INT, all_cell_starts.data(), ndim, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(cell_ends, ndim, MPI_INT, all_cell_ends.data(), ndim, MPI_INT, MPI_COMM_WORLD);

    if (OPP_rank == OPP_ROOT)
    {
        std::string log = "";
        for (int r = 0; r < OPP_comm_size; r++) {
            log += std::string("\nrank ") + std::to_string(r) + " start (";
            for (int d = 0; d < ndim; d++)
                log += std::to_string(all_cell_starts[r*ndim+d]) + ",";
            log += ") end (";
            for (int d = 0; d < ndim; d++)
                log += std::to_string(all_cell_ends[r*ndim+d]) + ",";
            log += ")";
        }
        opp_printf("cabana_color_pencil_x", "%s", log.c_str());
    }

#define CART_RANK_TO_INDEX(gcidx,ix,iy,iz,_x,_y) \
    int _ix, _iy, _iz;                                                    \
    _ix  = (gcidx);                        /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(_x);                   /* iy = iy+gpy*iz */            \
    _ix -= _iy*int(_x);                   /* ix = ix */                   \
    _iz  = _iy/int(_y);                   /* iz = iz */                   \
    _iy -= _iz*int(_y);                   /* iy = iy */                   \
    (ix) = _ix;                                                           \
    (iy) = _iy;                                                           \
    (iz) = _iz;                                                           \

    // used global id of cell and assign the color to the correct MPI rank
    const OPP_INT* gcid = ((OPP_INT*)cell_index->data);
    for (OPP_INT i = 0; i < cell_index->set->size; i++)
    {
        int coord[3] = {-1, -1 -1};
        CART_RANK_TO_INDEX(gcid[i], coord[0], coord[1], coord[2], deck.nx+2*NG, deck.ny+2*NG);

        if (coord[0] < NG) coord[0] = NG;
        else if (coord[0] >= (deck.nx + NG)) coord[0] = (deck.nx + NG - 1);
        if (coord[1] < NG) coord[1] = NG;
        else if (coord[1] >= (deck.ny + NG)) coord[1] = (deck.ny + NG - 1);
        if (coord[2] < NG) coord[2] = NG;
        else if (coord[2] >= (deck.nz + NG)) coord[2] = (deck.nz + NG - 1);

        for (int rank = 0; rank < OPP_comm_size; rank++) 
        {
            bool is_rank_suitable = true;

            if ((all_cell_starts[ndim*rank+0] > coord[1]) || (all_cell_ends[ndim*rank+0] <= coord[1]) || 
                (all_cell_starts[ndim*rank+1] > coord[2]) || (all_cell_ends[ndim*rank+1] <= coord[2]))
            {
                is_rank_suitable = false;
            }

            if (is_rank_suitable)
            {
                ((OPP_INT*)cell_colors->data)[i] = rank;
                break;
            }
        }    
    }

#undef CART_RANK_TO_INDEX
#endif
}

//*************************************************************************************************
/**
 * @brief This uses block colouring in YZ plane to colour the cell dats along the direction of X
 *        this directional partitioning minimize the particle MPI communications
 * @param deck - 
 * @param cell_index - cell index dat of the current rank, this includes the global cell indices
 * @param cell_index - cell colours to be enriched for partitioning
 * @return (void)
 */
void cabana_color_block_x(const Deck& deck, opp_dat cell_index, const opp_dat cell_colors) 
{
#ifdef USE_MPI

    const OPP_INT numClusters = deck.ny * deck.nz;
    std::vector<int> assignments;
    std::vector<int> clusterSizes(numClusters, 0);

    for (int r = 0; r < OPP_comm_size; r++) {
        int countForRank = opp_get_uniform_local_size(numClusters, r);
        for (int i = 0; i < countForRank; i++) {
            assignments.emplace_back(r);
            clusterSizes[r]++;
        }
    }

#define CART_RANK_TO_INDEX(gcidx,ix,iy,iz,_x,_y) \
    int _ix, _iy, _iz;                                                    \
    _ix  = (gcidx);                        /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(_x);                   /* iy = iy+gpy*iz */            \
    _ix -= _iy*int(_x);                   /* ix = ix */                   \
    _iz  = _iy/int(_y);                   /* iz = iz */                   \
    _iy -= _iz*int(_y);                   /* iy = iy */                   \
    (ix) = _ix;                                                           \
    (iy) = _iy;                                                           \
    (iz) = _iz;                                                           \

    // used global id of cell and assign the color to the correct MPI rank
    const OPP_INT* gcid = ((OPP_INT*)cell_index->data);
    for (OPP_INT i = 0; i < cell_index->set->size; i++)
    {
        int coord[3] = {-1, -1 -1};
        CART_RANK_TO_INDEX(gcid[i], coord[0], coord[1], coord[2], deck.nx+2*NG, deck.ny+2*NG);    
        
        if (coord[1] >= (deck.ny + NG)) coord[1] -= 2*NG;
        else if (coord[1] >= NG) coord[1] -= NG;
        if (coord[2] >= (deck.nz + NG)) coord[2] -= 2*NG;
        else if (coord[2] >= NG) coord[2] -= NG;

        const int idx = (coord[1] + deck.ny * coord[2]);
        ((OPP_INT*)cell_colors->data)[i] = assignments[idx]; 
    }

#undef CART_RANK_TO_INDEX

#endif
}
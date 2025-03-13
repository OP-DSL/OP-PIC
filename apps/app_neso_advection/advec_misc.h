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

class DataPointers;
void init_mesh(std::shared_ptr<DataPointers> m);
void distribute_data_over_ranks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m);

//*************************************************************************************************
/**
 * @brief Utility class to temporarily hold the mesh data until it is loaded by OP-PIC
 */
class DataPointers // This is just a placeholder for initializing // No use in DSL
{
    public:
        DataPointers() {}
        virtual ~DataPointers()
        {
            DeleteValues();   
        }

        inline void DeleteValues() {

            if (c_index) delete[] c_index;
            if (c_pos_ll) delete[] c_pos_ll;
            if (c_colors) delete[] c_colors;
            if (cell_cell_map) delete[] cell_cell_map;

            c_index = nullptr;
            c_pos_ll = nullptr;
            c_colors = nullptr;
            cell_cell_map = nullptr;
        }

        inline void CreateMeshArrays() {

            this->c_index       = new OPP_INT[this->n_cells * ONE];
            this->c_pos_ll      = new OPP_REAL[this->n_cells * DIM];  
            this->c_colors      = new OPP_INT[this->n_cells * ONE];        
            this->cell_cell_map = new OPP_INT[this->n_cells * NEIGHBOURS];

            for (int i = 0; i < this->n_cells * NEIGHBOURS; i++) this->cell_cell_map[i] = -1;
            
            for (int i = 0; i < this->n_cells; i++) this->c_colors[i] = 0;
        }

        int n_cells     = 0;

        OPP_INT* c_index       = nullptr;
        OPP_REAL* c_pos_ll   = nullptr;
        OPP_INT* c_colors      = nullptr;
        OPP_INT* cell_cell_map = nullptr;
};

//*************************************************************************************************
/**
 * @brief Initialize the rank specific mesh data to a DataPointers utility class shared pointer
 * @return std::shared_ptr<DataPointers>
 */
std::shared_ptr<DataPointers> LoadData() {

    std::shared_ptr<DataPointers> g_m(new DataPointers());

    if (OPP_rank == OPP_ROOT)      
        init_mesh(g_m);

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
void init_mesh(std::shared_ptr<DataPointers> m) {

    OPP_INT nx = opp_params->get<OPP_INT>("nx");
    OPP_INT ny = opp_params->get<OPP_INT>("ny");
    OPP_REAL cell_width = opp_params->get<OPP_REAL>("cell_width");

    m->n_cells = (nx * ny);

    m->CreateMeshArrays();

	opp_printf("Setup", "init_mesh n_cells=%d nx=%d ny=%d", m->n_cells, nx, ny);

	for (int x = 0; x < nx; x++) {
		for (int y = 0; y < ny; y++) {
 
            const int i = VOXEL(x,y, nx);
            
            m->c_index[i] = i;
            m->c_pos_ll[i*DIM + Dim::x] = x * cell_width;
            m->c_pos_ll[i*DIM + Dim::y] = y * cell_width;

            VOXEL_MAP(x-1, y,   nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::xd_y]);
            VOXEL_MAP(x+1, y,   nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::xu_y]);
            VOXEL_MAP(x,   y-1, nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_yd]);
            VOXEL_MAP(x,   y+1, nx, ny, m->cell_cell_map[i * NEIGHBOURS + CellMap::x_yu]);
        }
    }
}

//*************************************************************************************************
/**
 * @brief Initializes the particles in to the particle dats in the arguments, using the c_pos_ll dat
 *          Expect this to run on every MPI rank
 * @param p_id - opp_dat : Particle index relative to rank. TODO: make this global
 * @param p_pos - opp_dat : Particle 2D position (x,y)
 * @param p_vel - opp_dat : Particle 2D velocity (x,y)
 * @param p_cid - opp_map : Particle belonging cell index 
 * @param c_pos_ll - opp_dat : Lower left 2D position coordicate of the cell
 * @return (void)
 */
void init_particles(opp_dat p_id, opp_dat p_pos, opp_dat p_vel, opp_dat p_mdir, opp_map p_cid, 
                    opp_dat c_pos_ll, opp_dat gcid) 
{
    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles START");

    OPP_INT npart_per_cell = opp_params->get<OPP_INT>("npart_per_cell");

    std::mt19937 rng_pos(52234234 + OPP_rank);
    std::mt19937 rng_vel(52234231 + OPP_rank);

    const int cell_count = c_pos_ll->set->size;
    const int rank_npart = npart_per_cell * cell_count;

    if (rank_npart <= 0) {
        opp_printf("Setup", "Error No particles to add in rank %d", OPP_rank);
    }

    std::vector<std::vector<double>> positions;
    std::vector<int> cells;

    // Sample particles randomly in each local cell.
    opp_profiler->start("Init_get_pos");
    const double extents[] = { opp_params->get<OPP_REAL>("cell_width"), opp_params->get<OPP_REAL>("cell_width")};
    uniform_within_cartesian_cells(DIM, extents, (OPP_REAL*)c_pos_ll->data, cell_count, npart_per_cell, 
                                positions, cells, rng_pos);
    opp_profiler->end("Init_get_pos");

    // Sample some particle velocities.
    opp_profiler->start("Init_get_vel");
    auto velocities = get_normal_distribution(rank_npart, DIM, 0.0, 0.5, rng_vel);
    opp_profiler->end("Init_get_vel");

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles opp_increase_particle_count Start rank_npart=%d", rank_npart);

    // Host/Device space to store the particles.
    opp_profiler->start("Init_inc_parts");
    opp_increase_particle_count(p_id->set, rank_npart);
    opp_profiler->end("Init_inc_parts");

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles Load data to dats Start");

    opp_profiler->start("Init_assign");
    // Populate the host space with particle data.
    for (int px = 0; px < rank_npart; px++) {
        
        opp_get_data<OPP_REAL>(p_pos)[px * DIM + Dim::x] = positions.at(Dim::x).at(px);
        opp_get_data<OPP_REAL>(p_pos)[px * DIM + Dim::y] = positions.at(Dim::y).at(px);
        opp_get_data<OPP_REAL>(p_vel)[px * DIM + Dim::x] = velocities.at(Dim::x).at(px);
        opp_get_data<OPP_REAL>(p_vel)[px * DIM + Dim::y] = velocities.at(Dim::y).at(px);
        opp_get_data<OPP_INT>(p_mdir)[px * DIM + Dim::x] = 1;
        opp_get_data<OPP_INT>(p_mdir)[px * DIM + Dim::y] = 1;

        opp_get_data(p_cid)[px]       = cells.at(px);
        opp_get_data<OPP_INT>(p_id)[px] = px; //((int*)gcid->data)[cells.at(px)]; 
    }
    opp_profiler->end("Init_assign");

    if (OPP_rank == OPP_ROOT)
        opp_printf("Setup", "Init particles Uploading Start");

    opp_profiler->start("Init_upload");
    opp_upload_particle_set(p_id->set);
    opp_profiler->end("Init_upload");

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
void distribute_data_over_ranks(std::shared_ptr<DataPointers>& g_m, std::shared_ptr<DataPointers>& m)
{ 
#ifdef USE_MPI
    MPI_Bcast(&(g_m->n_cells), 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);

    m = std::make_shared<DataPointers>();
    m->n_cells     = opp_get_uniform_local_size(g_m->n_cells);
    m->CreateMeshArrays();

    opp_uniform_scatter_array(g_m->c_index      , m->c_index      , g_m->n_cells, m->n_cells, ONE);
    opp_uniform_scatter_array(g_m->c_pos_ll     , m->c_pos_ll     , g_m->n_cells, m->n_cells, DIM); 
    opp_uniform_scatter_array(g_m->c_colors     , m->c_colors     , g_m->n_cells, m->n_cells, ONE); 
    opp_uniform_scatter_array(g_m->cell_cell_map, m->cell_cell_map, g_m->n_cells, m->n_cells, NEIGHBOURS); 

    if (OPP_rank == OPP_ROOT)
        g_m->DeleteValues();
    
#else
    m = g_m;
#endif
}

/*************************************************************************************************
 * This is a utility function to get the total particles iterated during the simulation
*/
inline void get_global_values(const int64_t total_part_iter, int64_t& gbl_total_part_iter, 
    int64_t& gbl_max_iter, int64_t& gbl_min_iter) {

#ifdef USE_MPI
        MPI_Reduce(&total_part_iter, &gbl_total_part_iter, 1, MPI_INT64_T, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&total_part_iter, &gbl_max_iter, 1, MPI_INT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_part_iter, &gbl_min_iter, 1, MPI_INT64_T, MPI_MIN, 0, MPI_COMM_WORLD);
#else
        gbl_total_part_iter = total_part_iter;
        gbl_max_iter = total_part_iter;
        gbl_min_iter = total_part_iter;
#endif
}

/***************************************************************************************************
 * @brief This uses block colouring to colour the cell dats according to global index
 * @param cell_index - cell index dat of the current rank, this includes the global cell indices
 * @param cell_colors - cell colours to be enriched for partitioning
 * @return (void)
 */
inline void advec_color_block(const opp_dat cell_index, opp_dat cell_colors) 
{
#ifdef USE_MPI
    // used global id of cell and assign the color to the correct MPI rank
    // const OPP_INT* gcid = ((OPP_INT*)cell_index->data);
    for (OPP_INT i = 0; i < cell_index->set->size; i++)
    {
        ((OPP_INT*)cell_colors->data)[i] = OPP_rank; 
    }
#endif
}
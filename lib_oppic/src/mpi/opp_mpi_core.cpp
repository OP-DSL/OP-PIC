
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

#include <opp_mpi_core.h>

MPI_Comm OP_MPI_WORLD;
MPI_Comm OP_MPI_GLOBAL;

void opp_halo_create();
void opp_part_comm_init();

//*******************************************************************************
void opp_partition_core(std::string lib_name, op_set prime_set, op_map prime_map, op_dat data)
{
    if (lib_name == "PARMETIS_KWAY")
    {
#ifdef HAVE_PARMETIS
        if (prime_map != NULL)
        {
            opp_partition_kway(prime_map); // use parmetis kway partitioning
        }
        else
        {
            opp_abort("opp_partition PARMETIS_KWAY Error: Partitioning prime_map : NULL - UNSUPPORTED Partitioner Specification");  
        }
#else
        opp_abort("opp_partition_core PARMETIS_KWAY Error: Parmetis not installed or not defined");
#endif
    }
    else if (lib_name == "PARMETIS_GEOM")
    {
#ifdef HAVE_PARMETIS
        if (data != NULL)
        {
            opp_partition_geom(data); // use parmetis geometric partitioning
        }
        else
        {
            opp_abort("opp_partition PARMETIS_GEOM Error: Partitioning geom dat : NULL - UNSUPPORTED Partitioner Specification"); 
        }
#else
        opp_abort("opp_partition_core PARMETIS_GEOM Error: Parmetis not installed or not defined");
#endif
    }
    else if (lib_name == "EXTERNAL")
    {
        if (data != NULL)
        {
            opp_partition_external(prime_set, data); // use external partitioning dat
        }
        else
        {
            opp_abort("opp_partition EXTERNAL Error: Partitioning color dat : NULL - UNSUPPORTED Partitioner Specification"); 
        }
    }
    else if (lib_name != "")
    {
        opp_abort("opp_partition Error: Unsupported lib_name - UNSUPPORTED Partitioner Specification");
    }

    opp_halo_create();

    opp_part_comm_init(); 

    std::vector<std::vector<int>> set_sizes(oppic_sets.size());

    for (oppic_set set : oppic_sets)
    {
        std::vector<int>& recv_vec = set_sizes[set->index];
        recv_vec.resize(OPP_comm_size * 3);

        std::vector<int> sizes{ set->size, set->exec_size, set->nonexec_size };
        MPI_Gather(&(sizes[0]), 3, MPI_INT, &(recv_vec[0]), 3, MPI_INT, OPP_ROOT, OP_MPI_WORLD);
    }

    // print the set sizes of all ranks after partitioning
    if (OPP_rank == OPP_ROOT)
    {
        std::string log = "";

        for (oppic_set set : oppic_sets)
            log += "\t - " + std::string(set->name);

        opp_printf("opp_partition()", "(size|ieh|inh) %s", log.c_str());

        for (int i = 0; i < OPP_comm_size; i++)
        {
            log = "RANK [" + std::to_string(i) + "]";
            
            for (int j = 0; j < (int)oppic_sets.size(); j++)
                log += "\t- " + std::to_string(set_sizes[j][i * 3]) + "|" + 
                    std::to_string(set_sizes[j][i * 3 + 1]) + "|" + std::to_string(set_sizes[j][i * 3 + 2]);

            opp_printf("opp_partition()", "%s", log.c_str());
        }
    }
}

std::map<int, opp_dat> negative_mapping_indices;

//*******************************************************************************
void opp_sanitize_all_maps()
{
    for (int i = 0; i < (int)oppic_maps.size(); i++)
    {
        oppic_map map = oppic_maps[i];

        if (map->dim == 1) continue;

        if (OP_DEBUG) opp_printf("opp_sanitize_all_maps", " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

        std::string name = std::string("AUTO_DAT_") + map->name;
        oppic_dat dat = opp_decl_mesh_dat(map->from, map->dim, DT_INT, (char*)map->map, name.c_str());  
        negative_mapping_indices[map->index] = dat;

        memset(dat->data, 0, (map->from->size * dat->size));

        for (int n = 0; n < map->from->size; n++)
        {
            int positive_mapping = -1;
            std::vector<int> index;

            for (int d = 0; d < map->dim; d++)
            {
                if (map->map[n * map->dim + d] < 0)
                {
                    index.push_back(n * map->dim + d);
                }
                else
                {
                    positive_mapping = map->map[n * map->dim + d];
                }
            }

            if (positive_mapping >= 0)
            {
                for (int i : index)
                {
                    map->map[i] = positive_mapping;
                    ((int*)dat->data)[i] = 1;
                }
            }
            else
            {
                opp_printf("opp_sanitize_all_maps", "Error: No positive mapping found at %d in map: %s", n, map->name);
            }
        }
    }

    // if (OP_DEBUG)
    // {
    //     for (int i = 0; i < oppic_maps.size(); i++)
    //     {
    //         oppic_map map = oppic_maps[i];
            
    //         opp_printf("opp_sanitize_all_maps", OPP_rank, " map: %s | from->size: %d | dim: %d", map->name, map->from->size, map->dim);

    //         for (int n = 0; n < map->from->size; n++)
    //         {
    //             for (int d = 1; d < map->dim; d++)
    //             {
    //                 if (map->map[n * map->dim + d] < 0)
    //                 {
    //                     opp_printf("opp_sanitize_all_maps", OPP_rank, "Error: map: %s | ptr: %p | negative mapping at index: %d [%d]", map->name, map->map, n, n * map->dim + d);
    //                 }
    //             }
    //         }
    //     }
    // }
}

//*******************************************************************************
void opp_desanitize_all_maps()
{
    for (int i = 0; i < (int)oppic_maps.size(); i++)
    {
        oppic_map map = oppic_maps[i];

        if (map->dim == 1) continue;

        if (OP_DEBUG) opp_printf("opp_desanitize_all_maps", " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

        auto it = negative_mapping_indices.find(map->index);
        if (it == negative_mapping_indices.end())
        {
            opp_printf("opp_desanitize_all_maps", "Error: Negative mappings not found for map: %s", map->name);
            continue;
        }
            
        oppic_dat dat = it->second;

        for (int x = 0; x < (map->from->size + map->from->exec_size) * map->dim; x++)
        {
            if (((int*)dat->data)[x] == 1)
                map->map[x] = -1;
        }
    }

    // could realloc the dat to a lower size, if required
    negative_mapping_indices.clear();
}

//****************************************
void opp_get_start_end(opp_set set, opp_reset reset, int& start, int& end)
{
    switch (reset)
    {
        case OPP_Reset_Core:
            start = 0;
            end = set->core_size;
            break;
        case OPP_Reset_Set:
            start = 0;
            end = set->size;
            break;
        case OPP_Reset_All:
            start = 0;
            end = set->size + set->exec_size + set->nonexec_size;
            break;
        case OPP_Reset_ieh:
            start = set->size;
            end = set->size + set->exec_size;
            break;
        case OPP_Reset_inh:
            start = set->size + set->exec_size;
            end = set->size + set->exec_size + set->nonexec_size;
            break;
        default:
            opp_printf("opp_get_start_end", "Error: opp_reset failure");
    }
}

//*******************************************************************************
//*******************************************************************************

using namespace opp;

Comm::Comm(MPI_Comm comm_parent) {
    this->comm_parent = comm_parent;

    int rank_parent;
    CHECK(MPI_Comm_rank(comm_parent, &rank_parent))
    CHECK(MPI_Comm_split_type(comm_parent, MPI_COMM_TYPE_SHARED, 0,
                            MPI_INFO_NULL, &this->comm_intra))

    int rank_intra;
    CHECK(MPI_Comm_rank(this->comm_intra, &rank_intra))
    const int colour_intra = (rank_intra == 0) ? 1 : MPI_UNDEFINED;
    CHECK(MPI_Comm_split(comm_parent, colour_intra, 0, &this->comm_inter))

    CHECK(MPI_Comm_rank(this->comm_parent, &this->rank_parent))
    CHECK(MPI_Comm_rank(this->comm_intra, &this->rank_intra))
    CHECK(MPI_Comm_size(this->comm_parent, &this->size_parent))
    CHECK(MPI_Comm_size(this->comm_intra, &this->size_intra))
    
    if (comm_inter != MPI_COMM_NULL) {

        CHECK(MPI_Comm_rank(this->comm_inter, &this->rank_inter))
        CHECK(MPI_Comm_size(this->comm_inter, &this->size_inter))
    }

    if (OP_DEBUG)
        opp_printf("Comm", "rank_parent %d|s=%d rank_intra %d|s=%d rank_inter %d|s=%d",
            this->rank_parent, this->size_parent, this->rank_intra, this->size_intra, 
            this->rank_inter, this->size_inter);
};

Comm::~Comm() {

    if ((this->comm_intra != MPI_COMM_NULL) && (this->comm_intra != MPI_COMM_WORLD)) {
        
        CHECK(MPI_Comm_free(&this->comm_intra))
        this->comm_intra = MPI_COMM_NULL;
    }

    if ((this->comm_inter != MPI_COMM_NULL) && (this->comm_inter != MPI_COMM_WORLD)) {

        CHECK(MPI_Comm_free(&this->comm_inter))
        this->comm_intra = MPI_COMM_NULL;
    }
}

//*************************************************************************************************
// ndim : can be 2 or 3
// cell_counts : cell_counts in each direction
// cell_index : cell_index dat which holds global numbering
// cell_colors : local cell_colors dat to colour with most appropriate MPI rank
void opp_colour_cartesian_mesh(const int ndim, const std::vector<int> cell_counts, opp_dat cell_index, 
                            const opp_dat cell_colors)
{

    MPI_Comm comm_cart;
    int mpi_dims[3] = {0, 0, 0};
    int periods[3] = {1, 1, 1};
    int coords[3] = {0, 0, 0};
    int cell_starts[3] = {0, 0, 0}; // Holds the first cell this rank owns in each dimension.
    int cell_ends[3] = {1, 1, 1}; // Holds the last cell+1 this ranks owns in each dimension.

    MPI_Dims_create(OPP_comm_size, ndim, mpi_dims);

    std::vector<int> cell_count_ordering = reverse_argsort(cell_counts); // direction with most cells first to match mpi_dims order

    std::vector<int> mpi_dims_reordered(ndim);
    for (int dimx = 0; dimx < ndim; dimx++)  // reorder the mpi_dims to match the actual domain
        mpi_dims_reordered[cell_count_ordering[dimx]] = mpi_dims[dimx];

    MPI_Cart_create(MPI_COMM_WORLD, ndim, mpi_dims_reordered.data(), periods, 1, &comm_cart);
    MPI_Cart_get(comm_cart, ndim, mpi_dims, periods, coords);

    for (int dimx = 0; dimx < ndim; dimx++) 
        get_decomp_1d(mpi_dims[dimx], cell_counts[dimx], coords[dimx], &cell_starts[dimx], &cell_ends[dimx]);

    std::vector<int> all_cell_starts(OPP_comm_size * ndim);
    std::vector<int> all_cell_ends(OPP_comm_size * ndim);

    MPI_Allgather(cell_starts, ndim, MPI_INT, all_cell_starts.data(), ndim, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(cell_ends, ndim, MPI_INT, all_cell_ends.data(), ndim, MPI_INT, MPI_COMM_WORLD);

    if (OPP_rank == OPP_ROOT)
    {
        std::string log = "";
        for (int r = 0; r < OPP_comm_size; r++) {
            log += std::string("\nrank ") + std::to_string(r) + " start (" + 
                std::to_string(all_cell_starts[r*ndim+0]) + "," + std::to_string(all_cell_starts[r*ndim+1]) +
                ") end (" + std::to_string(all_cell_ends[r*ndim+0]) + "," + std::to_string(all_cell_ends[r*ndim+1]) + ")";
        }
        opp_printf("opp_colour_cartesian_mesh", "%s", log.c_str());
    }

#define CART_RANK_TO_INDEX(rank,ix,iy,iz,_x,_y) \
    int _ix, _iy, _iz;                                                    \
    _ix  = (rank);                        /* ix = ix+gpx*( iy+gpy*iz ) */ \
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
        CART_RANK_TO_INDEX(gcid[i], coord[0], coord[1], coord[2], cell_counts[0], cell_counts[1]);

        for (int rank = 0; rank < OPP_comm_size; rank++) 
        {
            bool is_rank_suitable[3] = { true, true, true };

            for (int dimx = 0; dimx < ndim; dimx++) 
            {
                if ((all_cell_starts[ndim*rank+dimx] > coord[dimx]) || 
                    (all_cell_ends[ndim*rank+dimx] <= coord[dimx]))
                {
                    is_rank_suitable[dimx] = false;
                    break;
                }
            }

            if (is_rank_suitable[0] && is_rank_suitable[1] && is_rank_suitable[2])
            {
                ((OPP_INT*)cell_colors->data)[i] = rank;
                break;
            }
        }    
    }

#undef CART_RANK_TO_INDEX
}
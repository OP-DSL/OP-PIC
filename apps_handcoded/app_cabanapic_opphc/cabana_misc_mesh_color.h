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

// include only to cabana_misc.h
#pragma once

/***************************************************************************************************
 * @brief This uses MPI routines to colour the cell dats like a bundle of pencils along the direction of X
 *        this directional partitioning minimize the particle MPI communications
 * @param deck - 
 * @param cell_index - cell index dat of the current rank, this includes the global cell indices
 * @param cell_index - cell colours to be enriched for partitioning
 * @return (void)
 */
inline void cabana_color_pencil(const Deck& deck, opp_dat cell_index, const opp_dat cell_colors) 
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
        opp_printf("cabana_color_pencil", "%s", log.c_str());
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

/***************************************************************************************************
 * @brief This uses block colouring in YZ plane to colour the cell dats along the direction of X
 *        this directional partitioning minimize the particle MPI communications
 * @param deck - 
 * @param cell_index - cell index dat of the current rank, this includes the global cell indices
 * @param cell_colors - cell colours to be enriched for partitioning
 * @return (void)
 */
inline void cabana_color_block(const Deck& deck, opp_dat cell_index, const opp_dat cell_colors) 
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


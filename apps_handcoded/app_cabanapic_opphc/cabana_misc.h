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
#include "cabana_misc_mesh_color.h"
#include "cabana_misc_mesh_loader.h"

int m0 = 0, m1 = 1, m2 = 2, m3 = 3, m4 = 4, m5 = 5;

void enrich_particles_random(const Deck& deck, const OPP_INT cell_count, OPP_REAL* pos, 
        OPP_REAL* vel, OPP_REAL* str_mid, OPP_INT* cid, OPP_REAL* weight, 
        const OPP_REAL* cell_pos_ll, const OPP_INT rank_part_start, const OPP_INT* ghost);
void enrich_particles_two_stream(const Deck& deck, const OPP_INT cell_count, 
        OPP_REAL* pos,  OPP_REAL* vel, OPP_REAL* str_mid, OPP_INT* cid, 
        OPP_REAL* weight, const OPP_INT* global_cids, const OPP_INT* ghost);

/***************************************************************************************************
 * @brief Initializes the particles in to the particle dats in the arguments, using the cell_pos_ll dat
 *          Expect this to run on every MPI rank
 * @param part_pos - opp_dat : Particle 3D position (x,y,z)
 * @param part_vel - opp_dat : Particle 3D velocity (x,y,z)
 * @param part_streak_mid - opp_dat : Particle 3D temporary position (x,y,z)
 * @param part_weight - opp_dat : Particle weight
 * @param part_mesh_rel - opp_dat : Particle belonging cell index 
 * @param cell_pos_ll - opp_dat : Lower left 2D position coordicate of the cell
 * @return (void)
 */
void init_particles(const Deck& deck,  opp_dat part_pos, opp_dat part_vel, opp_dat part_streak_mid,
                    opp_dat part_weight, opp_dat part_mesh_rel, opp_dat cell_pos_ll, opp_dat cell_cgid, opp_dat cell_ghost) 
{
    OPP_RUN_ON_ROOT() opp_printf("Setup", "Init particles START");

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
            (OPP_REAL*)part_streak_mid->data, (OPP_INT*)part_mesh_rel->data,
            (OPP_REAL*)part_weight->data, (OPP_INT*)cell_cgid->data, (OPP_INT*)cell_ghost->data);
    else
        enrich_particles_random(deck, all_cell_count, (OPP_REAL*)part_pos->data, (OPP_REAL*)part_vel->data, 
            (OPP_REAL*)part_streak_mid->data, (OPP_INT*)part_mesh_rel->data,
            (OPP_REAL*)part_weight->data, (OPP_REAL*)cell_pos_ll->data, rank_part_start, (OPP_INT*)cell_ghost->data);

    OPP_RUN_ON_ROOT() opp_printf("Setup", "Init particles Uploading Start");

    opp_upload_particle_set(part_pos->set);

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    OPP_RUN_ON_ROOT() opp_printf("Setup", "Init particles END");
}

/*************************************************************************************************
 * This function allocates data to particles randomly using a uniform distribution #similar to NESO-advection
*/
void enrich_particles_random(const Deck& deck, const OPP_INT cell_count, OPP_REAL* pos, 
        OPP_REAL* vel, OPP_REAL* str_mid, OPP_INT* cid, OPP_REAL* weight, 
        const OPP_REAL* cell_pos_ll, const OPP_INT rank_part_start, const OPP_INT* ghost) {
    
    OPP_RUN_ON_ROOT() opp_printf("Setup", "enrich_particles_random");

    const OPP_INT npart_per_cell = deck.nppc;
    const OPP_REAL const_weight  = deck.we;
    const OPP_REAL init_vel      = deck.v0;
    const OPP_REAL extents[DIM]  = { deck.dx, deck.dy, deck.dz };

    std::mt19937 rng(52234234);
    std::uniform_real_distribution<OPP_REAL> uniform_rng(0.0, 1.0);

    std::vector<OPP_REAL> uniform_dist_vec(npart_per_cell * DIM);
    for (int px = 0; px < npart_per_cell; px++) {
        for (int dimx = 0; dimx < DIM; dimx++) {
            uniform_dist_vec[px * DIM + dimx] = extents[dimx] * uniform_rng(rng);
        }
    }

    // Populate the host space with particle data.
    int p_idx = 0;
    for (int cx = 0; cx < cell_count; cx++) {
        if (ghost[cx] == 1) continue;
        for (int px = 0; px < npart_per_cell; px++) {

            pos[p_idx * DIM + Dim::x] = uniform_dist_vec[px * DIM + Dim::x];
            pos[p_idx * DIM + Dim::y] = uniform_dist_vec[px * DIM + Dim::y];
            pos[p_idx * DIM + Dim::z] = uniform_dist_vec[px * DIM + Dim::z];

            const int sign = (px % 2) ? 1 : -1; // For two stream
            // const int sign = 1; // For one stream

            vel[p_idx * DIM + Dim::x] = 0.0;
            vel[p_idx * DIM + Dim::y] = sign * init_vel;
            vel[p_idx * DIM + Dim::z] = 0.0;

            str_mid[p_idx * DIM + Dim::x] = 0.0;
            str_mid[p_idx * DIM + Dim::y] = 0.0;
            str_mid[p_idx * DIM + Dim::z] = 0.0;

            cid[p_idx]    = cx;
            weight[p_idx] = const_weight;

            p_idx++;
        }
    }
}

/*************************************************************************************************
 * This is similar to the original CabanaPIC particle enrichment
*/
void enrich_particles_two_stream(const Deck& deck, const OPP_INT cell_count, 
        OPP_REAL* pos,  OPP_REAL* vel, OPP_REAL* str_mid, OPP_INT* cid, 
        OPP_REAL* weight, const OPP_INT* global_cids, const OPP_INT* ghost) {
    
    OPP_RUN_ON_ROOT() opp_printf("Setup", "enrich_particles_two_stream");

    const OPP_INT npart_per_cell = deck.nppc;
    const OPP_REAL const_weight  = deck.we;
    const OPP_REAL v0            = deck.v0;
    const OPP_INT nx             = deck.nx;
    const OPP_INT ny             = deck.ny;
    const OPP_REAL dxp           = (2.0 / npart_per_cell);

    // Populate the host space with particle data.
    int part_idx = 0;
    for (int cx = 0; cx < cell_count; cx++) {
        
        if (ghost[cx] == 1) {
            continue;
        }
        const int global_cid = global_cids[cx];
        const int cell_particle_start = (global_cid * npart_per_cell);

        int ix, iy, iz;
        RANK_TO_INDEX(global_cid,ix,iy,iz,nx+2*NG,ny+2*NG);
        ix-=1; iy-=1; iz-=1;
        const int pre_ghost = (ix + nx * (iy + (ny * iz)));

        for (int p_idx = 0; p_idx < npart_per_cell; p_idx++) {
            
            int sign =  -1;
            const size_t pi2 = (cell_particle_start + p_idx); // TODO : looks like cell_particle_start is wrong due to ghosts
            const size_t pi = ((pi2) / 2);
            if (pi2 % 2 == 0)
                sign = 1;

            const int pic = (2 * pi) % npart_per_cell; //Every 2 particles have the same "pic".
            const double x = pic * dxp + 0.5 * dxp - 1.0;

            // Initialize velocity.(each cell length is 2)
            const double gam = 1.0 / sqrt(1.0 - v0 * v0);

            const double na = 0.0001 * sin(2.0 * 3.1415926 * ((x + 1.0 + pre_ghost * 2) / (2 * ny)));
            
            pos[part_idx * DIM + Dim::x] = 0.0;
            pos[part_idx * DIM + Dim::y] = x;
            pos[part_idx * DIM + Dim::z] = 0.0;

            vel[part_idx * DIM + Dim::x] = sign * v0 * gam * (1.0 + na * sign);
            vel[part_idx * DIM + Dim::y] = 0; 
            // vel[part_idx * DIM + Dim::y] = sign * 0.1 * v0 * gam * (1.0 + na * sign);
            vel[part_idx * DIM + Dim::z] = 0;

            str_mid[part_idx * DIM + Dim::x] = 0.0;
            str_mid[part_idx * DIM + Dim::y] = 0.0;
            str_mid[part_idx * DIM + Dim::z] = 0.0;

            cid[part_idx]    = cx;
            weight[part_idx] = const_weight;  

            part_idx++;
        }
    }
}

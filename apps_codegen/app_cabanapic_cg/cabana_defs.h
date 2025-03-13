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

#include <opp_lib.h>
#include <math.h>

#define NG                        1
#define ONE                       1
#define TWO                       2
#define DIM                       3
#define ACCUMULATOR_ARRAY_LENGTH  4
#define INTERP_LEN                18
#define ACC_LEN (ACCUMULATOR_ARRAY_LENGTH*DIM)

#define VOXEL(x,y,z, nx,ny,nz) ((x) + ((nx)+(NG*2))*((y) + ((ny)+(NG*2))*(z)))

#define VOXEL_MAP(_ix,_iy,_iz, nx,ny,nz, OUT) \
    { \
    int _x = _ix; \
    int _y = _iy; \
    int _z = _iz; \
    if ((_x < 0) || (_y < 0) || (_z < 0) || (_x > nx+NG) || (_y > ny+NG) || (_z > nz+NG)) OUT=-1; \
    else OUT = ((_x) + ((nx)+(NG*2))*((_y) + ((ny)+(NG*2))*(_z))); } \

#define VOXEL_MAP_NON_GHOST_PERIODIC(_ix,_iy,_iz, nx,ny,nz, OUT) \
    { \
    int _x = _ix; \
    int _y = _iy; \
    int _z = _iz; \
    if (_x == 0) _x = (nx); \
    if (_y == 0) _y = (ny); \
    if (_z == 0) _z = (nz); \
    if (_x == nx+1) _x = 1; \
    if (_y == ny+1) _y = 1; \
    if (_z == nz+1) _z = 1; \
    OUT = ((_x) + ((nx)+(NG*2))*((_y) + ((ny)+(NG*2))*(_z))); } \

#define RANK_TO_INDEX(rank,ix,iy,iz,_x,_y) \
    int _ix, _iy, _iz;                                                    \
    _ix  = (rank);                        /* ix = ix+gpx*( iy+gpy*iz ) */ \
    _iy  = _ix/int(_x);   /* iy = iy+gpy*iz */            \
    _ix -= _iy*int(_x);   /* ix = ix */                   \
    _iz  = _iy/int(_y);   /* iz = iz */                   \
    _iy -= _iz*int(_y);   /* iy = iy */                   \
    (ix) = _ix;                                                           \
    (iy) = _iy;                                                           \
    (iz) = _iz;                                                           \

enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

enum CellAcc {
    jfx = 0 * ACCUMULATOR_ARRAY_LENGTH,
    jfy = 1 * ACCUMULATOR_ARRAY_LENGTH,
    jfz = 2 * ACCUMULATOR_ARRAY_LENGTH,
};

enum CellInterp {
    ex = 0,
    dexdy,
    dexdz,
    d2exdydz,
    ey,
    deydz,
    deydx,
    d2eydzdx,
    ez,
    dezdx,
    dezdy,
    d2ezdxdy,
    cbx,
    dcbxdx,
    cby,
    dcbydy,
    cbz,
    dcbzdz,
};

#define NEIGHBOURS 12
#define FACES 6

enum CellMap {
    xd_yd_z = 0,    // USED
    xd_y_zd,        // USED
    xd_y_z,         // USED
    x_yd_zd,        // USED
    x_yd_z,         // USED
    x_y_zd,         // USED
    x_y_zu,         // USED
    x_yu_z,         // USED
    x_yu_zu,        // USED
    xu_y_z,         // USED
    xu_y_zu,        // USED
    xu_yu_z,        // USED
};

enum Face {
    xd = 0,   
    yd = 1,  
    zd = 2,
    xu = 3,
    yu = 4,
    zu = 5,
};

/***************************************************************************************************
 * @brief Utility class to temporarily hold the mesh data until it is loaded by OP-PIC
 */
class DataPointers
{
    public:
        DataPointers() {}
        virtual ~DataPointers()
        {
            DeleteValues();   
        }

        inline void DeleteValues() {

            if (c_index != nullptr) delete[] c_index;
            if (c_e != nullptr) delete[] c_e;
            if (c_b != nullptr) delete[] c_b;
            if (c_j != nullptr) delete[] c_j;
            if (c_acc != nullptr) delete[] c_acc;
            if (c_interp != nullptr) delete[] c_interp;
            if (c_pos_ll != nullptr) delete[] c_pos_ll;
            if (c_colors != nullptr) delete[] c_colors;
            if (c_ghost != nullptr) delete[] c_ghost;
            if (c_mask_right != nullptr) delete[] c_mask_right; 
            if (c_mask_ug != nullptr) delete[] c_mask_ug; 
            if (c_mask_ugb != nullptr) delete[] c_mask_ugb; 

            if (c2ngc_map != nullptr) delete[] c2ngc_map;
            if (c2c_map != nullptr) delete[] c2c_map;
            if (c2cug_map != nullptr) delete[] c2cug_map;
            if (c2cugb_map != nullptr) delete[] c2cugb_map;

            if (c2cug0_map != nullptr) delete[] c2cug0_map;
            if (c2cug1_map != nullptr) delete[] c2cug1_map;
            if (c2cug2_map != nullptr) delete[] c2cug2_map;
            if (c2cug3_map != nullptr) delete[] c2cug3_map;
            if (c2cug4_map != nullptr) delete[] c2cug4_map;
            if (c2cug5_map != nullptr) delete[] c2cug5_map;

            if (c2cugb0_map != nullptr) delete[] c2cugb0_map;
            if (c2cugb1_map != nullptr) delete[] c2cugb1_map;
            if (c2cugb2_map != nullptr) delete[] c2cugb2_map;
            if (c2cugb3_map != nullptr) delete[] c2cugb3_map;
            if (c2cugb4_map != nullptr) delete[] c2cugb4_map;
            if (c2cugb5_map != nullptr) delete[] c2cugb5_map;

            c_index = nullptr;
            c_e = nullptr;
            c_b = nullptr;
            c_j = nullptr;
            c_acc = nullptr;
            c_interp = nullptr;
            c_pos_ll = nullptr;
            c_colors = nullptr;
            c_ghost      = nullptr;
            c_mask_right = nullptr;
            c_mask_ug    = nullptr;
            c_mask_ugb   = nullptr;

            c2ngc_map = nullptr;
            c2c_map   = nullptr;
            c2cug_map   = nullptr;
            c2cugb_map  = nullptr;

            c2cug0_map  = nullptr;
            c2cug1_map  = nullptr;
            c2cug2_map  = nullptr;
            c2cug3_map  = nullptr;
            c2cug4_map  = nullptr;
            c2cug5_map  = nullptr;

            c2cugb0_map  = nullptr;
            c2cugb1_map  = nullptr;
            c2cugb2_map  = nullptr;
            c2cugb3_map  = nullptr;
            c2cugb4_map  = nullptr;
            c2cugb5_map  = nullptr;
        }

        inline void CreateMeshCommArrays() {

            this->c_index      = new OPP_INT[this->n_cells * ONE];
            this->c_pos_ll     = new OPP_REAL[this->n_cells * DIM];
            this->c_ghost      = new OPP_INT[this->n_cells * ONE];
            this->c_mask_right = new OPP_INT[this->n_cells * ONE];
            this->c_mask_ug    = new OPP_INT[this->n_cells * 6];
            this->c_mask_ugb   = new OPP_INT[this->n_cells * 6];

            this->c2c_map     = new OPP_INT[this->n_cells * NEIGHBOURS];
            this->c2ngc_map   = new OPP_INT[this->n_cells * FACES];
            this->c2cug_map   = new OPP_INT[this->n_cells * 6];
            this->c2cugb_map  = new OPP_INT[this->n_cells * 6];

            this->c2cug0_map  = new OPP_INT[this->n_cells];
            this->c2cug1_map  = new OPP_INT[this->n_cells];
            this->c2cug2_map  = new OPP_INT[this->n_cells];
            this->c2cug3_map  = new OPP_INT[this->n_cells];
            this->c2cug4_map  = new OPP_INT[this->n_cells];
            this->c2cug5_map  = new OPP_INT[this->n_cells];

            this->c2cugb0_map  = new OPP_INT[this->n_cells];
            this->c2cugb1_map  = new OPP_INT[this->n_cells];
            this->c2cugb2_map  = new OPP_INT[this->n_cells];
            this->c2cugb3_map  = new OPP_INT[this->n_cells];
            this->c2cugb4_map  = new OPP_INT[this->n_cells];
            this->c2cugb5_map  = new OPP_INT[this->n_cells];
        }

        inline void CreateMeshNonCommArrays(int cell_count) {

            this->c_e          = new OPP_REAL[cell_count * DIM];  
            this->c_b          = new OPP_REAL[cell_count * DIM];        
            this->c_j          = new OPP_REAL[cell_count * DIM];
            this->c_acc        = new OPP_REAL[cell_count * ACC_LEN]; 
            this->c_interp     = new OPP_REAL[cell_count * INTERP_LEN]; 
            this->c_colors    = new OPP_INT[cell_count];

            for (int n = 0; n < cell_count; n++) {

                for (int d = 0; d < DIM; d++) {
                    this->c_e[n * DIM + d] = 0.0;
                    this->c_b[n * DIM + d] = 0.0;
                    this->c_j[n * DIM + d] = 0.0;

                    for (int i = 0; i < ACCUMULATOR_ARRAY_LENGTH; i++)
                        this->c_acc[(n * DIM + d) * ACCUMULATOR_ARRAY_LENGTH + i] = 0.0;
                }

                for (int i = 0; i < INTERP_LEN; i++)
                    this->c_interp[n * INTERP_LEN + i] = 0.0;
                
                this->c_colors[n]    = 10000;
            }
        }

        inline void CreateMeshNonCommArrays() {
            CreateMeshNonCommArrays(this->n_cells);
        }
        int n_cells     = 0;
        int n_particles = 0;

        OPP_INT* c_index   = nullptr;
        OPP_REAL* c_e      = nullptr;
        OPP_REAL* c_b      = nullptr;
        OPP_REAL* c_j      = nullptr;
        OPP_REAL* c_acc    = nullptr;
        OPP_REAL* c_interp = nullptr;
        OPP_REAL* c_pos_ll = nullptr;
        OPP_INT* c_colors = nullptr;
        OPP_INT* c_ghost         = nullptr;
        OPP_INT* c_mask_right    = nullptr;
        OPP_INT* c_mask_ug       = nullptr;
        OPP_INT* c_mask_ugb      = nullptr;

        OPP_INT* c2ngc_map   = nullptr;
        OPP_INT* c2c_map     = nullptr;
        OPP_INT* c2cug_map   = nullptr;
        OPP_INT* c2cugb_map  = nullptr;
        
        OPP_INT* c2cug0_map  = nullptr;
        OPP_INT* c2cug1_map  = nullptr;
        OPP_INT* c2cug2_map  = nullptr;
        OPP_INT* c2cug3_map  = nullptr;
        OPP_INT* c2cug4_map  = nullptr;
        OPP_INT* c2cug5_map  = nullptr;

        OPP_INT* c2cugb0_map  = nullptr;
        OPP_INT* c2cugb1_map  = nullptr;
        OPP_INT* c2cugb2_map  = nullptr;
        OPP_INT* c2cugb3_map  = nullptr;
        OPP_INT* c2cugb4_map  = nullptr;
        OPP_INT* c2cugb5_map  = nullptr;
};

/*************************************************************************************************
 * Helper function for Deck constructor
*/
inline OPP_REAL courant_length(OPP_REAL lx, OPP_REAL ly, OPP_REAL lz, size_t nx, size_t ny, size_t nz ) {
    OPP_REAL w0, w1 = 0;
    if( nx>1 ) w0 = nx / lx, w1 += w0*w0;
    if( ny>1 ) w0 = ny / ly, w1 += w0*w0;
    if( nz>1 ) w0 = nz / lz, w1 += w0*w0;
    return sqrt(1 / w1);
}

/***************************************************************************************************
 * @brief Utility class to hold the simulation parameters similar to original
 */
class Deck {
    public:
    Deck(opp::Params& param) {
        num_steps = param.get<OPP_INT>("num_steps");
        v0 = param.get<OPP_REAL>("init_vel");
        nppc = param.get<OPP_INT>("num_part_per_cell");
        nx = param.get<OPP_INT>("nx");
        ny = param.get<OPP_INT>("ny");
        nz = param.get<OPP_INT>("nz");
        c = param.get<OPP_REAL>("c");
        qsp = param.get<OPP_REAL>("qsp");
        me = param.get<OPP_REAL>("me");

        const OPP_INT domain_expansion = param.get<OPP_INT>("domain_expansion");
        const OPP_REAL gam = 1.0 / sqrt(1.0 - v0 * v0);
        const OPP_REAL default_grid_len = param.get<OPP_REAL>("default_grid_len");
        len_x_global = default_grid_len;
        len_y_global = 0.628318530717959 * (gam * sqrt(gam));
        if (domain_expansion > 0) // Always expect to expand on z direction when weak scaling!   
            len_z_global = default_grid_len * domain_expansion;
        else
            len_z_global = default_grid_len * OPP_comm_size;    
        dx = len_x_global / nx;
        dy = len_y_global / ny;
        dz = len_z_global / nz;
        
        n0 = 2.0;
        dt = 0.99 * courant_length(len_x_global, len_y_global, len_z_global, nx, ny, nz) / c;
        
        Npe = n0 * len_x_global * len_y_global * len_z_global;
        Ne = (int64_t)nx * (int64_t)ny * (int64_t)nz * (int64_t)nppc;
        we = Npe / Ne;
        eps = param.get<OPP_REAL>("eps");

        cdt_dx = (c * dt / dx);
        cdt_dy = (c * dt / dy);
        cdt_dz = (c * dt / dz);
        qdt_2mc = qsp * dt / (2 * me * c);
        const OPP_REAL frac = 1.0;
        px = (nx>0) ? (frac * c * dt / dx) : 0;
        py = (ny>0) ? (frac * c * dt / dy) : 0;
        pz = (nz>0) ? (frac * c * dt / dz) : 0;
        acc_coefx = 0.25 / (dy * dz * dt);
        acc_coefy = 0.25 / (dz * dx * dt);
        acc_coefz = 0.25 / (dx * dy * dt);
        dt_eps0 = dt / eps;

        OPP_RUN_ON_ROOT() this->print();
    }

    inline void print() {
        opp_printf("Deck", "num_steps: %d ", num_steps);
        opp_printf("Deck", "v0: %2.25lE ", v0);
        opp_printf("Deck", "nppc: %d ", nppc);
        opp_printf("Deck", "nx: %d %d %d", nx, ny, nz);
        opp_printf("Deck", "len_x_global: %2.25lE %2.25lE %2.25lE", len_x_global, len_y_global, len_z_global);
        opp_printf("Deck", "dx: %2.25lE %2.25lE %2.25lE", dx, dy, dz);
        opp_printf("Deck", "n0: %2.25lE ", n0);
        opp_printf("Deck", "dt: %2.25lE ", dt);
        opp_printf("Deck", "c: %2.25lE ", c);
        opp_printf("Deck", "Npe: %2.25lE ", Npe);
        opp_printf("Deck", "Ne: %" PRId64 "", Ne);
        opp_printf("Deck", "we: %2.25lE ", we);
        opp_printf("Deck", "eps: %2.25lE ", eps);
        opp_printf("Deck", "qsp: %2.25lE ", qsp);
        opp_printf("Deck", "me: %2.25lE ", me);

        opp_printf("Deck", "cdt_dx: %2.25lE %2.25lE %2.25lE ", cdt_dx, cdt_dy, cdt_dz);
        opp_printf("Deck", "px: %2.25lE %2.25lE %2.25lE ", px, py, pz);
        opp_printf("Deck", "acc_coefx: %2.25lE %2.25lE %2.25lE ", acc_coefx, acc_coefy, acc_coefz);
        opp_printf("Deck", "qdt_2mc: %2.25lE ", qdt_2mc);
        opp_printf("Deck", "dt_eps0: %2.25lE ", dt_eps0);
    }

    OPP_INT num_steps = 0;
    OPP_REAL v0 = 0.0;
    OPP_INT nppc = 0;
    OPP_INT nx = 0;
    OPP_INT ny = 0;
    OPP_INT nz = 0;
    OPP_REAL len_x_global = 0;
    OPP_REAL len_y_global = 0;
    OPP_REAL len_z_global = 0;
    OPP_REAL dx = 0;
    OPP_REAL dy = 0;
    OPP_REAL dz = 0;
    OPP_REAL n0 = 0;
    OPP_REAL dt = 0;
    OPP_REAL c = 0;
    OPP_REAL Npe = 0;
    int64_t Ne = 0;
    OPP_REAL we = 0;
    OPP_REAL eps = 0;
    OPP_REAL qsp = 0;
    OPP_REAL me = 0;

    OPP_REAL cdt_dx = 0;
    OPP_REAL cdt_dy = 0;
    OPP_REAL cdt_dz = 0;
    OPP_REAL qdt_2mc = 0;
    OPP_REAL px = 0;
    OPP_REAL py = 0;
    OPP_REAL pz = 0;
    OPP_REAL acc_coefx = 0;
    OPP_REAL acc_coefy = 0;
    OPP_REAL acc_coefz = 0;
    OPP_REAL dt_eps0 = 0;
};

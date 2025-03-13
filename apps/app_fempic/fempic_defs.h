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

#include "opp_lib.h"
#ifdef USE_MPI
    #include "opp_mpi_core.h"
#endif
#include "opp_cluster.h"

#ifdef DEBUG_LOG
    #define FP_DEBUG true
#else
    #define FP_DEBUG false
#endif

#define INJ_EXCESS 100

#define ONE                1
#define DIM                3
#define N_PER_C            4
#define N_PER_IF           3
#define NEIGHB_C           4
#define DET_FIELDS         4
#define ALL_DET            (NEIGHB_C * DET_FIELDS)

const double EPS0 = 8.8541878e-12;      /*permittivity of free space*/
const double QE   = 1.602e-19;          /*elementary charge*/
const double AMU  = 1.660538921e-27;    /*atomic mass unit*/
const double Kb   = 8.617333262e-5;     /*Boltzmann's  constant*/

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

        inline void DeleteValues()
        {
            if (c_to_n) delete[] c_to_n;
            if (c_to_c) delete[] c_to_c;
            if (c_ef) delete[] c_ef;
            if (c_det) delete[] c_det;
            if (c_vol) delete[] c_vol;
            if (c_sd) delete[] c_sd;
            if (c_col) delete[] c_col;
            if (c_id) delete[] c_id;
            if (c_centroid) delete[] c_centroid;
            if (n_bnd_pot) delete[] n_bnd_pot;
            if (n_pot) delete[] n_pot;
            if (n_ion_den) delete[] n_ion_den;
            if (n_pos) delete[] n_pos;
            if (n_vol) delete[] n_vol;
            if (n_type) delete[] n_type;
            if (n_color) delete[] n_color;
            if (n_id) delete[] n_id;
            if (if_to_c) delete[] if_to_c;
            if (if_to_n) delete[] if_to_n;
            if (if_v_norm) delete[] if_v_norm;
            if (if_u_norm) delete[] if_u_norm;
            if (if_norm) delete[] if_norm;
            if (if_area) delete[] if_area;
            if (if_dist) delete[] if_dist;
            if (if_n_pos) delete[] if_n_pos;
            if (if_id) delete[] if_id;

            c_to_n = nullptr;
            c_to_c = nullptr;
            c_ef = nullptr;
            c_det = nullptr;
            c_vol = nullptr; 
            c_sd = nullptr;
            c_col = nullptr;
            c_id = nullptr;
            c_centroid = nullptr;

            n_bnd_pot = nullptr;
            n_pot = nullptr;
            n_ion_den = nullptr;
            n_pos = nullptr;
            n_vol = nullptr; 
            n_type = nullptr;
            n_color = nullptr;
            n_id = nullptr;

            if_to_c = nullptr;    
            if_to_n = nullptr;   
            if_v_norm = nullptr;
            if_u_norm = nullptr;
            if_norm = nullptr;  
            if_area = nullptr;    
            if_dist = nullptr;
            if_n_pos = nullptr;
            if_id = nullptr;
        }

        inline void CreateMeshArrays()
        {
            c_ef       = new double[n_cells * DIM];
            c_to_n     = new int[n_cells * N_PER_C];
            c_to_c     = new int[n_cells * NEIGHB_C];
            c_det      = new double[n_cells * ALL_DET]; // [alpha,beta,gamma,delta] * 4neighbours
            c_vol      = new double[n_cells];
            c_sd       = new double[n_cells * N_PER_C*DIM]; // arranged as [x,y,z] * 4 neighbours
            c_col      = new int[n_cells];
            c_id       = new int[n_cells];
            c_centroid = new double[n_cells * DIM];

            n_bnd_pot = new double[n_nodes];
            n_pot     = new double[n_nodes];
            n_ion_den = new double[n_nodes];
            n_pos     = new double[n_nodes * DIM];
            n_vol     = new double[n_nodes];
            n_type    = new int[n_nodes];
            n_color   = new int[n_nodes];
            n_id      = new int[n_nodes];

            if_to_c   = new int[n_ifaces];
            if_to_n   = new int[n_ifaces * N_PER_IF]; 
            if_v_norm = new double[n_ifaces * DIM]; 
            if_u_norm = new double[n_ifaces * DIM]; 
            if_norm   = new double[n_ifaces * DIM]; 
            if_area   = new double[n_ifaces]; 
            if_dist   = new int[n_ifaces];
            if_n_pos  = new double[n_ifaces * N_PER_IF * DIM];
            if_id     = new int[n_ifaces];  
        };

        int n_nodes  = 0;
        int n_cells  = 0;
        int n_ifaces = 0;

        int *c_to_n = nullptr;
        int *c_to_c = nullptr;
        double *c_ef = nullptr;
        double *c_det = nullptr;
        double *c_vol = nullptr; 
        double *c_sd = nullptr;
        int *c_col = nullptr; 
        int *c_id = nullptr;
        double *c_centroid = nullptr;

        double *n_bnd_pot = nullptr;
        double *n_pot = nullptr;
        double *n_ion_den = nullptr;
        double *n_pos = nullptr;
        double *n_vol = nullptr; 
        int *n_type = nullptr; 
        int *n_color = nullptr; 
        int *n_id = nullptr;

        int *if_to_c = nullptr;        // c_con
        int *if_to_n = nullptr;        // con[3]; 
        double *if_v_norm = nullptr;   // v[3]
        double *if_u_norm = nullptr;   // u[3]
        double *if_norm = nullptr;     // normal[3]
        double *if_area = nullptr;     // area
        int *if_dist = nullptr;
        double *if_n_pos = nullptr;
        int *if_id = nullptr;
};
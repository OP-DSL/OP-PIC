#pragma once

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

// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <string>
// #include <stdarg.h>
#include <opp_lib.h>

// #ifdef DEBUG_LOG
//     #define SP_DEBUG true
// #else
//     #define SP_DEBUG false
// #endif

// // global variables
// #define MPI_double MPI_DOUBLE
// #define STR_SIZE         501
// #define NG               100
// #define NP               10

// // definitions for diagnostics
// #define NDENSITY        0
// #define NEFIELD         1
// #define NVOLTAGE        2
// #define NVXPHASE        3
// #define NDIAGNOSTICS    4

// // variables for timing
// #define PARTICLES       0
// #define FIELDS          1
// #define MPITIME         2
// #define DIAGNOSTICS     3
// #define INIT            4
// #define NTIMES          5

// const double fourth         = 1.0 / 4.0;
// const double half           = 1.0 / 2.0;
// const double one            = 1.;
// const double one_third      = 1./3.;
// const double two_fifteenths = 2./15.;

extern double dtfactor;

extern double *pdata;
extern double *lhsbuf, *rhsbuf;
extern int npart, lpart, rpart;

// fields arrays
// extern double *Earray;
// extern double *phiarray;
// extern double *narray;
// extern double *nback;

extern double area;
extern double density;
extern double np2c;
extern double q, m, qm;
extern double qscale;
extern double epsilon;
extern double wp;
extern double El, Er;
extern double dx, dt, t;
extern double xl, xr, L, Ll;
extern int ntimesteps;

extern int lproc, rproc;
extern int nproc, rank;
extern int ver, subver;
extern int last;

extern int ng, nc, ngglobal;

extern int nl, nr;
extern int nlp, nrp;  // counters for sending particle buffers
extern int nlr, nrr;

extern double *tri_a, *tri_b, *tri_c, *gam; // TODO : Enrich properly

// extern bool diagnosticsflag;
extern double lhsvoltage;

// extern double tttimes[NTIMES];
// extern char names[NTIMES][STR_SIZE];
// extern double tstart, tend;  // global times
// extern double wtres;
extern double nparttot;

extern double tt;  // universal time

// variables for field communication
extern double *frbuffer;
extern double *fsbuffer;

// extern char diags[NDIAGNOSTICS][STR_SIZE];
// extern double * thediags[NDIAGNOSTICS];

// extern char procname[STR_SIZE];

//*************************************************************************************************
// global functions
// void Quit(void);

// inline double xcoord(int i);

// void gradient(double *grad, double *arr, int n, double scale);

enum Dir : char {
    Left = 0,
    Right,
};

void init(void);
void setcoeffs(double scale);

//*************************************************************************************************
void seq_field_solve_poissons_equation(opp_set set, opp_arg arg0, opp_arg arg1);
void seq_field_solve_get_potential_gradient(opp_set set, opp_arg arg0, opp_arg arg1);

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
            if (node_field_E_tmp) delete[] node_field_E_tmp;
            if (node_field_J_tmp) delete[] node_field_J_tmp;
            if (node_field_P_tmp) delete[] node_field_P_tmp;
            if (node_xlocal_tmp) delete[] node_xlocal_tmp;
            if (node_index_tmp) delete[] node_index_tmp;
            if (cell_to_nodes_tmp) delete[] cell_to_nodes_tmp;
            if (cell_to_cells_tmp) delete[] cell_to_cells_tmp;

            if (part_position_x_tmp) delete[] part_position_x_tmp;
            if (part_velocity_x_tmp) delete[] part_velocity_x_tmp;
            if (part_field_E_tmp) delete[] part_field_E_tmp;
            if (part_cell_index_tmp) delete[] part_cell_index_tmp;

            node_field_E_tmp = nullptr;
            node_field_J_tmp = nullptr;
            node_field_P_tmp = nullptr;
            node_xlocal_tmp = nullptr;
            node_index_tmp = nullptr;
            cell_to_nodes_tmp = nullptr;
            cell_to_cells_tmp = nullptr;

            part_position_x_tmp = nullptr;
            part_velocity_x_tmp = nullptr;
            part_field_E_tmp = nullptr;
            part_cell_index_tmp = nullptr;
        }

        int n_nodes;
        int n_cells;
        int n_particles;

        double *node_field_E_tmp = nullptr;
        double *node_field_J_tmp = nullptr;
        double *node_field_P_tmp = nullptr;
        double *node_xlocal_tmp = nullptr;
        int *node_index_tmp = nullptr;
        int *cell_to_nodes_tmp = nullptr;
        int *cell_to_cells_tmp = nullptr;

        double *part_position_x_tmp = nullptr;
        double *part_velocity_x_tmp = nullptr;
        double *part_field_E_tmp = nullptr;
        int *part_cell_index_tmp = nullptr;
};

inline std::unique_ptr<DataPointers> Load()
{ 
    ntimesteps = opp_params->get<OPP_INT>("num_steps");
    nc         = opp_params->get<OPP_INT>("ncpp");
    ng         = opp_params->get<OPP_INT>("ncpp") + 1;
    npart      = nc * opp_params->get<OPP_INT>("ppc");
    dtfactor   = opp_params->get<OPP_REAL>("dtfactor");
    lhsvoltage = opp_params->get<OPP_INT>("lhsv");;

    init();

    std::unique_ptr<DataPointers> m(new DataPointers());
    m->n_particles = npart;

    m->part_position_x_tmp = new double[npart];
    m->part_velocity_x_tmp = new double[npart];
    m->part_field_E_tmp    = new double[npart];
    m->part_cell_index_tmp = new int[npart];

    for (int i=0; i < npart; i++)
    {
        double position            = xl + ((i + 1)*(xr-xl))/double(npart + 1);
        m->part_position_x_tmp[i]  = position;
        m->part_velocity_x_tmp[i]  = 0.0;
        m->part_cell_index_tmp[i]  = ((position - xl)/dx);
        m->part_field_E_tmp[i]     = 0.0;
    }

    m->n_nodes = ng;
    m->n_cells = nc;

    m->node_field_E_tmp  = new double[ng];
    m->node_field_J_tmp  = new double[ng];
    m->node_field_P_tmp  = new double[ng];
    m->node_xlocal_tmp   = new double[ng];
    m->node_index_tmp    = new int[ng];
    m->cell_to_nodes_tmp = new int[nc * 2];
    m->cell_to_cells_tmp = new int[nc * 2];
    double xlocal = xl;

    for (int i=0; i < nc; i++)
    {
        m->cell_to_nodes_tmp[i * 2 + 0]   = i;
        m->cell_to_nodes_tmp[i * 2 + 1]   = (i+1);

        m->cell_to_cells_tmp[i * 2 + 0]   = (i-1);
        m->cell_to_cells_tmp[i * 2 + 1]   = (i+1);
    }
    m->cell_to_cells_tmp[0] = -1;              // outside the domain : Left
    m->cell_to_cells_tmp[(nc-1) * 2 + 1] = -1; // outside the domain : Right

    for (int i=0; i < ng; i++, xlocal +=dx)
    {
        m->node_field_E_tmp[i] = 0.0;      // Earray (Electric field)
        m->node_field_J_tmp[i] = 0.0;      // narray (Current density)
        m->node_field_P_tmp[i] = 0.0;      // phiarray (Potential)
        m->node_xlocal_tmp[i]  = xlocal;   // Local position of the node
        m->node_index_tmp[i]   = i;
    }

    {    
        tri_a = new double[ng];
        tri_b = new double[ng];
        tri_c = new double[ng];
        gam   = new double[ng];

        for(int i=0; i < ng; i++)
        {
            tri_a[i] = tri_b[i] = tri_c[i] = 0.0;
        }

        setcoeffs(-epsilon/(q*dx*dx));
    }

    return m;
} 
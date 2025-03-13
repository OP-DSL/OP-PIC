#pragma once

//*********************************************
// USER WRITTEN CODE
//*********************************************

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <stdarg.h>
#include <opp_lib.h>

#ifdef DEBUG_LOG
    #define SP_DEBUG true
#else
    #define SP_DEBUG false
#endif

// global variables
#define MPI_double MPI_DOUBLE
#define STR_SIZE         501
#define NG               100
#define NP               10

// definitions for diagnostics
#define NDENSITY        0
#define NEFIELD         1
#define NVOLTAGE        2
#define NVXPHASE        3
#define NDIAGNOSTICS    4

// variables for timing
#define PARTICLES       0
#define FIELDS          1
#define MPITIME         2
#define DIAGNOSTICS     3
#define INIT            4
#define NTIMES          5

const double fourth         = 1.0 / 4.0;
const double half           = 1.0 / 2.0;
const double one            = 1.;
const double one_third      = 1./3.;
const double two_fifteenths = 2./15.;

extern double dtfactor;

extern double *pdata;
extern double *lhsbuf, *rhsbuf;
extern int npart, lpart, rpart;

// fields arrays
extern double *Earray;
extern double *phiarray;
extern double *narray;
extern double *nback;

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

extern bool diagnosticsflag;
extern double lhsvoltage;

extern double tttimes[NTIMES];
extern char names[NTIMES][STR_SIZE];
extern double tstart, tend;  // global times
extern double wtres;
extern double nparttot;

extern double tt;  // universal time

// variables for field communication
extern double *frbuffer;
extern double *fsbuffer;

extern char diags[NDIAGNOSTICS][STR_SIZE];
extern double * thediags[NDIAGNOSTICS];

extern char procname[STR_SIZE];

//*************************************************************************************************
// global functions
void Quit(void);

inline double xcoord(int i);

void gradient(double *grad, double *arr, int n, double scale);

void init(void);

double rhsV(double t);

double lhsV(double t);

void setcoeffs(double scale);

//*************************************************************************************************
void seq_field_solve_poissons_equation(
    int set_size,
    double* node0_field_J,
    double* node0_field_P);

//*************************************************************************************************
void seq_field_solve_get_potential_gradient(
    int set_size,
    double* nodes_field_E,
    const double* nodes_field_P);

//*************************************************************************************************
void opp_par_loop_all__WeightFieldsToParticles(opp_set,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_par_loop_particle_all__PushParticles(opp_set,opp_arg,opp_arg,opp_arg,opp_arg);
void opp_par_loop_all__WeightParticlesToFields(opp_set,opp_arg,opp_arg,opp_arg);
void opp_par_loop_all__FieldSolveSumLaplace(opp_set,opp_arg,opp_arg);

//*************************************************************************************************


class Loader
{
    public:
        Loader() {}
        virtual ~Loader()
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

        double *part_position_x_tmp = nullptr;
        double *part_velocity_x_tmp = nullptr;
        double *part_field_E_tmp = nullptr;
        int *part_cell_index_tmp = nullptr;
};

inline std::unique_ptr<Loader> Load()
{ 
    ntimesteps = opp_params->get<OPP_INT>("nt");
    nc         = opp_params->get<OPP_INT>("ncpp");
    ng         = opp_params->get<OPP_INT>("ncpp") + 1;
    npart      = nc * opp_params->get<OPP_INT>("ppc");;
    dtfactor   = opp_params->get<OPP_REAL>("dtfactor");;
    lhsvoltage = opp_params->get<OPP_INT>("lhsv");;

    init();

    std::unique_ptr<Loader> loader(new Loader());
    loader->n_particles = npart;

    loader->part_position_x_tmp = new double[npart];
    loader->part_velocity_x_tmp = new double[npart];
    loader->part_field_E_tmp    = new double[npart];
    loader->part_cell_index_tmp = new int[npart];

    for (int i=0; i < npart; i++)
    {
        double position                 = xl + ((i + 1)*(xr-xl))/double(npart + 1);
        loader->part_position_x_tmp[i]  = position;
        loader->part_velocity_x_tmp[i]  = 0.0;
        loader->part_cell_index_tmp[i]  = ((position - xl)/dx);
        loader->part_field_E_tmp[i]     = 0.0;
    }

    loader->n_nodes = ng;
    loader->n_cells = nc;

    loader->node_field_E_tmp  = new double[ng];
    loader->node_field_J_tmp  = new double[ng];
    loader->node_field_P_tmp  = new double[ng];
    loader->node_xlocal_tmp   = new double[ng];
    loader->node_index_tmp    = new int[ng];
    loader->cell_to_nodes_tmp = new int[nc * 2];
    double xlocal = xl;

    for (int i=0, j=0; i < nc; i++, j+=2)
    {
        loader->cell_to_nodes_tmp[j]     = i;
        loader->cell_to_nodes_tmp[j+1]   = (i+1);
    }

    for (int i=0; i < ng; i++, xlocal +=dx)
    {
        loader->node_field_E_tmp[i] = 0.0;      // Earray (Electric field)
        loader->node_field_J_tmp[i] = 0.0;      // narray (Current density)
        loader->node_field_P_tmp[i] = 0.0;      // phiarray (Potential)
        loader->node_xlocal_tmp[i]  = xlocal;   // Local position of the node
        loader->node_index_tmp[i]   = i;
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

    return loader;
} 
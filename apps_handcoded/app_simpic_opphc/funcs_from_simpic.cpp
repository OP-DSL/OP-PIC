#include "simpic.h"
#include <cstring>

double dtfactor;

double *pdata;
double *lhsbuf, *rhsbuf;
int npart, lpart, rpart;

// fields arrays
double *Earray;
double *phiarray;
double *narray;
double *nback;

double area;
double density;
double np2c;
double q, m, qm;
double qscale;
double epsilon;
double wp;
double El, Er;
double dx, dt, t;
double xl, xr, L, Ll;
int ntimesteps;

int lproc, rproc;
int nproc, rank;
int ver, subver;
int last;

int ng, nc, ngglobal;

int nl, nr;
int nlp, nrp;  // counters for sending particle buffers
int nlr, nrr;

double *tri_a, *tri_b, *tri_c, *gam; // TODO : Enrich properly

bool diagnosticsflag;
double lhsvoltage;

double tttimes[NTIMES];
char names[NTIMES][STR_SIZE];
double tstart, tend;  // global times
double wtres;
double nparttot;

double tt;  // universal time

// variables for field communication
double *frbuffer;
double *fsbuffer;

char diags[NDIAGNOSTICS][STR_SIZE];
double * thediags[NDIAGNOSTICS];

char procname[STR_SIZE];


inline double xcoord(int i)
{
  return xl + i*dx;
}

void gradient(double *grad, double *arr, int n, double scale)
{
  int i;

  n--;
  for(i=1; i < n; i++)
    {
      grad[i] = scale*(arr[i+1] - arr[i-1]);
    }
}

void init(void)
{
    nparttot = 0.; // not used
    density = 1.E13;
    epsilon = 8.85E-12;
    area = 1.;
    L = 1.;
    q = 1.602E-19;
    m = 9.1E-31;

    rank = 0;
    nproc = 1;

    xl = rank/double(nproc);
    xr = (rank + 1.)/double(nproc);
    xl *= L;
    xr *= L;
    
    Ll = xr - xl;
    dx = Ll / double(nc);

    nl = rank - 1;
    nr = rank + 1;
    last = nproc - 1;
    if(rank == 0)
    {
        nl = last;
    }
    if(rank == last)
    {
        nr = 0;
    }

    np2c = density*area*Ll/double(npart);

    // calculate time step from plasma frequency
    wp = density*q*q/(epsilon*m);
    wp = sqrt(wp);

    dt = dtfactor/wp;

    qscale = np2c/(area*dx);

    t = ntimesteps*dt;
    qm = q*dt/m;

    diagnosticsflag = true;
}

double rhsV(double t)
{
    return 0.;
}

double lhsV(double t)
{
    return lhsvoltage;
}

void setcoeffs(double scale)
{
    if(rank == 0)
    {
        tri_a[0] = 0.0;
        tri_b[0] = 1.0;
        tri_c[0] = 0.0;
    }
    else
    {
        tri_b[0] = -scale*(1.0 + dx/xl);
        tri_c[0] = scale;
            
    }
    
    if(rank == last)
    {
        tri_a[nc] = 0.0;
        tri_b[nc] = 1.0;
        tri_c[nc] = 0.0;
    }
    else
    {
        tri_a[nc] = scale;
        tri_b[nc] = -scale*(1.0 + dx/(L - xr));
    }

    for(int i=1; i < nc; i++)
    {
        tri_a[i] = scale;
        tri_b[i] = (-2.0 * scale);
        tri_c[i] = scale;
    }
}



//*************************************************************************************************
void seq_field_solve_poissons_equation(
    int set_size,
    double* node0_field_J,
    double* node0_field_P)
{
    // modify density array
    double nlold, nrold;
    if(rank == 0)
    {
        nlold = node0_field_J[0];
        node0_field_J[0] = 0.;
    }
    else
    {
        node0_field_J[0] *= 2;
    }
    if(rank == last)
    {
        nrold = node0_field_J[nc];
        node0_field_J[nc] = 0.;
    }
    else
    {
        node0_field_J[nc] *= 2;
    }

    int nstrt = 0;

    // Tridiagonal matrix of Poisson equation 𝜙𝑗+1−2𝜙𝑗+𝜙𝑗−1=𝑝𝑗 is solved with Gaussian elimination by each processor
    {
        int j;
        double bet = tri_b[nstrt];
        node0_field_P[nstrt] = node0_field_J[nstrt]/bet;

        for(j = nstrt + 1; j < set_size; j++) 
        {
            gam[j]              = tri_c[j-1]/bet;
            bet                 = tri_b[j] - tri_a[j]*gam[j];
            node0_field_P[j]    = (node0_field_J[j] - tri_a[j]*node0_field_P[j-1])/bet;
        }

        for(j = nc - 1; j >= nstrt; j--)
        {
            node0_field_P[j] -= (gam[j+1] * node0_field_P[j+1]);
        }
    }

    // restore density array
    if(rank == 0)
    {
        node0_field_J[0] = 2*nlold;
    }
    if(rank == last)
    {
        node0_field_J[nc] = 2*nrold;
    }

}

//*************************************************************************************************
void seq_field_solve_get_potential_gradient(
    int set_size,
    double* nodes_field_E,
    const double* nodes_field_P)
{
    double scale = -0.5/dx;

    for(int i = 0; i < (set_size - 1); i++)
    {
        nodes_field_E[i] = scale * (nodes_field_P[i+1] - nodes_field_P[i-1]);
    }

    nodes_field_E[0]  = -(nodes_field_P[1] - nodes_field_P[0])/dx;
    nodes_field_E[nc] = -(nodes_field_P[nc] - nodes_field_P[nc-1])/dx;
}
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
#include "fempic_defs.h"
#include "minifempic_funcs.h"

#ifndef USE_PETSC
    #define Vec int
    #define Mat int
    #define KSP int
    #define KSPConvergedReason int
#endif

/*solver class*/
class FESolver {
public:
    enum Method { NonLinear, GaussSeidel, Lapack, Petsc };

    FESolver(opp_map c2n_map, opp_dat n_type, opp_dat n_pos, opp_dat n_bnd_pot,  
                int argc, char **argv);
    ~FESolver();

    void computePhi(opp_arg arg0, opp_arg arg1, opp_arg arg2);
    
    void preAssembly(opp_map c2n_map, opp_dat n_bnd_pot);
    void enrich_cell_shape_deriv(opp_dat cell_shape_deriv);

    bool linearSolve(double *ion_den); 
    void nonLinearSolve(double *ion_den); 
    void buildJmatrix();
    void buildF1Vector(double *ion_den);
    
    void summarize(std::ostream &out);  

protected:
    void addKe(std::map<int, std::map<int, double>>& sparse_K, int e, double ke[4][4]);
    void addFe(Vec *Fvec, int e, double fe[4]);
    double evalNa(int a, double xi, double eta, double zeta);
    void getNax(double nx[3], int e, int a);
    void initialzeMatrix(std::map<int, std::map<int, double>>& sparse_K);
    void computeNX(opp_dat n_pos, opp_map c2n_map);
    void sanityCheck();
    void init_node_to_eq_map(opp_dat n_type_dat);
    void initPetscStructures();

    Method fesolver_method = Method::Petsc;

    int *node_to_eq_map  = nullptr;        /*node_to_eq_map[n]=A*/
    double ***NX  = nullptr;                /*NX[e][a] is a dNa/dx [3] vector*/
    
    double *detJ = nullptr;     /* determinant of the jacobian x_xi */

    const opp_map c2n_map;

    const double n0 = 0.0;
    const double phi0 = 0.0;
    const double kTe = 0.0;
    const double wall_potential = 0.0;

    double *dLocal = nullptr;        /* d[neq] is the solution from the linear solver */  
    double *f1Local = nullptr;
    double *tempNEQ1 = nullptr, *tempNEQ2 = nullptr, *tempNEQ3 = nullptr;

    const int n_nodes_set = 0;
    const int n_nodes_inc_halo = 0; 
    const int n_elements_set = 0;
    const int n_elements_inc_halo = 0;

    int neq = 0;        /*number of unknowns/equations*/
    int global_neq = 0;
    int own_start = 0;
    int own_end = 0;

    /*quadrature points*/
    const double l[2] = { -sqrt(1.0/3.0), sqrt(1.0/3.0) };
    const double W[2] = { 1, 1 };
    const int n_int = 2;

    /* Petsc related variables */
    Vec         Bvec, F0vec, F1vec, Dvec, Yvec;       
    Mat         Jmat, Kmat;                                 
    KSP         ksp;            /* linear solver context */
    KSPConvergedReason reason;

    int *vecCol;                // in use - indices related to current MPI rank

// #ifdef USE_CUDA
    double *dLocal_d = nullptr;        /* d[neq] is the solution from the linear solver */  
    double *f1Local_d = nullptr;
    double *tempNEQ1_d = nullptr;
    double *tempNEQ2_d = nullptr;
    double *tempNEQ3_d = nullptr;
    double *detJ_d = nullptr; 

    int *node_to_eq_map_d= nullptr;
    
    int neqNBlocks = 0;
    int nodesNBlocks = 0;
    int nodes_inc_haloNBlocks = 0;
    int cells_inc_haloNBlocks = 0;
// #endif
};

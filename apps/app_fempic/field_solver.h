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

class FESolver {
public:
    FESolver(const opp_map c2n_map, const opp_dat n_type, const opp_dat n_pos, const opp_dat n_bnd_pot);

    ~FESolver();

    void compute_phi(opp_arg arg0, opp_arg arg1, opp_arg arg2);
    void enrich_cell_shape_deriv(opp_dat cell_shape_deriv);

protected:
    void init_f1_and_J(const opp_dat ion_den_dat);
    void build_f1_vector();
    void build_j_matrix();
    void compute_node_potential(const opp_dat n_bnd_pot_dat, opp_dat node_potential_dat);

    void pre_assembly(const opp_dat n_bnd_pot);
    void summarize(std::ostream &out);  

    void add_ke(std::map<int, std::map<int, double>>& sparse_K, int e, double ke[4][4]);
    void add_fe(Vec *Fvec, int e, double fe[4]);
    double evaluate_na(int a, double xi, double eta, double zeta);
    void get_nax(double nx[3], int e, int a);
    void initialze_matrix(std::map<int, std::map<int, double>>& sparse_K);
    void compute_nx(const opp_dat n_pos);
    void sanity_check();
    void init_node_to_eq_map(const opp_dat n_type_dat);
    void init_petsc_structures();
    void calculate_neq(const opp_dat n_type_dat);
    void init_device_variables();
    void destroy_device_variables();

    std::vector<int> node_to_eq_map;
    std::vector<std::array<std::array<double, 3>, 4>> NX;
    std::vector<double> detJ;     /* determinant of the jacobian x_xi */

    const opp_map c2n_map;

    const double CONST_n0 = 0.0;
    const double CONST_phi0 = 0.0;
    const double CONST_kTe = 0.0;
    const double CONST_EPS0 = 0.0; 
    const double CONST_QE = 0.0;     

    std::vector<double> dLocal;        /* d[neq] is the solution from the linear solver */  
    std::vector<double> f1Local;
    std::vector<double> tempNEQ1, tempNEQ2, tempNEQ3;

    const int n_nodes_set = 0;
    const int n_nodes_inc_halo = 0; 
    const int n_cells_set = 0;
    const int n_cells_inc_halo = 0;

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

    std::vector<int> vec_col;                // in use - indices related to current MPI rank

    double *dLocal_d = nullptr;
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

    double* l_DEV_CONST = nullptr;
    double* W_DEV_CONST = nullptr;
};

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

#include "../field_solver.h"

#define is_neq(a, b) ((a) != (b))
#define SCALLING 100

#if defined(USE_CUDA) || defined(USE_HIP) || defined(USE_SYCL)
#include "device.cpp"
#else
#include "host.cpp"
#endif

// #ifdef USE_PETSC
//*************************************************************************************************
FESolver::FESolver(
    const opp_map c2n_map, 
    const opp_dat n_type_dat, 
    const opp_dat n_pos_dat,
    const opp_dat n_bnd_pot_dat) :
        c2n_map(c2n_map),
        CONST_n0(opp_params->get<OPP_REAL>("plasma_den")),
        CONST_kTe(Kb * opp_params->get<OPP_REAL>("electron_temperature")),
        CONST_EPS0(EPS0),
        CONST_QE(QE),
        n_nodes_set(n_type_dat->set->size),
        n_nodes_inc_halo(n_type_dat->set->size + n_type_dat->set->exec_size + n_type_dat->set->nonexec_size),
        n_cells_set(c2n_map->from->size),
        n_cells_inc_halo(c2n_map->from->size + c2n_map->from->exec_size)
{
    if (OPP_DBG) opp_printf("FESolver", "FESolver");

    // opp_printf("FESolver", "n_nodes_set %d | n_nodes_inc_halo %d | n_cells_set %d | n_cells_inc_halo %d", 
    //     n_nodes_set, n_nodes_inc_halo, n_cells_set, n_cells_inc_halo);

    calculate_neq(n_type_dat);

    node_to_eq_map.resize(n_nodes_inc_halo);

    NX.resize(n_cells_inc_halo);
    detJ.resize(n_cells_inc_halo);
    
    dLocal.resize(neq, 0.0);
    f1Local.resize(neq, 0.0);
    tempNEQ1.resize(neq);
    tempNEQ2.resize(neq);
    tempNEQ3.resize(neq);

    init_petsc_structures();

    init_node_to_eq_map(n_type_dat);

    compute_nx(n_pos_dat);

    summarize(std::cout);

    pre_assembly(n_bnd_pot_dat);

    sanity_check();

    init_device_variables();
}

//*************************************************************************************************
FESolver::~FESolver() 
{    
    KSPDestroy(&ksp);
    VecDestroy(&Bvec);
    VecDestroy(&F0vec);
    VecDestroy(&F1vec);
    VecDestroy(&Dvec);
    VecDestroy(&Yvec);
    MatDestroy(&Jmat);
    MatDestroy(&Kmat);

    destroy_device_variables();
}

//*************************************************************************************************
/*
    wrapper for solving the Poisson's equation
    build F1 using ion_density and D
    B = K * D
    B = B - F0
    B = B - F1
    Build J using K and D
    Ksp Solve : J * Y = B (Y = solution)
    D = D - Y
*/
void FESolver::compute_phi(opp_arg arg0, opp_arg arg1, opp_arg arg2) 
{ 
    if (OPP_DBG) opp_printf("FESolver", "compute_phi START");

    opp_profiler->start("compute_phi");

    const int nargs = 3;
    opp_arg args[nargs];
    args[0] = arg0;  // node_potential
    args[1] = arg1;  // ion_den
    args[2] = arg2;  // n_bnd_pot
    args[1].idx = 2; // HACK to forcefully make halos to download

    opp_profiler->start("FSolve_halo");
    opp_mpi_halo_exchanges_grouped(arg1.dat->set, nargs, args, DEVICE_TYPE);
    opp_mpi_halo_wait_all(nargs, args);
    opp_profiler->end("FSolve_halo");

    opp_profiler->start("FSolve_solv_init");
    init_f1_and_J(args[1].dat); // ion_den
    opp_profiler->start("FSolve_solv_init");

    opp_profiler->start("FSolve_f1");
    build_f1_vector();
    opp_profiler->end("FSolve_f1");

    MatCopy(Kmat, Jmat, DIFFERENT_NONZERO_PATTERN);

    opp_profiler->start("FSolve_j");
    build_j_matrix();           // Build J using K and D
    opp_profiler->end("FSolve_j");

    MatMult(Kmat, Dvec, Bvec);  // B = K * D
    VecAXPY(Bvec, -1.0, F0vec); // B = B - F0
    VecAXPY(Bvec, -1.0, F1vec); // B = B - F1

    opp_profiler->start("FSolve_ksp");
    KSPSolve(ksp, Bvec, Yvec); // Jmat * Yvec = Bvec (Yvec = solution)
    opp_profiler->end("FSolve_ksp");

    VecAXPY(Dvec, -1.0, Yvec); // D = D - Y

    KSPGetConvergedReason(ksp, &reason);
    if (reason < 0) {
        char* str_reason;
        KSPGetConvergedReasonString(ksp, (const char**)(&str_reason));
        std::cerr << "linear_solve Petsc failed to converge : " << str_reason <<
            "; run with -ksp_converged_reason" << std::endl;           
    }

    opp_profiler->start("FSolve_npot");
    compute_node_potential(args[2].dat, args[0].dat); // (n_bnd_pot, node_potential)
    opp_profiler->end("FSolve_npot");

    opp_set_dirtybit_grouped(nargs, args, DEVICE_TYPE);

    opp_profiler->end("compute_phi");

    if (OPP_DBG) opp_printf("FESolver", "compute_phi DONE");
}

//*************************************************************************************************
void FESolver::calculate_neq(const opp_dat n_type_dat) 
{
    // TODO : DO NOT CALCULATE SOLUTION FOR IMPORT NON EXEC, the owning rank will do that
    for (int i = 0; i < n_nodes_set; i++) 
        if (opp_get_data<OPP_INT>(n_type_dat)[i] == NORMAL || opp_get_data<OPP_INT>(n_type_dat)[i] == OPEN) 
            neq++;
}

//*************************************************************************************************
void FESolver::init_node_to_eq_map(const opp_dat n_type_dat) // relation from node indices to equation indices
{
    int P=0;
    for (int n = 0; n < n_nodes_inc_halo; n++)
    {
        if (n < n_nodes_set && 
            (opp_get_data<OPP_INT>(n_type_dat)[n] == NORMAL || opp_get_data<OPP_INT>(n_type_dat)[n] == OPEN))
        {
            node_to_eq_map[n]=P;
            P++;
        }
        else 
        {
            node_to_eq_map[n]=INT_MIN;    /* dirichlet node or for halo region, until updated by MPI comm */
        }
    }

#ifdef USE_MPI
    // send all the node_to_eq_map values of the export non-exec halos
    // (send the values by adding by own_start and making the values negative)
    // receive all node_to_eq_map values of the import non exec halos and 
    // assign the values to the correct indices of node_to_eq_map array

    halo_list exp_nonexec_list = OPP_export_nonexec_list[n_type_dat->set->index];
    halo_list imp_nonexec_list = OPP_import_nonexec_list[n_type_dat->set->index];

    std::vector<int> send_buffer(exp_nonexec_list->size);
    std::vector<MPI_Request> send_req(exp_nonexec_list->ranks_size);

    for (int i = 0; i < exp_nonexec_list->ranks_size; i++) 
    {
        for (int j = 0; j < exp_nonexec_list->sizes[i]; j++) 
        {
            int index = exp_nonexec_list->disps[i] + j;           
            send_buffer[index] = node_to_eq_map[exp_nonexec_list->list[index]];

            if (send_buffer[index] != INT_MIN)
                send_buffer[index] = -1 * (send_buffer[index] + own_start);

            // opp_printf("SEND", "idx=%d val=%d to rank %d", exp_nonexec_list->list[index], 
            //     send_buffer[index], exp_nonexec_list->ranks[i]);
        }

        MPI_Isend(&send_buffer[exp_nonexec_list->disps[i]], exp_nonexec_list->sizes[i], MPI_INT,
                exp_nonexec_list->ranks[i], 11, OPP_MPI_WORLD, &(send_req[i]));
    }

    const int nonexec_init = n_type_dat->set->size + n_type_dat->set->exec_size;
    for (int i = 0; i < imp_nonexec_list->ranks_size; i++) 
    {
        // opp_printf("SEND", "count=%d from rank %d", imp_nonexec_list->sizes[i], imp_nonexec_list->ranks[i]);

        MPI_Recv(&(node_to_eq_map[nonexec_init + imp_nonexec_list->disps[i]]), imp_nonexec_list->sizes[i], MPI_INT,
            imp_nonexec_list->ranks[i], 11, OPP_MPI_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Waitall(exp_nonexec_list->ranks_size, send_req.data(), MPI_STATUSES_IGNORE);
#endif
}

//*************************************************************************************************
void FESolver::init_petsc_structures()
{
    MatCreate(PETSC_COMM_WORLD, &Jmat);
    MatCreate(PETSC_COMM_WORLD, &Kmat);

    MatSetSizes(Jmat, neq, neq, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetSizes(Kmat, neq, neq, PETSC_DETERMINE, PETSC_DETERMINE);

#ifndef USE_MPI
    MatSetType(Jmat, MATSEQAIJ);
    MatSetType(Kmat, MATSEQAIJ);

    VecCreateSeq(PETSC_COMM_WORLD, neq, &Bvec);
#else
    MatSetType(Jmat, MATMPIAIJ);
    MatSetType(Kmat, MATMPIAIJ);

    VecCreateMPI(PETSC_COMM_WORLD, neq, PETSC_DETERMINE, &Bvec);
#endif

    MatSetFromOptions(Jmat);
    MatSetFromOptions(Kmat);
    
    MatSetUp(Jmat);
    MatSetUp(Kmat);

    VecSetFromOptions(Bvec);
    
    VecDuplicate(Bvec, &F0vec);
    VecDuplicate(Bvec, &F1vec);
    VecDuplicate(Bvec, &Dvec);
    VecDuplicate(Bvec, &Yvec);

    VecSet(Bvec, 0.0);
    VecSet(F0vec, 0.0);
    VecSet(F1vec, 0.0);
    VecSet(Dvec, 0.0);
    VecSet(Yvec, 0.0);
    VecAssemblyBegin(Bvec); VecAssemblyEnd(Bvec);
    VecAssemblyBegin(F0vec); VecAssemblyEnd(F0vec);
    VecAssemblyBegin(F1vec); VecAssemblyEnd(F1vec);
    VecAssemblyBegin(Dvec); VecAssemblyEnd(Dvec);
    VecAssemblyBegin(Yvec); VecAssemblyEnd(Yvec);

    VecGetOwnershipRange(Bvec, &own_start, &own_end);
    VecGetSize(Bvec, &global_neq);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetOperators(ksp, Jmat, Jmat);
    // KSPSetTolerances(ksp, 1.e-2 / (neq * neq), 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT); 
    // KSPSetTolerances(ksp, 1.e-100 / (neq * neq), 1.e-100, PETSC_DEFAULT, PETSC_DEFAULT); 
    KSPSetTolerances(ksp, 1.e-100, 1.e-100, PETSC_DEFAULT , PETSC_DEFAULT); 
    KSPSetFromOptions(ksp); 

    //PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

    vec_col.resize(neq);
    for (int i = 0; i < neq; i++) 
        vec_col[i] = (i + own_start);
}

//*************************************************************************************************
/* preassembles the K matrix and "F0 vector */
void FESolver::pre_assembly(const opp_dat n_bnd_pot) 
{
    if (OPP_DBG)
        opp_printf("FESolver", "pre_assembly global_neq %d neq %d n_cells_inc_halo %d", 
            global_neq, neq, n_cells_inc_halo);

    std::map<int, std::map<int, double>> sparse_K;

    for (int e = 0; e < n_cells_inc_halo; e++) 
    {
        double ke[4][4];

        for (int a=0;a<4;a++)
            for (int b=0;b<4;b++) 
            {
                ke[a][b] = 0;    /*reset*/

                /*perform quadrature*/
                for (int k=0;k<n_int;k++)
                    for (int j=0;j<n_int;j++)
                        for (int i=0;i<n_int;i++) 
                        {
                            double nax[3],nbx[3];

                            get_nax(nax,e,a);
                            get_nax(nbx,e,b);

                            /*dot product*/
                            double dot=0;
                            for (int d=0;d<3;d++) 
                                dot+=nax[d]*nbx[d];
                            ke[a][b] += dot*detJ[e]*W[i]*W[j]*W[k];
                        }
            }

        add_ke(sparse_K, e,ke);
   
        double fe[4]; /*force vector*/

        for (int a=0;a<4;a++) 
        {  
            double fh=0;    /*second term int(na*h), always zero since support only h=0*/
            double fg = 0;  /*third term, -sum(kab*qb)*/

            for (int b=0;b<4;b++) 
            {
                // This can reach import halos
                int node_idx = c2n_map->map[e * c2n_map->dim + b]; 
                const double gb = opp_get_data<OPP_REAL>(n_bnd_pot)[node_idx];
                fg-=ke[a][b]*gb;
            }

            fe[a] = fh + fg; /*combine*/
        }

        add_fe(&F0vec, e,fe);

    }  /*end of element*/

    VecAssemblyBegin(F0vec); VecAssemblyEnd(F0vec);

    initialze_matrix(sparse_K);

    const double scalar2 = 1 / (SCALLING * SCALLING); 

    VecScale(F0vec, scalar2);   // downscalling since NX is scalled to avoid precision issues
    MatScale(Kmat, scalar2);    // downscalling since NX is scalled to avoid precision issues

    sparse_K.clear();
}

//*************************************************************************************************
void FESolver::initialze_matrix(std::map<int, std::map<int, double>>& sparse_K)
{
    if (OPP_DBG) opp_printf("FESolver", "initialze_matrix START");

    std::vector<int> matCol(neq);                         // Number of non zero columns per row
    std::vector<std::vector<int>> matIndex(neq);          // Non zero column indices per row
    std::vector<std::vector<double>> matIndexValues(neq); // Non zero column values per row (temp copy to enrich K matrix)

    int diag_max_fields = 0, off_diag_max_fields = 0;
    
    for (int i=own_start; i<own_end; i++) { // iterate own rows
    
        std::vector<int> tempVec;
        int local_diag_max_fields = 0, local_off_diag_max_fields = 0;
        
        std::map<int, double>& K_sparse_row = sparse_K[i-own_start];

        for (int j=0; j<global_neq; j++) { // iterate all columns
            
            double value = 0.0;
            auto it = K_sparse_row.find(j);
            if (it != K_sparse_row.end())
                value = it->second;
       
            if ((std::abs(value) > 1e-12) || (i == j)) { // significant ones and all diagonals  
                tempVec.push_back(j); 
                if (j >= own_start && j < own_end)
                    local_diag_max_fields++;
                else
                    local_off_diag_max_fields++;
            }         
        }

        matCol[i-own_start] = tempVec.size();
        matIndex[i-own_start].resize(tempVec.size());
        matIndexValues[i-own_start].resize(tempVec.size());

        std::copy(tempVec.begin(), tempVec.end(), matIndex[i-own_start].data());

        diag_max_fields = (diag_max_fields > local_diag_max_fields) ? 
                                diag_max_fields : local_diag_max_fields; 

        off_diag_max_fields = (off_diag_max_fields > local_off_diag_max_fields) ? 
                                off_diag_max_fields : local_off_diag_max_fields;  
    }

#ifndef USE_MPI
    MatSeqAIJSetPreallocation(Kmat, diag_max_fields, nullptr);
    MatSeqAIJSetPreallocation(Jmat, diag_max_fields, nullptr);
#else
    MatMPIAIJSetPreallocation(Kmat, diag_max_fields, NULL, off_diag_max_fields, NULL);
    MatMPIAIJSetPreallocation(Jmat, diag_max_fields, NULL, off_diag_max_fields, NULL);
#endif

    if (OPP_DBG) opp_printf("FESolver", "initialze_matrix diag_max_fields=%d off_diag_max_fields=%d", 
        diag_max_fields, off_diag_max_fields);

    for (int i=own_start; i<own_end; i++) {
    
        for (int j=0; j<matCol[i-own_start]; j++) 
            matIndexValues[i-own_start][j] = sparse_K[i-own_start][matIndex[i-own_start][j]];

        MatSetValues(Kmat, 1, &i, matCol[i-own_start], matIndex[i-own_start].data(), 
            matIndexValues[i-own_start].data(), INSERT_VALUES); 
    }

    MatAssemblyBegin(Kmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Kmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(Jmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Jmat, MAT_FINAL_ASSEMBLY);

    if (OPP_DBG) opp_printf("FESolver", "initialze_matrix END");
}

//*************************************************************************************************

// UTIL FUNCTIONS

//*************************************************************************************************
void FESolver::summarize(std::ostream &out) 
{
    if (OPP_DBG)
        opp_printf("FESolver", "own neq is %d, the global neq is %d, start_eq is %d, end_eq %d -------", 
            neq, global_neq, own_start, own_end);
}

//*************************************************************************************************
void FESolver::enrich_cell_shape_deriv(opp_dat cell_shape_deriv)
{
    if (OPP_DBG) opp_printf("FESolver", "enrich_cell_shape_deriv START");

    opp_profiler->start("EnrichCellShapeDeriv");

    // copying only up to set size, hence import exec halos will be zero
    for (int cellID = 0; cellID < n_cells_inc_halo; cellID++)
    {
        for (int nodeCon = 0; nodeCon < N_PER_C; nodeCon++)
        {
            opp_get_data<OPP_REAL>(cell_shape_deriv)[cellID * (N_PER_C*DIM) + nodeCon * DIM + 0 ] = 
                (NX[cellID][nodeCon][0] / SCALLING);
            
            opp_get_data<OPP_REAL>(cell_shape_deriv)[cellID * (N_PER_C*DIM) + nodeCon * DIM + 1 ] = 
                (NX[cellID][nodeCon][1] / SCALLING);
            
            opp_get_data<OPP_REAL>(cell_shape_deriv)[cellID * (N_PER_C*DIM) + nodeCon * DIM + 2 ] = 
                (NX[cellID][nodeCon][2] / SCALLING);
        }
    }

#ifdef USE_MPI
    cell_shape_deriv->dirtybit = 1;
    opp_arg arg0 = opp_arg_dat(cell_shape_deriv, OPP_RW);
    arg0.idx = 2; // HACK
    opp_mpi_halo_exchanges_grouped(cell_shape_deriv->set, 1, &arg0, Device_CPU);  
    opp_mpi_halo_wait_all(1, &arg0);
#endif

    cell_shape_deriv->dirty_hd = Dirty::Device; // To make GPU versions to download updated data!

    opp_profiler->end("EnrichCellShapeDeriv");

    if (OPP_DBG) opp_printf("FESolver", "enrich_cell_shape_deriv END");
}

/*adds contributions from element stiffness matrix*/
void FESolver::add_ke(std::map<int, std::map<int, double>>& sparse_K, int e, double ke[4][4]) // BUG : K is not created correctly
{ 
    if (neq <= 0) return;

    for (int a=0;a<4;a++)         /*tetrahedra*/
    {
        for (int b=0;b<4;b++) 
        {
            const int node_idx1 = c2n_map->map[e * c2n_map->dim + a];
            const int P = node_to_eq_map[node_idx1]; 

            if (P<0) continue;    /* skip g nodes or not in current ranks own row range */

            const int node_idx2 = c2n_map->map[e * c2n_map->dim + b];
            int Q = node_to_eq_map[node_idx2]; 

            if (Q<0)
            {
                if (Q == INT_MIN) continue; // skip g nodes
                else Q = -1 * Q;            // on a different MPI rank, get its correct column
            }
            else 
                Q += own_start;             // on the same MPI rank, get its correct column

            std::map<int, double>& K_sparse_row = sparse_K[P];
            auto it = K_sparse_row.find(Q);
            if (it == K_sparse_row.end()) {
                K_sparse_row.insert(std::pair<int,double>(Q, ke[a][b]));
            }
            else {
                it->second += ke[a][b];
            }
        }
    }
}

/*adds contributions from element force vector to a global F vector*/
void FESolver::add_fe(Vec *Fvec, int e, double fe[4]) 
{ 
    if (neq <= 0) return;
    
    for (int a=0;a<4;a++)    /*tetrahedra*/ 
    {
        const int node_idx = c2n_map->map[e * c2n_map->dim + a];
        const int P = node_to_eq_map[node_idx]; 
        if (P<0) continue;    /* skip g nodes or on a different rank (owner compute)*/
        
        VecSetValue(*Fvec, P + own_start, fe[a], ADD_VALUES); // F[P] += fe[a];
    }
}

/*evaluates shape function a at position (xi,eta,zeta)*/
double FESolver::evaluate_na(int a, double xi, double eta, double zeta) 
{
    switch(a) 
    {
        case 0: return xi;
        case 1: return eta;
        case 2: return zeta;
        case 3: return 1-xi-eta-zeta;       
    }

    return 0;    /*shouldn't happen*/
}

/*returns derivative of N[a] at some logical point since we are using linear elements, 
these are constant in each element*/
void FESolver::get_nax(double nx[3], int e, int a) 
{
    for (int d=0;d<3;d++)
        nx[d] = NX[e][a][d];
}

/*computes derivatives of the shape functions for all elements constants since 
using linear elements*/
void FESolver::compute_nx(const opp_dat n_pos) 
{ 
    /*derivatives of the shape functions vs. xi*/
    double na_xi[4][3] = {{1,0,0}, {0,1,0}, {0,0,1}, {-1,-1,-1}};

    for (int e = 0; e < n_cells_inc_halo; e++) 
    {       
        /*node indices*/
        int* map0idx = &((int*)c2n_map->map)[e * c2n_map->dim]; 

        double x[4][3];
        for (int a=0;a<4;a++) 
        {
            /*node positions*/
            const double *pos = opp_get_data<OPP_REAL>(n_pos) + (map0idx[a] * n_pos->dim); 
            for (int d=0;d<3;d++) 
                x[a][d] = pos[d];
        }

        double x_xi[3][3];

        for (int i=0;i<3;i++)    /*x/y/z*/
            for (int j=0;j<3;j++) /*xi/eta/zeta*/ 
            {
                x_xi[i][j] = 0;
                for (int a=0; a<4; a++)    /*tet node*/
                    x_xi[i][j] += na_xi[a][j]*x[a][i];
            }

        detJ[e] = det3(x_xi);

        double xi_x[3][3];
        inverse(x_xi,xi_x);
      
        for (int a=0;a<4;a++) /*evaluate na_x*/
            for (int d=0;d<3;d++) 
            {
                NX[e][a][d]=0;
                for (int k=0;k<3;k++)
                    // scalled to avoid precision issues
                    NX[e][a][d] += (na_xi[a][k]*xi_x[k][d] * SCALLING); 
            }
    }
}

//*************************************************************************************************
// Sanity check after assembling matrices
void FESolver::sanity_check()
{
    int j_gm, j_gn, k_gm, k_gn, j_row_start, j_row_end, k_row_start, k_row_end;

    MatGetSize(Jmat, &j_gm, &j_gn);
    MatGetSize(Kmat, &k_gm, &k_gn);

    MatGetOwnershipRange(Jmat, &j_row_start, &j_row_end);
    MatGetOwnershipRange(Kmat, &k_row_start, &k_row_end);

    if (is_neq(j_gm, j_gn) || is_neq(k_gm, k_gn) || is_neq(j_gm, global_neq) ||
        is_neq(j_row_start, k_row_start) || is_neq(j_row_end, k_row_end) || 
        is_neq(j_row_start, own_start) || is_neq(j_row_end, own_end))
    {
        opp_printf("FESolver::FESolver", "Error... Matrix vector sizes issue");
    }      
    
    if (OPP_DBG) 
        opp_printf("FESolver::FESolver", 
            "j_gm %d, j_gn %d, k_gm %d, k_gn %d, j_row_start %d, j_row_end %d, k_row_start %d, k_row_end %d", 
            j_gm, j_gn, k_gm, k_gn, j_row_start, j_row_end, k_row_start, k_row_end);
}

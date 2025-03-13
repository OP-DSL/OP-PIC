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

//*************************************************************************************************
FESolver::FESolver(
    const opp_map c2n_map, const opp_dat n_type_dat, const opp_dat n_pos_dat, const opp_dat n_bnd_pot_dat) :
        c2n_map(c2n_map),
        CONST_n0(opp_params->get<OPP_REAL>("plasma_den")),
        CONST_kTe(Kb * opp_params->get<OPP_REAL>("electron_temperature")),
        CONST_EPS0(EPS0),
        CONST_QE(QE)
{
    if (OPP_IS_VALID_PROCESS) {
        n_nodes_set = n_type_dat->set->size;
        n_nodes_inc_halo = n_type_dat->set->size + n_type_dat->set->exec_size + n_type_dat->set->nonexec_size;
        n_cells_set = c2n_map->from->size;
        n_cells_inc_halo = c2n_map->from->size + c2n_map->from->exec_size;

        // if (OPP_DBG) 
        OPP_RUN_ON_ROOT()
            opp_printf("FESolver", "n_nodes_set %d | n_nodes_inc_halo %d | n_cells_set %d | n_cells_inc_halo %d", 
                n_nodes_set, n_nodes_inc_halo, n_cells_set, n_cells_inc_halo);

        // this is required since n_bnd_pot_dat is updated using an opp_par_loop
        opp_mpi_force_halo_update_if_dirty(n_bnd_pot_dat->set, { n_bnd_pot_dat }, Device_CPU); 

        calculate_neq(n_type_dat);

        node_to_eq_map.resize(n_nodes_inc_halo);

        NX.resize(n_cells_inc_halo);
        detJ.resize(n_cells_inc_halo);
        
        dLocal.resize(neq, 0.0);
        f1Local.resize(neq, 0.0);
        tempNEQ1.resize(neq, 0.0);
        tempNEQ2.resize(neq, 0.0);
        tempNEQ3.resize(neq, 0.0);
    }

#ifdef USE_MPI
    MPI_Allreduce(&neq, &global_neq, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
#else
    global_neq = neq;
#endif

    init_petsc_structures();
    init_node_to_eq_map(n_type_dat);
    compute_nx(n_pos_dat);
    pre_assembly(n_bnd_pot_dat);
    sanity_check();
    init_device_variables();

    OPP_RUN_ON_ROOT()
            opp_printf("FESolver", "Init done");
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

/**************************************************************************************************
    wrapper for solving the Poisson's equation
    (1) build F1 using ion_density and D
    (2) B = K * D
    (3) B = B - F0
    (4) B = B - F1
    (5) Build J using K and D
    (6) Ksp Solve : J * Y = B (Y = solution)
    (7) D = D - Y
*/
void FESolver::compute_phi(opp_arg arg0, opp_arg arg1, opp_arg arg2) 
{ 
    if (OPP_DBG) opp_printf("FESolver", "compute_phi START");

    opp_profiler->start("compute_phi");

    // No need for halo exchanges since we only need own node ion densities

    opp_profiler->start("FSolve_solv_init");
    init_f1_and_J(arg1.dat); // ion_den
    opp_profiler->end("FSolve_solv_init");

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
    compute_node_potential(arg2.dat, arg0.dat); // (n_bnd_pot, node_potential)
    opp_profiler->end("FSolve_npot");

    opp_set_dirtybit_grouped(1, &arg0, DEVICE_TYPE);

    opp_profiler->end("compute_phi");

    if (OPP_DBG) opp_printf("FESolver", "compute_phi DONE");
}

//*************************************************************************************************
void FESolver::calculate_neq(const opp_dat n_type_dat) 
{
    for (int i = 0; i < n_nodes_set; i++) 
        if (opp_get_data<OPP_INT>(n_type_dat)[i] == NORMAL || opp_get_data<OPP_INT>(n_type_dat)[i] == OPEN) 
            neq++;
}

//*************************************************************************************************
void FESolver::init_node_to_eq_map(const opp_dat n_type_dat) // relation from node indices to equation indices
{ OPP_RETURN_IF_INVALID_PROCESS;

    int eq = 0;
    for (int n = 0; n < n_nodes_inc_halo; n++)
    {
        if (n < n_nodes_set && 
            (opp_get_data<OPP_INT>(n_type_dat)[n] == NORMAL || opp_get_data<OPP_INT>(n_type_dat)[n] == OPEN))
        {
            node_to_eq_map[n] = (eq++);
        }
        else 
        {
            node_to_eq_map[n] = INT_MIN;    /* dirichlet node or for halo region, until updated by MPI comm */
        }
    }

#ifdef USE_MPI
/*
Idea: with MPI there are halo regions that we need to take care of. Send the node_to_eq_map values of the 
export non-exec halos to correct ranks, by adding own_start+1 and making the value negative. The +1 is to 
avoid misunderstanding of mapping ZERO. When node_to_eq_map is accessed and if it is negative and greater
than INT_MIN, we know thats its from a halo, so do the reverse calculation to get the actual global mapping.
*/
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
                send_buffer[index] = -1 * (send_buffer[index] + own_start + 1);
        }

        MPI_Isend(&send_buffer[exp_nonexec_list->disps[i]], exp_nonexec_list->sizes[i], MPI_INT,
                exp_nonexec_list->ranks[i], 11, OPP_MPI_WORLD, &(send_req[i]));
    }

    const int nonexec_init = n_type_dat->set->size + n_type_dat->set->exec_size;
    for (int i = 0; i < imp_nonexec_list->ranks_size; i++) 
    {
        MPI_Recv(&(node_to_eq_map[nonexec_init + imp_nonexec_list->disps[i]]), imp_nonexec_list->sizes[i], 
            MPI_INT, imp_nonexec_list->ranks[i], 11, OPP_MPI_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Waitall(exp_nonexec_list->ranks_size, send_req.data(), MPI_STATUSES_IGNORE);
#endif
}

//*************************************************************************************************
void FESolver::init_petsc_structures()
{
    MatCreate(PETSC_COMM_WORLD, &Jmat);
    MatCreate(PETSC_COMM_WORLD, &Kmat);

    MatSetSizes(Jmat, PETSC_DECIDE, PETSC_DECIDE, global_neq, global_neq);
    MatSetSizes(Kmat, PETSC_DECIDE, PETSC_DECIDE, global_neq, global_neq);

#ifndef USE_MPI
    MatSetType(Jmat, MATSEQAIJ); MatSetType(Kmat, MATSEQAIJ);
    VecCreateSeq(PETSC_COMM_WORLD, global_neq, &Bvec);
#else
    MatSetType(Jmat, MATMPIAIJ); MatSetType(Kmat, MATMPIAIJ);
    VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, global_neq, &Bvec); 
#endif

    MatSetUp(Jmat); MatSetUp(Kmat);

    VecSet(Bvec, 0.0);
    VecAssemblyBegin(Bvec); VecAssemblyEnd(Bvec);

    duplicate_vec(&Bvec, &F0vec);
    duplicate_vec(&Bvec, &F1vec);
    // duplicate_vec(&Bvec, &Dvec); // Dvec is created as a ghost vector below
    duplicate_vec(&Bvec, &Yvec);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetOperators(ksp, Jmat, Jmat);
    KSPSetTolerances(ksp, 1.e-10, 1.e-10, PETSC_DEFAULT , PETSC_DEFAULT); 
    KSPSetFromOptions(ksp); 
    //PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_MATLAB);

    vec_col.resize(neq);
    std::vector<int> ghost_vec_indices, ghost_indices;

    int pet_own_start = 0, pet_own_end = 0;
    VecGetOwnershipRange(Bvec, &pet_own_start, &pet_own_end);
    
    if (OPP_IS_VALID_PROCESS) {

        std::vector<int> neq_per_rank(OPP_comm_size, 0);
#ifdef USE_MPI       
        MPI_Allgather(&neq, 1, MPI_INT, neq_per_rank.data(), 1, MPI_INT, OPP_MPI_WORLD);
#endif

        for (int i = 0; i < OPP_rank; ++i)
            own_start += neq_per_rank[i];
        own_end = own_start + neq;

        for (int i = 0; i < neq; i++) 
            vec_col[i] = (i + own_start);

        for (int i = pet_own_start; i < pet_own_end; ++i) {
            ghost_vec_indices.push_back(i);
        }
        for (auto gbl_idx : vec_col) {
            if (gbl_idx < pet_own_start || gbl_idx >= pet_own_end) {
                ghost_vec_indices.push_back(gbl_idx);
                ghost_indices.push_back(gbl_idx);
            }
        }

        for (size_t i = 0; i < vec_col.size(); ++i) {
            auto it = std::find(ghost_vec_indices.begin(), ghost_vec_indices.end(), vec_col[i]);
            
            if (it != ghost_vec_indices.end()) {
                ex_indices.push_back(std::distance(ghost_vec_indices.begin(), it));
            } 
            else {
                opp_printf("ERROR", 
                    "init_petsc_structures - vec_col value %d not in ghost_vec_indices", i);
            }
        }
    }

    // if (OPP_DBG)
        opp_printf("init_petsc_structures", 
            "pet_own_start %d pet_own_end %d | neq %d | own_start %d own_end %d | ghost_indices %zu", 
            pet_own_start, pet_own_end, neq, own_start, own_end, ghost_indices.size());

    VecCreateGhost(PETSC_COMM_WORLD, (pet_own_end - pet_own_start), global_neq, 
                    ghost_indices.size(), ghost_indices.data(), &Dvec);
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
        double ke[4][4], fe[4]; /*force vector*/

        for (int a = 0; a < 4; a++)
        {
            for (int b = 0; b < 4; b++) 
            {
                ke[a][b] = 0.0;    /*reset*/

                /*perform quadrature*/
                for (int k = 0; k < n_int; k++)
                    for (int j = 0; j < n_int; j++)
                        for (int i = 0; i < n_int; i++) 
                        {
                            double dot = 0.0; /*dot product*/
                            for (int d = 0; d < DIM; d++) 
                                dot += (NX[e][a][d] * NX[e][b][d]);
                            ke[a][b] += (dot * detJ[e] * W[i] * W[j] * W[k]);
                        }
            }
        }
        add_ke(sparse_K, e, ke);

        for (int a = 0; a < 4; a++) 
        {  
            double fh = 0;  /*second term int(na*h), always zero since support only h=0*/
            double fg = 0;  /*third term, -sum(kab*qb)*/

            for (int b = 0; b < 4; b++) 
            {
                const int node_idx = c2n_map->map[e * c2n_map->dim + b]; // This can reach import halos
                fg -= (ke[a][b] * opp_get_data<OPP_REAL>(n_bnd_pot)[node_idx]);
            }

            fe[a] = fh + fg; /*combine*/
        }

        add_fe(&F0vec, e, fe);

    }  /*end of element*/

    VecAssemblyBegin(F0vec); VecAssemblyEnd(F0vec);

    initialze_matrix(sparse_K);

    const double scalar2 = 1.0 / (SCALLING * SCALLING); 

    VecScale(F0vec, scalar2);   // downscalling since NX is scalled to avoid precision issues
    MatScale(Kmat, scalar2);    // downscalling since NX is scalled to avoid precision issues

    sparse_K.clear();
}

//*************************************************************************************************
void FESolver::initialze_matrix(std::map<int, std::map<int, double>>& sparse_K)
{
    if (OPP_DBG) 
        opp_printf("FESolver", "initialze_matrix START own_start %d own_end %d", own_start, own_end);

    std::vector<int> matCol(neq);                         // Number of non zero columns per row
    std::vector<std::vector<int>> matIndex(neq);          // Non zero column indices per row
    std::vector<std::vector<double>> matIndexValues(neq); // Non zero column values per row (temp copy to enrich K matrix)

    int diag_max_fields = 0, off_diag_max_fields = 0;
    
    for (int i=own_start; i<own_end; i++) { // iterate own rows
    
        std::vector<int> tempVec;
        int diag = 0, off_diag = 0;
        
        std::map<int, double>& K_sparse_row = sparse_K[i-own_start];

        for (int j=0; j<global_neq; j++) { // iterate all columns
            
            double value = 0.0;
            auto it = K_sparse_row.find(j);
            if (it != K_sparse_row.end())
                value = it->second;
       
            if ((std::abs(value) > 1e-12) || (i == j)) { // significant ones and all diagonals  
                tempVec.push_back(j); 
                if (j >= own_start && j < own_end) diag++;
                else off_diag++;
            }         
        }

        matCol[i-own_start] = tempVec.size();
        matIndex[i-own_start].resize(tempVec.size());
        matIndexValues[i-own_start].resize(tempVec.size());

        std::copy(tempVec.begin(), tempVec.end(), matIndex[i-own_start].data());

        diag_max_fields = (diag_max_fields > diag) ? diag_max_fields : diag; 
        off_diag_max_fields = (off_diag_max_fields > off_diag) ? off_diag_max_fields : off_diag;  
    }

// #ifndef USE_MPI
//     MatSeqAIJSetPreallocation(Kmat, diag_max_fields, nullptr);
//     MatSeqAIJSetPreallocation(Jmat, diag_max_fields, nullptr);
// #else
//     MatMPIAIJSetPreallocation(Kmat, diag_max_fields, NULL, off_diag_max_fields, NULL);
//     MatMPIAIJSetPreallocation(Jmat, diag_max_fields, NULL, off_diag_max_fields, NULL);
// #endif
//     MatSetOption(Kmat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
//     MatSetOption(Jmat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

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
void FESolver::duplicate_vec(Vec* vec_mimic, Vec* vec_new)
{
    VecDuplicate(*vec_mimic, vec_new);
    VecSet(*vec_new, 0.0);
    VecAssemblyBegin(*vec_new); VecAssemblyEnd(*vec_new);
}

//*************************************************************************************************
void FESolver::enrich_cell_shape_deriv(opp_dat cell_shape_deriv) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    if (OPP_DBG) opp_printf("FESolver", "enrich_cell_shape_deriv START");

    opp_profiler->start("EnrichCellShapeDeriv");

    OPP_REAL* sd = opp_get_data<OPP_REAL>(cell_shape_deriv);
    for (int cellID = 0; cellID < n_cells_set; cellID++)
    {
        for (int nodeCon = 0; nodeCon < N_PER_C; nodeCon++)
        {
            for (int d = 0; d < DIM; d++)
                sd[cellID * (N_PER_C*DIM) + nodeCon * DIM + d] = (NX[cellID][nodeCon][d] / SCALLING);
        }
    }

    cell_shape_deriv->dirtybit = 1; // make it dirty to force update
    opp_mpi_force_halo_update_if_dirty(cell_shape_deriv->set, {cell_shape_deriv}, Device_CPU);
    cell_shape_deriv->dirty_hd = Dirty::Device; // To make GPU versions to download updated data!

    opp_profiler->end("EnrichCellShapeDeriv");

    if (OPP_DBG) opp_printf("FESolver", "enrich_cell_shape_deriv END");
}

/*adds contributions from element stiffness matrix*/
void FESolver::add_ke(std::map<int, std::map<int, double>>& sparse_K, int e, double ke[4][4])
{ 
    if (neq <= 0) return;

    for (int a=0; a<4; a++)         /*tetrahedra*/
    {
        for (int b=0; b<4; b++) 
        {
            const int node_idx1 = c2n_map->map[e * c2n_map->dim + a];
            const int P = node_to_eq_map[node_idx1]; 

            if (P < 0) continue;    /* skip g nodes or not in current ranks own row range */

            const int node_idx2 = c2n_map->map[e * c2n_map->dim + b];
            int Q = node_to_eq_map[node_idx2]; 

            if (Q < 0)
            {
                if (Q == INT_MIN) continue; // skip g nodes
                else Q = (-1 * Q - 1);      // on a different MPI rank, get its correct column
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
    
    for (int a=0; a<4; a++)    /*tetrahedra*/ 
    {
        const int node_idx = c2n_map->map[e * c2n_map->dim + a];
        const int P = node_to_eq_map[node_idx]; 
        if (P < 0) continue;    /* skip g nodes or on a different rank (owner compute)*/
        
        VecSetValue(*Fvec, P + own_start, fe[a], ADD_VALUES); // F[P] += fe[a];
    }
}

/*computes derivatives of the shape functions for all elements constants since using linear elements*/
void FESolver::compute_nx(const opp_dat n_pos) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    /*derivatives of the shape functions vs. xi*/
    double na_xi[4][3] = {{1,0,0}, {0,1,0}, {0,0,1}, {-1,-1,-1}};

    for (int e = 0; e < n_cells_inc_halo; e++) // for all cells inc exec halos
    {       
        const int* n_indices = &((int*)c2n_map->map)[e * N_PER_C]; 
        double x_xi[3][3], xi_x[3][3];

        for (int dim = 0; dim < 3; dim++)
        {
            for (int j = 0; j < 3; j++) /*xi/eta/zeta*/ 
            {
                x_xi[dim][j] = 0;
                for (int a = 0; a < 4; a++)    /*tet node*/
                {
                    x_xi[dim][j] += (na_xi[a][j] * 
                        opp_get_data<OPP_REAL>(n_pos)[n_indices[a] * DIM + dim]);
                }
            }
        }

        detJ[e] = det3(x_xi);

        inverse(x_xi,xi_x);
      
        for (int a = 0; a < 4; a++) // for all nodes in tetra
        {
            for (int d = 0; d < 3; d++) // for all dims
            {
                NX[e][a][d] = 0.0;
                for (int k = 0; k < 3; k++)
                { // scalled to avoid precision issues
                    NX[e][a][d] += (na_xi[a][k] * xi_x[k][d] * SCALLING); 
                }
            }
        }
    }
}

void FESolver::sanity_check()
{
    if (OPP_DBG)
        opp_printf("FESolver", "own neq is %d, the global neq is %d, start_eq is %d, end_eq %d -------", 
            neq, global_neq, own_start, own_end);

    int j_gm, j_gn, k_gm, k_gn, j_row_start, j_row_end, k_row_start, k_row_end;

    MatGetSize(Jmat, &j_gm, &j_gn);
    MatGetSize(Kmat, &k_gm, &k_gn);

    MatGetOwnershipRange(Jmat, &j_row_start, &j_row_end);
    MatGetOwnershipRange(Kmat, &k_row_start, &k_row_end);

    if (is_neq(j_gm, j_gn) || is_neq(k_gm, k_gn) || is_neq(j_gm, global_neq) ||
        is_neq(j_row_start, k_row_start) || is_neq(j_row_end, k_row_end))
    {
        opp_printf("FESolver::FESolver", "Error... Matrix vector sizes issue");
    }      
    
    if (OPP_DBG) 
        opp_printf("FESolver::FESolver", 
            "j_gm %d, j_gn %d, k_gm %d, k_gn %d, j_row_start %d, j_row_end %d, k_row_start %d, k_row_end %d", 
            j_gm, j_gn, k_gm, k_gn, j_row_start, j_row_end, k_row_start, k_row_end);
}


#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdio.h>

#include "FESolver.h"
#include "fempic.h"

//*************************************************************************************************
FESolver::FESolver(
    opp::Params* params, 
    oppic_set cells_set, 
    oppic_map cell_to_nodes_map, 
    oppic_dat node_type_dat, 
    oppic_dat node_pos_dat,
    oppic_dat node_bnd_pot_dat,
    int argc, char **argv) 
{
    phi0 = 0;
    n0 = params->get<OPP_REAL>("plasma_den");
    kTe = Kb * params->get<OPP_REAL>("electron_temperature");
    wall_potential = -(params->get<OPP_REAL>("wall_potential"));

    n_nodes    = node_type_dat->set->size + node_type_dat->set->exec_size + node_type_dat->set->nonexec_size; // DO NOT CALCULATE SOLUTION FOR IMPORT NON EXEC 
    n_elements = cells_set->size + node_type_dat->set->exec_size;
    neq        = 0; /*count number of unknowns*/
    
    /*OPEN nodes are "h" nodes*/
    for (std::size_t i = 0; i < n_nodes; i++) // TODO : DO NOT CALCULATE SOLUTION FOR IMPORT NON EXEC
        if (((int*)node_type_dat->data)[i] == NORMAL || ((int*)node_type_dat->data)[i] == OPEN) 
            neq++;

    summarize(std::cout);

    ID = new int[n_nodes];
   
    LM = new int*[n_elements]; /*allocate location matrix, n_elements*4 */
    for (int e=0;e<n_elements;e++) 
        LM[e] = new int[4];

    /*allocate NX matrix*/
    NX = new double**[n_elements];
    for (int e=0;e<n_elements;e++) 
    {
        NX[e] = new double*[4];
        for (int a=0;a<4;a++) 
            NX[e][a] = new double[3];
    }
    
    d = new double[neq]; /*solution array*/
    for (int n=0;n<neq;n++) 
        d[n]=0;    /*initial guess*/

    detJ = new double[n_elements];

    /*set up the ID array note valid values are 0 to neq-1 and -1 indicates "g" boundary node*/
    int P=0;
    for (int n=0;n<n_nodes;n++)
    {
        if (((int*)node_type_dat->data)[n] == NORMAL || ((int*)node_type_dat->data)[n] == OPEN)
        {
            ID[n]=P;
            P++;
        }
        else 
            ID[n]=-1;    /*dirichlet node*/
    }

    /*now set up the LM matrix*/
    for (int e=0;e<n_elements;e++)
        for (int a=0;a<4;a++)    /*tetrahedra*/ 
        {
            int node_idx = cell_to_nodes_map->map[e * cell_to_nodes_map->dim + a]; // TODO_PET : This can reach the halo region :()
            LM[e][a] = ID[node_idx]; 
        }

    /*set quadrature points*/
    l[0]=-sqrt(1.0/3.0); 
    l[1]=sqrt(1.0/3.0);
    W[0]=1; 
    W[1]=1;
    n_int = 2;

    computeNX(node_pos_dat, cell_to_nodes_map); /*compute NX matrix*/

    /* Petsc */
    PetscInitialize(&argc, &argv, PETSC_NULL, "FESolver::Petsc");

    vecCol           = new int[neq];
    matCol           = new int[neq];
    matIndex         = new int*[neq];
    matIndexValues   = new double*[neq];

    for (int i = 0; i < neq; i++)
    {
        vecCol[i] = i;
        matCol[i] = 0;
    }

    MatCreate(PETSC_COMM_WORLD, &Jmat);
    MatSetType(Jmat, MATSEQAIJ);
    MatSetSizes(Jmat, PETSC_DECIDE, PETSC_DECIDE, neq, neq);
    MatSetFromOptions(Jmat);
    MatSetUp(Jmat);

    MatCreate(PETSC_COMM_WORLD, &Kmat);
    MatSetType(Kmat, MATSEQAIJ);
    MatSetSizes(Kmat, PETSC_DECIDE, PETSC_DECIDE, neq, neq);
    MatSetFromOptions(Kmat);
    MatSetUp(Kmat);

    VecCreateSeq(PETSC_COMM_WORLD, neq, &Bvec);
    VecSetFromOptions(Bvec);
    
    VecDuplicate(Bvec, &Xvec);
    VecDuplicate(Bvec, &F0vec);
    VecDuplicate(Bvec, &F1vec);
    VecDuplicate(Bvec, &Dvec);
    VecDuplicate(Bvec, &Yvec);

    VecSet(F0vec, 0.0);
    VecSet(F1vec, 0.0);
    VecSet(Dvec, 0.0);
    VecSet(Yvec, 0.0);
    VecAssemblyBegin(Dvec); VecAssemblyEnd(Dvec);
    VecAssemblyBegin(Yvec); VecAssemblyEnd(Yvec);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetType(ksp, KSPCG);
    KSPSetOperators(ksp, Jmat, Jmat);
    // KSPSetTolerances(ksp, 1.e-2 / (neq * neq), 1.e-50, PETSC_DEFAULT, PETSC_DEFAULT); 
    KSPSetTolerances(ksp, 1.e-100 / (neq * neq), 1.e-100, PETSC_DEFAULT, PETSC_DEFAULT); 
    KSPSetFromOptions(ksp); 

    preAssembly(cell_to_nodes_map, node_bnd_pot_dat);
}

//*************************************************************************************************
FESolver::~FESolver() 
{ //TRACE_ME;

    for (int e=0;e<n_elements;e++) 
    {
        for (int a=0;a<4;a++) 
            delete[] NX[e][a];
        delete[] NX[e];
        delete[] LM[e];
    }

    delete[] LM;
    delete[] NX;
    delete[] ID;
    delete[] d;
    delete[] detJ;

    /* Petsc */
    delete[] vecCol;
    delete[] matCol;
    for (int i=0; i<neq; i++) { delete[] matIndex[i]; delete[] matIndexValues[i]; }
    delete[] matIndex;
    delete[] matIndexValues;
    
    KSPDestroy(&ksp);
    VecDestroy(&Xvec);
    VecDestroy(&Bvec);
    VecDestroy(&F0vec);
    VecDestroy(&F1vec);
    VecDestroy(&Dvec);
    VecDestroy(&Yvec);
    MatDestroy(&Jmat);
    MatDestroy(&Kmat);

    PetscFinalize();
}

//*************************************************************************************************
/*preassembles the K matrix and "h" and "g" parts of the force vector*/
void FESolver::preAssembly(oppic_map cell_to_nodes_map, oppic_dat node_bnd_pot) 
{ //TRACE_ME;

    double **K = new double*[neq];
    for (int i=0;i<neq;i++)
    {
        K[i] = new double[neq];
        
        for (int j=0;j<neq;j++) 
            K[i][j] = 0;
    }
        
    for (int e=0;e<n_elements;e++) 
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

                            getNax(nax,e,a);
                            getNax(nbx,e,b);

                            /*dot product*/
                            double dot=0;
                            for (int d=0;d<3;d++) dot+=nax[d]*nbx[d];
                            ke[a][b] += dot*detJ[e]*W[i]*W[j]*W[k];
                        }
            }

        addKe(K, e,ke);
   
        double fe[4]; /*force vector*/

        for (int a=0;a<4;a++) 
        {
            /*second term int(na*h), always zero since support only h=0*/
            double fh=0;

            /*third term, -sum(kab*qb)*/
            double fg = 0;
            for (int b=0;b<4;b++) 
            {
                int node_idx = cell_to_nodes_map->map[e * cell_to_nodes_map->dim + b]; // TODO_PET : This can reach import halos, g include halos
                double gb = ((double*)node_bnd_pot->data)[node_idx];
                fg-=ke[a][b]*gb;
            }

            /*combine*/
            fe[a] = fh + fg;
        }

        addFe(&F0vec, e,fe);

    }  /*end of element*/

    VecAssemblyBegin(F0vec); VecAssemblyEnd(F0vec);

    initialzeMatrix(K);

    for (int i=0;i<neq;i++) 
        delete[] K[i];
    delete[] K;
}

//*************************************************************************************************
/*computes "ff" part of F*/
void FESolver::buildF1Vector(double *ion_den) 
{ //TRACE_ME;

    double *f = new double[neq];

    for (int n=0;n<n_nodes;n++) 
    {
        int A = ID[n];
        if (A<0) continue;    /*skip known nodes*/
        f[A] = (QE/EPS0)*(ion_den[n]+n0*exp((d[A]-phi0)/kTe));
    }

    for (int e=0;e<n_elements;e++) 
    {
        double fe[4];
        for (int a=0;a<4;a++) 
        {
            /*first term is int(na*f), set to zero for now*/
            double ff=0;
            int A = LM[e][a];
            if (A>=0)    /*if unknown node*/ 
            {
                /*perform quadrature*/
                for (int k=0;k<n_int;k++)
                    for (int j=0;j<n_int;j++)
                        for (int i=0;i<n_int;i++) 
                        {
                            /*change of limits*/
                            double xi = 0.5*(l[i]+1);
                            double eta = 0.5*(l[j]+1);
                            double zeta = 0.5*(l[k]+1);

                            double Na=evalNa(a,xi,eta,zeta);
                            ff += f[A]*Na*detJ[e]*W[i]*W[j]*W[k];
                        }

                ff*=(1.0/8.0);    /*change of limits*/
                fe[a] = ff;
            }
        }

        addFe(&F1vec, e,fe);
    }

    VecAssemblyBegin(F1vec); VecAssemblyEnd(F1vec);

    delete[] f;
}

//*************************************************************************************************
/*builds J matrix for solver*/
void FESolver::buildJmatrix() 
{ //TRACE_ME;

    double *fp_term = new double[neq];
    double *FP = new double[neq];

    for (int n=0;n<neq;n++) 
        FP[n] = 0;

    for (int n=0;n<neq;n++) 
        fp_term[n] = -QE/EPS0*n0*exp((d[n]-phi0)/kTe)*(1/kTe);

    MatCopy(Kmat, Jmat, DIFFERENT_NONZERO_PATTERN);

    double fe[4];  /*build fprime vector*/

    for (int e=0;e<n_elements;e++) 
    {
        for (int a=0;a<4;a++) 
        {
            double ff=0;
            int A = LM[e][a];
            if (A>=0)    /*if unknown node*/ 
            {
                /*perform quadrature*/
                for (int k=0;k<n_int;k++)
                    for (int j=0;j<n_int;j++)
                        for (int i=0;i<n_int;i++) 
                        {
                            /*change of limits*/
                            double xi = 0.5*(l[i]+1);
                            double eta = 0.5*(l[j]+1);
                            double zeta = 0.5*(l[k]+1);

                            double Na=evalNa(a,xi,eta,zeta);
                            ff += fp_term[A]*Na*detJ[e]*W[i]*W[j]*W[k];
                        }

                ff*=(1.0/8.0);    /*change of limits*/
            }

            fe[a] = ff;
        }

        /*assembly*/
        for (int a=0;a<4;a++)    /*tetrahedra*/ 
        {
            int P = LM[e][a];
            if (P<0) continue;    /*skip g nodes*/

            FP[P] += fe[a];
        }
    }

    /*subtract diagonal term*/
    for (int u=0;u<neq;u++) 
    {
        MatSetValue(Jmat, u, u, (-FP[u]), ADD_VALUES);     // J[u][u]-=FP[u];
    }

    MatAssemblyBegin(Jmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Jmat, MAT_FINAL_ASSEMBLY);

    delete[] fp_term;
    delete[] FP;
}


//*************************************************************************************************
/*wrapper for solving the non-linear Poisson's equation*/
void FESolver::computePhi(oppic_arg arg0, oppic_arg arg1, oppic_arg arg2) 
{ TRACE_ME;

    int nargs = 3;
    oppic_arg args[nargs];
    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    int set_size = oppic_mpi_halo_exchanges(arg1.dat->set, nargs, args);

    double *ion_den = (double*)arg1.dat->data;
    double *node_potential = (double*)arg0.dat->data;
    double *node_bnd_pot = (double*)arg2.dat->data;

    solve(ion_den);

    /* combine d (solution from linear solver) and boundary potential to get node_potential */
    for (int n=0;n<n_nodes;n++) 
    {
        node_potential[n] = node_bnd_pot[n]; /*zero on non-g nodes*/

        int A=ID[n];
        if (A>=0) /*is this a non-boundary node?*/
            node_potential[n] += d[A];
    }

    oppic_mpi_set_dirtybit(nargs, args);
}

//*************************************************************************************************
void FESolver::solve(double *ion_den) 
{ TRACE_ME;

    buildF1Vector(ion_den);

    MatMult(Kmat, Dvec, Bvec); // B = K * D
    
    VecAXPY(Bvec, -1.0, F0vec); // B = B - F0
    
    VecAXPY(Bvec, -1.0, F1vec); // B = B - F1

    buildJmatrix();

    KSPSolve(ksp, Bvec, Yvec);

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);
    if (reason < 0)
    {
        char* str_reason = "\0";
        KSPGetConvergedReasonString(ksp, (const char**)(&str_reason));
        std::cerr << "Petsc failed to converge : " << str_reason << "; run with -ksp_converged_reason"<< std::endl;
    }

    VecAXPY(Dvec, -1.0, Yvec); // D = D - Y

    VecGetValues(Dvec, neq, vecCol, d); // For the calculation at computePhi()
}

//*************************************************************************************************
void FESolver::summarize(std::ostream &out) 
{
    out << "FE SOLVER INFORMATION" << std::endl
        << "---------------------" << std::endl;
    out << "  Number of unkowns: " << neq << std::endl;
    out << std::endl;
}

//*************************************************************************************************
void FESolver::initialzeMatrix(double **p_A)
{
    if (!matIndexCreated)
    {
        int max_fields = 0;
        for (int i=0; i<neq; i++)
        {
            std::vector<int> tempVec;    
            for (int j=0; j<neq; j++)
            {    
                // significant ones and all diagonals            
                if ((std::abs(p_A[i][j]) > 1e-12) || (i == j)) tempVec.push_back(j);            
            }

            matCol[i] = tempVec.size();
            matIndex[i] = new int[tempVec.size()];
            matIndexValues[i] = new double[tempVec.size()];

            std::copy(tempVec.begin(), tempVec.end(), matIndex[i]);

            max_fields = (max_fields > tempVec.size()) ? max_fields : tempVec.size();
        }

        MatSeqAIJSetPreallocation(Kmat, max_fields, nullptr);
        MatSeqAIJSetPreallocation(Jmat, max_fields, nullptr);

        printf("FESolver::initialzeMatrix max_fields=%d\n\n", max_fields);
        matIndexCreated = true;
    }

    for (int i=0; i<neq; i++)
    {
        for (int j=0; j<matCol[i]; j++) matIndexValues[i][j] = p_A[i][matIndex[i][j]];

        MatSetValues(Kmat, 1, &i, matCol[i], matIndex[i], matIndexValues[i], INSERT_VALUES); 
    }

    MatAssemblyBegin(Kmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Kmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(Jmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Jmat, MAT_FINAL_ASSEMBLY);
}

//*************************************************************************************************

// UTIL FUNCTIONS

void FESolver::enrich_cell_shape_deriv(oppic_dat cell_shape_deriv)
{
    // copying only up to set size, hence import exec halos will be zero
    for (int cellID = 0; cellID < cell_shape_deriv->set->size; cellID++)
    {
        for (int nodeCon = 0; nodeCon < NODES_PER_CELL; nodeCon++)
        {
            ((double*)cell_shape_deriv->data)[cellID * (NODES_PER_CELL*DIMENSIONS) + nodeCon * DIMENSIONS + 0 ] = NX[cellID][nodeCon][0];
            
            ((double*)cell_shape_deriv->data)[cellID * (NODES_PER_CELL*DIMENSIONS) + nodeCon * DIMENSIONS + 1 ] = NX[cellID][nodeCon][1];
            
            ((double*)cell_shape_deriv->data)[cellID * (NODES_PER_CELL*DIMENSIONS) + nodeCon * DIMENSIONS + 2 ] = NX[cellID][nodeCon][2];
        }
    }
}

/*adds contributions from element stiffness matrix*/
void FESolver::addKe(double** K, int e, double ke[4][4]) 
{ //TRACE_ME;

    for (int a=0;a<4;a++)    /*tetrahedra*/
        for (int b=0;b<4;b++) 
        {
            int P = LM[e][a];
            int Q = LM[e][b];
            if (P<0 || Q<0) continue;    /*skip g nodes*/

            K[P][Q] += ke[a][b];
        }
}

/*adds contributions from element force vector to a global F vector*/
void FESolver::addFe(Vec *Fvec, int e, double fe[4]) 
{ //TRACE_ME;

    for (int a=0;a<4;a++)    /*tetrahedra*/ 
    {
        int P = LM[e][a];
        if (P<0) continue;    /*skip g nodes*/
        
        VecSetValue(*Fvec, P, fe[a], ADD_VALUES); // F[P] += fe[a];
    }
}

/*evaluates shape function a at position (xi,eta,zeta)*/
double FESolver::evalNa(int a, double xi, double eta, double zeta) 
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

/*returns derivative of N[a] at some logical point since we are using linear elements, these are constant in each element*/
void FESolver::getNax(double nx[3], int e, int a) 
{
    for (int d=0;d<3;d++)
        nx[d] = NX[e][a][d];
}

/*computes derivatives of the shape functions for all elements constants since using linear elements*/
void FESolver::computeNX(oppic_dat node_pos, oppic_map cell_to_nodes_map) 
{ //TRACE_ME;

    /*derivatives of the shape functions vs. xi*/
    double na_xi[4][3] = {{1,0,0}, {0,1,0}, {0,0,1}, {-1,-1,-1}};

    for (int e=0;e<n_elements;e++) 
    {       
        int* map0idx = &((int*)cell_to_nodes_map->map)[e * cell_to_nodes_map->dim]; /*node positions*/

        double x[4][3];
        for (int a=0;a<4;a++) 
        {
            double *pos = &((double*)node_pos->data)[map0idx[a] * node_pos->dim];
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
                    NX[e][a][d]+=na_xi[a][k]*xi_x[k][d];
            }
    }
}

//*************************************************************************************************
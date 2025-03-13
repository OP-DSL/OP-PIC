#include "hip/hip_runtime.h"
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

#include "FESolver.h"
#include "opp_hip.h"

#define is_neq(a, b) ((a) != (b))

#define SCALLING 100
#define USE_HIP

bool print_petsc = false;

__constant__ double l_DEVICE_CONST[2];
__constant__ double W_DEVICE_CONST[2];
__constant__ int n_int_DEVICE_CONST;

__constant__ double QE_DEVICE_CONST;
__constant__ double EPS0_DEVICE_CONST;
__constant__ double n0_DEVICE_CONST;
__constant__ double phi0_DEVICE_CONST;
__constant__ double kTe_DEVICE_CONST;

#ifdef USE_PETSC
//*************************************************************************************************
FESolver::FESolver(
    opp_map c2n_map, 
    opp_dat n_type_dat, 
    opp_dat n_pos_dat,
    opp_dat n_bnd_pot_dat,
    int argc, char **argv) :
        c2n_map(c2n_map),
        n0(opp_params->get<OPP_REAL>("plasma_den")),
        kTe(Kb * opp_params->get<OPP_REAL>("electron_temperature")),
        wall_potential(-(opp_params->get<OPP_REAL>("wall_potential"))),
        n_nodes_set(n_type_dat->set->size),
        n_nodes_inc_halo(n_type_dat->set->size + n_type_dat->set->exec_size + n_type_dat->set->nonexec_size),
        n_elements_set(c2n_map->from->size),
        n_elements_inc_halo(c2n_map->from->size + c2n_map->from->exec_size)
{
    if (OPP_DBG) opp_printf("FESolver", "FESolver");

    // opp_printf("FESolver", "n_nodes_set %d | n_nodes_inc_halo %d | n_elements_set %d | n_elements_inc_halo %d", 
    //     n_nodes_set, n_nodes_inc_halo, n_elements_set, n_elements_inc_halo);

    // TODO : DO NOT CALCULATE SOLUTION FOR IMPORT NON EXEC, the owning rank will do that
    for (int i = 0; i < n_nodes_set; i++) 
        if (((int*)n_type_dat->data)[i] == NORMAL || ((int*)n_type_dat->data)[i] == OPEN) 
            neq++;

    node_to_eq_map = new int[n_nodes_inc_halo];

    NX = new double**[n_elements_inc_halo]; /*allocate NX matrix*/
    detJ = new double[n_elements_inc_halo];

    for (int e = 0; e < n_elements_inc_halo; e++) 
    {
        NX[e] = new double*[4];
        for (int a=0;a<4;a++) 
            NX[e][a] = new double[3];
    }
    
    dLocal = new double[neq];        /*solution array*/
    f1Local = new double[neq];
    tempNEQ1 = new double[neq];
    tempNEQ2 = new double[neq];
    tempNEQ3 = new double[neq];
    
    for (int n=0;n<neq;n++) 
    {
        dLocal[n] = 0;                 /*initial guess*/
        f1Local[n] = 0.0;
    }

    initPetscStructures();

    init_node_to_eq_map(n_type_dat);

    computeNX(n_pos_dat, c2n_map);

    summarize(std::cout);

    preAssembly(c2n_map, n_bnd_pot_dat);

    if (print_petsc) 
    {
        MatView(Kmat, PETSC_VIEWER_STDOUT_WORLD);
        VecView(F0vec, PETSC_VIEWER_STDOUT_WORLD);
    }

    sanityCheck();

#ifdef USE_HIP
    neqNBlocks            = (neq - 1) / OPP_gpu_threads_per_block + 1;
    nodesNBlocks          = (n_nodes_set - 1) / OPP_gpu_threads_per_block + 1;
    nodes_inc_haloNBlocks = (n_nodes_inc_halo - 1) / OPP_gpu_threads_per_block + 1;
    cells_inc_haloNBlocks = (n_elements_inc_halo - 1) / OPP_gpu_threads_per_block + 1;

    cutilSafeCall(hipMalloc(&dLocal_d, neq * sizeof(double)));
    cutilSafeCall(hipMalloc(&f1Local_d, neq * sizeof(double)));
    cutilSafeCall(hipMalloc(&tempNEQ1_d, neq * sizeof(double))); // Assigned, no need to init here
    cutilSafeCall(hipMalloc(&tempNEQ2_d, neq * sizeof(double))); // Assigned, no need to init here
    cutilSafeCall(hipMalloc(&tempNEQ3_d, neq * sizeof(double))); // Assigned, no need to init here
    cutilSafeCall(hipMalloc(&detJ_d, n_elements_inc_halo * sizeof(double)));
    cutilSafeCall(hipMalloc(&node_to_eq_map_d, n_nodes_inc_halo * sizeof(int)));

    cutilSafeCall(hipMemcpy(dLocal_d, dLocal, neq * sizeof(double), hipMemcpyHostToDevice));
    cutilSafeCall(hipMemcpy(f1Local_d, f1Local, neq * sizeof(double), hipMemcpyHostToDevice));
    cutilSafeCall(hipMemcpy(detJ_d, detJ, n_elements_inc_halo * sizeof(double), hipMemcpyHostToDevice));
    cutilSafeCall(hipMemcpy(node_to_eq_map_d, node_to_eq_map, n_nodes_inc_halo * sizeof(int), 
                                hipMemcpyHostToDevice));

    cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(l_DEVICE_CONST), l, sizeof(double) * 2));
    cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(W_DEVICE_CONST), W, sizeof(double) * 2));
    cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(n_int_DEVICE_CONST), &n_int, sizeof(int)));

    cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(QE_DEVICE_CONST), &QE, sizeof(double)));
    cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(EPS0_DEVICE_CONST), &EPS0, sizeof(double)));
    cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(n0_DEVICE_CONST), &n0, sizeof(double)));
    cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(phi0_DEVICE_CONST), &phi0, sizeof(double)));
    cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(kTe_DEVICE_CONST), &kTe, sizeof(double)));
#endif
}

//*************************************************************************************************
FESolver::~FESolver() 
{
#ifdef USE_HIP
    if (dLocal_d != NULL) cutilSafeCall(hipFree(dLocal_d));
    if (f1Local_d != NULL) cutilSafeCall(hipFree(f1Local_d));
    if (tempNEQ1_d != nullptr) cutilSafeCall(hipFree(tempNEQ1_d));
    if (tempNEQ2_d != nullptr) cutilSafeCall(hipFree(tempNEQ2_d));
    if (tempNEQ3_d != nullptr) cutilSafeCall(hipFree(tempNEQ3_d));
    if (detJ_d != nullptr) cutilSafeCall(hipFree(detJ_d));
    if (node_to_eq_map_d != nullptr) cutilSafeCall(hipFree(node_to_eq_map_d));
#endif

    for (int e = 0; e < n_elements_inc_halo; e++) 
    {
        for (int a=0;a<4;a++) 
            delete[] NX[e][a];
        delete[] NX[e];
    }

    delete[] NX;
    delete[] node_to_eq_map;
    delete[] dLocal;
    delete[] detJ;
    
    if (f1Local != nullptr) delete[] f1Local;

    if (tempNEQ1 != nullptr) delete[] tempNEQ1;
    if (tempNEQ2 != nullptr) delete[] tempNEQ2;
    if (tempNEQ3 != nullptr) delete[] tempNEQ3;

    /* Petsc */
    delete[] vecCol;
    
    KSPDestroy(&ksp);
    VecDestroy(&Bvec);
    VecDestroy(&F0vec);
    VecDestroy(&F1vec);
    VecDestroy(&Dvec);
    VecDestroy(&Yvec);
    MatDestroy(&Jmat);
    MatDestroy(&Kmat);
}

//*************************************************************************************************
void FESolver::init_node_to_eq_map(opp_dat n_type_dat) // relation from node indices to equation indices
{
    int P=0;
    for (int n = 0; n < n_nodes_inc_halo; n++)
    {
        if (n < n_nodes_set && 
            (((int*)n_type_dat->data)[n] == NORMAL || ((int*)n_type_dat->data)[n] == OPEN))
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

    int* send_buffer = new int[exp_nonexec_list->size];
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

    int nonexec_init = n_type_dat->set->size + n_type_dat->set->exec_size;
    for (int i = 0; i < imp_nonexec_list->ranks_size; i++) 
    {
        // opp_printf("SEND", "count=%d from rank %d", imp_nonexec_list->sizes[i], imp_nonexec_list->ranks[i]);

        MPI_Recv(&(node_to_eq_map[nonexec_init + imp_nonexec_list->disps[i]]), imp_nonexec_list->sizes[i], MPI_INT,
            imp_nonexec_list->ranks[i], 11, OPP_MPI_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Waitall(exp_nonexec_list->ranks_size, send_req.data(), MPI_STATUSES_IGNORE);

    send_req.clear();
    delete[] send_buffer;
#endif
}

//*************************************************************************************************
void FESolver::initPetscStructures()
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

    vecCol = new int[neq];
    for (int i = 0; i < neq; i++) vecCol[i] = i + own_start;
}

//*************************************************************************************************
/* preassembles the K matrix and "F0 vector */
void FESolver::preAssembly(opp_map c2n_map, opp_dat n_bnd_pot) 
{
    if (OPP_DBG)
        opp_printf("FESolver", "preAssembly global_neq %d neq %d n_elements_inc_halo %d", 
            global_neq, neq, n_elements_inc_halo);

    std::map<int, std::map<int, double>> sparse_K;

    for (int e = 0; e < n_elements_inc_halo; e++) 
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

        addKe(sparse_K, e,ke);
   
        double fe[4]; /*force vector*/

        for (int a=0;a<4;a++) 
        {  
            double fh=0;    /*second term int(na*h), always zero since support only h=0*/
            double fg = 0;  /*third term, -sum(kab*qb)*/

            for (int b=0;b<4;b++) 
            {
                // This can reach import halos
                int node_idx = c2n_map->map[e * c2n_map->dim + b]; 
                double gb = ((double*)n_bnd_pot->data)[node_idx];
                fg-=ke[a][b]*gb;
            }

            fe[a] = fh + fg; /*combine*/
        }

        addFe(&F0vec, e,fe);

    }  /*end of element*/

    VecAssemblyBegin(F0vec); VecAssemblyEnd(F0vec);

    initialzeMatrix(sparse_K);

    // double scalar1 = 1 / SCALLING;
    double scalar2 = 1 / (SCALLING * SCALLING); 

    VecScale(F0vec, scalar2);   // downscalling since NX is scalled to avoid precision issues
    MatScale(Kmat, scalar2);    // downscalling since NX is scalled to avoid precision issues

    sparse_K.clear();
}

//*************************************************************************************************
__global__ void initBuildF1VectorKernel(const int* node_to_eq_map_d, const double* dLocal_d, 
                                            const double* ion_den_d, double* tempNEQ1_d, const int n_nodes_inc_halo) 
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < n_nodes_inc_halo) {
        const int A = node_to_eq_map_d[tid];
        if (A >= 0) {
            tempNEQ1_d[A] = (QE_DEVICE_CONST / EPS0_DEVICE_CONST) * 
                             (ion_den_d[tid] + n0_DEVICE_CONST * exp((dLocal_d[A] - phi0_DEVICE_CONST) / kTe_DEVICE_CONST));
        }
    }
}

//*************************************************************************************************
__global__ void computeF1VectorValuesKernel(const int n_elements_inc_halo, const int* c2n_map_d,  
                    const int* node_to_eq_map_d, const double* tempNEQ1_d, double* f1Local_d, const double* detJ_d, const int stride) 
{
    int e = blockIdx.x * blockDim.x + threadIdx.x;

    if (e < n_elements_inc_halo) 
    {
        double Na = 0.0;
        double fe[4];
        double ff = 0.0;

        for (int a = 0; a < 4; a++) {
            const int node_idx = c2n_map_d[e + stride * a];  // Assuming 4 nodes per element
            const int A = node_to_eq_map_d[node_idx];

            if (A >= 0) {
                ff = 0.0;

                for (int k = 0; k < n_int_DEVICE_CONST; k++) {
                    for (int j = 0; j < n_int_DEVICE_CONST; j++) {
                        for (int i = 0; i < n_int_DEVICE_CONST; i++) {
                            switch (a) {
                                case 0: Na = 0.5 * (l_DEVICE_CONST[i] + 1); break;
                                case 1: Na = 0.5 * (l_DEVICE_CONST[j] + 1); break;
                                case 2: Na = 0.5 * (l_DEVICE_CONST[k] + 1); break;
                                case 3: Na = 1 - 0.5 * (l_DEVICE_CONST[i] + 1) - 0.5 * (l_DEVICE_CONST[j] + 1) 
                                                    - 0.5 * (l_DEVICE_CONST[k] + 1); break;
                                default: Na = 0;
                            }

                            ff += tempNEQ1_d[A] * Na * detJ_d[e] * W_DEVICE_CONST[i] * W_DEVICE_CONST[j] * W_DEVICE_CONST[k];
                        }
                    }
                }

                ff *= (1.0 / 8.0);
                fe[a] = ff;
            }
        }

        for (int a = 0; a < 4; a++) {
            const int node_idx = c2n_map_d[e + stride * a];
            const int P = node_to_eq_map_d[node_idx];

            if (P >= 0) {
                atomicAdd(&f1Local_d[P], fe[a]);
            }
        }
    }
}

//*************************************************************************************************
void FESolver::buildF1Vector(double *ion_den_d) 
{

}

//*************************************************************************************************
__global__ void initBuildJmatrixKernel(double *f2_d, double *f3_d, const double *d_d, const int neq)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    if (tid < neq)
    {
        f2_d[tid] = 0.0;
        f3_d[tid] = -QE_DEVICE_CONST / EPS0_DEVICE_CONST * n0_DEVICE_CONST * 
                        exp((d_d[tid] - phi0_DEVICE_CONST) / kTe_DEVICE_CONST) * (1 / kTe_DEVICE_CONST);
    }
}

//*************************************************************************************************
__global__ void computeJmatrixValuesKernel(const int n_elements_inc_halo, const int* c2n_map_d, 
                const int* node_to_eq_map_d, const double* detJ_d, double* tempNEQ3_d, double* tempNEQ2_d, const int stride) 
{
    int e = blockIdx.x * blockDim.x + threadIdx.x;

    if (e < n_elements_inc_halo) 
    {
        double fe[4];
        double Na = 0.0;

        for (int a = 0; a < 4; a++) {
            double ff = 0;

            //const int node_idx = c2n_map_d[e * 4 + a];  // Assuming 4 nodes per element
            const int node_idx = c2n_map_d[e + stride * a];  // Assuming 4 nodes per element
            const int A = node_to_eq_map_d[node_idx];

            if (A >= 0) {
                for (int k = 0; k < n_int_DEVICE_CONST; k++) {
                    for (int j = 0; j < n_int_DEVICE_CONST; j++) {
                        for (int i = 0; i < n_int_DEVICE_CONST; i++) {
                            switch (a) {
                                case 0: Na = 0.5 * (l_DEVICE_CONST[i] + 1); break;
                                case 1: Na = 0.5 * (l_DEVICE_CONST[j] + 1); break;
                                case 2: Na = 0.5 * (l_DEVICE_CONST[k] + 1); break;
                                case 3: Na = 1 - 0.5 * (l_DEVICE_CONST[i] + 1) - 0.5 * (l_DEVICE_CONST[j] + 1) 
                                                - 0.5 * (l_DEVICE_CONST[k] + 1); break;
                                default: Na = 0;
                            }

                            ff += tempNEQ3_d[A] * Na * detJ_d[e] * W_DEVICE_CONST[i] * W_DEVICE_CONST[j] * W_DEVICE_CONST[k];
                        }
                    }
                }
                ff *= (1.0 / 8.0);
            }

            fe[a] = ff;
        }

        for (int a = 0; a < 4; a++) {
            const int node_idx = c2n_map_d[e + stride * a];
            const int P = node_to_eq_map_d[node_idx];
            if (P >= 0) {  // Skip g nodes
                atomicAdd(&tempNEQ2_d[P], fe[a]);
            }
        }
    }
}

//*************************************************************************************************
void FESolver::buildJmatrix() 
{ 

}

//*************************************************************************************************
__global__ void computeNodePotentialKernel(const int n_nodes_set_d, const int* node_to_eq_map_d, const double* dLocal_d, 
                                            const double* n_bnd_pot_d, double* node_potential_d) 
{
    int n = blockIdx.x * blockDim.x + threadIdx.x;

    if (n < n_nodes_set_d) {
        node_potential_d[n] = n_bnd_pot_d[n]; /*zero on non-g nodes*/

        if (!(n_bnd_pot_d[n] == -10000000 || n_bnd_pot_d[n] == 0)) {
            printf("WRONG BND POT %2.20lE\n", n_bnd_pot_d[n]);
        }
        int A = node_to_eq_map_d[n];
        if (A >= 0) /*is this a non-boundary node?*/
            node_potential_d[n] += dLocal_d[A];
    }
}

//*************************************************************************************************
/*wrapper for solving the non-linear Poisson's equation*/
void FESolver::computePhi(opp_arg arg0, opp_arg arg1, opp_arg arg2) 
{ 
    if (OPP_DBG) opp_printf("FESolver", "computePhi START");

    opp_profiler->start("ComputePhi");

    int nargs = 3;
    opp_arg args[nargs];
    args[0] = arg0;  // node_potential
    args[1] = arg1;  // ion_den
    args[2] = arg2;  // n_bnd_pot
    args[1].idx = 2; // HACK to forcefully make halos to download

    opp_profiler->start("FSolve_halo");
    opp_mpi_halo_exchanges_grouped(arg1.dat->set, nargs, args, Device_GPU);
    opp_mpi_halo_wait_all(nargs, args);
    opp_profiler->end("FSolve_halo");

    if (linearSolve((double*)arg1.dat->data_d) == false)
    {
        char* str_reason;
        KSPGetConvergedReasonString(ksp, (const char**)(&str_reason));
        std::cerr << "linearSolve Petsc failed to converge : " << str_reason <<
            "; run with -ksp_converged_reason" << std::endl;           
    }

    opp_profiler->start("FSolve_npot");
    computeNodePotentialKernel <<<nodesNBlocks, OPP_gpu_threads_per_block>>> (n_nodes_set, node_to_eq_map_d, dLocal_d, 
                                                        (double*)arg2.dat->data_d, (double*)arg0.dat->data_d);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(hipDeviceSynchronize());
    opp_profiler->end("FSolve_npot");

    opp_profiler->end("ComputePhi");

    if (OPP_DBG) opp_printf("FESolver", "computePhi DONE");
}

/*
    build F1 using ion_density and D
    B = K * D
    B = B - F0
    B = B - F1
    Build J using K and D
    Ksp Solve : J * Y = B (Y = solution)
    D = D - Y
*/
//*************************************************************************************************
bool FESolver::linearSolve(double *ion_den_d) 
{ 

    if (OPP_DBG) opp_printf("FESolver", "linearSolve START");

    opp_profiler->start("FSolve_solv_init");
    initBuildF1VectorKernel <<<nodes_inc_haloNBlocks, OPP_gpu_threads_per_block>>> (node_to_eq_map_d, dLocal_d, 
                                                                        ion_den_d, tempNEQ1_d, n_nodes_inc_halo);

    initBuildJmatrixKernel <<<neqNBlocks, OPP_gpu_threads_per_block>>> (tempNEQ2_d, tempNEQ3_d, dLocal_d, neq);

    cutilSafeCall(hipDeviceSynchronize()); // Sync since we need tempNEQ1_d and tempNEQ3_d
    opp_profiler->end("FSolve_solv_init");

    opp_profiler->start("FSolve_j_and_f1");
    const int stride = c2n_map->from->size + c2n_map->from->exec_size + c2n_map->from->nonexec_size;

    computeF1VectorValuesKernel <<<cells_inc_haloNBlocks, OPP_gpu_threads_per_block>>> (n_elements_inc_halo, 
                                            c2n_map->map_d, node_to_eq_map_d, tempNEQ1_d, f1Local_d, detJ_d, stride);

    computeJmatrixValuesKernel <<<cells_inc_haloNBlocks, OPP_gpu_threads_per_block>>> (n_elements_inc_halo, 
                    c2n_map->map_d, node_to_eq_map_d, detJ_d, tempNEQ3_d, tempNEQ2_d, stride);

    MatCopy(Kmat, Jmat, DIFFERENT_NONZERO_PATTERN);

    MatMult(Kmat, Dvec, Bvec);  // B = K * D
    
    VecAXPY(Bvec, -1.0, F0vec); // B = B - F0

    cutilSafeCall(hipDeviceSynchronize()); // Sync since we need f1Local_d and tempNEQ2_d
    opp_profiler->end("FSolve_j_and_f1");

    cutilSafeCall(hipMemcpy(f1Local, f1Local_d, neq * sizeof(double), hipMemcpyDeviceToHost));
    cutilSafeCall(hipMemcpy(tempNEQ2, tempNEQ2_d, neq * sizeof(double), hipMemcpyDeviceToHost));

    VecSetValues(F1vec, neq, vecCol, f1Local, INSERT_VALUES);
    VecAssemblyBegin(F1vec); VecAssemblyEnd(F1vec);

    for (int u=0;u<neq;u++) {   /*subtract diagonal term*/
        MatSetValue(Jmat, (u + own_start), (u + own_start), (-tempNEQ2[u]), ADD_VALUES); // J[u][u]-=tempNEQ2[u];
    }
    MatAssemblyBegin(Jmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Jmat, MAT_FINAL_ASSEMBLY);

    VecAXPY(Bvec, -1.0, F1vec); // B = B - F1

    opp_profiler->start("FSolve_ksp");
    KSPSolve(ksp, Bvec, Yvec); // Jmat * Yvec = Bvec (Yvec = solution)
    opp_profiler->end("FSolve_ksp");

    VecAXPY(Dvec, -1.0, Yvec); // D = D - Y

    opp_profiler->start("FSolve_d");
    VecGetValues(Dvec, neq, vecCol, dLocal); // For the calculation at computePhi()

    cutilSafeCall(hipMemcpy(dLocal_d, dLocal, neq * sizeof(double), hipMemcpyHostToDevice));
    opp_profiler->end("FSolve_d");

    KSPGetConvergedReason(ksp, &reason);

    if (OPP_DBG) opp_printf("FESolver", "linearSolve DONE");

    return (reason >= 0);
}

//*************************************************************************************************
void FESolver::initialzeMatrix(std::map<int, std::map<int, double>>& sparse_K)
{
    if (OPP_DBG) opp_printf("FESolver", "initialzeMatrix START");

    int *matCol             = new int[neq];     // Number of non zero columns per row
    int **matIndex          = new int*[neq];    // Non zero column indices per row
    double **matIndexValues = new double*[neq]; // Non zero column values per row (temp copy to enrich K matrix)

    int diag_max_fields = 0, off_diag_max_fields = 0;
    
    for (int i=own_start; i<own_end; i++) // iterate own rows
    {
        std::vector<int> tempVec;
        int local_diag_max_fields = 0, local_off_diag_max_fields = 0;
        
        std::map<int, double>& K_sparse_row = sparse_K[i-own_start];

        for (int j=0; j<global_neq; j++) // iterate all columns
        {   
            double value = 0.0;
            auto it = K_sparse_row.find(j);
            if (it != K_sparse_row.end()) {
                value = it->second;
            }

            // significant ones and all diagonals         
            if ((std::abs(value) > 1e-12) || (i == j)) 
            {
                tempVec.push_back(j); 

                if (j >= own_start && j < own_end)
                    local_diag_max_fields++;
                else
                    local_off_diag_max_fields++;
            }         
        }

        matCol[i-own_start] = tempVec.size();
        matIndex[i-own_start] = new int[tempVec.size()];
        matIndexValues[i-own_start] = new double[tempVec.size()];

        std::copy(tempVec.begin(), tempVec.end(), matIndex[i-own_start]);

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

    if (OPP_DBG) opp_printf("FESolver", "initialzeMatrix diag_max_fields=%d off_diag_max_fields=%d", 
        diag_max_fields, off_diag_max_fields);

    for (int i=own_start; i<own_end; i++)
    {
        for (int j=0; j<matCol[i-own_start]; j++) 
            matIndexValues[i-own_start][j] = sparse_K[i-own_start][matIndex[i-own_start][j]];

        MatSetValues(Kmat, 1, &i, matCol[i-own_start], matIndex[i-own_start], 
            matIndexValues[i-own_start], INSERT_VALUES); 
    }

    MatAssemblyBegin(Kmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Kmat, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(Jmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Jmat, MAT_FINAL_ASSEMBLY);

    delete[] matCol;
    for (int i=0; i<neq; i++) { delete[] matIndex[i]; delete[] matIndexValues[i]; }
    delete[] matIndex;
    delete[] matIndexValues;

    if (OPP_DBG) opp_printf("FESolver", "initialzeMatrix END");
}

//*************************************************************************************************

// UTIL FUNCTIONS

//*************************************************************************************************
void FESolver::summarize(std::ostream &out) 
{
    // opp_printf("FESolver", "FE SOLVER INFORMATION");
    // opp_printf("FESolver", "---------------------");
    if (OPP_DBG)
        opp_printf("FESolver", "own neq is %d, the global neq is %d, start_eq is %d, end_eq %d -------", 
            neq, global_neq, own_start, own_end);
    // opp_printf("FESolver", "---------------------");
}

void FESolver::enrich_cell_shape_deriv(opp_dat cell_shape_deriv)
{
    if (OPP_DBG) opp_printf("FESolver", "enrich_cell_shape_deriv START");

    opp_profiler->start("EnrichCellShapeDeriv");

    // copying only up to set size, hence import exec halos will be zero
    for (int cellID = 0; cellID < n_elements_inc_halo; cellID++)
    {
        for (int nodeCon = 0; nodeCon < N_PER_C; nodeCon++)
        {
            ((double*)cell_shape_deriv->data)[cellID * (N_PER_C*DIM) + nodeCon * DIM + 0 ] = 
                (NX[cellID][nodeCon][0] / SCALLING);
            
            ((double*)cell_shape_deriv->data)[cellID * (N_PER_C*DIM) + nodeCon * DIM + 1 ] = 
                (NX[cellID][nodeCon][1] / SCALLING);
            
            ((double*)cell_shape_deriv->data)[cellID * (N_PER_C*DIM) + nodeCon * DIM + 2 ] = 
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
void FESolver::addKe(std::map<int, std::map<int, double>>& sparse_K, int e, double ke[4][4]) // BUG : K is not created correctly
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
void FESolver::addFe(Vec *Fvec, int e, double fe[4]) 
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

/*returns derivative of N[a] at some logical point since we are using linear elements, 
these are constant in each element*/
void FESolver::getNax(double nx[3], int e, int a) 
{
    for (int d=0;d<3;d++)
        nx[d] = NX[e][a][d];
}

/*computes derivatives of the shape functions for all elements constants since 
using linear elements*/
void FESolver::computeNX(opp_dat n_pos, opp_map c2n_map) 
{ 

    /*derivatives of the shape functions vs. xi*/
    double na_xi[4][3] = {{1,0,0}, {0,1,0}, {0,0,1}, {-1,-1,-1}};

    for (int e = 0; e < n_elements_inc_halo; e++) 
    {       
        /*node indices*/
        int* map0idx = &((int*)c2n_map->map)[e * c2n_map->dim]; 

        double x[4][3];
        for (int a=0;a<4;a++) 
        {
            /*node positions*/
            double *pos = &((double*)n_pos->data)[map0idx[a] * n_pos->dim]; 
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
void FESolver::sanityCheck()
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

//*************************************************************************************************

#else

    FESolver::FESolver(
        opp_map c2n_map, 
        opp_dat n_type, 
        opp_dat n_pos,  
        opp_dat n_bnd_pot,
        int argc, char **argv) : c2n_map(c2n_map) {};
    FESolver::~FESolver() {};

    void FESolver::computePhi(opp_arg arg0, opp_arg arg1, opp_arg arg2) {};
    
    void FESolver::preAssembly(opp_map c2n_map, opp_dat n_bnd_pot) {};
    void FESolver::enrich_cell_shape_deriv(opp_dat cell_shape_deriv) {};

    bool FESolver::linearSolve(double *ion_den) { return true; }; 
    void FESolver::buildJmatrix() {};
    void FESolver::buildF1Vector(double *ion_den) {};
    
    void FESolver::summarize(std::ostream &out) {};  

    void FESolver::addKe(std::map<int, std::map<int, double>>& sparse_K, int e, double ke[4][4]) {};
    void FESolver::addFe(Vec *Fvec, int e, double fe[4]) {};
    double FESolver::evalNa(int a, double xi, double eta, double zeta) { return -1.0; };
    void FESolver::getNax(double nx[3], int e, int a) {};
    void FESolver::initialzeMatrix(std::map<int, std::map<int, double>>& sparse_K) {};
    void FESolver::computeNX(opp_dat n_pos, opp_map c2n_map) {};
    void FESolver::sanityCheck() {};
    void FESolver::init_node_to_eq_map(opp_dat n_type_dat) {};
    void FESolver::initPetscStructures() {};
    
#endif
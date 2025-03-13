
#define DEVICE_TYPE Device_GPU

#if defined(USE_CUDA)
    #include "opp_cuda.h"
#elif defined(USE_HIP)
    #include "opp_hip.h"
#elif defined(USE_SYCL)
    #include "opp_sycl.h"
#endif

//*************************************************************************************************
void FESolver::init_device_variables() 
{ OPP_RETURN_IF_INVALID_PROCESS;

    OPP_RUN_ON_ROOT() opp_printf("FESolver", "Running on Device");

    neqNBlocks            = (neq - 1) / OPP_gpu_threads_per_block + 1;
    nodesNBlocks          = (n_nodes_set - 1) / OPP_gpu_threads_per_block + 1;
    nodes_inc_haloNBlocks = (n_nodes_inc_halo - 1) / OPP_gpu_threads_per_block + 1;
    cells_inc_haloNBlocks = (n_cells_inc_halo - 1) / OPP_gpu_threads_per_block + 1;

    opp_mem::copy_host_to_dev<double>(dLocal_d, dLocal.data(), neq, true, true, neq);
    opp_mem::copy_host_to_dev<double>(f1Local_d, f1Local.data(), neq, true, true, neq);
    opp_mem::copy_host_to_dev<double>(detJ_d, detJ.data(), n_cells_inc_halo, true, true, n_cells_inc_halo);
    opp_mem::copy_host_to_dev<int>(node_to_eq_map_d, node_to_eq_map.data(), n_nodes_inc_halo, 
                                    true, true, n_nodes_inc_halo);

    tempNEQ1_d = opp_mem::dev_malloc<double>(neq);
    tempNEQ2_d = opp_mem::dev_malloc<double>(neq);
    tempNEQ3_d = opp_mem::dev_malloc<double>(neq);

    opp_mem::copy_host_to_dev<double>(l_DEV_CONST, l, 2, true, true, 2);
    opp_mem::copy_host_to_dev<double>(W_DEV_CONST, W, 2, false, true, 2);
}

//*************************************************************************************************
void FESolver::destroy_device_variables() 
{ OPP_RETURN_IF_INVALID_PROCESS;

    opp_mem::dev_free(dLocal_d);
    opp_mem::dev_free(f1Local_d);
    opp_mem::dev_free(detJ_d);
    opp_mem::dev_free(node_to_eq_map_d);

    opp_mem::dev_free(tempNEQ1_d);
    opp_mem::dev_free(tempNEQ2_d);
    opp_mem::dev_free(tempNEQ3_d);
}

//*************************************************************************************************
OPP_GLOBAL_FUNCTION
void init_f1_kernel(const int* n2eq_map_l, const double* local_l, const double* ion_den_l, double* tempNEQ1_l, 
        const int n_nodes_inc_halo, double qe, double eps0, double n0, double phi0, double kTe ADDITIONAL_PARAMETERS) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;
    if (tid < n_nodes_inc_halo) {
        const int eq_idx = n2eq_map_l[tid];
        if (eq_idx >= 0) {
            tempNEQ1_l[eq_idx] = (qe / eps0) * (ion_den_l[tid] + 
                    n0 * exp((local_l[eq_idx] - phi0) / kTe));
        }
    }
}

OPP_GLOBAL_FUNCTION
void init_J_kernel(double *tempNEQ2_l, double *tempNEQ3_l, const double *dLocal_l, const int neq, 
        double qe, double eps0, double n0, double phi0, double kTe ADDITIONAL_PARAMETERS)
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;
    if (tid < neq) {
        tempNEQ2_l[tid] = 0.0;
        tempNEQ3_l[tid] = -qe / eps0 * n0 * 
                        exp((dLocal_l[tid] - phi0) / kTe) * (1 / kTe);
    }
}

OPP_GLOBAL_FUNCTION
void compute_f1_kernel(const int* c2n_map_l, const int* n2eq_map_l, const double* tempNEQ1_l,
        double* f1Local_l, const double* detJ_l, const int n_cells_inc_halo, const int cells_stride, 
        int n_int, const double* l_d, const double* w_d ADDITIONAL_PARAMETERS) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;
    if (tid < n_cells_inc_halo) {
        double Na = 0.0, ff = 0.0, fe[4];
        for (int a = 0; a < 4; a++) {
            const int node_idx = c2n_map_l[tid + cells_stride * a];  // Assuming 4 nodes per element
            const int eq_idx = n2eq_map_l[node_idx];
            if (eq_idx >= 0) {
                ff = 0.0;
                for (int k = 0; k < n_int; k++) {
                    for (int j = 0; j < n_int; j++) {
                        for (int i = 0; i < n_int; i++) {
                            switch (a) {
                                case 0: Na = 0.5 * (l_d[i] + 1); break;
                                case 1: Na = 0.5 * (l_d[j] + 1); break;
                                case 2: Na = 0.5 * (l_d[k] + 1); break;
                                case 3: Na = 1 - 0.5 * (l_d[i] + 1) - 0.5 * (l_d[j] + 1) 
                                                    - 0.5 * (l_d[k] + 1); break;
                                default: Na = 0;
                            }
                            ff += tempNEQ1_l[eq_idx] * Na * detJ_l[tid] * w_d[i] * w_d[j] * w_d[k];
                        }
                    }
                }
                ff *= (1.0 / 8.0);
                fe[a] = ff;
            }
        }
        for (int a = 0; a < 4; a++) {
            const int node_idx = c2n_map_l[tid + cells_stride * a];
            const int eq_idx = n2eq_map_l[node_idx];
            if (eq_idx >= 0) {
                OPP_ATOMIC_FETCH_ADD(&f1Local_l[eq_idx], fe[a]);
            }
        }
    }
}

OPP_GLOBAL_FUNCTION
void compute_J_kernel(const int* c2n_map_l, const int* n2eq_map_l, const double* detJ_l, 
        double* tempNEQ3_l, double* tempNEQ2_l, const int n_cells_inc_halo, const int cells_stride, 
        int n_int, const double* l_d, const double* w_d ADDITIONAL_PARAMETERS) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;
    if (tid < n_cells_inc_halo) { 
        double Na = 0.0, fe[4];
        for (int a = 0; a < 4; a++) {
            double ff = 0;
            const int node_idx = c2n_map_l[tid + cells_stride * a];  // Assuming 4 nodes per element
            const int eq_idx = n2eq_map_l[node_idx];
            if (eq_idx >= 0) {
                for (int k = 0; k < n_int; k++) {
                    for (int j = 0; j < n_int; j++) {
                        for (int i = 0; i < n_int; i++) {
                            switch (a) {
                                case 0: Na = 0.5 * (l_d[i] + 1); break;
                                case 1: Na = 0.5 * (l_d[j] + 1); break;
                                case 2: Na = 0.5 * (l_d[k] + 1); break;
                                case 3: Na = 1 - 0.5 * (l_d[i] + 1) - 0.5 * (l_d[j] + 1) 
                                                - 0.5 * (l_d[k] + 1); break;
                                default: Na = 0;
                            }
                            ff += tempNEQ3_l[eq_idx] * Na * detJ_l[tid] * 
                                            w_d[i] *  w_d[j] * w_d[k];
                        }
                    }
                }
                ff *= (1.0 / 8.0);
            }
            fe[a] = ff;
        }
        for (int a = 0; a < 4; a++) {
            const int node_idx = c2n_map_l[tid + cells_stride * a];
            const int eq_idx = n2eq_map_l[node_idx];
            if (eq_idx >= 0) {  // Skip g nodes
                OPP_ATOMIC_FETCH_ADD(&tempNEQ2_l[eq_idx], fe[a]);
            }
        }
    }
}

OPP_GLOBAL_FUNCTION 
void compute_node_potential_kernel(const int n_nodes_set_d, const int* node_to_eq_map_d, 
        const double* dLocal_d, const double* n_bnd_pot_d, double* node_potential_d ADDITIONAL_PARAMETERS) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;
    if (tid < n_nodes_set_d) {
        const int eq_idx = node_to_eq_map_d[tid];
        const double val = ((eq_idx >= 0) ? dLocal_d[eq_idx] : 0); /*is this a non-boundary node?*/
        node_potential_d[tid] = n_bnd_pot_d[tid] + val;
    }
}

//*************************************************************************************************
void FESolver::init_f1_and_J(const opp_dat ion_den_dat) 
{ OPP_RETURN_IF_INVALID_PROCESS;

    const double* ion_den_d = (double*)(ion_den_dat->data_d);

#ifdef USE_SYCL
    opp_queue->submit([&](sycl::handler &cgh) {
        const int* n2eq_map_l = node_to_eq_map_d;
        const double* local_l = dLocal_d;
        const double* ion_den_l = ion_den_d;
        double* tempNEQ1_l = tempNEQ1_d;
        const int n_nodes_inc_halo_l = n_nodes_inc_halo;
        const double qe = CONST_QE;
        const double eps0 = CONST_EPS0;
        const double n0 = CONST_n0;
        const double phi0 = CONST_phi0;
        const double kTe = CONST_kTe;

        cgh.parallel_for(sycl::nd_range<1>(OPP_gpu_threads_per_block*nodes_inc_haloNBlocks,OPP_gpu_threads_per_block), 
            [=](sycl::nd_item<1> item) {
                init_f1_kernel(n2eq_map_l, local_l, ion_den_l, tempNEQ1_l, n_nodes_inc_halo_l, 
                    qe, eps0, n0, phi0, kTe, item);
            });
    });  

    opp_queue->submit([&](sycl::handler &cgh) {
        const double* local_l = dLocal_d;
        double *tempNEQ2_l = tempNEQ2_d;
        double *tempNEQ3_l = tempNEQ3_d;
        const double qe = CONST_QE;
        const double eps0 = CONST_EPS0;
        const double n0 = CONST_n0;
        const double phi0 = CONST_phi0;
        const double kTe = CONST_kTe;
        const int neq_l = neq;
        
        cgh.parallel_for(sycl::nd_range<1>(OPP_gpu_threads_per_block*neqNBlocks,OPP_gpu_threads_per_block), 
            [=](sycl::nd_item<1> item) {
                init_J_kernel(tempNEQ2_l, tempNEQ3_l, local_l, neq_l, qe, eps0, n0, phi0, kTe, item);
            });
    });    
#else
    init_f1_kernel <<<nodes_inc_haloNBlocks, OPP_gpu_threads_per_block>>> (
        node_to_eq_map_d, 
        dLocal_d,
        ion_den_d, 
        tempNEQ1_d, 
        n_nodes_inc_halo, CONST_QE, CONST_EPS0, CONST_n0, CONST_phi0, CONST_kTe
    );
    init_J_kernel <<<neqNBlocks, OPP_gpu_threads_per_block>>> (
        tempNEQ2_d, 
        tempNEQ3_d, 
        dLocal_d, 
        neq, CONST_QE, CONST_EPS0, CONST_n0, CONST_phi0, CONST_kTe
    );
#endif

    OPP_DEVICE_SYNCHRONIZE();
}

//*************************************************************************************************
void FESolver::build_f1_vector() 
{ 
    if (OPP_IS_VALID_PROCESS) {
        const int cells_stride = c2n_map->from->size + c2n_map->from->exec_size + c2n_map->from->nonexec_size;

#ifdef USE_SYCL
        opp_queue->submit([&](sycl::handler &cgh) {
            const int* c2n_map_l = c2n_map->map_d;
            const int* n2eq_map_l = node_to_eq_map_d;
            const double* tempNEQ1_l = tempNEQ1_d;
            double* f1Local_l = f1Local_d;
            const double* detJ_l = detJ_d;
            const double *l_d = l_DEV_CONST;
            const double *w_d = W_DEV_CONST;
            const int n_cells_inc_halo_l = n_cells_inc_halo;
            const int n_int_l = n_int;

            cgh.parallel_for(sycl::nd_range<1>(OPP_gpu_threads_per_block*cells_inc_haloNBlocks,OPP_gpu_threads_per_block), 
                [=](sycl::nd_item<1> item) {
                    compute_f1_kernel(c2n_map_l, n2eq_map_l, tempNEQ1_l, f1Local_l, detJ_l, 
                        n_cells_inc_halo_l, cells_stride, n_int_l, l_d, w_d, item);
                });
        });
#else
        compute_f1_kernel <<<cells_inc_haloNBlocks, OPP_gpu_threads_per_block>>> (
            c2n_map->map_d, 
            node_to_eq_map_d, 
            tempNEQ1_d, 
            f1Local_d, 
            detJ_d, 
            n_cells_inc_halo, cells_stride, n_int, l_DEV_CONST, W_DEV_CONST
        );
#endif

        OPP_DEVICE_SYNCHRONIZE();

        opp_mem::copy_dev_to_host<double>(f1Local.data(), f1Local_d, neq);
    }

    VecSetValues(F1vec, neq, vec_col.data(), f1Local.data(), INSERT_VALUES);
    VecAssemblyBegin(F1vec); VecAssemblyEnd(F1vec);
}

//*************************************************************************************************
void FESolver::build_j_matrix() 
{   
    if (OPP_IS_VALID_PROCESS) {

        const int cells_stride = c2n_map->from->size + c2n_map->from->exec_size + c2n_map->from->nonexec_size;
#ifdef USE_SYCL
        opp_queue->submit([&](sycl::handler &cgh) {
            const int* c2n_map_l = c2n_map->map_d;
            const int* n2eq_map_l = node_to_eq_map_d;
            const double* detJ_l = detJ_d;
            double* tempNEQ3_l = tempNEQ3_d;
            double* tempNEQ2_l = tempNEQ2_d;
            const double *l_d = l_DEV_CONST;
            const double *w_d = W_DEV_CONST;
            const int n_cells_inc_halo_l = n_cells_inc_halo;
            const int n_int_l = n_int;

            cgh.parallel_for(sycl::nd_range<1>(OPP_gpu_threads_per_block*cells_inc_haloNBlocks,OPP_gpu_threads_per_block), 
                [=](sycl::nd_item<1> item) {
                    compute_J_kernel(c2n_map_l, n2eq_map_l, detJ_l, tempNEQ3_l, tempNEQ2_l,
                        n_cells_inc_halo_l, cells_stride, n_int_l, l_d, w_d, item);
                });
        });
#else
        compute_J_kernel <<<cells_inc_haloNBlocks, OPP_gpu_threads_per_block>>> (
            c2n_map->map_d, 
            node_to_eq_map_d, 
            detJ_d, 
            tempNEQ3_d, 
            tempNEQ2_d, 
            n_cells_inc_halo, cells_stride, n_int, l_DEV_CONST, W_DEV_CONST
        );
#endif

        OPP_DEVICE_SYNCHRONIZE();

        opp_mem::copy_dev_to_host<double>(tempNEQ2.data(), tempNEQ2_d, neq);
    }

    for (int u=0;u<neq;u++)   /* subtract diagonal term - J[u][u]-=tempNEQ2[u]; */
        MatSetValue(Jmat, (u + own_start), (u + own_start), (-tempNEQ2[u]), ADD_VALUES); 
    MatAssemblyBegin(Jmat, MAT_FINAL_ASSEMBLY); MatAssemblyEnd(Jmat, MAT_FINAL_ASSEMBLY);
}

//*************************************************************************************************
void FESolver::compute_node_potential(const opp_dat n_bnd_pot_dat, opp_dat node_potential_dat) {

    // VecGetValues(Dvec, neq, vec_col.data(), dLocal.data()); 
    VecGhostUpdateBegin(Dvec, INSERT_VALUES, SCATTER_FORWARD);
    VecGhostUpdateEnd(Dvec, INSERT_VALUES, SCATTER_FORWARD);
    VecGetArray(Dvec, &tmpDptr);
    for (int i = 0; i < neq; i++) {
        dLocal[i] = tmpDptr[ex_indices[i]];
    }
    VecRestoreArray(Dvec, &tmpDptr);

    if (OPP_IS_VALID_PROCESS) {
        const double* n_bnd_pot = (double*)(n_bnd_pot_dat->data_d);
        double* node_potential = (double*)(node_potential_dat->data_d);

        opp_mem::copy_host_to_dev<double>(dLocal_d, dLocal.data(), neq);

#ifdef USE_SYCL
        opp_queue->submit([&](sycl::handler &cgh) {
            const int n_nodes_set_l = n_nodes_set;
            const int* n2eq_map_l = node_to_eq_map_d;
            const double* dLocal_l = dLocal_d;
            const double *n_bnd_pot_l = n_bnd_pot;
            double *node_potential_l = node_potential;

            cgh.parallel_for(sycl::nd_range<1>(OPP_gpu_threads_per_block*nodesNBlocks,OPP_gpu_threads_per_block), 
                [=](sycl::nd_item<1> item) {
                    compute_node_potential_kernel(n_nodes_set_l, n2eq_map_l, dLocal_l, 
                        n_bnd_pot_l, node_potential_l, item);
                });
        });
#else
        compute_node_potential_kernel <<<nodesNBlocks, OPP_gpu_threads_per_block>>> (
            n_nodes_set, 
            node_to_eq_map_d, 
            dLocal_d,
            n_bnd_pot, 
            node_potential
        );
#endif

        OPP_DEVICE_SYNCHRONIZE();
    }
}

//*************************************************************************************************
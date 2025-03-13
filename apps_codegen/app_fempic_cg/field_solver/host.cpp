
#define DEVICE_TYPE Device_CPU

//*************************************************************************************************
void FESolver::init_device_variables() { 
    OPP_RUN_ON_ROOT() opp_printf("FESolver", "Running on Host");
}

//*************************************************************************************************
void FESolver::destroy_device_variables() {

}

//*************************************************************************************************
void FESolver::init_f1_and_J(const opp_dat ion_den_dat) 
{   
    const OPP_REAL* ion_den = opp_get_data<OPP_REAL>(ion_den_dat);
    for (int n = 0; n < n_nodes_inc_halo; n++) 
    {
        const int eq_idx = node_to_eq_map[n];
        if (eq_idx < 0) continue;    /*is this a non-boundary node?*/
        
        tempNEQ1[eq_idx] = (CONST_QE / CONST_EPS0) * ((ion_den[n]) + 
                        CONST_n0 * exp((dLocal[eq_idx] - CONST_phi0) / CONST_kTe));
    }

    for (int n = 0; n < neq; n++) {
        tempNEQ2[n] = 0;
        tempNEQ3[n] = -CONST_QE / CONST_EPS0 * 
                        CONST_n0 * exp((dLocal[n] - CONST_phi0) / CONST_kTe) * (1 / CONST_kTe);
    }
}

//*************************************************************************************************
void FESolver::build_f1_vector() 
{
    double Na = 0.0, ff = 0.0; 

    for (int e = 0; e < n_cells_inc_halo; e++) 
    {
        for (int a = 0; a < 4; a++) // get 4 nodes in tetrahedra
        {
            // IDEA : F1vec will be enriched only for the own range (owner compute),
            // hence import halo region calculation is not required here. is it correct?
            const int node_idx = c2n_map->map[e * N_PER_C + a];
            const int A = node_to_eq_map[node_idx]; 
            if (A >= 0)    /*if unknown node*/ 
            {
                ff = 0.0;

                /*perform quadrature*/
                for (int k = 0; k < n_int; k++) 
                {
                    for (int j = 0; j < n_int; j++)
                    {
                        for (int i = 0; i < n_int; i++) 
                        {
                            switch(a) // evaluate_na
                            {
                                case 0: Na = 0.5*(l[i]+1); break;
                                case 1: Na = 0.5*(l[j]+1); break;
                                case 2: Na = 0.5*(l[k]+1); break;
                                case 3: Na = 1 - 0.5*(l[i]+1) - 0.5*(l[j]+1) - 0.5*(l[k]+1); break; 
                                default: Na = 0;  
                            }

                            ff += (tempNEQ1[A] * Na * detJ[e] * W[i] * W[j] * W[k]);
                        }
                    }
                }

                ff *= (1.0/8.0);    /*change of limits*/
                f1Local[A] += ff;
            }
        }
    }

    VecSetValues(F1vec, neq, vec_col.data(), f1Local.data(), INSERT_VALUES);
    VecAssemblyBegin(F1vec); VecAssemblyEnd(F1vec);
}

//*************************************************************************************************
void FESolver::build_j_matrix() 
{     
    double Na = 0.0;

    for (int e = 0; e < n_cells_inc_halo; e++) 
    {
        for (int a = 0; a < 4; a++) // get 4 nodes in tetrahedra
        {
            double ff=0;

            const int node_idx = c2n_map->map[e * N_PER_C + a];
            const int A = node_to_eq_map[node_idx]; 
            if (A >= 0)    /*if unknown node*/ 
            {
                for (int k=0;k<n_int;k++) /*perform quadrature*/
                    for (int j=0;j<n_int;j++)
                        for (int i=0;i<n_int;i++) 
                        {   
                            switch(a) // evaluate_na
                            {
                                case 0: Na = 0.5*(l[i]+1); break;
                                case 1: Na = 0.5*(l[j]+1); break;
                                case 2: Na = 0.5*(l[k]+1); break;
                                case 3: Na = 1 - 0.5*(l[i]+1) - 0.5*(l[j]+1) - 0.5*(l[k]+1); break; 
                                default: Na = 0;  
                            }

                            ff += (tempNEQ3[A] * Na * detJ[e] * W[i] * W[j] * W[k]);
                        }

                ff *= (1.0 / 8.0);    /*change of limits*/
                tempNEQ2[A] += ff;
            }
        }
    }

    for (int u = 0; u < neq; u++) /*subtract diagonal term*/
        MatSetValue(Jmat, (u + own_start), (u + own_start), (-tempNEQ2[u]), ADD_VALUES); // J[u][u]-=tempNEQ2[u];

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

    const double* n_bnd_pot = (double*)(n_bnd_pot_dat->data);
    double* node_potential = (double*)(node_potential_dat->data);

    /* combine d (solution from linear solver) and boundary potential to get node_potential */
    for (int n = 0; n < n_nodes_set; n++) // owner compute, hence can have only upto node set size
    {
        node_potential[n] = n_bnd_pot[n]; /*zero on non-g nodes*/

        const int A = node_to_eq_map[n];
        if (A >= 0)           /*is this a non-boundary node?*/
            node_potential[n] += dLocal[A];
    }
}

//*************************************************************************************************
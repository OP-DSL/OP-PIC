
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k6 {
inline void compute_node_charge_density_kernel(
    double *node_charge_den,
    const double *node_volume
) {
    (*node_charge_den) *= (CONST_spwt[0] / node_volume[0]);
}
}

void opp_par_loop_all__compute_node_charge_density_kernel(opp_set set,
    opp_arg arg0, // n_charge_den | OPP_RW
    opp_arg arg1 // n_volume | OPP_READ
) 
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("compute_node_charge_density_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__compute_node_charge_density_kernel set_size %d", set->size);

    const int nthreads = omp_get_max_threads(); 

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
 
    opp_mpi_halo_wait_all(nargs, args);    
  
    
  
    #pragma omp parallel for
    for (int n = 0; n < iter_size; ++n) 
    {
        {
            opp_k6::compute_node_charge_density_kernel(
                (OPP_REAL *)args[0].data + (n * 1),
                (const OPP_REAL *)args[1].data + (n * 1)
            );
        }
    }
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("compute_node_charge_density_kernel");
}

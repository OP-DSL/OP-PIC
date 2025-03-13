
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k5 {
inline void deposit_charge_on_nodes_kernel(
    const double *part_lc,
    double *node_charge_den0,
    double *node_charge_den1,
    double *node_charge_den2,
    double *node_charge_den3
) {
    node_charge_den0[0] += part_lc[0];
    node_charge_den1[0] += part_lc[1];
    node_charge_den2[0] += part_lc[2];
    node_charge_den3[0] += part_lc[3];
}
}

void opp_par_loop_all__deposit_charge_on_nodes_kernel(opp_set set,
    opp_arg arg0, // p_lc | OPP_READ
    opp_arg arg1, // n_charge_den | OPP_INC
    opp_arg arg2, // n_charge_den | OPP_INC
    opp_arg arg3, // n_charge_den | OPP_INC
    opp_arg arg4 // n_charge_den | OPP_INC
) 
{
    const int nargs = 5;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = arg4;

    opp_profiler->start("deposit_charge_on_nodes_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__deposit_charge_on_nodes_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
    opp_init_double_indirect_reductions(nargs, args);
    opp_mpi_halo_wait_all(nargs, args);    
 
    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 
    for (int n = 0; n < iter_size; ++n) 
    {
        opp_p2c = OPP_mesh_relation_data + n;
        const OPP_INT *map0 = args[1].map_data + (opp_p2c[0] * 4);
   
        opp_k5::deposit_charge_on_nodes_kernel(
            (const OPP_REAL *)args[0].data + (n * 4),
            (OPP_REAL *)args[1].data + (map0[0] * 1),
            (OPP_REAL *)args[2].data + (map0[1] * 1),
            (OPP_REAL *)args[3].data + (map0[2] * 1),
            (OPP_REAL *)args[4].data + (map0[3] * 1)
        );
    }
    opp_exchange_double_indirect_reductions(nargs, args);
    opp_complete_double_indirect_reductions(nargs, args);
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("deposit_charge_on_nodes_kernel");
}

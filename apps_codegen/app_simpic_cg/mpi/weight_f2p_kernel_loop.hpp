
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k1 {
void weight_f2p_kernel(
        const double* node0_field_E,  //LHS
        const double* node1_field_E,  //RHS
        const double* particle0_position_x,
        double* particle0_field_E
    )
{
    double xx = ((particle0_position_x[0] - CONST_xl[0]) / CONST_dx[0]); // Makes Global position to local position comapared to the cell
    int n = int(xx);
    double frac = (xx - n);

    particle0_field_E[0] = ((frac * node1_field_E[0]) + ((1.0 - frac) * node0_field_E[0]));
}
}

void opp_par_loop_all__weight_f2p_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // n_field_e | OPP_READ
    opp_arg arg1, // n_field_e | OPP_READ
    opp_arg arg2, // p_pos_x | OPP_READ
    opp_arg arg3 // p_field_e | OPP_WRITE
) 
{
    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    opp_profiler->start("weight_f2p_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__weight_f2p_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);

    opp_mpi_halo_wait_all(nargs, args);    
 
    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 
    for (int n = 0; n < iter_size; ++n) 
    {
        opp_p2c = OPP_mesh_relation_data + n;
        const OPP_INT *map0 = args[0].map_data + (opp_p2c[0] * 2);
   
        opp_k1::weight_f2p_kernel(
            (const OPP_REAL *)args[0].data + (map0[0] * 1),
            (const OPP_REAL *)args[1].data + (map0[1] * 1),
            (const OPP_REAL *)args[2].data + (n * 1),
            (OPP_REAL *)args[3].data + (n * 1)
        );
    }
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("weight_f2p_kernel");
}


//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k3 {
inline void calculate_new_pos_vel_kernel(
    const double *cell_ef,
    double *part_pos,
    double *part_vel
) {
    const double coefficient1 = CONST_charge[0] / CONST_mass[0] * (CONST_dt[0]);
    for (int i = 0; i < 3; i++) {
        part_vel[i] += (coefficient1 * cell_ef[i]);
        part_pos[i] += part_vel[i] * (CONST_dt[0]);
    }
}
}

void opp_par_loop_all__calculate_new_pos_vel_kernel(opp_set set,
    opp_arg arg0, // c_ef | OPP_READ
    opp_arg arg1, // p_pos | OPP_WRITE
    opp_arg arg2 // p_vel | OPP_WRITE
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("calculate_new_pos_vel_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__calculate_new_pos_vel_kernel set_size %d", set->size);

    const int iter_size = set->size;
    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 

    for (int n = 0; n < iter_size; ++n) 
    {
        opp_p2c = OPP_mesh_relation_data + n;
   
        opp_k3::calculate_new_pos_vel_kernel(
            (const OPP_REAL *)args[0].data + (opp_p2c[0] * 3),
            (OPP_REAL *)args[1].data + (n * 3),
            (OPP_REAL *)args[2].data + (n * 3)
        );
    }

    opp_profiler->end("calculate_new_pos_vel_kernel");
}

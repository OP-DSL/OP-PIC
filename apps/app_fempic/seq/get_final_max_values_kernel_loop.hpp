
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k9 {
inline void get_final_max_values_kernel(
    const double* n_charge_den,
    double* max_n_charge_den,
    const double* n_pot,
    double* max_n_pot
) {
    *max_n_charge_den = ((abs(*n_charge_den) > *max_n_charge_den) ? (abs(*n_charge_den)) : (*max_n_charge_den));
    *max_n_pot = ((*n_pot > *max_n_pot) ? (*n_pot) : (*max_n_pot));
}
}

void opp_par_loop_all__get_final_max_values_kernel(opp_set set,
    opp_arg arg0, // n_charge_den | OPP_READ
    opp_arg arg1, // | OPP_MAX
    opp_arg arg2, // n_potential | OPP_READ
    opp_arg arg3 // | OPP_MAX
) 
{
    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    opp_profiler->start("get_final_max_values_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__get_final_max_values_kernel set_size %d", set->size);

    const int iter_size = set->size;

    for (int n = 0; n < iter_size; ++n) 
    {
        opp_k9::get_final_max_values_kernel(
            (const OPP_REAL *)args[0].data + (n * 1),
            (OPP_REAL *)args[1].data,
            (const OPP_REAL *)args[2].data + (n * 1),
            (OPP_REAL *)args[3].data
        );
    }

    opp_profiler->end("get_final_max_values_kernel");
}

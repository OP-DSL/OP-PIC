
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k4 {
void sum_laplace_kernel(
        const double* node0_xlocal,
        double* node0_field_P
    )
{
    double rv = 0.0;
    double lv = CONST_lhs_voltage[0];

    double frac = ((*node0_xlocal) / CONST_L[0]);
    (*node0_field_P) += (frac * rv + (1. - frac) * lv);
}
}

void opp_par_loop_all__sum_laplace_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // n_xlocal | OPP_READ
    opp_arg arg1 // n_field_p | OPP_RW
) 
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("sum_laplace_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__sum_laplace_kernel set_size %d", set->size);

    const int iter_size = set->size;

    for (int n = 0; n < iter_size; ++n) 
    {
        opp_k4::sum_laplace_kernel(
            (const OPP_REAL *)args[0].data + (n * 1),
            (OPP_REAL *)args[1].data + (n * 1)
        );
    }

    opp_profiler->end("sum_laplace_kernel");
}

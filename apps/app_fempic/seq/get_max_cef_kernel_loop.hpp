
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k8 {
inline void get_max_cef_kernel(
    const double* val,
    double* max_val
) {
    for (int dim = 0; dim < 3; ++dim) {
        *max_val = ((val[dim] > *max_val) ? (val[dim]) : (*max_val));
    }
}
}

void opp_par_loop_all__get_max_cef_kernel(opp_set set,
    opp_arg arg0, // c_ef | OPP_READ
    opp_arg arg1 // | OPP_MAX
) 
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("get_max_cef_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__get_max_cef_kernel set_size %d", set->size);

    const int iter_size = set->size;

    for (int n = 0; n < iter_size; ++n) 
    {
        opp_k8::get_max_cef_kernel(
            (const OPP_REAL *)args[0].data + (n * 3),
            (OPP_REAL *)args[1].data
        );
    }

    opp_profiler->end("get_max_cef_kernel");
}

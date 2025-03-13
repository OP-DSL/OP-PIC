
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k1 {
inline void update_pos_kernel(const double* p_vel, double* p_pos, int* p_mdir)
{
    for (int dm = 0; dm < 2; dm++) {

        const double offset = p_vel[dm] * CONST_dt[0];
        p_pos[dm] += offset; // s1 = s0 + ut

        // correct for periodic boundary conditions
        const int n_extent_offset_int = std::abs(p_pos[dm]) + 2.0;
        const double temp_pos = p_pos[dm] + n_extent_offset_int * CONST_extents[dm];
        p_pos[dm] = std::fmod(temp_pos, CONST_extents[dm]);

        p_mdir[dm] = (offset > 0) ? 1 : -1;
    }
}
}

void opp_par_loop_all__update_pos_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // p_vel | OPP_READ
    opp_arg arg1, // p_pos | OPP_RW
    opp_arg arg2 // p_mdir | OPP_WRITE
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("update_pos_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__update_pos_kernel set_size %d", set->size);

    const int iter_size = set->size;

    for (int n = 0; n < iter_size; ++n) 
    {
        opp_k1::update_pos_kernel(
            (const OPP_REAL *)args[0].data + (n * 2),
            (OPP_REAL *)args[1].data + (n * 2),
            (OPP_INT *)args[2].data + (n * 2)
        );
    }

    opp_profiler->end("update_pos_kernel");
}

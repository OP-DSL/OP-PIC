
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k1 {
inline void update_pos_kernel(const double* part_vel, double* p_pos)
{
    for (int dm = 0; dm < 2; dm++) {

        p_pos[dm] += part_vel[dm] * CONST_dt[0]; // s1 = s0 + ut

        // correct for periodic boundary conditions
        const int n_extent_offset_int = std::abs(p_pos[dm]) + 2.0;
        const double temp_pos = p_pos[dm] + n_extent_offset_int * CONST_extents[dm];
        p_pos[dm] = std::fmod(temp_pos, CONST_extents[dm]);
    }
}
}

void opp_par_loop_all__update_pos_kernel(opp_set set,
    opp_arg arg0, // p_vel | OPP_READ
    opp_arg arg1 // p_pos | OPP_RW
) 
{
    const int nargs = 2;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;

    opp_profiler->start("update_pos_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__update_pos_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);

    opp_mpi_halo_wait_all(nargs, args);    
 
    for (int n = 0; n < iter_size; ++n) 
    {
        opp_k1::update_pos_kernel(
            (const OPP_REAL *)args[0].data + (n * 2),
            (OPP_REAL *)args[1].data + (n * 2)
        );
    }
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("update_pos_kernel");
}

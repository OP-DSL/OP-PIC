
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k8 {
enum Dim {
    x = 0,
    y = 1,
    z = 2,
};

inline void compute_energy_kernel(
    const int* cell0_ghost,
    const double* cell_field,
    double* energy)
{
    if (cell0_ghost[0] == 0)
    {
        energy[0] += cell_field[Dim::x] * cell_field[Dim::x] +
                    cell_field[Dim::y] * cell_field[Dim::y] +
                    cell_field[Dim::z] * cell_field[Dim::z];
    }
}
}

void opp_par_loop_all__compute_energy_kernel(opp_set set,
    opp_arg arg0, // c_mask_ghost | OPP_READ
    opp_arg arg1, // c_e | OPP_READ
    opp_arg arg2 // | OPP_INC
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("compute_energy_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__compute_energy_kernel set_size %d", set->size);

    const int iter_size = set->size;

    for (int n = 0; n < iter_size; ++n) 
    {
        opp_k8::compute_energy_kernel(
            (const OPP_INT *)args[0].data + (n * 1),
            (const OPP_REAL *)args[1].data + (n * 3),
            (OPP_REAL *)args[2].data
        );
    }

    opp_profiler->end("compute_energy_kernel");
}

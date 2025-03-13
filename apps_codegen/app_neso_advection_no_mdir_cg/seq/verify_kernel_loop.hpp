
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k3 {
enum Dim {
    x = 0,
    y = 1,
};

inline void verify_kernel(
        const double* p_pos,
        const int* c_gbl_idx,
        int* incorrect_count)
{
    // get the cell boundaries for the current cell_index - using global cell index
    int ix = -1, iy = -1;
    int _ix, _iy; _ix  = ((*c_gbl_idx)); _iy  = _ix/int(CONST_ndimcells[Dim::x]); _ix -= _iy*int(CONST_ndimcells[Dim::x]); (ix) = _ix; (iy) = _iy;;

    if (ix < 0 || iy < 0)
    {
        // opp_printf("VERIFY", "Incorrect ix[%d] iy[%d] for global cell[%d] nx[%d]",
        //     ix, iy, (*c_gbl_idx), CONST_ndimcells[Dim::x]);
        (*incorrect_count)++;
        return;
    }

    // get the boundaries of that cell
    const double boundary_ll[2] = { (ix * CONST_cell_width[0]), (iy * CONST_cell_width[0]) };

    // check whether the current particle is within those boundaries or not!
    const double p_pos_x = p_pos[Dim::x];
    if (p_pos_x < boundary_ll[Dim::x] ||
            p_pos_x > (boundary_ll[Dim::x] + CONST_cell_width[0])) {

        (*incorrect_count)++;
        return;
    }

    const double p_pos_y = p_pos[Dim::y];
    if (p_pos_y < boundary_ll[Dim::y] ||
            p_pos_y > (boundary_ll[Dim::y] + CONST_cell_width[0])) {

        (*incorrect_count)++;
        return;
    }
}
}

void opp_par_loop_all__verify_kernel(opp_set set,
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1, // c_idx | OPP_READ
    opp_arg arg2 // | OPP_INC
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("verify_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__verify_kernel set_size %d", set->size);

    const int iter_size = set->size;
    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 
    OPP_INT arg2_local[1] = {0};


    for (int n = 0; n < iter_size; ++n) 
    {
        opp_p2c = OPP_mesh_relation_data + n;
   
        if (n == set->size) {
            memcpy(arg2.data, arg2_local, 1 * sizeof(OPP_INT));
        }

        opp_k3::verify_kernel(
            (const OPP_REAL *)args[0].data + (n * 2),
            (const OPP_INT *)args[1].data + (opp_p2c[0] * 1),
            (OPP_INT *)args[2].data
        );
    }

    opp_profiler->end("verify_kernel");
}

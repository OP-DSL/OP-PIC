
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

    const int nthreads = omp_get_max_threads(); 

    const int iter_size = opp_mpi_halo_exchanges(set, nargs, args);
 
    opp_mpi_halo_wait_all(nargs, args);    
    OPP_mesh_relation_data = ((OPP_INT *)set->mesh_relation_dat->data); 

    OPP_INT arg2_l[nthreads * 1];
    for (int thr = 0; thr < nthreads; thr++)
        for (int d = 0; d < 1; d++)
            arg2_l[1 * thr + d] = OPP_INT_ZERO;

  
    
  
 
    #pragma omp parallel for 
    for (int thr = 0; thr < nthreads; thr++)
    {
        const size_t start  = ((size_t)iter_size * thr) / nthreads;
        const size_t finish = ((size_t)iter_size * (thr+1)) / nthreads;
      
        for (size_t n = start; n < finish; n++)
        { 
            OPP_INT* opp_p2c = OPP_mesh_relation_data + n;
   
            opp_k3::verify_kernel(
                (const OPP_REAL *)args[0].data + (n * 2),
                (const OPP_INT *)args[1].data + (opp_p2c[0] * 1),
                arg2_l + (1 * thr)
            );
        }
    }
    
    for (int thr = 0; thr < nthreads; thr++) 
        for (int d = 0; d < 1; d++)
            ((OPP_INT*)args[2].data)[d] += (arg2_l[1 * thr + d]);
#ifdef USE_MPI 
    opp_mpi_reduce(&args[2], (OPP_INT *)args[2].data);
#endif
 
    opp_set_dirtybit(nargs, args);

    opp_profiler->end("verify_kernel");
}


//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k4_dat0_stride = -1;
OPP_INT opp_k4_dat1_stride = -1;
OPP_INT opp_k4_dat2_stride = -1;
OPP_INT opp_k4_dat3_stride = -1;
OPP_INT opp_k4_c2c_map_stride = -1;

OPP_INT* opp_k4_dat0_stride_s = nullptr;
OPP_INT* opp_k4_dat1_stride_s = nullptr;
OPP_INT* opp_k4_dat2_stride_s = nullptr;
OPP_INT* opp_k4_dat3_stride_s = nullptr;
OPP_INT* opp_k4_c2c_map_stride_s = nullptr;

//--------------------------------------------------------------
void opp_particle_move__move_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0,   // p_pos | OPP_READ
    opp_arg arg1,   // p_lc | OPP_WRITE
    opp_arg arg2,   // c_volume | OPP_READ
    opp_arg arg3   // c_det | OPP_READ
) 
{
    if (OPP_DBG) opp_printf("APP", "opp_particle_move__move_kernel set_size %d", set->size);

    opp_profiler->start("move_kernel");

    const int nargs = 5;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

    opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
    const opp_set c_set = c2c_map->from;
    const OPP_INT c2c_stride = c_set->size + c_set->exec_size + c_set->nonexec_size;

    opp_set_stride(cells_set_size_s, cells_set_size, set->cells_set->size);
    opp_set_stride(opp_k4_c2c_map_stride_s, opp_k4_c2c_map_stride, c2c_stride);

    opp_mpi_halo_wait_all(nargs, args);

#ifdef OPP_BLOCK_SIZE_4
    const int block_size = OPP_BLOCK_SIZE_4;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    int num_blocks = 200;

    do {
        opp_set_stride(comm_iteration_s, comm_iteration, OPP_comm_iteration);
        opp_set_stride(opp_k4_dat0_stride_s, opp_k4_dat0_stride, args[0].dat->set->set_capacity);
        opp_set_stride(opp_k4_dat1_stride_s, opp_k4_dat1_stride, args[1].dat->set->set_capacity);
        opp_set_stride(opp_k4_dat2_stride_s, opp_k4_dat2_stride, args[2].dat->set->set_capacity);
        opp_set_stride(opp_k4_dat3_stride_s, opp_k4_dat3_stride, args[3].dat->set->set_capacity);

        opp_init_particle_move(set, nargs, args);

        num_blocks = (OPP_iter_end - OPP_iter_start - 1) / block_size + 1;

        opp_queue->submit([&](sycl::handler &cgh) {
            
            const OPP_INT* comm_iteration_sycl = comm_iteration_s;
            const OPP_INT* opp_cell_set_size_sycl = cells_set_size_s;

            OPP_INT *remove_count = (OPP_INT *)set->particle_remove_count_d;
            OPP_INT *remove_part_indices = (OPP_INT *)OPP_remove_particle_indices_d;
            OPP_INT *move_part_indices = (OPP_INT *)OPP_move_particle_indices_d;
            OPP_INT *move_cell_indices = (OPP_INT *)OPP_move_cell_indices_d;
            OPP_INT *move_count = (OPP_INT *)OPP_move_count_d;

            const OPP_INT* opp_k4_c2c_map_stride_sycl = opp_k4_c2c_map_stride_s;
            const OPP_INT* opp_k4_dat0_stride_sycl = opp_k4_dat0_stride_s;
            const OPP_INT* opp_k4_dat1_stride_sycl = opp_k4_dat1_stride_s;
            const OPP_INT* opp_k4_dat2_stride_sycl = opp_k4_dat2_stride_s;
            const OPP_INT* opp_k4_dat3_stride_sycl = opp_k4_dat3_stride_s;
   

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // p_pos
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // p_lc
            OPP_REAL* dat2_sycl = (OPP_REAL*)args[2].data_d;     // c_volume
            OPP_REAL* dat3_sycl = (OPP_REAL*)args[3].data_d;     // c_det
            
            OPP_INT *p2c_map_sycl = (OPP_INT *)p2c_map->p2c_dat->data_d;
            const OPP_INT *c2c_map_sycl = (OPP_INT *)c2c_map->map_d; 

            const OPP_INT iter_start = OPP_iter_start;
            const OPP_INT iter_end = OPP_iter_end; 

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            auto  move_kernel_sycl = [=](
                char& opp_move_status_flag, const bool opp_move_hop_iter_one_flag, // Added by code-gen
                const OPP_INT* opp_c2c, OPP_INT* opp_p2c, // Added by code-gen
                const double *point_pos, double* point_lc,
                const double *cell_volume, const double *cell_det
            ) {

                const double coefficient2 = (1.0 / 6.0) / (*cell_volume);

                for (int i=0; i<4; i++) { /*loop over vertices*/

                    point_lc[(i) * opp_k4_dat1_stride_sycl[0]] = coefficient2 * (
                        cell_det[(i * 4 + 0) * opp_k4_dat3_stride_sycl[0]] -
                        cell_det[(i * 4 + 1) * opp_k4_dat3_stride_sycl[0]] * point_pos[(0) * opp_k4_dat0_stride_sycl[0]] +
                        cell_det[(i * 4 + 2) * opp_k4_dat3_stride_sycl[0]] * point_pos[(1) * opp_k4_dat0_stride_sycl[0]] -
                        cell_det[(i * 4 + 3) * opp_k4_dat3_stride_sycl[0]] * point_pos[(2) * opp_k4_dat0_stride_sycl[0]]);
                }

                if (!(point_lc[(0) * opp_k4_dat1_stride_sycl[0]] < 0.0 ||
                    point_lc[(0) * opp_k4_dat1_stride_sycl[0]] > 1.0 ||
                    point_lc[(1) * opp_k4_dat1_stride_sycl[0]] < 0.0 ||
                    point_lc[(1) * opp_k4_dat1_stride_sycl[0]] > 1.0 ||
                    point_lc[(2) * opp_k4_dat1_stride_sycl[0]] < 0.0 ||
                    point_lc[(2) * opp_k4_dat1_stride_sycl[0]] > 1.0 ||
                    point_lc[(3) * opp_k4_dat1_stride_sycl[0]] < 0.0 ||
                    point_lc[(3) * opp_k4_dat1_stride_sycl[0]] > 1.0)) {

                    { opp_move_status_flag = OPP_MOVE_DONE; };
                    return;
                }

                // outside the last known cell, find most negative weight and
                // use that cell_index to reduce computations
                int min_i = 0;
                double min_lc = point_lc[(0) * opp_k4_dat1_stride_sycl[0]];

                for (int i=1; i<4; i++) {
                    if (point_lc[(i) * opp_k4_dat1_stride_sycl[0]] < min_lc) {
                        min_lc = point_lc[(i) * opp_k4_dat1_stride_sycl[0]];
                        min_i = i;
                    }
                }

                if (opp_c2c[(min_i) * opp_k4_c2c_map_stride_sycl[0]] >= 0) { // is there a neighbor in this direction?
                    (*opp_p2c) = opp_c2c[(min_i) * opp_k4_c2c_map_stride_sycl[0]];
                    { opp_move_status_flag = OPP_NEED_MOVE; };
                }
                else {
                    (*opp_p2c) = INT_MAX;
                    { opp_move_status_flag = OPP_NEED_REMOVE; };
                }
            };

            // -----------------------------------------------------------------------------------------
            auto opp_part_check_status = 
                [=](char& move_flag, bool& iter_flag, int* c_idx, int p_idx) -> bool {
                
                iter_flag = false;
                if (move_flag == OPP_MOVE_DONE) {
                    return false;
                }
                else if (move_flag == OPP_NEED_REMOVE) {
                    c_idx[0] = MAX_CELL_INDEX;
                    const int removeIdx = opp_atomic_fetch_add(remove_count, 1);
                    remove_part_indices[removeIdx] = p_idx;

                    return false;
                }
                else if (c_idx[0] >= opp_cell_set_size_sycl[0]) {
                    // cell_id is not owned by the current mpi rank, need to communicate
                    const int moveIdx = opp_atomic_fetch_add(move_count, 1);
                    move_part_indices[moveIdx] = p_idx;
                    move_cell_indices[moveIdx] = c_idx[0];

                    // To be removed from the current rank, packing will be done prior exchange & removal
                    const int removeIdx = opp_atomic_fetch_add(remove_count, 1);
                    remove_part_indices[removeIdx] = p_idx;

                    return false;
                }
                return true; // cell_id is an own cell and move_flag == OPP_NEED_MOVE
            };

            // -----------------------------------------------------------------------------------------
            auto opp_move_kernel = [=](sycl::nd_item<1> item) {
                
                const int tid = item.get_global_linear_id();
                const int n = tid + iter_start;

                if (n < iter_end) {
                    OPP_INT *opp_p2c = (p2c_map_sycl + n);
                    char move_flag = OPP_NEED_MOVE;
                    bool iter_one_flag = (comm_iteration_sycl[0] > 0) ? false : true;

                    
                    do {
                        const OPP_INT p2c = opp_p2c[0];
                        const OPP_INT* opp_c2c = c2c_map_sycl + p2c;

                        move_kernel_sycl(
                            move_flag, iter_one_flag, opp_c2c, opp_p2c,
                            dat0_sycl + n, // p_pos 
                            dat1_sycl + n, // p_lc 
                            dat2_sycl + p2c, // c_volume 
                            dat3_sycl + p2c // c_det 
             
                        ); 
                  
                    } while (opp_part_check_status(move_flag, iter_one_flag, opp_p2c, n));
                }        
            };

            // -----------------------------------------------------------------------------------------
            cgh.parallel_for<class opp_particle_move>(
                    sycl::nd_range<1>(block_size * num_blocks, block_size), opp_move_kernel);
        });
    
    } while (opp_finalize_particle_move(set)); // MPI communication iteration

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    opp_queue->wait();
 
    opp_profiler->end("move_kernel");
}
void opp_init_direct_hop_cg(double grid_spacing, int dim, const opp_dat c_gbl_id, const opp::BoundingBox& b_box, 
    opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1, // p_lc | OPP_WRITE
    opp_arg arg2, // c_volume | OPP_READ
    opp_arg arg3 // c_det | OPP_READ
) {
    opp_profiler->start("Setup_Mover");
    
    opp_profiler->end("Setup_Mover");
}

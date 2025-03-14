
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k2_dat0_stride = -1;
OPP_INT opp_k2_dat1_stride = -1;
OPP_INT opp_k2_dat2_stride = -1;
OPP_INT opp_k2_c2c_map_stride = -1;

OPP_INT* opp_k2_dat0_stride_s = nullptr;
OPP_INT* opp_k2_dat1_stride_s = nullptr;
OPP_INT* opp_k2_dat2_stride_s = nullptr;
OPP_INT* opp_k2_c2c_map_stride_s = nullptr;

//--------------------------------------------------------------
void opp_particle_move__move_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0,   // p_field_e | OPP_READ
    opp_arg arg1,   // p_vel_x | OPP_RW
    opp_arg arg2   // p_pos_x | OPP_RW
) 
{
    if (OPP_DBG) opp_printf("APP", "opp_particle_move__move_kernel set_size %d", set->size);

    opp_profiler->start("move_kernel");

    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

    opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
    const opp_set c_set = c2c_map->from;
    const OPP_INT c2c_stride = c_set->size + c_set->exec_size + c_set->nonexec_size;

    opp_set_stride(cells_set_size_s, cells_set_size, set->cells_set->size);
    opp_set_stride(opp_k2_c2c_map_stride_s, opp_k2_c2c_map_stride, c2c_stride);

    opp_mpi_halo_wait_all(nargs, args);

#ifdef OPP_BLOCK_SIZE_2
    const int block_size = OPP_BLOCK_SIZE_2;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    int num_blocks = 200;

    do {
        opp_set_stride(comm_iteration_s, comm_iteration, OPP_comm_iteration);
        opp_set_stride(opp_k2_dat0_stride_s, opp_k2_dat0_stride, args[0].dat->set->set_capacity);
        opp_set_stride(opp_k2_dat1_stride_s, opp_k2_dat1_stride, args[1].dat->set->set_capacity);
        opp_set_stride(opp_k2_dat2_stride_s, opp_k2_dat2_stride, args[2].dat->set->set_capacity);

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

            const OPP_INT* opp_k2_c2c_map_stride_sycl = opp_k2_c2c_map_stride_s;
            const OPP_INT* opp_k2_dat0_stride_sycl = opp_k2_dat0_stride_s;
            const OPP_INT* opp_k2_dat1_stride_sycl = opp_k2_dat1_stride_s;
            const OPP_INT* opp_k2_dat2_stride_sycl = opp_k2_dat2_stride_s;
   
            const OPP_REAL* CONST_xl_sycl = CONST_xl_s;
            const OPP_INT* CONST_neighbour_cell_sycl = CONST_neighbour_cell_s;
            const OPP_REAL* CONST_xr_sycl = CONST_xr_s;
            const OPP_REAL* CONST_qm_sycl = CONST_qm_s;
            const OPP_REAL* CONST_dt_sycl = CONST_dt_s;
            const OPP_INT* CONST_rank_sycl = CONST_rank_s;
            const OPP_INT* CONST_comm_size_sycl = CONST_comm_size_s;
            const OPP_REAL* CONST_dx_sycl = CONST_dx_s;

            OPP_REAL* dat0_sycl = (OPP_REAL*)args[0].data_d;     // p_field_e
            OPP_REAL* dat1_sycl = (OPP_REAL*)args[1].data_d;     // p_vel_x
            OPP_REAL* dat2_sycl = (OPP_REAL*)args[2].data_d;     // p_pos_x
            
            OPP_INT *p2c_map_sycl = (OPP_INT *)p2c_map->p2c_dat->data_d;
            const OPP_INT *c2c_map_sycl = (OPP_INT *)c2c_map->map_d; 

            const OPP_INT iter_start = OPP_iter_start;
            const OPP_INT iter_end = OPP_iter_end; 

            // user provided elemental kernel
            // -----------------------------------------------------------------------------------------
            enum Dir : char {
                Left = 0,
                Right,
            };

             move_kernel_sycl = [=](
                    char& opp_move_status_flag, const bool opp_move_hop_iter_one_flag, // Added by code-gen
                const OPP_INT* opp_c2c, OPP_INT* opp_p2c, // Added by code-gen
                const double* part_field_E,
                    double* part_velocity_x,
                    double* part_position_x
                )
            {
                if ((opp_move_hop_iter_one_flag))
                {
                    part_velocity_x[(0) * opp_k2_dat1_stride_sycl[0]] += (CONST_qm_sycl[0] * part_field_E[(0) * opp_k2_dat0_stride_sycl[0]]);
                    part_position_x[(0) * opp_k2_dat2_stride_sycl[0]] += (part_velocity_x[(0) * opp_k2_dat1_stride_sycl[0]] * CONST_dt_sycl[0]);
                }

                // since particle cell index can be easily calculated with global positioning, no need to search by iterating
                if ((part_position_x[(0) * opp_k2_dat2_stride_sycl[0]] > CONST_xl_sycl[0]) && (part_position_x[(0) * opp_k2_dat2_stride_sycl[0]] < CONST_xr_sycl[0]))
                {
                    double xx = ((part_position_x[(0) * opp_k2_dat2_stride_sycl[0]] - CONST_xl_sycl[0]) / CONST_dx_sycl[0]); // Makes Global position to local position comapared to the cell
                    opp_p2c[0] = int(xx);
                    { opp_move_status_flag = OPP_MOVE_DONE; };
                }
                else if ((part_position_x[(0) * opp_k2_dat2_stride_sycl[0]] >= CONST_xl_sycl[0]) && (CONST_rank_sycl[0] == 0) ||
                         (part_position_x[(0) * opp_k2_dat2_stride_sycl[0]] <= CONST_xr_sycl[0]) && (CONST_rank_sycl[0] == (CONST_comm_size_sycl[0]-1)))
                {
                    opp_p2c[0] = INT_MAX;
                    { opp_move_status_flag = OPP_NEED_REMOVE; };
                }
                else
                {
                    opp_p2c[0] = (part_position_x[(0) * opp_k2_dat2_stride_sycl[0]] < CONST_xl_sycl[0]) ?
                                    CONST_neighbour_cell_sycl[Dir::Left] : CONST_neighbour_cell_sycl[Dir::Right];
                    { opp_move_status_flag = OPP_NEED_MOVE; };
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
                            dat0_sycl + n, // p_field_e 
                            dat1_sycl + n, // p_vel_x 
                            dat2_sycl + n // p_pos_x 
             
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
    OPP_DEVICE_SYNCHRONIZE();
 
    opp_profiler->end("move_kernel");
}

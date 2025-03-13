
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k2 {
enum Dir : char {
    Left = 0,
    Right,
};

void move_kernel(
        const double* part_field_E,
        double* part_velocity_x,
        double* part_position_x
    )
{
    if ((opp_move_hop_iter_one_flag))
    {
        part_velocity_x[0] += (CONST_qm[0] * part_field_E[0]);
        part_position_x[0] += (part_velocity_x[0] * CONST_dt[0]);
    }

    // since particle cell index can be easily calculated with global positioning, no need to search by iterating
    if ((part_position_x[0] > CONST_xl[0]) && (part_position_x[0] < CONST_xr[0]))
    {
        double xx = ((part_position_x[0] - CONST_xl[0]) / CONST_dx[0]); // Makes Global position to local position comapared to the cell
        opp_p2c[0] = int(xx);
        { opp_move_status_flag = OPP_MOVE_DONE; };
    }
    else if ((part_position_x[0] >= CONST_xl[0]) && (CONST_rank[0] == 0) ||
             (part_position_x[0] <= CONST_xr[0]) && (CONST_rank[0] == (CONST_comm_size[0]-1)))
    {
        opp_p2c[0] = INT_MAX;
        { opp_move_status_flag = OPP_NEED_REMOVE; };
    }
    else
    {
        opp_p2c[0] = (part_position_x[0] < CONST_xl[0]) ?
                        CONST_neighbour_cell[Dir::Left] : CONST_neighbour_cell[Dir::Right];
        { opp_move_status_flag = OPP_NEED_MOVE; };
    }
}
}

void opp_particle_move__move_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // p_field_e | OPP_READ
    opp_arg arg1, // p_vel_x | OPP_RW
    opp_arg arg2 // p_pos_x | OPP_RW
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

    OPP_mesh_relation_data = (OPP_INT*)p2c_map->p2c_dat->data;

    opp_mpi_halo_exchanges(set, nargs, args);
        
    opp_mpi_halo_wait_all(nargs, args);

    // lambda function for multi hop particle movement
    auto multihop_mover = [&](const int n) {

        opp_p2c = OPP_mesh_relation_data + n;

        if (opp_p2c[0] == MAX_CELL_INDEX) {
            return;
        }

        opp_move_status_flag = OPPX_MOVE_DONE; 
        opp_move_hop_iter_one_flag = true;

        do {
            opp_c2c = &((c2c_map->map)[opp_p2c[0] * 2]);

            opp_k2::move_kernel(
                (const OPP_REAL *)args[0].data + (n * 1),
                (OPP_REAL *)args[1].data + (n * 1),
                (OPP_REAL *)args[2].data + (n * 1)
            );

        } while (opp_check_part_move_status(opp_p2c[0], n, set->particle_remove_count));
    };

    // ----------------------------------------------------------------------------
    opp_init_particle_move(set, nargs, args);


    opp_profiler->start("Mv_AllMv0");

    // ----------------------------------------------------------------------------
    // check whether all particles not marked for global comm is within cell, 
    // and if not mark to move between cells within the MPI rank, mark for neighbour comm
    opp_profiler->start("move_kernel_only");
    const int start1 = OPP_iter_start;
    const int end1 = OPP_iter_end;
    for (int i = start1; i < end1; i++) { 
        
        multihop_mover(i);
    }
    opp_profiler->end("move_kernel_only");

    opp_profiler->end("Mv_AllMv0");


    // ----------------------------------------------------------------------------
    // Do neighbour communication and if atleast one particle is received by the currect rank, 
    // then iterate over the newly added particles
    while (opp_finalize_particle_move(set)) {
        
        const std::string profName = std::string("Mv_AllMv") + std::to_string(OPP_comm_iteration);
        opp_profiler->start(profName);
        
        opp_init_particle_move(set, nargs, args);

        // check whether particle is within cell, and if not move between cells within the MPI rank, mark for neighbour comm
        opp_profiler->start("move_kernel_only");
        const int start3 = OPP_iter_start;
        const int end3 = OPP_iter_end;
        for (int i = start3; i < end3; i++) {      

            multihop_mover(i);
        }
        opp_profiler->end("move_kernel_only");

        opp_profiler->end(profName);
    }

    opp_set_dirtybit(nargs, args);

    opp_profiler->end("move_kernel");
}

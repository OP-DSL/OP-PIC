
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k2 {
enum CellMap {
    xd_y = 0,
    xu_y,
    x_yd,
    x_yu
};

enum Dim {
    x = 0,
    y = 1,
};

inline void move_kernel(const double* p_pos, const double* c_pos_ll)
{
    // check for x direction movement
    const double p_pos_x = p_pos[Dim::x];
    if (p_pos_x < c_pos_ll[Dim::x]) {
        opp_p2c[0] = opp_c2c[CellMap::xd_y];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }
    if (p_pos_x > (c_pos_ll[Dim::x] + CONST_cell_width[0])) {
        opp_p2c[0] = opp_c2c[CellMap::xu_y];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }

    // check for y direction movement
    const double p_pos_y = p_pos[Dim::y];
    if (p_pos_y < c_pos_ll[Dim::y]) {
        opp_p2c[0] = opp_c2c[CellMap::x_yd];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }
    if (p_pos_y > (c_pos_ll[Dim::y] + CONST_cell_width[0])) {
        opp_p2c[0] = opp_c2c[CellMap::x_yu];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }

    { opp_move_status_flag = OPP_MOVE_DONE; };
}
}

void opp_particle_move__move_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1 // c_pos_ll | OPP_READ
) 
{
    if (OPP_DBG) opp_printf("APP", "opp_particle_move__move_kernel set_size %d", set->size);

    opp_profiler->start("move_kernel");

    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

    OPP_mesh_relation_data = (OPP_INT*)p2c_map->p2c_dat->data;
#ifdef LOG_HOPS
    OPP_move_max_hops = 0;
#endif

    opp_mpi_halo_exchanges(set, nargs, args);
        
    opp_mpi_halo_wait_all(nargs, args);

    // lambda function for multi hop particle movement
    auto multihop_mover = [&](const int n) {

        opp_p2c = OPP_mesh_relation_data + n;

        if (opp_p2c[0] == MAX_CELL_INDEX) {
            return;
        }

        opp_move_status_flag = OPP_MOVE_DONE; 
        opp_move_hop_iter_one_flag = true;

#ifdef LOG_HOPS
        int hops = 0;
#endif
        do {
            opp_c2c = &((c2c_map->map)[opp_p2c[0] * 4]);

            opp_k2::move_kernel(
                (const OPP_REAL *)args[0].data + (n * 2),
                (const OPP_REAL *)args[1].data + (opp_p2c[0] * 2)
            );
#ifdef LOG_HOPS
            hops++;
#endif
        } while (opp_check_part_move_status(opp_p2c[0], n, set->particle_remove_count));

#ifdef LOG_HOPS
        OPP_move_max_hops = (OPP_move_max_hops < hops) ? hops : OPP_move_max_hops;
#endif    
    };

    // ----------------------------------------------------------------------------
    opp_init_particle_move(set, nargs, args);

    if (useGlobalMove) {
        
        globalMover->initGlobalMove();

        opp_profiler->start("GblMv_Move");

        // check whether particles needs to be moved over global move routine
        const int start = OPP_iter_start;
        const int end = OPP_iter_end;
        for (int i = start; i < end; i++) {       

            // TODO: we assume pos is in arg 0, Is there a better way?
            const opp_point* point = (const opp_point*)&(((OPP_REAL*)args[0].data)[i * 2]); 

            // check for global move, and if satisfy global move criteria, then remove the particle from current rank
            opp_part_checkForGlobalMove2D(set, *point, i, OPP_mesh_relation_data[i]);
        }

        opp_profiler->end("GblMv_Move");

        globalMover->communicate(set);
    }

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
    // finalize the global move routine and iterate over newly added particles and check whether they need neighbour comm
    if (useGlobalMove && globalMover->finalize(set) > 0) {
        
        opp_profiler->start("GblMv_AllMv");

        // need to change arg data since particle resize in globalMover::finalize could change the pointer in dat->data 
        for (int i = 0; i < nargs; i++)
            if (args[i].argtype == OPP_ARG_DAT && args[i].dat->set->is_particle)
                args[i].data = args[i].dat->data;

        // check whether the new particle is within cell, and if not move between cells within the MPI rank, 
        // mark for neighbour comm. Do only for the globally moved particles 
        opp_profiler->start("move_kernel_only");
        const int start2 = (set->size - set->diff);
        const int end2 = set->size;
        for (int i = start2; i < end2; i++) { 
                
            multihop_mover(i);                 
        }
        opp_profiler->end("move_kernel_only");

        opp_profiler->end("GblMv_AllMv");
    }

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

void opp_init_direct_hop_cg(double grid_spacing, const opp_dat c_gbl_id, const opp::BoundingBox& b_box, 
    opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1 // c_pos_ll | OPP_READ
) {
    opp_profiler->start("Setup_Mover");
    
    useGlobalMove = opp_params->get<OPP_BOOL>("opp_global_move");

    if (useGlobalMove) {

        const int nargs = 3;
        opp_arg args[nargs];

        args[0] = arg0;
        args[1] = arg1;
        args[2] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

#ifdef USE_MPI
        opp_mpi_halo_exchanges(c_gbl_id->set, nargs, args);

        comm = std::make_shared<opp::Comm>(MPI_COMM_WORLD);
        globalMover = std::make_unique<opp::GlobalParticleMover>(comm->comm_parent);

        opp_mpi_halo_wait_all(nargs, args);
#endif

        boundingBox = std::make_shared<opp::BoundingBox>(b_box);
        cellMapper = std::make_shared<opp::CellMapper>(boundingBox, grid_spacing, comm);

        const int c_set_size = c_gbl_id->set->size;

        // lambda function for dh mesh search loop
        auto all_cell_checker = [&](const opp_point& point, int& cid) {          

            for (int ci = 0; ci < c_set_size; ++ci) {
                opp_move_status_flag = OPP_NEED_MOVE;  
                opp_move_hop_iter_one_flag = true;

                int temp_ci = ci; // we dont want to get iterating ci changed within the kernel, hence get a copy
                
                opp_p2c = &(temp_ci);           
                opp_c2c = &((c2c_map->map)[temp_ci * 4]);

                opp_k2::move_kernel(
                    (const OPP_REAL*)&point,
                    (const OPP_REAL *)args[1].data + (temp_ci * 2) // c_pos_ll| OPP_READ
                );
                if (opp_move_status_flag == OPP_MOVE_DONE) {       
                    cid = temp_ci;
                    break;
                }
            }
        };
        
        if (opp_params->get<OPP_BOOL>("opp_dh_data_generate")) {
            cellMapper->generateStructuredMesh(c_gbl_id->set, c_gbl_id, all_cell_checker);
        }
        else {
            cellMapper->generateStructuredMeshFromFile(c_gbl_id->set, c_gbl_id);  
        }
        
        opp_profiler->reg("GlbToLocal");
        opp_profiler->reg("GblMv_Move");
        opp_profiler->reg("GblMv_AllMv");
    }

    opp_profiler->end("Setup_Mover");
}


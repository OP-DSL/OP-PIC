
//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k2 {
enum Dim {
    x = 0,
    y = 1,
};

enum CellMap {
    xd_y = 0,
    xu_y,
    x_yd,
    x_yu
};

inline void move_kernel(char& opp_move_status_flag, const bool opp_move_hop_iter_one_flag, // Added by code-gen
    const OPP_INT* opp_c2c, OPP_INT* opp_p2c, // Added by code-gen
    const double* p_pos, int* p_mdir, const double* c_pos_ll)
{
    // check for x direction movement
    const double p_pos_x_diff = (p_pos[Dim::x] - c_pos_ll[Dim::x]);
    if ((p_pos_x_diff >= 0.0) && (p_pos_x_diff <= CONST_cell_width[0])) {
        p_mdir[Dim::x] = 0; // within cell in x direction
    }
    else if (p_mdir[Dim::x] > 0) {
        opp_p2c[0] = opp_c2c[CellMap::xu_y];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }
    else if (p_mdir[Dim::x] < 0) {
        opp_p2c[0] = opp_c2c[CellMap::xd_y];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }

    // check for y direction movement
    const double p_pos_y_diff = (p_pos[Dim::y] - c_pos_ll[Dim::y]);
    if ((p_pos_y_diff >= 0.0) && (p_pos_y_diff <= CONST_cell_width[0])) {
        p_mdir[Dim::y] = 0; // within cell in y direction
    }
    else if (p_mdir[Dim::y] > 0) {
        opp_p2c[0] = opp_c2c[CellMap::x_yu];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }
    else if (p_mdir[Dim::y] < 0) {
        opp_p2c[0] = opp_c2c[CellMap::x_yd];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }

    { opp_move_status_flag = OPP_MOVE_DONE; };
}
}

void opp_particle_move__move_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1, // p_mdir | OPP_RW
    opp_arg arg2 // c_pos_ll | OPP_READ
) 
{
    if (OPP_DBG) opp_printf("APP", "opp_particle_move__move_kernel set_size %d", set->size);

    opp_profiler->start("move_kernel");

    const int nthreads = omp_get_max_threads();

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
    auto multihop_mover = [&](const int n, const int thread) {

        OPP_INT *opp_p2c = OPP_mesh_relation_data + n;

        if (opp_p2c[0] == MAX_CELL_INDEX) {
            return;
        }

        OPP_INT *opp_c2c = nullptr;
        char move_flag = OPP_MOVE_DONE;
        bool iter_one_flag = true;

        do {
            move_flag = OPP_MOVE_DONE;
            opp_c2c = c2c_map->map + (opp_p2c[0] * 4);

            opp_k2::move_kernel(
                move_flag, iter_one_flag, opp_c2c, opp_p2c, 
                (const OPP_REAL *)args[0].data + (n * 2),
                (OPP_INT *)args[1].data + (n * 2),
                (const OPP_REAL *)args[2].data + (opp_p2c[0] * 2)
            );

        } while (opp_check_part_move_status(move_flag, iter_one_flag, opp_p2c[0], n, thread));
    };

    // ----------------------------------------------------------------------------
    opp_init_particle_move(set, nargs, args);
    const int total_count = OPP_iter_end - OPP_iter_start;

    if (useGlobalMove) { // For now Global Move will not work with MPI

#ifdef USE_MPI         
        globalMover->initGlobalMove();
#endif

        opp_profiler->start("GblMv_Move");

        // check whether particles needs to be moved over global move routine
        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++)
        {
            const size_t start  = ((size_t)total_count * thr) / nthreads;
            const size_t finish = ((size_t)total_count * (thr+1)) / nthreads;

            for (size_t i = start; i < finish; i++)
            {
                OPP_INT* opp_p2c = OPP_mesh_relation_data + i;
                // TODO: we assume pos is in arg 0, Change this!
                const opp_point* point = (const opp_point*)&(((OPP_REAL*)args[0].data)[i * 2]); 

                // check for global move, and if satisfy global move criteria, then remove the particle from current rank
                if (opp_part_checkForGlobalMove2D(set, *point, i, opp_p2c[0], thr)) { 
                    
                    part_remove_count_per_thr[thr] += 1;
                    continue;  
                }
            }
        }

        opp_profiler->end("GblMv_Move");

#ifdef USE_MPI 
        opp_gather_gbl_move_indices();
        globalMover->communicate(set);
#endif
    }

    opp_profiler->start("Mv_AllMv0");

    // ----------------------------------------------------------------------------
    // check whether all particles not marked for global comm is within cell, 
    // and if not mark to move between cells within the MPI rank, mark for neighbour comm
    opp_profiler->start("move_kernel_only");
    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        const size_t start  = OPP_iter_start + ((size_t)total_count * thr) / nthreads;
        const size_t finish = OPP_iter_start + ((size_t)total_count * (thr+1)) / nthreads;

        for (size_t i = start; i < finish; i++)
        {   
            multihop_mover(i, thr);
        }
    }
    opp_profiler->end("move_kernel_only");

    opp_profiler->end("Mv_AllMv0");

#ifdef USE_MPI 
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
        const int iter_start = (set->size - set->diff);
        const int iter_count = set->diff;
        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++)
        {
            const size_t start  = iter_start + ((size_t)iter_count * thr) / nthreads;
            const size_t finish = iter_start + ((size_t)iter_count * (thr+1)) / nthreads;

            for (size_t i = start; i < finish; i++)
            {   
                multihop_mover(i, thr);
            }        
        }
        opp_profiler->end("move_kernel_only");

        opp_profiler->end("GblMv_AllMv");
    }
#endif

    // ----------------------------------------------------------------------------
    // Do neighbour communication and if atleast one particle is received by the currect rank, 
    // then iterate over the newly added particles
    while (opp_finalize_particle_move(set)) {
        
        const std::string profName = std::string("Mv_AllMv") + std::to_string(OPP_comm_iteration);
        opp_profiler->start(profName);
        
        opp_init_particle_move(set, nargs, args);

        // check whether particle is within cell, and if not move between cells within the MPI rank, mark for neighbour comm
        opp_profiler->start("move_kernel_only");
        const int iter_count = (OPP_iter_end - OPP_iter_start);
        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++)
        {
            const size_t start  = OPP_iter_start + ((size_t)iter_count * thr) / nthreads;
            const size_t finish = OPP_iter_start + ((size_t)iter_count * (thr+1)) / nthreads;

            for (size_t i = start; i < finish; i++)
            {   
                multihop_mover(i, thr);
            }
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
    opp_arg arg1, // p_mdir | OPP_RW
    opp_arg arg2 // c_pos_ll | OPP_READ
) {
    opp_profiler->start("Setup_Mover");
    
    useGlobalMove = opp_params->get<OPP_BOOL>("opp_global_move");

    if (useGlobalMove) {

        const int nargs = 4;
        opp_arg args[nargs];

        args[0] = arg0;
        args[1] = arg1;
        args[2] = arg2;
        args[3] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

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
 
            // we dont want to change the original arrays during dh mesh generation, hence duplicate except OPP_READ
            OPP_INT arg1_temp[2];
            for (int ci = 0; ci < c_set_size; ++ci) {
                char opp_move_status_flag = OPP_NEED_MOVE;  
                bool opp_move_hop_iter_one_flag = true;

                int temp_ci = ci; // we dont want to get iterating ci changed within the kernel, hence get a copy
                
                OPP_INT* opp_p2c = &(temp_ci);           
                OPP_INT* opp_c2c = &((c2c_map->map)[temp_ci * 4]);
                // arg1 is OPP_RW, hence get a copy just incase
                std::memcpy(&arg1_temp, (OPP_INT *)args[1].data, (sizeof(OPP_INT) * 2));

                opp_k2::move_kernel(
                    opp_move_status_flag, opp_move_hop_iter_one_flag, opp_c2c, opp_p2c, 
                    (const OPP_REAL*)&point,
                    arg1_temp, // p_mdir| OPP_RW
                    (const OPP_REAL *)args[2].data + (temp_ci * 2) // c_pos_ll| OPP_READ
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



//*********************************************
// AUTO GENERATED CODE
//*********************************************

namespace opp_k4 {
inline void move_kernel(
    char& opp_move_status_flag, const bool opp_move_hop_iter_one_flag, // Added by code-gen
    const OPP_INT* opp_c2c, OPP_INT* opp_p2c, // Added by code-gen
    const double *point_pos, double* point_lc,
    const double *cell_volume, const double *cell_det
) {

    const double coefficient2 = (1.0 / 6.0) / (*cell_volume);

    for (int i=0; i<4; i++) { /*loop over vertices*/

        point_lc[i] = coefficient2 * (
            cell_det[i * 4 + 0] -
            cell_det[i * 4 + 1] * point_pos[0] +
            cell_det[i * 4 + 2] * point_pos[1] -
            cell_det[i * 4 + 3] * point_pos[2]);
    }

    if (!(point_lc[0] < 0.0 ||
        point_lc[0] > 1.0 ||
        point_lc[1] < 0.0 ||
        point_lc[1] > 1.0 ||
        point_lc[2] < 0.0 ||
        point_lc[2] > 1.0 ||
        point_lc[3] < 0.0 ||
        point_lc[3] > 1.0)) {

        { opp_move_status_flag = OPP_MOVE_DONE; };
        return;
    }

    // outside the last known cell, find most negative weight and
    // use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = point_lc[0];

    for (int i=1; i<4; i++) {
        if (point_lc[i] < min_lc) {
            min_lc = point_lc[i];
            min_i = i;
        }
    }

    if (opp_c2c[min_i] >= 0) { // is there a neighbor in this direction?
        (*opp_p2c) = opp_c2c[min_i];
        { opp_move_status_flag = OPP_NEED_MOVE; };
    }
    else {
        (*opp_p2c) = INT_MAX;
        { opp_move_status_flag = OPP_NEED_REMOVE; };
    }
}
}

void opp_particle_move__move_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1, // p_lc | OPP_WRITE
    opp_arg arg2, // c_volume | OPP_READ
    opp_arg arg3 // c_det | OPP_READ
) 
{
    if (OPP_DBG) opp_printf("APP", "opp_particle_move__move_kernel set_size %d", set->size);

    opp_profiler->start("move_kernel");

    const int nthreads = omp_get_max_threads();

    const int nargs = 5;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;
    args[4] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

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

            opp_k4::move_kernel(
                move_flag, iter_one_flag, opp_c2c, opp_p2c, 
                (const OPP_REAL *)args[0].data + (n * 3),
                (OPP_REAL *)args[1].data + (n * 4),
                (const OPP_REAL *)args[2].data + (opp_p2c[0] * 1),
                (const OPP_REAL *)args[3].data + (opp_p2c[0] * 16)
            );

        } while (opp_check_part_move_status(move_flag, iter_one_flag, opp_p2c[0], n, thread));
    };

    // ----------------------------------------------------------------------------
    opp_init_particle_move(set, 0, nullptr);

#ifdef USE_MPI 
    if (useGlobalMove) { // For now Global Move will not work with MPI
        
        globalMover->initGlobalMove();

        opp_profiler->start("GblMv_Move");

        // check whether particles needs to be moved over global move routine
        #pragma omp parallel for
        for (int thr = 0; thr < nthreads; thr++)
        {
            const size_t start  = ((size_t)size * thr) / nthreads;
            const size_t finish = ((size_t)size * (thr+1)) / nthreads;

            for (size_t i = start; i < finish; i++)
            {
                OPP_INT* opp_p2c = OPP_mesh_relation_data + i;
                // TODO: we assume pos is in arg 0, Change this!
                const opp_point* point = (const opp_point*)&(((OPP_REAL*)args[0].data)[i * 3]); 

                // check for global move, and if satisfy global move criteria, then remove the particle from current rank
                if (opp_part_checkForGlobalMove3D(set, *point, i, opp_p2c[0])) { 
                    
                    set->particle_remove_count++;
                    continue;  
                }
            }
        }

        opp_profiler->end("GblMv_Move");

        globalMover->communicate(set);
    }
#endif

    opp_profiler->start("Mv_AllMv0");

    // ----------------------------------------------------------------------------
    // check whether all particles not marked for global comm is within cell, 
    // and if not mark to move between cells within the MPI rank, mark for neighbour comm
    opp_profiler->start("move_kernel_only");
    const int total_count = OPP_iter_end - OPP_iter_start;
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
        
        opp_init_particle_move(set, 0, nullptr);

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

// ----------------------------------------------------------------------------
inline void gen_dh_structured_mesh(opp_set set, const opp_dat c_gbl_id, opp_map c2c_map,
        const int nargs, opp_arg* args) { 

    if (OPP_rank == 0)            
        opp_printf("APP", "gen_dh_structured_mesh cells [%s] global grid dims %d %d %d",
            set->name, cellMapper->globalGridDims.x, cellMapper->globalGridDims.y, cellMapper->globalGridDims.z);

    const int set_size_inc_halo = set->size + set->exec_size + set->nonexec_size;
    if (set_size_inc_halo <= 0) {
        opp_printf("APP", "Error... set_size_inc_halo <= 0 for set %s", set->name);
        opp_abort("Error... APP set_size_inc_halo <= 0");
    }


    std::map<size_t, opp_point> removed_coords;
    const opp_point& min_glb_coords = boundingBox->getGlobalMin();
    const opp_point& maxCoordinate = boundingBox->getLocalMax(); // required for GET_VERT define
    bool opp_move_hop_iter_one_flag = true;

    // lambda function for dh mesh search loop
    auto all_cell_checker = [&](const opp_point& point, int& cid) {          
 
        // we dont want to change the original arrays during dh mesh generation, hence duplicate except OPP_READ
        OPP_REAL arg1_temp[4];

        for (int ci = 0; ci < set->size; ++ci) {
            char opp_move_status_flag = OPP_MOVE_DONE;  

            int temp_ci = ci; // we dont want to get iterating ci changed within the kernel, hence get a copy
            
            OPP_INT* opp_p2c = &(temp_ci);           
            OPP_INT* opp_c2c = &((c2c_map->map)[temp_ci * 4]);
    

            opp_k4::move_kernel(
                opp_move_status_flag, opp_move_hop_iter_one_flag, opp_c2c, opp_p2c, 
                (const OPP_REAL*)&point,
                arg1_temp, // p_lc| OPP_WRITE
                (const OPP_REAL *)args[2].data + (temp_ci * 1), // c_volume| OPP_READ
                (const OPP_REAL *)args[3].data + (temp_ci * 16) // c_det| OPP_READ
            );
            if (opp_move_status_flag == OPP_MOVE_DONE) {       
                cid = ci;
                break;
            }
        }
    };

    cellMapper->createStructMeshMappingArrays();

    // Step 1 : Get centroids of the structured mesh cells and try to relate them to unstructured mesh indices
    if (OPP_rank == 0) opp_printf("APP", "gen_dh_structured_mesh Step 1 Start");
    opp_profiler->start("Setup_Mover_s1");
    double x = 0.0, y = 0.0, z = 0.0;
    
    #pragma omp parallel for
    for (int dz = cellMapper->localGridStart.z; dz < cellMapper->localGridEnd.z; dz++) {       
        z = min_glb_coords.z + dz * cellMapper->gridSpacing;        
        for (int dy = cellMapper->localGridStart.y; dy < cellMapper->localGridEnd.y; dy++) {            
            y = min_glb_coords.y + dy * cellMapper->gridSpacing;           
            for (int dx = cellMapper->localGridStart.x; dx < cellMapper->localGridEnd.x; dx++) {                
                x = min_glb_coords.x + dx * cellMapper->gridSpacing;               
                
                size_t index = (size_t)(dx + dy * cellMapper->globalGridDims.x + 
                            dz * cellMapper->globalGridDims.x * cellMapper->globalGridDims.y);                
                const opp_point centroid = cellMapper->getCentroidOfBox(opp_point(x, y ,z));
                int cid = MAX_CELL_INDEX;

                all_cell_checker(centroid, cid); // Find in which cell this centroid lies

                if (cid == MAX_CELL_INDEX) {
                    #pragma omp critical
                    {
                        removed_coords.insert(std::make_pair(index, opp_point(x, y ,z)));
                    }
                }
                else if (cid < set->size) // write only if the structured cell belong to the current MPI rank                    
                    cellMapper->enrichStructuredMesh(index, ((int*)c_gbl_id->data)[cid], OPP_rank);
            }
        }
    }
    opp_profiler->end("Setup_Mover_s1");

    // Step 2 : For MPI, get the inter-node values reduced to the structured mesh
    if (OPP_rank == 0) opp_printf("APP", "gen_dh_structured_mesh Step 2 Start");
    opp_profiler->start("Setup_Mover_s2");
#ifdef USE_MPI
    cellMapper->reduceInterNodeMappings(1);

    // The marked structured cells from this rank might be filled by another rank, so if already filled, 
    // no need to recalculate from current rank
    for (auto it = removed_coords.begin(); it != removed_coords.end(); ) {
        size_t removed_idx = it->first;
        if (cellMapper->structMeshToRankMapping[removed_idx] != MAX_CELL_INDEX) {
            it = removed_coords.erase(it); // This structured index is already written by another rank
            // opp_printf("APP", "index %zu already in %d", this->structMeshToRankMapping[removed_idx], removed_idx);
        } 
        else
            ++it; 
    }

    cellMapper->waitBarrier();    
#endif
    opp_profiler->end("Setup_Mover_s2");

    // Step 3 : Iterate over NEED_REMOVE points, Check whether atleast one vertex of the structured mesh is within 
    //          an unstructured mesh cell. If multiple are found, get the minimum cell index to match with MPI
    if (OPP_rank == 0) opp_printf("APP", "gen_dh_structured_mesh Step 3 Start");
    opp_profiler->start("Setup_Mover_s3");
    std::vector<size_t> removed_coords_keys;
    removed_coords_keys.reserve(removed_coords.size());
    for (const auto& pair : removed_coords)
        removed_coords_keys.push_back(pair.first);

    #pragma omp parallel for
    for (size_t i = 0; i < removed_coords_keys.size(); ++i) {
        const size_t index = removed_coords_keys[i];
        opp_point& p = removed_coords[index];
        double &x = p.x, &y = p.y, &z = p.z;
        
        const double gs = cellMapper->gridSpacing;
        int most_suitable_cid = MAX_CELL_INDEX, most_suitable_gbl_cid = MAX_CELL_INDEX;

        std::array<opp_point,8> vertices = {
            opp_point(GET_VERT(x,x),    GET_VERT(y,y),    GET_VERT(z,z)),
            opp_point(GET_VERT(x,x),    GET_VERT(y,y+gs), GET_VERT(z,z)),
            opp_point(GET_VERT(x,x),    GET_VERT(y,y+gs), GET_VERT(z,z+gs)),
            opp_point(GET_VERT(x,x),    GET_VERT(y,y),    GET_VERT(z,z+gs)),
            opp_point(GET_VERT(x,x+gs), GET_VERT(y,y),    GET_VERT(z,z)),
            opp_point(GET_VERT(x,x+gs), GET_VERT(y,y+gs), GET_VERT(z,z)),
            opp_point(GET_VERT(x,x+gs), GET_VERT(y,y),    GET_VERT(z,z+gs)),
            opp_point(GET_VERT(x,x+gs), GET_VERT(y,y+gs), GET_VERT(z,z+gs)),
        };

        for (const auto& point : vertices) {
            int cid = MAX_CELL_INDEX;

            all_cell_checker(point, cid);

            if ((cid != MAX_CELL_INDEX) && (cid < set_size_inc_halo)) { 
                const int gbl_cid = ((OPP_INT*)c_gbl_id->data)[cid];
                if (most_suitable_gbl_cid > gbl_cid) {
                    most_suitable_gbl_cid = gbl_cid;
                    most_suitable_cid = cid;
                }
            }
        }    

        // Allow neighbours to write on-behalf of the current rank, to reduce issues
        cellMapper->lockWindows();
        int avail_gbl_cid = cellMapper->structMeshToCellMapping[index]; 
        if ((most_suitable_gbl_cid != MAX_CELL_INDEX) && (most_suitable_gbl_cid < avail_gbl_cid) && 
                    (most_suitable_cid < set->size)) {          
            cellMapper->enrichStructuredMesh(index, most_suitable_gbl_cid, OPP_rank);       
        }
        cellMapper->unlockWindows();
    }
    opp_profiler->end("Setup_Mover_s3");

    // Step 4 : For MPI, get the inter-node values reduced to the structured mesh
    if (OPP_rank == 0) opp_printf("APP", "gen_dh_structured_mesh Step 4 Start");
    opp_profiler->start("Setup_Mover_s4");
#ifdef USE_MPI
    cellMapper->reduceInterNodeMappings(2);
#endif
    opp_profiler->end("Setup_Mover_s4");

    // Step 5 : For MPI, convert the global cell coordinates to rank local coordinates for increased performance
    if (OPP_rank == 0) opp_printf("APP", "gen_dh_structured_mesh Step 5 Start");
    opp_profiler->start("Setup_Mover_s5");
#ifdef USE_MPI
    cellMapper->convertToLocalMappings(c_gbl_id);
#endif
    opp_profiler->end("Setup_Mover_s5");

    if (OPP_rank == 0) opp_printf("APP", "gen_dh_structured_mesh DONE");
}

void opp_init_direct_hop_cg(double grid_spacing, int dim, const opp_dat c_gbl_id, const opp::BoundingBox& b_box, 
    opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1, // p_lc | OPP_WRITE
    opp_arg arg2, // c_volume | OPP_READ
    opp_arg arg3 // c_det | OPP_READ
) {
    opp_profiler->start("Setup_Mover");
    
    useGlobalMove = opp_params->get<OPP_BOOL>("opp_global_move");

    if (useGlobalMove) {

        const int nargs = 4;
        opp_arg args[nargs];

        args[0] = arg0;
        args[1] = arg1;
        args[2] = arg2;
        args[3] = arg3;

#ifdef USE_MPI
        opp_mpi_halo_exchanges(c_gbl_id->set, nargs, args);

        comm = std::make_shared<opp::Comm>(MPI_COMM_WORLD);
        globalMover = std::make_unique<opp::GlobalParticleMover>(comm->comm_parent);

        opp_mpi_halo_wait_all(nargs, args);
#endif

        boundingBox = std::make_shared<opp::BoundingBox>(b_box);
        cellMapper = std::make_shared<opp::CellMapper>(boundingBox, grid_spacing, comm);

        gen_dh_structured_mesh(c_gbl_id->set, c_gbl_id, c2c_map, nargs, args);

        opp_profiler->reg("GlbToLocal");
        opp_profiler->reg("GblMv_Move");
        opp_profiler->reg("GblMv_AllMv");
    }

    opp_profiler->end("Setup_Mover");
}


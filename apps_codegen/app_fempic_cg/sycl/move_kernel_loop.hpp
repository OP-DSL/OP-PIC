
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

namespace opp_k4 {

namespace host {
inline void move_kernel(
    const double *point_pos, double* point_lc,
    const double *cell_volume, const double *cell_det
) {

    // Step (1) : Specifying computations to be carried out for each mesh element, until its final destination cell
    const double coefficient2 = (1.0 / 6.0) / (*cell_volume);
    for (int i=0; i<4; i++) { /*loop over vertices*/

        point_lc[i] = coefficient2 * (
            cell_det[i * 4 + 0] -
            cell_det[i * 4 + 1] * point_pos[0] +
            cell_det[i * 4 + 2] * point_pos[1] -
            cell_det[i * 4 + 3] * point_pos[2]);
    }
    // Step (1) : end

    // Step (2) : A method to identify if the particle has reached its final mesh cell
    if ((point_lc[0] > (-1) * 1e-10) && ((point_lc[0] - 1.0) < 1e-10) &&
        (point_lc[1] > (-1) * 1e-10) && ((point_lc[1] - 1.0) < 1e-10) &&
        (point_lc[2] > (-1) * 1e-10) && ((point_lc[2] - 1.0) < 1e-10) &&
        (point_lc[3] > (-1) * 1e-10) && ((point_lc[3] - 1.0) < 1e-10)) {

        // Step (3) : Computations to be carried out at the final destination mesh cell
        /* No computations in Mini-FEM-PIC */

        { opp_move_status_flag = OPP_MOVE_DONE; };
        return;
    }

    // Outside the last known cell, find most negative weight and move to it if valid

    // Step (5) : Calculate the next most probable cell index to search
    int min_i = 0;
    double min_lc = point_lc[0];
    for (int i=1; i<4; i++) {
        if (point_lc[i] < min_lc) {
            min_lc = point_lc[i];
            min_i = i;
        }
    }

    if (opp_c2c[min_i] >= 0) { // is there a neighbor in this direction?

        // Step X : Any calculations required on all hopped cells
        /* No computations in Mini-FEM-PIC */

        (*opp_p2c) = opp_c2c[min_i]; // Assign new p2c mapping
        { opp_move_status_flag = OPP_NEED_MOVE; };
    }
    else {
        // Step (4) : Actions to be carried out if the particle has moved out of the mesh domain
        (*opp_p2c) = INT_MAX;
        { opp_move_status_flag = OPP_NEED_REMOVE; };
    }
}
}

void opp_dev_move_kernel_sycl(opp_set set, const int nargs, opp_arg *args, opp_map p2c_map, opp_map c2c_map) {

    opp_set_stride(OPP_comm_iteration_d, OPP_comm_iteration_h, OPP_comm_iteration);
    opp_set_stride(opp_k4_dat0_stride_s, opp_k4_dat0_stride, args[0].dat->set->set_capacity);
    opp_set_stride(opp_k4_dat1_stride_s, opp_k4_dat1_stride, args[1].dat->set->set_capacity);
    opp_set_stride(opp_k4_dat2_stride_s, opp_k4_dat2_stride, args[2].dat->set->set_capacity);
    opp_set_stride(opp_k4_dat3_stride_s, opp_k4_dat3_stride, args[3].dat->set->set_capacity);

#ifdef OPP_BLOCK_SIZE_4
    const int block_size = OPP_BLOCK_SIZE_4;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif
    
    const int num_blocks = (OPP_iter_end - OPP_iter_start - 1) / block_size + 1;

    opp_profiler->start("move_kernel_only");

    opp_queue->submit([&](sycl::handler &cgh) {
        
        const OPP_INT* comm_iteration = OPP_comm_iteration_d;
        const OPP_INT* cell_set_size = OPP_cells_set_size_d;

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

            // Step (1) : Specifying computations to be carried out for each mesh element, until its final destination cell
            const double coefficient2 = (1.0 / 6.0) / (*cell_volume);
            for (int i=0; i<4; i++) { /*loop over vertices*/

                point_lc[(i) * opp_k4_dat1_stride_sycl[0]] = coefficient2 * (
                    cell_det[(i * 4 + 0) * opp_k4_dat3_stride_sycl[0]] -
                    cell_det[(i * 4 + 1) * opp_k4_dat3_stride_sycl[0]] * point_pos[(0) * opp_k4_dat0_stride_sycl[0]] +
                    cell_det[(i * 4 + 2) * opp_k4_dat3_stride_sycl[0]] * point_pos[(1) * opp_k4_dat0_stride_sycl[0]] -
                    cell_det[(i * 4 + 3) * opp_k4_dat3_stride_sycl[0]] * point_pos[(2) * opp_k4_dat0_stride_sycl[0]]);
            }
            // Step (1) : end

            // Step (2) : A method to identify if the particle has reached its final mesh cell
            if ((point_lc[(0) * opp_k4_dat1_stride_sycl[0]] > (-1) * 1e-10) && ((point_lc[(0) * opp_k4_dat1_stride_sycl[0]] - 1.0) < 1e-10) &&
                (point_lc[(1) * opp_k4_dat1_stride_sycl[0]] > (-1) * 1e-10) && ((point_lc[(1) * opp_k4_dat1_stride_sycl[0]] - 1.0) < 1e-10) &&
                (point_lc[(2) * opp_k4_dat1_stride_sycl[0]] > (-1) * 1e-10) && ((point_lc[(2) * opp_k4_dat1_stride_sycl[0]] - 1.0) < 1e-10) &&
                (point_lc[(3) * opp_k4_dat1_stride_sycl[0]] > (-1) * 1e-10) && ((point_lc[(3) * opp_k4_dat1_stride_sycl[0]] - 1.0) < 1e-10)) {

                // Step (3) : Computations to be carried out at the final destination mesh cell
                /* No computations in Mini-FEM-PIC */

                { opp_move_status_flag = OPP_MOVE_DONE; };
                return;
            }

            // Outside the last known cell, find most negative weight and move to it if valid

            // Step (5) : Calculate the next most probable cell index to search
            int min_i = 0;
            double min_lc = point_lc[(0) * opp_k4_dat1_stride_sycl[0]];
            for (int i=1; i<4; i++) {
                if (point_lc[(i) * opp_k4_dat1_stride_sycl[0]] < min_lc) {
                    min_lc = point_lc[(i) * opp_k4_dat1_stride_sycl[0]];
                    min_i = i;
                }
            }

            if (opp_c2c[(min_i) * opp_k4_c2c_map_stride_sycl[0]] >= 0) { // is there a neighbor in this direction?

                // Step X : Any calculations required on all hopped cells
                /* No computations in Mini-FEM-PIC */

                (*opp_p2c) = opp_c2c[(min_i) * opp_k4_c2c_map_stride_sycl[0]]; // Assign new p2c mapping
                { opp_move_status_flag = OPP_NEED_MOVE; };
            }
            else {
                // Step (4) : Actions to be carried out if the particle has moved out of the mesh domain
                (*opp_p2c) = INT_MAX;
                { opp_move_status_flag = OPP_NEED_REMOVE; };
            }
        };

        // -----------------------------------------------------------------------------------------
        auto opp_move_kernel = [=](sycl::nd_item<1> item) {
            
            const int tid = item.get_global_linear_id();
            const int n = tid + iter_start;

            if (n < iter_end) {
                OPP_INT *opp_p2c = (p2c_map_sycl + n);
                if (opp_p2c[0] == MAX_CELL_INDEX) {
                    return;
                }
                
                char move_flag = OPP_NEED_MOVE;
                bool iter_one_flag = (comm_iteration[0] > 0) ? false : true;

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
            
                } while (opp_part_check_status_device(move_flag, iter_one_flag, opp_p2c, n, 
                            remove_count, remove_part_indices, 
                            move_part_indices, move_cell_indices, move_count, cell_set_size));
            }        
        };

        // -----------------------------------------------------------------------------------------
        cgh.parallel_for<class opp_particle_move>(
                sycl::nd_range<1>(block_size * num_blocks, block_size), opp_move_kernel);
    });

    OPP_DEVICE_SYNCHRONIZE();

    opp_profiler->end("move_kernel_only");
}
} // end of namespace

//--------------------------------------------------------------
void opp_particle_move__move_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0,   // p_pos | OPP_READ
    opp_arg arg1,   // p_lc | OPP_WRITE
    opp_arg arg2,   // c_volume | OPP_READ
    opp_arg arg3   // c_det | OPP_READ
) 
{ OPP_RETURN_IF_INVALID_PROCESS;

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
 
    const OPP_INT c2c_stride = c2c_map->from->size + c2c_map->from->exec_size + c2c_map->from->nonexec_size;

    opp_set_stride(OPP_cells_set_size_d, OPP_cells_set_size, set->cells_set->size);
    opp_set_stride(opp_k4_c2c_map_stride_s, opp_k4_c2c_map_stride, c2c_stride);

    opp_mpi_halo_wait_all(nargs, args);

    opp_init_particle_move(set, nargs, args);

    if (useGlobalMove) {
           
#ifdef USE_MPI         
        globalMover->initGlobalMove();
        opp_init_dh_device(set);
#endif
        opp_profiler->start("GblMv_Move");
        // check whether particles need to be moved via the global move routine
        opp_dev_checkForGlobalMove_sycl(set,
            (OPP_REAL*)args[0].data_d,    // p_pos 
            (OPP_INT *)args[4].data_d     // p2c_map
        );
        opp_profiler->end("GblMv_Move");

#ifdef USE_MPI 
        opp_gather_dh_move_indices(set);
        globalMover->communicate(set);
#endif
    }

    opp_k4::opp_dev_move_kernel_sycl(set, nargs, args, p2c_map, c2c_map);

#ifdef USE_MPI 
    // ----------------------------------------------------------------------------
    // finalize the global move routine and iterate over newly added particles and check whether they need neighbour comm
    if (useGlobalMove) { 
        
        opp_profiler->start("GblMv_finalize");
        const int finalized = globalMover->finalize(set);
        opp_profiler->end("GblMv_finalize");

        if (finalized > 0) {
            opp_profiler->start("GblMv_AllMv");

            // need to change arg data since particle resize in globalMover::finalize could change the pointer in dat->data 
            for (int i = 0; i < nargs; i++)
                if (args[i].argtype == OPP_ARG_DAT && args[i].dat->set->is_particle)
                    args[i].data_d = args[i].dat->data_d;

            // check whether the new particle is within cell, and if not move between cells within the MPI rank, 
            // mark for neighbour comm. Do only for the globally moved particles 
            OPP_iter_start = (set->size - set->diff);
            OPP_iter_end = set->size;
            
            opp_k4::opp_dev_move_kernel_sycl(set, nargs, args, p2c_map, c2c_map); 

            opp_profiler->end("GblMv_AllMv");
        }
    }
#endif

    while (opp_finalize_particle_move(set)) {
        
        opp_init_particle_move(set, nargs, args);

        opp_k4::opp_dev_move_kernel_sycl(set, nargs, args, p2c_map, c2c_map);
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();
 
    opp_profiler->end("move_kernel");
}

void opp_init_direct_hop_cg(double grid_spacing, const opp_dat c_gbl_id, const opp::BoundingBox& b_box, 
    opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1, // p_lc | OPP_WRITE
    opp_arg arg2, // c_volume | OPP_READ
    opp_arg arg3 // c_det | OPP_READ
) { OPP_RETURN_IF_INVALID_PROCESS;

    opp_profiler->start("Setup_Mover");

    useGlobalMove = opp_params->get<OPP_BOOL>("opp_global_move");

    if (OPP_DBG) opp_printf("opp_init_direct_hop_cg", "START useGlobalMove=%s", useGlobalMove ? "YES" : "NO");

    if (useGlobalMove) {

        const int nargs = 5;
        opp_arg args[nargs];

        args[0] = arg0;
        args[1] = arg1;
        args[2] = arg2;
        args[3] = arg3;
        args[4] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

#ifdef USE_MPI
        opp_mpi_halo_exchanges_grouped(c_gbl_id->set, nargs, args, Device_CPU);

        comm = std::make_shared<opp::Comm>(OPP_MPI_WORLD);
        globalMover = std::make_unique<opp::GlobalParticleMover>(comm->comm_parent);

        opp_mpi_halo_wait_all(nargs, args);
#endif

        boundingBox = std::make_shared<opp::BoundingBox>(b_box);
        cellMapper = std::make_shared<opp::CellMapper>(boundingBox, grid_spacing, comm);
        
        const int c_set_size = c_gbl_id->set->size;

        // lambda function for dh mesh search loop
        auto all_cell_checker = [&](const opp_point& point, int& cid) {          
 
            // we dont want to change the original arrays during dh mesh generation, hence duplicate except OPP_READ
            OPP_REAL arg1_temp[4];
            for (int ci = 0; ci < c_set_size; ++ci) {
                opp_move_status_flag = OPP_NEED_MOVE;  
                opp_move_hop_iter_one_flag = true;
                
                int temp_ci = ci; // we dont want to get iterating ci changed within the kernel, hence get a copy
                
                opp_p2c = &(temp_ci);           
                opp_c2c = &((c2c_map->map)[temp_ci * 4]);

                opp_k4::host::move_kernel(
                    (const OPP_REAL*)&point,
                    arg1_temp, // p_lc| OPP_WRITE
                    (const OPP_REAL *)args[2].data + (temp_ci * 1), // c_volume| OPP_READ
                    (const OPP_REAL *)args[3].data + (temp_ci * 16) // c_det| OPP_READ
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
    }

    opp_profiler->end("Setup_Mover");
}

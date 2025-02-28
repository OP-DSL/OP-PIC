
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k4_dat0_stride = -1;
OPP_INT opp_k4_dat1_stride = -1;
OPP_INT opp_k4_dat2_stride = -1;
OPP_INT opp_k4_dat3_stride = -1;
OPP_INT opp_k4_c2c_map_stride = -1;

__constant__ OPP_INT opp_k4_dat0_stride_d;
__constant__ OPP_INT opp_k4_dat1_stride_d;
__constant__ OPP_INT opp_k4_dat2_stride_d;
__constant__ OPP_INT opp_k4_dat3_stride_d;
__constant__ OPP_INT opp_k4_c2c_map_stride_d;



namespace opp_k4 {

namespace host {
inline void move_kernel(
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

    if ((point_lc[0] > -1e-10) && ((point_lc[0] - 1.0) < 1e-10) &&
        (point_lc[1] > -1e-10) && ((point_lc[1] - 1.0) < 1e-10) &&
        (point_lc[2] > -1e-10) && ((point_lc[2] - 1.0) < 1e-10) &&
        (point_lc[3] > -1e-10) && ((point_lc[3] - 1.0) < 1e-10)) {

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

__device__ inline void move_kernel(
    char& opp_move_status_flag, const bool opp_move_hop_iter_one_flag, // Added by code-gen
    const OPP_INT* opp_c2c, OPP_INT* opp_p2c, // Added by code-gen
    const double *point_pos, double* point_lc,
    const double *cell_volume, const double *cell_det
) {
    const double coefficient2 = (1.0 / 6.0) / (*cell_volume);

    for (int i=0; i<4; i++) { /*loop over vertices*/

        point_lc[(i) * opp_k4_dat1_stride_d] = coefficient2 * (
            cell_det[(i * 4 + 0) * opp_k4_dat3_stride_d] -
            cell_det[(i * 4 + 1) * opp_k4_dat3_stride_d] * point_pos[(0) * opp_k4_dat0_stride_d] +
            cell_det[(i * 4 + 2) * opp_k4_dat3_stride_d] * point_pos[(1) * opp_k4_dat0_stride_d] -
            cell_det[(i * 4 + 3) * opp_k4_dat3_stride_d] * point_pos[(2) * opp_k4_dat0_stride_d]);
    }

    if ((point_lc[(0) * opp_k4_dat1_stride_d] > -1e-18) && ((point_lc[(0) * opp_k4_dat1_stride_d] - 1.0) < -1e-18) &&
        (point_lc[(1) * opp_k4_dat1_stride_d] > -1e-18) && ((point_lc[(1) * opp_k4_dat1_stride_d] - 1.0) < -1e-18) &&
        (point_lc[(2) * opp_k4_dat1_stride_d] > -1e-18) && ((point_lc[(2) * opp_k4_dat1_stride_d] - 1.0) < -1e-18) &&
        (point_lc[(3) * opp_k4_dat1_stride_d] > -1e-18) && ((point_lc[(3) * opp_k4_dat1_stride_d] - 1.0) < -1e-18)) {

        { opp_move_status_flag = OPP_MOVE_DONE; };
        return;
    }

    // outside the last known cell, find most negative weight and
    // use that cell_index to reduce computations
    int min_i = 0;
    double min_lc = point_lc[(0) * opp_k4_dat1_stride_d];

    for (int i=1; i<4; i++) {
        if (point_lc[(i) * opp_k4_dat1_stride_d] < min_lc) {
            min_lc = point_lc[(i) * opp_k4_dat1_stride_d];
            min_i = i;
        }
    }

    if (opp_c2c[(min_i) * opp_k4_c2c_map_stride_d] >= 0) { // is there a neighbor in this direction?
        (*opp_p2c) = opp_c2c[(min_i) * opp_k4_c2c_map_stride_d];
        { opp_move_status_flag = OPP_NEED_MOVE; };
    }
    else {
        (*opp_p2c) = INT_MAX;
        { opp_move_status_flag = OPP_NEED_REMOVE; };
    }
}
}

//--------------------------------------------------------------
__global__ void opp_dev_move_kernel(
    const OPP_REAL *__restrict__ dat0,     // p_pos
    OPP_REAL *__restrict__ dat1,     // p_lc
    const OPP_REAL *__restrict__ dat2,     // c_volume
    const OPP_REAL *__restrict__ dat3,     // c_det
    OPP_INT *__restrict__ p2c_map,
    const OPP_INT *__restrict__ c2c_map,
    OPP_INT *__restrict__ particle_remove_count,
    OPP_INT *__restrict__ particle_remove_indices,
    OPP_INT *__restrict__ move_particle_indices,
    OPP_INT *__restrict__ move_cell_indices,
    OPP_INT *__restrict__ move_count,
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int n = OPP_DEVICE_GLOBAL_LINEAR_ID + start;

    if (n < end) {

        OPP_INT *opp_p2c = (p2c_map + n);
        if (opp_p2c[0] == MAX_CELL_INDEX) {
            return;
        }

        char move_flag = OPP_NEED_MOVE;
        bool iter_one_flag = (OPP_comm_iteration_d > 0) ? false : true;

        do {
            const OPP_INT p2c = opp_p2c[0]; // get the value here, since the kernel might change it
            const OPP_INT* opp_c2c = c2c_map + p2c;           

            opp_k4::move_kernel(
                move_flag, iter_one_flag, opp_c2c, opp_p2c,
                dat0 + n, // p_pos 
                dat1 + n, // p_lc 
                dat2 + p2c, // c_volume 
                dat3 + p2c // c_det 
          
            );

        } while (opp_part_check_status_device(move_flag, iter_one_flag, opp_p2c, n, 
            *particle_remove_count, particle_remove_indices, move_particle_indices, 
            move_cell_indices, move_count));        
    }
    
}


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

    opp_mem::dev_copy_to_symbol<OPP_INT>(OPP_cells_set_size_d, &OPP_cells_set_size, &(set->cells_set->size), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_c2c_map_stride_d, &opp_k4_c2c_map_stride, &c2c_stride, 1);

    opp_mpi_halo_wait_all(nargs, args);

#ifdef OPP_BLOCK_SIZE_4
    const int block_size = OPP_BLOCK_SIZE_4;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    int num_blocks = 200;

    opp_init_particle_move(set, nargs, args);
    if (useGlobalMove) {
           
#ifdef USE_MPI         
        globalMover->initGlobalMove();
        opp_init_dh_device(set);
#endif
        opp_profiler->start("GblMv_Move");

        opp_mem::dev_copy_to_symbol<OPP_INT>(cellMapper_pos_stride_d, &cellMapper_pos_stride, &(args[0].dat->set->set_capacity), 1);
        opp_mem::dev_copy_to_symbol<OPP_INT>(OPP_rank_d, &OPP_rank, 1);

        OPP_DEV_CHECK(hipMemcpyToSymbol(HIP_SYMBOL(opp_minSavedDHGrid_d), opp_minSavedDHGrid, 3 * sizeof(size_t)));
        OPP_DEV_CHECK(hipMemcpyToSymbol(HIP_SYMBOL(opp_maxSavedDHGrid_d), opp_maxSavedDHGrid, 3 * sizeof(size_t)));

        // check whether particles need to be moved via the global move routine
        num_blocks = (OPP_iter_end - OPP_iter_start - 1) / block_size + 1;
        opp_dev_checkForGlobalMove3D_kernel<<<num_blocks, block_size>>>(
            (OPP_REAL*)args[0].data_d,    // p_pos 
            (OPP_INT *)args[4].data_d,    // p2c_map
            cellMapper->structMeshToCellMapping_d, 
            cellMapper->structMeshToRankMapping_d,
            cellMapper->oneOverGridSpacing_d, 
            cellMapper->minGlbCoordinate_d, 
            cellMapper->globalGridDims_d, 
            cellMapper->globalGridSize_d,
            set->particle_remove_count_d, 
            OPP_remove_particle_indices_d, 
            dh_indices_d.part_indices, 
            dh_indices_d.cell_indices, 
            dh_indices_d.rank_indices, 
            dh_indices_d.move_count,
            OPP_iter_start, OPP_iter_end
        );
        OPP_DEVICE_SYNCHRONIZE();

        opp_profiler->end("GblMv_Move");

#ifdef USE_MPI 
        opp_gather_dh_move_indices(set);
        globalMover->communicate(set);
#endif
    }




    opp_profiler->start("Mv_AllMv0");
    // ----------------------------------------------------------------------------
    // check whether all particles not marked for global comm is within cell, 
    // and if not mark to move between cells within the MPI rank, mark for neighbour comm

    opp_mem::dev_copy_to_symbol<OPP_INT>(OPP_comm_iteration_d, &OPP_comm_iteration, 1);
    num_blocks = (OPP_iter_end - OPP_iter_start - 1) / block_size + 1;
    
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat0_stride_d, &opp_k4_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat1_stride_d, &opp_k4_dat1_stride, &(args[1].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat2_stride_d, &opp_k4_dat2_stride, &(args[2].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat3_stride_d, &opp_k4_dat3_stride, &(args[3].dat->set->set_capacity), 1);

    {
        opp_profiler->start("move_kernel_only");
        opp_dev_move_kernel<<<num_blocks, block_size>>>(
            (OPP_REAL *)args[0].data_d,    // p_pos
            (OPP_REAL *)args[1].data_d,    // p_lc
            (OPP_REAL *)args[2].data_d,    // c_volume
            (OPP_REAL *)args[3].data_d,    // c_det
            (OPP_INT *)args[4].data_d,    // p2c_map
            (OPP_INT *)c2c_map->map_d,    // c2c_map
            (OPP_INT *)set->particle_remove_count_d,
            (OPP_INT *)OPP_remove_particle_indices_d,
            (OPP_INT *)OPP_move_particle_indices_d,
            (OPP_INT *)OPP_move_cell_indices_d,
            (OPP_INT *)OPP_move_count_d,
            OPP_iter_start,
            OPP_iter_end
        );
        OPP_DEVICE_SYNCHRONIZE(); 
        opp_profiler->end("move_kernel_only");
    }


    opp_profiler->end("Mv_AllMv0");
    
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
            const int start2 = (set->size - set->diff);
            const int end2 = set->size;
            num_blocks = (end2 - start2 - 1) / block_size + 1;

            opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat0_stride_d, &opp_k4_dat0_stride, &(args[0].dat->set->set_capacity), 1);
            opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat1_stride_d, &opp_k4_dat1_stride, &(args[1].dat->set->set_capacity), 1);
            opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat2_stride_d, &opp_k4_dat2_stride, &(args[2].dat->set->set_capacity), 1);
            opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat3_stride_d, &opp_k4_dat3_stride, &(args[3].dat->set->set_capacity), 1);

            opp_profiler->start("move_kernel_only");      
            opp_dev_move_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,    // p_pos
                (OPP_REAL *)args[1].data_d,    // p_lc
                (OPP_REAL *)args[2].data_d,    // c_volume
                (OPP_REAL *)args[3].data_d,    // c_det
                (OPP_INT *)args[4].data_d,    // p2c_map
                (OPP_INT *)c2c_map->map_d,    // c2c_map
                (OPP_INT *)set->particle_remove_count_d,
                (OPP_INT *)OPP_remove_particle_indices_d,
                (OPP_INT *)OPP_move_particle_indices_d,
                (OPP_INT *)OPP_move_cell_indices_d,
                (OPP_INT *)OPP_move_count_d,
                start2,
                end2
            );
            OPP_DEVICE_SYNCHRONIZE(); 
            opp_profiler->end("move_kernel_only");

            opp_profiler->end("GblMv_AllMv");
        }
    }
#endif

    // ----------------------------------------------------------------------------
    // Do neighbour communication and if atleast one particle is received by the currect rank, 
    // then iterate over the newly added particles
    while (opp_finalize_particle_move(set)) {

        opp_init_particle_move(set, nargs, args);

        opp_mem::dev_copy_to_symbol<OPP_INT>(OPP_comm_iteration_d, &OPP_comm_iteration, 1);
        num_blocks = (OPP_iter_end - OPP_iter_start - 1) / block_size + 1;
        
        opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat0_stride_d, &opp_k4_dat0_stride, &(args[0].dat->set->set_capacity), 1);
        opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat1_stride_d, &opp_k4_dat1_stride, &(args[1].dat->set->set_capacity), 1);
        opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat2_stride_d, &opp_k4_dat2_stride, &(args[2].dat->set->set_capacity), 1);
        opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k4_dat3_stride_d, &opp_k4_dat3_stride, &(args[3].dat->set->set_capacity), 1);

        {
            opp_profiler->start("move_kernel_only");
            opp_dev_move_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,    // p_pos
                (OPP_REAL *)args[1].data_d,    // p_lc
                (OPP_REAL *)args[2].data_d,    // c_volume
                (OPP_REAL *)args[3].data_d,    // c_det
                (OPP_INT *)args[4].data_d,    // p2c_map
                (OPP_INT *)c2c_map->map_d,    // c2c_map
                (OPP_INT *)set->particle_remove_count_d,
                (OPP_INT *)OPP_remove_particle_indices_d,
                (OPP_INT *)OPP_move_particle_indices_d,
                (OPP_INT *)OPP_move_cell_indices_d,
                (OPP_INT *)OPP_move_count_d,
                OPP_iter_start,
                OPP_iter_end
            );
            OPP_DEVICE_SYNCHRONIZE(); 
            opp_profiler->end("move_kernel_only");
        }

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


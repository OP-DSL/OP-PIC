
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k2_dat0_stride = -1;
OPP_INT opp_k2_dat1_stride = -1;
OPP_INT opp_k2_c2c_map_stride = -1;

__constant__ OPP_INT opp_k2_dat0_stride_d;
__constant__ OPP_INT opp_k2_dat1_stride_d;
__constant__ OPP_INT opp_k2_c2c_map_stride_d;



namespace opp_k2 {

namespace host {
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

__device__ inline void move_kernel(char& opp_move_status_flag, const bool opp_move_hop_iter_one_flag, // Added by code-gen
    const OPP_INT* opp_c2c, OPP_INT* opp_p2c, // Added by code-gen
    const double* p_pos, const double* c_pos_ll)
{
    // check for x direction movement
    const double p_pos_x = p_pos[(Dim::x) * opp_k2_dat0_stride_d];
    if (p_pos_x < c_pos_ll[(Dim::x) * opp_k2_dat1_stride_d]) {
        opp_p2c[0] = opp_c2c[(CellMap::xd_y) * opp_k2_c2c_map_stride_d];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }
    if (p_pos_x > (c_pos_ll[(Dim::x) * opp_k2_dat1_stride_d] + CONST_cell_width_d[0])) {
        opp_p2c[0] = opp_c2c[(CellMap::xu_y) * opp_k2_c2c_map_stride_d];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }

    // check for y direction movement
    const double p_pos_y = p_pos[(Dim::y) * opp_k2_dat0_stride_d];
    if (p_pos_y < c_pos_ll[(Dim::y) * opp_k2_dat1_stride_d]) {
        opp_p2c[0] = opp_c2c[(CellMap::x_yd) * opp_k2_c2c_map_stride_d];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }
    if (p_pos_y > (c_pos_ll[(Dim::y) * opp_k2_dat1_stride_d] + CONST_cell_width_d[0])) {
        opp_p2c[0] = opp_c2c[(CellMap::x_yu) * opp_k2_c2c_map_stride_d];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }

    { opp_move_status_flag = OPP_MOVE_DONE; };
}
}

//--------------------------------------------------------------
__global__ void opp_dev_move_kernel(
    const OPP_REAL *__restrict__ dat0,     // p_pos
    const OPP_REAL *__restrict__ dat1,     // c_pos_ll
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

            opp_k2::move_kernel(
                move_flag, iter_one_flag, opp_c2c, opp_p2c,
                dat0 + n, // p_pos 
                dat1 + p2c // c_pos_ll 
          
            );

        } while (opp_part_check_status_device(move_flag, iter_one_flag, opp_p2c, n, 
            *particle_remove_count, particle_remove_indices, move_particle_indices, 
            move_cell_indices, move_count));        
    }
    
}


void opp_particle_move__move_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0,   // p_pos | OPP_READ
    opp_arg arg1   // c_pos_ll | OPP_READ
) 
{
    if (OPP_DBG) opp_printf("APP", "opp_particle_move__move_kernel set_size %d", set->size);

    opp_profiler->start("move_kernel");

    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

    opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);

    const OPP_INT c2c_stride = c2c_map->from->size + c2c_map->from->exec_size + c2c_map->from->nonexec_size;

    opp_mem::dev_copy_to_symbol<OPP_INT>(OPP_cells_set_size_d, &OPP_cells_set_size, &(set->cells_set->size), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_c2c_map_stride_d, &opp_k2_c2c_map_stride, &c2c_stride, 1);

    opp_mpi_halo_wait_all(nargs, args);

#ifdef OPP_BLOCK_SIZE_2
    const int block_size = OPP_BLOCK_SIZE_2;
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

        hipError_t err1 = hipMemcpyToSymbol(HIP_SYMBOL(opp_minSavedDHGrid_d), opp_minSavedDHGrid, 3 * sizeof(size_t));
        if (err1 != hipSuccess) {
            throw std::runtime_error(std::string("dev_copy_to_symbol: ") + hipGetErrorString(err1));
        }
        hipError_t err2 = hipMemcpyToSymbol(HIP_SYMBOL(opp_maxSavedDHGrid_d), opp_maxSavedDHGrid, 3 * sizeof(size_t));
        if (err2 != hipSuccess) {
            throw std::runtime_error(std::string("dev_copy_to_symbol: ") + hipGetErrorString(err2));
        }

        // check whether particles need to be moved via the global move routine
        num_blocks = (OPP_iter_end - OPP_iter_start - 1) / block_size + 1;
        opp_dev_checkForGlobalMove2D_kernel<<<num_blocks, block_size>>>(
            (OPP_REAL*)args[0].data_d,    // p_pos 
            (OPP_INT *)args[2].data_d,    // p2c_map
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
    
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat0_stride_d, &opp_k2_dat0_stride, &(args[0].dat->set->set_capacity), 1);
    opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat1_stride_d, &opp_k2_dat1_stride, &(args[1].dat->set->set_capacity), 1);

    {
        opp_profiler->start("move_kernel_only");
        opp_dev_move_kernel<<<num_blocks, block_size>>>(
            (OPP_REAL *)args[0].data_d,    // p_pos
            (OPP_REAL *)args[1].data_d,    // c_pos_ll
            (OPP_INT *)args[2].data_d,    // p2c_map
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

            opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat0_stride_d, &opp_k2_dat0_stride, &(args[0].dat->set->set_capacity), 1);
            opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat1_stride_d, &opp_k2_dat1_stride, &(args[1].dat->set->set_capacity), 1);

            opp_profiler->start("move_kernel_only");      
            opp_dev_move_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,    // p_pos
                (OPP_REAL *)args[1].data_d,    // c_pos_ll
                (OPP_INT *)args[2].data_d,    // p2c_map
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
        
        opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat0_stride_d, &opp_k2_dat0_stride, &(args[0].dat->set->set_capacity), 1);
        opp_mem::dev_copy_to_symbol<OPP_INT>(opp_k2_dat1_stride_d, &opp_k2_dat1_stride, &(args[1].dat->set->set_capacity), 1);

        {
            opp_profiler->start("move_kernel_only");
            opp_dev_move_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,    // p_pos
                (OPP_REAL *)args[1].data_d,    // c_pos_ll
                (OPP_INT *)args[2].data_d,    // p2c_map
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
    opp_arg arg1 // c_pos_ll | OPP_READ
) {
    opp_profiler->start("Setup_Mover");

    useGlobalMove = opp_params->get<OPP_BOOL>("opp_global_move");

    if (OPP_DBG) opp_printf("opp_init_direct_hop_cg", "START useGlobalMove=%s", useGlobalMove ? "YES" : "NO");

    if (useGlobalMove) {

        const int nargs = 3;
        opp_arg args[nargs];

        args[0] = arg0;
        args[1] = arg1;
        args[2] = opp_arg_dat(p2c_map->p2c_dat, OPP_RW); // required to make dirty or should manually make it dirty

#ifdef USE_MPI
        opp_mpi_halo_exchanges_grouped(c_gbl_id->set, nargs, args, Device_CPU);

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

                opp_k2::host::move_kernel(
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
    }

    opp_profiler->end("Setup_Mover");
}


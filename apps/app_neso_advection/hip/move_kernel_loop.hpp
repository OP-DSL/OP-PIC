
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k2_dat0_stride = -1;
OPP_INT opp_k2_dat1_stride = -1;
OPP_INT opp_k2_dat2_stride = -1;
OPP_INT opp_k2_c2c_map_stride = -1;

__constant__ OPP_INT opp_k2_dat0_stride_d;
__constant__ OPP_INT opp_k2_dat1_stride_d;
__constant__ OPP_INT opp_k2_dat2_stride_d;
__constant__ OPP_INT opp_k2_c2c_map_stride_d;



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

__device__ inline void move_kernel(char& opp_move_status_flag, const bool opp_move_hop_iter_one_flag, // Added by code-gen
    const OPP_INT* opp_c2c, OPP_INT* opp_p2c, // Added by code-gen
    const double* p_pos, int* p_mdir, const double* c_pos_ll)
{
    // check for x direction movement
    const double p_pos_x_diff = (p_pos[(Dim::x) * opp_k2_dat0_stride_d] - c_pos_ll[(Dim::x) * opp_k2_dat2_stride_d]);
    if ((p_pos_x_diff >= 0.0) && (p_pos_x_diff <= CONST_cell_width_d[0])) {
        p_mdir[(Dim::x) * opp_k2_dat1_stride_d] = 0; // within cell in x direction
    }
    else if (p_mdir[(Dim::x) * opp_k2_dat1_stride_d] > 0) {
        opp_p2c[0] = opp_c2c[(CellMap::xu_y) * opp_k2_c2c_map_stride_d];

        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }
    else if (p_mdir[(Dim::x) * opp_k2_dat1_stride_d] < 0) {
        opp_p2c[0] = opp_c2c[(CellMap::xd_y) * opp_k2_c2c_map_stride_d];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }

    // check for y direction movement
    const double p_pos_y_diff = (p_pos[(Dim::y) * opp_k2_dat0_stride_d] - c_pos_ll[(Dim::y) * opp_k2_dat2_stride_d]);
    if ((p_pos_y_diff >= 0.0) && (p_pos_y_diff <= CONST_cell_width_d[0])) {
        p_mdir[(Dim::y) * opp_k2_dat1_stride_d] = 0; // within cell in y direction
    }
    else if (p_mdir[(Dim::y) * opp_k2_dat1_stride_d] > 0) {
        opp_p2c[0] = opp_c2c[(CellMap::x_yu) * opp_k2_c2c_map_stride_d];

        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }
    else if (p_mdir[(Dim::y) * opp_k2_dat1_stride_d] < 0) {
        opp_p2c[0] = opp_c2c[(CellMap::x_yd) * opp_k2_c2c_map_stride_d];
        { opp_move_status_flag = OPP_NEED_MOVE; }; return;
    }

    { opp_move_status_flag = OPP_MOVE_DONE; };
}

//*******************************************************************************
// Returns true only if another hop is required by the current rank
__device__ inline bool opp_part_check_status_hip(char& move_flag, bool& iter_one_flag, 
        int* cell_id, int particle_index, int& remove_count, int *remove_particle_indices, 
        int *move_particle_indices, int *move_cell_indices, int *move_count) 
{
    iter_one_flag = false;

    if (move_flag == OPP_MOVE_DONE)
    {
        return false;
    }
    else if (move_flag == OPP_NEED_REMOVE)
    {
        *cell_id = MAX_CELL_INDEX;
        const int removeArrayIndex = atomicAdd(&remove_count, 1);
        remove_particle_indices[removeArrayIndex] = particle_index;

        return false;
    }
    else if (*cell_id >= OPP_cells_set_size_d)
    {
        // cell_id cell is not owned by the current mpi rank, need to communicate
        int moveArrayIndex = atomicAdd(move_count, 1);
        move_particle_indices[moveArrayIndex] = particle_index;
        move_cell_indices[moveArrayIndex] = *cell_id;

        // Needs to be removed from the current rank, 
        // particle packing will be done just prior exchange and removal
        move_flag = OPP_NEED_REMOVE; 
        const int removeArrayIndex = atomicAdd(&remove_count, 1);
        remove_particle_indices[removeArrayIndex] = particle_index;

        return false;
    }

    // cell_id is an own cell and move_flag == OPP_NEED_MOVE
    return true;
}

}

//--------------------------------------------------------------
__global__ void opp_dev_move_kernel(
    const OPP_REAL *__restrict__ dat0,     // p_pos
    OPP_INT *__restrict__ dat1,     // p_mdir
    const OPP_REAL *__restrict__ dat2,     // c_pos_ll
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
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;

        OPP_INT *opp_p2c = (p2c_map + n);
        char move_flag = OPP_NEED_MOVE;
        bool iter_one_flag = (OPP_comm_iteration_d > 0) ? false : true;

        do
        {
            const OPP_INT p2c = opp_p2c[0]; // get the value here, since the kernel might change it
            const OPP_INT* opp_c2c = c2c_map + p2c;           

            opp_k2::move_kernel(
                move_flag, iter_one_flag, opp_c2c, opp_p2c,
                dat0 + n, // p_pos 
                dat1 + n, // p_mdir 
                dat2 + p2c // c_pos_ll 
          
            );

        } while (opp_k2::opp_part_check_status_hip(move_flag, iter_one_flag, opp_p2c, n, 
            *particle_remove_count, particle_remove_indices, move_particle_indices, 
            move_cell_indices, move_count));        
    }
    
}

//--------------------------------------------------------------

void opp_particle_move__move_kernel(opp_set set, opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0,   // p_pos | OPP_READ
    opp_arg arg1,   // p_mdir | OPP_RW
    opp_arg arg2   // c_pos_ll | OPP_READ
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

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
    if (OPP_cells_set_size != set->cells_set->size) {
        OPP_cells_set_size = set->cells_set->size; 
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(OPP_cells_set_size_d), &OPP_cells_set_size, sizeof(int)));
    }
    const OPP_INT c2c_stride = c2c_map->from->size + c2c_map->from->exec_size + c2c_map->from->nonexec_size;
    if (opp_k2_c2c_map_stride != c2c_stride) {
        opp_k2_c2c_map_stride = c2c_stride;
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k2_c2c_map_stride_d), &opp_k2_c2c_map_stride, sizeof(OPP_INT)));
    }

    opp_mpi_halo_wait_all(nargs, args);

#ifdef OPP_BLOCK_SIZE_2
    const int block_size = OPP_BLOCK_SIZE_2;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    int num_blocks = 200;

    do 
    {
        if (opp_k2_dat0_stride != args[0].dat->set->set_capacity) {
            opp_k2_dat0_stride = args[0].dat->set->set_capacity;
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k2_dat0_stride_d), &opp_k2_dat0_stride, sizeof(OPP_INT)));
        }
        if (opp_k2_dat1_stride != args[1].dat->set->set_capacity) {
            opp_k2_dat1_stride = args[1].dat->set->set_capacity;
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k2_dat1_stride_d), &opp_k2_dat1_stride, sizeof(OPP_INT)));
        }
        if (opp_k2_dat2_stride != args[2].dat->set->set_capacity) {
            opp_k2_dat2_stride = args[2].dat->set->set_capacity;
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k2_dat2_stride_d), &opp_k2_dat2_stride, sizeof(OPP_INT)));
        }

        opp_init_particle_move(set, nargs, args);
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(OPP_comm_iteration_d), &OPP_comm_iteration, sizeof(int)));

        num_blocks = (OPP_iter_end - OPP_iter_start - 1) / block_size + 1;

        {
            opp_dev_move_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,    // p_pos
                (OPP_INT *)args[1].data_d,    // p_mdir
                (OPP_REAL *)args[2].data_d,    // c_pos_ll
                (OPP_INT *)args[3].data_d,    // p2c_map
                (OPP_INT *)c2c_map->map_d,    // c2c_map
                (OPP_INT *)set->particle_remove_count_d,
                (OPP_INT *)OPP_remove_particle_indices_d,
                (OPP_INT *)OPP_move_particle_indices_d,
                (OPP_INT *)OPP_move_cell_indices_d,
                (OPP_INT *)OPP_move_count_d,
                OPP_iter_start,
                OPP_iter_end
            );
        }

    } while (opp_finalize_particle_move(set)); 

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("move_kernel");
}
void opp_init_direct_hop_cg(double grid_spacing, int dim, const opp_dat c_gbl_id, const opp::BoundingBox& b_box, 
    opp_map c2c_map, opp_map p2c_map,
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1, // p_mdir | OPP_RW
    opp_arg arg2 // c_pos_ll | OPP_READ
) {
    opp_profiler->start("Setup_Mover");
    
    opp_profiler->end("Setup_Mover");
}

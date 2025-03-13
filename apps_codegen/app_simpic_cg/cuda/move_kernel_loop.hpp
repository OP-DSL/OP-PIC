
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
enum Dir : char {
    Left = 0,
    Right,
};

__device__ void move_kernel(
        char& opp_move_status_flag, const bool opp_move_hop_iter_one_flag, // Added by code-gen
    const OPP_INT* opp_c2c, OPP_INT* opp_p2c, // Added by code-gen
    const double* part_field_E,
        double* part_velocity_x,
        double* part_position_x
    )
{
    if ((opp_move_hop_iter_one_flag))
    {
        part_velocity_x[(0) * opp_k2_dat1_stride_d] += (CONST_qm_d[0] * part_field_E[(0) * opp_k2_dat0_stride_d]);
        part_position_x[(0) * opp_k2_dat2_stride_d] += (part_velocity_x[(0) * opp_k2_dat1_stride_d] * CONST_dt_d[0]);
    }

    // since particle cell index can be easily calculated with global positioning, no need to search by iterating
    if ((part_position_x[(0) * opp_k2_dat2_stride_d] > CONST_xl_d[0]) && (part_position_x[(0) * opp_k2_dat2_stride_d] < CONST_xr_d[0]))
    {
        double xx = ((part_position_x[(0) * opp_k2_dat2_stride_d] - CONST_xl_d[0]) / CONST_dx_d[0]); // Makes Global position to local position comapared to the cell
        opp_p2c[0] = int(xx);
        { opp_move_status_flag = OPP_MOVE_DONE; };
    }
    else if ((part_position_x[(0) * opp_k2_dat2_stride_d] >= CONST_xl_d[0]) && (CONST_rank_d[0] == 0) ||
             (part_position_x[(0) * opp_k2_dat2_stride_d] <= CONST_xr_d[0]) && (CONST_rank_d[0] == (CONST_comm_size_d[0]-1)))
    {
        opp_p2c[0] = INT_MAX;
        { opp_move_status_flag = OPP_NEED_REMOVE; };
    }
    else
    {
        opp_p2c[0] = (part_position_x[(0) * opp_k2_dat2_stride_d] < CONST_xl_d[0]) ?
                        CONST_neighbour_cell_d[Dir::Left] : CONST_neighbour_cell_d[Dir::Right];
        { opp_move_status_flag = OPP_NEED_MOVE; };
    }
}

//*******************************************************************************
// Returns true only if another hop is required by the current rank
__device__ inline bool opp_part_check_status_cuda(char& move_flag, bool& iter_one_flag, 
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

__global__ void opp_dev_move_kernel(
    const OPP_REAL *__restrict__ dat0,     // p_field_e
    OPP_REAL *__restrict__ dat1,     // p_vel_x
    OPP_REAL *__restrict__ dat2,     // p_pos_x
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
                dat0 + n, // p_field_e 
                dat1 + n, // p_vel_x 
                dat2 + n // p_pos_x 
          
            );

        } while (opp_k2::opp_part_check_status_cuda(move_flag, iter_one_flag, opp_p2c, n, 
            *particle_remove_count, particle_remove_indices, move_particle_indices, 
            move_cell_indices, move_count));        
    }
}

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

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
    if (OPP_cells_set_size != set->cells_set->size) {
        OPP_cells_set_size = set->cells_set->size; 
        cutilSafeCall(cudaMemcpyToSymbol(OPP_cells_set_size_d, &OPP_cells_set_size, sizeof(int)));
    }
    const OPP_INT c2c_stride = c2c_map->from->size + c2c_map->from->exec_size + c2c_map->from->nonexec_size;
    if (opp_k2_c2c_map_stride != c2c_stride) {
        opp_k2_c2c_map_stride = c2c_stride;
        cutilSafeCall(cudaMemcpyToSymbol(opp_k2_c2c_map_stride_d, &opp_k2_c2c_map_stride, sizeof(OPP_INT)));
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
            cutilSafeCall(cudaMemcpyToSymbol(opp_k2_dat0_stride_d, &opp_k2_dat0_stride, sizeof(OPP_INT)));
        }
        if (opp_k2_dat1_stride != args[1].dat->set->set_capacity) {
            opp_k2_dat1_stride = args[1].dat->set->set_capacity;
            cutilSafeCall(cudaMemcpyToSymbol(opp_k2_dat1_stride_d, &opp_k2_dat1_stride, sizeof(OPP_INT)));
        }
        if (opp_k2_dat2_stride != args[2].dat->set->set_capacity) {
            opp_k2_dat2_stride = args[2].dat->set->set_capacity;
            cutilSafeCall(cudaMemcpyToSymbol(opp_k2_dat2_stride_d, &opp_k2_dat2_stride, sizeof(OPP_INT)));
        }

        opp_init_particle_move(set, nargs, args);
        cutilSafeCall(cudaMemcpyToSymbol(OPP_comm_iteration_d, &OPP_comm_iteration, sizeof(int)));

        num_blocks = (OPP_iter_end - OPP_iter_start - 1) / block_size + 1;

        opp_dev_move_kernel<<<num_blocks, block_size>>>(
            (OPP_REAL *)args[0].data_d,    // p_field_e
            (OPP_REAL *)args[1].data_d,    // p_vel_x
            (OPP_REAL *)args[2].data_d,    // p_pos_x
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

    } while (opp_finalize_particle_move(set)); 

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("move_kernel");
}

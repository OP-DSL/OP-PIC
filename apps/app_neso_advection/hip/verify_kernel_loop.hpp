
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k3_dat0_stride = -1;
OPP_INT opp_k3_dat1_stride = -1;

__constant__ OPP_INT opp_k3_dat0_stride_d;
__constant__ OPP_INT opp_k3_dat1_stride_d;



namespace opp_k3 {
enum Dim {
    x = 0,
    y = 1,
};

__device__ inline void verify_kernel(
        const double* p_pos,
        const int* c_gbl_idx,
        int* incorrect_count)
{
    // get the cell boundaries for the current cell_index - using global cell index
    int ix = -1, iy = -1;
    int _ix, _iy; _ix  = ((*c_gbl_idx)); _iy  = _ix/int(CONST_ndimcells_d[Dim::x]); _ix -= _iy*int(CONST_ndimcells_d[Dim::x]); (ix) = _ix; (iy) = _iy;;

    if (ix < 0 || iy < 0)
    {
        // opp_printf("VERIFY", "Incorrect ix[%d] iy[%d] for global cell[%d] nx[%d]",
        //     ix, iy, (*c_gbl_idx), CONST_ndimcells[Dim::x]);
        (*incorrect_count)++;
        return;
    }

    // get the boundaries of that cell
    const double boundary_ll[2] = { (ix * CONST_cell_width_d[0]), (iy * CONST_cell_width_d[0]) };

    // check whether the current particle is within those boundaries or not!
    const double p_pos_x = p_pos[(Dim::x) * opp_k3_dat0_stride_d];
    if (p_pos_x < boundary_ll[Dim::x] ||
            p_pos_x > (boundary_ll[Dim::x] + CONST_cell_width_d[0])) {

        (*incorrect_count)++;
        return;
    }

    const double p_pos_y = p_pos[(Dim::y) * opp_k3_dat0_stride_d];
    if (p_pos_y < boundary_ll[Dim::y] ||
            p_pos_y > (boundary_ll[Dim::y] + CONST_cell_width_d[0])) {

        (*incorrect_count)++;
        return;
    }
}

}

//--------------------------------------------------------------
__global__ void opp_dev_verify_kernel(
    const OPP_REAL *__restrict__ dat0,  // p_pos
    const OPP_INT *__restrict__ dat1,  // c_idx
    const OPP_INT *__restrict__ p2c_map,
    OPP_INT *gbl2,
    const OPP_INT start,
    const OPP_INT end
) 
{
    OPP_INT gbl2_local[1];
    for (int d = 0; d < 1; ++d)
        gbl2_local[d] = OPP_INT_ZERO;

    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    for (int n = thread_id; n < (end - start); n += blockDim.x * gridDim.x) {
        
        const OPP_INT* opp_p2c = p2c_map + n;
        
        opp_k3::verify_kernel(
            dat0 + n, // p_pos 
            dat1 + opp_p2c[0], // c_idx 
            gbl2_local // 
        );
    }

    for (int d = 0; d < 1; ++d)
        opp_reduction<OPP_INC>(gbl2 + blockIdx.x * 1 + d, gbl2_local[d]);
    
}

//--------------------------------------------------------------
void opp_par_loop_all__verify_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // p_pos | OPP_READ
    opp_arg arg1, // c_idx | OPP_READ
    opp_arg arg2 // | OPP_INC
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("verify_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__verify_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    OPP_INT *arg2_host_data = (OPP_INT *)args[2].data;

    if (opp_k3_dat0_stride != args[0].dat->set->set_capacity) {
        opp_k3_dat0_stride = args[0].dat->set->set_capacity;
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k3_dat0_stride_d), &opp_k3_dat0_stride, sizeof(OPP_INT)));
    }
    if (opp_k3_dat1_stride != args[1].dat->set->set_capacity) {
        opp_k3_dat1_stride = args[1].dat->set->set_capacity;
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k3_dat1_stride_d), &opp_k3_dat1_stride, sizeof(OPP_INT)));
    }

#ifdef OPP_BLOCK_SIZE_3
    const int block_size = OPP_BLOCK_SIZE_3;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    int num_blocks = 200;

    int max_blocks = (MAX(set->core_size, set->size + set->exec_size - set->core_size) - 1) / block_size + 1;

    int reduction_bytes = 0;
    int reduction_size = 0;

    reduction_bytes += ROUND_UP(max_blocks * 1 * sizeof(OPP_INT));
    reduction_size   = MAX(reduction_size, sizeof(OPP_INT));

    opp_reallocReductArrays(reduction_bytes);
    reduction_bytes = 0;

    args[2].data   = OPP_reduct_h + reduction_bytes;
    args[2].data_d = OPP_reduct_d + reduction_bytes;

    for (int b = 0; b < max_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            ((OPP_INT *)args[2].data)[b * 1 + d] = OPP_INT_ZERO;
    }

    reduction_bytes += ROUND_UP(max_blocks * 1 * sizeof(OPP_INT));

    opp_mvReductArraysToDevice(reduction_bytes);
    
    if (iter_size > 0) 
    {
        const OPP_INT start = 0;
        const OPP_INT end = iter_size;
        {
            opp_dev_verify_kernel<<<num_blocks, block_size, (reduction_size * block_size)>>>(
                (OPP_REAL *)args[0].data_d,     // p_pos
                (OPP_INT *)args[1].data_d,     // c_idx
                (OPP_INT *)set->mesh_relation_dat->data_d,
                (OPP_INT *)args[2].data_d,
                start,
                end
            );
        }
    }

    opp_mvReductArraysToHost(reduction_bytes);

    for (int b = 0; b < max_blocks; ++b) {
        for (int d = 0; d < 1; ++d)
            arg2_host_data[d] += ((OPP_INT *)args[2].data)[b * 1 + d];
    }

    args[2].data = (char *)arg2_host_data;
    opp_mpi_reduce(&args[2], arg2_host_data);

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   
 
    opp_profiler->end("verify_kernel");
}

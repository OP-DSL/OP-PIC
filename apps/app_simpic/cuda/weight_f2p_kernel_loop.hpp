
//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k1_dat0_stride = -1;
OPP_INT opp_k1_dat1_stride = -1;
OPP_INT opp_k1_dat2_stride = -1;
OPP_INT opp_k1_map0_stride = -1;

__constant__ OPP_INT opp_k1_dat0_stride_d;
__constant__ OPP_INT opp_k1_dat1_stride_d;
__constant__ OPP_INT opp_k1_dat2_stride_d;
__constant__ OPP_INT opp_k1_map0_stride_d;



namespace opp_k1 {
__device__ void weight_f2p_kernel(
        const double* node0_field_E,  //LHS
        const double* node1_field_E,  //RHS
        const double* particle0_position_x,
        double* particle0_field_E
    )
{
    double xx = ((particle0_position_x[(0) * opp_k1_dat1_stride_d] - CONST_xl_d[0]) / CONST_dx_d[0]); // Makes Global position to local position comapared to the cell
    int n = int(xx);
    double frac = (xx - n);

    particle0_field_E[(0) * opp_k1_dat2_stride_d] = ((frac * node1_field_E[(0) * opp_k1_dat0_stride_d]) + ((1.0 - frac) * node0_field_E[(0) * opp_k1_dat0_stride_d]));
}

}

//--------------------------------------------------------------
__global__ void opp_dev_weight_f2p_kernel(
    const OPP_REAL *__restrict__ dat0,  // n_field_e
    const OPP_REAL *__restrict__ dat1,  // p_pos_x
    OPP_REAL *__restrict__ dat2,  // p_field_e
    const OPP_INT *__restrict__ p2c_map,
    const OPP_INT *__restrict__ map0,  // c2n_map
    const OPP_INT start,
    const OPP_INT end
) 
{
    const int thread_id = threadIdx.x + blockIdx.x * blockDim.x;

    if (thread_id + start < end) {

        const int n = thread_id + start;
        
        const OPP_INT* opp_p2c = p2c_map + n;
        
        opp_k1::weight_f2p_kernel(
            dat0 + map0[opp_p2c[0] + opp_k1_map0_stride_d * 0], // n_field_e 
            dat0 + map0[opp_p2c[0] + opp_k1_map0_stride_d * 1], // n_field_e 
            dat1 + n, // p_pos_x 
            dat2 + n // p_field_e 
        );
    }
    
}

//--------------------------------------------------------------
void opp_par_loop_all__weight_f2p_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // n_field_e | OPP_READ
    opp_arg arg1, // n_field_e | OPP_READ
    opp_arg arg2, // p_pos_x | OPP_READ
    opp_arg arg3 // p_field_e | OPP_WRITE
) 
{
    const int nargs = 4;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;
    args[3] = arg3;

    opp_profiler->start("weight_f2p_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__weight_f2p_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
 
 
    if (opp_k1_dat0_stride != args[0].dat->set->set_capacity) {
        opp_k1_dat0_stride = args[0].dat->set->set_capacity;
        cutilSafeCall(cudaMemcpyToSymbol(opp_k1_dat0_stride_d, &opp_k1_dat0_stride, sizeof(OPP_INT)));
    }
    if (opp_k1_dat1_stride != args[2].dat->set->set_capacity) {
        opp_k1_dat1_stride = args[2].dat->set->set_capacity;
        cutilSafeCall(cudaMemcpyToSymbol(opp_k1_dat1_stride_d, &opp_k1_dat1_stride, sizeof(OPP_INT)));
    }
    if (opp_k1_dat2_stride != args[3].dat->set->set_capacity) {
        opp_k1_dat2_stride = args[3].dat->set->set_capacity;
        cutilSafeCall(cudaMemcpyToSymbol(opp_k1_dat2_stride_d, &opp_k1_dat2_stride, sizeof(OPP_INT)));
    }
    if (opp_k1_map0_stride != args[0].size) {
        opp_k1_map0_stride = args[0].size;
        cutilSafeCall(cudaMemcpyToSymbol(opp_k1_map0_stride_d, &opp_k1_map0_stride, sizeof(OPP_INT)));
    }

#ifdef OPP_BLOCK_SIZE_1
    const int block_size = OPP_BLOCK_SIZE_1;
#else
    const int block_size = OPP_gpu_threads_per_block;
#endif

    opp_mpi_halo_wait_all(nargs, args);

    int num_blocks = 200;
    if (iter_size > 0) 
    {
        const OPP_INT start = 0;
        const OPP_INT end = iter_size;
        num_blocks = (end - start - 1) / block_size + 1;

        {
            opp_dev_weight_f2p_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // n_field_e
                (OPP_REAL *)args[2].data_d,     // p_pos_x
                (OPP_REAL *)args[3].data_d,     // p_field_e
                (OPP_INT *)set->mesh_relation_dat->data_d,
                args[0].map_data_d,     // c2n_map
                start,
                end
            );
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    cutilSafeCall(cudaDeviceSynchronize());   
 
    opp_profiler->end("weight_f2p_kernel");
}


//*********************************************
// AUTO GENERATED CODE
//*********************************************

OPP_INT opp_k3_dat0_stride = -1;
OPP_INT opp_k3_dat1_stride = -1;
OPP_INT opp_k3_map0_stride = -1;

__constant__ OPP_INT opp_k3_dat0_stride_d;
__constant__ OPP_INT opp_k3_dat1_stride_d;
__constant__ OPP_INT opp_k3_map0_stride_d;

OPP_INT opp_k3_sr_set_stride = -1;
__constant__ OPP_INT opp_k3_sr_set_stride_d;

thrust::device_vector<OPP_INT> sr_dat0_keys_dv;
thrust::device_vector<OPP_REAL> sr_dat0_values_dv;

namespace opp_k3 {
__device__ void weight_p2f_kernel(
        double* node0_field_J,
        double* node1_field_J,
        const double* particle0_position_x
    )
{
    double xx = ((particle0_position_x[(0) * opp_k3_dat1_stride_d] - CONST_xl_d[0]) / CONST_dx_d[0]); // Makes Global position to local position comapared to the cell
    int n = int(xx);
    double frac = (xx - n);

    (*node0_field_J) += (CONST_qscale_d[0] * (1.0 - frac));  // Can change qscale to be calculated from particle data
    (*node1_field_J) += (CONST_qscale_d[0] * frac);
}

//--------------------------------------------------------------
__global__ void assign_values( // Used for Segmented Reductions
    const OPP_INT *__restrict keys,
    const OPP_REAL *__restrict values,
    OPP_REAL *__restrict dat,
    const int start,
    const int end) 
{
    const int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        const int n = tid + start;
        const int mapping = keys[n];  
        dat[mapping] += values[n];
    }
}
}

//--------------------------------------------------------------
__global__ void opp_dev_weight_p2f_kernel(
    OPP_REAL *__restrict__ dat0,  // n_field_j
    const OPP_REAL *__restrict__ dat1,  // p_pos_x
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
        
        OPP_REAL arg0_0_local[1];
        for (int d = 0; d < 1; ++d)
            arg0_0_local[d] = OPP_REAL_ZERO;

        OPP_REAL arg1_1_local[1];
        for (int d = 0; d < 1; ++d)
            arg1_1_local[d] = OPP_REAL_ZERO;

        opp_k3::weight_p2f_kernel(
            arg0_0_local, // n_field_j 
            arg1_1_local, // n_field_j 
            dat1 + n // p_pos_x 
        );

        for (int d = 0; d < 1; ++d)
            atomicAdd(dat0 + map0[opp_k3_map0_stride_d * 0 + opp_p2c[0]] + (d * opp_k3_dat0_stride_d), arg0_0_local[d]);

        for (int d = 0; d < 1; ++d)
            atomicAdd(dat0 + map0[opp_k3_map0_stride_d * 1 + opp_p2c[0]] + (d * opp_k3_dat0_stride_d), arg1_1_local[d]);
    }
    
}

//--------------------------------------------------------------
__global__ void opp_dev_sr_weight_p2f_kernel( // Used for Segmented Reductions
    OPP_INT *__restrict__ sr_dat0_keys,     // sr keys for n_field_j
    OPP_REAL *__restrict__ sr_dat0_values,     // sr values for n_field_j
    const OPP_REAL *__restrict__ dat1,  // p_pos_x
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

        OPP_REAL arg0_0_local[1];
        for (int d = 0; d < 1; ++d)
            arg0_0_local[d] = OPP_REAL_ZERO;

        OPP_REAL arg1_1_local[1];
        for (int d = 0; d < 1; ++d)
            arg1_1_local[d] = OPP_REAL_ZERO;

        opp_k3::weight_p2f_kernel(
            arg0_0_local, // n_field_j 
            arg1_1_local, // n_field_j 
            dat1 + n // p_pos_x 
        );

        int offset = 0;
        for (int d = 0; d < 1; ++d, ++offset) {
            sr_dat0_values[n + opp_k3_sr_set_stride_d * offset] = arg0_0_local[d];
            sr_dat0_keys[n + opp_k3_sr_set_stride_d * offset] = map0[opp_k3_map0_stride_d * 0 + opp_p2c[0]] + (d * opp_k3_dat0_stride_d);
        }
        for (int d = 0; d < 1; ++d, ++offset) {
            sr_dat0_values[n + opp_k3_sr_set_stride_d * offset] = arg1_1_local[d];
            sr_dat0_keys[n + opp_k3_sr_set_stride_d * offset] = map0[opp_k3_map0_stride_d * 1 + opp_p2c[0]] + (d * opp_k3_dat0_stride_d);
        }
    }
    
}

void opp_par_loop_all__weight_p2f_kernel(opp_set set, opp_iterate_type, 
    opp_arg arg0, // n_field_j | OPP_INC
    opp_arg arg1, // n_field_j | OPP_INC
    opp_arg arg2 // p_pos_x | OPP_READ
) 
{
    const int nargs = 3;
    opp_arg args[nargs];

    args[0] = arg0;
    args[1] = arg1;
    args[2] = arg2;

    opp_profiler->start("weight_p2f_kernel");

    if (OPP_DBG) opp_printf("APP", "opp_par_loop_all__weight_p2f_kernel set_size %d", set->size);

    const int iter_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);

#ifdef USE_MPI
    opp_init_double_indirect_reductions_device(nargs, args);
#endif
 
 
    if (opp_k3_dat0_stride != args[0].dat->set->set_capacity) {
        opp_k3_dat0_stride = args[0].dat->set->set_capacity;
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k3_dat0_stride_d), &opp_k3_dat0_stride, sizeof(OPP_INT)));
    }
    if (opp_k3_dat1_stride != args[2].dat->set->set_capacity) {
        opp_k3_dat1_stride = args[2].dat->set->set_capacity;
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k3_dat1_stride_d), &opp_k3_dat1_stride, sizeof(OPP_INT)));
    }
    if (opp_k3_map0_stride != args[0].size) {
        opp_k3_map0_stride = args[0].size;
        cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k3_map0_stride_d), &opp_k3_map0_stride, sizeof(OPP_INT)));
    }

#ifdef OPP_BLOCK_SIZE_3
    const int block_size = OPP_BLOCK_SIZE_3;
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

        if (!opp_params->get<OPP_BOOL>("use_reg_red")) // Do atomics ----------       
        {
            opp_dev_weight_p2f_kernel<<<num_blocks, block_size>>>(
                (OPP_REAL *)args[0].data_d,     // n_field_j
                (OPP_REAL *)args[2].data_d,     // p_pos_x
                (OPP_INT *)set->mesh_relation_dat->data_d,
                args[0].map_data_d,     // c2n_map
                start,
                end
            );
        }
     
        else // Do segmented reductions ----------       
        {
            if (opp_k3_sr_set_stride != set->size) {
                opp_k3_sr_set_stride = set->size;
                cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(opp_k3_sr_set_stride_d), &opp_k3_sr_set_stride, sizeof(OPP_INT)));
            }

            size_t operating_size_dat0 = 0, resize_size_dat0 = 0;

            operating_size_dat0 += (size_t)(args[0].dat->dim);
            resize_size_dat0 += (size_t)(args[0].dat->dim);
            operating_size_dat0 += (size_t)(args[1].dat->dim);
            resize_size_dat0 += (size_t)(args[1].dat->dim);

            operating_size_dat0 *= (size_t)(set->size);
            resize_size_dat0 *= (size_t)(set->set_capacity);

            if (resize_size_dat0 > sr_dat0_keys_dv.size()) { // resize only if current vector is small        
                sr_dat0_keys_dv.resize(resize_size_dat0, 0);
                sr_dat0_values_dv.resize(resize_size_dat0, 0);
            }
        
            // Create key/value pairs
            opp_profiler->start("SR_CrKeyVal");
            opp_dev_sr_weight_p2f_kernel<<<num_blocks, block_size>>>( 
                opp_get_dev_raw_ptr<OPP_INT>(sr_dat0_keys_dv),     // sr keys for n_field_j
                opp_get_dev_raw_ptr<OPP_REAL>(sr_dat0_values_dv),     // sr values for n_field_j
                (OPP_REAL *)args[2].data_d,     // p_pos_x
                (OPP_INT *)set->mesh_relation_dat->data_d,
                args[0].map_data_d,     // c2n_map
                start,
                end
            );
            OPP_DEVICE_SYNCHRONIZE();
            opp_profiler->end("SR_CrKeyVal");

            // Sort by keys to bring the identical keys together
            opp_profiler->start("SR_SortByKey");
            thrust::sort_by_key(sr_dat0_keys_dv.begin(), sr_dat0_keys_dv.begin() + operating_size_dat0, 
                sr_dat0_values_dv.begin());
            opp_profiler->end("SR_SortByKey");

            // Compute the unique keys and their corresponding values
            opp_profiler->start("SR_RedByKey");
            auto new_end = thrust::reduce_by_key(
                sr_dat0_keys_dv.begin(), sr_dat0_keys_dv.begin() + operating_size_dat0,
                sr_dat0_values_dv.begin(),
                sr_dat0_keys_dv.begin(),
                sr_dat0_values_dv.begin());        
            opp_profiler->end("SR_RedByKey");

            const size_t reduced_size = (new_end.first - sr_dat0_keys_dv.begin());
            
            // Assign reduced values to the nodes using keys/values
            opp_profiler->start("SR_Assign");                
            opp_k3::assign_values<<<num_blocks, block_size>>> ( // TODO : check whether num_blocks is correct
                opp_get_dev_raw_ptr<OPP_INT>(sr_dat0_keys_dv),
                opp_get_dev_raw_ptr<OPP_REAL>(sr_dat0_values_dv),
                (OPP_REAL *) args[0].data_d,
                0, reduced_size);
            OPP_DEVICE_SYNCHRONIZE();
            opp_profiler->end("SR_Assign");

            // Last: clear the thrust vectors if this is the last iteration (avoid crash)
            if (opp_params->get<OPP_INT>("num_steps") == (OPP_main_loop_iter + 1)) {
                OPP_DEVICE_SYNCHRONIZE();
                sr_dat0_values_dv.clear(); sr_dat0_values_dv.shrink_to_fit();
                sr_dat0_keys_dv.clear(); sr_dat0_keys_dv.shrink_to_fit();
            }        
        }
    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();   

#ifdef USE_MPI    
    opp_exchange_double_indirect_reductions_device(nargs, args);
    opp_complete_double_indirect_reductions_device(nargs, args);
#endif
 
    opp_profiler->end("weight_p2f_kernel");
}

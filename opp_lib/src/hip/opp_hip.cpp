/* 
BSD 3-Clause License

Copyright (c) 2022, OP-DSL

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "opp_hip.h"
#include <unistd.h>
#include "opp_particle_comm.cpp"
#include "opp_increase_part_count.cpp"

// arrays for global constants and reductions
int OPP_consts_bytes = 0, OPP_reduct_bytes = 0;
char *OPP_reduct_h = nullptr;
char *OPP_reduct_d = nullptr;
char *OPP_consts_h = nullptr;
char *OPP_consts_d = nullptr;

//****************************************
void opp_init(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << std::endl;
        exit(-1);
    }

#ifdef USE_PETSC
    PetscInitialize(&argc, &argv, PETSC_NULLPTR, "opp::PetscHIP");
#else
    #ifdef USE_MPI
        MPI_Init(&argc, &argv);
    #endif
#endif

#ifdef USE_MPI
    OPP_MPI_WORLD = MPI_COMM_WORLD;
    
    MPI_Comm_rank(OPP_MPI_WORLD, &OPP_rank);
    MPI_Comm_size(OPP_MPI_WORLD, &OPP_comm_size);
#endif

    if (OPP_rank == OPP_ROOT)
    {
        std::string log = "Running on HIP";
#ifdef USE_MPI
        log += "+MPI with " + std::to_string(OPP_comm_size) + " ranks";
#endif        
        opp_printf("OP-PIC", "%s", log.c_str());
        opp_printf("OP-PIC", "---------------------------------------------");
    }

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif  

    opp_init_core(argc, argv);
    opp_params->write(std::cout);
    cutilDeviceInit(argc, argv);

    // cutilSafeCall(hipDeviceSetCacheConfig(hipFuncCachePreferL1));
    // cutilSafeCall(hipDeviceSetCacheConfig(hipFuncCachePreferShared));
    // cutilSafeCall(hipDeviceSetSharedMemConfig(hipSharedMemBankSizeEightByte));

    OPP_auto_soa = 1; // TODO : Make this configurable with args

    const int threads_per_block = opp_params->get<OPP_INT>("opp_threads_per_block");   
    if (threads_per_block > 0 && threads_per_block < INT_MAX)
        OPP_gpu_threads_per_block = threads_per_block;

    const int gpu_direct = opp_params->get<OPP_INT>("opp_gpu_direct");   
    if (gpu_direct > 0 && gpu_direct < INT_MAX)
        OPP_gpu_direct = gpu_direct;

    int deviceId = -1;
    cutilSafeCall(hipGetDevice(&deviceId));
    hipDeviceProp_t prop;
    cutilSafeCall(hipGetDeviceProperties(&prop, deviceId));

    OPP_gpu_shared_mem_per_block = prop.sharedMemPerBlock;

    char hostname[256];
    if (gethostname(hostname, sizeof(hostname)) != 0) 
    {
        opp_printf("cutilDeviceInit", "Failed to get hostname of MPI rank %d", OPP_rank);
        opp_abort();
    }

    opp_printf("opp_init", 
        "Device: %d [%s] on Host [%s] threads=%d Shared_Mem=%lubytes GPU_Direct=%d", deviceId, 
        prop.name, hostname, OPP_gpu_threads_per_block, prop.sharedMemPerBlock, OPP_gpu_direct);

    opp_profiler->reg("PSS_Resize");
    opp_profiler->reg("PSS_Copy");
    opp_profiler->reg("PSS_HoleFill");
    opp_profiler->reg("PSS_Sequence");
    opp_profiler->reg("PSS_SortKey0");
    opp_profiler->reg("PSS_SortKey");
    opp_profiler->reg("PSS_Dats");
    opp_profiler->reg("PHF_Sort");
    opp_profiler->reg("PHF_CopyIf");
    opp_profiler->reg("PHF_Dats");

    opp_profiler->reg("Inc_part_count");
    opp_profiler->reg("Inc_part_count_DWN");
    opp_profiler->reg("Inc_part_count_INC");
    opp_profiler->reg("Inc_part_count_UPL");
    opp_profiler->reg("Mv_Pack");
    opp_profiler->reg("Mv_Unpack");
    opp_profiler->reg("Mv_PackExDir");
    opp_profiler->reg("Mv_UnpackDir");
    opp_profiler->reg("Mv_Finalize"); 
    opp_profiler->reg("Mv_fill");
}

//****************************************
void opp_exit()
{
    globalMover.reset();
    cellMapper.reset();
    boundingBox.reset();
    comm.reset();

    if (OPP_reduct_h) free(OPP_reduct_h);
    if (OPP_reduct_d) cutilSafeCall(hipFree(OPP_reduct_d));
    if (OPP_consts_h) free(OPP_consts_h);
    if (OPP_consts_d) cutilSafeCall(hipFree(OPP_consts_d));

#ifdef USE_MPI 
        opp_halo_destroy(); // free memory allocated to halos and mpi_buffers 
        opp_partition_destroy(); // free memory used for holding partition information
        opp_part_comm_destroy(); // free memory allocated for particle communication
#endif

    opp_hip_exit();
    opp_exit_core();

#ifdef USE_PETSC
    PetscFinalize();
#endif
}

//****************************************
void opp_abort(std::string s)
{
    opp_printf("opp_abort", "%s", s.c_str());
#ifdef USE_MPI 
    MPI_Abort(OPP_MPI_WORLD, 2);
#else
    exit(-1);
#endif
}

//****************************************
void opp_hip_exit() 
{
    if (!OPP_hybrid_gpu)
        return;

    for (auto& a : opp_maps) 
    {
        cutilSafeCall(hipFree(a->map_d));
    }

    for (auto& a : opp_dats) 
    {
        // cutilSafeCall(hipFree(a->data_d));
        if (a->thrust_int) delete a->thrust_int;
        if (a->thrust_real) delete a->thrust_real;
        if (a->thrust_int_sort) delete a->thrust_int_sort;
        if (a->thrust_real_sort) delete a->thrust_real_sort;
    }

    if (opp_saved_mesh_relation_d != nullptr)
    {
        cutilSafeCall(hipFree(opp_saved_mesh_relation_d));
        // opp_host_free(opp_saved_mesh_relation_d);
    }

    cellIdx_dv.clear();
    cellIdx_dv.shrink_to_fit();

    i_dv.clear();
    i_dv.shrink_to_fit();

    send_part_cell_idx_dv.clear();
    send_part_cell_idx_dv.shrink_to_fit();

    temp_int_dv.clear();
    temp_int_dv.shrink_to_fit();

    temp_real_dv.clear();
    temp_real_dv.shrink_to_fit();

    OPP_thrust_move_particle_indices_d.clear();
    OPP_thrust_move_particle_indices_d.shrink_to_fit();

    OPP_thrust_move_cell_indices_d.clear();
    OPP_thrust_move_cell_indices_d.shrink_to_fit();

    OPP_thrust_remove_particle_indices_d.clear();
    OPP_thrust_remove_particle_indices_d.shrink_to_fit();

    ps_to_indices_dv.clear(); ps_to_indices_dv.shrink_to_fit(); 
    ps_from_indices_dv.clear(); ps_from_indices_dv.shrink_to_fit(); 
    ps_sequence_dv.clear(); ps_sequence_dv.shrink_to_fit(); 

    for (auto it = particle_indices_hv.begin(); it != particle_indices_hv.end(); it++) it->second.clear();
    for (auto it = cell_indices_hv.begin(); it != cell_indices_hv.end(); it++) it->second.clear();
    for (auto it = particle_indices_dv.begin(); it != particle_indices_dv.end(); it++)
    {
        it->second.clear();
        it->second.shrink_to_fit();
    }
    for (auto it = send_data.begin(); it != send_data.end(); it++)
    {
        it->second.clear();
        it->second.shrink_to_fit();
    }
    for (auto it = recv_data.begin(); it != recv_data.end(); it++)
    {
        it->second.clear();
        it->second.shrink_to_fit();
    }

    if (OPP_need_remove_flags_d != nullptr)
    {
        cutilSafeCall(hipFree(OPP_need_remove_flags_d));
    }

    if (OPP_move_count_d != nullptr)
    {
        cutilSafeCall(hipFree(OPP_move_count_d));
    } 
}

//****************************************
opp_set opp_decl_set(int size, char const *name)
{
    return opp_decl_set_core(size, name);
}

//****************************************
opp_map opp_decl_map(opp_set from, opp_set to, int dim, int *imap, char const *name)
{
    opp_map map = opp_decl_map_core(from, to, dim, imap, name);

    if (from->is_particle)
    {
        opp_create_device_arrays(map->p2c_dat);
        opp_upload_dat(map->p2c_dat);
    }
    else
    {
        opp_upload_map(map, true);
    }

    return map;
}

//****************************************
opp_dat opp_decl_mesh_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name)
{
    std::string type = "";
    int size = -1;

    getDatTypeSize(dtype, type, size);
    
    opp_dat dat = opp_decl_dat_core(set, dim, type.c_str(), size, (char*)data, name);

    opp_create_device_arrays(dat);

    opp_upload_dat(dat);

    return dat;
}

opp_dat opp_decl_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name, bool cell_index)
{
    if (set->is_particle) 
        return opp_decl_part_dat(set, dim, dtype, data, name, cell_index);
    else
        return opp_decl_mesh_dat(set, dim, dtype, data, name);
}

//****************************************
opp_map opp_decl_map_txt(opp_set from, opp_set to, int dim, const char* file_name, 
                                char const *name)
{
    int* map_data = (int*)opp_load_from_file_core(file_name, from->size, dim, "int", sizeof(int));

    opp_map map = opp_decl_map(from, to, dim, map_data, name);

    opp_host_free(map_data);

    return map;
}

//****************************************
opp_dat opp_decl_dat_txt(opp_set set, int dim, opp_data_type dtype, const char* file_name, 
                                char const *name)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    char* dat_data = (char*)opp_load_from_file_core(file_name, set->size, dim, type.c_str(), size);

    opp_dat dat = opp_decl_dat_core(set, dim, type.c_str(), size, dat_data, name);

    opp_host_free(dat_data);

    return dat;
}

//****************************************
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, int dim, const char *typ, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat1");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, nullptr, acc);
}
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, int dim, const char *typ, opp_map p2c_map, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat2");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, p2c_map->p2c_dat, acc);
}

//****************************************
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, opp_map p2c_map, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat3");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, p2c_map->p2c_dat, acc);
}
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat4");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, nullptr, acc);
}
opp_arg opp_arg_dat(opp_dat dat, opp_map p2c_map, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat5");
    return opp_arg_dat_core(dat, -1, NULL, dat->dim, dat->type, p2c_map->p2c_dat, acc);
}
opp_arg opp_arg_dat(opp_dat dat, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat6");
    return opp_arg_dat_core(dat, -1, NULL, dat->dim, dat->type, nullptr, acc);
}

//****************************************
opp_arg opp_arg_dat(opp_map data_map, opp_access acc, bool offset)
{
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat7");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, -1, nullptr, 
                data_map->p2c_dat->dim, data_map->p2c_dat->type, nullptr, acc);
    return opp_arg_dat_core(data_map, -1, nullptr, nullptr, acc);
}
opp_arg opp_arg_dat(opp_map data_map, opp_map p2c_map, opp_access acc, bool offset)
{
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat8");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, -1, nullptr, 
            data_map->p2c_dat->dim, data_map->p2c_dat->type, p2c_map->p2c_dat, acc);
    return opp_arg_dat_core(data_map, -1, nullptr, p2c_map->p2c_dat, acc);
}
opp_arg opp_arg_dat(opp_map data_map, int idx, opp_map map, opp_access acc, bool offset)
{
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat9");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, idx, map, 
            data_map->p2c_dat->dim, data_map->p2c_dat->type, nullptr, acc);
    return opp_arg_dat_core(data_map, idx, map, nullptr, acc);
}
opp_arg opp_arg_dat(opp_map data_map, int idx, opp_map map, opp_map p2c_map, opp_access acc, bool offset)
{
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat10");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, idx, map, 
            data_map->p2c_dat->dim, data_map->p2c_dat->type, p2c_map->p2c_dat, acc);
    return opp_arg_dat_core(data_map, idx, map, p2c_map->p2c_dat, acc);
}


//****************************************
// template <class T> opp_arg opp_arg_gbl(T *data, int dim, char const *typ, opp_access acc);
opp_arg opp_arg_gbl(double *data, int dim, char const *typ, opp_access acc)
{
    return opp_arg_gbl_core(data, dim, typ, acc);
}
opp_arg opp_arg_gbl(int *data, int dim, char const *typ, opp_access acc)
{
    return opp_arg_gbl_core(data, dim, typ, acc);
}
opp_arg opp_arg_gbl(const bool *data, int dim, char const *typ, opp_access acc)
{
    return opp_arg_gbl_core(data, dim, typ, acc);
}

//****************************************
opp_set opp_decl_particle_set(int size, char const *name, opp_set cells_set)
{
    opp_set set = opp_decl_particle_set_core(size, name, cells_set);

    cutilSafeCall(hipMalloc((void**)&(set->particle_remove_count_d), sizeof(int)));
    cutilSafeCall(hipDeviceSynchronize());
  
    return set;
}
opp_set opp_decl_particle_set(char const *name, opp_set cells_set)
{
    return opp_decl_particle_set(0, name, cells_set);
}

//****************************************
opp_dat opp_decl_part_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name, 
                            bool cell_index)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    opp_dat dat = opp_decl_particle_dat_core(set, dim, type.c_str(), size, (char*)data, name, cell_index);

    opp_create_device_arrays(dat);

    opp_upload_dat(dat);

    return dat;
}

//****************************************
opp_dat opp_decl_particle_dat_txt(opp_set set, int dim, opp_data_type dtype, const char* file_name, 
                                        char const *name, bool cell_index)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    char* dat_data = (char*)opp_load_from_file_core(file_name, set->size, dim, type.c_str(), size);

    opp_dat dat = opp_decl_particle_dat_core(set, dim, type.c_str(), size, dat_data, name, cell_index);

    opp_host_free(dat_data);

    return dat;
}

//****************************************
void opp_increase_particle_count(opp_set part_set, const int num_particles_to_insert)
{ 
    opp_profiler->start("Inc_part_count");

    bool need_resizing = (part_set->set_capacity < (part_set->size + num_particles_to_insert)) ? true : false;

    if (OPP_DBG) 
        opp_printf("opp_increase_particle_count", "need_resizing %s", need_resizing ? "YES" : "NO");

    // TODO : We should be able to do a device to device copy instead of getting to host

    if (need_resizing) 
    {
        opp_profiler->start("Inc_part_count_DWN");
        opp_download_particle_set(part_set, true); 
        opp_profiler->end("Inc_part_count_DWN");
    }

    opp_profiler->start("Inc_part_count_INC");
    if (!opp_increase_particle_count_core(part_set, num_particles_to_insert))
    {
        opp_printf("opp_increase_particle_count", "Error at opp_increase_particle_count_core");
        opp_abort();
    }
    opp_profiler->end("Inc_part_count_INC");

    if (need_resizing)
    {
        opp_profiler->start("Inc_part_count_UPL");
        for (opp_dat& current_dat : *(part_set->particle_dats))
        {
            if (OPP_DBG) opp_printf("opp_increase_particle_count", "resizing dat [%s] set_capacity [%d]", 
                            current_dat->name, part_set->set_capacity);

            // TODO : We might be able to copy only the old data from device to device!

            opp_create_device_arrays(current_dat, true);

            opp_upload_dat(current_dat);

            current_dat->dirty_hd = Dirty::NotDirty;
        }   
        opp_profiler->end("Inc_part_count_UPL");     
    } 

    opp_profiler->end("Inc_part_count");
}

// opp_inc_part_count_with_distribution() is in opp_increase_part_count.cu

//****************************************
void opp_particle_sort(opp_set set)
{ 
    particle_sort_device(set, false);
}

//****************************************
void opp_print_dat_to_txtfile(opp_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{
    if (OPP_DBG) opp_printf("opp_print_dat_to_txtfile", "writing file [%s]", file_name_suffix);

    if (dat->dirty_hd == Dirty::Host) 
        opp_download_dat(dat);

    std::string prefix = std::string(file_name_prefix) + "_h";
    opp_print_dat_to_txtfile_core(dat, prefix.c_str(), file_name_suffix);
}

//****************************************
void opp_print_map_to_txtfile(opp_map map, const char *file_name_prefix, const char *file_name_suffix)
{
    if (OPP_DBG) opp_printf("opp_print_map_to_txtfile", "writing file [%s]", file_name_suffix);
    
    std::string prefix = std::string(file_name_prefix) + "_h";

    opp_print_map_to_txtfile_core(map, file_name_prefix, file_name_suffix);
}

//****************************************
// Set the complete dat to zero (upto array capacity)
void opp_reset_dat_impl(opp_dat dat, char* val, opp_reset reset)
{
    if (OPP_DBG) 
        opp_printf("opp_reset_dat_impl", "dat [%s] dim [%d] dat size [%d] set size [%d] set capacity [%d]", 
            dat->name, dat->dim, dat->size, dat->set->size, dat->set->set_capacity);

#ifdef USE_MPI
    int start = 0, end = dat->set->size;
    opp_get_start_end(dat->set, reset, start, end);

    size_t element_size = dat->size / dat->dim;

    for (int64_t i = 0; i < dat->dim; i++) 
    { 
        size_t data_d_offset = (i * dat->set->set_capacity + start) * element_size;

        if (OPP_DBG) 
            opp_printf("opp_reset_dat_impl", "dat %s dim %lld bytes_to_copy_per_dim %zu %p offset %zu", 
                dat->name, i, (end - start) * element_size, dat->data_d, data_d_offset);

        cutilSafeCall(hipMemset((dat->data_d + data_d_offset), 0, (end - start) * element_size));  
    }
#else
    cutilSafeCall(hipMemset((double*)(dat->data_d), 0, dat->size * dat->set->set_capacity));
#endif

    cutilSafeCall(hipDeviceSynchronize());

    dat->dirty_hd = Dirty::Host;
    dat->dirtybit = 1;
}


//*******************************************************************************
// This routine assumes that the host data structures are clean
void opp_partition(std::string lib_name, opp_set prime_set, opp_map prime_map, opp_dat data)
{
#ifdef USE_MPI
    // remove all negative mappings and copy the first mapping of the current element for all negative mappings
    opp_sanitize_all_maps();

    opp_partition_core(lib_name, prime_set, prime_map, data);

    opp_desanitize_all_maps();

    for (opp_dat dat : opp_dats) 
    {
        // no halos for particle dats
        if (dat->set->is_particle) 
            continue;

        opp_create_device_arrays(dat, true);
        opp_upload_dat(dat);
    }

    for (opp_map map : opp_maps) 
    {
        opp_upload_map(map, true);
    }
#endif
}

//*******************************************************************************
opp_dat opp_fetch_data(opp_dat dat) {

    if (dat->set->is_particle) {
        opp_printf("opp_fetch_data", "Error Cannot rearrange particle dats");
        opp_abort();
    }

    if (dat->dirty_hd == Dirty::Host) 
        opp_download_dat(dat);

#ifdef USE_MPI
    // rearrange data backe to original order in mpi
    return opp_mpi_get_data(dat);
#else
    return dat;
#endif
}


//*******************************************************************************
void opp_mpi_print_dat_to_txtfile(opp_dat dat, const char *file_name) 
{
    cutilSafeCall(hipDeviceSynchronize());
#ifdef USE_MPI
    const std::string prefixed_file_name = std::string("mpi_files/MPI_") + 
                                            std::to_string(OPP_comm_size) + std::string("_") + file_name;

    // rearrange data back to original order in mpi
    opp_dat temp = opp_fetch_data(dat);
    
    print_dat_to_txtfile_mpi(temp, prefixed_file_name.c_str());

    opp_host_free(temp->data);
    opp_host_free(temp->set);
    opp_host_free(temp);
#endif
}




//****************************************
void opp_dump_dat(opp_dat dat)
{
    cutilSafeCall(hipDeviceSynchronize());

    if (dat->dirty_hd == Dirty::Host) 
        opp_download_dat(dat);

    opp_dump_dat_core(dat);
}

//****************************************
// DEVICE->HOST | this invalidates what is in the HOST
void opp_download_dat(opp_dat dat) 
{ 
    cutilSafeCall(hipDeviceSynchronize());

    size_t set_size = dat->set->set_capacity;
    if (strstr(dat->type, ":soa") != NULL || (OPP_auto_soa && dat->dim > 1)) 
    {
        if (OPP_DBG) opp_printf("opp_download_dat", "GPU->CPU SOA | %s", dat->name);

        char *temp_data = (char *)opp_host_malloc(dat->size * set_size * sizeof(char));
        cutilSafeCall(hipMemcpy(temp_data, dat->data_d, set_size * dat->size, hipMemcpyDeviceToHost));
        
        int element_size = dat->size / dat->dim;
        for (int i = 0; i < dat->dim; i++) 
        {
            for (int j = 0; j < set_size; j++) 
            {
                for (int c = 0; c < element_size; c++) 
                {
                    dat->data[dat->size * j + element_size * i + c] = 
                        temp_data[element_size * i * set_size + element_size * j + c];
                }
            }
        }
        opp_host_free(temp_data);
    } 
    else 
    {
        if (OPP_DBG) opp_printf("opp_download_dat", "GPU->CPU NON-SOA| %s", dat->name);

        cutilSafeCall(hipMemcpy(dat->data, dat->data_d, set_size * dat->size, hipMemcpyDeviceToHost));
    }

    dat->dirty_hd = Dirty::NotDirty;
}

//****************************************
// HOST->DEVICE | this invalidates what is in the DEVICE
void opp_upload_dat(opp_dat dat)
{ 

    size_t set_capacity = dat->set->set_capacity;

    if (strstr(dat->type, ":soa") != NULL || (OPP_auto_soa && dat->dim > 1)) 
    {
        if (OPP_DBG) opp_printf("opp_upload_dat","CPU->GPU SOA | %s", dat->name);

        char *temp_data = (char *)opp_host_malloc(dat->size * set_capacity * sizeof(char));
        int element_size = dat->size / dat->dim;

        for (int i = 0; i < dat->dim; i++) 
        {
            for (int j = 0; j < set_capacity; j++) 
            {
                for (int c = 0; c < element_size; c++) 
                {
                    temp_data[element_size * i * set_capacity + element_size * j + c] = 
                        dat->data[dat->size * j + element_size * i + c];
                }
            }
        }

        opp_cpHostToDevice((void **)&(dat->data_d), (void **)&(temp_data), (dat->size * set_capacity), 
                                (dat->size * set_capacity), false);
        opp_host_free(temp_data);
    } 
    else 
    {
        if (OPP_DBG) opp_printf("opp_upload_dat", "CPU->GPU NON-SOA| %s", dat->name);
        opp_cpHostToDevice((void **)&(dat->data_d), (void **)&(dat->data), (dat->size * set_capacity), 
                                (dat->size * set_capacity), false);
    }

    dat->dirty_hd = Dirty::NotDirty;
}

//****************************************
// HOST->DEVICE | this invalidates what is in the DEVICE
void opp_upload_map(opp_map map, bool create_new) 
{
    if (OPP_DBG) opp_printf("opp_upload_map", "CPU->GPU | %s %s", map->name, create_new ? "NEW" : "COPY");
    int set_size = map->from->size + map->from->exec_size + map->from->nonexec_size;
    int *temp_map = (int *)opp_host_malloc(map->dim * set_size * sizeof(int));

    const int set_size_plus_exec = map->from->size + map->from->exec_size;
    
    for (int i = 0; i < map->dim; i++) 
    {
        for (int j = 0; j < set_size; j++) 
        {
            if (j >= set_size_plus_exec)
                temp_map[i * set_size + j] = -10;
            else
                temp_map[i * set_size + j] = map->map[map->dim * j + i];
        }
    }

    size_t copy_size = map->dim * set_size * sizeof(int);
    opp_cpHostToDevice((void **)&(map->map_d), (void **)&(temp_map), copy_size, copy_size, create_new);
    opp_host_free(temp_map);
}

//****************************************
// DEVICE -> HOST
void opp_download_particle_set(opp_set particles_set, bool force_download)
{

    if (OPP_DBG) opp_printf("opp_download_particle_set", "set [%s]", particles_set->name);

    cutilSafeCall(hipDeviceSynchronize());

    for (opp_dat& current_dat : *(particles_set->particle_dats))
    {
        if (current_dat->data_d == NULL)
        {
            if (OPP_DBG) opp_printf("opp_download_particle_set", "device pointer is NULL in dat [%s]", 
                            current_dat->name);
            continue;
        }
        if (current_dat->dirty_hd != Dirty::Host && !force_download)
        {
            if (OPP_DBG) opp_printf("opp_download_particle_set", "host is not dirty in dat [%s]", 
                            current_dat->name);
            continue;
        }
        opp_download_dat(current_dat);
    }  
}

//****************************************
// HOST->DEVICE
void opp_upload_particle_set(opp_set particles_set, bool realloc)
{ 

    if (OPP_DBG) opp_printf("opp_upload_particle_set", "set [%s]", particles_set->name);

    for (opp_dat& current_dat : *(particles_set->particle_dats))
    {
        if (realloc)
        {
            opp_create_device_arrays(current_dat, realloc);
        }

        opp_upload_dat(current_dat);
    }  
}

//*******************************************************************************

//****************************************
void opp_init_double_indirect_reductions_hip(int nargs, opp_arg *args)
{
#ifdef USE_MPI
    opp_init_double_indirect_reductions(nargs, args);
#endif
}

//****************************************
void opp_exchange_double_indirect_reductions_hip(int nargs, opp_arg *args)
{
#ifdef USE_MPI

    if (OPP_DBG) opp_printf("opp_exchange_double_indirect_reductions_hip", "ALL START");

    cutilSafeCall(hipDeviceSynchronize());

    for (int n = 0; n < nargs; n++) 
    {
        bool already_done = false;

        // check if the dat is mapped with double indirect mapping
        if (is_double_indirect_reduction(args[n]))
        {
            // Check if dat reduction was already done within these args
            for (int m = 0; m < n; m++) 
            {
                if (args[n].dat == args[m].dat)
                    already_done = true;
            }

            if (!already_done)
            {
                opp_download_dat(args[n].dat);
            }
        }
    }

    opp_exchange_double_indirect_reductions(nargs, args);

    if (OPP_DBG) opp_printf("opp_exchange_double_indirect_reductions_hip", "ALL END");
#endif
}

//****************************************
void opp_complete_double_indirect_reductions_hip(int nargs, opp_arg *args)
{
#ifdef USE_MPI
    if (OPP_DBG) opp_printf("opp_complete_double_indirect_reductions_hip", "ALL START");

    opp_complete_double_indirect_reductions(nargs, args);

    for (int n = 0; n < nargs; n++) 
    {
        bool already_done = false;

        // check if the dat is mapped with double indirect mapping
        if (is_double_indirect_reduction(args[n]))
        {
            // Check if dat reduction was already done within these args
            for (int m = 0; m < n; m++) 
            {
                if (args[n].dat == args[m].dat)
                    already_done = true;
            }

            if (!already_done)
            {
                opp_upload_dat(args[n].dat);
            }
        }
    }  

    if (OPP_DBG) opp_printf("opp_complete_double_indirect_reductions_hip", "ALL END");  
#endif
}

//****************************************

// **************************************** UTILITY FUNCTIONS ****************************************

//****************************************
void __hipSafeCall(hipError_t err, const char *file, const int line) 
{
    if (hipSuccess != err) 
    {
        // fprintf(stderr, "%s(%i) : cutilSafeCall() Runtime API error : %s.\n", file, line, 
        //     hipGetErrorString(err));
        std::string log = std::string(file) + "(" + std::to_string(line);
        log += std::string(") cutilSafeCall() Runtime API error : ") + hipGetErrorString(err);
        opp_abort(log.c_str());
    }
}

//****************************************
void __cutilCheckMsg(const char *errorMessage, const char *file, const int line) 
{
    hipError_t err = hipGetLastError();
    if (hipSuccess != err) 
    {
        fprintf(stderr, "%s(%i) : cutilCheckMsg() error : %s : %s.\n", file, line, errorMessage, 
            hipGetErrorString(err));
        opp_abort();
    }
}

//****************************************
void opp_cpHostToDevice(void **data_d, void **data_h, size_t copy_size, size_t alloc_size, bool create_new) 
{
    if (create_new)
    {
        if (*data_d != NULL) cutilSafeCall(hipFree(*data_d));
        cutilSafeCall(hipMalloc(data_d, alloc_size));
    }

    cutilSafeCall(hipMemcpy(*data_d, *data_h, copy_size, hipMemcpyHostToDevice));
    cutilSafeCall(hipDeviceSynchronize());
}

//****************************************
void opp_create_device_arrays(opp_dat dat, bool create_new)
{
    if (OPP_DBG) opp_printf("opp_create_device_arrays", "%s %s", dat->name, dat->type);

    char* temp_char_d = nullptr;

    if (strcmp(dat->type, "double") == 0)
    {
        if (dat->set->size > 0)
        {
            if (create_new && dat->thrust_real)
            {
                delete dat->thrust_real;
                delete dat->thrust_real_sort;
            }

            dat->thrust_real = new thrust::device_vector<double>(dat->set->set_capacity * dat->dim);
            dat->data_d = (char*)thrust::raw_pointer_cast(dat->thrust_real->data());

            dat->thrust_real_sort = new thrust::device_vector<double>(dat->set->set_capacity * dat->dim);
            temp_char_d = (char*)thrust::raw_pointer_cast(dat->thrust_real_sort->data());
        } 
    } 
    else if (strcmp(dat->type, "int") == 0 )
    {
        if (dat->set->size > 0)
        {
            if (create_new && dat->thrust_int)
            {
                delete dat->thrust_int;
                delete dat->thrust_int_sort;
            }

            dat->thrust_int = new thrust::device_vector<int>(dat->set->set_capacity * dat->dim);
            dat->data_d = (char*)thrust::raw_pointer_cast(dat->thrust_int->data());

            dat->thrust_int_sort = new thrust::device_vector<int>(dat->set->set_capacity * dat->dim);
            temp_char_d = (char*)thrust::raw_pointer_cast(dat->thrust_int_sort->data());
        } 
    }
    else
    {
        std::cerr << "opp_create_device_arrays DEVICE not implemented for type: " << dat->type << " dat name: " << 
            dat->name << std::endl;
        opp_abort();
    }

    if (OPP_DBG) opp_printf("opp_create_device_arrays", "Device array of dat [%s][%p][%p] Capacity [%d]", 
                        dat->name, dat->data_d, temp_char_d, dat->set->set_capacity * dat->dim);
}

//****************************************
void cutilDeviceInit(int argc, char **argv) 
{
    (void)argc;
    (void)argv;
    int deviceCount;

    cutilSafeCall(hipGetDeviceCount(&deviceCount));
    if (deviceCount == 0) 
    {
        opp_printf("cutilDeviceInit", "cutil error: no devices supporting DEVICE");
        opp_abort();
    }
  
    int int_rank = OPP_rank;
#ifdef USE_MPI
        opp::Comm comm(MPI_COMM_WORLD);
        int_rank = comm.rank_intra;
#endif

    // Test we have access to a device
    hipError_t err = hipSetDevice(int_rank % deviceCount);
    if (err == hipSuccess) 
    {
        float *test;
        if (hipMalloc((void **)&test, sizeof(float)) != hipSuccess) {
            OPP_hybrid_gpu = 0;
        }
        else {
            cutilSafeCall(hipFree(test));
            OPP_hybrid_gpu = 1;
        }
    }

    if (OPP_hybrid_gpu == 0) 
    {
        opp_printf("cutilDeviceInit", "Error... Init device Device Failed");
        opp_abort();
    }
}

//****************************************
void print_last_hip_error()
{
    printf("ANY hip ERRORS? %s\n", hipGetErrorString(hipGetLastError()));
}

//*******************************************************************************
void opp_colour_cartesian_mesh(const int ndim, std::vector<int> cell_counts, opp_dat cell_index, 
                            const opp_dat cell_colors, const int cell_ghosts) 
{
#ifdef USE_MPI  
    __opp_colour_cartesian_mesh(ndim, cell_counts, cell_index, cell_colors, cell_ghosts);
#endif
}



//*******************************************************************************
void* opp_host_malloc(size_t size)
{
    return malloc(size);
}

//*******************************************************************************
void* opp_host_realloc(void* ptr, size_t new_size)
{
    return realloc(ptr, new_size);
}

//*******************************************************************************
void opp_host_free(void* ptr)
{
    free(ptr);
}







//*******************************************************************************
void opp_reallocReductArrays(int reduct_bytes) 
{
    if (reduct_bytes > OPP_reduct_bytes) 
    {
        if (OPP_reduct_bytes > 0) 
        {
            free(OPP_reduct_h);
            cutilSafeCall(hipFree(OPP_reduct_d));
        }
        OPP_reduct_bytes = 4 * reduct_bytes; // 4 is arbitrary, more than needed
        OPP_reduct_h = (char *)malloc(OPP_reduct_bytes);
        cutilSafeCall(hipMalloc((void **)&OPP_reduct_d, OPP_reduct_bytes));
    }
}

void opp_mvReductArraysToDevice(int reduct_bytes) 
{
    cutilSafeCall(hipMemcpy(OPP_reduct_d, OPP_reduct_h, reduct_bytes, hipMemcpyHostToDevice));
    cutilSafeCall(hipDeviceSynchronize());
}

void opp_mvReductArraysToHost(int reduct_bytes) 
{
    cutilSafeCall(hipMemcpy(OPP_reduct_h, OPP_reduct_d, reduct_bytes, hipMemcpyDeviceToHost));
    cutilSafeCall(hipDeviceSynchronize());
}

void opp_reallocConstArrays(int consts_bytes) 
{
    if (consts_bytes > OPP_consts_bytes) 
    {
        if (OPP_consts_bytes > 0) 
        {
            free(OPP_consts_h);
            cutilSafeCall(hipFree(OPP_consts_d));
        }
        OPP_consts_bytes = 4 * consts_bytes; // 4 is arbitrary, more than needed
        OPP_consts_h = (char *)malloc(OPP_consts_bytes);
        cutilSafeCall(hipMalloc((void **)&OPP_consts_d, OPP_consts_bytes));
    }
}

void opp_mvConstArraysToDevice(int consts_bytes) 
{
    cutilSafeCall(hipMemcpy(OPP_consts_d, OPP_consts_h, consts_bytes, hipMemcpyHostToDevice));
    cutilSafeCall(hipDeviceSynchronize());
}

void opp_mvConstArraysToHost(int consts_bytes) 
{
    cutilSafeCall(hipMemcpy(OPP_consts_h, OPP_consts_d, consts_bytes, hipMemcpyDeviceToHost));
    cutilSafeCall(hipDeviceSynchronize());
}



// **************************************** REMOVED FUNCTIONS ****************************************

// // TODO: Try to do this in hip
// // make cell index of the particle to be removed as int_max,
// // sort all arrays
// // make sizes changed instead of resizing
// // Could do this inside device itself
// //****************************************
// void opp_finalize_particle_move_hip(opp_set set)
// {
//     cutilSafeCall(hipMemcpy(set->particle_statuses, set->particle_statuses_d, set->size * sizeof(int), 
//                                 hipMemcpyDeviceToHost));

//     set->particle_remove_count = 0;

//     for (int j = 0; j < set->size; j++)
//     {
//         if (set->particle_statuses[j] == OPP_NEED_REMOVE)
//         {
//             (set->particle_remove_count)++; // Could use a device variable and use atomicAdds inside device code
//         }
//     }   

//     if (OPP_DBG) printf("\topp_finalize_particle_move_hip set [%s] with particle_remove_count [%d]\n", 
//                     set->name, set->particle_remove_count);

//     if (set->particle_remove_count <= 0)
//     {
//         cutilSafeCall(hipFree(set->particle_statuses_d));
//         set->particle_statuses_d = NULL;
//         return;
//     }

//     opp_download_particle_set(set); // TODO : Find a better way

//     opp_finalize_particle_move_core(set);

//     opp_upload_particle_set(set, true /* realloc the device pointers */);

//     cutilSafeCall(hipFree(set->particle_statuses_d));
//     set->particle_statuses_d = NULL;

//     if (OPP_DBG) printf("\topp_finalize_particle_move_hip set [%s] with new size [%d]\n", 
//                     set->name, set->size);
// }

// //****************************************
// void opp_pack_marked_particles_to_move(opp_set set)
// {
// #ifdef USE_MPI

//     OPP_need_remove_flags.resize(OPP_need_remove_flags_size, 0);
// if (OPP_DBG) opp_printf("CHECK", "allocating OPP_need_remove_flags with size %d %zu", 
//                 OPP_need_remove_flags_size, OPP_need_remove_flags.size());

//     cutilSafeCall(hipMemcpy((char*)&(OPP_need_remove_flags[0]), OPP_need_remove_flags_d, 
//         sizeof(char) * OPP_need_remove_flags_size, hipMemcpyDeviceToHost));

// if (OPP_DBG) opp_printf("PPPPPPPPP", "SSSSSS");
// int p = 0;
// for (size_t particle_index = 0; particle_index < OPP_need_remove_flags.size(); particle_index++)
// {
//     if (OPP_need_remove_flags[particle_index] == 1)
//         p++;
// }
// if (OPP_DBG) opp_printf("PPPPPPPPP", 
//                 "FLAG All %d Move %d", OPP_need_remove_flags.size(), p);

//     int flagged = 0;

//     for (size_t particle_index = 0; particle_index < OPP_need_remove_flags.size(); particle_index++)
//     {
//         if (OPP_need_remove_flags[particle_index] == 1)
//         { flagged++;
//             std::map<int, opp_particle_comm_data>& set_part_com_data = opp_part_comm_neighbour_data[set];
            
//             int map0idx = OPP_mesh_relation_data[particle_index];

//             auto it = set_part_com_data.find(map0idx);
//             if (it == set_part_com_data.end())
//             {
//                 opp_printf("opp_pack_marked_particles_to_move", 
//                     "Error: cell %d cannot be found in opp_part_comm_neighbour_data map", map0idx);
//                 return; // unlikely, need exit(-1) to abort instead!
//             }

//             opp_particle_comm_data& comm_data = it->second;

//             opp_part_mark_move(set, particle_index, comm_data);

//             // This particle is already packed, hence needs to be removed from the current rank
//             OPP_mesh_relation_data[particle_index] = MAX_CELL_INDEX;
//         }
//     }

// if (OPP_DBG) opp_printf("opp_pack_marked_particles_to_move", 
//                 "FLAG All %d Move %d", OPP_need_remove_flags.size(), flagged);

//     // This is because hole filling/ sorting is done on the device arrays
//     opp_upload_dat(set->mesh_relation_dat);
// #endif
// }

// **************************************** ***************** ****************************************
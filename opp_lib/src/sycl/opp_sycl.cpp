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

#include <opp_sycl.h>
#include "opp_particle_mover.cpp"
#include "opp_increase_part_count.cpp"
#include "opp_sycl_utils.cpp"

sycl::queue *opp_queue = nullptr;

// arrays for global constants and reductions
int OPP_consts_bytes = 0, OPP_reduct_bytes = 0;
char *OPP_reduct_h = nullptr, *OPP_reduct_d = nullptr;
char *OPP_consts_h = nullptr, *OPP_consts_d = nullptr;

std::vector<char*> opp_consts;

char opp_move_status_flag = OPP_MOVE_DONE;
bool opp_move_hop_iter_one_flag = true;
OPP_INT* opp_p2c = nullptr;
OPP_INT* opp_c2c = nullptr;

opp_dh_indices dh_indices_h;
opp_dh_indices dh_indices_d;

//****************************************
void opp_init(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << std::endl;
        opp_abort();
    }

    if (OPP_DBG)
        opp_printf("opp_init", "Starting init");

#if defined(USE_PETSC)
    PetscInitialize(&argc, &argv, PETSC_NULLPTR, "opp::PetscSYCL");
#elif defined(USE_MPI)
        MPI_Init(&argc, &argv);
#endif

    opp_init_core(argc, argv);

#ifdef USE_MPI
    int size = -1, rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int opp_partitions = opp_params->get<OPP_INT>("opp_partitions", false); 
    if (opp_partitions != INT_MAX) {
        if (size % opp_partitions != 0) {
            opp_abort("The MPI ranks is not divisible by partition count");
        }

        const int gap = size / opp_partitions;
        int color = (rank % gap == 0) ? 1 : MPI_UNDEFINED;
        MPI_Comm_split(MPI_COMM_WORLD, color, rank, &OPP_MPI_WORLD);

        if (color != MPI_UNDEFINED) {
            MPI_CHECK(MPI_Comm_rank(OPP_MPI_WORLD, &OPP_rank))
            MPI_CHECK(MPI_Comm_size(OPP_MPI_WORLD, &OPP_comm_size))
        } 
        else {
            // Handle case where the rank is not part of the new communicator
            OPP_rank = -1 * rank;
            OPP_comm_size = 0;
        }
    }
    else {
        OPP_MPI_WORLD = MPI_COMM_WORLD;
        OPP_rank = rank;
        OPP_comm_size = size;
    }
#endif

    if (OPP_rank == OPP_ROOT) {
        std::string log = "Running on SYCL";
#ifdef USE_MPI
        log += "+MPI with " + std::to_string(OPP_comm_size) + " GPU ranks, but with " + std::to_string(size) + " total ranks";
#endif        
        opp_printf("OP-PIC", "%s", log.c_str());
        opp_printf("OP-PIC", "---------------------------------------------");
    }
    // opp_printf("OP-PIC", "old rank %d new rank %d", rank, OPP_rank);

#ifdef USE_MPI
    MPI_Barrier(OPP_MPI_WORLD);
#endif  

    OPP_RETURN_IF_INVALID_PROCESS;

    opp_params->write(std::cout);
    opp_device_init(argc, argv);

    OPP_auto_soa = 1; // TODO : Make this configurable with args

    const int threads_per_block = opp_params->get<OPP_INT>("opp_threads_per_block");   
    if (threads_per_block > 0 && threads_per_block < INT_MAX)
        OPP_gpu_threads_per_block = threads_per_block;

    const int gpu_direct = opp_params->get<OPP_INT>("opp_gpu_direct");   
    if (gpu_direct > 0 && gpu_direct < INT_MAX)
        OPP_gpu_direct = gpu_direct;

    const int deviceId = dpct::dev_mgr::instance().current_device_id();
    dpct::device_info prop;
    dpct::get_device_info(prop, dpct::dev_mgr::instance().get_device(deviceId));

    // DPCT1019:13: local_mem_size in SYCL is not a complete equivalent of sharedMemPerBlock in CUDA.
    // OPP_gpu_shared_mem_per_block = prop.get_local_mem_size();

    char hostname[256];
    if (gethostname(hostname, sizeof(hostname)) != 0) {
        opp_printf("opp_init", "Failed to get hostname of MPI rank %d", OPP_rank);
        opp_abort();
    }

    opp_printf("opp_init",  "Device: %d [%s] on Host [%s] threads=%d Shared_Mem=%lubytes GPU_Direct=%d",
        deviceId, prop.get_name(), hostname, OPP_gpu_threads_per_block, prop.get_local_mem_size(), OPP_gpu_direct);
}

//****************************************
void opp_exit()
{
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif  

    if (OPP_IS_VALID_PROCESS) {

        globalMover.reset();
        cellMapper.reset();
        boundingBox.reset();
        comm.reset();

        opp_mem::host_free(OPP_reduct_h);
        opp_mem::host_free(OPP_consts_h);

#ifdef USE_MPI 
        opp_halo_destroy(); // free memory allocated to halos and mpi_buffers 
        opp_partition_destroy(); // free memory used for holding partition information
        opp_part_comm_destroy(); // free memory allocated for particle communication
#endif

        opp_device_exit();
    }

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
void opp_device_init(int argc, char **argv) {

    char temp[64];
    bool run_on_gpu = true;
    for (int n = 1; n < argc; n++) {
        char *pch = strstr(argv[n], "-cpu");
        if (pch != NULL)
            run_on_gpu = false;
    }

    if (OPP_DBG) opp_printf("opp_device_init", "run_on_gpu=%s", run_on_gpu ? "TRUE" : "FALSE");

    try {
//         std::vector<sycl::device> selected_devices;
//         std::vector<sycl::device> all_devices = sycl::device::get_devices();   
//         if (all_devices.size() == 0) {
//             opp_printf("opp_device_init", "Error: no supporting devices found");
//             opp_abort();
//         }
//         for (const auto& dev : all_devices) {
//             if (OPP_DBG) {
//                 std::string device_type;  
//                 if (dev.is_cpu()) device_type = "CPU";
//                 else if (dev.is_gpu()) device_type = "GPU";
//                 else if (dev.is_accelerator()) device_type = "Accelerator";
//                 else device_type = "Unknown";
//                 std::cout << "opp_device_init[" << OPP_rank << "] - Detected device [ " 
//                     << dev.get_info<sycl::info::device::name>() << " (" 
//                     << device_type << ") ]" << std::endl;
//             }
//             // Should we include CPU and GPU loops as well?
//             if ((dev.is_cpu() && !run_on_gpu) || (dev.is_gpu() && run_on_gpu)) {
//                 selected_devices.push_back(dev);
//             }
//         }
//         if (selected_devices.size() == 0) {
//             opp_printf("opp_device_init", "Requested %s but none found", run_on_gpu ? "GPUs" : "CPUs");
//             opp_abort();
//         }
// #ifdef USE_MPI
//         opp::Comm comm(OPP_MPI_WORLD);
//         const int int_rank = comm.rank_intra;
// #else
//         const int int_rank = OPP_rank;
// #endif
//         // Select a GPU based on the MPI rank
//         sycl::device selected_device = selected_devices[int_rank % selected_devices.size()];
//         opp_queue = new sycl::queue(selected_device);

        // Instead of above commented code, use a wrapper to set device selectors!
        // Run using, mpirun -np 8 ./wrap_sycl.sh bin/sycl_mpi congigs/xx.param
        opp_queue = new sycl::queue(sycl::default_selector_v);
        
        // Test we have access to a device
        float *test = opp_mem::dev_malloc<float>(1);
        opp_mem::dev_free(test);
    }
    catch (sycl::exception const &exc) {
        
        std::cerr << exc.what() << "Exception caught at file:" << __FILE__
                    << ", line:" << __LINE__ << std::endl;
        opp_abort("opp_device_init - Error... Init device Device Failed");
    }
}

//****************************************
void opp_device_exit() 
{
    for (auto& set : opp_sets) {
        if (set->is_particle) 
            opp_mem::dev_free(set->particle_remove_count_d);
    }

    for (auto& map : opp_maps) {
        opp_mem::dev_free(map->map_d);
    }

    for (auto& dat : opp_dats) {
        opp_mem::dev_free(dat->data_d);
        opp_mem::dev_free(dat->data_swap_d);
    }

    opp_mem::dev_free(opp_saved_mesh_relation_d);

    opp_mem::dev_free(OPP_move_particle_indices_d);
    opp_mem::dev_free(OPP_move_cell_indices_d);
    opp_mem::dev_free(OPP_remove_particle_indices_d);
    opp_mem::dev_free(OPP_move_count_d);

    opp_mem::dev_free(hf_from_indices_dp);
    opp_mem::dev_free(ps_from_indices_dp);

    opp_mem::dev_free(OPP_reduct_d);
    opp_mem::dev_free(OPP_consts_d);

    for (auto& a : opp_consts)
        opp_mem::dev_free(a);

    // below are for GPU direct particle communication
    for (auto it = particle_indices_hv.begin(); it != particle_indices_hv.end(); it++) 
        it->second.clear();
    for (auto it = cell_indices_hv.begin(); it != cell_indices_hv.end(); it++) 
        it->second.clear();
    for (auto it = particle_indices_dv.begin(); it != particle_indices_dv.end(); it++) {
        it->second.clear(); it->second.shrink_to_fit();
    }
    for (auto it = send_data.begin(); it != send_data.end(); it++) {
        it->second.clear(); it->second.shrink_to_fit();
    }
    for (auto it = recv_data.begin(); it != recv_data.end(); it++) {
        it->second.clear(); it->second.shrink_to_fit();
    }
}

//********************************************************************************
opp_set opp_decl_set(int size, char const *name)
{
    return opp_decl_set_core(size, name);
}

//****************************************
opp_set opp_decl_particle_set(int size, char const *name, opp_set cells_set)
{
    opp_set set = opp_decl_particle_set_core(size, name, cells_set);
    if (OPP_IS_VALID_PROCESS)
        set->particle_remove_count_d = opp_mem::dev_malloc<OPP_INT>(1);
    return set;
}
opp_set opp_decl_particle_set(char const *name, opp_set cells_set)
{
    return opp_decl_particle_set(0, name, cells_set);
}

//****************************************
opp_map opp_decl_map(opp_set from, opp_set to, int dim, int *imap, char const *name)
{ OPP_RETURN_NULL_IF_INVALID_PROCESS;
    opp_map map = opp_decl_map_core(from, to, dim, imap, name);

    if (from->is_particle) {
        opp_create_dat_device_arrays(map->p2c_dat);
        opp_upload_dat(map->p2c_dat);
    }
    else {
        opp_upload_map(map, true);
    }

    return map;
}

//****************************************
opp_dat opp_decl_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name)
{ OPP_RETURN_NULL_IF_INVALID_PROCESS;
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    opp_dat dat = opp_decl_dat_core(set, dim, type.c_str(), size, (char*)data, name);

    opp_create_dat_device_arrays(dat);
    opp_upload_dat(dat);

    return dat;
}

//********************************************************************************
opp_map opp_decl_map_txt(opp_set from, opp_set to, int dim, const char* file_name, char const *name)
{ OPP_RETURN_NULL_IF_INVALID_PROCESS;
    int* map_data = (int*)opp_load_from_file_core(file_name, from->size, dim, "int", sizeof(int));
    opp_map map = opp_decl_map(from, to, dim, map_data, name);
    opp_mem::host_free(map_data);

    return map;
}

//****************************************
opp_dat opp_decl_dat_txt(opp_set set, int dim, opp_data_type dtype, const char* file_name, char const *name)
{ OPP_RETURN_NULL_IF_INVALID_PROCESS;
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    char* dat_data = (char*)opp_load_from_file_core(file_name, set->size, dim, type.c_str(), size);

    opp_dat dat = opp_decl_dat_core(set, dim, type.c_str(), size, dat_data, name);

    opp_mem::host_free(dat_data);

    return dat;
}

//********************************************************************************
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, int dim, const char *typ, opp_access acc, 
                        bool offset)
{ OPP_RETURN_ARG_IF_INVALID_PROCESS;
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat1");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, nullptr, acc);
}
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, int dim, const char *typ, opp_map p2c_map, 
                        opp_access acc, bool offset)
{ OPP_RETURN_ARG_IF_INVALID_PROCESS;
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat2");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, p2c_map->p2c_dat, acc);
}

//****************************************
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, opp_map p2c_map, opp_access acc, bool offset)
{ OPP_RETURN_ARG_IF_INVALID_PROCESS;
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat3");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, p2c_map->p2c_dat, acc);
}
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, opp_access acc, bool offset)
{ OPP_RETURN_ARG_IF_INVALID_PROCESS;
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat4");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, nullptr, acc);
}
opp_arg opp_arg_dat(opp_dat dat, opp_map p2c_map, opp_access acc, bool offset)
{ OPP_RETURN_ARG_IF_INVALID_PROCESS;
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat5");
    return opp_arg_dat_core(dat, -1, NULL, dat->dim, dat->type, p2c_map->p2c_dat, acc);
}
opp_arg opp_arg_dat(opp_dat dat, opp_access acc, bool offset)
{ OPP_RETURN_ARG_IF_INVALID_PROCESS;
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat6");
    return opp_arg_dat_core(dat, -1, NULL, dat->dim, dat->type, nullptr, acc);
}

//****************************************
opp_arg opp_arg_dat(opp_map data_map, opp_access acc, bool offset)
{ OPP_RETURN_ARG_IF_INVALID_PROCESS;
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat7");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, -1, nullptr, 
                data_map->p2c_dat->dim, data_map->p2c_dat->type, nullptr, acc);
    return opp_arg_dat_core(data_map, -1, nullptr, nullptr, acc);
}
opp_arg opp_arg_dat(opp_map data_map, opp_map p2c_map, opp_access acc, bool offset)
{ OPP_RETURN_ARG_IF_INVALID_PROCESS;
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat8");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, -1, nullptr, 
            data_map->p2c_dat->dim, data_map->p2c_dat->type, p2c_map->p2c_dat, acc);
    return opp_arg_dat_core(data_map, -1, nullptr, p2c_map->p2c_dat, acc);
}
opp_arg opp_arg_dat(opp_map data_map, int idx, opp_map map, opp_access acc, bool offset)
{ OPP_RETURN_ARG_IF_INVALID_PROCESS;
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat9");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, idx, map, 
            data_map->p2c_dat->dim, data_map->p2c_dat->type, nullptr, acc);
    return opp_arg_dat_core(data_map, idx, map, nullptr, acc);
}
opp_arg opp_arg_dat(opp_map data_map, int idx, opp_map map, opp_map p2c_map, opp_access acc, bool offset)
{ OPP_RETURN_ARG_IF_INVALID_PROCESS;
    if (data_map == nullptr) opp_abort("dat is NULL at opp_arg_dat10");
    if (data_map->from->is_particle)
        return opp_arg_dat_core(data_map->p2c_dat, idx, map, 
            data_map->p2c_dat->dim, data_map->p2c_dat->type, p2c_map->p2c_dat, acc);
    return opp_arg_dat_core(data_map, idx, map, p2c_map->p2c_dat, acc);
}

//********************************************************************************
void opp_print_dat_to_txtfile(opp_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{ OPP_RETURN_IF_INVALID_PROCESS;
    if (OPP_DBG) opp_printf("opp_print_dat_to_txtfile", "writing file [%s %s] data %p data_d %p", 
                    file_name_prefix, file_name_suffix, dat->data, dat->data_d);

    if (dat->dirty_hd == Dirty::Host) 
        opp_download_dat(dat);

    std::string prefix = std::string(file_name_prefix) + "_sycl";
    opp_print_dat_to_txtfile_core(dat, prefix.c_str(), file_name_suffix);
}

//****************************************
void opp_print_map_to_txtfile(opp_map map, const char *file_name_prefix, const char *file_name_suffix)
{ OPP_RETURN_IF_INVALID_PROCESS;
    if (OPP_DBG) opp_printf("opp_print_map_to_txtfile", "writing file [%s]", file_name_suffix);
    
    std::string prefix = std::string(file_name_prefix) + "_sycl" + std::to_string(OPP_rank) + "_";

    opp_print_map_to_txtfile_core(map, prefix.c_str(), file_name_suffix);
}

//****************************************
opp_dat opp_fetch_data(opp_dat dat) 
{ OPP_RETURN_NULL_IF_INVALID_PROCESS;
    if (dat->set->is_particle) {
        opp_printf("opp_fetch_data", "Error Cannot rearrange particle dats");
        opp_abort();
    }

    if (dat->dirty_hd == Dirty::Host) 
        opp_download_dat(dat);

#ifdef USE_MPI
    return opp_mpi_get_data(dat); // rearrange data back to original order in mpi
#else
    return dat;
#endif
}

//****************************************
void opp_mpi_print_dat_to_txtfile(opp_dat dat, const char *file_name) 
{ OPP_RETURN_IF_INVALID_PROCESS;
    OPP_DEVICE_SYNCHRONIZE();

#ifdef USE_MPI
    const std::string prefixed_file_name = std::string("mpi_files/SYCL_") + 
                std::to_string(OPP_comm_size) + std::string("_iter") + 
                std::to_string(OPP_main_loop_iter) + std::string("_") + file_name;

    opp_dat temp = opp_fetch_data(dat); // rearrange data back to original order in mpi
    
    print_dat_to_txtfile_mpi(temp, prefixed_file_name.c_str());

    opp_mem::host_free(temp->data);
    opp_mem::host_free(temp->set);
    opp_mem::host_free(temp);
#endif
}

//****************************************
void opp_dump_dat(opp_dat dat)
{ OPP_RETURN_IF_INVALID_PROCESS;
    OPP_DEVICE_SYNCHRONIZE();

    if (dat->dirty_hd == Dirty::Host) opp_download_dat(dat);

    opp_dump_dat_core(dat);
}

//********************************************************************************
// Set the complete dat to zero (upto array capacity)
void opp_reset_dat_impl(opp_dat dat, char* val, opp_reset reset)
{ OPP_RETURN_IF_INVALID_PROCESS;
    if (OPP_DBG) 
        opp_printf("opp_reset_dat_impl", "dat [%s] dim [%d] dat size [%d] set size [%d] set capacity [%d]", 
            dat->name, dat->dim, dat->size, dat->set->size, dat->set->set_capacity);

#ifdef USE_MPI
    int start = 0, end = dat->set->size;
    opp_get_start_end(dat->set, reset, start, end);

    const size_t element_size = dat->size / dat->dim;

    for (int64_t i = 0; i < dat->dim; i++) { 
        size_t data_d_offset = (i * dat->set->set_capacity + start) * element_size;

        if (OPP_DBG) 
            opp_printf("opp_reset_dat_impl", "dat %s dim %lld bytes_to_copy_per_dim %zu %p offset %zu", 
                dat->name, i, (end - start) * element_size, dat->data_d, data_d_offset);

        opp_mem::dev_memset<char>(dat->data_d + data_d_offset, (end - start) * element_size, 0);
    }
#else
    opp_mem::dev_memset<char>(dat->data_d, dat->size * dat->set->set_capacity, 0);
#endif

    dat->dirty_hd = Dirty::Host;
    dat->dirtybit = 1;
}

//*******************************************************************************
// This routine assumes that the host data structures are clean
void opp_partition(std::string lib_name, opp_set prime_set, opp_map prime_map, opp_dat data)
{ OPP_RETURN_IF_INVALID_PROCESS;
#ifdef USE_MPI
    // remove all negative mappings and copy the first mapping of the current element for all negative mappings
    opp_sanitize_all_maps();

    opp_partition_core(lib_name, prime_set, prime_map, data);

    opp_desanitize_all_maps();

    for (opp_dat dat : opp_dats) {

        if (dat->set->is_particle) 
            continue; // no halos for particle dats

        opp_create_dat_device_arrays(dat, true);
        opp_upload_dat(dat);
    }

    for (opp_map map : opp_maps) 
        opp_upload_map(map, true);
#endif
}


//*******************************************************************************
void opp_init_double_indirect_reductions_device(int nargs, opp_arg *args)
{
#ifdef USE_MPI
    opp_init_double_indirect_reductions(nargs, args);
#endif
}

//****************************************
void opp_exchange_double_indirect_reductions_device(int nargs, opp_arg *args)
{
#ifdef USE_MPI

    if (OPP_DBG) opp_printf("opp_exchange_double_indirect_reductions_device", "ALL START");

    OPP_DEVICE_SYNCHRONIZE();

    for (int n = 0; n < nargs; n++) {
        bool already_done = false;

        // check if the dat is mapped with double indirect mapping
        if (is_double_indirect_reduction(args[n])) {
            
            // Check if dat reduction was already done within these args
            for (int m = 0; m < n; m++) {
                if (args[n].dat == args[m].dat)
                    already_done = true;
            }

            if (!already_done) {
                opp_download_dat(args[n].dat);
            }
        }
    }

    opp_exchange_double_indirect_reductions(nargs, args);

    if (OPP_DBG) opp_printf("opp_exchange_double_indirect_reductions_device", "ALL END");
#endif
}

//****************************************
void opp_complete_double_indirect_reductions_device(int nargs, opp_arg *args)
{
#ifdef USE_MPI
    if (OPP_DBG) opp_printf("opp_complete_double_indirect_reductions_device", "ALL START");

    opp_complete_double_indirect_reductions(nargs, args);

    for (int n = 0; n < nargs; n++) {
        bool already_done = false;

        // check if the dat is mapped with double indirect mapping
        if (is_double_indirect_reduction(args[n])) {
            
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

    if (OPP_DBG) opp_printf("opp_complete_double_indirect_reductions_device", "ALL END");  
#endif
}

//*******************************************************************************
void opp_reallocReductArrays(int reduct_bytes) 
{
    if (reduct_bytes > OPP_reduct_bytes) {
        if (OPP_reduct_bytes > 0) {
            opp_mem::host_free(OPP_reduct_h);
            opp_mem::dev_free(OPP_reduct_d);
        }
        OPP_reduct_bytes = 4 * reduct_bytes; // 4 is arbitrary, more than needed
        OPP_reduct_h = opp_mem::host_malloc<char>(OPP_reduct_bytes);
        OPP_reduct_d = opp_mem::dev_malloc<char>(OPP_reduct_bytes); // use opp_mem::dev_realloc instead   
    }
}

//****************************************
void opp_mvReductArraysToDevice(int reduct_bytes) 
{
    opp_mem::copy_host_to_dev<char>(OPP_reduct_d, OPP_reduct_h, reduct_bytes);
}

//****************************************
void opp_mvReductArraysToHost(int reduct_bytes) 
{
    OPP_DEVICE_SYNCHRONIZE();
    opp_mem::copy_dev_to_host<char>(OPP_reduct_h, OPP_reduct_d, reduct_bytes);
}

//*******************************************************************************
void opp_reallocConstArrays(int consts_bytes) 
{
    if (consts_bytes > OPP_consts_bytes) {
        if (OPP_consts_bytes > 0) {
            opp_mem::host_free(OPP_consts_h);
            opp_mem::dev_free(OPP_consts_d);
        }
        OPP_consts_bytes = 4 * consts_bytes; // 4 is arbitrary, more than needed
        OPP_consts_h = opp_mem::host_malloc<char>(OPP_consts_bytes);
        OPP_consts_d = opp_mem::dev_malloc<char>(OPP_consts_bytes); // use opp_mem::dev_realloc instead
    }
}

//****************************************
void opp_mvConstArraysToDevice(int consts_bytes) 
{
    opp_mem::copy_host_to_dev<char>(OPP_consts_d, OPP_consts_h, consts_bytes);
}

//****************************************
void opp_mvConstArraysToHost(int consts_bytes) 
{
    OPP_DEVICE_SYNCHRONIZE();
    opp_mem::copy_dev_to_host<char>(OPP_consts_h, OPP_consts_d, consts_bytes);
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
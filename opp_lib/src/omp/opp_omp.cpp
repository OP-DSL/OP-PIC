
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

#include <opp_omp.h>
#include "opp_increase_part_count.cpp"
#include "opp_particle_mover.cpp"

int opp_nthreads = 1;
std::vector<int> part_remove_count_per_thr;
std::vector<std::vector<int>> move_part_indices_per_thr;
std::vector<std::vector<opp_tmp_gbl_move_info>> gbl_move_indices_per_thr;

//****************************************
void opp_init(int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << std::endl;
        exit(-1);
    }

#ifdef USE_PETSC
    PetscInitialize(&argc, &argv, PETSC_NULLPTR, "opp::PetscOMP");
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

    opp_nthreads = omp_get_max_threads();

    if (OPP_rank == OPP_ROOT) {
        std::string log = "Running on OMP";
#ifdef USE_MPI
        log += "+MPI with " + std::to_string(OPP_comm_size) + " ranks";
#endif    
        log += " and " + std::to_string(opp_nthreads) + " threads per rank";
        opp_printf("OP-PIC", "%s", log.c_str());
        opp_printf("OP-PIC", "---------------------------------------------");
    }

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif  

    opp_init_core(argc, argv);
    opp_params->write(std::cout);
}

//****************************************
void opp_exit()
{
    if (OPP_DBG) opp_printf("opp_exit", "");
    
    globalMover.reset();
    cellMapper.reset();
    boundingBox.reset();
    comm.reset();

#ifdef USE_MPI 
        opp_halo_destroy(); // free memory allocated to halos and mpi_buffers 
        opp_partition_destroy(); // free memory used for holding partition information
        opp_part_comm_destroy(); // free memory allocated for particle communication
#endif

    for (auto& dat : opp_dats) {
        for (size_t thr = 1; thr < dat->thread_data->size(); thr++) { 
            opp_host_free(dat->thread_data->at(thr)); 
        }
    }
    opp_exit_core();

#ifdef USE_PETSC
    PetscFinalize();
#elif defined(USE_MPI)
    MPI_Finalize();
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
opp_set opp_decl_set(int size, char const *name)
{
    return opp_decl_set_core(size, name);
}

//****************************************
opp_map opp_decl_map(opp_set from, opp_set to, int dim, int *imap, char const *name)
{
    return opp_decl_map_core(from, to, dim, imap, name);
}

//****************************************
opp_dat opp_decl_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    return opp_decl_dat_core(set, dim, type.c_str(), size, (char*)data, name);
}

//****************************************
opp_map opp_decl_map_txt(opp_set from, opp_set to, int dim, const char* file_name, char const *name)
{
    int* map_data = (int*)opp_load_from_file_core(file_name, from->size, dim, "int", sizeof(int));

    opp_map map = opp_decl_map(from, to, dim, map_data, name);

    opp_host_free(map_data);

    return map;
}

//****************************************
opp_dat opp_decl_dat_txt(opp_set set, int dim, opp_data_type dtype, const char* file_name, char const *name)
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
opp_set opp_decl_particle_set(int size, char const *name, opp_set cells_set)
{
    return opp_decl_particle_set_core(size, name, cells_set);
}
opp_set opp_decl_particle_set(char const *name, opp_set cells_set)
{
    return opp_decl_particle_set_core(name, cells_set);
}

//****************************************
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, int dim, const char *typ, opp_access acc, bool offset)
{
    if (dat == nullptr) opp_abort("dat is NULL at opp_arg_dat1");
    return opp_arg_dat_core(dat, idx, map, dat->dim, dat->type, nullptr, acc);
}
opp_arg opp_arg_dat(opp_dat dat, int idx, opp_map map, int dim, const char *typ, opp_map p2c_map, 
                        opp_access acc, bool offset)
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
void opp_print_dat_to_txtfile(opp_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{
    std::string prefix = std::string(file_name_prefix) + "_omp";
    opp_print_dat_to_txtfile_core(dat, prefix.c_str(), file_name_suffix);
}

//****************************************
void opp_print_map_to_txtfile(opp_map map, const char *file_name_prefix, const char *file_name_suffix)
{
    std::string prefix = std::string(file_name_prefix) + "_omp";
    opp_print_map_to_txtfile_core(map, prefix.c_str(), file_name_suffix);
}

//****************************************
void opp_dump_dat(opp_dat dat)
{
    opp_dump_dat_core(dat);
}

//****************************************
void opp_reset_dat_impl(opp_dat dat, char* val, opp_reset reset)
{
    if (OPP_DBG) 
        opp_printf("opp_reset_dat_impl", "dat [%s] dim [%d] dat size [%d] set size [%d] set capacity [%d]", 
            dat->name, dat->dim, dat->size, dat->set->size, dat->set->set_capacity);

    int start = 0, end = dat->set->size;

#ifdef USE_MPI
    opp_get_start_end(dat->set, reset, start, end);
#endif

    // TODO : reset pragma omp
    for (size_t i = (size_t)start; i < (size_t)end; i++) {
        memcpy(dat->data + i * dat->size, val, dat->size);
    }
}

#ifndef USE_MPI // if USE_MPI is defined, the functions in opp_mpi_particle_comm.cpp and opp_mpi_halo.cpp is used
//****************************************
int opp_mpi_halo_exchanges_grouped(opp_set set, int nargs, opp_arg *args, DeviceType device) 
{
    return set->size;
}

//****************************************
int opp_mpi_halo_exchanges(opp_set set, int nargs, opp_arg *args) 
{
    return set->size;
}

//****************************************
void opp_mpi_force_halo_update_if_dirty(opp_set set, std::vector<opp_dat> dats, DeviceType device) 
{
    // Nothing to execute here
}

//****************************************
void opp_mpi_halo_wait_all(int nargs, opp_arg *args)
{
    // Nothing to execute here
}
#endif // if USE_MPI is defined, the functions in opp_mpi_particle_comm.cpp and opp_mpi_halo.cpp is used

//*******************************************************************************
void opp_partition(std::string lib_name, opp_set prime_set, opp_map prime_map, opp_dat data)
{
#ifdef USE_MPI
    opp_profiler->start("opp_partition");

    // remove all negative mappings and copy the first mapping of the current element for all negative mappings
    opp_sanitize_all_maps();

    opp_partition_core(lib_name, prime_set, prime_map, data);

    opp_desanitize_all_maps();
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    opp_profiler->end("opp_partition");
#endif
}

//*******************************************************************************
void opp_mpi_print_dat_to_txtfile(opp_dat dat, const char *file_name) 
{
#ifdef USE_MPI
    const std::string prefixed_file_name = std::string("mpi_files/MPI_") + 
                std::to_string(OPP_comm_size) + std::string("_iter") + 
                std::to_string(OPP_main_loop_iter) + std::string("_") + file_name;

    if (dat->set->is_particle) {
        opp_printf("opp_mpi_print_dat_to_txtfile", "Error Cannot rearrange particle dats");
        opp_abort();
    }

    // rearrange data back to original order in mpi
    opp_dat temp = opp_mpi_get_data(dat);
    
    print_dat_to_txtfile_mpi(temp, prefixed_file_name.c_str());

    opp_host_free(temp->data);
    opp_host_free(temp->set);
    opp_host_free(temp);
#endif
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
opp_dat opp_fetch_data(opp_dat dat) 
{
#ifdef USE_MPI
    if (dat->set->is_particle) {
        opp_printf("opp_fetch_data", "Error Cannot rearrange particle dats");
        opp_abort();
    }

    // rearrange data backe to original order in mpi
    return opp_mpi_get_data(dat);
#else
    return nullptr;
#endif
}

//*******************************************************************************
// Below API only for GPU backends
void opp_upload_dat(opp_dat dat) {}
void opp_download_dat(opp_dat dat) {}
void opp_download_particle_set(opp_set particles_set, bool force_download) {}
void opp_upload_particle_set(opp_set particles_set, bool realloc) {}

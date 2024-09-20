
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

#include <opp_mpi.h>

opp_move_var move_var;
void opp_part_pack(opp_set set);
void opp_part_unpack(opp_set set);

char opp_move_status_flag = OPPX_MOVE_DONE;
bool opp_move_hop_iter_one_flag = true;
OPP_INT* opp_p2c = nullptr;
OPP_INT* opp_c2c = nullptr;

//*******************************************************************************
void opp_init(int argc, char **argv) 
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config_file> ..." << std::endl;
        exit(-1);
    }

#ifdef USE_PETSC
    PetscInitialize(&argc, &argv, PETSC_NULLPTR, "opp::PetscMPI");
#else
    MPI_Init(&argc, &argv);
#endif

    OPP_MPI_WORLD = MPI_COMM_WORLD;

    MPI_Comm_rank(MPI_COMM_WORLD, &OPP_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &OPP_comm_size);

    if (OPP_rank == OPP_ROOT)
    {
        std::string log = "Running on MPI with " + std::to_string(OPP_comm_size) + " ranks";
        opp_printf("OP-PIC", "%s", log.c_str());
        opp_printf("OP-PIC", "---------------------------------------------");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    opp_init_core(argc, argv);
    opp_params->write(std::cout);

    opp_profiler->reg("Mv_Finalize");
}

//*******************************************************************************
void opp_exit() 
{
    if (OPP_DBG) opp_printf("opp_exit", "");

    globalMover.reset();
    cellMapper.reset();
    boundingBox.reset();
    comm.reset();
    
    {   
        opp_halo_destroy(); // free memory allocated to halos and mpi_buffers 
        opp_partition_destroy(); // free memory used for holding partition information
        opp_part_comm_destroy(); // free memory allocated for particle communication
        
        // opp_host_free(set_import_buffer_size);

        // for (int i = 0; i < OPP_import_index; i++)
        //     opp_host_free(OPP_import_list[i]);
        // if (OPP_import_list)
        //     opp_host_free(OPP_import_list);
        
        // for (int i = 0; i < OPP_export_index; i++)
        //     opp_host_free(OPP_export_list[i]);
        // if (OPP_export_list)
        //     opp_host_free(OPP_export_list);
    }

    opp_exit_core();

#ifdef USE_PETSC
    PetscFinalize();
#else
    MPI_Finalize();
#endif
    
}

//****************************************
void opp_abort(std::string s)
{
    opp_printf("opp_abort", "%s", s.c_str());
    MPI_Abort(MPI_COMM_WORLD, 2);
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

    if (OPP_DBG) opp_printf("opp_decl_map", OPP_rank, " map: %s | ptr: %p | dim: %d", map->name, map->map, map->dim);

    return map;
}

//****************************************
opp_dat opp_decl_mesh_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    return opp_decl_dat_core(set, dim, type.c_str(), size, (char*)data, name);
}

//****************************************
opp_dat opp_decl_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name, bool cell_index)
{
    if (set->is_particle) 
        return opp_decl_part_dat(set, dim, dtype, data, name, cell_index);
    else
        return opp_decl_mesh_dat(set, dim, dtype, data, name);
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
    return opp_decl_particle_set_core(size, name, cells_set);
}
opp_set opp_decl_particle_set(char const *name, opp_set cells_set)
{
    return opp_decl_particle_set_core(name, cells_set);
}

//****************************************
opp_dat opp_decl_part_dat(opp_set set, int dim, opp_data_type dtype, void *data, char const *name, bool cell_index)
{
    std::string type = "";
    int size = -1;
    getDatTypeSize(dtype, type, size);

    return opp_decl_particle_dat_core(set, dim, type.c_str(), size, (char*)data, name, cell_index);
}

//****************************************
opp_dat opp_decl_particle_dat_txt(opp_set set, int dim, opp_data_type dtype, const char* file_name, char const *name, bool cell_index)
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
void opp_increase_particle_count(opp_set particles_set, const int num_particles_to_insert)
{
    if (!opp_increase_particle_count_core(particles_set, num_particles_to_insert))
    {
        opp_printf("opp_increase_particle_count", "Error: opp_increase_particle_count_core failed for particle set [%s]", particles_set->name);
        opp_abort("opp_increase_particle_count");        
    }
}

//****************************************
void opp_inc_part_count_with_distribution(opp_set particles_set, int num_particles_to_insert, opp_dat part_dist, bool calc_new)
{
    if (OPP_DBG) opp_printf("opp_inc_part_count_with_distribution", "num_particles_to_insert [%d]", num_particles_to_insert);

    if (!opp_inc_part_count_with_distribution_core(particles_set, num_particles_to_insert, part_dist))
    {
        opp_printf("opp_inc_part_count_with_distribution", "Error: opp_inc_part_count_with_distribution_core failed for particle set [%s]", particles_set->name);
        opp_abort("opp_inc_part_count_with_distribution_core");        
    }
}

//****************************************
void opp_reset_num_particles_to_insert(opp_set set) // unused
{
    opp_reset_num_particles_to_insert_core(set);
}

//****************************************
void opp_particle_sort(opp_set set) // unused
{ 
    opp_particle_sort_core(set);
}

//****************************************
void opp_print_dat_to_txtfile(opp_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{
    std::string prefix = std::string(file_name_prefix) + "_r" + std::to_string(OPP_rank) + "_m" + std::to_string(OPP_comm_size);
    opp_print_dat_to_txtfile_core(dat, prefix.c_str(), file_name_suffix);
}

//****************************************
void opp_print_map_to_txtfile(opp_map map, const char *file_name_prefix, const char *file_name_suffix)
{
    std::string prefix = std::string(file_name_prefix) + "_r" + std::to_string(OPP_rank) + "_m" + std::to_string(OPP_comm_size);
    opp_print_map_to_txtfile_core(map, prefix.c_str(), file_name_suffix);
}

//****************************************
void opp_init_particle_move(opp_set set, int nargs, opp_arg *args)
{ 

    opp_init_particle_move_core(set);

    opp_move_part_indices.clear();
    opp_move_part_indices.reserve(20000);

    if (OPP_comm_iteration == 0)
    {
        OPP_iter_start = 0;
        OPP_iter_end   = set->size;          
        OPP_part_comm_count_per_iter = 0; 
    }
    else
    {
        // need to change the arg data since particle communication could change the pointer in realloc dat->data
        for (int i = 0; i < nargs; i++)
        {
            if (args[i].argtype == OPP_ARG_DAT && args[i].dat->set->is_particle)
            {
                args[i].data = args[i].dat->data;
                if (OPP_DBG) opp_printf("SSSS", "dat %s", args[i].dat->name);
            }
        }
    }

    if (OPP_DBG) opp_printf("opp_init_particle_move", "comm_iter=%d start=%d end=%d", 
                OPP_comm_iteration, OPP_iter_start, OPP_iter_end);

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 
}

//****************************************
bool opp_finalize_particle_move(opp_set set)
{ 

    if (OPP_DBG) opp_printf("opp_finalize_particle_move", "Start particle set [%s]", set->name);

    opp_profiler->start("Mv_Finalize");

    OPP_part_comm_count_per_iter += (int)opp_move_part_indices.size();

    opp_process_marked_particles(set); 

    opp_part_pack(set);
    
    // send the counts and send the particles  
    opp_part_exchange(set);  

    // Can fill the holes here, since the communicated particles will be added at the end
    opp_finalize_particle_move_core(set);

    if (OPP_auto_sort == 1)
    {
        if (OPP_DBG) opp_printf("opp_finalize_particle_move", "auto sorting particle set [%s]", set->name);
        opp_particle_sort(set);
    }

    if (opp_part_check_all_done(set))
    {
        if (OPP_max_comm_iteration < OPP_comm_iteration)
            OPP_max_comm_iteration = OPP_comm_iteration;

        OPP_comm_iteration = 0; // reset for the next par loop
        
        opp_profiler->end("Mv_Finalize");

        return false; // all mpi ranks do not have anything to communicate to any rank
    }
        
    opp_part_wait_all(set); // wait till all the particles are communicated
    
    // increase the particle count if required and unpack the communicated particle buffer 
    // in to separate particle dats
    opp_part_unpack(set);

    // cleanSendRecvBuffers(set);

    OPP_iter_start = set->size - set->diff;
    OPP_iter_end   = set->size;  

    OPP_comm_iteration++;  

    opp_profiler->end("Mv_Finalize");

    return true;
}

//****************************************
void opp_reset_dat_impl(opp_dat dat, char* val, opp_reset reset)
{
    if (!val)
    {
        opp_printf("opp_reset_dat_impl", "Error: val is NULL");
        return;
    }

    int start = -1;
    int end = -1;

    opp_get_start_end(dat->set, reset, start, end);

    for (int i = start; i < end; i++)
    {
        memcpy(dat->data + i * dat->size, val, dat->size);
    }

    // TODO : Check whether this is OK for all the reset options!
    dat->dirtybit = 0;
}

//*******************************************************************************
void opp_partition(std::string lib_name, opp_set prime_set, opp_map prime_map, opp_dat data)
{
    opp_profiler->start("opp_partition");

    // remove all negative mappings and copy the first mapping of the current element for all negative mappings
    opp_sanitize_all_maps();

    opp_partition_core(lib_name, prime_set, prime_map, data);

    opp_desanitize_all_maps();

    MPI_Barrier(MPI_COMM_WORLD);
    
    opp_profiler->end("opp_partition");
}

//*******************************************************************************
void opp_mpi_print_dat_to_txtfile(opp_dat dat, const char *file_name) 
{
    const std::string prefixed_file_name = std::string("mpi_files/MPI_") + 
                std::to_string(OPP_comm_size) + std::string("_iter") + 
                std::to_string(OPP_main_loop_iter) + std::string("_") + file_name;
    // rearrange data back to original order in mpi
    opp_dat temp = opp_mpi_get_data(dat);
    
    print_dat_to_txtfile_mpi(temp, prefixed_file_name.c_str());

    opp_host_free(temp->data);
    opp_host_free(temp->set);
    opp_host_free(temp);
}

//*******************************************************************************
opp_move_var opp_get_move_var(int thread)
{
    // no perf improvement by using a buffered move var, could create a new here instead
    
    move_var.move_status = OPP_MOVE_DONE;

    if (OPP_comm_iteration != 0) // TRUE means communicated particles, no need to do the iteration one calculations
        move_var.iteration_one = false;
    else
        move_var.iteration_one = true;

    return move_var; // passing the object for now :(
}

//*******************************************************************************
opp_dat opp_fetch_data(opp_dat dat) {
    if (dat->set->is_particle) {
        opp_printf("opp_fetch_data", "Error Cannot rearrange particle dats");
        opp_abort();
    }

    // rearrange data backe to original order in mpi
    return opp_mpi_get_data(dat);
}

//*******************************************************************************
// Copy a dat from host to device
void opp_upload_dat(opp_dat dat) {}

//*******************************************************************************
// Copy a dat from device to host
void opp_download_dat(opp_dat dat) {}

//*******************************************************************************
// Copy all dats of the set from device to host
void opp_download_particle_set(opp_set particles_set, bool force_download) {}

//*******************************************************************************
// Copy all dats of the set from host to device
void opp_upload_particle_set(opp_set particles_set, bool realloc) {}

//*******************************************************************************
void opp_colour_cartesian_mesh(const int ndim, std::vector<int> cell_counts, opp_dat cell_index, 
                            const opp_dat cell_colors, const int cell_ghosts) 
{  
    __opp_colour_cartesian_mesh(ndim, cell_counts, cell_index, cell_colors, cell_ghosts);
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
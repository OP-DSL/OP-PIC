
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

#include <opp_lib_core.h>

#include "opp_particle_organize_core.cpp"
#include "opp_increase_part_count_core.cpp"
#include "opp_lib_file_core.cpp"

//****************************************
std::vector<opp_set> opp_sets;
std::vector<opp_map> opp_maps;
std::vector<opp_dat> opp_dats;

// int OPP_hybrid_gpu                   = 0;
int OPP_maps_base_index              = 0;
int OPP_auto_soa                     = 0;
int OPP_gpu_direct                   = 0;
double OPP_part_alloc_mult           = 1;
int OPP_fill_period                 = 2;
opp_fill_type OPP_fill_type         = OPP_HoleFill_All;
int OPP_mpi_part_alloc_mult         = 1;
int OPP_rank                        = 0;
int OPP_comm_size                   = 1;
int OPP_comm_iteration              = 0;
int OPP_max_comm_iteration          = 0;
int OPP_iter_start                  = 0;
int OPP_iter_end                    = 0;
int OPP_main_loop_iter              = 0;
int OPP_gpu_threads_per_block       = OPP_DEFAULT_GPU_THREADS_PER_BLOCK;
size_t OPP_gpu_shared_mem_per_block = 0;
int *OPP_mesh_relation_data         = nullptr;
int *OPP_mesh_relation_data_d       = nullptr;
int OPP_part_cells_set_size         = 0;
int OPP_part_comm_count_per_iter    = 0;
int OPP_move_max_hops               = -1;
int OPP_move_moreX_hops             = -1;
bool OPP_dh_data_dump               = false;
size_t opp_maxSavedDHGrid[3];
size_t opp_minSavedDHGrid[3];

std::unique_ptr<opp::Params> opp_params;
std::unique_ptr<opp::Profiler> opp_profiler;

//****************************************
void opp_init_core(int argc, char **argv) 
{
    opp_sets.clear();
    opp_maps.clear();
    opp_dats.clear();

    opp_params = std::make_unique<opp::Params>(argv[1]);
    opp_profiler = std::make_unique<opp::Profiler>();

    // these will be overidden by args
    OPP_part_alloc_mult = opp_params->get<OPP_REAL>("opp_allocation_multiple");

    if (opp_params->get<OPP_STRING>("opp_fill") == "HoleFill_All")
        OPP_fill_type = OPP_HoleFill_All;
    else if (opp_params->get<OPP_STRING>("opp_fill") == "Sort_All")
        OPP_fill_type = OPP_Sort_All;
    else if (opp_params->get<OPP_STRING>("opp_fill") == "Shuffle_All")
        OPP_fill_type = OPP_Shuffle_All;
    else if (opp_params->get<OPP_STRING>("opp_fill") == "Sort_Periodic") {
        OPP_fill_type = OPP_Sort_Periodic;
        if (opp_params->get<OPP_INT>("opp_fill_period") >= 0)
            OPP_fill_period = opp_params->get<OPP_INT>("opp_fill_period");
    }
    else if (opp_params->get<OPP_STRING>("opp_fill") == "Shuffle_Periodic") {
        OPP_fill_type = OPP_Shuffle_Periodic;
        if (opp_params->get<OPP_INT>("opp_fill_period") >= 0)
            OPP_fill_period = opp_params->get<OPP_INT>("opp_fill_period");
    }

    for (int n = 1; n < argc; n++) {
        opp_set_args_core(argv[n]);
    }

    opp_profiler->reg("Mv_holefill_core");
    opp_profiler->reg("Mv_sort_core");
    opp_profiler->reg("Mv_shuffle_core");
}

//****************************************
void opp_exit_core() 
{ 
    if (opp_profiler.get()) {    
        if (opp_params->get<OPP_BOOL>("opp_profile_all"))
            opp_profiler->printProfile(true);
        opp_profiler->printProfile();
    }

    for (auto& a : opp_maps) {
        opp_host_free(a->map);
        opp_host_free((char*)a->name);
        opp_host_free(a);
    }
    opp_maps.clear();

    for (auto& dat : opp_dats) {
        opp_host_free(dat->data);
        delete dat->thread_data;
        opp_host_free((char*)dat->name);
        opp_host_free((char*)dat->type);
        opp_host_free(dat);
    }
    opp_dats.clear();

    for (auto& a : opp_sets) {
        delete a->indexes_to_remove;
        delete a->particle_dats;
        delete a->cell_index_v_part_index_map;
        opp_host_free((char*)a->name);
        opp_host_free(a);
    }
    opp_sets.clear();
}

//****************************************
void opp_set_args_core(char *argv) 
{
    char temp[64];
    char *pch;

    pch = strstr(argv, "OPP_ALLOC_MULT=");
    if (pch != NULL) {
        strncpy(temp, pch, 17);
        OPP_part_alloc_mult = atof(temp + 15);
        
        if (OPP_rank == OPP_ROOT)
            printf("\topp_set_args_core OPP_part_alloc_mult = %lf\n", OPP_part_alloc_mult);
        return;
    }

    pch = strstr(argv, "OPP_MPI_ALLOC_MULT=");
    if (pch != NULL) {
        strncpy(temp, pch, 21);
        OPP_mpi_part_alloc_mult = atoi(temp + 19);
        
        if (OPP_rank == OPP_ROOT)
            printf("\topp_set_args_core OPP_mpi_part_alloc_mult = %d\n", OPP_mpi_part_alloc_mult);
        return;
    }

    pch = strstr(argv, "OPP_DH_DATA_DUMP=");
    if (pch != NULL) {
        strncpy(temp, pch, 18);
        int dump = atoi(temp + 17);
        OPP_dh_data_dump = (dump != 0);
        
        if (OPP_rank == OPP_ROOT)
            printf("\topp_set_args_core OPP_dh_data_dump = %s %d\n", OPP_dh_data_dump ? "TRUE" : "FALSE", dump);
        return;
    }
}

//****************************************
opp_set opp_decl_set_core(int size, char const *name) 
{
    opp_set set          = (opp_set)opp_host_malloc(sizeof(opp_set_core));
    set->index             = opp_sets.size();
    set->size              = size;
    set->name              = copy_str(name);
    set->core_size         = size;
    set->exec_size         = 0;
    set->nonexec_size      = 0;

    set->is_particle       = false;
    set->set_capacity      = size;
    set->diff              = 0;
    set->mesh_relation_dat = NULL;
    set->cells_set         = NULL;
    set->particle_size     = 0;

    set->indexes_to_remove           = new std::vector<int>();
    set->particle_dats               = new std::vector<opp_dat>();
    set->cell_index_v_part_index_map = new std::map<int, part_index>();

    opp_sets.push_back(set);
    return set;
}

//****************************************
// opp_decl_particle_set should make an opp_set which changes size during the loop
opp_set opp_decl_particle_set_core(int size, char const *name, opp_set cells_set) 
{   
    opp_set set = opp_decl_set_core(size, name);

    set->is_particle = true;
    set->cells_set = cells_set;

    return set;
}

opp_set opp_decl_particle_set_core(char const *name, opp_set cells_set) 
{
    return opp_decl_particle_set_core(0, name, cells_set);
}

//****************************************
opp_map opp_decl_map_core(opp_set from, opp_set to, int dim, OPP_INT *imap, char const *name) 
{
    if (from == NULL) {
        opp_printf("opp_decl_map", "Error -- invalid 'from' set for map %s", name);
        exit(-1);
    }

    if (to == NULL) {
        opp_printf("opp_decl_map", "Error -- invalid 'to' set for map %s", name);
        exit(-1);
    }

    if (dim <= 0) {
        opp_printf("opp_decl_map", "Error -- negative/zero dimension for map %s", name);
        exit(-1);
    }

    if (to->is_particle) {
        opp_printf("opp_decl_map", "Error -- cannot have mappings to a dynamic (particle) set [%s to %s]", 
            from->name, to->name);
        exit(-1);
    }

    opp_map map  = (opp_map)opp_host_malloc(sizeof(opp_map_core));
    map->index    = opp_maps.size();
    map->from     = from;
    map->to       = to;
    map->dim      = dim;

    map->map      = (OPP_INT *)opp_host_malloc((size_t)from->size * (size_t)dim * sizeof(OPP_INT));
    memcpy(map->map, imap, sizeof(OPP_INT) * from->size * dim);  

    if (OPP_maps_base_index == 1) { // convert map to 0 based indexing -- i.e. reduce each map value by 1
        for (size_t i = 0; i < (size_t)from->size * dim; i++)
            (map->map[i])--;
    }

    map->map_d        = NULL;
    map->name         = copy_str(name);
    map->user_managed = 1;
    map->p2c_dat      = NULL;

    if (!from->is_particle) {
        opp_maps.push_back(map);
    }
    else {
        map->p2c_dat = opp_decl_dat_core(from, dim, "int", sizeof(OPP_INT), (char*)map->map, name);

        from->mesh_relation_dat = map->p2c_dat;
        map->p2c_dat->is_cell_index = true;

        if (map->p2c_dat->data != NULL) {
            opp_host_free(map->map);
        }
        map->map = NULL;
        map->p2c_dat->p2c_map = map;
    }
    
    return map;
}

//****************************************
opp_dat opp_decl_dat_core(opp_set set, int dim, char const *type, int size, char *data, char const *name) 
{
    if (set == NULL) {
        opp_printf("opp_decl_dat", "Error -- invalid set for data: %s\n", name);
        exit(-1);
    }

    if (dim <= 0) {
        opp_printf("opp_decl_dat", "Error -- negative/zero dimension for data: %s\n", name);
        exit(-1);
    }

    opp_dat dat = (opp_dat)opp_host_malloc(sizeof(opp_dat_core));
    dat->index    = opp_dats.size();
    dat->set      = set;
    dat->dim      = dim;

    if (set->size > 0) {
        const size_t bytes = (size_t)dim * (size_t)size * 
                            (size_t)(set->size + set->exec_size + set->nonexec_size) * sizeof(char);
        dat->data = (char *)opp_host_malloc(bytes);
        memcpy(dat->data, data, (size_t)dim * (size_t)size * (size_t)set->size * sizeof(char));
    }
    else {
        dat->data = NULL;
    }

    dat->data_d        = NULL;
    dat->data_swap_d   = NULL;
    dat->name          = copy_str(name);
    dat->type          = copy_str(type);
    dat->size          = dim * size;
    dat->user_managed  = 1;
    dat->mpi_buffer    = NULL;
    dat->mpi_reduc_buffer = NULL;
    dat->reduc_comm    = OPP_Reduc_NO_Comm;
    dat->buffer_d      = NULL;
    dat->buffer_d_r    = NULL;
    dat->dirty_hd      = Dirty::NotDirty;
    dat->dirtybit      = 0;

    dat->thread_data        = new std::vector<char*>();
    dat->thread_data_d      = NULL;
    dat->is_cell_index      = false;
    dat->thrust_int         = NULL;
    dat->thrust_real        = NULL;
    dat->thrust_int_sort    = NULL;
    dat->thrust_real_sort   = NULL;
    dat->p2c_map            = NULL;

    set->particle_size += dat->size;

    opp_dats.push_back(dat);

    if (set->is_particle) {
        set->particle_dats->push_back(dat); 
    }

    return dat;
}

//****************************************
opp_arg opp_arg_dat_core(opp_dat dat, int idx, opp_map map, int dim, const char *typ, 
                                opp_dat p2c_map, opp_access acc) 
{
    opp_arg arg;
    arg.index       = -1;
    arg.argtype     = OPP_ARG_DAT;

    arg.dat         = dat;
    arg.map         = map;
    arg.p2c_map     = p2c_map;
    arg.dim         = dim;
    arg.idx         = idx;
    
    arg.size        = ((map == NULL) ? 
                        -1 : (map->from->size + map->from->exec_size + map->from->nonexec_size));
    arg.data        = ((dat == NULL) ? NULL : dat->data);
    arg.data_d      = ((dat == NULL) ? NULL : dat->data_d);
    arg.map_data    = ((idx == -1 || idx == -2) ? NULL : map->map);
    arg.map_data_d  = ((idx == -1 || idx == -2) ? NULL : map->map_d);

    arg.type        = typ;     // Not used
    arg.acc         = acc;
    arg.opt         = 1;
    arg.sent        = 0;
    
    return arg;
}

// arg.map has the map, can change to mapping data map if required
opp_arg opp_arg_dat_core(opp_map data_map, int idx, opp_map map, opp_dat p2c_map, opp_access acc)
{
    opp_arg arg;
    arg.argtype     = OPP_ARG_MAP;

    arg.dat         = NULL;
    arg.map         = map;
    arg.p2c_map     = p2c_map;
    arg.dim         = data_map->dim;
    arg.idx         = idx;
    
    arg.size        = (data_map->from->size + data_map->from->exec_size + data_map->from->nonexec_size);
    arg.data        = (char*)data_map->map;
    arg.data_d      = (char*)data_map->map_d;
    arg.map_data    = ((map == NULL) ? NULL : map->map);
    arg.map_data_d  = ((map == NULL) ? NULL : map->map_d);

    arg.type        = "int";
    arg.acc         = acc;
    arg.opt         = 1;

    return arg;
}

//****************************************
// Set Dirty::Device, making HOST data to be clean
void opp_set_dirtybit(int nargs, opp_arg *args) 
{
    for (int n = 0; n < nargs; n++) {
        if ((args[n].opt == 1) && (args[n].argtype == OPP_ARG_DAT) &&
                (args[n].acc == OPP_INC || args[n].acc == OPP_WRITE || args[n].acc == OPP_RW)) {
            if (OPP_DBG) 
                opp_printf("opp_set_dirtybit", "Setting Dirty::Device| %s", args[n].dat->name);
            args[n].dat->dirty_hd = Dirty::Device;
            if (!args[n].dat->set->is_particle)
                args[n].dat->dirtybit = 1;
        }
    }
}

//****************************************
// Set Dirty::Host, making DEVICE data to be clean
void opp_set_dirtybit_device(int nargs, opp_arg *args) 
{
    for (int n = 0; n < nargs; n++) {
        if ((args[n].opt == 1) && (args[n].argtype == OPP_ARG_DAT) &&
                (args[n].acc == OPP_INC || args[n].acc == OPP_WRITE || args[n].acc == OPP_RW)) {
            if (OPP_DBG) 
                opp_printf("opp_set_dirtybit_device", "Setting Dirty::Host| %s", args[n].dat->name);          
            args[n].dat->dirty_hd = Dirty::Host;
            if (!args[n].dat->set->is_particle)
                args[n].dat->dirtybit = 1;
        }
    }
}

//****************************************
void opp_set_dirtybit_grouped(int nargs, opp_arg *args, DeviceType device)
{
    return (device == Device_CPU) ? 
                opp_set_dirtybit(nargs, args) : opp_set_dirtybit_device(nargs, args);
}

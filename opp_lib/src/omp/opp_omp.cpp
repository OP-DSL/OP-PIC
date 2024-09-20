
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

int OPP_nthreads = 1;
std::vector<int> part_remove_count_per_thr;
std::vector<std::vector<int>> move_part_indices_per_thr;
std::vector<std::vector<opp_tmp_gbl_move_info>> gbl_move_indices_per_thr;
std::vector<opp_move_var> move_var_per_thr;

void opp_part_pack(opp_set set);
void opp_part_unpack(opp_set set);

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

    OPP_nthreads = omp_get_max_threads();

    if (OPP_rank == OPP_ROOT)
    {
        std::string log = "Running on OMP";
#ifdef USE_MPI
        log += "+MPI with " + std::to_string(OPP_comm_size) + " ranks";
#endif    
        log += " and " + std::to_string(OPP_nthreads) + " threads per rank";
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
        exit(-1);        
    }
}

//****************************************
bool opp_inc_part_count_with_distribution_omp(opp_set particles_set, int num_particles_to_insert, opp_dat part_dist)
{
    opp_profiler->start("IncPartCountWithDistribution");

    if (!opp_increase_particle_count_core(particles_set, num_particles_to_insert))
    {
        opp_profiler->end("IncPartCountWithDistribution");
        return false;
    }

    int* part_mesh_connectivity = (int *)particles_set->mesh_relation_dat->data;
    int* distribution           = (int *)part_dist->data;

    int startX = (particles_set->size - particles_set->diff);
    int distribution_set_size = part_dist->set->size;

    int nthreads = omp_get_max_threads();

    int set_size = particles_set->diff;

    #pragma omp parallel for
    for (int thr = 0; thr < nthreads; thr++)
    {
        size_t start  = ((size_t)set_size * thr) / nthreads;
        size_t finish = ((size_t)set_size * (thr+1)) / nthreads;

        for (size_t i = start; i < finish; i++)
        { 
            // iterating for all particles here, makes the performance hit.
            // can we change here?
            for (int j = 0; j < distribution_set_size; j++) 
            {
                if (i < distribution[j]) 
                {
                    part_mesh_connectivity[startX + i] = j;
                    break;
                }  
            }
        }
    } 

    opp_profiler->end("IncPartCountWithDistribution");

    return true;
}

//****************************************
void opp_inc_part_count_with_distribution(opp_set particles_set, int num_particles_to_insert, opp_dat part_dist, bool calc_new)
{
    if (OPP_DBG) opp_printf("opp_inc_part_count_with_distribution", "num_particles_to_insert [%d]", num_particles_to_insert);

    // if (!opp_inc_part_count_with_distribution_omp(particles_set, num_particles_to_insert, part_dist))
    // {
    //     opp_printf("opp_inc_part_count_with_distribution", "Error: opp_inc_part_count_with_distribution_omp failed for particle set [%s]", particles_set->name);
    //     exit(-1);        
    // }

    // perf of core is better than omp :(

    if (!opp_inc_part_count_with_distribution_core(particles_set, num_particles_to_insert, part_dist))
    {
        opp_printf("opp_inc_part_count_with_distribution", "Error: opp_inc_part_count_with_distribution_core failed for particle set [%s]", particles_set->name);
        exit(-1);        
    }
}

//****************************************
void opp_reset_num_particles_to_insert(opp_set set)
{
    opp_reset_num_particles_to_insert_core(set);
}

//****************************************
void opp_mark_particle_to_remove(opp_set set, int particle_index)
{
    opp_mark_particle_to_remove_core(set, particle_index);
}

//****************************************
void opp_remove_marked_particles_from_set(opp_set set)
{
    opp_remove_marked_particles_from_set_core(set);
}
void opp_remove_marked_particles_from_set(opp_set set, std::vector<int>& idx_to_remove)
{
    opp_remove_marked_particles_from_set_core(set, idx_to_remove);
}

//****************************************
void opp_init_particle_move(opp_set set, int nargs, opp_arg *args)
{
    part_remove_count_per_thr.resize(OPP_nthreads);
    std::fill(part_remove_count_per_thr.begin(), part_remove_count_per_thr.end(), 0);

    move_var_per_thr.clear();
    move_var_per_thr.resize(OPP_nthreads);

#ifdef USE_MPI
    opp_move_part_indices.clear();
    opp_move_part_indices.reserve(20000);

    move_part_indices_per_thr.resize(OPP_nthreads);
    for (int t = 0; t < OPP_nthreads; t++) 
    {
        move_part_indices_per_thr[t].clear();
        move_part_indices_per_thr[t].reserve(1000);
    }
    // std::fill(move_part_indices_per_thr.begin(), move_part_indices_per_thr.end(), std::vector<int>());

    gbl_move_indices_per_thr.resize(OPP_nthreads);
    for (int t = 0; t < OPP_nthreads; t++) 
    {
        gbl_move_indices_per_thr[t].clear();
        gbl_move_indices_per_thr[t].reserve(1000);
    }    
#endif

    if (OPP_comm_iteration == 0)
    {
        OPP_iter_start = 0;
        OPP_iter_end   = set->size;          
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

    OPP_mesh_relation_data = ((int *)set->mesh_relation_dat->data); 

    opp_init_particle_move_core(set);
}

//****************************************
void opp_mark_particle_to_move(opp_set set, int particle_index, int move_status)
{
    opp_mark_particle_to_move_core(set, particle_index, move_status);
}

//****************************************
bool opp_finalize_particle_move(opp_set set)
{ 

    if (OPP_DBG) opp_printf("opp_finalize_particle_move", "Start particle set [%s]", set->name);

    opp_profiler->start("Mv_Finalize");

    for (int i = 0; i < OPP_nthreads; i++) 
        set->particle_remove_count += part_remove_count_per_thr[i];

#ifdef USE_MPI
    // process indices of each thread now
    for (int i = 0; i < OPP_nthreads; i++)
    {
        opp_move_part_indices.insert(opp_move_part_indices.end(), 
            move_part_indices_per_thr[i].begin(), move_part_indices_per_thr[i].end());
    }

    opp_process_marked_particles(set);

    opp_part_pack(set);

    // send the counts and send the particles  
    opp_part_exchange(set); 
#endif

    opp_finalize_particle_move_omp(set);

    if (OPP_auto_sort == 1)
    {
        if (OPP_DBG) printf("\topp_finalize_particle_move auto sorting particle set [%s]\n", set->name);
        opp_particle_sort(set);
    }

#ifdef USE_MPI
    if (opp_part_check_all_done(set))
    {
        if (OPP_max_comm_iteration < OPP_comm_iteration)
            OPP_max_comm_iteration = OPP_comm_iteration;

        OPP_comm_iteration = 0; // reset for the next par loop
        
        opp_profiler->end("Mv_Finalize");
        return false; // all mpi ranks do not have anything to communicate to any rank
    }

    opp_part_wait_all(set); // wait till all the particles are communicated

    if (OPP_DBG)
        opp_printf("opp_finalize_particle_move", "set [%s] size prior unpack %d", set->name, set->size);

    // increase the particle count if required and unpack the communicated particle buffer 
    // in to separate particle dats
    opp_part_unpack(set);    

    OPP_iter_start = set->size - set->diff;
    OPP_iter_end   = set->size;  

    OPP_comm_iteration++;  

    opp_profiler->end("Mv_Finalize");
    return true;
#else

    opp_profiler->end("Mv_Finalize");
    return false;
#endif
}

#ifdef USE_OLD_FINALIZE

//****************************************
void opp_finalize_particle_move_omp(opp_set set)
{
    if (OPP_DBG) opp_printf("opp_finalize_particle_move_omp", "set [%s] with particle_remove_count [%d]\n", 
        set->name, set->particle_remove_count);

    if (set->particle_remove_count <= 0) return;

    if (OPP_auto_sort == 0) // if not auto sorting, fill the holes
    {
        // getting a backup of cell index since it will also be rearranged using a random OMP thread
        int *mesh_relation_data = (int *)opp_host_malloc(set->set_capacity * set->mesh_relation_dat->size); 
        
        memcpy((char*)mesh_relation_data, set->mesh_relation_dat->data, 
            set->set_capacity * set->mesh_relation_dat->size);

        #pragma omp parallel for
        for (int i = 0; i < (int)set->particle_dats->size(); i++)
        {
            opp_dat current_opp_dat = set->particle_dats->at(i);
            int removed_count = 0;
            int skip_count = 0;

            for (size_t j = 0; j < set->size; j++)
            {
                if (mesh_relation_data[j] != MAX_CELL_INDEX) continue;

                char* dat_removed_ptr = (char *)(current_opp_dat->data + (size_t)(j * current_opp_dat->size));

                // BUG_FIX: (set->size - removed_count - 1) This index could marked to be removed, and if marked, 
                // then there could be an array index out of bounds access error in the future
                while (mesh_relation_data[set->size - removed_count - skip_count - 1] == MAX_CELL_INDEX)
                {
                    skip_count++;
                }
                if (j >= (set->size - removed_count - skip_count - 1)) 
                {
                    if (OPP_DBG) 
                        opp_printf("opp_finalize_particle_move_omp", 
                            "Current Iteration index [%d] and replacement index %d; hence breaking\n", 
                            j, (set->size - removed_count - skip_count - 1));
                    break;
                }

                size_t offset_byte = (size_t)(set->size - removed_count - skip_count - 1) * current_opp_dat->size;
                char* dat_to_replace_ptr = (char *)(current_opp_dat->data + offset_byte);
                
                // Get the last element and replace the hole // Not the Optimum!!!
                // TODO : Can we make NULL data and handle it in sort?
                memcpy(dat_removed_ptr, dat_to_replace_ptr, current_opp_dat->size); 

                removed_count++;
            }

            // current_opp_dat->data = (char *)opp_host_realloc(current_opp_dat->data, 
            //       (size_t)(set->size - removed_count) * (size_t)current_opp_dat->size);
        }

        opp_host_free(mesh_relation_data);
    }

    set->size -= set->particle_remove_count;
    set->particle_remove_count = 0;
}

#else

//****************************************
void opp_finalize_particle_move_omp(opp_set set)
{
    if (OPP_DBG) opp_printf("opp_finalize_particle_move_omp", "set [%s] with particle_remove_count [%d]\n", 
        set->name, set->particle_remove_count);

    if (set->particle_remove_count <= 0) return;

    if (OPP_auto_sort == 0) // if not auto sorting, fill the holes
    {
        int *mesh_relation_data = (int *)set->mesh_relation_dat->data;
        std::vector<std::pair<size_t, size_t>> swap_indices;    // contain hole index and the index from back to swap

        // Idea: The last available element should be copied to the hole
        // In the below scope we try to calculate the element to be swapped with the hole
        {
            // set->particle_remove_count   // the particle count that should be removed
            int removed_count = 0;          // how many elements currently being removed
            int skip_count = 0;             // how many elements from the back is skipped ..
                                            // .. due to that element is also to be removed

            for (size_t j = 0; j < (size_t)set->size; j++)
            {
                // skip if the current index is not to be removed
                if (mesh_relation_data[j] != MAX_CELL_INDEX) 
                    continue;

                // handle if the element from the back is also to be removed
                while ((set->size - removed_count - skip_count - 1 >= 0) && 
                    (mesh_relation_data[set->size - removed_count - skip_count - 1] == MAX_CELL_INDEX))
                {
                    skip_count++;
                }

                // check whether the holes are at the back!
                if ((set->size - removed_count - skip_count - 1 < 0) ||
                    (j >= (size_t)(set->size - removed_count - skip_count - 1))) 
                {
                    if (OPP_DBG) 
                        opp_printf("opp_finalize_particle_move_core", 
                        "Current Iteration index [%d] and replacement index %d; hence breaking [%s]", 
                        j, (set->size - removed_count - skip_count - 1), set->name);
                    break;
                }

                swap_indices.push_back(std::make_pair(j, (size_t)(set->size - removed_count - skip_count - 1)));

                removed_count++;
            }
        }

        // For all the dats, fill the holes using the swap_indices
        #pragma omp parallel for
        for (size_t i = 0; i < set->particle_dats->size(); i++)
        {
            opp_dat dat = set->particle_dats->at(i);

            for (const auto& x : swap_indices)
            {
                char* dat_removed_ptr = (char *)(dat->data + (x.first * dat->size));

                size_t offset_byte = x.second * dat->size;
                char* dat_to_replace_ptr = (char *)(dat->data + offset_byte);
                
                memcpy(dat_removed_ptr, dat_to_replace_ptr, dat->size); 
            }
        }
    }

    set->size -= set->particle_remove_count;
    set->particle_remove_count = 0;
}

#endif

//****************************************
void opp_particle_sort(opp_set set)
{ 
    
    if (OPP_DBG) printf("\topp_particle_sort set [%s]\n", set->name);
    
    int* mesh_relation_data = (int*)set->mesh_relation_dat->data;

    std::vector<size_t> idx_before_sort = sort_indexes(mesh_relation_data, set->size);

    #pragma omp parallel for
    for (int i = 0; i < (int)set->particle_dats->size(); i++)
    {    
        auto& dat = set->particle_dats->at(i);
        char *new_data = (char *)opp_host_malloc(set->set_capacity * dat->size);
        char *old_data = (char*)dat->data;
        
        for (int j = 0; j < set->size; j++)
        {
            memcpy(new_data + j * dat->size, old_data + idx_before_sort[j] * dat->size, dat->size);
        }

        opp_host_free(dat->data);
        dat->data = new_data;

        if (dat->is_cell_index)
        { 
            int* cell_index_array = (int*)dat->data;
            int current_cell_index = -1, previous_cell_index = -1;
            std::map<int, part_index>& map = *(set->cell_index_v_part_index_map);
            map.clear();

            for (int j = 0; j < set->size; j++)
            {    
                current_cell_index = cell_index_array[j];
            
                if ((current_cell_index != previous_cell_index) && (current_cell_index >= 0))
                {
                    part_index& pi = map[current_cell_index];
                    pi.start = j;

                    if (previous_cell_index >= 0) map[previous_cell_index].end = (j - 1);
                }
                previous_cell_index = current_cell_index;
            }
            map[previous_cell_index].end = (set->size - 1);
        }
    }
}

//****************************************
void opp_print_dat_to_txtfile(opp_dat dat, const char *file_name_prefix, const char *file_name_suffix)
{
    std::string prefix = std::string(file_name_prefix) + "_o";
    opp_print_dat_to_txtfile_core(dat, prefix.c_str(), file_name_suffix);
}

//****************************************
void opp_print_map_to_txtfile(opp_map map, const char *file_name_prefix, const char *file_name_suffix)
{
    opp_print_map_to_txtfile_core(map, file_name_prefix, file_name_suffix);
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
    for (int i = start; i < end; i++)
    {
        memcpy(dat->data + i * dat->size, val, dat->size);
    }
}

//****************************************
bool opp_part_check_status_omp(opp_move_var& m, int map0idx, opp_set set, 
    int particle_index, int& remove_count, int thread)
{
    m.iteration_one = false;

    if (m.move_status == OPP_MOVE_DONE)
    {
        return false;
    }
    else if (m.move_status == OPP_NEED_REMOVE)
    {
        part_remove_count_per_thr[thread] += 1;
        OPP_mesh_relation_data[particle_index] = MAX_CELL_INDEX;

        return false;
    }
#ifdef USE_MPI
    else if (map0idx >= OPP_part_cells_set_size)
    {
        // map0idx cell is not owned by the current mpi rank (it is in the import exec halo region), need to communicate
        move_part_indices_per_thr[thread].push_back(particle_index);
        return false;
    }
#endif

    return true;
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
void opp_mpi_halo_wait_all(int nargs, opp_arg *args)
{
    // Nothing to execute here
}
#endif // if USE_MPI is defined, the functions in opp_mpi_particle_comm.cpp and opp_mpi_halo.cpp is used

//****************************************
opp_move_var opp_get_move_var(int thread)
{
    // opp_move_var& move_var = move_var_per_thr[thread];

    // move_var.move_status = OPP_MOVE_DONE;
    // move_var.iteration_one = true;

    opp_move_var move_var;

    if (OPP_comm_iteration != 0) // TRUE means communicated particles, no need to do the iteration one calculations
        move_var.iteration_one = false;
    else
        move_var.iteration_one = true;

    return move_var;
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
                                            std::to_string(OPP_comm_size) + std::string("_") + file_name;

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

#ifdef USE_MPI
//*******************************************************************************
opp_dat opp_fetch_data(opp_dat dat) {
    if (dat->set->is_particle) {
        opp_printf("opp_fetch_data", "Error Cannot rearrange particle dats");
        opp_abort();
    }

    // rearrange data backe to original order in mpi
    return opp_mpi_get_data(dat);
}
#endif
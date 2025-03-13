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

// MPI Halo related global variables

halo_list *OPP_export_exec_list; // EEH list
halo_list *OPP_import_exec_list; // IEH list

halo_list *OPP_import_nonexec_list; // INH list
halo_list *OPP_export_nonexec_list; // ENH list

// global variables to hold partition information on an MPI rank
int OPP_part_index = 0;
part *OPP_part_list;

// Save original partition ranges
int **orig_part_range = NULL;

// Sliding planes data structures
int OPP_import_index = 0, OPP_import_max = 0;
int OPP_export_index = 0, OPP_export_max = 0;
op_import_handle *OPP_import_list = NULL;
op_export_handle *OPP_export_list = NULL;

int **export_exec_list_d = NULL;
int **export_nonexec_list_d = NULL;
int **export_exec_list_disps_d = NULL;
int **export_nonexec_list_disps_d = NULL;
int **import_exec_list_disps_d = NULL;
int **import_nonexec_list_disps_d = NULL;

// Timing
double t1, t2, c1, c2;

/*******************************************************************************
 * Main MPI halo creation routine
 *******************************************************************************/
void __opp_halo_create() 
{
    if (OPP_DBG) opp_printf("opp_halo_create", "start");

    // declare timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    double time;
    double max_time;
    op_timers(&cpu_t1, &wall_t1); // timer start for list create

    // create new communicator for OP mpi operation
    int my_rank, comm_size;
    // MPI_Comm_dup(OPP_MPI_WORLD, &OPP_MPI_WORLD);
    MPI_Comm_rank(OPP_MPI_WORLD, &my_rank);
    MPI_Comm_size(OPP_MPI_WORLD, &comm_size);

    /* Compute global partition range information for each set*/
    int **part_range = (int **)opp_host_malloc((int)opp_sets.size() * sizeof(int *)); // part_range[set][ start and end ]
    get_part_range(part_range, my_rank, comm_size, OPP_MPI_WORLD);

    // save this partition range information if it is not already saved during
    // a call to some partitioning routine
    if (orig_part_range == NULL) 
    {
        orig_part_range = (int **)opp_host_malloc((int)opp_sets.size() * sizeof(int *));
        for (int s = 0; s < (int)opp_sets.size(); s++) 
        {
            opp_set set = opp_sets[s];
            orig_part_range[set->index] = (int *)opp_host_malloc(2 * comm_size * sizeof(int));
            for (int j = 0; j < comm_size; j++) 
            {
                orig_part_range[set->index][2 * j] = part_range[set->index][2 * j];
                orig_part_range[set->index][2 * j + 1] = part_range[set->index][2 * j + 1];
            }
        }
    }

    OPP_export_exec_list = (halo_list *)opp_host_malloc((int)opp_sets.size() * sizeof(halo_list));

    /*----- STEP 1 - Construct export lists for execute set elements and related mapping table entries -----*/

    // declare temporaty scratch variables to hold set export lists and mapping table export lists
    int s_i;
    int *set_list;

    int cap_s = 1000; // keep track of the temp array capacities

    for (size_t s = 0; s < opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        // create a temporaty scratch space to hold export list for this set
        s_i = 0;
        cap_s = 1000;
        set_list = (int *)opp_host_malloc(cap_s * sizeof(int));

        for (int e = 0; e < set->size; e++) // for each elment of this set
        {      
            for (size_t m = 0; m < opp_maps.size(); m++) // for each maping table
            { 
                opp_map map = opp_maps[m];

                if (compare_sets(map->from, set) == 1) // need to select mappings FROM this set
                { 
                    int part, local_index;
                    for (int j = 0; j < map->dim; j++) // for each element pointed at by this entry
                    { 
                        part = get_partition(map->map[e * map->dim + j],
                                            part_range[map->to->index], &local_index,
                                            comm_size);
                        if (s_i >= cap_s) 
                        {
                            cap_s = cap_s * 2;
                            set_list = (int *)opp_host_realloc(set_list, cap_s * sizeof(int));
                        }

                        if (part != my_rank) 
                        {
                            set_list[s_i++] = part; // add to set export list
                            set_list[s_i++] = e;
                        }
                    }
                }
            }
        }

        // create set export list
        halo_list h_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));
        create_export_list(set, set_list, h_list, s_i, comm_size, my_rank);
        OPP_export_exec_list[set->index] = h_list;
        opp_host_free(set_list); // free temp list
    }

    /*---- STEP 2 - construct import lists for mappings and execute sets------*/

    OPP_import_exec_list = (halo_list *)opp_host_malloc((int)opp_sets.size() * sizeof(halo_list));

    int *neighbors, *sizes;
    int ranks_size;

    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        //-----Discover neighbors-----
        ranks_size = 0;
        neighbors = (int *)opp_host_malloc(comm_size * sizeof(int));
        sizes = (int *)opp_host_malloc(comm_size * sizeof(int));

        halo_list list = OPP_export_exec_list[set->index];

        find_neighbors_set(list, neighbors, sizes, &ranks_size, my_rank, comm_size, OPP_MPI_WORLD);
        std::vector<MPI_Request> request_send(list->ranks_size);

        int *rbuf, cap = 0, index = 0;

        for (int i = 0; i < list->ranks_size; i++) 
        {
            // printf("export from %d to %d set %10s, list of size %d \n", 
            //  my_rank,list->ranks[i],set->name,list->sizes[i]);
            int *sbuf = &list->list[list->disps[i]];
            MPI_Isend(sbuf, list->sizes[i], MPI_INT, list->ranks[i], s, OPP_MPI_WORLD, &request_send[i]);
        }

        for (int i = 0; i < ranks_size; i++)
            cap = cap + sizes[i];
        int *temp = (int *)opp_host_malloc(cap * sizeof(int));

        // import this list from those neighbors
        for (int i = 0; i < ranks_size; i++) 
        {
            // printf("import from %d to %d set %10s, list of size %d\n", 
            //  neighbors[i], my_rank, set->name, sizes[i]);
            rbuf = (int *)opp_host_malloc(sizes[i] * sizeof(int));
            MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], s, OPP_MPI_WORLD, MPI_STATUS_IGNORE);
            memcpy(&temp[index], (void *)&rbuf[0], sizes[i] * sizeof(int));
            index = index + sizes[i];
            opp_host_free(rbuf);
        }

        MPI_Waitall(list->ranks_size, request_send.data(), MPI_STATUSES_IGNORE);

        // create import lists
        halo_list h_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));
        create_import_list(set, temp, h_list, index, neighbors, sizes, ranks_size, comm_size, my_rank);
        OPP_import_exec_list[set->index] = h_list;
    }

    /*--STEP 3 -Exchange mapping table entries using the import/export lists--*/

    for (int m = 0; m < (int)opp_maps.size(); m++) // for each maping table
    { 
        opp_map map = opp_maps[m];
        halo_list i_list = OPP_import_exec_list[map->from->index];
        halo_list e_list = OPP_export_exec_list[map->from->index];

        std::vector<MPI_Request> request_send(e_list->ranks_size);

        // prepare bits of the mapping tables to be exported
        int **sbuf = (int **)opp_host_malloc(e_list->ranks_size * sizeof(int *));

        for (int i = 0; i < e_list->ranks_size; i++) 
        {
            sbuf[i] = (int *)opp_host_malloc((size_t)e_list->sizes[i] * map->dim * sizeof(int));
            for (int j = 0; j < e_list->sizes[i]; j++) 
            {
                for (int p = 0; p < map->dim; p++) 
                {
                    sbuf[i][j * map->dim + p] = map->map[map->dim * (e_list->list[e_list->disps[i] + j]) + p];
                }
            }

            MPI_Isend(sbuf[i], map->dim * e_list->sizes[i], MPI_INT, e_list->ranks[i],
                        m, OPP_MPI_WORLD, &request_send[i]);
        }

        // prepare space for the incomming mapping tables - realloc each mapping tables in each mpi process
        opp_maps[map->index]->map = (int *)opp_host_realloc(opp_maps[map->index]->map,
            (map->dim * (size_t)(map->from->size + i_list->size)) * sizeof(int));

        int init = map->dim * (map->from->size);
        for (int i = 0; i < i_list->ranks_size; i++) 
        {
            MPI_Recv( &(opp_maps[map->index]->map[init + i_list->disps[i] * map->dim]),
                map->dim * i_list->sizes[i], MPI_INT, i_list->ranks[i], m, OPP_MPI_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Waitall(e_list->ranks_size, request_send.data(), MPI_STATUSES_IGNORE);
        
        for (int i = 0; i < e_list->ranks_size; i++)
            opp_host_free(sbuf[i]);
        opp_host_free(sbuf);
    }

    /*-- STEP 4 - Create import lists for non-execute set elements using mapping
        table entries including the additional mapping table entries --*/

    OPP_import_nonexec_list = (halo_list *)opp_host_malloc((int)opp_sets.size() * sizeof(halo_list));
    OPP_export_nonexec_list = (halo_list *)opp_host_malloc((int)opp_sets.size() * sizeof(halo_list));

    // declare temporaty scratch variables to hold non-exec set export lists
    s_i = 0;
    set_list = NULL;
    cap_s = 1000; // keep track of the temp array capacity

    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;   

        halo_list exec_set_list = OPP_import_exec_list[set->index];

        // create a temporaty scratch space to hold nonexec export list for this set
        s_i = 0;
        set_list = (int *)opp_host_malloc(cap_s * sizeof(int));

        for (int m = 0; m < (int)opp_maps.size(); m++) // for each maping table
        { 
            opp_map map = opp_maps[m];
            halo_list exec_map_list = OPP_import_exec_list[map->from->index];

            // need to select  mappings TO this set
            if (compare_sets(map->to, set) == 1 
#ifdef NO_EXTENDED_HALOS // this will not create non exec halos for same set mappings (cell to cell)
                && compare_sets(map->from, map->to) != 1
#endif
            )
            { 
                // for each entry in this mapping table: original+execlist
                int len = map->from->size + exec_map_list->size;
                for (int e = 0; e < len; e++) 
                {
                    int part;
                    int local_index;
                    for (int j = 0; j < map->dim; j++) // for each element pointed at by this entry
                    { 
                        part = get_partition(map->map[e * map->dim + j],
                                            part_range[map->to->index], &local_index, comm_size);

                        if (s_i >= cap_s) 
                        {
                            cap_s = cap_s * 2;
                            set_list = (int *)opp_host_realloc(set_list, cap_s * sizeof(int));
                        }

                        if (part != my_rank) 
                        {
                            int found = -1;
                            // check in exec list
                            int rank = binary_search(exec_set_list->ranks, part, 0, exec_set_list->ranks_size - 1);

                            if (rank >= 0) 
                            {
                                found = binary_search(exec_set_list->list, local_index,
                                                    exec_set_list->disps[rank],
                                                    exec_set_list->disps[rank] + exec_set_list->sizes[rank] - 1);
                            }

                            if (found < 0) 
                            {
                                // not in this partition and not found in  exec list
                                // add to non-execute set_list
                                set_list[s_i++] = part;
                                set_list[s_i++] = local_index;
                            }
                        }
                    }
                }
            }
        }

        // create non-exec set import list
        halo_list h_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));
        create_nonexec_import_list(set, set_list, h_list, s_i, comm_size, my_rank);
        opp_host_free(set_list); // free temp list
        OPP_import_nonexec_list[set->index] = h_list;
    }

    /*----------- STEP 5 - construct non-execute set export lists -------------*/

    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        //-----Discover neighbors-----
        ranks_size = 0;
        neighbors = (int *)opp_host_malloc(comm_size * sizeof(int));
        sizes = (int *)opp_host_malloc(comm_size * sizeof(int));

        halo_list list = OPP_import_nonexec_list[set->index];
        find_neighbors_set(list, neighbors, sizes, &ranks_size, my_rank, comm_size,
                        OPP_MPI_WORLD);

        std::vector<MPI_Request> request_send(list->ranks_size);
        int *rbuf, cap = 0, index = 0;

        for (int i = 0; i < list->ranks_size; i++) 
        {
            int *sbuf = &list->list[list->disps[i]];
            MPI_Isend(sbuf, list->sizes[i], MPI_INT, list->ranks[i], s, OPP_MPI_WORLD, &request_send[i]);
        }

        for (int i = 0; i < ranks_size; i++)
            cap = cap + sizes[i];
        int *temp = (int *)opp_host_malloc(cap * sizeof(int));

        // export this list to those neighbors
        for (int i = 0; i < ranks_size; i++) 
        {
            rbuf = (int *)opp_host_malloc(sizes[i] * sizeof(int));
            MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], s, OPP_MPI_WORLD, MPI_STATUS_IGNORE);
            memcpy(&temp[index], (void *)&rbuf[0], sizes[i] * sizeof(int));
            index = index + sizes[i];
            opp_host_free(rbuf);
        }

        MPI_Waitall(list->ranks_size, request_send.data(), MPI_STATUSES_IGNORE);

        halo_list h_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));
        create_nonexec_export_list(set, temp, h_list, index, neighbors, sizes, ranks_size, comm_size, my_rank);
        OPP_export_nonexec_list[set->index] = h_list;
    }

    /*-STEP 6 - Exchange execute set elements/data using the import/export lists--*/

    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];

        if (set->is_particle) continue; 

        halo_list i_list = OPP_import_exec_list[set->index];
        halo_list e_list = OPP_export_exec_list[set->index];

        int d = -1; // d is just simply the tag for mpi comms
        for (int k = 0; k < (int)opp_dats.size(); k++) // for each dat
        {
            d++; // increase tag to do mpi comm for the next opp_dat
            opp_dat dat = opp_dats[k];

            if (compare_sets(set, dat->set) == 1) // if this data array is defined on this set
            { 
                std::vector<MPI_Request> request_send(e_list->ranks_size);

                // prepare execute set element data to be exported
                char **sbuf = (char **)opp_host_malloc(e_list->ranks_size * sizeof(char *));

                for (int i = 0; i < e_list->ranks_size; i++) 
                {
                    sbuf[i] = (char *)opp_host_malloc((size_t)e_list->sizes[i] * (size_t)dat->size);
                    for (int j = 0; j < e_list->sizes[i]; j++) 
                    {
                        int set_elem_index = e_list->list[e_list->disps[i] + j];
                        memcpy(&sbuf[i][j * (size_t)dat->size],
                            (void *)&dat->data[(size_t)dat->size * (set_elem_index)], dat->size);
                    }

                    MPI_Isend(sbuf[i], (size_t)dat->size * e_list->sizes[i], MPI_CHAR,
                                e_list->ranks[i], d, OPP_MPI_WORLD, &request_send[i]);
                }

                // prepare space for the incomming data - realloc each
                // data array in each mpi process
                dat->data =
                    (char *)opp_host_realloc(dat->data, (size_t)(set->size + i_list->size) * (size_t)dat->size);

                int init = set->size * (size_t)dat->size;
                for (int i = 0; i < i_list->ranks_size; i++) 
                {
                    MPI_Recv(&(dat->data[init + i_list->disps[i] * (size_t)dat->size]),
                            (size_t)dat->size * i_list->sizes[i], MPI_CHAR, i_list->ranks[i], d,
                            OPP_MPI_WORLD, MPI_STATUS_IGNORE);
                }

                MPI_Waitall(e_list->ranks_size, request_send.data(), MPI_STATUSES_IGNORE);
                
                for (int i = 0; i < e_list->ranks_size; i++)
                    opp_host_free(sbuf[i]);
                opp_host_free(sbuf);
            }
        }
    }

    /*-STEP 7 - Exchange non-execute set elements/data using the import/export lists--*/

    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;  

        halo_list i_list = OPP_import_nonexec_list[set->index];
        halo_list e_list = OPP_export_nonexec_list[set->index];

        int d = -1; // d is just simply the tag for mpi comms

        for (int k = 0; k < (int)opp_dats.size(); k++) // for each dat
        {
            d++; // increase tag to do mpi comm for the next opp_dat
            opp_dat dat = opp_dats[k];

            if (compare_sets(set, dat->set) == 1)  // if this data array is defined on this set
            {
                // printf("on rank %d, The data array is %10s\n",my_rank,dat->name);
                std::vector<MPI_Request> request_send(e_list->ranks_size);

                // prepare non-execute set element data to be exported
                char **sbuf = (char **)opp_host_malloc(e_list->ranks_size * sizeof(char *));

                for (int i = 0; i < e_list->ranks_size; i++) 
                {
                    sbuf[i] = (char *)opp_host_malloc(e_list->sizes[i] * (size_t)dat->size);
                    for (int j = 0; j < e_list->sizes[i]; j++) 
                    {
                        int set_elem_index = e_list->list[e_list->disps[i] + j];
                        memcpy(&sbuf[i][j * (size_t)dat->size], 
                            (void *)&dat->data[(size_t)dat->size * (set_elem_index)], dat->size);
                    }
                    MPI_Isend(sbuf[i], (size_t)dat->size * e_list->sizes[i], MPI_CHAR,
                                e_list->ranks[i], d, OPP_MPI_WORLD, &request_send[i]);
                }

                // prepare space for the incomming nonexec-data - realloc each data array in each mpi process
                halo_list exec_i_list = OPP_import_exec_list[set->index];

                dat->data = (char *)opp_host_realloc(dat->data,
                    (size_t)(set->size + exec_i_list->size + i_list->size) * (size_t)dat->size);

                size_t init = (size_t)(set->size + exec_i_list->size) * (size_t)dat->size;
                for (int i = 0; i < i_list->ranks_size; i++) {
                MPI_Recv(&(dat->data[init + i_list->disps[i] * (size_t)dat->size]),
                        (size_t)dat->size * i_list->sizes[i], MPI_CHAR, i_list->ranks[i], d,
                        OPP_MPI_WORLD, MPI_STATUS_IGNORE);
                }

                MPI_Waitall(e_list->ranks_size, request_send.data(), MPI_STATUSES_IGNORE);
                for (int i = 0; i < e_list->ranks_size; i++)
                    opp_host_free(sbuf[i]);
                opp_host_free(sbuf);
            }
        }
    }

    /*-STEP 8 ----------------- Renumber Mapping tables-----------------------*/

    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        for (int m = 0; m < (int)opp_maps.size(); m++) // for each maping table
        { 
            opp_map map = opp_maps[m];

            if (compare_sets(map->to, set) == 1) { // need to select
                                                    // mappings TO this set

                halo_list exec_set_list = OPP_import_exec_list[set->index];
                halo_list nonexec_set_list = OPP_import_nonexec_list[set->index];

                halo_list exec_map_list = OPP_import_exec_list[map->from->index];

                // for each entry in this mapping table: original+execlist
                int len = map->from->size + exec_map_list->size;
                for (int e = 0; e < len; e++) 
                {
                    for (int j = 0; j < map->dim; j++) // for each element  pointed at by this entry
                    {     
                        int part;
                        int local_index = 0;
                        part = get_partition(map->map[e * map->dim + j],
                                            part_range[map->to->index], &local_index,
                                            comm_size);

                        if (part == my_rank) 
                        {
                            opp_maps[map->index]->map[e * map->dim + j] = local_index;
                        } 
                        else 
                        {
                            int found = -1;
                            // check in exec list
                            int rank1 = binary_search(exec_set_list->ranks, part, 0,
                                                        exec_set_list->ranks_size - 1);
                            // check in nonexec list
                            int rank2 = binary_search(nonexec_set_list->ranks, part, 0,
                                                        nonexec_set_list->ranks_size - 1);

                            if (rank1 >= 0) {
                                found = binary_search(exec_set_list->list, local_index,
                                                    exec_set_list->disps[rank1],
                                                    exec_set_list->disps[rank1] +
                                                        exec_set_list->sizes[rank1] - 1);
                                if (found >= 0) 
                                {
                                    opp_maps[map->index]->map[e * map->dim + j] = found + map->to->size;
                                }
                            }

                            if (rank2 >= 0 && found < 0) 
                            {
                                found = binary_search(nonexec_set_list->list, local_index,
                                                    nonexec_set_list->disps[rank2],
                                                    nonexec_set_list->disps[rank2] +
                                                        nonexec_set_list->sizes[rank2] - 1);
                                if (found >= 0) 
                                {
                                    opp_maps[map->index]->map[e * map->dim + j] = found + set->size + 
                                                                                        exec_set_list->size;
                                }
                            }

                            if (found < 0)
                            {
                                if (compare_sets(map->from, map->to) == 1 && (e >= map->from->size))
                                {
                                    opp_maps[map->index]->map[e * map->dim + j] = -2;
                                }
                                else
                                {
                                    printf("ERROR: Set %10s Element %d needed on rank %d from partition %d\n",
                                        set->name, local_index, my_rank, part);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /*-STEP 9 ---------------- Create MPI send Buffers-----------------------*/

    for (int k = 0; k < (int)opp_dats.size(); k++) // for each dat
    {
        opp_dat dat = opp_dats[k];

        if (dat->set->is_particle) continue;

        {
            op_mpi_buffer mpi_buf = (op_mpi_buffer)opp_host_malloc(sizeof(op_mpi_buffer_core));

            halo_list exec_e_list = OPP_export_exec_list[dat->set->index];
            halo_list nonexec_e_list = OPP_export_nonexec_list[dat->set->index];

            mpi_buf->buf_exec = (char *)opp_host_malloc((size_t)(exec_e_list->size) * (size_t)dat->size);
            mpi_buf->buf_nonexec = (char *)opp_host_malloc((size_t)(nonexec_e_list->size) * (size_t)dat->size);

            halo_list exec_i_list = OPP_import_exec_list[dat->set->index];
            halo_list nonexec_i_list = OPP_import_nonexec_list[dat->set->index];

            mpi_buf->s_req = (MPI_Request *)opp_host_malloc(
                sizeof(MPI_Request) * (exec_e_list->ranks_size + nonexec_e_list->ranks_size));
            mpi_buf->r_req = (MPI_Request *)opp_host_malloc(
                sizeof(MPI_Request) * (exec_i_list->ranks_size + nonexec_i_list->ranks_size));

            mpi_buf->s_num_req = 0;
            mpi_buf->r_num_req = 0;
            dat->mpi_buffer = mpi_buf;
        }

        // This is specifically for double indirected reduction operations
        // for now, sending only non exec halos, check whether we need exec halos too!
        {
            op_mpi_buffer mpi_buf = (op_mpi_buffer)opp_host_malloc(sizeof(op_mpi_buffer_core));

            // halo_list exec_i_list = OPP_import_exec_list[dat->set->index];
            halo_list nonexec_i_list = OPP_import_nonexec_list[dat->set->index];
            
            // halo_list exec_e_list = OPP_export_exec_list[dat->set->index];
            halo_list nonexec_e_list = OPP_export_nonexec_list[dat->set->index];

            int recv_buf_size = (int)(nonexec_e_list->size) * (int)dat->size;           
            mpi_buf->buf_nonexec = (char *)opp_host_malloc(recv_buf_size);

            if (OPP_DBG) opp_printf("opp_halo_create", "STEP 9 - dat [%s] mpi_reduc_buffer buf_nonexec size [%d]", 
                            dat->name, recv_buf_size);

            //(exec_i_list->ranks_size + nonexec_i_list->ranks_size));
            //(exec_e_list->ranks_size + nonexec_e_list->ranks_size));

            mpi_buf->s_req = (MPI_Request *)opp_host_malloc(sizeof(MPI_Request) * (nonexec_i_list->ranks_size)); 
            mpi_buf->r_req = (MPI_Request *)opp_host_malloc(sizeof(MPI_Request) * (nonexec_e_list->ranks_size)); 

            mpi_buf->s_num_req = 0;
            mpi_buf->r_num_req = 0;
            dat->mpi_reduc_buffer = mpi_buf;
        }
    }

    // set dirty bits of all data arrays to 0 for each data array
    for (int k = 0; k < (int)opp_dats.size(); k++) // for each dat
    {
        opp_dat dat = opp_dats[k];
        dat->dirtybit = 0;
    }

    /*-STEP 10 -------------------- Separate core elements------------------------*/

    int **core_elems = (int **)opp_host_malloc((int)opp_sets.size() * sizeof(int *));
    int **exp_elems = (int **)opp_host_malloc((int)opp_sets.size() * sizeof(int *));

    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        halo_list exec = OPP_export_exec_list[set->index];
        halo_list nonexec = OPP_export_nonexec_list[set->index];

        if (exec->size > 0) 
        {
            exp_elems[set->index] = (int *)opp_host_malloc(exec->size * sizeof(int));
            memcpy(exp_elems[set->index], exec->list, exec->size * sizeof(int));
            quickSort(exp_elems[set->index], 0, exec->size - 1);

            int num_exp = removeDups(exp_elems[set->index], exec->size);
            core_elems[set->index] = (int *)opp_host_malloc(set->size * sizeof(int));
            int count = 0;
            for (int e = 0; e < set->size; e++) // for each elment of this set
            {
                if ((binary_search(exp_elems[set->index], e, 0, num_exp - 1) < 0)) 
                {
                    core_elems[set->index][count++] = e;
                }
            }
            quickSort(core_elems[set->index], 0, count - 1);

            if (count + num_exp != set->size)
                printf("sizes not equal\n");
            set->core_size = count;

            // for each data array defined on this set seperate its elements
            for (int k = 0; k < (int)opp_dats.size(); k++) // for each dat
            {
                opp_dat dat = opp_dats[k];

                if (compare_sets(set, dat->set) == 1) // if this data array is defined on this set
                {
                    char *new_dat = (char *)opp_host_malloc((size_t)set->size * (size_t)dat->size);
                    for (int i = 0; i < count; i++) {
                        memcpy(&new_dat[i * (size_t)dat->size], 
                            &dat->data[core_elems[set->index][i] * (size_t)dat->size], dat->size);
                    }
                    for (int i = 0; i < num_exp; i++) {
                        memcpy(&new_dat[(count + i) * (size_t)dat->size], 
                            &dat->data[exp_elems[set->index][i] * (size_t)dat->size], dat->size);
                    }
                    memcpy(&dat->data[0], &new_dat[0], set->size * (size_t)dat->size);
                    opp_host_free(new_dat);
                }
            }

            // for each mapping defined from this set seperate its elements
            for (int m = 0; m < (int)opp_maps.size(); m++) // for each set
            { 
                opp_map map = opp_maps[m];

                if (compare_sets(map->from, set) == 1)  // if this mapping is defined from this set
                { 
                    int *new_map = (int *)opp_host_malloc((size_t)set->size * map->dim * sizeof(int));
                    for (int i = 0; i < count; i++) {
                        memcpy(&new_map[i * (size_t)map->dim],
                            &map->map[core_elems[set->index][i] * (size_t)map->dim],
                            map->dim * sizeof(int));
                    }
                    for (int i = 0; i < num_exp; i++) 
                    {
                        memcpy(&new_map[(count + i) * (size_t)map->dim],
                            &map->map[exp_elems[set->index][i] * (size_t)map->dim],
                            map->dim * sizeof(int));
                    }
                    memcpy(&map->map[0], &new_map[0], set->size * (size_t)map->dim * sizeof(int));
                    opp_host_free(new_map);
                }
            }

            for (int i = 0; i < exec->size; i++) 
            {
                int index = binary_search(exp_elems[set->index], exec->list[i], 0, num_exp - 1);
                if (index < 0)
                    printf("Problem in seperating core elements - exec list\n");
                else
                    exec->list[i] = count + index;
            }

            for (int i = 0; i < nonexec->size; i++) 
            {
                int index = binary_search(core_elems[set->index], nonexec->list[i], 0, count - 1);
                if (index < 0) 
                {
                    index = binary_search(exp_elems[set->index], nonexec->list[i], 0, num_exp - 1);
                    if (index < 0)
                        printf("Problem in seperating core elements - nonexec list\n");
                    else
                        nonexec->list[i] = count + index;
                } 
                else
                    nonexec->list[i] = index;
            }
        } 
        else 
        {
            core_elems[set->index] = (int *)opp_host_malloc(set->size * sizeof(int));
            exp_elems[set->index] = (int *)opp_host_malloc(0 * sizeof(int));
            for (int e = 0; e < set->size; e++)  // for each elment of this set
            {
                core_elems[set->index][e] = e;
            }
            set->core_size = set->size;
        }
    }

    // now need to renumber mapping tables as the elements are seperated
    for (int m = 0; m < (int)opp_maps.size(); m++) { // for each set
        opp_map map = opp_maps[m];

        halo_list exec_map_list = OPP_import_exec_list[map->from->index];
        // for each entry in this mapping table: original+execlist
        int len = map->from->size + exec_map_list->size;
        for (int e = 0; e < len; e++) 
        {
            for (int j = 0; j < map->dim; j++) // for each element pointed  at by this entry 
            { 
                if (map->map[e * map->dim + j] < map->to->size) 
                {
                    int index = binary_search(core_elems[map->to->index],
                                                map->map[e * map->dim + j], 0,
                                                map->to->core_size - 1);
                    if (index < 0) 
                    {
                        index = binary_search(exp_elems[map->to->index], map->map[e * map->dim + j], 0,
                                            (map->to->size) - (map->to->core_size) - 1);
                        if (index < 0)
                        {
                            if (compare_sets(map->from, map->to) == 1)
                                opp_maps[map->index]->map[e * (size_t)map->dim + j] = -3;
                            else
                                printf("Problem in seperating core elements - renumbering map\n");
                        }
                        else
                            opp_maps[map->index]->map[e * (size_t)map->dim + j] = map->to->core_size + index;
                    } 
                    else
                        opp_maps[map->index]->map[e * (size_t)map->dim + j] = index;
                }
            }
        }
    }

    /*-STEP 11 ----------- Save the original set element indexes------------------*/

    // if OPP_part_list is empty, (i.e. no previous partitioning done) then create it and store the seperation of 
    // elements using core_elems and exp_elems
    if (OPP_part_index != (int)opp_sets.size()) 
    {
        // allocate memory for list
        OPP_part_list = (part *)opp_host_malloc((int)opp_sets.size() * sizeof(part));

        for (int s = 0; s < (int)opp_sets.size(); s++)  // for each set
        {
            opp_set set = opp_sets[s];

            if (set->is_particle) continue;

            int *g_index = (int *)opp_host_malloc(sizeof(int) * set->size);
            int *partition = (int *)opp_host_malloc(sizeof(int) * set->size);
            for (int i = 0; i < set->size; i++) 
            {
                g_index[i] = get_global_index(i, my_rank, part_range[set->index], comm_size);
                partition[i] = my_rank;
            }
            decl_partition(set, g_index, partition);

            // combine core_elems and exp_elems to one memory block
            int *temp = (int *)opp_host_malloc(sizeof(int) * set->size);
            memcpy(&temp[0], core_elems[set->index], set->core_size * sizeof(int));
            memcpy(&temp[set->core_size], exp_elems[set->index],
                    (set->size - set->core_size) * sizeof(int));

            // update OPP_part_list[set->index]->g_index
            for (int i = 0; i < set->size; i++) 
            {
                temp[i] = OPP_part_list[set->index]->g_index[temp[i]];
            }
            opp_host_free(OPP_part_list[set->index]->g_index);
            OPP_part_list[set->index]->g_index = temp;
        }
    } 
    else  // OPP_part_list exists (i.e. a partitioning has been done) update the seperation of elements
    {        

        for (int s = 0; s < (int)opp_sets.size(); s++)  // for each set
        { 
            opp_set set = opp_sets[s];

            if (set->is_particle) continue;

            // combine core_elems and exp_elems to one memory block
            int *temp = (int *)opp_host_malloc(sizeof(int) * set->size);
            memcpy(&temp[0], core_elems[set->index], set->core_size * sizeof(int));
            memcpy(&temp[set->core_size], exp_elems[set->index], (set->size - set->core_size) * sizeof(int));

            // update OPP_part_list[set->index]->g_index
            for (int i = 0; i < set->size; i++) 
            {
                temp[i] = OPP_part_list[set->index]->g_index[temp[i]];
            }
            opp_host_free(OPP_part_list[set->index]->g_index);
            OPP_part_list[set->index]->g_index = temp;
        }
    }

    // set up exec and nonexec sizes
    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        set->exec_size = OPP_import_exec_list[set->index]->size;
        set->nonexec_size = OPP_import_nonexec_list[set->index]->size;

        set->set_capacity = set->size + set->exec_size + set->nonexec_size;
    }

    /*-STEP 12 ---------- Clean up and Compute rough halo size numbers------------*/

    for (int i = 0; i < (int)opp_sets.size(); i++) 
    {
        opp_host_free(part_range[i]);

        if (opp_sets[i]->is_particle) continue;
        
        opp_host_free(core_elems[i]);
        opp_host_free(exp_elems[i]);
    }
    opp_host_free(part_range);
    opp_host_free(exp_elems);
    opp_host_free(core_elems);

    op_timers(&cpu_t2, &wall_t2); // timer stop for list create
    // compute import/export lists creation time
    time = wall_t2 - wall_t1;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, OPP_ROOT, OPP_MPI_WORLD);

    // compute avg/min/max set sizes and exec sizes accross the MPI universe
    int avg_size = 0, min_size = 0, max_size = 0;
    for (int s = 0; s < (int)opp_sets.size(); s++) 
    {
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        // number of set elements first
        MPI_Reduce(&set->size, &avg_size, 1, MPI_INT, MPI_SUM, OPP_ROOT,
                OPP_MPI_WORLD);
        MPI_Reduce(&set->size, &min_size, 1, MPI_INT, MPI_MIN, OPP_ROOT,
                OPP_MPI_WORLD);
        MPI_Reduce(&set->size, &max_size, 1, MPI_INT, MPI_MAX, OPP_ROOT,
                OPP_MPI_WORLD);

        if (OPP_DBG && my_rank == OPP_ROOT) 
        {
            opp_printf("opp_halo_create", "Num of %8s (avg | min | max)", set->name);
            opp_printf("opp_halo_create", "total elems         %10d %10d %10d", 
                avg_size / comm_size, min_size, max_size);
        }

        avg_size = 0;
        min_size = 0;
        max_size = 0;

        // number of OWNED elements second
        MPI_Reduce(&set->core_size, &avg_size, 1, MPI_INT, MPI_SUM, OPP_ROOT,
                OPP_MPI_WORLD);
        MPI_Reduce(&set->core_size, &min_size, 1, MPI_INT, MPI_MIN, OPP_ROOT,
                OPP_MPI_WORLD);
        MPI_Reduce(&set->core_size, &max_size, 1, MPI_INT, MPI_MAX, OPP_ROOT,
                OPP_MPI_WORLD);

        if (OPP_DBG && my_rank == OPP_ROOT) 
        {
            opp_printf("opp_halo_create", "core elems         %10d %10d %10d", 
                avg_size / comm_size, min_size, max_size);
        }
        avg_size = 0;
        min_size = 0;
        max_size = 0;

        // number of exec halo elements third
        MPI_Reduce(&OPP_import_exec_list[set->index]->size, &avg_size, 1, MPI_INT,
                MPI_SUM, OPP_ROOT, OPP_MPI_WORLD);
        MPI_Reduce(&OPP_import_exec_list[set->index]->size, &min_size, 1, MPI_INT,
                MPI_MIN, OPP_ROOT, OPP_MPI_WORLD);
        MPI_Reduce(&OPP_import_exec_list[set->index]->size, &max_size, 1, MPI_INT,
                MPI_MAX, OPP_ROOT, OPP_MPI_WORLD);

        if (OPP_DBG && my_rank == OPP_ROOT) 
        {
            opp_printf("opp_halo_create", "exec halo elems     %10d %10d %10d", 
                avg_size / comm_size, min_size, max_size);
        }
        avg_size = 0;
        min_size = 0;
        max_size = 0;

        // number of non-exec halo elements fourth
        MPI_Reduce(&OPP_import_nonexec_list[set->index]->size, &avg_size, 1, MPI_INT,
                MPI_SUM, OPP_ROOT, OPP_MPI_WORLD);
        MPI_Reduce(&OPP_import_nonexec_list[set->index]->size, &min_size, 1, MPI_INT,
                MPI_MIN, OPP_ROOT, OPP_MPI_WORLD);
        MPI_Reduce(&OPP_import_nonexec_list[set->index]->size, &max_size, 1, MPI_INT,
                MPI_MAX, OPP_ROOT, OPP_MPI_WORLD);

        if (OPP_DBG && my_rank == OPP_ROOT) 
        {
            opp_printf("opp_halo_create", "non-exec halo elems %10d %10d %10d", 
                avg_size / comm_size, min_size, max_size);
        }
        avg_size = 0;
        min_size = 0;
        max_size = 0;
        if (OPP_DBG && my_rank == OPP_ROOT) 
        {
            opp_printf("opp_halo_create", "-------------------------------------");
        }
    }

    // compute avg/min/max number of MPI neighbors per process accross the MPI
    // universe
    avg_size = 0, min_size = 0, max_size = 0;
    for (int s = 0; s < (int)opp_sets.size(); s++) 
    {
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        // number of exec halo neighbors first
        MPI_Reduce(&OPP_import_exec_list[set->index]->ranks_size, &avg_size, 1,
                MPI_INT, MPI_SUM, OPP_ROOT, OPP_MPI_WORLD);
        MPI_Reduce(&OPP_import_exec_list[set->index]->ranks_size, &min_size, 1,
                MPI_INT, MPI_MIN, OPP_ROOT, OPP_MPI_WORLD);
        MPI_Reduce(&OPP_import_exec_list[set->index]->ranks_size, &max_size, 1,
                MPI_INT, MPI_MAX, OPP_ROOT, OPP_MPI_WORLD);

        if (OPP_DBG && my_rank == OPP_ROOT) 
        {
            opp_printf("opp_halo_create", "MPI neighbors for exchanging %8s (avg | min | max)", 
                set->name);
            opp_printf("opp_halo_create", "exec halo elems     %4d %4d %4d", 
                avg_size / comm_size, min_size, max_size);
        }
        avg_size = 0;
        min_size = 0;
        max_size = 0;

        // number of non-exec halo neighbors second
        MPI_Reduce(&OPP_import_nonexec_list[set->index]->ranks_size, &avg_size, 1,
                MPI_INT, MPI_SUM, OPP_ROOT, OPP_MPI_WORLD);
        MPI_Reduce(&OPP_import_nonexec_list[set->index]->ranks_size, &min_size, 1,
                MPI_INT, MPI_MIN, OPP_ROOT, OPP_MPI_WORLD);
        MPI_Reduce(&OPP_import_nonexec_list[set->index]->ranks_size, &max_size, 1,
                MPI_INT, MPI_MAX, OPP_ROOT, OPP_MPI_WORLD);

        if (OPP_DBG && my_rank == OPP_ROOT) 
        {
            opp_printf("opp_halo_create", "non-exec halo elems %4d %4d %4d", 
                avg_size / comm_size, min_size, max_size);
        }
        avg_size = 0;
        min_size = 0;
        max_size = 0;
        if (OPP_DBG && my_rank == OPP_ROOT) {
            opp_printf("opp_halo_create", "-------------------------------------");
        }
    }

    // compute average worst case halo size in Bytes
    int tot_halo_size = 0;
    for (int s = 0; s < (int)opp_sets.size(); s++) 
    {
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        for (int k = 0; k < (int)opp_dats.size(); k++) // for each dat
        {
            opp_dat dat = opp_dats[k];

            if (compare_sets(dat->set, set) == 1) 
            {
                halo_list exec_imp = OPP_import_exec_list[set->index];
                halo_list nonexec_imp = OPP_import_nonexec_list[set->index];
                tot_halo_size = tot_halo_size + exec_imp->size * (size_t)dat->size +
                                nonexec_imp->size * (size_t)dat->size;
            }
        }
    }
    int avg_halo_size;
    MPI_Reduce(&tot_halo_size, &avg_halo_size, 1, MPI_INT, MPI_SUM, OPP_ROOT, OPP_MPI_WORLD);

    // print performance results
    if (OPP_DBG && my_rank == OPP_ROOT) 
    {
        opp_printf("opp_halo_create", "Max total halo creation time = %lf", max_time);
        opp_printf("opp_halo_create", "Average (worst case) Halo size = %d Bytes", 
            avg_halo_size / comm_size);
    }

    if (OPP_DBG) opp_printf("opp_halo_create", "end");
}

/*******************************************************************************
 * Routine to Clean-up all MPI halos(called at the end of an OP2 MPI application)
*******************************************************************************/
void __opp_halo_destroy() 
{
    // remove halos from opp_dats
    for (int k = 0; k < (int)opp_dats.size(); k++) // for each dat
    {
        opp_dat dat = opp_dats[k];
        dat->data = (char *)opp_host_realloc(dat->data, (size_t)dat->set->size * dat->size);
    }

    //free lists
    for (int s = 0; s < (int)opp_sets.size(); s++) 
    {
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        opp_host_free(OPP_import_exec_list[set->index]->ranks);
        opp_host_free(OPP_import_exec_list[set->index]->disps);
        opp_host_free(OPP_import_exec_list[set->index]->sizes);
        opp_host_free(OPP_import_exec_list[set->index]->list);
        opp_host_free(OPP_import_exec_list[set->index]);

        opp_host_free(OPP_import_nonexec_list[set->index]->ranks);
        opp_host_free(OPP_import_nonexec_list[set->index]->disps);
        opp_host_free(OPP_import_nonexec_list[set->index]->sizes);
        opp_host_free(OPP_import_nonexec_list[set->index]->list);
        opp_host_free(OPP_import_nonexec_list[set->index]);

        opp_host_free(OPP_export_exec_list[set->index]->ranks);
        opp_host_free(OPP_export_exec_list[set->index]->disps);
        opp_host_free(OPP_export_exec_list[set->index]->sizes);
        opp_host_free(OPP_export_exec_list[set->index]->list);
        opp_host_free(OPP_export_exec_list[set->index]);

        opp_host_free(OPP_export_nonexec_list[set->index]->ranks);
        opp_host_free(OPP_export_nonexec_list[set->index]->disps);
        opp_host_free(OPP_export_nonexec_list[set->index]->sizes);
        opp_host_free(OPP_export_nonexec_list[set->index]->list);
        opp_host_free(OPP_export_nonexec_list[set->index]);
    }
    opp_host_free(OPP_import_exec_list);
    opp_host_free(OPP_import_nonexec_list);
    opp_host_free(OPP_export_exec_list);
    opp_host_free(OPP_export_nonexec_list);

    for (int k = 0; k < (int)opp_dats.size(); k++) // for each dat
    {
        opp_dat dat = opp_dats[k];

        if (dat->set->is_particle) continue;

        opp_host_free(((op_mpi_buffer)(dat->mpi_buffer))->buf_exec);
        opp_host_free(((op_mpi_buffer)(dat->mpi_buffer))->buf_nonexec);
        opp_host_free(((op_mpi_buffer)(dat->mpi_buffer))->s_req);
        opp_host_free(((op_mpi_buffer)(dat->mpi_buffer))->r_req);
        opp_host_free((op_mpi_buffer)(dat->mpi_buffer));

        //opp_host_free(((op_mpi_buffer)(dat->mpi_reduc_buffer))->buf_exec);
        opp_host_free(((op_mpi_buffer)(dat->mpi_reduc_buffer))->buf_nonexec);
        opp_host_free(((op_mpi_buffer)(dat->mpi_reduc_buffer))->s_req);
        opp_host_free(((op_mpi_buffer)(dat->mpi_reduc_buffer))->r_req);
        opp_host_free((op_mpi_buffer)(dat->mpi_reduc_buffer));
    }
}

/*******************************************************************************
 * Routine to exchange MPI halos of the arg using the host data structures
*******************************************************************************/
void __opp_mpi_host_halo_exchange(opp_arg *arg, int exec_flag)
{
    opp_dat dat = arg->dat;

    if (arg->sent == 1) 
    {
        opp_printf("opp_mpi_halo_exchange", "Error: Halo exchange already in flight for dat %s", dat->name);
        fflush(stdout);
        opp_abort("opp_mpi_halo_exchange");
    }

    arg->sent = 0; // reset flag

    if (OPP_DBG) 
        opp_printf("opp_mpi_halo_exchange", "Exchanging Halo of data array [%s]", dat->name);

    halo_list imp_exec_list    = OPP_import_exec_list[dat->set->index];
    halo_list imp_nonexec_list = OPP_import_nonexec_list[dat->set->index];

    halo_list exp_exec_list    = OPP_export_exec_list[dat->set->index];
    halo_list exp_nonexec_list = OPP_export_nonexec_list[dat->set->index];

    //-------first exchange exec elements related to this data array--------

    // sanity checks
    if (compare_sets(imp_exec_list->set, dat->set) == 0) 
    {
        opp_printf("opp_mpi_halo_exchange", "Error: Import list and set mismatch\n");
        opp_abort("opp_mpi_halo_exchange Error: Import list and set mismatch");
    }
    if (compare_sets(exp_exec_list->set, dat->set) == 0) 
    {
        opp_printf("opp_mpi_halo_exchange", "Error: Export list and set mismatch\n");
        opp_abort("opp_mpi_halo_exchange Error: Export list and set mismatch");
    }

    op_mpi_buffer mpi_buffer = (op_mpi_buffer)(dat->mpi_buffer);

    opp_profiler->startMpiComm("", opp::OPP_Mesh);

    double total_send_size = 0.0;

    int set_elem_index = -1;
    for (int i = 0; i < exp_exec_list->ranks_size; i++) 
    {
        for (int j = 0; j < exp_exec_list->sizes[i]; j++) 
        {
            set_elem_index = exp_exec_list->list[exp_exec_list->disps[i] + j];
            memcpy(&mpi_buffer->buf_exec[exp_exec_list->disps[i] * dat->size + j * dat->size],
                (void *)&dat->data[dat->size * (set_elem_index)], dat->size);
        }
        
        MPI_Isend(&mpi_buffer->buf_exec[exp_exec_list->disps[i] * dat->size], 
            dat->size * exp_exec_list->sizes[i], MPI_CHAR, exp_exec_list->ranks[i], dat->index, 
            OPP_MPI_WORLD, &mpi_buffer->s_req[mpi_buffer->s_num_req++]);
        
        total_send_size += (dat->size * exp_exec_list->sizes[i] * 1.0f);
    }

    int init = dat->set->size * dat->size;
    for (int i = 0; i < imp_exec_list->ranks_size; i++) 
    {
        MPI_Irecv(&(dat->data[init + imp_exec_list->disps[i] * dat->size]), 
            dat->size * imp_exec_list->sizes[i], MPI_CHAR, imp_exec_list->ranks[i], dat->index, 
            OPP_MPI_WORLD, &mpi_buffer->r_req[mpi_buffer->r_num_req++]);
    }

    //-----second exchange nonexec elements related to this data array------
    
    // sanity checks
    if (compare_sets(imp_nonexec_list->set, dat->set) == 0) 
    {
        opp_printf("opp_mpi_halo_exchange", "Error: Non-Import list and set mismatch");
        opp_abort("opp_mpi_halo_exchange Error: Non-Import list and set mismatch");
    }
    if (compare_sets(exp_nonexec_list->set, dat->set) == 0) 
    {
        opp_printf("opp_mpi_halo_exchange", "Error: Non-Export list and set mismatch");
        opp_abort("opp_mpi_halo_exchange Error: Non-Export list and set mismatch");
    }

    for (int i = 0; i < exp_nonexec_list->ranks_size; i++) 
    {
        for (int j = 0; j < exp_nonexec_list->sizes[i]; j++) 
        {
            set_elem_index = exp_nonexec_list->list[exp_nonexec_list->disps[i] + j];
            memcpy(&mpi_buffer->buf_nonexec[exp_nonexec_list->disps[i] * dat->size + j * dat->size],
                (void *)&dat->data[dat->size * (set_elem_index)], dat->size);
        }

        MPI_Isend(&mpi_buffer->buf_nonexec[exp_nonexec_list->disps[i] * dat->size], 
            dat->size * exp_nonexec_list->sizes[i], MPI_CHAR, exp_nonexec_list->ranks[i], dat->index, 
            OPP_MPI_WORLD, &mpi_buffer->s_req[mpi_buffer->s_num_req++]);
        
        total_send_size += (dat->size * exp_nonexec_list->sizes[i] * 1.0f);
    }

    int nonexec_init = (dat->set->size + imp_exec_list->size) * dat->size;
    for (int i = 0; i < imp_nonexec_list->ranks_size; i++) 
    {
        MPI_Irecv(&(dat->data[nonexec_init + imp_nonexec_list->disps[i] * dat->size]), 
            dat->size * imp_nonexec_list->sizes[i], MPI_CHAR, imp_nonexec_list->ranks[i], dat->index, 
            OPP_MPI_WORLD, &mpi_buffer->r_req[mpi_buffer->r_num_req++]);
    }
    
    opp_profiler->addTransferSize("", opp::OPP_Mesh, total_send_size, 1);

    opp_profiler->endMpiComm("", opp::OPP_Mesh);

    // clear dirty bit
    dat->dirtybit = 0;
    arg->sent = 1;
}

/*******************************************************************************
 * Routine to wait and complete all MPI halo exchanges done with host data structures
*******************************************************************************/
void __opp_mpi_host_halo_wait_all(int nargs, opp_arg *args)
{
    if (OPP_DBG) opp_printf("__opp_mpi_host_halo_wait_all", "START");

    opp_profiler->startMpiComm("", opp::OPP_Mesh);

    for (int n = 0; n < nargs; n++) 
    {
        opp_arg *arg = &args[n];

        if (arg->opt && arg->argtype == OPP_ARG_DAT && arg->sent == 1) 
        {
            opp_dat dat = arg->dat;

            op_mpi_buffer mpi_buffer = (op_mpi_buffer)(dat->mpi_buffer);

            MPI_Waitall(mpi_buffer->s_num_req, mpi_buffer->s_req, MPI_STATUSES_IGNORE);
            MPI_Waitall(mpi_buffer->r_num_req, mpi_buffer->r_req, MPI_STATUSES_IGNORE);

            // The halos are arranged in order, after the set size -- import exec halo -- import non exec halo
            // hence rearranging is not required

            mpi_buffer->s_num_req = 0;
            mpi_buffer->r_num_req = 0;

            arg->sent = 2; // set flag to indicate completed comm
        }
    }

    opp_profiler->endMpiComm("", opp::OPP_Mesh);

    if (OPP_DBG) opp_printf("__opp_mpi_host_halo_wait_all", "END");
}



// Utility functions

/*******************************************************************************
 * Routine to create a generic halo list
 * (used in both import and export list creation)
 *******************************************************************************/
void create_list(int *list, int *ranks, int *disps, int *sizes, int *ranks_size, int *total, 
        int *temp_list, int size, int comm_size, int my_rank) 
{
    (void)my_rank;
    int index = 0;
    int total_size = 0;
    if (size < 0) printf("problem\n");

    // negative values set as an initialisation
    for (int r = 0; r < comm_size; r++) 
    {
        disps[r] = ranks[r] = -99;
        sizes[r] = 0;
    }

    for (int r = 0; r < comm_size; r++) 
    {
        sizes[index] = disps[index] = 0;

        int *temp = (int *)opp_host_malloc((size / 2) * sizeof(int));
        for (int i = 0; i < size; i = i + 2) 
        {
            if (temp_list[i] == r)
                temp[sizes[index]++] = temp_list[i + 1];
        }
        if (sizes[index] > 0) 
        {
            ranks[index] = r;
            // sort temp,
            quickSort(temp, 0, sizes[index] - 1);
            // eliminate duplicates in temp
            sizes[index] = removeDups(temp, sizes[index]);
            total_size = total_size + sizes[index];

            if (index > 0)
                disps[index] = disps[index - 1] + sizes[index - 1];
            // add to end of exp_list
            for (int e = 0; e < sizes[index]; e++)
                list[disps[index] + e] = temp[e];

            index++;
        }
        opp_host_free(temp);
    }

    *total = total_size;
    *ranks_size = index;
}

/*******************************************************************************
 * Routine to create an export list
 *******************************************************************************/
void create_export_list(opp_set set, int *temp_list, halo_list h_list, int size, int comm_size, int my_rank) 
{
    int *ranks = (int *)opp_host_malloc(comm_size * sizeof(int));
    int *list = (int *)opp_host_malloc((size / 2) * sizeof(int));
    int *disps = (int *)opp_host_malloc(comm_size * sizeof(int));
    int *sizes = (int *)opp_host_malloc(comm_size * sizeof(int));

    int ranks_size = 0;
    int total_size = 0;

    create_list(list, ranks, disps, sizes, &ranks_size, &total_size, temp_list, size, comm_size, my_rank);

    h_list->set = set;
    h_list->size = total_size;
    h_list->ranks = ranks;
    h_list->ranks_size = ranks_size;
    h_list->disps = disps;
    h_list->sizes = sizes;
    h_list->list = list;
}

/*******************************************************************************
 * Routine to create an import list
 *******************************************************************************/
void create_import_list(opp_set set, int *temp_list, halo_list h_list, int total_size, int *ranks, int *sizes, 
                        int ranks_size, int comm_size, int my_rank) 
{
    (void)my_rank;
    int *disps = (int *)opp_host_malloc(comm_size * sizeof(int));
    disps[0] = 0;
    for (int i = 0; i < ranks_size; i++) 
    {
        if (i > 0)
            disps[i] = disps[i - 1] + sizes[i - 1];
    }

    h_list->set = set;
    h_list->size = total_size;
    h_list->ranks = ranks;
    h_list->ranks_size = ranks_size;
    h_list->disps = disps;
    h_list->sizes = sizes;
    h_list->list = temp_list;
}

/*******************************************************************************
 * Routine to create an nonexec-import list (only a wrapper)
 *******************************************************************************/
void create_nonexec_import_list(opp_set set, int *temp_list, halo_list h_list, int size, 
    int comm_size, int my_rank) 
{
    create_export_list(set, temp_list, h_list, size, comm_size, my_rank);
}

/*******************************************************************************
 * Routine to create an nonexec-export list (only a wrapper)
 *******************************************************************************/
void create_nonexec_export_list(opp_set set, int *temp_list, halo_list h_list, int total_size,
                                    int *ranks, int *sizes, int ranks_size, int comm_size, int my_rank) 
{
    create_import_list(set, temp_list, h_list, total_size, ranks, sizes, ranks_size, comm_size, my_rank);
}
/*
* Open source copyright declaration based on BSD open source template:
* http://www.opensource.org/licenses/bsd-license.php
*
* This file is part of the OP2 distribution.
*
* Copyright (c) 2011, Mike Giles and others. Please see the AUTHORS file in
* the main source directory for a full list of copyright holders.
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * The name of Mike Giles may not be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY Mike Giles ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL Mike Giles BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
* op_mpi_core.c
*
* Implements the OP2 Distributed memory (MPI) halo creation, halo exchange and
* support utility routines/functions
*
* written by: Gihan R. Mudalige, (Started 01-03-2011)
*/

// mpi header
#include <mpi.h>
#include <opp_mpi.h>

// MPI Halo related global variables

halo_list *OP_export_exec_list; // EEH list
halo_list *OP_import_exec_list; // IEH list

halo_list *OP_import_nonexec_list; // INH list
halo_list *OP_export_nonexec_list; // ENH list

// global variables to hold partition information on an MPI rank
int OP_part_index = 0;
part *OP_part_list;

// Save original partition ranges
int **orig_part_range = NULL;

// Sliding planes data structures
int OP_import_index = 0, OP_import_max = 0;
int OP_export_index = 0, OP_export_max = 0;
op_import_handle *OP_import_list = NULL;
op_export_handle *OP_export_list = NULL;

// Timing
double t1, t2, c1, c2;


/*******************************************************************************
 * Routine to declare partition information for a given set
 *******************************************************************************/
void decl_partition(op_set set, int *g_index, int *partition) 
{
    part p = (part)malloc(sizeof(part_core));
    p->set = set;
    p->g_index = g_index;
    p->elem_part = partition;
    p->is_partitioned = 0;
    OP_part_list[set->index] = p;
    OP_part_index++;
}

/*******************************************************************************
 * Routine to get partition range on all mpi ranks for all sets
 *******************************************************************************/
void get_part_range(int **part_range, int my_rank, int comm_size, MPI_Comm Comm) 
{
    (void)my_rank;
    for (int s = 0; s < OP_set_index; s++) 
    {
        op_set set = OP_set_list[s];

        int *sizes = (int *)malloc(sizeof(int) * comm_size);
        MPI_Allgather(&set->size, 1, MPI_INT, sizes, 1, MPI_INT, Comm);

        part_range[set->index] = (int *)malloc(2 * comm_size * sizeof(int));

        int disp = 0;
        for (int i = 0; i < comm_size; i++) 
        {
            part_range[set->index][2 * i] = disp;
            disp = disp + sizes[i] - 1;
            part_range[set->index][2 * i + 1] = disp;
            disp++;
        #ifdef DEBUG
            if (my_rank == OPP_MPI_ROOT && OP_diags > 5)
                printf("range of %10s in rank %d: %d-%d\n", set->name, i,
                    part_range[set->index][2 * i],
                    part_range[set->index][2 * i + 1]);
        #endif
        }
        free(sizes);
    }
}

/*******************************************************************************
 * Routine to get partition (i.e. mpi rank) where global_index is located and
 * its local index
 *******************************************************************************/
int get_partition(int global_index, int *part_range, int *local_index, int comm_size) 
{
    for (int i = 0; i < comm_size; i++) 
    {
        if (global_index >= part_range[2 * i] && global_index <= part_range[2 * i + 1]) 
        {
            *local_index = global_index - part_range[2 * i];
            return i;
        }
    }

    if (OP_DEBUG) 
    {    
        std::string log = "";
        for (int i = 0; i < comm_size; i++) 
        {
            log += std::to_string(i) + "|" + std::to_string(part_range[2 * i]) + "|" + std::to_string(part_range[2 * i + 1]) + "\n";
        }

        opp_printf("get_partition()", OPP_my_rank, "Error: orphan global index %d part_range->\n%s", global_index, log.c_str());
    }

    MPI_Abort(OP_MPI_WORLD, 2);
    return 1;
}

/*******************************************************************************
 * Routine to convert a local index in to a global index
 *******************************************************************************/
int get_global_index(int local_index, int partition, int *part_range, int comm_size) 
{
    (void)comm_size;
    int g_index = part_range[2 * partition] + local_index;
#ifdef DEBUG
    if (g_index > part_range[2 * (comm_size - 1) + 1] && OP_diags > 2)
        printf("Global index larger than set size\n");
#endif
    return g_index;
}

/*******************************************************************************
 * Routine to find the MPI neighbors given a halo list
 *******************************************************************************/
void find_neighbors_set(halo_list List, int *neighbors, int *sizes, int *ranks_size, int my_rank, int comm_size, MPI_Comm Comm) 
{
    int *temp = (int *)malloc(comm_size * sizeof(int));
    int *r_temp = (int *)malloc(comm_size * comm_size * sizeof(int));

    for (int r = 0; r < comm_size * comm_size; r++)
        r_temp[r] = -99;
    for (int r = 0; r < comm_size; r++)
        temp[r] = -99;

    int n = 0;

    for (int r = 0; r < comm_size; r++) 
    {
        if (List->ranks[r] >= 0)
        temp[List->ranks[r]] = List->sizes[r];
    }

    MPI_Allgather(temp, comm_size, MPI_INT, r_temp, comm_size, MPI_INT, Comm);

    for (int i = 0; i < comm_size; i++) 
    {
        if (i != my_rank)
        {
            if (r_temp[i * comm_size + my_rank] > 0) 
            {
                neighbors[n] = i;
                sizes[n] = r_temp[i * comm_size + my_rank];
                n++;
            }
        }
    }
    *ranks_size = n;
    free(temp);
    free(r_temp);
}

/*******************************************************************************
 * Routine to create a generic halo list
 * (used in both import and export list creation)
 *******************************************************************************/
void create_list(int *list, int *ranks, int *disps, int *sizes, int *ranks_size, int *total, int *temp_list, 
        int size, int comm_size, int my_rank) 
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

        int *temp = (int *)malloc((size / 2) * sizeof(int));
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
        free(temp);
    }

    *total = total_size;
    *ranks_size = index;
}

/*******************************************************************************
 * Routine to create an export list
 *******************************************************************************/
void create_export_list(op_set set, int *temp_list, halo_list h_list, int size, int comm_size, int my_rank) 
{
    int *ranks = (int *)malloc(comm_size * sizeof(int));
    int *list = (int *)malloc((size / 2) * sizeof(int));
    int *disps = (int *)malloc(comm_size * sizeof(int));
    int *sizes = (int *)malloc(comm_size * sizeof(int));

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
void create_import_list(op_set set, int *temp_list, halo_list h_list, int total_size, int *ranks, int *sizes, 
                        int ranks_size, int comm_size, int my_rank) 
{
    (void)my_rank;
    int *disps = (int *)malloc(comm_size * sizeof(int));
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
/* static */  void create_nonexec_import_list(op_set set, int *temp_list, halo_list h_list, int size, int comm_size, int my_rank) 
{
    create_export_list(set, temp_list, h_list, size, comm_size, my_rank);
}

/*******************************************************************************
 * Routine to create an nonexec-export list (only a wrapper)
 *******************************************************************************/
/* static */  void create_nonexec_export_list(op_set set, int *temp_list, halo_list h_list, int total_size,
                                    int *ranks, int *sizes, int ranks_size, int comm_size, int my_rank) 
{
    create_import_list(set, temp_list, h_list, total_size, ranks, sizes, ranks_size, comm_size, my_rank);
}

/*******************************************************************************
 * Main MPI halo creation routine
 *******************************************************************************/
void opp_halo_create() 
{
    opp_printf("opp_halo_create", OPP_my_rank, "start");

    // declare timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    double time;
    double max_time;
    op_timers(&cpu_t1, &wall_t1); // timer start for list create

    // create new communicator for OP mpi operation
    int my_rank, comm_size;
    // MPI_Comm_dup(OP_MPI_WORLD, &OP_MPI_WORLD);
    MPI_Comm_rank(OP_MPI_WORLD, &my_rank);
    MPI_Comm_size(OP_MPI_WORLD, &comm_size);

    /* Compute global partition range information for each set*/
    int **part_range = (int **)malloc(OP_set_index * sizeof(int *)); // part_range[set][ start and end ]
    get_part_range(part_range, my_rank, comm_size, OP_MPI_WORLD);

    // save this partition range information if it is not already saved during
    // a call to some partitioning routine
    if (orig_part_range == NULL) 
    {
        orig_part_range = (int **)malloc(OP_set_index * sizeof(int *));
        for (int s = 0; s < OP_set_index; s++) 
        {
            op_set set = OP_set_list[s];
            orig_part_range[set->index] = (int *)malloc(2 * comm_size * sizeof(int));
            for (int j = 0; j < comm_size; j++) 
            {
                orig_part_range[set->index][2 * j] = part_range[set->index][2 * j];
                orig_part_range[set->index][2 * j + 1] = part_range[set->index][2 * j + 1];
            }
        }
    }

    OP_export_exec_list = (halo_list *)malloc(OP_set_index * sizeof(halo_list));

    /*----- STEP 1 - Construct export lists for execute set elements and related mapping table entries -----*/

    // declare temporaty scratch variables to hold set export lists and mapping table export lists
    int s_i;
    int *set_list;

    int cap_s = 1000; // keep track of the temp array capacities

    for (int s = 0; s < OP_set_index; s++) // for each set
    { 
        op_set set = OP_set_list[s];

        // create a temporaty scratch space to hold export list for this set
        s_i = 0;
        cap_s = 1000;
        set_list = (int *)malloc(cap_s * sizeof(int));

        for (int e = 0; e < set->size; e++) // for each elment of this set
        {      
            for (int m = 0; m < OP_map_index; m++) // for each maping table
            { 
                op_map map = OP_map_list[m];

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
                            set_list = (int *)realloc(set_list, cap_s * sizeof(int));
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
        halo_list h_list = (halo_list)malloc(sizeof(halo_list_core));
        create_export_list(set, set_list, h_list, s_i, comm_size, my_rank);
        OP_export_exec_list[set->index] = h_list;
        free(set_list); // free temp list
    }

    /*---- STEP 2 - construct import lists for mappings and execute sets------*/

    OP_import_exec_list = (halo_list *)malloc(OP_set_index * sizeof(halo_list));

    int *neighbors, *sizes;
    int ranks_size;

    for (int s = 0; s < OP_set_index; s++) // for each set
    { 
        op_set set = OP_set_list[s];

        //-----Discover neighbors-----
        ranks_size = 0;
        neighbors = (int *)malloc(comm_size * sizeof(int));
        sizes = (int *)malloc(comm_size * sizeof(int));

        halo_list list = OP_export_exec_list[set->index];

        find_neighbors_set(list, neighbors, sizes, &ranks_size, my_rank, comm_size, OP_MPI_WORLD);
        MPI_Request request_send[list->ranks_size];

        int *rbuf, cap = 0, index = 0;

        for (int i = 0; i < list->ranks_size; i++) 
        {
            // printf("export from %d to %d set %10s, list of size %d \n", my_rank,list->ranks[i],set->name,list->sizes[i]);
            int *sbuf = &list->list[list->disps[i]];
            MPI_Isend(sbuf, list->sizes[i], MPI_INT, list->ranks[i], s, OP_MPI_WORLD, &request_send[i]);
        }

        for (int i = 0; i < ranks_size; i++)
            cap = cap + sizes[i];
        int *temp = (int *)malloc(cap * sizeof(int));

        // import this list from those neighbors
        for (int i = 0; i < ranks_size; i++) 
        {
            // printf("import from %d to %d set %10s, list of size %d\n", neighbors[i], my_rank, set->name, sizes[i]);
            rbuf = (int *)malloc(sizes[i] * sizeof(int));
            MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], s, OP_MPI_WORLD, MPI_STATUS_IGNORE);
            memcpy(&temp[index], (void *)&rbuf[0], sizes[i] * sizeof(int));
            index = index + sizes[i];
            free(rbuf);
        }

        MPI_Waitall(list->ranks_size, request_send, MPI_STATUSES_IGNORE);

        // create import lists
        halo_list h_list = (halo_list)malloc(sizeof(halo_list_core));
        create_import_list(set, temp, h_list, index, neighbors, sizes, ranks_size, comm_size, my_rank);
        OP_import_exec_list[set->index] = h_list;
    }

    /*--STEP 3 -Exchange mapping table entries using the import/export lists--*/

    for (int m = 0; m < OP_map_index; m++) // for each maping table
    { 
        op_map map = OP_map_list[m];
        halo_list i_list = OP_import_exec_list[map->from->index];
        halo_list e_list = OP_export_exec_list[map->from->index];

        MPI_Request request_send[e_list->ranks_size];

        // prepare bits of the mapping tables to be exported
        int **sbuf = (int **)malloc(e_list->ranks_size * sizeof(int *));

        for (int i = 0; i < e_list->ranks_size; i++) 
        {
            sbuf[i] = (int *)malloc((size_t)e_list->sizes[i] * map->dim * sizeof(int));
            for (int j = 0; j < e_list->sizes[i]; j++) 
            {
                for (int p = 0; p < map->dim; p++) 
                {
                    sbuf[i][j * map->dim + p] = map->map[map->dim * (e_list->list[e_list->disps[i] + j]) + p];
                }
            }

            MPI_Isend(sbuf[i], map->dim * e_list->sizes[i], MPI_INT, e_list->ranks[i],
                        m, OP_MPI_WORLD, &request_send[i]);
        }

        // prepare space for the incomming mapping tables - realloc each mapping tables in each mpi process
        OP_map_list[map->index]->map = (int *)realloc(OP_map_list[map->index]->map,
            (map->dim * (size_t)(map->from->size + i_list->size)) * sizeof(int));

        int init = map->dim * (map->from->size);
        for (int i = 0; i < i_list->ranks_size; i++) 
        {
            MPI_Recv( &(OP_map_list[map->index]->map[init + i_list->disps[i] * map->dim]),
                map->dim * i_list->sizes[i], MPI_INT, i_list->ranks[i], m, OP_MPI_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Waitall(e_list->ranks_size, request_send, MPI_STATUSES_IGNORE);
        
        for (int i = 0; i < e_list->ranks_size; i++)
            free(sbuf[i]);
        free(sbuf);
    }

    /*-- STEP 4 - Create import lists for non-execute set elements using mapping
        table entries including the additional mapping table entries --*/

    OP_import_nonexec_list = (halo_list *)malloc(OP_set_index * sizeof(halo_list));
    OP_export_nonexec_list = (halo_list *)malloc(OP_set_index * sizeof(halo_list));

    // declare temporaty scratch variables to hold non-exec set export lists
    s_i = 0;
    set_list = NULL;
    cap_s = 1000; // keep track of the temp array capacity

    for (int s = 0; s < OP_set_index; s++) // for each set
    { 
        op_set set = OP_set_list[s];
        halo_list exec_set_list = OP_import_exec_list[set->index];

        // create a temporaty scratch space to hold nonexec export list for this set
        s_i = 0;
        set_list = (int *)malloc(cap_s * sizeof(int));

        for (int m = 0; m < OP_map_index; m++) // for each maping table
        { 
            op_map map = OP_map_list[m];
            halo_list exec_map_list = OP_import_exec_list[map->from->index];

            if (compare_sets(map->to, set) == 1) // need to select  mappings TO this set
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
                            set_list = (int *)realloc(set_list, cap_s * sizeof(int));
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
        halo_list h_list = (halo_list)malloc(sizeof(halo_list_core));
        create_nonexec_import_list(set, set_list, h_list, s_i, comm_size, my_rank);
        free(set_list); // free temp list
        OP_import_nonexec_list[set->index] = h_list;
    }

    /*----------- STEP 5 - construct non-execute set export lists -------------*/

    for (int s = 0; s < OP_set_index; s++) // for each set
    { 
        op_set set = OP_set_list[s];

        //-----Discover neighbors-----
        ranks_size = 0;
        neighbors = (int *)malloc(comm_size * sizeof(int));
        sizes = (int *)malloc(comm_size * sizeof(int));

        halo_list list = OP_import_nonexec_list[set->index];
        find_neighbors_set(list, neighbors, sizes, &ranks_size, my_rank, comm_size,
                        OP_MPI_WORLD);

        MPI_Request request_send[list->ranks_size];
        int *rbuf, cap = 0, index = 0;

        for (int i = 0; i < list->ranks_size; i++) 
        {
            int *sbuf = &list->list[list->disps[i]];
            MPI_Isend(sbuf, list->sizes[i], MPI_INT, list->ranks[i], s, OP_MPI_WORLD, &request_send[i]);
        }

        for (int i = 0; i < ranks_size; i++)
            cap = cap + sizes[i];
        int *temp = (int *)malloc(cap * sizeof(int));

        // export this list to those neighbors
        for (int i = 0; i < ranks_size; i++) 
        {
            rbuf = (int *)malloc(sizes[i] * sizeof(int));
            MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], s, OP_MPI_WORLD, MPI_STATUS_IGNORE);
            memcpy(&temp[index], (void *)&rbuf[0], sizes[i] * sizeof(int));
            index = index + sizes[i];
            free(rbuf);
        }

        MPI_Waitall(list->ranks_size, request_send, MPI_STATUSES_IGNORE);

        halo_list h_list = (halo_list)malloc(sizeof(halo_list_core));
        create_nonexec_export_list(set, temp, h_list, index, neighbors, sizes, ranks_size, comm_size, my_rank);
        OP_export_nonexec_list[set->index] = h_list;
    }

    /*-STEP 6 - Exchange execute set elements/data using the import/export lists--*/

    for (int s = 0; s < OP_set_index; s++) // for each set
    { 
        op_set set = OP_set_list[s];
        halo_list i_list = OP_import_exec_list[set->index];
        halo_list e_list = OP_export_exec_list[set->index];

        int d = -1; // d is just simply the tag for mpi comms
        for (int k = 0; k < OP_dat_index; k++) // for each dat
        {
            d++; // increase tag to do mpi comm for the next op_dat
            op_dat dat = OP_dat_list[k];

            if (compare_sets(set, dat->set) == 1) // if this data array is defined on this set
            { 
                MPI_Request request_send[e_list->ranks_size];

                // prepare execute set element data to be exported
                char **sbuf = (char **)malloc(e_list->ranks_size * sizeof(char *));

                for (int i = 0; i < e_list->ranks_size; i++) 
                {
                    sbuf[i] = (char *)malloc((size_t)e_list->sizes[i] * (size_t)dat->size);
                    for (int j = 0; j < e_list->sizes[i]; j++) 
                    {
                        int set_elem_index = e_list->list[e_list->disps[i] + j];
                        memcpy(&sbuf[i][j * (size_t)dat->size],
                            (void *)&dat->data[(size_t)dat->size * (set_elem_index)], dat->size);
                    }

                    MPI_Isend(sbuf[i], (size_t)dat->size * e_list->sizes[i], MPI_CHAR,
                                e_list->ranks[i], d, OP_MPI_WORLD, &request_send[i]);
                }

                // prepare space for the incomming data - realloc each
                // data array in each mpi process
                dat->data =
                    (char *)realloc(dat->data, (size_t)(set->size + i_list->size) * (size_t)dat->size);

                int init = set->size * (size_t)dat->size;
                for (int i = 0; i < i_list->ranks_size; i++) 
                {
                    MPI_Recv(&(dat->data[init + i_list->disps[i] * (size_t)dat->size]),
                            (size_t)dat->size * i_list->sizes[i], MPI_CHAR, i_list->ranks[i], d,
                            OP_MPI_WORLD, MPI_STATUS_IGNORE);
                }

                MPI_Waitall(e_list->ranks_size, request_send, MPI_STATUSES_IGNORE);
                
                for (int i = 0; i < e_list->ranks_size; i++)
                    free(sbuf[i]);
                free(sbuf);
            }
        }
    }

    /*-STEP 7 - Exchange non-execute set elements/data using the import/export lists--*/

    for (int s = 0; s < OP_set_index; s++) // for each set
    { 
        op_set set = OP_set_list[s];
        halo_list i_list = OP_import_nonexec_list[set->index];
        halo_list e_list = OP_export_nonexec_list[set->index];

        int d = -1; // d is just simply the tag for mpi comms

        for (int k = 0; k < OP_dat_index; k++) // for each dat
        {
            d++; // increase tag to do mpi comm for the next op_dat
            op_dat dat = OP_dat_list[k];

            if (compare_sets(set, dat->set) == 1)  // if this data array is defined on this set
            {
                // printf("on rank %d, The data array is %10s\n",my_rank,dat->name);
                MPI_Request request_send[e_list->ranks_size];

                // prepare non-execute set element data to be exported
                char **sbuf = (char **)malloc(e_list->ranks_size * sizeof(char *));

                for (int i = 0; i < e_list->ranks_size; i++) 
                {
                    sbuf[i] = (char *)malloc(e_list->sizes[i] * (size_t)dat->size);
                    for (int j = 0; j < e_list->sizes[i]; j++) 
                    {
                        int set_elem_index = e_list->list[e_list->disps[i] + j];
                        memcpy(&sbuf[i][j * (size_t)dat->size], (void *)&dat->data[(size_t)dat->size * (set_elem_index)], dat->size);
                    }
                    MPI_Isend(sbuf[i], (size_t)dat->size * e_list->sizes[i], MPI_CHAR,
                                e_list->ranks[i], d, OP_MPI_WORLD, &request_send[i]);
                }

                // prepare space for the incomming nonexec-data - realloc each data array in each mpi process
                halo_list exec_i_list = OP_import_exec_list[set->index];

                dat->data = (char *)realloc(dat->data,
                    (size_t)(set->size + exec_i_list->size + i_list->size) * (size_t)dat->size);

                size_t init = (size_t)(set->size + exec_i_list->size) * (size_t)dat->size;
                for (int i = 0; i < i_list->ranks_size; i++) {
                MPI_Recv(&(dat->data[init + i_list->disps[i] * (size_t)dat->size]),
                        (size_t)dat->size * i_list->sizes[i], MPI_CHAR, i_list->ranks[i], d,
                        OP_MPI_WORLD, MPI_STATUS_IGNORE);
                }

                MPI_Waitall(e_list->ranks_size, request_send, MPI_STATUSES_IGNORE);
                for (int i = 0; i < e_list->ranks_size; i++)
                    free(sbuf[i]);
                free(sbuf);
            }
        }
    }

    /*-STEP 8 ----------------- Renumber Mapping tables-----------------------*/

    for (int s = 0; s < OP_set_index; s++) // for each set
    { 
        op_set set = OP_set_list[s];

        for (int m = 0; m < OP_map_index; m++) // for each maping table
        { 
            op_map map = OP_map_list[m];

            if (compare_sets(map->to, set) == 1) { // need to select
                                                    // mappings TO this set

                halo_list exec_set_list = OP_import_exec_list[set->index];
                halo_list nonexec_set_list = OP_import_nonexec_list[set->index];

                halo_list exec_map_list = OP_import_exec_list[map->from->index];

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
                            OP_map_list[map->index]->map[e * map->dim + j] = local_index;
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
                                    OP_map_list[map->index]->map[e * map->dim + j] = found + map->to->size;
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
                                    OP_map_list[map->index]->map[e * map->dim + j] = found + set->size + exec_set_list->size;
                                }
                            }

                            if (found < 0)
                                printf("ERROR: Set %10s Element %d needed on rank %d from partition %d\n",
                                    set->name, local_index, my_rank, part);
                        }
                    }
                }
            }
        }
    }

    /*-STEP 9 ---------------- Create MPI send Buffers-----------------------*/

    for (int k = 0; k < OP_dat_index; k++) // for each dat
    {
        op_dat dat = OP_dat_list[k];

        op_mpi_buffer mpi_buf = (op_mpi_buffer)malloc(sizeof(op_mpi_buffer_core));

        halo_list exec_e_list = OP_export_exec_list[dat->set->index];
        halo_list nonexec_e_list = OP_export_nonexec_list[dat->set->index];

        mpi_buf->buf_exec = (char *)malloc((size_t)(exec_e_list->size) * (size_t)dat->size);
        mpi_buf->buf_nonexec = (char *)malloc((size_t)(nonexec_e_list->size) * (size_t)dat->size);

        halo_list exec_i_list = OP_import_exec_list[dat->set->index];
        halo_list nonexec_i_list = OP_import_nonexec_list[dat->set->index];

        mpi_buf->s_req = (MPI_Request *)malloc(
            sizeof(MPI_Request) * (exec_e_list->ranks_size + nonexec_e_list->ranks_size));
        mpi_buf->r_req = (MPI_Request *)malloc(
            sizeof(MPI_Request) * (exec_i_list->ranks_size + nonexec_i_list->ranks_size));

        mpi_buf->s_num_req = 0;
        mpi_buf->r_num_req = 0;
        dat->mpi_buffer = mpi_buf;
    }

    // set dirty bits of all data arrays to 0 for each data array
    for (int k = 0; k < OP_dat_index; k++) // for each dat
    {
        op_dat dat = OP_dat_list[k];
        dat->dirtybit = 0;
    }

    /*-STEP 10 -------------------- Separate core elements------------------------*/

    int **core_elems = (int **)malloc(OP_set_index * sizeof(int *));
    int **exp_elems = (int **)malloc(OP_set_index * sizeof(int *));

    for (int s = 0; s < OP_set_index; s++) // for each set
    { 
        op_set set = OP_set_list[s];

        halo_list exec = OP_export_exec_list[set->index];
        halo_list nonexec = OP_export_nonexec_list[set->index];

        if (exec->size > 0) 
        {
            exp_elems[set->index] = (int *)malloc(exec->size * sizeof(int));
            memcpy(exp_elems[set->index], exec->list, exec->size * sizeof(int));
            quickSort(exp_elems[set->index], 0, exec->size - 1);

            int num_exp = removeDups(exp_elems[set->index], exec->size);
            core_elems[set->index] = (int *)malloc(set->size * sizeof(int));
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
            for (int k = 0; k < OP_dat_index; k++) // for each dat
            {
                op_dat dat = OP_dat_list[k];

                if (compare_sets(set, dat->set) == 1) // if this data array is defined on this set
                {
                    char *new_dat = (char *)malloc((size_t)set->size * (size_t)dat->size);
                    for (int i = 0; i < count; i++) {
                        memcpy(&new_dat[i * (size_t)dat->size], &dat->data[core_elems[set->index][i] * (size_t)dat->size], dat->size);
                    }
                    for (int i = 0; i < num_exp; i++) {
                        memcpy(&new_dat[(count + i) * (size_t)dat->size], &dat->data[exp_elems[set->index][i] * (size_t)dat->size], dat->size);
                    }
                    memcpy(&dat->data[0], &new_dat[0], set->size * (size_t)dat->size);
                    free(new_dat);
                }
            }

            // for each mapping defined from this set seperate its elements
            for (int m = 0; m < OP_map_index; m++) // for each set
            { 
                op_map map = OP_map_list[m];

                if (compare_sets(map->from, set) == 1)  // if this mapping is defined from this set
                { 
                    int *new_map = (int *)malloc((size_t)set->size * map->dim * sizeof(int));
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
                    free(new_map);
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
            core_elems[set->index] = (int *)malloc(set->size * sizeof(int));
            exp_elems[set->index] = (int *)malloc(0 * sizeof(int));
            for (int e = 0; e < set->size; e++)  // for each elment of this set
            {
                core_elems[set->index][e] = e;
            }
            set->core_size = set->size;
        }
    }

    // now need to renumber mapping tables as the elements are seperated
    for (int m = 0; m < OP_map_index; m++) { // for each set
        op_map map = OP_map_list[m];

        halo_list exec_map_list = OP_import_exec_list[map->from->index];
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
                            printf("Problem in seperating core elements - renumbering map\n");
                        else
                            OP_map_list[map->index]->map[e * (size_t)map->dim + j] = map->to->core_size + index;
                    } 
                    else
                        OP_map_list[map->index]->map[e * (size_t)map->dim + j] = index;
                }
            }
        }
    }

    /*-STEP 11 ----------- Save the original set element indexes------------------*/

    // if OP_part_list is empty, (i.e. no previous partitioning done) then create it and store the seperation of 
    // elements using core_elems and exp_elems
    if (OP_part_index != OP_set_index) 
    {
        // allocate memory for list
        OP_part_list = (part *)malloc(OP_set_index * sizeof(part));

        for (int s = 0; s < OP_set_index; s++)  // for each set
        {
            op_set set = OP_set_list[s];

            int *g_index = (int *)malloc(sizeof(int) * set->size);
            int *partition = (int *)malloc(sizeof(int) * set->size);
            for (int i = 0; i < set->size; i++) 
            {
                g_index[i] = get_global_index(i, my_rank, part_range[set->index], comm_size);
                partition[i] = my_rank;
            }
            decl_partition(set, g_index, partition);

            // combine core_elems and exp_elems to one memory block
            int *temp = (int *)malloc(sizeof(int) * set->size);
            memcpy(&temp[0], core_elems[set->index], set->core_size * sizeof(int));
            memcpy(&temp[set->core_size], exp_elems[set->index],
                    (set->size - set->core_size) * sizeof(int));

            // update OP_part_list[set->index]->g_index
            for (int i = 0; i < set->size; i++) 
            {
                temp[i] = OP_part_list[set->index]->g_index[temp[i]];
            }
            free(OP_part_list[set->index]->g_index);
            OP_part_list[set->index]->g_index = temp;
        }
    } 
    else  // OP_part_list exists (i.e. a partitioning has been done) update the seperation of elements
    {        

        for (int s = 0; s < OP_set_index; s++)  // for each set
        { 
            op_set set = OP_set_list[s];

            // combine core_elems and exp_elems to one memory block
            int *temp = (int *)malloc(sizeof(int) * set->size);
            memcpy(&temp[0], core_elems[set->index], set->core_size * sizeof(int));
            memcpy(&temp[set->core_size], exp_elems[set->index], (set->size - set->core_size) * sizeof(int));

            // update OP_part_list[set->index]->g_index
            for (int i = 0; i < set->size; i++) 
            {
                temp[i] = OP_part_list[set->index]->g_index[temp[i]];
            }
            free(OP_part_list[set->index]->g_index);
            OP_part_list[set->index]->g_index = temp;
        }
    }

    // set up exec and nonexec sizes
    for (int s = 0; s < OP_set_index; s++) // for each set
    { 
        op_set set = OP_set_list[s];
        set->exec_size = OP_import_exec_list[set->index]->size;
        set->nonexec_size = OP_import_nonexec_list[set->index]->size;
    }

    /*-STEP 12 ---------- Clean up and Compute rough halo size numbers------------*/

    for (int i = 0; i < OP_set_index; i++) 
    {
        free(part_range[i]);
        free(core_elems[i]);
        free(exp_elems[i]);
    }
    free(part_range);
    free(exp_elems);
    free(core_elems);

    op_timers(&cpu_t2, &wall_t2); // timer stop for list create
    // compute import/export lists creation time
    time = wall_t2 - wall_t1;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, OPP_MPI_ROOT, OP_MPI_WORLD);

    // compute avg/min/max set sizes and exec sizes accross the MPI universe
    int avg_size = 0, min_size = 0, max_size = 0;
    for (int s = 0; s < OP_set_index; s++) 
    {
        op_set set = OP_set_list[s];

        // number of set elements first
        MPI_Reduce(&set->size, &avg_size, 1, MPI_INT, MPI_SUM, OPP_MPI_ROOT,
                OP_MPI_WORLD);
        MPI_Reduce(&set->size, &min_size, 1, MPI_INT, MPI_MIN, OPP_MPI_ROOT,
                OP_MPI_WORLD);
        MPI_Reduce(&set->size, &max_size, 1, MPI_INT, MPI_MAX, OPP_MPI_ROOT,
                OP_MPI_WORLD);

        if (my_rank == OPP_MPI_ROOT) 
        {
            printf("Num of %8s (avg | min | max)\n", set->name);
            printf("total elems         %10d %10d %10d\n", avg_size / comm_size, min_size, max_size);
        }

        avg_size = 0;
        min_size = 0;
        max_size = 0;

        // number of OWNED elements second
        MPI_Reduce(&set->core_size, &avg_size, 1, MPI_INT, MPI_SUM, OPP_MPI_ROOT,
                OP_MPI_WORLD);
        MPI_Reduce(&set->core_size, &min_size, 1, MPI_INT, MPI_MIN, OPP_MPI_ROOT,
                OP_MPI_WORLD);
        MPI_Reduce(&set->core_size, &max_size, 1, MPI_INT, MPI_MAX, OPP_MPI_ROOT,
                OP_MPI_WORLD);

        if (my_rank == OPP_MPI_ROOT) 
        {
            printf("core elems         %10d %10d %10d \n", avg_size / comm_size, min_size, max_size);
        }
        avg_size = 0;
        min_size = 0;
        max_size = 0;

        // number of exec halo elements third
        MPI_Reduce(&OP_import_exec_list[set->index]->size, &avg_size, 1, MPI_INT,
                MPI_SUM, OPP_MPI_ROOT, OP_MPI_WORLD);
        MPI_Reduce(&OP_import_exec_list[set->index]->size, &min_size, 1, MPI_INT,
                MPI_MIN, OPP_MPI_ROOT, OP_MPI_WORLD);
        MPI_Reduce(&OP_import_exec_list[set->index]->size, &max_size, 1, MPI_INT,
                MPI_MAX, OPP_MPI_ROOT, OP_MPI_WORLD);

        if (my_rank == OPP_MPI_ROOT) 
        {
            printf("exec halo elems     %10d %10d %10d \n", avg_size / comm_size, min_size, max_size);
        }
        avg_size = 0;
        min_size = 0;
        max_size = 0;

        // number of non-exec halo elements fourth
        MPI_Reduce(&OP_import_nonexec_list[set->index]->size, &avg_size, 1, MPI_INT,
                MPI_SUM, OPP_MPI_ROOT, OP_MPI_WORLD);
        MPI_Reduce(&OP_import_nonexec_list[set->index]->size, &min_size, 1, MPI_INT,
                MPI_MIN, OPP_MPI_ROOT, OP_MPI_WORLD);
        MPI_Reduce(&OP_import_nonexec_list[set->index]->size, &max_size, 1, MPI_INT,
                MPI_MAX, OPP_MPI_ROOT, OP_MPI_WORLD);

        if (my_rank == OPP_MPI_ROOT) 
        {
            printf("non-exec halo elems %10d %10d %10d \n", avg_size / comm_size, min_size, max_size);
        }
        avg_size = 0;
        min_size = 0;
        max_size = 0;
        if (my_rank == OPP_MPI_ROOT) 
        {
            printf("-----------------------------------------------------\n");
        }
    }

    if (my_rank == OPP_MPI_ROOT) 
    {
        printf("\n\n");
    }

    // compute avg/min/max number of MPI neighbors per process accross the MPI
    // universe
    avg_size = 0, min_size = 0, max_size = 0;
    for (int s = 0; s < OP_set_index; s++) 
    {
        op_set set = OP_set_list[s];

        // number of exec halo neighbors first
        MPI_Reduce(&OP_import_exec_list[set->index]->ranks_size, &avg_size, 1,
                MPI_INT, MPI_SUM, OPP_MPI_ROOT, OP_MPI_WORLD);
        MPI_Reduce(&OP_import_exec_list[set->index]->ranks_size, &min_size, 1,
                MPI_INT, MPI_MIN, OPP_MPI_ROOT, OP_MPI_WORLD);
        MPI_Reduce(&OP_import_exec_list[set->index]->ranks_size, &max_size, 1,
                MPI_INT, MPI_MAX, OPP_MPI_ROOT, OP_MPI_WORLD);

        if (my_rank == OPP_MPI_ROOT) 
        {
            printf("MPI neighbors for exchanging %8s (avg | min | max)\n", set->name);
            printf("exec halo elems     %4d %4d %4d\n", avg_size / comm_size, min_size, max_size);
        }
        avg_size = 0;
        min_size = 0;
        max_size = 0;

        // number of non-exec halo neighbors second
        MPI_Reduce(&OP_import_nonexec_list[set->index]->ranks_size, &avg_size, 1,
                MPI_INT, MPI_SUM, OPP_MPI_ROOT, OP_MPI_WORLD);
        MPI_Reduce(&OP_import_nonexec_list[set->index]->ranks_size, &min_size, 1,
                MPI_INT, MPI_MIN, OPP_MPI_ROOT, OP_MPI_WORLD);
        MPI_Reduce(&OP_import_nonexec_list[set->index]->ranks_size, &max_size, 1,
                MPI_INT, MPI_MAX, OPP_MPI_ROOT, OP_MPI_WORLD);

        if (my_rank == OPP_MPI_ROOT) 
        {
            printf("non-exec halo elems %4d %4d %4d\n", avg_size / comm_size, min_size, max_size);
        }
        avg_size = 0;
        min_size = 0;
        max_size = 0;
        if (my_rank == OPP_MPI_ROOT) {
            printf("-----------------------------------------------------\n");
        }
    }

    // compute average worst case halo size in Bytes
    int tot_halo_size = 0;
    for (int s = 0; s < OP_set_index; s++) 
    {
        op_set set = OP_set_list[s];

        for (int k = 0; k < OP_dat_index; k++) // for each dat
        {
            op_dat dat = OP_dat_list[k];

            if (compare_sets(dat->set, set) == 1) 
            {
                halo_list exec_imp = OP_import_exec_list[set->index];
                halo_list nonexec_imp = OP_import_nonexec_list[set->index];
                tot_halo_size = tot_halo_size + exec_imp->size * (size_t)dat->size +
                                nonexec_imp->size * (size_t)dat->size;
            }
        }
    }
    int avg_halo_size;
    MPI_Reduce(&tot_halo_size, &avg_halo_size, 1, MPI_INT, MPI_SUM, OPP_MPI_ROOT, OP_MPI_WORLD);

    // print performance results
    if (my_rank == OPP_MPI_ROOT) 
    {
        printf("Max total halo creation time = %lf\n", max_time);
        printf("Average (worst case) Halo size = %d Bytes\n", avg_halo_size / comm_size);
    }

    opp_printf("opp_halo_create", OPP_my_rank, "end");
}

/*******************************************************************************
 * Routine to Clean-up all MPI halos(called at the end of an OP2 MPI application)
*******************************************************************************/
void opp_halo_destroy() 
{
    // remove halos from op_dats
    for (int k = 0; k < OP_dat_index; k++) // for each dat
    {
        op_dat dat = OP_dat_list[k];
        dat->data = (char *)realloc(dat->data, (size_t)dat->set->size * dat->size);
    }

    //free lists
    // for (int s = 0; s < OP_set_index; s++) 
    // {
    //     op_set set = OP_set_list[s];

    //     if (set->is_particle) continue;

    //     free(OP_import_exec_list[set->index]->ranks);
    //     free(OP_import_exec_list[set->index]->disps);
    //     free(OP_import_exec_list[set->index]->sizes);
    //     free(OP_import_exec_list[set->index]->list);
    //     free(OP_import_exec_list[set->index]);

    //     free(OP_import_nonexec_list[set->index]->ranks);
    //     free(OP_import_nonexec_list[set->index]->disps);
    //     free(OP_import_nonexec_list[set->index]->sizes);
    //     free(OP_import_nonexec_list[set->index]->list);
    //     free(OP_import_nonexec_list[set->index]);

    //     free(OP_export_exec_list[set->index]->ranks);
    //     free(OP_export_exec_list[set->index]->disps);
    //     free(OP_export_exec_list[set->index]->sizes);
    //     free(OP_export_exec_list[set->index]->list);
    //     free(OP_export_exec_list[set->index]);

    //     free(OP_export_nonexec_list[set->index]->ranks);
    //     free(OP_export_nonexec_list[set->index]->disps);
    //     free(OP_export_nonexec_list[set->index]->sizes);
    //     free(OP_export_nonexec_list[set->index]->list);
    //     free(OP_export_nonexec_list[set->index]);
    // }
    // free(OP_import_exec_list);
    // free(OP_import_nonexec_list);
    // free(OP_export_exec_list);
    // free(OP_export_nonexec_list);

    // for (int k = 0; k < OP_dat_index; k++) // for each dat
    // {
    //     op_dat dat = OP_dat_list[k];

    //     if (dat->set->is_particle) continue;

    //     free(((op_mpi_buffer)(dat->mpi_buffer))->buf_exec);
    //     free(((op_mpi_buffer)(dat->mpi_buffer))->buf_nonexec);
    //     free(((op_mpi_buffer)(dat->mpi_buffer))->s_req);
    //     free(((op_mpi_buffer)(dat->mpi_buffer))->r_req);
    // }
}
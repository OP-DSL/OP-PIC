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
 * op_mpi_part_core.c
 *
 * Implements the OP2 Distributed memory (MPI) Partitioning wrapper routines,
 * data migration and support utility functions
 *
 * written by: Gihan R. Mudalige, (Started 07-04-2011)
 */

#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <float.h> // double min/max

#include <opp_mpi.h>

#ifdef HAVE_PARMETIS
#include <parmetis.h>

// #ifdef PARMETIS_VER_4
    typedef idx_t idxtype;
// #else
//     typedef int idx_t;
// #endif
#endif

//**********************************
// TODO : Remove these
double OPP_hybrid_balance = 1.0;
#define real_t float
//**********************************


#ifdef __cplusplus
extern "C" {
#endif

MPI_Comm OPP_PART_WORLD;

/*******************************************************************************
 * Utility function to find the number of times a value appears in an array
 *******************************************************************************/
/* static */ int frequencyof(int value, int *array, int size) 
{
    int frequency = 0;
    for (int i = 0; i < size; i++) 
    {
        if (array[i] == value)
            frequency++;
    }
    return frequency;
}

/*******************************************************************************
 * Utility function to find the mode of a set of numbers in an array
 *******************************************************************************/
/* static */ int find_mode(int *array, int size) 
{
    int count = 0, mode = array[0], current;
    for (int i = 0; i < size; i++) 
    {
        if (i > 0 && array[i] == array[i - 1])
            continue;
        current = frequencyof(array[i], array, size);
        if (count < current) 
        {
            count = current;
            mode = array[i];
        }
    }
    return mode;
}

/*******************************************************************************
 * Utility function to see if a target opp_set is held in an opp_set array
 *******************************************************************************/
/* static */ int compare_all_sets(opp_set target_set, opp_set other_sets[], int size) 
{
    for (int i = 0; i < size; i++) 
    {
        if (compare_sets(target_set, other_sets[i]) == 1)
            return i;
    }
    return -1;
}

/*******************************************************************************
 * Special routine to create export list during partitioning map->to set
 * from map_>from set in partition_to_set()
 *******************************************************************************/
/* static */ int *create_exp_list_2(opp_set set, int *temp_list, halo_list h_list, int *part_list, int size, 
                            int comm_size, int my_rank) 
{
    (void)my_rank;
    int *ranks = (int *)opp_host_malloc(comm_size * sizeof(int));
    int *to_list = (int *)opp_host_malloc((size / 3) * sizeof(int));
    part_list = (int *)opp_host_malloc((size / 3) * sizeof(int));
    int *disps = (int *)opp_host_malloc(comm_size * sizeof(int));
    int *sizes = (int *)opp_host_malloc(comm_size * sizeof(int));

    int index = 0;
    size_t total_size = 0;

    // negative values set as an initialisation
    for (int r = 0; r < comm_size; r++) 
    {
        disps[r] = ranks[r] = -99;
        sizes[r] = 0;
    }

    for (int r = 0; r < comm_size; r++) 
    {
        sizes[index] = 0;
        disps[index] = 0;
        int *temp_to = (int *)opp_host_malloc((size / 3) * sizeof(int));
        int *temp_part = (int *)opp_host_malloc((size / 3) * sizeof(int));

        for (int i = 0; i < size; i = i + 3) 
        {
            if (temp_list[i] == r) 
            {
                temp_to[sizes[index]] = temp_list[i + 1];
                temp_part[sizes[index]] = temp_list[i + 2];
                sizes[index]++;
            }
        }

        if (sizes[index] > 0) 
        {
            ranks[index] = r;
            // no sorting
            total_size = total_size + sizes[index];
            // no eleminating duplicates
            if (index > 0)
                disps[index] = disps[index - 1] + sizes[index - 1];

            // add to end of t_list and p_list
            for (int e = 0; e < sizes[index]; e++) 
            {
                to_list[disps[index] + e] = temp_to[e];
                part_list[disps[index] + e] = temp_part[e];
            }
            index++;
        }
        opp_host_free(temp_to);
        opp_host_free(temp_part);
    }

    h_list->set = set;
    h_list->size = total_size;
    h_list->ranks = ranks;
    h_list->ranks_size = index;
    h_list->disps = disps;
    h_list->sizes = sizes;
    h_list->list = to_list;

    return part_list;
}

/*******************************************************************************
 * Special routine to create import list during partitioning map->to set
 * from map_>from set in partition_to_set()
 *******************************************************************************/
/* static */ void create_imp_list_2(opp_set set, int *temp_list, halo_list h_list, int total_size, int *ranks, int *sizes,
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
 * Routine to use a partitioned map->to set to partition the map->from set
 *******************************************************************************/
/* static */ int partition_from_set(opp_map map, int my_rank, int comm_size, int **part_range) 
{
    (void)my_rank;
    part p_set = OPP_part_list[map->to->index];

    size_t cap = 100;
    size_t count = 0;
    int *temp_list = (int *)opp_host_malloc(cap * sizeof(int));

    halo_list pi_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));

    // go through the map and build an import list of the non-local "to" elements
    for (int i = 0; i < map->from->size; i++) 
    {
        int part, local_index;
        for (int j = 0; j < map->dim; j++) 
        {
            part = get_partition(map->map[i * map->dim + j], part_range[map->to->index], &local_index, comm_size);
            if (count >= cap) 
            {
                cap = cap * 2;
                temp_list = (int *)opp_host_realloc(temp_list, cap * sizeof(int));
            }

            if (part != my_rank) 
            {
                temp_list[count++] = part;
                temp_list[count++] = local_index;
            }
        }
    }
    create_export_list(map->to, temp_list, pi_list, count, comm_size, my_rank);
    opp_host_free(temp_list);

    // now, discover neighbors and create export list of "to" elements
    int ranks_size = 0;

    int *neighbors = (int *)opp_host_malloc(comm_size * sizeof(int));
    int *sizes = (int *)opp_host_malloc(comm_size * sizeof(int));

    halo_list pe_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));
    find_neighbors_set(pi_list, neighbors, sizes, &ranks_size, my_rank, comm_size, OPP_PART_WORLD);

    //  MPI_Request request_send[pi_list->ranks_size];
    std::vector<MPI_Request> request_send(pi_list->ranks_size);

    int *rbuf;
    cap = 0;
    count = 0;

    for (int i = 0; i < pi_list->ranks_size; i++) 
    {
        int *sbuf = &pi_list->list[pi_list->disps[i]];
        MPI_Isend(sbuf, pi_list->sizes[i], MPI_INT, pi_list->ranks[i], 1, OPP_PART_WORLD, &request_send[i]);
    }

    for (int i = 0; i < ranks_size; i++)
        cap = cap + sizes[i];
    temp_list = (int *)opp_host_malloc(cap * sizeof(int));

    for (int i = 0; i < ranks_size; i++) 
    {
        rbuf = (int *)opp_host_malloc(sizes[i] * sizeof(int));
        MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], 1, OPP_PART_WORLD, MPI_STATUS_IGNORE);
        memcpy(&temp_list[count], (void *)&rbuf[0], sizes[i] * sizeof(int));
        count = count + sizes[i];
        opp_host_free(rbuf);
    }
    MPI_Waitall(pi_list->ranks_size, request_send.data(), MPI_STATUSES_IGNORE);
    create_import_list(map->to, temp_list, pe_list, count, neighbors, sizes, ranks_size, comm_size, my_rank);

    // use the import and export lists to exchange partition information of this "to" set
    std::vector<MPI_Request> request_send_p(pe_list->ranks_size);

    // first - prepare partition information of the "to" set element to be exported
    int **sbuf = (int **)opp_host_malloc(pe_list->ranks_size * sizeof(int *));
    for (int i = 0; i < pe_list->ranks_size; i++) 
    {
        sbuf[i] = (int *)opp_host_malloc(pe_list->sizes[i] * sizeof(int));
        for (int j = 0; j < pe_list->sizes[i]; j++) 
        {
            int elem = pe_list->list[pe_list->disps[i] + j];
            sbuf[i][j] = p_set->elem_part[elem];
        }
        MPI_Isend(sbuf[i], pe_list->sizes[i], MPI_INT, pe_list->ranks[i], 2, OPP_PART_WORLD, &request_send_p[i]);
    }

    // second - prepare space for the incomming partition information of the "to" set
    int *imp_part = (int *)opp_host_malloc(sizeof(int) * pi_list->size);

    // third - receive
    for (int i = 0; i < pi_list->ranks_size; i++) 
    {
        MPI_Recv(&imp_part[pi_list->disps[i]], pi_list->sizes[i], MPI_INT, pi_list->ranks[i], 2, OPP_PART_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Waitall(pe_list->ranks_size, request_send_p.data(), MPI_STATUSES_IGNORE);
    for (int i = 0; i < pe_list->ranks_size; i++)
        opp_host_free(sbuf[i]);
    opp_host_free(sbuf);

    // allocate memory to hold the partition details for the set thats going to be partitioned
    int *partition = (int *)opp_host_malloc(sizeof(int) * map->from->size);

    // go through the mapping table and the imported partition information and partition the "from" set
    for (int i = 0; i < map->from->size; i++) 
    {
        int part, local_index;
        int *found_parts = (int*)opp_host_malloc(sizeof(int)*map->dim);
        for (int j = 0; j < map->dim; j++) 
        {
            part = get_partition(map->map[i * map->dim + j], part_range[map->to->index], &local_index, comm_size);

            if (part == my_rank)
                found_parts[j] = p_set->elem_part[local_index];
            else // get partition information from imported data
            {
                int r = binary_search(pi_list->ranks, part, 0, pi_list->ranks_size - 1);
                if (r >= 0) 
                {
                    int elem = binary_search(&pi_list->list[pi_list->disps[r]], local_index, 0, pi_list->sizes[r] - 1);
                    if (elem >= 0)
                        found_parts[j] = imp_part[pi_list->disps[r] + elem];
                    else 
                    {
                        printf("Element %d not found in partition import list\n", local_index);
                        opp_abort("Element not found in partition import list");
                    }
                } 
                else 
                {
                    printf("Rank %d not found in partition import list\n", part);
                    opp_abort("Rank not found in partition import list");
                }
            }
        }
        partition[i] = find_mode(found_parts, map->dim);
        opp_host_free(found_parts);
    }

    OPP_part_list[map->from->index]->elem_part = partition;
    OPP_part_list[map->from->index]->is_partitioned = 1;

    // cleanup
    opp_host_free(imp_part);
    opp_host_free(pi_list->list);
    opp_host_free(pi_list->ranks);
    opp_host_free(pi_list->sizes);
    opp_host_free(pi_list->disps);
    opp_host_free(pi_list);
    opp_host_free(pe_list->list);
    opp_host_free(pe_list->ranks);
    opp_host_free(pe_list->sizes);
    opp_host_free(pe_list->disps);
    opp_host_free(pe_list);

    return 1;
}

/*******************************************************************************
 * Routine to use the partitioned map->from set to partition the map->to set
 *******************************************************************************/
/* static */ int partition_to_set(opp_map map, int my_rank, int comm_size, int **part_range) 
{
    part p_set = OPP_part_list[map->from->index];

    int cap = 300;
    int count = 0;
    int *temp_list = (int *)opp_host_malloc(cap * sizeof(int));

    halo_list pe_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));
    int *part_list_e = NULL; // corresponding "to" element's partition infomation
    // exported to an mpi rank

    // go through the map and if any element pointed to by a mapping table entry
    //(i.e. a "from" set element) is in a foreign partition, add the partition
    // of the from element to be exported to that mpi foreign process
    // also collect information about the local "to" elements
    for (int i = 0; i < map->from->size; i++) 
    {
        int part;
        int local_index;

        for (int j = 0; j < map->dim; j++) 
        {
            part = get_partition(map->map[i * map->dim + j], part_range[map->to->index], &local_index, comm_size);

            if (part != my_rank) 
            {
                if (count >= cap) 
                {
                    cap = cap * 3;
                    temp_list = (int *)opp_host_realloc(temp_list, cap * sizeof(int));
                }

                temp_list[count++] = part; // curent partition (i.e. mpi rank)
                temp_list[count++] = local_index; // map->map[i*map->dim+j];//global index
                temp_list[count++] = p_set->elem_part[i]; // new partition
            }
        }
    }

    part_list_e = create_exp_list_2(map->to, temp_list, pe_list, part_list_e, count, comm_size, my_rank);
    opp_host_free(temp_list);

    int ranks_size = 0;
    int *neighbors = (int *)opp_host_malloc(comm_size * sizeof(int));
    int *sizes = (int *)opp_host_malloc(comm_size * sizeof(int));

    // to_part_list tpi_list;
    halo_list pi_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));
    int *part_list_i = NULL; // corresponding "to" element's partition infomation imported from an mpi rank

    find_neighbors_set(pe_list, neighbors, sizes, &ranks_size, my_rank, comm_size, OPP_PART_WORLD);

    //  MPI_Request request_send_t[pe_list->ranks_size];
    MPI_Request *request_send_t = (MPI_Request *)opp_host_malloc(pe_list->ranks_size * sizeof(MPI_Request)); 
    MPI_Request *request_send_p = (MPI_Request *)opp_host_malloc(pe_list->ranks_size * sizeof(MPI_Request));

    int *rbuf_t, *rbuf_p;
    cap = 0;
    count = 0;

    for (int i = 0; i < pe_list->ranks_size; i++) 
    {
        int *sbuf_t = &pe_list->list[pe_list->disps[i]];
        int *sbuf_p = &part_list_e[pe_list->disps[i]];
        MPI_Isend(sbuf_t, pe_list->sizes[i], MPI_INT, pe_list->ranks[i], 1, OPP_PART_WORLD, &request_send_t[i]);
        MPI_Isend(sbuf_p, pe_list->sizes[i], MPI_INT, pe_list->ranks[i], 2, OPP_PART_WORLD, &request_send_p[i]);
    }

    for (int i = 0; i < ranks_size; i++)
        cap = cap + sizes[i];
    int *temp_list_t = (int *)opp_host_malloc(cap * sizeof(int));
    part_list_i = (int *)opp_host_malloc(cap * sizeof(int));

    for (int i = 0; i < ranks_size; i++) 
    {
        rbuf_t = (int *)opp_host_malloc(sizes[i] * sizeof(int));
        rbuf_p = (int *)opp_host_malloc(sizes[i] * sizeof(int));

        MPI_Recv(rbuf_t, sizes[i], MPI_INT, neighbors[i], 1, OPP_PART_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(rbuf_p, sizes[i], MPI_INT, neighbors[i], 2, OPP_PART_WORLD, MPI_STATUS_IGNORE);

        memcpy(&temp_list_t[count], (void *)&rbuf_t[0], sizes[i] * sizeof(int));
        memcpy(&part_list_i[count], (void *)&rbuf_p[0], sizes[i] * sizeof(int));

        count = count + sizes[i];
        opp_host_free(rbuf_t);
        opp_host_free(rbuf_p);
    }
    MPI_Waitall(pe_list->ranks_size, request_send_t, MPI_STATUSES_IGNORE);
    MPI_Waitall(pe_list->ranks_size, request_send_p, MPI_STATUSES_IGNORE);

    create_imp_list_2(map->to, temp_list_t, pi_list, count, neighbors, sizes, ranks_size, comm_size, my_rank);

    //-----go through local mapping table as well as the imported information
    // and partition the "to" set
    cap = map->to->size;
    count = 0;
    int *to_elems = (int *)opp_host_malloc(sizeof(int) * cap);
    int *parts = (int *)opp_host_malloc(sizeof(int) * cap);

    //--first the local mapping table
    int local_index;
    int part;
    for (int i = 0; i < map->from->size; i++) 
    {
        for (int j = 0; j < map->dim; j++) 
        {
            part = get_partition(map->map[i * map->dim + j], part_range[map->to->index], &local_index, comm_size);
            if (part == my_rank) 
            {
                if (count >= cap) 
                {
                    cap = cap * 2;
                    parts = (int *)opp_host_realloc(parts, sizeof(int) * cap);
                    to_elems = (int *)opp_host_realloc(to_elems, sizeof(int) * cap);
                }
                to_elems[count] = local_index;
                parts[count++] = p_set->elem_part[i];
            }
        }
    }

    // copy pi_list.list and part_list_i to to_elems and parts
    if (count + pi_list->size > 0) {
        to_elems = (int *)opp_host_realloc(to_elems, sizeof(int) * (count + pi_list->size));
        parts = (int *)opp_host_realloc(parts, sizeof(int) * (count + pi_list->size));
    }

    memcpy(&to_elems[count], (void *)&pi_list->list[0], pi_list->size * sizeof(int));
    memcpy(&parts[count], (void *)&part_list_i[0], pi_list->size * sizeof(int));

    int *partition = (int *)opp_host_malloc(sizeof(int) * map->to->size);
    for (int i = 0; i < map->to->size; i++) {
        partition[i] = -99;
    }

    count = count + pi_list->size;

    // sort both to_elems[] and correspondingly parts[] arrays
    if (count > 0)
        quickSort_2(to_elems, parts, 0, count - 1);

    if (count > comm_size * 10) 
    {
        int *part_counter = (int *)opp_host_malloc(comm_size * sizeof(int));
        for (int i = 0; i < count;) 
        {
            memset(part_counter, 0, comm_size * sizeof(int));
            int curr = to_elems[i];
            do 
            {
                part_counter[parts[i]]++;
                i++;

                if (i >= count)
                    break;
            } while (curr == to_elems[i]);
            
            int maxpos = 0;
            for (int j = 0; j < comm_size; j++)
                if (part_counter[maxpos] < part_counter[j])
                    maxpos = j;
            partition[curr] = maxpos;
        }
        opp_host_free(part_counter);
    } 
    else 
    {
        int *found_parts;
        for (int i = 0; i < count;) 
        {
            int curr = to_elems[i];
            int c = 0;
            cap = map->dim;
            found_parts = (int *)opp_host_malloc(sizeof(int) * cap);

            do 
            {
                if (c >= cap) 
                {
                    cap = cap * 2;
                    found_parts = (int *)opp_host_realloc(found_parts, sizeof(int) * cap);
                }
                found_parts[c++] = parts[i];
                i++;
                
                if (i >= count)
                    break;
            } while (curr == to_elems[i]);

            partition[curr] = find_mode(found_parts, c);
            opp_host_free(found_parts);
        }
    }

    if (count + pi_list->size > 0) 
    {
        opp_host_free(to_elems);
        opp_host_free(parts);
    }

    // check if this "from" set is an "on to" set
    // need to check this globally on all processors
    int ok = 1;
    for (int i = 0; i < map->to->size; i++) 
    {
        if (partition[i] < 0) 
        {
            printf("on rank %d: Map %s is not an an on-to mapping from set %s to set %s\n",
                    my_rank, map->name, map->from->name, map->to->name);

            // return -1;
            ok = -1;
            break;
        }
    }

    // check if globally this map was giving us an on-to set mapping
    int *global_ok_array = (int *)opp_host_malloc(comm_size * sizeof(int));
    for (int r = 0; r < comm_size; r++)
        global_ok_array[r] = 1;

    MPI_Allgather(&ok, 1, MPI_INT, global_ok_array, 1, MPI_INT, OPP_PART_WORLD);
    int result = 1;
    for (int r = 0; r < comm_size; r++) 
    {
        if (global_ok_array[r] < 0) 
        {
            // printf("Rank %d reported problem partitioning\n",r);
            result = -1;
        }
    }
    opp_host_free(global_ok_array);

    if (result == 1) 
    {
        OPP_part_list[map->to->index]->elem_part = partition;
        OPP_part_list[map->to->index]->is_partitioned = 1;
    } 
    else 
    {
        opp_host_free(partition);
    }

    // cleanup
    opp_host_free(pi_list->list);
    opp_host_free(pi_list->ranks);
    opp_host_free(pi_list->sizes);
    opp_host_free(pi_list->disps);
    opp_host_free(pi_list);
    opp_host_free(pe_list->list);
    opp_host_free(pe_list->ranks);
    opp_host_free(pe_list->sizes);
    opp_host_free(pe_list->disps);
    opp_host_free(pe_list);
    opp_host_free(part_list_i);
    opp_host_free(part_list_e);

    opp_host_free(request_send_p);
    opp_host_free(request_send_t);

    return result;
}

/*******************************************************************************
 * Routine to partition all secondary sets using primary set partition
 *******************************************************************************/

/* static */ void partition_all(opp_set primary_set, int my_rank, int comm_size) 
{
    if (OPP_DBG) opp_printf("partition_all", "start");

    // Compute global partition range information for each set
    int **part_range = (int **)opp_host_malloc(opp_sets.size() * sizeof(int *));
    get_part_range(part_range, my_rank, comm_size, OPP_PART_WORLD);

    int sets_partitioned = 1;
    int maps_used = 0;

    std::vector<opp_set> all_partitioned_sets(opp_sets.size());
    std::vector<int> all_used_maps(opp_maps.size(), -1);

    // add all particle sets to avoid partitioning errors
    for (size_t i = 0; i < opp_sets.size(); i++) 
    {
        opp_set set = opp_sets[i];
        if (set->is_particle)
            all_partitioned_sets[sets_partitioned++] = set;
    }

    // begin with the partitioned primary set
    all_partitioned_sets[0] = opp_sets[primary_set->index];

    int error = 0;
    while (sets_partitioned < (int)opp_sets.size() && error == 0) 
    {
        std::vector<int> cost(opp_maps.size());
        for (size_t i = 0; i < opp_maps.size(); i++)
            cost[i] = 99;

        // compute a "cost" associated with using each mapping table
        for (size_t m = 0; m < opp_maps.size(); m++) 
        {
            opp_map map = opp_maps[m];

            if (linear_search(all_used_maps.data(), map->index, 0, maps_used - 1) < 0) // if not used before
            {
                part to_set = OPP_part_list[map->to->index];
                part from_set = OPP_part_list[map->from->index];

                // partitioning a set using a mapping from a partitioned set costs
                // more than partitioning a set using a mapping to a partitioned set
                // i.e. preferance is given to the latter over the former
                if (from_set->is_partitioned == 1 && 
                    compare_all_sets(map->from, all_partitioned_sets.data(), sets_partitioned) >= 0)
                {
                    cost[map->index] = 2;
                }
                else if (to_set->is_partitioned == 1 &&
                        compare_all_sets(map->to, all_partitioned_sets.data(), sets_partitioned) >= 0)
                {
                    cost[map->index] = (map->dim == 1 ? 0 : 1);
                }
            }
        }

        while (1) 
        {
            int selected = min(cost.data(), (int)opp_maps.size());

            if (selected >= 0) 
            {
                opp_map map = opp_maps[selected];

                // partition using this map
                part to_set = OPP_part_list[map->to->index];
                part from_set = OPP_part_list[map->from->index];

                if (to_set->is_partitioned == 1 && from_set->is_partitioned == 0) 
                {
                    if (partition_from_set(map, my_rank, comm_size, part_range) > 0) 
                    {
                        all_partitioned_sets[sets_partitioned++] = map->from;
                        all_used_maps[maps_used++] = map->index;
                        break;
                    } 
                    else // partitioning unsuccessful with this map- find another map
                        cost[selected] = 99;
                } 
                else if (from_set->is_partitioned == 1 && to_set->is_partitioned == 0) 
                {
                    if (partition_to_set(map, my_rank, comm_size, part_range) > 0) 
                    {
                        all_partitioned_sets[sets_partitioned++] = map->to;
                        all_used_maps[maps_used++] = map->index;
                        break;
                    } 
                    else // partitioning unsuccessful with this map - find another map
                        cost[selected] = 99;
                } 
                else 
                {
                    cost[selected] = 99;
                }
            } 
            else // partitioning error;
            {
                opp_printf("partition_all", "On rank %d: Partitioning error", my_rank);
                error = 1;
                break;
            }
        }
    }

    if (my_rank == OPP_ROOT) 
    {
        int die = 0;
        opp_printf("partition_all", "Sets partitioned = %d", sets_partitioned);
        
        if (sets_partitioned != (int)opp_sets.size()) 
        {
            for (size_t s = 0; s < opp_sets.size(); s++) 
            { // for each set
                opp_set set = opp_sets[s];
                part P = OPP_part_list[set->index];
                if (P->is_partitioned != 1) 
                {
                    opp_printf("partition_all", "Unable to find mapping between primary set and %s", P->set->name);
                    if (P->set->size != 0)
                        die = 1;
                }
            }
            
            if (die) 
            {
                opp_printf("partition_all", "Partitioning aborted !");
                opp_abort("partition_all Partitioning aborted !");
            }
        }
    }

    for (size_t i = 0; i < opp_sets.size(); i++)
        opp_host_free(part_range[i]);
    opp_host_free(part_range);
}

/*******************************************************************************
 * Routine to renumber mapping table entries with new partition's indexes
 *******************************************************************************/
/* static */ void renumber_maps(int my_rank, int comm_size) 
{
    // get partition rage information
    int **part_range = (int **)opp_host_malloc((int)opp_sets.size() * sizeof(int *));
    get_part_range(part_range, my_rank, comm_size, OPP_PART_WORLD);

    // find elements of the "to" set thats not in this local process
    for (int m = 0; m < (int)opp_maps.size(); m++)  // for each maping table
    { 
        opp_map map = opp_maps[m];

        int cap = 1000;
        int count = 0;
        int *req_list = (int *)opp_host_malloc(cap * sizeof(int));

        for (int i = 0; i < map->from->size; i++) 
        {
            int local_index;
            for (int j = 0; j < map->dim; j++) 
            {
                local_index = binary_search(OPP_part_list[map->to->index]->g_index,
                                map->map[i * map->dim + j], 0, map->to->size - 1);

                if (count >= cap) 
                {
                    cap = cap * 2;
                    req_list = (int *)opp_host_realloc(req_list, cap * sizeof(int));
                }

                if (local_index < 0) // not in this partition
                {
                    // store the global index of the element
                    req_list[count++] = map->map[i * map->dim + j];
                }
            }
        }
        // sort and remove duplicates
        if (count > 0) 
        {
            quickSort(req_list, 0, count - 1);
            count = removeDups(req_list, count);
            req_list = (int *)opp_host_realloc(req_list, count * sizeof(int));
        }

        // do an allgather to findout how many elements that each process will
        // be requesting partition information about
        std::vector<int> recv_count(comm_size);
        MPI_Allgather(&count, 1, MPI_INT, recv_count.data(), 1, MPI_INT, OPP_PART_WORLD);

        // discover global size of these required elements
        int g_count = 0;
        for (int i = 0; i < comm_size; i++)
            g_count += recv_count[i];

        // prepare for an allgatherv
        int disp = 0;
        int *displs = (int *)opp_host_malloc(comm_size * sizeof(int));
        for (int i = 0; i < comm_size; i++) 
        {
            displs[i] = disp;
            disp = disp + recv_count[i];
        }

        // allocate memory to hold the global indexes of elements requiring
        // partition details
        int *g_index = (int *)opp_host_malloc(sizeof(int) * g_count);

        MPI_Allgatherv(req_list, count, MPI_INT, g_index, recv_count.data(), displs, MPI_INT, OPP_PART_WORLD);
        opp_host_free(req_list);

        if (g_count > 0) 
        {
            quickSort(g_index, 0, g_count - 1);
            g_count = removeDups(g_index, g_count);
            g_index = (int *)opp_host_realloc(g_index, g_count * sizeof(int));
        }

        // printf("on rank %d map %s needs set %s : before g_count = %d\n",
        //    my_rank, map->name, map->to->name, g_count);

        // go through the recieved global g_index array and see if any local
        // element's
        // partition details are requested by some foreign process
        int *exp_index = (int *)opp_host_malloc(sizeof(int) * g_count);
        int *exp_g_index = (int *)opp_host_malloc(sizeof(int) * g_count);

        int exp_count = 0;
        for (int i = 0; i < g_count; i++) 
        {
            int local_index = binary_search(OPP_part_list[map->to->index]->g_index, g_index[i], 0, map->to->size - 1);
            int global_index;
            if (local_index >= 0) 
            {
                exp_g_index[exp_count] = g_index[i];

                global_index = get_global_index(local_index, my_rank,
                                                part_range[map->to->index], comm_size);
                exp_index[exp_count++] = global_index;
            }
        }
        opp_host_free(g_index);

        // realloc exp_index, exp_g_index
        exp_index = (int *)opp_host_realloc(exp_index, sizeof(int) * exp_count);
        exp_g_index = (int *)opp_host_realloc(exp_g_index, sizeof(int) * exp_count);

        // now export to every MPI rank, these partition info with an all-to-all
        MPI_Allgather(&exp_count, 1, MPI_INT, recv_count.data(), 1, MPI_INT, OPP_PART_WORLD);
        disp = 0;
        opp_host_free(displs);
        displs = (int *)opp_host_malloc(comm_size * sizeof(int));

        for (int i = 0; i < comm_size; i++) 
        {
            displs[i] = disp;
            disp = disp + recv_count[i];
        }

        // allocate memory to hold the incomming partition details and allgatherv
        g_count = 0;
        for (int i = 0; i < comm_size; i++)
            g_count += recv_count[i];
        
        int *all_imp_index = (int *)opp_host_malloc(sizeof(int) * g_count);
        g_index = (int *)opp_host_malloc(sizeof(int) * g_count);


        MPI_Allgatherv(exp_g_index, exp_count, MPI_INT, g_index, recv_count.data(), displs, MPI_INT, OPP_PART_WORLD);

        MPI_Allgatherv(exp_index, exp_count, MPI_INT, all_imp_index, recv_count.data(), displs, MPI_INT, OPP_PART_WORLD);

        opp_host_free(exp_index);
        opp_host_free(exp_g_index);

        // sort all_imp_index according to g_index array
        if (g_count > 0)
            quickSort_2(g_index, all_imp_index, 0, g_count - 1);

        // now we hopefully have all the informattion required to renumber this map
        // so now, again go through each entry of this mapping table and renumber
        for (int i = 0; i < map->from->size; i++) 
        {
            int local_index, global_index;
            for (int j = 0; j < map->dim; j++) 
            {
                local_index = binary_search(OPP_part_list[map->to->index]->g_index,
                                map->map[i * map->dim + j], 0, map->to->size - 1);

                if (local_index < 0) // not in this partition
                {
                    // need to search through g_index array
                    int found = binary_search(g_index, map->map[i * map->dim + j], 0,
                                                g_count - 1);
                    if (found < 0)
                        printf("Problem in renumbering rank %d map %s index %d d %d\n", 
                            OPP_rank, map->name, i, j);
                    else 
                    {
                        opp_maps[map->index]->map[i * map->dim + j] = all_imp_index[found];
                    }
                } 
                else // in this partition
                {
                    global_index = get_global_index(local_index, my_rank, part_range[map->to->index], comm_size);
                    opp_maps[map->index]->map[i * map->dim + j] = global_index;
                }
            }
        }

        opp_host_free(g_index);
        opp_host_free(displs);
        opp_host_free(all_imp_index);
    }
    for (int i = 0; i < (int)opp_sets.size(); i++)
        opp_host_free(part_range[i]);
    opp_host_free(part_range);

    if (OPP_DBG) opp_printf("partition_all", "end");
}

/*******************************************************************************
 * Routine to perform data migration to new partitions
 *******************************************************************************/
/* static */ void migrate_all(int my_rank, int comm_size) 
{
    if (OPP_DBG) opp_printf("migrate_all", "start");

    /*--STEP 1 - Create Imp/Export Lists for reverse migrating elements ----------*/

    // create imp/exp lists for reverse migration
    std::vector<halo_list> pe_list(opp_sets.size()); // export list for each set
    std::vector<halo_list> pi_list(opp_sets.size()); // import list for each set

    // create partition export lists
    int *temp_list;
    size_t count, cap;

    for (int s = 0; s < (int)opp_sets.size(); s++) 
    { // for each set
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        part p = OPP_part_list[set->index];

        // create a temporaty scratch space to hold export list for this set's partition information
        count = 0;
        cap = 1000;
        temp_list = (int *)opp_host_malloc(cap * sizeof(int));

        for (int i = 0; i < set->size; i++) 
        {
            if (p->elem_part[i] != my_rank) 
            {
                if (count >= cap) 
                {
                    cap = cap * 2;
                    temp_list = (int *)opp_host_realloc(temp_list, cap * sizeof(int));
                }
                temp_list[count++] = p->elem_part[i];
                temp_list[count++] = i; // part.g_index[i];
            }
        }
        // create partition export list
        pe_list[set->index] = (halo_list)opp_host_malloc(sizeof(halo_list_core));
        create_export_list(set, temp_list, pe_list[set->index], count, comm_size, my_rank);
        opp_host_free(temp_list);
    }

    // create partition import lists
    int *neighbors, *sizes;
    int ranks_size;

    for (int s = 0; s < (int)opp_sets.size(); s++)  // for each set
    { 
        opp_set set = opp_sets[s];
        
        if (set->is_particle) continue;
        
        halo_list exp = pe_list[set->index];

        //-----Discover neighbors-----
        ranks_size = 0;
        neighbors = (int *)opp_host_malloc(comm_size * sizeof(int));
        sizes = (int *)opp_host_malloc(comm_size * sizeof(int));

        find_neighbors_set(exp, neighbors, sizes, &ranks_size, my_rank, comm_size, OPP_PART_WORLD);

        //    MPI_Request request_send[exp->ranks_size];
        MPI_Request *request_send = (MPI_Request *)opp_host_malloc(exp->ranks_size * sizeof(MPI_Request));

        int *rbuf;
        cap = 0;
        count = 0;

        for (int i = 0; i < exp->ranks_size; i++) 
        {
            // printf("export from %d to %d set %10s, list of size %d \n",
            // my_rank,exp->ranks[i],set->name,exp->sizes[i]);
            int *sbuf = &exp->list[exp->disps[i]];
            MPI_Isend(sbuf, exp->sizes[i], MPI_INT, exp->ranks[i], 1, OPP_PART_WORLD, &request_send[i]);
        }

        for (int i = 0; i < ranks_size; i++)
            cap = cap + sizes[i];
        temp_list = (int *)opp_host_malloc(cap * sizeof(int));

        for (int i = 0; i < ranks_size; i++) 
        {
            // printf("import from %d to %d set %10s, list of size %d\n",
            // neighbors[i], my_rank, set->name, sizes[i]);
            rbuf = (int *)opp_host_malloc(sizes[i] * sizeof(int));

            MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], 1, OPP_PART_WORLD, MPI_STATUS_IGNORE);
            memcpy(&temp_list[count], (void *)&rbuf[0], sizes[i] * sizeof(int));
            count = count + sizes[i];
            opp_host_free(rbuf);
        }

        MPI_Waitall(exp->ranks_size, request_send, MPI_STATUSES_IGNORE);
        pi_list[set->index] = (halo_list)opp_host_malloc(sizeof(halo_list_core));
        create_import_list(set, temp_list, pi_list[set->index], count, neighbors, sizes, ranks_size, comm_size, my_rank);

        opp_host_free(request_send);
    }

    /*--STEP 2 - Perform Partitioning Data migration -----------------------------*/

    // data migration first ......
    for (int s = 0; s < (int)opp_sets.size(); s++)  // for each set
    {    
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        halo_list imp = pi_list[set->index];
        halo_list exp = pe_list[set->index];

        //    MPI_Request request_send[exp->ranks_size];
        MPI_Request *request_send = (MPI_Request *)opp_host_malloc(exp->ranks_size * sizeof(MPI_Request));

        // migrate data defined on this set
        int d = -1; // d is just simply the tag for mpi comms
        for (int k = 0; k < (int)opp_dats.size(); k++)
        {
            d++; // increase tag to do mpi comm for the next opp_dat
            opp_dat dat = opp_dats[k];

            if (compare_sets(dat->set, set) == 1)  // this data array is defined on this set
            {
                // prepare bits of the data array to be exported
                char **sbuf = (char **)opp_host_malloc(exp->ranks_size * sizeof(char *));

                for (int i = 0; i < exp->ranks_size; i++) 
                {
                    sbuf[i] = (char *)opp_host_malloc(exp->sizes[i] * (size_t)dat->size);
                    for (int j = 0; j < exp->sizes[i]; j++) 
                    {
                        int index = exp->list[exp->disps[i] + j];
                        memcpy(&sbuf[i][j * (size_t)dat->size], (void *)&dat->data[(size_t)dat->size * (index)], dat->size);
                    }

                    if ((size_t)dat->size * exp->sizes[i] > (size_t)INT_MAX) 
                        printf("Integer overflow at %s: %d\n",__FILE__,__LINE__);
                    MPI_Isend(sbuf[i], (size_t)dat->size * exp->sizes[i], MPI_CHAR, exp->ranks[i], d, OPP_PART_WORLD, &request_send[i]);
                }

                char *rbuf = (char *)opp_host_malloc((size_t)dat->size * imp->size);
                for (int i = 0; i < imp->ranks_size; i++) 
                {
                    if ((size_t)dat->size * imp->sizes[i] > (size_t)INT_MAX) 
                        printf("Integer overflow at %s: %d\n",__FILE__,__LINE__);
                    MPI_Recv(&rbuf[(size_t)imp->disps[i] * (size_t)dat->size], (size_t)dat->size * imp->sizes[i],MPI_CHAR, imp->ranks[i], d, 
                        OPP_PART_WORLD, MPI_STATUS_IGNORE);
                }

                MPI_Waitall(exp->ranks_size, request_send, MPI_STATUSES_IGNORE);
                for (int i = 0; i < exp->ranks_size; i++)
                    opp_host_free(sbuf[i]);
                opp_host_free(sbuf);

                // delete the data entirs that has been sent and create a
                // modified data array
                char *new_dat = (char *)opp_host_malloc((size_t)dat->size * (set->size + imp->size));

                count = 0;
                for (int i = 0; i < dat->set->size; i++) // iterate over old set size
                {
                    if (OPP_part_list[set->index]->elem_part[i] == my_rank) 
                    {
                        memcpy(&new_dat[count * (size_t)dat->size], (void *)&dat->data[(size_t)dat->size * i], dat->size);
                        count++;
                    }
                }

                memcpy(&new_dat[count * (size_t)dat->size], (void *)rbuf, (size_t)dat->size * (size_t)imp->size);
                count = count + imp->size;
                new_dat = (char *)opp_host_realloc(new_dat, (size_t)dat->size * count);
                opp_host_free(rbuf);

                opp_host_free(dat->data);
                dat->data = new_dat;
            }
        }

        opp_host_free(request_send);
    }

    // mapping tables second ......
    for (int s = 0; s < (int)opp_sets.size(); s++)  // for each set
    {
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        halo_list imp = pi_list[set->index];
        halo_list exp = pe_list[set->index];

        //    MPI_Request request_send[exp->ranks_size];
        MPI_Request *request_send = (MPI_Request *)opp_host_malloc(exp->ranks_size * sizeof(MPI_Request));

        // migrate mapping tables from this set
        for (int m = 0; m < (int)opp_maps.size(); m++)  // for each maping table
        {
            opp_map map = opp_maps[m];

            if (compare_sets(map->from, set) == 1) // need to select  mappings FROM this set
            { 
                // prepare bits of the mapping tables to be exported
                int **sbuf = (int **)opp_host_malloc(exp->ranks_size * sizeof(int *));

                // send mapping table entirs to relevant mpi processes
                for (int i = 0; i < exp->ranks_size; i++) 
                {
                    sbuf[i] = (int *)opp_host_malloc(exp->sizes[i] * map->dim * sizeof(int));
                    for (int j = 0; j < exp->sizes[i]; j++) 
                    {
                        for (int p = 0; p < map->dim; p++) 
                        {
                            sbuf[i][j * map->dim + p] = map->map[map->dim * (exp->list[exp->disps[i] + j]) + p];
                        }
                    }
                    MPI_Isend(sbuf[i], map->dim * exp->sizes[i], MPI_INT, exp->ranks[i], m, OPP_PART_WORLD, &request_send[i]);
                }

                int *rbuf = (int *)opp_host_malloc(map->dim * sizeof(int) * imp->size);

                // receive mapping table entirs from relevant mpi processes
                for (int i = 0; i < imp->ranks_size; i++) 
                {
                    MPI_Recv(&rbuf[(size_t)imp->disps[i] * map->dim], map->dim * imp->sizes[i], MPI_INT, imp->ranks[i], m, OPP_PART_WORLD, MPI_STATUS_IGNORE);
                }

                MPI_Waitall(exp->ranks_size, request_send, MPI_STATUSES_IGNORE);
                for (int i = 0; i < exp->ranks_size; i++)
                    opp_host_free(sbuf[i]);
                opp_host_free(sbuf);

                // delete the mapping table entirs that has been sent and create a
                // modified mapping table
                int *new_map = (int *)opp_host_malloc(sizeof(int) * (set->size + imp->size) * map->dim);

                count = 0;
                for (int i = 0; i < map->from->size; i++)  // iterate over old size of the maping table
                { 
                    if (OPP_part_list[map->from->index]->elem_part[i] == my_rank) 
                    {
                        memcpy(&new_map[count * map->dim], (void *)&opp_maps[map->index]->map[map->dim * i], map->dim * sizeof(int));
                        count++;
                    }
                }
                memcpy(&new_map[count * map->dim], (void *)rbuf,
                    map->dim * sizeof(int) * imp->size);
                count = count + imp->size;
                new_map = (int *)opp_host_realloc(new_map, sizeof(int) * count * map->dim);

                opp_host_free(rbuf);
                opp_host_free(opp_maps[map->index]->map);
                opp_maps[map->index]->map = new_map;
            }
        }

        opp_host_free(request_send);
    }

    /*--STEP 3 - Update Partitioning Information and Sort Set Elements------------*/

    // need to exchange the original g_index
    for (int s = 0; s < (int)opp_sets.size(); s++) 
    { // for each set
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        halo_list imp = pi_list[set->index];
        halo_list exp = pe_list[set->index];

        //    MPI_Request request_send[exp->ranks_size];
        MPI_Request *request_send = (MPI_Request *)opp_host_malloc(exp->ranks_size * sizeof(MPI_Request));

        // prepare bits of the original g_index array to be exported
        int **sbuf = (int **)opp_host_malloc(exp->ranks_size * sizeof(int *));

        // send original g_index values to relevant mpi processes
        for (int i = 0; i < exp->ranks_size; i++) 
        {
            sbuf[i] = (int *)opp_host_malloc(exp->sizes[i] * sizeof(int));
            for (int j = 0; j < exp->sizes[i]; j++) 
            {
                sbuf[i][j] = OPP_part_list[set->index]->g_index[exp->list[exp->disps[i] + j]];
            }
            MPI_Isend(sbuf[i], exp->sizes[i], MPI_INT, exp->ranks[i], s, OPP_PART_WORLD, &request_send[i]);
        }

        int *rbuf = (int *)opp_host_malloc(sizeof(int) * imp->size);

        // receive original g_index values from relevant mpi processes
        for (int i = 0; i < imp->ranks_size; i++) 
        {
            MPI_Recv(&rbuf[imp->disps[i]], imp->sizes[i], MPI_INT, imp->ranks[i], s, OPP_PART_WORLD, MPI_STATUS_IGNORE);
        }

        MPI_Waitall(exp->ranks_size, request_send, MPI_STATUSES_IGNORE);
        
        for (int i = 0; i < exp->ranks_size; i++)
            opp_host_free(sbuf[i]);
        opp_host_free(sbuf);

        // delete the g_index entirs that has been sent and create a modified g_index
        int *new_g_index = (int *)opp_host_malloc(sizeof(int) * (set->size + imp->size));

        count = 0;
        for (int i = 0; i < set->size; i++) 
        { // iterate over old
                                                // size of the g_index array
            if (OPP_part_list[set->index]->elem_part[i] == my_rank) 
            {
                new_g_index[count] = OPP_part_list[set->index]->g_index[i];
                count++;
            }
        }

        memcpy(&new_g_index[count], (void *)rbuf, sizeof(int) * imp->size);
        count = count + imp->size;
        new_g_index = (int *)opp_host_realloc(new_g_index, sizeof(int) * count);
        int *new_part = (int *)opp_host_malloc(sizeof(int) * count);
        
        for (size_t i = 0; i < count; i++)
            new_part[i] = my_rank;

        opp_host_free(rbuf);
        opp_host_free(OPP_part_list[set->index]->g_index);
        opp_host_free(OPP_part_list[set->index]->elem_part);

        OPP_part_list[set->index]->elem_part = new_part;
        OPP_part_list[set->index]->g_index = new_g_index;

        opp_sets[set->index]->size = count;
        OPP_part_list[set->index]->set = opp_sets[set->index];

        opp_host_free(request_send);
    }

    // re-set values in mapping tables
    for (int m = 0; m < (int)opp_maps.size(); m++) 
    { // for each maping table
        opp_map map = opp_maps[m];

        opp_maps[map->index]->from = opp_sets[map->from->index];
        opp_maps[map->index]->to = opp_sets[map->to->index];
    }

    // re-set values in data arrays
    for (int k = 0; k < (int)opp_dats.size(); k++)
    {
        opp_dat dat = opp_dats[k];

        if (dat->set->is_particle) continue;

        dat->set = opp_sets[dat->set->index];
    }

    // finally .... need to sort for each set, data on the set and mapping tables from this set accordiing 
    // to the OPP_part_list[set.index]->g_index array values.
    for (int s = 0; s < (int)opp_sets.size(); s++) // for each set
    { 
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        // first ... data on this set
        for (int k = 0; k < (int)opp_dats.size(); k++)
        {
            opp_dat dat = opp_dats[k];

            if (compare_sets(dat->set, set) == 1) 
            {
                if (set->size > 0) 
                {
                    int *temp = (int *)opp_host_malloc(sizeof(int) * set->size);
                    memcpy(temp, (void *)OPP_part_list[set->index]->g_index, sizeof(int) * set->size);
                    quickSort_dat(temp, dat->data, 0, set->size - 1, dat->size);
                    opp_host_free(temp);
                }
            }
        }

        // second ... mapping tables
        for (int m = 0; m < (int)opp_maps.size(); m++)  // for each maping table
        { 
            opp_map map = opp_maps[m];

            if (compare_sets(map->from, set) == 1) 
            {
                if (set->size > 0) 
                {
                    int *temp = (int *)opp_host_malloc(sizeof(int) * set->size);
                    memcpy(temp, (void *)OPP_part_list[set->index]->g_index, sizeof(int) * set->size);
                    quickSort_map(temp, opp_maps[map->index]->map, 0, set->size - 1, map->dim);
                    opp_host_free(temp);
                }
            }
        }
        if (set->size > 0)
            quickSort(OPP_part_list[set->index]->g_index, 0, set->size - 1);
    }

    // cleanup : destroy pe_list, pi_list
    for (size_t s = 0; s < opp_sets.size(); s++)  // for each set
    { 
        opp_set set = opp_sets[s];

        if (set->is_particle) continue;

        opp_host_free(pe_list[set->index]->ranks);
        opp_host_free(pe_list[set->index]->disps);
        opp_host_free(pe_list[set->index]->sizes);
        opp_host_free(pe_list[set->index]->list);

        opp_host_free(pi_list[set->index]->ranks);
        opp_host_free(pi_list[set->index]->disps);
        opp_host_free(pi_list[set->index]->sizes);
        opp_host_free(pi_list[set->index]->list);
        opp_host_free(pe_list[set->index]);
        opp_host_free(pi_list[set->index]);
    }

    if (OPP_DBG) opp_printf("migrate_all", "end");
}

/*******************************************************************************
 * Routine to revert back to the original partitioning
 *******************************************************************************/
void opp_partition_destroy() 
{
    if (OPP_DBG) opp_printf("opp_partition_destroy", "start");

    // destroy OPP_part_list[]
    for (int s = 0; s < (int)opp_sets.size(); s++)  // for each set
    { 
        opp_set set = opp_sets[s];
        // if (set->is_particle) continue;
        opp_host_free(OPP_part_list[set->index]->g_index);
        opp_host_free(OPP_part_list[set->index]->elem_part);
        opp_host_free(OPP_part_list[set->index]);
    }
    opp_host_free(OPP_part_list);
    for (int i = 0; i < (int)opp_sets.size(); i++)
    {
        // if (opp_sets[i]->is_particle) continue;
        opp_host_free(orig_part_range[i]);
    }

    opp_host_free(orig_part_range);

    if (OPP_DBG) opp_printf("opp_partition_destroy", "end");
}

/*******************************************************************************
 * Wrapper routine to partition a given set Using ParMETIS PartKway()
 *******************************************************************************/
void opp_partition_kway(opp_map primary_map) 
{
    if (primary_map->from == primary_map->to)
    {
        std::cerr << "from set == to set in provided primary map: " << primary_map->name << std::endl;
        return;
    }

    if (OPP_DBG) opp_printf("opp_partition_kway", "START");

#ifdef HAVE_PARMETIS

    // declare timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2, time, max_time;

    op_timers(&cpu_t1, &wall_t1); // timer start for partitioning

    // create new communicator for partitioning
    int my_rank, comm_size;
    MPI_Comm_dup(OPP_MPI_WORLD, &OPP_PART_WORLD);
    MPI_Comm_rank(OPP_PART_WORLD, &my_rank);
    MPI_Comm_size(OPP_PART_WORLD, &comm_size);

#ifdef DEBUG
    // check if the  primary_map is an on to map from the from-set to the to-set
    if (is_onto_map(primary_map) != 1) {
        opp_printf("opp_partition_kway", "Map %s is not an onto map from set %s to set %s",
            primary_map->name, primary_map->from->name, primary_map->to->name);
        opp_abort("opp_partition_kway Map is not an onto map");
    }
#endif

    /*--STEP 0 - initialise partitioning data stauctures with the current (block) partitioning information */

    // Compute global partition range information for each set
    int **part_range = (int **)opp_host_malloc((int)opp_sets.size() * sizeof(int *));
    get_part_range(part_range, my_rank, comm_size, OPP_PART_WORLD);

    // save the original part_range for future partition reversing
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

    // allocate memory for list
    OPP_part_list = (part *)opp_host_malloc((int)opp_sets.size() * sizeof(part));

    for (int s = 0; s < (int)opp_sets.size(); s++) { // for each set
        opp_set set = opp_sets[s];
        // printf("set %s size = %d\n", set.name, set.size);
        int *g_index = (int *)opp_host_malloc(sizeof(int) * set->size);
        for (int i = 0; i < set->size; i++)
        g_index[i] =
            get_global_index(i, my_rank, part_range[set->index], comm_size);
        decl_partition(set, g_index, NULL);
    }

    /*--STEP 1 - Construct adjacency list of the to-set of the primary_map  -------*/

    // create export list
    int c = 0;
    int cap = 1000;
    int *list = (int *)opp_host_malloc(cap * sizeof(int)); // temp list

    for (int e = 0; e < primary_map->from->size; e++)  // for each maping table entry
    { 
        int part, local_index;
        for (int j = 0; j < primary_map->dim; j++) // for each element pointed at by this entry
        { 
            part = get_partition(primary_map->map[e * primary_map->dim + j], part_range[primary_map->to->index], &local_index, comm_size);
            if (c >= cap) 
            {
                cap = cap * 2;
                list = (int *)opp_host_realloc(list, cap * sizeof(int));
            }

            if (part != my_rank) 
            {
                list[c++] = part; // add to export list
                list[c++] = e;
            }
        }
    }
    halo_list exp_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));
    create_export_list(primary_map->from, list, exp_list, c, comm_size, my_rank);
    opp_host_free(list); // free temp list

    // create import list
    int *neighbors, *sizes;
    int ranks_size;

    //-----Discover neighbors-----
    ranks_size = 0;
    neighbors = (int *)opp_host_malloc(comm_size * sizeof(int));
    sizes = (int *)opp_host_malloc(comm_size * sizeof(int));

    // halo_list list = OPP_export_exec_list[set->index];
    find_neighbors_set(exp_list, neighbors, sizes, &ranks_size, my_rank,
                        comm_size, OPP_PART_WORLD);

    MPI_Request *request_send = (MPI_Request *)opp_host_malloc(exp_list->ranks_size * sizeof(MPI_Request));

    int *rbuf, index = 0;
    cap = 0;

    for (int i = 0; i < exp_list->ranks_size; i++) 
    {
        int *sbuf = &exp_list->list[exp_list->disps[i]];
        MPI_Isend(sbuf, exp_list->sizes[i], MPI_INT, exp_list->ranks[i], primary_map->index, OPP_PART_WORLD, &request_send[i]);
    }

    for (int i = 0; i < ranks_size; i++)
        cap = cap + sizes[i];
    int *temp = (int *)opp_host_malloc(cap * sizeof(int));

    // import this list from those neighbors
    for (int i = 0; i < ranks_size; i++) 
    {
        rbuf = (int *)opp_host_malloc(sizes[i] * sizeof(int));
        MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], primary_map->index, OPP_PART_WORLD, MPI_STATUS_IGNORE);
        memcpy(&temp[index], (void *)&rbuf[0], sizes[i] * sizeof(int));
        index = index + sizes[i];
        opp_host_free(rbuf);
    }
    MPI_Waitall(exp_list->ranks_size, request_send, MPI_STATUSES_IGNORE);

    halo_list imp_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));
    create_import_list(primary_map->from, temp, imp_list, index, neighbors, sizes,
                        ranks_size, comm_size, my_rank);

    //
    // Exchange mapping table entries using the import/export lists
    //

    // prepare bits of the mapping tables to be exported
    int **sbuf = (int **)opp_host_malloc(exp_list->ranks_size * sizeof(int *));

    for (int i = 0; i < exp_list->ranks_size; i++) 
    {
        sbuf[i] = (int *)opp_host_malloc(exp_list->sizes[i] * primary_map->dim * sizeof(int));
        for (int j = 0; j < exp_list->sizes[i]; j++) 
        {
            for (int p = 0; p < primary_map->dim; p++) 
            {
                sbuf[i][j * primary_map->dim + p] = primary_map->map[primary_map->dim * (exp_list->list[exp_list->disps[i] + j]) + p];
            }
        }
        MPI_Isend(sbuf[i], primary_map->dim * exp_list->sizes[i], MPI_INT, exp_list->ranks[i], primary_map->index, 
            OPP_PART_WORLD, &request_send[i]);
    }

    // prepare space for the incomming mapping tables
    int *foreign_maps = (int *)opp_host_malloc(primary_map->dim * (imp_list->size) * sizeof(int));

    for (int i = 0; i < imp_list->ranks_size; i++) {
        MPI_Recv(&foreign_maps[(size_t)imp_list->disps[i] * primary_map->dim], primary_map->dim * imp_list->sizes[i], MPI_INT, 
                imp_list->ranks[i], primary_map->index, OPP_PART_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Waitall(exp_list->ranks_size, request_send, MPI_STATUSES_IGNORE);
    for (int i = 0; i < exp_list->ranks_size; i++)
        opp_host_free(sbuf[i]);
    opp_host_free(sbuf);

    int **adj = (int **)opp_host_malloc(primary_map->to->size * sizeof(int *));
    int *adj_i = (int *)opp_host_malloc(primary_map->to->size * sizeof(int));
    int *adj_cap = (int *)opp_host_malloc(primary_map->to->size * sizeof(int));

    for (int i = 0; i < primary_map->to->size; i++)
        adj_i[i] = 0;
    for (int i = 0; i < primary_map->to->size; i++)
        adj_cap[i] = primary_map->dim;
    for (int i = 0; i < primary_map->to->size; i++)
        adj[i] = (int *)opp_host_malloc(adj_cap[i] * sizeof(int));

    // go through each from-element of local primary_map and construct adjacency list
    for (int i = 0; i < primary_map->from->size; i++)
    {
        int part, local_index;
        for (int j = 0; j < primary_map->dim; j++)  // for each element pointed at by this entry
        { 
            part = get_partition(primary_map->map[i * primary_map->dim + j],
                                part_range[primary_map->to->index], &local_index,
                                comm_size);

            if (part == my_rank) 
            {
                for (int k = 0; k < primary_map->dim; k++) 
                {
                    if (adj_i[local_index] >= adj_cap[local_index]) 
                    {
                        adj_cap[local_index] = adj_cap[local_index] * 2;
                        adj[local_index] = (int *)opp_host_realloc(adj[local_index], adj_cap[local_index] * sizeof(int));
                    }
                    adj[local_index][adj_i[local_index]++] = primary_map->map[i * primary_map->dim + k];
                }
            }
        }
    }

    // go through each from-element of foreign primary_map and add to adjacency list
    for (int i = 0; i < imp_list->size; i++) 
    {
        int part, local_index;
        for (int j = 0; j < primary_map->dim; j++) 
        { // for each element pointed at by this entry
            part = get_partition(foreign_maps[i * primary_map->dim + j], part_range[primary_map->to->index], &local_index, comm_size);

            if (part == my_rank) 
            {
                for (int k = 0; k < primary_map->dim; k++) 
                {
                    if (adj_i[local_index] >= adj_cap[local_index]) 
                    {
                        adj_cap[local_index] = adj_cap[local_index] * 2;
                        adj[local_index] = (int *)opp_host_realloc(adj[local_index], adj_cap[local_index] * sizeof(int));
                    }
                    adj[local_index][adj_i[local_index]++] = foreign_maps[i * primary_map->dim + k];
                }
            }
        }
    }
    opp_host_free(foreign_maps);

    //
    // Setup data structures for ParMetis PartKway
    //
    idx_t comm_size_pm = comm_size;

    idx_t *vtxdist = (idx_t *)opp_host_malloc(sizeof(idx_t) * (comm_size + 1));
    for (int i = 0; i < comm_size; i++) 
    {
        vtxdist[i] = part_range[primary_map->to->index][2 * i];
    }
    vtxdist[comm_size] = part_range[primary_map->to->index][2 * (comm_size - 1) + 1] + 1;

    idx_t *xadj = (idx_t *)opp_host_malloc(sizeof(idx_t) * (primary_map->to->size + 1));
    cap = (primary_map->to->size) * primary_map->dim;

    idx_t *adjncy = (idx_t *)opp_host_malloc(sizeof(idx_t) * cap);
    int count = 0;
    int prev_count = 0;
    for (int i = 0; i < primary_map->to->size; i++) 
    {
        int g_index = get_global_index(i, my_rank, part_range[primary_map->to->index], comm_size);
        quickSort(adj[i], 0, adj_i[i] - 1);
        adj_i[i] = removeDups(adj[i], adj_i[i]);

        if (adj_i[i] < 2) 
        {
            opp_printf("opp_partition_kway", "The from set: %s of primary map: %s is not an on to set of to-set: %s",
                    primary_map->from->name, primary_map->name, primary_map->to->name);
            opp_printf("opp_partition_kway", "Need to select a different primary map");
            opp_abort("opp_partition_kway Need to select a different primary map");
        }

        adj[i] = (int *)opp_host_realloc(adj[i], adj_i[i] * sizeof(int));
        for (int j = 0; j < adj_i[i]; j++) 
        {
            if (adj[i][j] != g_index) 
            {
                if (count >= cap) 
                {
                    cap = cap * 2;
                    adjncy = (idx_t *)opp_host_realloc(adjncy, sizeof(idx_t) * cap);
                }
                adjncy[count++] = (idx_t)adj[i][j];
            }
        }
        if (i != 0) 
        {
            xadj[i] = prev_count;
            prev_count = count;
        } 
        else 
        {
            xadj[i] = 0;
            prev_count = count;
        }
    }
    xadj[primary_map->to->size] = count;

    for (int i = 0; i < primary_map->to->size; i++)
        opp_host_free(adj[i]);
    opp_host_free(adj_i);
    opp_host_free(adj_cap);
    opp_host_free(adj);

    idx_t *partition_pm = (idx_t *)opp_host_malloc(sizeof(idx_t) * primary_map->to->size);
    for (int i = 0; i < primary_map->to->size; i++) 
    {
        partition_pm[i] = -99;
    }

    idx_t edge_cut = 0;
    idx_t numflag = 0;
    idx_t wgtflag = 0;
    idx_t options[3] = {1, 3, 15};

    int *hybrid_flags = (int *)opp_host_malloc(comm_size * sizeof(int));
    for (int i = 0; i < comm_size; i++) {
        hybrid_flags[i] = 0;
    }
    // MPI_Allgather(&OPP_hybrid_gpu, 1, MPI_INT, hybrid_flags, 1, MPI_INT, OPP_PART_WORLD);
    double total = 0;
    for (int i = 0; i < comm_size; i++)
        total += hybrid_flags[i] == 1 ? OPP_hybrid_balance : 1.0;

    idx_t ncon = 1;
    real_t *tpwgts = (real_t *)opp_host_malloc(comm_size * sizeof(real_t) * ncon);
    for (int i = 0; i < comm_size * ncon; i++)
        tpwgts[i] = hybrid_flags[i] == 1 ? OPP_hybrid_balance / total : 1.0 / total;

    real_t *ubvec = (real_t *)opp_host_malloc(sizeof(real_t) * ncon);
    *ubvec = 1.05;

    // clean up before calling ParMetis
    for (int i = 0; i < (int)opp_sets.size(); i++)
        opp_host_free(part_range[i]);
    opp_host_free(part_range);
    opp_host_free(imp_list->list);
    opp_host_free(imp_list->disps);
    opp_host_free(imp_list->ranks);
    opp_host_free(imp_list->sizes);
    opp_host_free(exp_list->list);
    opp_host_free(exp_list->disps);
    opp_host_free(exp_list->ranks);
    opp_host_free(exp_list->sizes);
    opp_host_free(imp_list);
    opp_host_free(exp_list);

    if (my_rank == OPP_ROOT) 
    {
        opp_printf("opp_partition_kway", "-----------------------------------------------------------");
        opp_printf("opp_partition_kway", "ParMETIS_V3_PartKway Output");
        opp_printf("opp_partition_kway", "-----------------------------------------------------------");
    }

    ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, NULL, NULL, &wgtflag, &numflag, &ncon, &comm_size_pm, tpwgts, ubvec, 
        options, &edge_cut, partition_pm, &OPP_PART_WORLD);
        
    if (my_rank == OPP_ROOT)
        opp_printf("opp_partition_kway", "-----------------------------------------------------------");
    opp_host_free(vtxdist);
    opp_host_free(xadj);
    opp_host_free(adjncy);
    opp_host_free(ubvec);
    opp_host_free(tpwgts);

    int *partition = (int *)opp_host_malloc(sizeof(int) * primary_map->to->size);
    for (int i = 0; i < primary_map->to->size; i++) 
    {
        // sanity check to see if all elements were partitioned
        if (partition_pm[i] < 0 || partition_pm[i] >= comm_size) 
        {
            opp_printf("opp_partition_kway", "Partitioning problem: on rank %d, set %s element %d not assigned a partition",
                    my_rank, primary_map->to->name, i);
            opp_abort("opp_partition_kway Partitioning problem");
        }
        partition[i] = partition_pm[i];
    }
    opp_host_free(partition_pm);

    // initialise primary set as partitioned
    OPP_part_list[primary_map->to->index]->elem_part = partition;
    OPP_part_list[primary_map->to->index]->is_partitioned = 1;

    /*-STEP 2 - Partition all other sets,migrate data and renumber mapping tables-*/

    partition_all(primary_map->to, my_rank, comm_size);

    migrate_all(my_rank, comm_size);

    renumber_maps(my_rank, comm_size);

    op_timers(&cpu_t2, &wall_t2); // timer stop for partitioning

    time = wall_t2 - wall_t1;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, OPP_ROOT, OPP_PART_WORLD);
    MPI_Comm_free(&OPP_PART_WORLD);
  
    if (my_rank == OPP_ROOT) opp_printf("opp_partition_kway", "Max total Kway partitioning time = %lf", max_time);

    opp_host_free(request_send);
#endif

    if (OPP_DBG) opp_printf("opp_partition_kway", "end");
}

/*****************************************************************************************************************************************
 * This routine partitions based on information contained in an opp_dat called
 *partvecXXXX (number of total partitions, padded with 0s)
 *****************************************************************************************************************************************/
void opp_partition_external(opp_set primary_set, opp_dat partvec) 
{
    if (OPP_DBG) opp_printf("opp_partition_external", "start");

    // declare timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    double time;
    double max_time;

    op_timers(&cpu_t1, &wall_t1); // timer start for partitioning

    // create new communicator for partitioning
    int my_rank, comm_size;
    MPI_Comm_dup(OPP_MPI_WORLD, &OPP_PART_WORLD);
    MPI_Comm_rank(OPP_PART_WORLD, &my_rank);
    MPI_Comm_size(OPP_PART_WORLD, &comm_size);

    /*--STEP 0 - initialise partitioning data stauctures with the current (block) partitioning information */

    // Compute global partition range information for each set
    int **part_range = (int **)opp_host_malloc((int)opp_sets.size() * sizeof(int *));
    get_part_range(part_range, my_rank, comm_size, OPP_PART_WORLD);

    // save the original part_range for future partition reversing
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

    // allocate memory for list
    OPP_part_list = (part *)opp_host_malloc((int)opp_sets.size() * sizeof(part));

    for (int s = 0; s < (int)opp_sets.size(); s++) { // for each set
        opp_set set = opp_sets[s];
        int *g_index = (int *)opp_host_malloc(sizeof(int) * set->size);
        for (int i = 0; i < set->size; i++)
        g_index[i] =
            get_global_index(i, my_rank, part_range[set->index], comm_size);
        decl_partition(set, g_index, NULL);
    }

    int *partition = (int *)opp_host_malloc(sizeof(int) * primary_set->size);
    memcpy(partition, partvec->data, sizeof(int) * primary_set->size);

    // initialise primary set as partitioned
    OPP_part_list[primary_set->index]->elem_part = partition;
    OPP_part_list[primary_set->index]->is_partitioned = 1;

    // free part range
    for (int i = 0; i < (int)opp_sets.size(); i++)
        opp_host_free(part_range[i]);
    opp_host_free(part_range);

    /*-STEP 2 - Partition all other sets,migrate data and renumber mapping tables-*/

    // partition all other sets
    partition_all(primary_set, my_rank, comm_size);

    // migrate data, sort elements
    migrate_all(my_rank, comm_size);

    // renumber mapping tables
    renumber_maps(my_rank, comm_size);

    op_timers(&cpu_t2, &wall_t2); // timer stop for partitioning

    if (OPP_DBG) opp_printf("opp_partition_external", "done - calculating time");

    // print time for partitioning
    time = wall_t2 - wall_t1;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, OPP_ROOT, OPP_PART_WORLD);
    MPI_Comm_free(&OPP_PART_WORLD);
    if (my_rank == OPP_ROOT)
       opp_printf("opp_partition_external", "Max total random partitioning time = %lf", max_time);

    if (OPP_DBG) opp_printf("opp_partition_external", "end");
}


/*******************************************************************************
 * Wrapper routine to use ParMETIS_V3_PartGeom() which partitions a set
 * Using its XYZ Geometry Data
 *******************************************************************************/
void opp_partition_geom(opp_dat coords) 
{
    opp_printf("opp_partition_geom", "start");

#ifdef HAVE_PARMETIS

    // declare timers
    double cpu_t1, cpu_t2, wall_t1, wall_t2;
    double time;
    double max_time;

    op_timers(&cpu_t1, &wall_t1); // timer start for partitioning

    // create new communicator for partitioning
    int my_rank, comm_size;
    MPI_Comm_dup(OPP_MPI_WORLD, &OPP_PART_WORLD);
    MPI_Comm_rank(OPP_PART_WORLD, &my_rank);
    MPI_Comm_size(OPP_PART_WORLD, &comm_size);

    /*--STEP 0 - initialise partitioning data stauctures with the current (block)
        partitioning information */

    // Compute global partition range information for each set
    int **part_range = (int **)opp_host_malloc((int)opp_sets.size() * sizeof(int *));
    get_part_range(part_range, my_rank, comm_size, OPP_PART_WORLD);

    // save the original part_range for future partition reversing
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

    // allocate memory for list
    OPP_part_list = (part *)opp_host_malloc((int)opp_sets.size() * sizeof(part));

    for (int s = 0; s < (int)opp_sets.size(); s++)  // for each set
    { 
        opp_set set = opp_sets[s];
        // printf("set %s size = %d\n", set.name, set.size);
        int *g_index = (int *)opp_host_malloc(sizeof(int) * set->size);
        for (int i = 0; i < set->size; i++)
            g_index[i] = get_global_index(i, my_rank, part_range[set->index], comm_size);
        decl_partition(set, g_index, NULL);
    }

    /*--- STEP 1 - Partition primary set using its coordinates (1D,2D or 3D) -----*/

    // Setup data structures for ParMetis PartGeom
    idx_t *vtxdist = (idx_t *)opp_host_malloc(sizeof(idx_t) * (comm_size + 1));
    idx_t *partition = (idx_t *)opp_host_malloc(sizeof(idx_t) * coords->set->size);

    idx_t ndims = coords->dim;
    real_t *xyz = 0;

    // Create ParMetis compatible coordinates array
    //- i.e. coordinates should be floats
    if (ndims == 3 || ndims == 2 || ndims == 1) 
    {
        xyz = (real_t *)opp_host_malloc(coords->set->size * coords->dim * sizeof(real_t));
        size_t mult = coords->size / coords->dim;
        for (int i = 0; i < coords->set->size; i++) 
        {
            double temp;
            for (int e = 0; e < coords->dim; e++) 
            {
                memcpy(&temp, (void *)&(coords->data[(i * coords->dim + e) * mult]), mult);
                xyz[i * coords->dim + e] = (real_t)temp;
            }
        }
    } 
    else 
    {
        printf("Dimensions of Coordinate array not one of 3D,2D or 1D\n");
        printf("Not supported by ParMetis - Indicate correct coordinates array\n");
        opp_abort("opp_partition_geom Not supported by ParMetis");
    }

    for (int i = 0; i < comm_size; i++) {
        vtxdist[i] = part_range[coords->set->index][2 * i];
    }
    vtxdist[comm_size] =
        part_range[coords->set->index][2 * (comm_size - 1) + 1] + 1;

    // use xyz coordinates to feed into ParMETIS_V3_PartGeom
    ParMETIS_V3_PartGeom(vtxdist, &ndims, xyz, partition, &OPP_PART_WORLD);
    opp_host_free(xyz);
    opp_host_free(vtxdist);

    // free part range
    for (int i = 0; i < (int)opp_sets.size(); i++)
        opp_host_free(part_range[i]);
    opp_host_free(part_range);

    // saniti check to see if all elements were partitioned
    for (int i = 0; i < coords->size; i++) 
    {
        if (partition[i] < 0) 
        {
            printf("Partitioning problem: on rank %d, set %s element %d not assigned a partition\n",
                    my_rank, coords->name, i);
            opp_abort("opp_partition_geom Partitioning problem");
        }
    }

    // initialise primary set as partitioned
    OPP_part_list[coords->set->index]->elem_part = (int *)partition;
    OPP_part_list[coords->set->index]->is_partitioned = 1;

    /*-STEP 2 - Partition all other sets,migrate data and renumber mapping tables-*/

    // partition all other sets
    partition_all(coords->set, my_rank, comm_size);

    // migrate data, sort elements
    migrate_all(my_rank, comm_size);

    // renumber mapping tables
    renumber_maps(my_rank, comm_size);

    op_timers(&cpu_t2, &wall_t2); // timer stop for partitioning

    // printf time for partitioning
    time = wall_t2 - wall_t1;
    MPI_Reduce(&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, OPP_ROOT, OPP_PART_WORLD);
    MPI_Comm_free(&OPP_PART_WORLD);
    if (my_rank == OPP_ROOT)
        printf("Max total geometric partitioning time = %lf\n", max_time);

#endif
    
    opp_printf("opp_partition_geom", "end");
}

#ifdef __cplusplus
}
#endif
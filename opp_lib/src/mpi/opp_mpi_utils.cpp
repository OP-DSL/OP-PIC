
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

MPI_Comm OPP_MPI_IO_WORLD;

extern const char fmt_double[] = " %2.25lE";
extern const char fmt_float[] = " %f";
extern const char fmt_int[] = " %d";

void _mpi_allgather(int *l, int *g, int size, int *recevcnts, int *displs, MPI_Comm comm) 
{
    MPI_Allgatherv(l, size, MPI_INT, g, recevcnts, displs, MPI_INT, comm);
}

void _mpi_allgather(float *l, float *g, int size, int *recevcnts, int *displs, MPI_Comm comm) 
{
    MPI_Allgatherv(l, size, MPI_FLOAT, g, recevcnts, displs, MPI_FLOAT, comm);
}

void _mpi_allgather(double *l, double *g, int size, int *recevcnts, int *displs, MPI_Comm comm) 
{
    MPI_Allgatherv(l, size, MPI_DOUBLE, g, recevcnts, displs, MPI_DOUBLE, comm);
}

void _mpi_gather(int *l, int *g, int size, int *recevcnts, int *displs, MPI_Comm comm) 
{
    MPI_Gatherv(l, size, MPI_INT, g, recevcnts, displs, MPI_INT, OPP_ROOT, comm);
}

void _mpi_gather(float *l, float *g, int size, int *recevcnts, int *displs, MPI_Comm comm) 
{
    MPI_Gatherv(l, size, MPI_FLOAT, g, recevcnts, displs, MPI_FLOAT, OPP_ROOT, comm);
}

void _mpi_gather(double *l, double *g, int size, int *recevcnts, int *displs, MPI_Comm comm) 
{
    MPI_Gatherv(l, size, MPI_DOUBLE, g, recevcnts, displs, MPI_DOUBLE, OPP_ROOT, comm);
}

void checked_write(int v, const char *file_name) 
{
    if (v) 
    {
        opp_printf("checked_write", "Error: error writing to %s\n", file_name);
        opp_abort("checked_write");
    }
}

template <typename T, const char *fmt>
void write_txt(FILE *fp, int g_size, int elem_size, T *g_array, const char *file_name) 
{
    checked_write(fprintf(fp, "%d %d\n", g_size, elem_size) < 0, file_name);

    for (int i = 0; i < g_size; i++) 
    {
        for (int j = 0; j < elem_size; j++)
        checked_write(fprintf(fp, fmt, g_array[i * elem_size + j]) < 0, file_name);
                        
        fprintf(fp, "\n");
    }
}

template <typename T, void (*F)(FILE *, int, int, T *, const char *)>
void write_file(opp_dat dat, const char *file_name) 
{
    // create new communicator for output
    int rank, comm_size;
    MPI_Comm_dup(OPP_MPI_WORLD, &OPP_MPI_IO_WORLD);
    MPI_Comm_rank(OPP_MPI_IO_WORLD, &rank);
    MPI_Comm_size(OPP_MPI_IO_WORLD, &comm_size);

    // compute local number of elements in dat
    int count = dat->set->size;

    T *l_array = (T *)opp_host_malloc(dat->dim * (count) * sizeof(T));
    memcpy(l_array, (void *)&(dat->data[0]), (size_t)dat->size * count);

    int l_size = count;
    int elem_size = dat->dim;
    int *recevcnts = (int *)opp_host_malloc(comm_size * sizeof(int));
    int *displs = (int *)opp_host_malloc(comm_size * sizeof(int));
    int disp = 0;
    T *g_array = 0;

    MPI_Allgather(&l_size, 1, MPI_INT, recevcnts, 1, MPI_INT, OPP_MPI_IO_WORLD);

    int g_size = 0;
    for (int i = 0; i < comm_size; i++) 
    {
        g_size += recevcnts[i];
        recevcnts[i] = elem_size * recevcnts[i];
    }

    for (int i = 0; i < comm_size; i++) 
    {
        displs[i] = disp;
        disp = disp + recevcnts[i];
    }

    if (rank == OPP_ROOT)
        g_array = (T *)opp_host_malloc(elem_size * g_size * sizeof(T));

    _mpi_gather(l_array, g_array, l_size * elem_size, recevcnts, displs, OPP_MPI_IO_WORLD);

    if (rank == OPP_ROOT) 
    {
        FILE *fp;
        if ((fp = fopen(file_name, "w")) == NULL) 
        {
            printf("can't open file %s\n", file_name);
            opp_abort("write_file can't open file");
        }

        // Write binary or text as requested by the caller
        F(fp, g_size, elem_size, g_array, file_name);

        fclose(fp);
        opp_host_free(g_array);
    }

    opp_host_free(l_array);
    opp_host_free(recevcnts);
    opp_host_free(displs);

    MPI_Comm_free(&OPP_MPI_IO_WORLD);
}

void print_dat_to_txtfile_mpi(opp_dat dat, const char *file_name) 
{
    if (strcmp(dat->type, "double") == 0)
        write_file<double, write_txt<double, fmt_double> >(dat, file_name);
    
    else if (strcmp(dat->type, "float") == 0)
        write_file<float, write_txt<float, fmt_float> >(dat, file_name);
    
    else if (strcmp(dat->type, "int") == 0)
        write_file<int, write_txt<int, fmt_int> >(dat, file_name);
    
    else
        printf("Unknown type %s, cannot be written to file %s\n", dat->type,
            file_name);
}


opp_dat opp_mpi_get_data(opp_dat dat) 
{
    // create new communicator for fetching
    int my_rank, comm_size;
    MPI_Comm_rank(OPP_MPI_WORLD, &my_rank);
    MPI_Comm_size(OPP_MPI_WORLD, &comm_size);

    // make a copy of the distributed opp_dat on to a distributed temporary opp_dat
    opp_dat temp_dat = (opp_dat)opp_host_malloc(sizeof(opp_dat_core));
    char *data = (char *)opp_host_malloc((size_t)dat->set->size * dat->size);
    memcpy(data, dat->data, dat->set->size * (size_t)dat->size);

    // use orig_part_range to fill in OPP_part_list[set->index]->elem_part with
    // original partitioning information
    for (int i = 0; i < dat->set->size; i++)
    {
        int local_index;
        OPP_part_list[dat->set->index]->elem_part[i] = get_partition(OPP_part_list[dat->set->index]->g_index[i],
            orig_part_range[dat->set->index], &local_index, comm_size);
    }

    halo_list pe_list;
    halo_list pi_list;

    // create export list
    part p = OPP_part_list[dat->set->index];
    int count = 0;
    int cap = 1000;
    int *temp_list = (int *)opp_host_malloc(cap * sizeof(int));

    for (int i = 0; i < dat->set->size; i++) 
    {
        if (p->elem_part[i] != my_rank) 
        {
            if (count >= cap) 
            {
                cap = cap * 2;
                temp_list = (int *)opp_host_realloc(temp_list, cap * sizeof(int));
            }
            temp_list[count++] = p->elem_part[i];
            temp_list[count++] = i;
        }
    }

    pe_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));
    
    create_export_list(dat->set, temp_list, pe_list, count, comm_size, my_rank);
    
    opp_host_free(temp_list);

    // create import list

    int *neighbors, *sizes;
    int ranks_size;

    //-----Discover neighbors-----
    ranks_size = 0;
    neighbors = (int *)opp_host_malloc(comm_size * sizeof(int));
    sizes = (int *)opp_host_malloc(comm_size * sizeof(int));

    find_neighbors_set(pe_list, neighbors, sizes, &ranks_size, my_rank, comm_size, OPP_MPI_WORLD);
    
    std::vector<MPI_Request> request_send(pe_list->ranks_size);

    int *rbuf;
    cap = 0;
    count = 0;

    for (int i = 0; i < pe_list->ranks_size; i++) 
    {
        int *sbuf = &pe_list->list[pe_list->disps[i]];
        MPI_Isend(sbuf, pe_list->sizes[i], MPI_INT, pe_list->ranks[i], 1, OPP_MPI_WORLD, &request_send[i]);
    }

    for (int i = 0; i < ranks_size; i++)
        cap = cap + sizes[i];
    temp_list = (int *)opp_host_malloc(cap * sizeof(int));

    for (int i = 0; i < ranks_size; i++) 
    {
        rbuf = (int *)opp_host_malloc(sizes[i] * sizeof(int));
        MPI_Recv(rbuf, sizes[i], MPI_INT, neighbors[i], 1, OPP_MPI_WORLD, MPI_STATUS_IGNORE);

        memcpy(&temp_list[count], (void *)&rbuf[0], sizes[i] * sizeof(int));
        count = count + sizes[i];
        
        opp_host_free(rbuf);
    }

    MPI_Waitall(pe_list->ranks_size, request_send.data(), MPI_STATUSES_IGNORE);
    
    pi_list = (halo_list)opp_host_malloc(sizeof(halo_list_core));
    
    create_import_list(dat->set, temp_list, pi_list, count, neighbors, sizes,
                        ranks_size, comm_size, my_rank);

    // migrate the temp "data" array to the original MPI ranks

    // prepare bits of the data array to be exported
    char **sbuf_char = (char **)opp_host_malloc(pe_list->ranks_size * sizeof(char *));

    for (int i = 0; i < pe_list->ranks_size; i++) 
    {
        sbuf_char[i] = (char *)opp_host_malloc((size_t)pe_list->sizes[i] * (size_t)dat->size);
        
        for (int j = 0; j < pe_list->sizes[i]; j++) 
        {
            int index = pe_list->list[pe_list->disps[i] + j];
            memcpy(&sbuf_char[i][j * (size_t)dat->size], (void *)&data[(size_t)dat->size * (index)],
                    dat->size);
        }

        MPI_Isend(sbuf_char[i], (size_t)dat->size * pe_list->sizes[i], MPI_CHAR,
                pe_list->ranks[i], dat->index, OPP_MPI_WORLD, &request_send[i]);
    }

    char *rbuf_char = (char *)opp_host_malloc((size_t)dat->size * pi_list->size);
    for (int i = 0; i < pi_list->ranks_size; i++) 
    {
        MPI_Recv(&rbuf_char[pi_list->disps[i] * (size_t)dat->size],
                (size_t)dat->size * pi_list->sizes[i], MPI_CHAR, pi_list->ranks[i],
                dat->index, OPP_MPI_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Waitall(pe_list->ranks_size, request_send.data(), MPI_STATUSES_IGNORE);
    
    for (int i = 0; i < pe_list->ranks_size; i++)
        opp_host_free(sbuf_char[i]);
    opp_host_free(sbuf_char);

    // delete the data entirs that has been sent and create a modified data array
    char *new_dat = (char *)opp_host_malloc((size_t)dat->size * (dat->set->size + pi_list->size));

    count = 0;
    for (int i = 0; i < dat->set->size; i++) // iterate over old set size
    {
        if (OPP_part_list[dat->set->index]->elem_part[i] == my_rank) 
        {
            memcpy(&new_dat[count * (size_t)dat->size], (void *)&data[(size_t)dat->size * i], dat->size);
            count++;
        }
    }

    memcpy(&new_dat[count * (size_t)dat->size], (void *)rbuf_char, (size_t)dat->size * pi_list->size);
    count = count + pi_list->size;
    new_dat = (char *)opp_host_realloc(new_dat, (size_t)dat->size * count);
    
    opp_host_free(rbuf_char);
    opp_host_free(data);
    data = new_dat;

    // make a copy of the original g_index and migrate that also to the original MPI process
    // prepare bits of the original g_index array to be exported
    int **sbuf = (int **)opp_host_malloc(pe_list->ranks_size * sizeof(int *));

    // send original g_index values to relevant mpi processes
    for (int i = 0; i < pe_list->ranks_size; i++) 
    {
        sbuf[i] = (int *)opp_host_malloc(pe_list->sizes[i] * sizeof(int));
        for (int j = 0; j < pe_list->sizes[i]; j++) 
        {
            sbuf[i][j] = OPP_part_list[dat->set->index]->g_index[pe_list->list[pe_list->disps[i] + j]];
        }

        MPI_Isend(sbuf[i], pe_list->sizes[i], MPI_INT, pe_list->ranks[i], dat->index, 
            OPP_MPI_WORLD, &request_send[i]);
    }

    rbuf = (int *)opp_host_malloc(sizeof(int) * pi_list->size);

    // receive original g_index values from relevant mpi processes
    for (int i = 0; i < pi_list->ranks_size; i++) 
    {
        MPI_Recv(&rbuf[pi_list->disps[i]], pi_list->sizes[i], MPI_INT, pi_list->ranks[i], dat->index, 
            OPP_MPI_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Waitall(pe_list->ranks_size, request_send.data(), MPI_STATUSES_IGNORE);
    
    for (int i = 0; i < pe_list->ranks_size; i++)
        opp_host_free(sbuf[i]);
    opp_host_free(sbuf);

    // delete the g_index entirs that has been sent and create a modified g_index
    int *new_g_index = (int *)opp_host_malloc(sizeof(int) * (dat->set->size + pi_list->size));

    count = 0;
    for (int i = 0; i < dat->set->size; i++)  // iterate over old size of the g_index array
    { 
        if (OPP_part_list[dat->set->index]->elem_part[i] == my_rank) 
        {
            new_g_index[count] = OPP_part_list[dat->set->index]->g_index[i];
            count++;
        }
    }

    memcpy(&new_g_index[count], (void *)rbuf, sizeof(int) * pi_list->size);
    count = count + pi_list->size;
    new_g_index = (int *)opp_host_realloc(new_g_index, sizeof(int) * count);
    opp_host_free(rbuf);

    // sort elements in temporaty data according to new_g_index
    quickSort_dat(new_g_index, data, 0, count - 1, dat->size);

    // cleanup
    opp_host_free(pe_list->ranks);
    opp_host_free(pe_list->disps);
    opp_host_free(pe_list->sizes);
    opp_host_free(pe_list->list);
    opp_host_free(pe_list);
    opp_host_free(pi_list->ranks);
    opp_host_free(pi_list->disps);
    opp_host_free(pi_list->sizes);
    opp_host_free(pi_list->list);
    opp_host_free(pi_list);
    opp_host_free(new_g_index);

    // remember that the original set size is now given by count
    opp_set set = (opp_set)opp_host_malloc(sizeof(opp_set_core));
    set->index = dat->set->index;
    set->size = count;
    set->name = dat->set->name;

    temp_dat->index = dat->index;
    temp_dat->set = set;
    temp_dat->dim = dat->dim;
    temp_dat->data = data;
    temp_dat->data_d = NULL;
    temp_dat->name = dat->name;
    temp_dat->type = dat->type;
    temp_dat->size = dat->size;

    return temp_dat;
}


/*******************************************************************************
 * Routine to declare partition information for a given set
 *******************************************************************************/
void decl_partition(opp_set set, int *g_index, int *partition) 
{
    part p = (part)opp_host_malloc(sizeof(part_core));
    p->set = set;
    p->g_index = g_index;
    p->elem_part = partition;
    p->is_partitioned = 0;
    OPP_part_list[set->index] = p;
    OPP_part_index++;
}

/*******************************************************************************
 * Routine to get partition range on all mpi ranks for all sets
 *******************************************************************************/
void get_part_range(int **part_range, int my_rank, int comm_size, MPI_Comm Comm) 
{
    (void)my_rank;
    for (int s = 0; s < (int)opp_sets.size(); s++) 
    {
        opp_set set = opp_sets[s];

        int *sizes = (int *)opp_host_malloc(sizeof(int) * comm_size);
        MPI_Allgather(&set->size, 1, MPI_INT, sizes, 1, MPI_INT, Comm);

        part_range[set->index] = (int *)opp_host_malloc(2 * comm_size * sizeof(int));

        int disp = 0;
        for (int i = 0; i < comm_size; i++) 
        {
            part_range[set->index][2 * i] = disp;
            disp = disp + sizes[i] - 1;
            part_range[set->index][2 * i + 1] = disp;
            disp++;
        #ifdef DEBUG
            if (my_rank == OPP_ROOT && OPP_diags > 5)
                printf("range of %10s in rank %d: %d-%d\n", set->name, i,
                    part_range[set->index][2 * i],
                    part_range[set->index][2 * i + 1]);
        #endif
        }
        opp_host_free(sizes);
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

    if (OPP_DBG) 
    {    
        std::string log = "";
        for (int i = 0; i < comm_size; i++) 
        {
            log += std::to_string(i) + "|" + std::to_string(part_range[2 * i]) + "|" + 
                    std::to_string(part_range[2 * i + 1]) + "\n";
        }

        opp_printf("get_partition()", "Error: orphan global index %d part_range->\n%s", 
            global_index, log.c_str());
    }

    opp_abort("get_partition");
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
    if (g_index > part_range[2 * (comm_size - 1) + 1] && OPP_diags > 2)
        printf("Global index larger than set size\n");
#endif
    return g_index;
}

/*******************************************************************************
 * Routine to find the MPI neighbors given a halo list
 *******************************************************************************/
void find_neighbors_set(halo_list List, int *neighbors, int *sizes, int *ranks_size, int my_rank, int comm_size, MPI_Comm Comm) 
{
    std::vector<int> temp(comm_size, -99);
    std::vector<int> r_temp(comm_size * comm_size, -99);

    for (int r = 0; r < comm_size; r++) 
    {
        if (List->ranks[r] >= 0)
            temp[List->ranks[r]] = List->sizes[r];
    }

    MPI_Allgather(temp.data(), comm_size, MPI_INT, r_temp.data(), comm_size, MPI_INT, Comm);

    int n = 0;
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
}

bool is_double_indirect_reduction(opp_arg& arg)
{
    if (arg.argtype == OPP_ARG_DAT && arg.p2c_map != nullptr && 
            arg.idx != -1 && arg.acc == OPP_INC)
        return true;
    
    return false;
}
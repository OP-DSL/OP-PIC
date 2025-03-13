
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

//*******************************************************************************
void opp_init_double_indirect_reductions(int nargs, opp_arg *args)
{
    if (OPP_DBG) opp_printf("opp_init_double_indirect_reductions", "ALL START");

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
                char* reset_values = nullptr;
                switch (args[n].acc)
                {
                    case OPP_INC:
                    case OPP_MAX:
                        if (strcmp(args[n].dat->type, "double") == 0)
                            reset_values = (char*)opp_zero_double16;
                        else if (strcmp(args[n].dat->type, "int") == 0)
                            reset_values = (char*)opp_zero_int16;
                        else
                            opp_printf("opp_init_double_indirect_reductions", 
                                "Error: Datatype [%s] dat [%s] not implemented",
                                args[n].dat->type, args[n].dat->name);
                        break;
                    default:
                        opp_printf("opp_init_double_indirect_reductions", 
                            "Error: Reduction of dat [%s] operator not implemented",
                            args[n].dat->name);
                        return;
                }

                // reset the import non execute halos
                opp_reset_dat_impl(args[n].dat, reset_values, OPP_Reset_inh);
            }
        }
    }

    if (OPP_DBG) opp_printf("opp_init_double_indirect_reductions", "ALL END");
}

//*******************************************************************************
void opp_exchange_double_indirect_reductions(int nargs, opp_arg *args) 
{
    if (OPP_DBG) opp_printf("opp_exchange_double_indirect_reductions", "ALL START");

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
                opp_reduc_comm reduc_comm = OPP_Reduc_NO_Comm;
                switch (args[n].acc)
                {
                    case OPP_INC: reduc_comm = OPP_Reduc_SUM_Comm; break;
                    case OPP_MAX: reduc_comm = OPP_Reduc_MAX_Comm; break;
                    case OPP_MIN: reduc_comm = OPP_Reduc_MIN_Comm; break;
                }

                opp_exchange_double_indirect_reductions(args[n].dat, reduc_comm);
            }
        }
    }

    if (OPP_DBG) opp_printf("opp_exchange_double_indirect_reductions", "ALL END");
}

//*******************************************************************************
void opp_exchange_double_indirect_reductions(opp_dat dat, opp_reduc_comm reduc_comm) 
{
    if (OPP_DBG) opp_printf("opp_exchange_double_indirect_reductions", "START dat %s", dat->name);

    halo_list imp_nonexec_list = OPP_import_nonexec_list[dat->set->index];
    halo_list exp_nonexec_list = OPP_export_nonexec_list[dat->set->index];
    op_mpi_buffer reduc_buf    = (op_mpi_buffer)(dat->mpi_reduc_buffer);

    opp_profiler->startMpiComm("", opp::OPP_Mesh);

    double total_send_size = 0.0;

    // TODO : Include execute as well

    // send reduction contributions related to import non-execute lists
    int init = (dat->set->size + dat->set->exec_size) * dat->size;
    for (int i = 0; i < imp_nonexec_list->ranks_size; i++) 
    { 
        int send_rank    = imp_nonexec_list->ranks[i];
        char* send_buf   = &(dat->data[init + imp_nonexec_list->disps[i] * dat->size]);
        size_t send_size    = (size_t)dat->size * imp_nonexec_list->sizes[i];
        MPI_Request* req = &(reduc_buf->s_req[reduc_buf->s_num_req++]);

        if (OPP_DBG) 
            opp_printf("opp_exchange_double_indirect_reductions", "SEND SIZE %zu bytes to %d", 
                send_size, send_rank);

        MPI_Isend(send_buf, send_size, MPI_CHAR, send_rank, dat->index, OPP_MPI_WORLD, req);

        total_send_size += (send_size * 1.0f);
    }

    opp_profiler->addTransferSize("", opp::OPP_Mesh, total_send_size);

    // receive reduction contributions related to export non-execute lists
    for (int i = 0; i < exp_nonexec_list->ranks_size; i++) 
    {
        int recv_rank    = exp_nonexec_list->ranks[i];
        char* recv_buf   = &(reduc_buf->buf_nonexec[exp_nonexec_list->disps[i] * dat->size]);
        size_t recv_size    = (size_t)dat->size * exp_nonexec_list->sizes[i];
        MPI_Request* req = &(reduc_buf->r_req[reduc_buf->r_num_req++]);

        if (OPP_DBG) 
            opp_printf("opp_exchange_double_indirect_reductions", "RECEIVE SIZE %zu bytes from %d", 
                recv_size, recv_rank);

        MPI_Irecv(recv_buf, recv_size, MPI_CHAR, recv_rank, dat->index, OPP_MPI_WORLD, req);
    }

    dat->reduc_comm = reduc_comm;

    opp_profiler->endMpiComm("", opp::OPP_Mesh);

    if (OPP_DBG) opp_printf("opp_exchange_double_indirect_reductions", "END dat %s", dat->name);
}

//*******************************************************************************
void opp_complete_double_indirect_reductions(int nargs, opp_arg *args) 
{
    if (OPP_DBG) opp_printf("opp_complete_double_indirect_reductions", "ALL START");

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
                opp_complete_double_indirect_reductions(args[n].dat);
        }
    }

    if (OPP_DBG) opp_printf("opp_complete_double_indirect_reductions", "ALL END");
}

//*******************************************************************************
void opp_complete_double_indirect_reductions(opp_dat dat) 
{
    if (dat->reduc_comm == OPP_Reduc_NO_Comm)
    {
        //if (OPP_DBG) 
        opp_printf("opp_complete_double_indirect_reductions", 
            "No reduction communication in flight for dat %s", dat->name);
        return;
    }

    if (OPP_DBG) opp_printf("opp_complete_double_indirect_reductions", "START dat %s", dat->name);

    opp_profiler->startMpiComm("", opp::OPP_Mesh);

    op_mpi_buffer reduc_buf = (op_mpi_buffer)(dat->mpi_reduc_buffer);

    MPI_Waitall(reduc_buf->s_num_req, reduc_buf->s_req, MPI_STATUSES_IGNORE);
    MPI_Waitall(reduc_buf->r_num_req, reduc_buf->r_req, MPI_STATUSES_IGNORE);

    opp_profiler->endMpiComm("", opp::OPP_Mesh);

    reduc_buf->s_num_req = 0;
    reduc_buf->r_num_req = 0;

    halo_list exp_nonexec_list = OPP_export_nonexec_list[dat->set->index];

    opp_reduc_comm reduc_comm = dat->reduc_comm;
    int dim = dat->dim;
    double *d_recv_buf = nullptr, *d_dat_buf = nullptr;
    int *i_recv_buf = nullptr, *i_dat_buf = nullptr;

    for (int i = 0; i < exp_nonexec_list->ranks_size; i++) 
    {
        int *set_indices = &(exp_nonexec_list->list[exp_nonexec_list->disps[i]]);

        if (strcmp(dat->type, "double") == 0) 
        {
            d_recv_buf = (double*)(&(reduc_buf->buf_nonexec[exp_nonexec_list->disps[i] * dat->size]));         
            d_dat_buf = (double*)dat->data;

            for (int j = 0; j < exp_nonexec_list->sizes[i]; j++)
            {
                // opp_printf("opp_complete_double_indirect_reductions", "Index %d", set_indices[j]);

                opp_reduce_dat_element<OPP_REAL>(&(d_dat_buf[dim * set_indices[j]]), &(d_recv_buf[dim * j]), 
                    dim, reduc_comm);
            }
        } 
        else if (strcmp(dat->type, "int") == 0) 
        {
            i_recv_buf = (int*)&(reduc_buf->buf_nonexec[exp_nonexec_list->disps[i] * dat->size]);
            i_dat_buf = (int*)dat->data;

            for (int j = 0; j < exp_nonexec_list->sizes[i]; j++)
            {
                opp_reduce_dat_element<OPP_INT>(&(i_dat_buf[dim * set_indices[j]]), &(i_recv_buf[j * dim]), 
                    dim, reduc_comm);
            }
        } 
        else
        {
            opp_printf("opp_complete_double_indirect_reductions", 
                "Error: Reduction for data type [%s] not implemented", dat->type);
            return;
        } 

        if (OPP_DBG)
        {
            opp_printf("opp_complete_double_indirect_reductions", "Dat %s received %d elements from rank %d", 
                dat->name, exp_nonexec_list->sizes[i], exp_nonexec_list->ranks[i]);
        }
    }

    dat->reduc_comm = OPP_Reduc_NO_Comm;

    if (OPP_DBG) opp_printf("opp_complete_double_indirect_reductions", "END dat %s", dat->name);
}


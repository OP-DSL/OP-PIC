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

#include "opp_sycl.h"

void export_halo_gather(int *list, char *dat, int copy_size,
                                   int elem_size, char *export_buffer,
                                   const sycl::nd_item<1> &item_ct1) 
{
    const int id = item_ct1.get_global_linear_id();
    if (id < copy_size) 
    {
        int off = 0;
        if (elem_size % 16 == 0) 
        {
            off += 16 * (elem_size / 16);
            for (int i = 0; i < elem_size / 16; i++) 
            {
                ((sycl::double2 *)(export_buffer + id * elem_size))[i] =
                    ((sycl::double2 *)(dat + list[id] * elem_size))[i];
            }
        } 
        else if (elem_size % 8 == 0) 
        {
            off += 8 * (elem_size / 8);
            for (int i = 0; i < elem_size / 8; i++) 
            {
                ((double *)(export_buffer + id * elem_size))[i] =
                    ((double *)(dat + list[id] * elem_size))[i];
            }
        }
        for (int i = off; i < elem_size; i++) 
        {
            export_buffer[id * elem_size + i] = dat[list[id] * elem_size + i];
        }
    }
}

void export_halo_gather_soa(int *list, char *dat, int copy_size,
                                       int elem_size, char *export_buffer,
                                       int set_size, int dim,
                                       const sycl::nd_item<1> &item_ct1) 
{
    const int id = item_ct1.get_global_linear_id();
    const int size_of = elem_size / dim;
    
    if (id < copy_size) 
    {
        if (size_of == 8) 
        {
            for (int i = 0; i < dim; i++) 
            {
                ((double *)(export_buffer + id * elem_size))[i] =
                    ((double *)(dat + list[id] * size_of))[i * set_size];
            }
        } 
        else 
        {
            for (int i = 0; i < dim; i++) 
            {
                for (int j = 0; j < size_of; j++) 
                {
                    export_buffer[id * elem_size + i * size_of + j] =
                        dat[list[id] * size_of + i * set_size * size_of + j];
                }
            }
        }
    }
}

void gather_data_to_buffer(opp_arg arg, halo_list exp_exec_list,
                           halo_list exp_nonexec_list) 
{
    if (OPP_DBG) opp_printf("opp_helper", "Running gather_data_to_buffer");

    const int threads = 192;
    const int blocks = 1 + ((exp_exec_list->size - 1) / threads);
    const int blocks2 = 1 + ((exp_nonexec_list->size - 1) / threads);

    if (strstr(arg.dat->type, ":soa") != NULL || (OPP_auto_soa && arg.dat->dim > 1)) 
    {
        const int set_size = arg.dat->set->size + arg.dat->set->exec_size +
                    arg.dat->set->nonexec_size;

        {
            dpct::has_capability_or_fail(
                opp_queue->get_device(), {sycl::aspect::fp64});

            opp_queue->submit([&](sycl::handler &cgh) {
                int *export_exec_list_d_arg_dat_set_index_ct0 =
                    export_exec_list_d[arg.dat->set->index];
                int exp_exec_list_size_ct2 = exp_exec_list->size;
                int arg_dat_size_ct3 = arg.dat->size;
                char *arg_dat_buffer_d_ct4 = arg.dat->buffer_d;
                int arg_dat_dim_ct6 = arg.dat->dim;

                cgh.parallel_for(
                    sycl::nd_range<1>(threads * blocks, threads),
                    [=](sycl::nd_item<1> item_ct1) {
                        export_halo_gather_soa(
                            export_exec_list_d_arg_dat_set_index_ct0,
                            arg.data_d, exp_exec_list_size_ct2,
                            arg_dat_size_ct3, arg_dat_buffer_d_ct4, set_size,
                            arg_dat_dim_ct6, item_ct1);
                    });
            });
        }

        {
            dpct::has_capability_or_fail(
                opp_queue->get_device(), {sycl::aspect::fp64});

            opp_queue->submit([&](sycl::handler &cgh) {
                int *export_nonexec_list_d_arg_dat_set_index_ct0 =
                    export_nonexec_list_d[arg.dat->set->index];
                int exp_nonexec_list_size_ct2 = exp_nonexec_list->size;
                int arg_dat_size_ct3 = arg.dat->size;
                char *arg_dat_buffer_d_exp_exec_list_size_arg_dat_size_ct4 =
                    arg.dat->buffer_d + exp_exec_list->size * arg.dat->size;
                int arg_dat_dim_ct6 = arg.dat->dim;

                cgh.parallel_for(
                    sycl::nd_range<1>(threads * blocks2, threads),
                    [=](sycl::nd_item<1> item_ct1) {
                        export_halo_gather_soa(
                            export_nonexec_list_d_arg_dat_set_index_ct0,
                            arg.data_d, exp_nonexec_list_size_ct2,
                            arg_dat_size_ct3,
                            arg_dat_buffer_d_exp_exec_list_size_arg_dat_size_ct4,
                            set_size, arg_dat_dim_ct6, item_ct1);
                    });
            });
        }

    } 
    else 
    {
        {
            dpct::has_capability_or_fail(
                opp_queue->get_device(), {sycl::aspect::fp64});

            opp_queue->submit([&](sycl::handler &cgh) {
                int *export_exec_list_d_arg_dat_set_index_ct0 =
                    export_exec_list_d[arg.dat->set->index];
                int exp_exec_list_size_ct2 = exp_exec_list->size;
                int arg_dat_size_ct3 = arg.dat->size;
                char *arg_dat_buffer_d_ct4 = arg.dat->buffer_d;

                cgh.parallel_for(
                    sycl::nd_range<1>(threads * blocks, threads),
                    [=](sycl::nd_item<1> item_ct1) {
                        export_halo_gather(
                            export_exec_list_d_arg_dat_set_index_ct0,
                            arg.data_d, exp_exec_list_size_ct2,
                            arg_dat_size_ct3, arg_dat_buffer_d_ct4, item_ct1);
                    });
            });
        }

        {
            dpct::has_capability_or_fail(
                opp_queue->get_device(), {sycl::aspect::fp64});

            opp_queue->submit([&](sycl::handler &cgh) {
                int *export_nonexec_list_d_arg_dat_set_index_ct0 =
                    export_nonexec_list_d[arg.dat->set->index];
                int exp_nonexec_list_size_ct2 = exp_nonexec_list->size;
                int arg_dat_size_ct3 = arg.dat->size;
                char *arg_dat_buffer_d_exp_exec_list_size_arg_dat_size_ct4 =
                    arg.dat->buffer_d + exp_exec_list->size * arg.dat->size;

                cgh.parallel_for(
                    sycl::nd_range<1>(threads * blocks2, threads),
                    [=](sycl::nd_item<1> item_ct1) {
                        export_halo_gather(
                            export_nonexec_list_d_arg_dat_set_index_ct0,
                            arg.data_d, exp_nonexec_list_size_ct2,
                            arg_dat_size_ct3,
                            arg_dat_buffer_d_exp_exec_list_size_arg_dat_size_ct4,
                            item_ct1);
                    });
            });
        }
    }

    OPP_DEVICE_SYNCHRONIZE();
}




void import_halo_scatter_soa(int offset, char *dat, int copy_size,
                                        int elem_size, char *import_buffer,
                                        int set_size, int dim,
                                        const sycl::nd_item<1> &item_ct1) 
{
    const int id = item_ct1.get_global_linear_id();
    const int size_of = elem_size / dim;
    
    if (id < copy_size) 
    {
        if (size_of == 8) 
        {
            for (int i = 0; i < dim; i++) 
            {
                ((double *)(dat + (offset + id) * size_of))[i * set_size] =
                    ((double *)(import_buffer + id * elem_size))[i];
            }
        } 
        else 
        {
            for (int i = 0; i < dim; i++) 
            {
                for (int j = 0; j < size_of; j++) 
                {
                    dat[(offset + id) * size_of + i * set_size * size_of + j] =
                        import_buffer[id * elem_size + i * size_of + j];
                }
            }
        }
    }
}

void scatter_data_from_buffer(opp_arg arg) 
{
    if (OPP_DBG) opp_printf("opp_helper", "Running scatter_data_from_buffer");
    
    const int threads = 192;
    const int blocks = 1 + ((arg.dat->set->exec_size - 1) / threads);
    const int blocks2 = 1 + ((arg.dat->set->nonexec_size - 1) / threads);

    if (strstr(arg.dat->type, ":soa") != NULL || (OPP_auto_soa && arg.dat->dim > 1)) 
    {
        const int set_size = arg.dat->set->size + arg.dat->set->exec_size +
                    arg.dat->set->nonexec_size;
        int offset = arg.dat->set->size;
        int copy_size = arg.dat->set->exec_size;

        {
            dpct::has_capability_or_fail(
                opp_queue->get_device(), {sycl::aspect::fp64});

            opp_queue->submit([&](sycl::handler &cgh) {
                int arg_dat_size_ct3 = arg.dat->size;
                char *arg_dat_buffer_d_r_ct4 = arg.dat->buffer_d_r;
                int arg_dat_dim_ct6 = arg.dat->dim;

                cgh.parallel_for(
                    sycl::nd_range<1>(threads * blocks, threads),
                    [=](sycl::nd_item<1> item_ct1) {
                        import_halo_scatter_soa(
                            offset, arg.data_d, copy_size, arg_dat_size_ct3,
                            arg_dat_buffer_d_r_ct4, set_size, arg_dat_dim_ct6,
                            item_ct1);
                    });
            });
        }

        offset += arg.dat->set->exec_size;
        copy_size = arg.dat->set->nonexec_size;

        {
            dpct::has_capability_or_fail(
                opp_queue->get_device(), {sycl::aspect::fp64});

            opp_queue->submit([&](sycl::handler &cgh) {
                int arg_dat_size_ct3 = arg.dat->size;
                char
                    *arg_dat_buffer_d_r_arg_dat_set_exec_size_arg_dat_size_ct4 =
                        arg.dat->buffer_d_r +
                        arg.dat->set->exec_size * arg.dat->size;
                int arg_dat_dim_ct6 = arg.dat->dim;

                cgh.parallel_for(
                    sycl::nd_range<1>(threads * blocks, threads),
                    [=](sycl::nd_item<1> item_ct1) {
                        import_halo_scatter_soa(
                            offset, arg.data_d, copy_size, arg_dat_size_ct3,
                            arg_dat_buffer_d_r_arg_dat_set_exec_size_arg_dat_size_ct4,
                            set_size, arg_dat_dim_ct6, item_ct1);
                    });
            });
        }
    }

    OPP_DEVICE_SYNCHRONIZE();
}
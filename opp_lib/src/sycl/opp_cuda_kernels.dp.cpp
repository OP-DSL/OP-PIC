#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <sycl/sycl.hpp>
#include <dpct/dpct.hpp>
#include "opp_sycl.h"

void export_halo_gather(int *list, char *dat, int copy_size,
                                   int elem_size, char *export_buffer,
                                   const sycl::nd_item<3> &item_ct1) 
{
    const int id = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                   item_ct1.get_local_id(2);

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
                                       const sycl::nd_item<3> &item_ct1) 
{
    const int id = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
                   item_ct1.get_local_id(2);
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
    const int threads = 192;
    const int blocks = 1 + ((exp_exec_list->size - 1) / threads);
    const int blocks2 = 1 + ((exp_nonexec_list->size - 1) / threads);

    if (strstr(arg.dat->type, ":soa") != NULL || (OPP_auto_soa && arg.dat->dim > 1)) 
    {
        const int set_size = arg.dat->set->size + arg.dat->set->exec_size +
                    arg.dat->set->nonexec_size;

        {
            dpct::has_capability_or_fail(
                dpct::get_in_order_queue().get_device(), {sycl::aspect::fp64});

            dpct::get_in_order_queue().submit([&](sycl::handler &cgh) {
                int *export_exec_list_d_arg_dat_set_index_ct0 =
                    export_exec_list_d[arg.dat->set->index];
                int exp_exec_list_size_ct2 = exp_exec_list->size;
                int arg_dat_size_ct3 = arg.dat->size;
                char *arg_dat_buffer_d_ct4 = arg.dat->buffer_d;
                int arg_dat_dim_ct6 = arg.dat->dim;

                cgh.parallel_for(
                    sycl::nd_range<3>(sycl::range<3>(1, 1, blocks) *
                                          sycl::range<3>(1, 1, threads),
                                      sycl::range<3>(1, 1, threads)),
                    [=](sycl::nd_item<3> item_ct1) {
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
                dpct::get_in_order_queue().get_device(), {sycl::aspect::fp64});

            dpct::get_in_order_queue().submit([&](sycl::handler &cgh) {
                int *export_nonexec_list_d_arg_dat_set_index_ct0 =
                    export_nonexec_list_d[arg.dat->set->index];
                int exp_nonexec_list_size_ct2 = exp_nonexec_list->size;
                int arg_dat_size_ct3 = arg.dat->size;
                char *arg_dat_buffer_d_exp_exec_list_size_arg_dat_size_ct4 =
                    arg.dat->buffer_d + exp_exec_list->size * arg.dat->size;
                int arg_dat_dim_ct6 = arg.dat->dim;

                cgh.parallel_for(
                    sycl::nd_range<3>(sycl::range<3>(1, 1, blocks2) *
                                          sycl::range<3>(1, 1, threads),
                                      sycl::range<3>(1, 1, threads)),
                    [=](sycl::nd_item<3> item_ct1) {
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
                dpct::get_in_order_queue().get_device(), {sycl::aspect::fp64});

            dpct::get_in_order_queue().submit([&](sycl::handler &cgh) {
                int *export_exec_list_d_arg_dat_set_index_ct0 =
                    export_exec_list_d[arg.dat->set->index];
                int exp_exec_list_size_ct2 = exp_exec_list->size;
                int arg_dat_size_ct3 = arg.dat->size;
                char *arg_dat_buffer_d_ct4 = arg.dat->buffer_d;

                cgh.parallel_for(
                    sycl::nd_range<3>(sycl::range<3>(1, 1, blocks) *
                                          sycl::range<3>(1, 1, threads),
                                      sycl::range<3>(1, 1, threads)),
                    [=](sycl::nd_item<3> item_ct1) {
                        export_halo_gather(
                            export_exec_list_d_arg_dat_set_index_ct0,
                            arg.data_d, exp_exec_list_size_ct2,
                            arg_dat_size_ct3, arg_dat_buffer_d_ct4, item_ct1);
                    });
            });
        }

        {
            dpct::has_capability_or_fail(
                dpct::get_in_order_queue().get_device(), {sycl::aspect::fp64});

            dpct::get_in_order_queue().submit([&](sycl::handler &cgh) {
                int *export_nonexec_list_d_arg_dat_set_index_ct0 =
                    export_nonexec_list_d[arg.dat->set->index];
                int exp_nonexec_list_size_ct2 = exp_nonexec_list->size;
                int arg_dat_size_ct3 = arg.dat->size;
                char *arg_dat_buffer_d_exp_exec_list_size_arg_dat_size_ct4 =
                    arg.dat->buffer_d + exp_exec_list->size * arg.dat->size;

                cgh.parallel_for(
                    sycl::nd_range<3>(sycl::range<3>(1, 1, blocks2) *
                                          sycl::range<3>(1, 1, threads),
                                      sycl::range<3>(1, 1, threads)),
                    [=](sycl::nd_item<3> item_ct1) {
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

    cutilSafeCall(
        DPCT_CHECK_ERROR(dpct::get_current_device().queues_wait_and_throw()));
}




void import_halo_scatter_soa(int offset, char *dat, int copy_size,
                                        int elem_size, char *import_buffer,
                                        int set_size, int dim,
                                        const sycl::nd_item<3> &item_ct1) 
{
    int id = item_ct1.get_group(2) * item_ct1.get_local_range(2) +
             item_ct1.get_local_id(2);
    int size_of = elem_size / dim;
    
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
    const int threads = 192;
    const int blocks = 1 + ((arg.dat->set->exec_size - 1) / 192);
    const int blocks2 = 1 + ((arg.dat->set->nonexec_size - 1) / 192);

    if (strstr(arg.dat->type, ":soa") != NULL || (OPP_auto_soa && arg.dat->dim > 1)) 
    {
        const int set_size = arg.dat->set->size + arg.dat->set->exec_size +
                    arg.dat->set->nonexec_size;
        int offset = arg.dat->set->size;
        int copy_size = arg.dat->set->exec_size;

        {
            dpct::has_capability_or_fail(
                dpct::get_in_order_queue().get_device(), {sycl::aspect::fp64});

            dpct::get_in_order_queue().submit([&](sycl::handler &cgh) {
                int arg_dat_size_ct3 = arg.dat->size;
                char *arg_dat_buffer_d_r_ct4 = arg.dat->buffer_d_r;
                int arg_dat_dim_ct6 = arg.dat->dim;

                cgh.parallel_for(
                    sycl::nd_range<3>(sycl::range<3>(1, 1, blocks) *
                                          sycl::range<3>(1, 1, threads),
                                      sycl::range<3>(1, 1, threads)),
                    [=](sycl::nd_item<3> item_ct1) {
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
                dpct::get_in_order_queue().get_device(), {sycl::aspect::fp64});

            dpct::get_in_order_queue().submit([&](sycl::handler &cgh) {
                int arg_dat_size_ct3 = arg.dat->size;
                char
                    *arg_dat_buffer_d_r_arg_dat_set_exec_size_arg_dat_size_ct4 =
                        arg.dat->buffer_d_r +
                        arg.dat->set->exec_size * arg.dat->size;
                int arg_dat_dim_ct6 = arg.dat->dim;

                cgh.parallel_for(
                    sycl::nd_range<3>(sycl::range<3>(1, 1, blocks2) *
                                          sycl::range<3>(1, 1, threads),
                                      sycl::range<3>(1, 1, threads)),
                    [=](sycl::nd_item<3> item_ct1) {
                        import_halo_scatter_soa(
                            offset, arg.data_d, copy_size, arg_dat_size_ct3,
                            arg_dat_buffer_d_r_arg_dat_set_exec_size_arg_dat_size_ct4,
                            set_size, arg_dat_dim_ct6, item_ct1);
                    });
            });
        }
    }
}
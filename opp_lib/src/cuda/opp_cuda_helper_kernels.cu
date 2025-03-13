#include "opp_cuda.h"

/*******************************************************************************/
__global__ void export_halo_gather(int *list, char *dat, int copy_size,
                                   int elem_size, char *export_buffer) 
{
    const int id = OPP_DEVICE_GLOBAL_LINEAR_ID;

    if (id < copy_size) {
        int off = 0;
        if (elem_size % 16 == 0) {
            const int elem_size_over_16 = (elem_size / 16);
            off += elem_size;
            for (int i = 0; i < elem_size_over_16; i++) {
                ((double2 *)(export_buffer + id * elem_size))[i] =
                    ((double2 *)(dat + list[id] * elem_size))[i];
            }
        } 
        else if (elem_size % 8 == 0) {
            const int elem_size_over_8 = (elem_size / 8);
            off += elem_size;
            for (int i = 0; i < elem_size_over_8; i++) {
                ((double *)(export_buffer + id * elem_size))[i] =
                    ((double *)(dat + list[id] * elem_size))[i];
            }
        }
        for (int i = off; i < elem_size; i++) {
            export_buffer[id * elem_size + i] = dat[list[id] * elem_size + i];
        }
    }
}
__global__ void export_halo_gather_soa(int *list, char *dat, int copy_size,
                        int elem_size, char *export_buffer, int set_size, int dim) 
{
    const int id = OPP_DEVICE_GLOBAL_LINEAR_ID;
    const int size_of = elem_size / dim;
    
    if (id < copy_size) {
        if (size_of == 8) {
            for (int i = 0; i < dim; i++) {
                ((double *)(export_buffer + id * elem_size))[i] =
                    ((double *)(dat + list[id] * size_of))[i * set_size];
            }
        } 
        else {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < size_of; j++) {
                    export_buffer[id * elem_size + i * size_of + j] =
                        dat[list[id] * size_of + i * set_size * size_of + j];
                }
            }
        }
    }
}

void gather_data_to_buffer(opp_arg arg, halo_list exp_exec_list, halo_list exp_nonexec_list) 
{
    const int blocks = 1 + ((exp_exec_list->size - 1) / opp_const_threads_per_block);
    const int blocks2 = 1 + ((exp_nonexec_list->size - 1) / opp_const_threads_per_block);

    if (strstr(arg.dat->type, ":soa") != NULL || (OPP_auto_soa && arg.dat->dim > 1)) {
        const int set_size = arg.dat->set->size + arg.dat->set->exec_size +
                    arg.dat->set->nonexec_size;

        export_halo_gather_soa<<<blocks, opp_const_threads_per_block>>>(
            export_exec_list_d[arg.dat->set->index], arg.data_d,
            exp_exec_list->size, arg.dat->size, arg.dat->buffer_d, set_size,
            arg.dat->dim);

        export_halo_gather_soa<<<blocks2, opp_const_threads_per_block>>>(
            export_nonexec_list_d[arg.dat->set->index], arg.data_d,
            exp_nonexec_list->size, arg.dat->size,
            arg.dat->buffer_d + exp_exec_list->size * arg.dat->size, set_size,
            arg.dat->dim);

    } 
    else {
        export_halo_gather<<<blocks, opp_const_threads_per_block>>>(
            export_exec_list_d[arg.dat->set->index], arg.data_d,
            exp_exec_list->size, arg.dat->size, arg.dat->buffer_d);

        export_halo_gather<<<blocks2, opp_const_threads_per_block>>>(
            export_nonexec_list_d[arg.dat->set->index], arg.data_d,
            exp_nonexec_list->size, arg.dat->size,
            arg.dat->buffer_d + exp_exec_list->size * arg.dat->size);
    }

    OPP_DEVICE_SYNCHRONIZE();
}

/*******************************************************************************/
__global__ void import_halo_scatter_soa(int offset, char *dat, int copy_size,
                        int elem_size, char *import_buffer, int set_size, int dim) 
{
    const int id = OPP_DEVICE_GLOBAL_LINEAR_ID;
    const int size_of = elem_size / dim;
    
    if (id < copy_size) {
        if (size_of == 8) {
            for (int i = 0; i < dim; i++) {
                ((double *)(dat + (offset + id) * size_of))[i * set_size] =
                    ((double *)(import_buffer + id * elem_size))[i];
            }
        } 
        else {
            for (int i = 0; i < dim; i++) {
                for (int j = 0; j < size_of; j++) {
                    dat[(offset + id) * size_of + i * set_size * size_of + j] =
                        import_buffer[id * elem_size + i * size_of + j];
                }
            }
        }
    }
}

void scatter_data_from_buffer(opp_arg arg) 
{
    const int blocks = 1 + ((arg.dat->set->exec_size - 1) / opp_const_threads_per_block);
    const int blocks2 = 1 + ((arg.dat->set->nonexec_size - 1) / opp_const_threads_per_block);

    if (strstr(arg.dat->type, ":soa") != NULL || (OPP_auto_soa && arg.dat->dim > 1)) {
        const int set_size = arg.dat->set->size + arg.dat->set->exec_size +
                    arg.dat->set->nonexec_size;
        int offset = arg.dat->set->size;
        int copy_size = arg.dat->set->exec_size;

        import_halo_scatter_soa<<<blocks, opp_const_threads_per_block>>>(
            offset, arg.data_d, copy_size, arg.dat->size, arg.dat->buffer_d_r,
            set_size, arg.dat->dim);

        offset += arg.dat->set->exec_size;
        copy_size = arg.dat->set->nonexec_size;

        import_halo_scatter_soa<<<blocks2, opp_const_threads_per_block>>>(
            offset, arg.data_d, copy_size, arg.dat->size,
            arg.dat->buffer_d_r + arg.dat->set->exec_size * arg.dat->size, set_size,
            arg.dat->dim);
    }
}

/*******************************************************************************/
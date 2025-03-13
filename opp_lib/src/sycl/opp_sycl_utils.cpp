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

#include <opp_sycl.h>

//****************************************
void opp_create_dat_device_arrays(opp_dat dat, bool create_new)
{
    if (OPP_DBG) 
        opp_printf("opp_create_dat_device_arrays", "%s %s", dat->name, dat->type);

    if (create_new) {
        opp_mem::dev_free(dat->data_d);
        opp_mem::dev_free(dat->data_swap_d);
    }

    const size_t alloc_count = (size_t)(dat->set->set_capacity * dat->size);
    dat->data_d = opp_mem::dev_malloc<char>(alloc_count);
    // if (OPP_fill_type != OPP_HoleFill_All)
        dat->data_swap_d = opp_mem::dev_malloc<char>(alloc_count);

    if (OPP_DBG) 
        opp_printf("opp_create_dat_device_arrays", "Device array of dat [%s][%p][%p] [%zubytes]", 
                dat->name, dat->data_d, dat->data_swap_d, alloc_count);
}

//********************************************************************************
// DEVICE->HOST | this invalidates what is in the HOST
void opp_download_dat(opp_dat dat) 
{
    OPP_DEVICE_SYNCHRONIZE();

    const size_t set_size = dat->set->set_capacity;

    if (strstr(dat->type, ":soa") != NULL || (OPP_auto_soa && dat->dim > 1)) {
        if (OPP_DBG) 
            opp_printf("opp_download_dat", "GPU->CPU SOA | %s to %p from %p", 
                dat->name, dat->data, dat->data_d);

        std::vector<char> tmp_data(dat->size * set_size);
        opp_mem::copy_dev_to_host<char>(tmp_data.data(), dat->data_d, set_size * dat->size);

        const int element_size = dat->size / dat->dim;
        for (int i = 0; i < dat->dim; i++) {
            for (int j = 0; j < set_size; j++) {
                for (int c = 0; c < element_size; c++) {
                    dat->data[dat->size * j + element_size * i + c] = 
                        tmp_data[element_size * i * set_size + element_size * j + c];
                }
            }
        }
    } 
    else {
        if (OPP_DBG) 
            opp_printf("opp_download_dat", "GPU->CPU NON-SOA| %s to %p from %p", 
                dat->name, dat->data, dat->data_d);
        opp_mem::copy_dev_to_host<char>(dat->data, dat->data_d, set_size * dat->size);
    }

    dat->dirty_hd = Dirty::NotDirty;
}

//*******************************************************************************
// HOST->DEVICE | this invalidates what is in the DEVICE
void opp_upload_dat(opp_dat dat)
{ 
    const size_t capacity = dat->set->set_capacity;
    const size_t dat_size = dat->size;
    const size_t dat_dim = dat->dim;

    std::vector<char> tmp_data;
    char* host_data = dat->data;

    if (strstr(dat->type, ":soa") != NULL || (OPP_auto_soa && dat_dim > 1)) {

        if (OPP_DBG) opp_printf("opp_upload_dat","CPU->GPU SOA | %s", dat->name);
        
        tmp_data.resize(dat_size * capacity);
        const size_t element_size = dat_size / dat_dim;
        
        for (size_t i = 0; i < dat_dim; i++) {
            for (size_t j = 0; j < capacity; j++) {
                for (size_t c = 0; c < element_size; c++) {
                    tmp_data[element_size * i * capacity + element_size * j + c] = 
                        dat->data[dat_size * j + element_size * i + c];         
                }
            }
        }
        host_data = tmp_data.data();
    } 
    else {
        if (OPP_DBG) opp_printf("opp_upload_dat", "CPU->GPU NON-SOA| %s", dat->name);
    }

    opp_mem::copy_host_to_dev<char>(dat->data_d, host_data, (dat_size * capacity));

    dat->dirty_hd = Dirty::NotDirty;
}

//*******************************************************************************
// HOST->DEVICE | this invalidates what is in the DEVICE
void opp_upload_map(opp_map map, bool create_new) 
{
    if (OPP_DBG) 
        opp_printf("opp_upload_map", "CPU->GPU | %s %s", map->name, create_new ? "NEW" : "COPY");

    const int set_size = map->from->size + map->from->exec_size + map->from->nonexec_size;
    std::vector<OPP_INT> tmp_map(map->dim * set_size);
    const int set_size_plus_exec = map->from->size + map->from->exec_size;
    
    for (int i = 0; i < map->dim; i++) {
        for (int j = 0; j < set_size; j++) {
            if (j >= set_size_plus_exec)
                tmp_map[i * set_size + j] = -10;
            else
                tmp_map[i * set_size + j] = map->map[map->dim * j + i];
        }
    }

    if (create_new)
        opp_mem::copy_host_to_dev<OPP_INT>(map->map_d, tmp_map.data(), (map->dim * set_size), 
                                        false, true, (map->dim * set_size));
    else
        opp_mem::copy_host_to_dev<OPP_INT>(map->map_d, tmp_map.data(), (map->dim * set_size));        
}

//*******************************************************************************
// DEVICE -> HOST
void opp_download_particle_set(opp_set particles_set, bool force_download)
{
    if (OPP_DBG) 
        opp_printf("opp_download_particle_set", "set [%s]", particles_set->name);

    OPP_DEVICE_SYNCHRONIZE();

    for (opp_dat& current_dat : *(particles_set->particle_dats)) {
        if (current_dat->data_d == NULL) {
            if (OPP_DBG) 
                opp_printf("opp_download_particle_set", "device ptr is NULL in dat [%s]", 
                            current_dat->name);
            continue;
        }
        if (current_dat->dirty_hd != Dirty::Host && !force_download) {
            if (OPP_DBG) 
                opp_printf("opp_download_particle_set", "host is not dirty in dat [%s]", 
                            current_dat->name);
            continue;
        }
        opp_download_dat(current_dat);
    }  
}

//*******************************************************************************
// HOST->DEVICE
void opp_upload_particle_set(opp_set particles_set, bool realloc)
{ 
    if (OPP_DBG) opp_printf("opp_upload_particle_set", "set [%s]", particles_set->name);

    for (opp_dat& current_dat : *(particles_set->particle_dats)) {
        if (realloc) {
            opp_create_dat_device_arrays(current_dat, realloc);
        }

        opp_upload_dat(current_dat);
    }  
}

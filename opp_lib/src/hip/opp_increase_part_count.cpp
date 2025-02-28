
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

#include "opp_hip.h"

OPP_INT* opp_saved_mesh_relation_d = nullptr;
size_t opp_saved_mesh_relation_size = 0;

//****************************************
void opp_increase_particle_count(opp_set set, const OPP_INT insert_count)
{ OPP_RETURN_IF_INVALID_PROCESS;

    opp_profiler->start("opp_inc_part_count");

    bool need_resizing = (set->set_capacity < (set->size + insert_count)) ? true : false;

    if (OPP_DBG) 
        opp_printf("opp_increase_particle_count", "need_resizing %s", need_resizing ? "YES" : "NO");

    // TODO : We should be able to do a device to device copy instead of getting to host

    if (need_resizing) {
        opp_profiler->start("opp_inc_part_count_DWN");
        opp_download_particle_set(set, true); 
        opp_profiler->end("opp_inc_part_count_DWN");
    }

    opp_profiler->start("opp_inc_part_count_INC");
    if (!opp_increase_particle_count_core(set, insert_count)) {
        opp_printf("opp_increase_particle_count", "Error at opp_increase_particle_count_core");
        opp_abort();
    }
    opp_profiler->end("opp_inc_part_count_INC");

    if (need_resizing) {

        opp_profiler->start("opp_inc_part_count_UPL");
        for (opp_dat& current_dat : *(set->particle_dats)) {
            if (OPP_DBG) 
                opp_printf("opp_increase_particle_count", "resizing dat [%s] set_capacity [%d]", 
                            current_dat->name, set->set_capacity);

            // TODO : We might be able to copy only the old data from device to device!

            opp_create_dat_device_arrays(current_dat, true);
            opp_upload_dat(current_dat);

            current_dat->dirty_hd = Dirty::NotDirty;
        }   
        opp_profiler->end("opp_inc_part_count_UPL");     
    } 

    opp_profiler->end("opp_inc_part_count");
}

//****************************************
__global__ void opp_dev_assign_mesh_relation_kernel(
    OPP_INT *__restrict mesh_relation,
    const OPP_INT *__restrict distribution,
    OPP_INT start,
    OPP_INT end,
    OPP_INT inj_start,
    OPP_INT inlet_size) 
{
    const int tid = OPP_DEVICE_GLOBAL_LINEAR_ID;

    if (tid + start < end) {    
        const int n = tid + start;
        for (int i = 0; i < inlet_size; i++) {
            if (tid < distribution[i]) {
                mesh_relation[n + inj_start] = i; 
                break;
            } 
        }  
    }
}

//****************************************
void opp_inc_part_count_with_distribution(opp_set set, OPP_INT insert_count, 
                                             opp_dat iface_dist, bool calc_new)
{ OPP_RETURN_IF_INVALID_PROCESS;

    if (OPP_DBG) 
        opp_printf("opp_inc_part_count_with_distribution", "insert_count [%d] %s", 
        insert_count, (calc_new ? "NEW" : "COPY"));

    opp_profiler->start("opp_inc_parts_with_distr");

    opp_dat mesh_rel_dat  = set->mesh_relation_dat;

    // TODO : BUG What happens if the complete particle is dirty in device?

    opp_increase_particle_count(set, insert_count);

    const int nargs1 = 2;
    opp_arg args1[nargs1];

    // if iface particle distribution is dirty in device, get it to the device
    args1[0] = opp_arg_dat(iface_dist, OPP_READ);
    args1[1] = opp_arg_dat(mesh_rel_dat, OPP_WRITE);

    const OPP_INT set_size = opp_mpi_halo_exchanges_grouped(set, nargs1, args1, Device_GPU);
    opp_mpi_halo_wait_all(nargs1, args1);

    if (set_size > 0) {
        const OPP_INT start     = 0;
        const OPP_INT end       = set->diff;
        const OPP_INT inj_start = (set->size - set->diff);

        if (end - start > 0) {
            if (calc_new) {
                if (OPP_DBG) 
                    opp_printf("opp_inc_part_count_with_distribution", 
                        "Calculating all from new");

                const int nthread = OPP_gpu_threads_per_block;
                const int nblocks = (end - start - 1) / nthread + 1;

                opp_dev_assign_mesh_relation_kernel<<<nblocks, nthread>>>(
                    (OPP_INT *) mesh_rel_dat->data_d,
                    (OPP_INT *) iface_dist->data_d,
                    start, 
                    end, 
                    inj_start,
                    iface_dist->set->size);                
            }
            else {
                const size_t copy_size = (end - start);
                OPP_INT* inj_mesh_relations_d = (OPP_INT *)mesh_rel_dat->data_d + inj_start;

                if (opp_saved_mesh_relation_d == nullptr) {
                    if (OPP_DBG) 
                        opp_printf("opp_inc_part_count_with_distribution", 
                            "Allocating saved_mesh_relation_d with size [%zu]", copy_size);

                    opp_saved_mesh_relation_size = copy_size;
                    opp_saved_mesh_relation_d = opp_mem::dev_malloc<OPP_INT>(copy_size);     

                    const int nthread = OPP_gpu_threads_per_block;
                    const int nblocks = (end - start - 1) / nthread + 1;

                    opp_dev_assign_mesh_relation_kernel<<<nblocks, nthread>>>(
                        (OPP_INT *) mesh_rel_dat->data_d,
                        (OPP_INT *) iface_dist->data_d,
                        start, 
                        end, 
                        inj_start,
                        iface_dist->set->size);

                    OPP_DEVICE_SYNCHRONIZE();

                    // save the mesh relation data for next iteration
                    opp_mem::copy_dev_to_dev<OPP_INT>(opp_saved_mesh_relation_d, 
                                                    inj_mesh_relations_d, copy_size);
                }
                else {
                    if (OPP_DBG) 
                        opp_printf("opp_inc_part_count_with_distribution", 
                            "Copying saved_mesh_relation_d with size [%zu]", copy_size);

                    if (opp_saved_mesh_relation_size != copy_size) {
                        opp_printf("opp_inc_part_count_with_distribution", 
                            "ERROR... saved_mesh_relation_size [%d] does not match with new copy size [%d]", 
                            opp_saved_mesh_relation_size, copy_size);
                    }

                    // Copy from the saved mesh relation data
                    opp_mem::copy_dev_to_dev<OPP_INT>(inj_mesh_relations_d, 
                                                opp_saved_mesh_relation_d, copy_size);
                }
            }
        }
    }

    opp_set_dirtybit_grouped(nargs1, args1, Device_GPU);
    OPP_DEVICE_SYNCHRONIZE();

    opp_profiler->end("opp_inc_parts_with_distr");
}


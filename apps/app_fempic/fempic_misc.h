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

//*********************************************
// USER WRITTEN CODE
//*********************************************

#include "fempic_misc_mesh_colour.h"

/*************************************************************************************************
 * This function will enrich the particle distribution values into if_dist_dat
 * In addition, dummy_random dat will get enriched with random values in "rand_file"
*/
inline int init_inject_distributions(opp_dat if_dist_dat, opp_dat if_area_dat, opp_dat dummy_random)
{
    if (OPP_DBG) opp_printf("init_inject_distributions", "START");

    const double plasma_den   = opp_params->get<OPP_REAL>("plasma_den");
    const double dt           = opp_params->get<OPP_REAL>("dt");
    const double ion_velocity = opp_params->get<OPP_REAL>("ion_velocity");
    const double spwt         = opp_params->get<OPP_REAL>("spwt");

    int total_inject_count = 0, max_inject_count_per_face = 0;
    double remainder = 0.0;

    OPP_RUN_ON_ROOT()
        opp_printf("init_inject_distributions", "plasma_den=%2.15lE dt=%2.25lE ion_vel=%2.15lE spwt=%2.15lE", 
            plasma_den, dt, ion_velocity, spwt);

    // find the number of particles to be injected through each inlet face and 
    // get the max injected particle count per face
    for (int faceID=0; faceID<if_area_dat->set->size; faceID++)
    {   
        remainder = 0.0; // always make remainder to zero, and if not the MPI results will change

        double num_per_sec = plasma_den * ion_velocity * ((double*)if_area_dat->data)[faceID];
        double num_real = num_per_sec * dt;
        double fnum_mp = num_real / spwt + remainder;
        int num_mp = (int)fnum_mp;
        // remainder = fnum_mp - num_mp;

        total_inject_count += num_mp;

        ((int*)if_dist_dat->data)[faceID] = total_inject_count;

        if (max_inject_count_per_face < num_mp)
            max_inject_count_per_face = num_mp;
    }

    if (OPP_DBG)
        opp_printf("init_inject_distributions", "RAND_FILE inj_count %d max_inj_count_per_face %d", 
            total_inject_count, max_inject_count_per_face);  
            
    // increase dummy random particle set size, to load the random numbers for particle injections
    opp_increase_particle_count(dummy_random->set, (total_inject_count + INJ_EXCESS));  

    int total_size = -1, fsize = -1, fdim = -1;
    FILE *fp = NULL;
    const std::string rand_file_path = opp_params->get<OPP_STRING>("rand_file");

    // read from MPI ROOT and broadcast // alternatively could use HDF5 files
    if (OPP_rank == OPP_ROOT)
    {       
        if ((fp = fopen(rand_file_path.c_str(), "r")) == NULL)
        {
            opp_printf("init_inject_distributions", "Unable to open file %s\n", 
                rand_file_path.c_str());
            opp_abort();
        }

        if (fscanf(fp, "%d %d\n", &fsize, &fdim) != 2)
        {
            opp_printf("init_inject_distributions", "Error reading file data from %s\n", 
                rand_file_path.c_str());
            opp_abort();
        }

        total_size = fsize * fdim;

        if (max_inject_count_per_face * dummy_random->dim > total_size)
        {
            opp_printf("init_inject_distributions", "dim and/or set_size issue in file %s\n", 
                rand_file_path.c_str());
            opp_abort();     
        }
    }

#ifdef USE_MPI
    // Load the whole file and bradcast
    // TODO : We can reduce communications by sending only the required size
    MPI_Bcast(&total_size, 1, MPI_INT, OPP_ROOT, MPI_COMM_WORLD);
#endif

    double* dist = new double[total_size];
    if (OPP_rank == OPP_ROOT)
    {
        for (int n = 0; n < fsize; n++)
        {
            if (fscanf(fp, " %lf %lf\n", &dist[n * 2 + 0], &dist[n * 2 + 1]) != 2) 
            {
                opp_printf("init_inject_distributions", "Error reading from %s at index %d\n", 
                    rand_file_path.c_str(), n);
                opp_abort();
            }
        }

        fclose(fp);
    }

#ifdef USE_MPI
    // Load the whole file and bradcast
    // TODO : We can reduce communications by sending only the required size
    MPI_Bcast(dist, total_size, MPI_DOUBLE, OPP_ROOT, MPI_COMM_WORLD);
#endif

    if (if_area_dat->set->size > 0) 
    {
        double* random_dat = (double *)dummy_random->data;
        int* distribution  = (int *)if_dist_dat->data;
        int iface_index = 0, part_in_face_index = 0;

        // This trouble is only because we need mpi results to match with seq and others
        for (int i = 0; i < dummy_random->set->size; i++)
        {
            if (iface_index < if_dist_dat->set->size && i >= distribution[iface_index])
            {
                iface_index++; // check whether it is j or j-1
                part_in_face_index = 0; 
            }

            random_dat[i * 2 + 0] = dist[part_in_face_index * 2 + 0];
            random_dat[i * 2 + 1] = dist[part_in_face_index * 2 + 1];

            part_in_face_index++;
        }
    }
    
    delete[] dist;

    dummy_random->dirty_hd = Dirty::Device; // To make GPU versions to download updated data!
    if_dist_dat->dirty_hd = Dirty::Device; // To make GPU versions to download updated data!
    
    if (OPP_DBG) opp_printf("init_inject_distributions", "total_inject_count %d", total_inject_count);

    return total_inject_count;
}

/*************************************************************************************************
 * This is a utility function to count the particles per cell and c_part_count dat will get enriched
*/
inline void print_per_cell_particle_counts(opp_dat c_part_count, opp_dat part_mesh_relation)     
{
    // validate the final particle counts in each cell
    opp_reset_dat(c_part_count, (char*)opp_zero_int16);
    for (int p = 0; p < part_mesh_relation->set->size; p++)
    {
        int c_index    = ((int *)part_mesh_relation->data)[p];
        ((int *)c_part_count->data)[c_index] += 1;
    }

#ifdef USE_MPI
    opp_mpi_print_dat_to_txtfile(c_part_count, "c_part_count.dat");
#else
    opp_print_dat_to_txtfile(c_part_count, "", "c_part_count.dat");
#endif
}

/*************************************************************************************************
 * This is a utility function create a global level log string using local values
 * @param max_c_ef - this is already reduced glocally by the opp_par_loop
 * @param max_n_potential - this is already reduced glocally by the opp_par_loop
 * @param local_part_count - these are local values
 * @param local_parts_injected - these are local values
 * @param local_part_removed - these are local values
 * @return std::string
*/
inline std::string get_global_level_log(double max_c_ef, double max_n_potential, 
    int local_part_count, int local_parts_injected, int local_part_removed)
{
    std::string log = "";
    int64_t global_part_size = 0, global_inj_size = 0, global_removed = 0;
    int64_t glb_parts, gbl_max_parts, gbl_min_parts;
    int64_t glb_part_comms, gbl_max_part_comms, gbl_min_part_comms;
    int64_t glb_sum_max_hops, gbl_max_max_hops, gbl_min_max_hops;
    int global_max_comm_iteration = 0, gbl_move_moreX_hops = 0;

#ifdef USE_MPI
    MPI_Reduce(&OPP_max_comm_iteration, &global_max_comm_iteration, 1, MPI_INT, MPI_MAX, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&OPP_move_moreX_hops, &gbl_move_moreX_hops, 1, MPI_INT, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);

    int64_t temp_local_part_count     = (int64_t)local_part_count;
    int64_t temp_local_parts_injected = (int64_t)local_parts_injected;
    int64_t temp_local_part_removed   = (int64_t)local_part_removed;
    MPI_Reduce(&temp_local_part_count, &global_part_size, 1, MPI_INT64_T, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&temp_local_parts_injected, &global_inj_size, 1, MPI_INT64_T, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(&temp_local_part_removed, &global_removed, 1, MPI_INT64_T, MPI_SUM, OPP_ROOT, MPI_COMM_WORLD);

#else
    global_part_size = local_part_count;
    global_inj_size = local_parts_injected;
    global_removed = local_part_removed;
    global_max_comm_iteration = OPP_max_comm_iteration;
    gbl_max_max_hops = 0;
    gbl_min_max_hops = OPP_move_max_hops;
    gbl_move_moreX_hops = OPP_move_moreX_hops;
#endif

    opp_get_global_values(local_part_count, glb_parts, gbl_max_parts, gbl_min_parts);   
    opp_get_global_values(OPP_part_comm_count_per_iter, glb_part_comms, gbl_max_part_comms, gbl_min_part_comms);
    opp_get_global_values(OPP_move_max_hops, glb_sum_max_hops, gbl_max_max_hops, gbl_min_max_hops);

    log += std::string("\t np: ") + str(global_part_size, "%" PRId64);
    log += std::string(" (") + str(global_inj_size, "%" PRId64);
    log += std::string(" added, ") + str(global_removed, "%" PRId64);
    log += std::string(" removed)\t max c_ef: ") + str(max_c_ef, "%2.15lE");
    log += std::string(" max n_pot: ") + str(max_n_potential, "%2.15lE");
    log += std::string(" | max_comm_iteration: ") + str(global_max_comm_iteration, "%d");

    log += std::string(" | Gbl parts: ") + str(glb_parts, "%" PRId64);
    log += std::string(" Min: ") + str(gbl_min_parts, "%" PRId64);
    log += std::string(" Max: ") + str(gbl_max_parts, "%" PRId64);
    log += std::string(" | Gbl comms: ") + str(glb_part_comms, "%" PRId64);
    log += std::string(" Min: ") + str(gbl_min_part_comms, "%" PRId64);
    log += std::string(" Max: ") + str(gbl_max_part_comms, "%" PRId64);
#ifdef LOG_HOPS
    log += std::string(" | Hops: Min: ") + str(gbl_min_max_hops, "%" PRId64);
    log += std::string(" Max: ") + str(gbl_max_max_hops, "%" PRId64);
    log += std::string(" | more") + std::to_string(X_HOPS) + "_hops: " + str(gbl_move_moreX_hops, "%d");
#endif

    return log;
}


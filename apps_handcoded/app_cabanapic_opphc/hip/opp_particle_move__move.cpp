
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
// AUTO GENERATED CODE
//*********************************************

int m_OPP_HOST_1 = -1;
int m_OPP_HOST_2 = -1;
int m_OPP_HOST_3 = -1;
int m_OPP_HOST_5 = -1;
int m_OPP_HOST_6 = -1;
int m_OPP_HOST_7 = -1;

__constant__ int m_OPP_DEVICE_1;
__constant__ int m_OPP_DEVICE_2;
__constant__ int m_OPP_DEVICE_3;
__constant__ int m_OPP_DEVICE_5;
__constant__ int m_OPP_DEVICE_6;
__constant__ int m_OPP_DEVICE_7;

//user function
//*************************************************************************************************
__device__ void dev_weight_current_to_accumulator__kernel(
        OPP_REAL* cell0_acc,
        const OPP_REAL* q,
        const OPP_REAL& ux, const OPP_REAL& uy, const OPP_REAL& uz,
        const OPP_REAL& dx, const OPP_REAL& dy, const OPP_REAL& dz)
{
    OPP_REAL v0 = 0.0f, v1 = 0.0f, v2 = 0.0f, v3 = 0.0f, v4 = 0.0f, v5 = 0.0f;

    v5 = (*q) * ux * uy * uz * CONST_DEV_one_third;              // Compute correction
 
    #define CALC_J(X,Y,Z)                                                     \
    v4  = (*q)*u##X;             /* v2 = q ux                            */   \
    v1  = v4*d##Y;               /* v1 = q ux dy                         */   \
    v0  = v4-v1;                 /* v0 = q ux (1-dy)                     */   \
    v1 += v4;                    /* v1 = q ux (1+dy)                     */   \
    v4  = CONST_DEV_one+d##Z;    /* v4 = 1+dz                            */   \
    v2  = v0*v4;                 /* v2 = q ux (1-dy)(1+dz)               */   \
    v3  = v1*v4;                 /* v3 = q ux (1+dy)(1+dz)               */   \
    v4  = CONST_DEV_one-d##Z;    /* v4 = 1-dz                            */   \
    v0 *= v4;                    /* v0 = q ux (1-dy)(1-dz)               */   \
    v1 *= v4;                    /* v1 = q ux (1+dy)(1-dz)               */   \
    v0 += v5;                    /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */   \
    v1 -= v5;                    /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */   \
    v2 -= v5;                    /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */   \
    v3 += v5;                    /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

    CALC_J( x,y,z );
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfx + 0)]), v0); 
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfx + 1)]), v1); 
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfx + 2)]), v2); 
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfx + 3)]), v3); 

    CALC_J( y,z,x );
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfy + 0)]), v0); 
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfy + 1)]), v1); 
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfy + 2)]), v2); 
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfy + 3)]), v3); 

    CALC_J( z,x,y );
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfz + 0)]), v0); 
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfz + 1)]), v1); 
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfz + 2)]), v2); 
    atomicAdd(&(cell0_acc[m_OPP_DEVICE_6 * (CellAcc::jfz + 3)]), v3); 

    #undef CALC_J
}

__device__ void dev_push_particles__kernel(opp_move_var& m, 
    OPP_INT* part_cid, 
    OPP_REAL* part_vel, 
    OPP_REAL* part_pos, 
    OPP_REAL* part_streak_mid, 
    const OPP_REAL* part_weight, 
    const OPP_REAL* cell_interp, 
    OPP_REAL* cell_acc,
    const OPP_INT* cell_cell_map)
{
    if (m.iteration_one)
    {
        const OPP_REAL& ex       = cell_interp[m_OPP_DEVICE_5 * CellInterp::ex];
        const OPP_REAL& dexdy    = cell_interp[m_OPP_DEVICE_5 * CellInterp::dexdy];
        const OPP_REAL& dexdz    = cell_interp[m_OPP_DEVICE_5 * CellInterp::dexdz];
        const OPP_REAL& d2exdydz = cell_interp[m_OPP_DEVICE_5 * CellInterp::d2exdydz];
        const OPP_REAL& ey       = cell_interp[m_OPP_DEVICE_5 * CellInterp::ey];
        const OPP_REAL& deydz    = cell_interp[m_OPP_DEVICE_5 * CellInterp::deydz];
        const OPP_REAL& deydx    = cell_interp[m_OPP_DEVICE_5 * CellInterp::deydx];
        const OPP_REAL& d2eydzdx = cell_interp[m_OPP_DEVICE_5 * CellInterp::d2eydzdx];
        const OPP_REAL& ez       = cell_interp[m_OPP_DEVICE_5 * CellInterp::ez];
        const OPP_REAL& dezdx    = cell_interp[m_OPP_DEVICE_5 * CellInterp::dezdx];
        const OPP_REAL& dezdy    = cell_interp[m_OPP_DEVICE_5 * CellInterp::dezdy];
        const OPP_REAL& d2ezdxdy = cell_interp[m_OPP_DEVICE_5 * CellInterp::d2ezdxdy];
        OPP_REAL cbx             = cell_interp[m_OPP_DEVICE_5 * CellInterp::cbx];
        const OPP_REAL& dcbxdx   = cell_interp[m_OPP_DEVICE_5 * CellInterp::dcbxdx];
        OPP_REAL cby             = cell_interp[m_OPP_DEVICE_5 * CellInterp::cby];
        const OPP_REAL& dcbydy   = cell_interp[m_OPP_DEVICE_5 * CellInterp::dcbydy];
        OPP_REAL cbz             = cell_interp[m_OPP_DEVICE_5 * CellInterp::cbz];
        const OPP_REAL& dcbzdz   = cell_interp[m_OPP_DEVICE_5 * CellInterp::dcbzdz];

        const OPP_REAL& dx = part_pos[m_OPP_DEVICE_2 * Dim::x];             // Load position
        const OPP_REAL& dy = part_pos[m_OPP_DEVICE_2 * Dim::y];             // Load position
        const OPP_REAL& dz = part_pos[m_OPP_DEVICE_2 * Dim::z];             // Load position

        const OPP_REAL hax  = CONST_DEV_qdt_2mc * ( ( ex + dy*dexdy ) + dz * ( dexdz + dy*d2exdydz ) );
        const OPP_REAL hay  = CONST_DEV_qdt_2mc * ( ( ey + dz*deydz ) + dx * ( deydx + dz*d2eydzdx ) );
        const OPP_REAL haz  = CONST_DEV_qdt_2mc * ( ( ez + dx*dezdx ) + dy * ( dezdy + dx*d2ezdxdy ) );

        cbx  = cbx + dx*dcbxdx;                     // Interpolate B
        cby  = cby + dy*dcbydy;
        cbz  = cbz + dz*dcbzdz;

        OPP_REAL ux = part_vel[m_OPP_DEVICE_1 * Dim::x];             // Load velocity
        OPP_REAL uy = part_vel[m_OPP_DEVICE_1 * Dim::y];             // Load velocity
        OPP_REAL uz = part_vel[m_OPP_DEVICE_1 * Dim::z];             // Load velocity

        ux  += hax;                                 // Half advance E
        uy  += hay;
        uz  += haz;

        OPP_REAL v0   = CONST_DEV_qdt_2mc/sqrt(1 + (ux*ux + (uy*uy + uz*uz)));
                                                    // Boris - scalars
        OPP_REAL v1   = cbx*cbx + (cby*cby + cbz*cbz);
        OPP_REAL v2   = (v0*v0)*v1;
        OPP_REAL v3   = v0*(1+v2*(CONST_DEV_one_third+v2*CONST_DEV_two_fifteenths));
        OPP_REAL v4   = v3/(1+v1*(v3*v3));
        v4  += v4;
    
        v0   = ux + v3*( uy*cbz - uz*cby );         // Boris - uprime
        v1   = uy + v3*( uz*cbx - ux*cbz );
        v2   = uz + v3*( ux*cby - uy*cbx );
        ux  += v4*( v1*cbz - v2*cby );              // Boris - rotation
        uy  += v4*( v2*cbx - v0*cbz );
        uz  += v4*( v0*cby - v1*cbx );
        ux  += hax;                                 // Half advance E
        uy  += hay;
        uz  += haz;

        part_vel[m_OPP_DEVICE_1 * Dim::x] = ux;                      // save new velocity
        part_vel[m_OPP_DEVICE_1 * Dim::y] = uy;                      // save new velocity
        part_vel[m_OPP_DEVICE_1 * Dim::z] = uz;                      // save new velocity

        /**/                                        // Get norm displacement
        v0   = 1/sqrt(1 + (ux*ux+ (uy*uy + uz*uz)));
        ux  *= CONST_DEV_cdt_d[Dim::x];
        uy  *= CONST_DEV_cdt_d[Dim::y];
        uz  *= CONST_DEV_cdt_d[Dim::z];
    
        ux  *= v0;
        uy  *= v0;
        uz  *= v0;
        
        v0   = dx + ux;                             // Streak midpoint (inbnds)
        v1   = dy + uy;
        v2   = dz + uz;
        v3   = v0 + ux;                             // New position
        v4   = v1 + uy;
        const OPP_REAL v5   = v2 + uz;

        // moving within the cell // Likely
        if (  v3<=1 &&  v4<=1 &&  v5<=1 && -v3<=1 && -v4<=1 && -v5<=1 ) 
        {
            part_pos[m_OPP_DEVICE_2 * Dim::x] = v3;            // save new position
            part_pos[m_OPP_DEVICE_2 * Dim::y] = v4;            // save new position
            part_pos[m_OPP_DEVICE_2 * Dim::z] = v5;            // save new position

            const OPP_REAL q = part_weight[0] * CONST_DEV_qsp;

            dev_weight_current_to_accumulator__kernel(
                cell_acc,
                &q,
                ux, uy, uz,
                v0, v1, v2);

            m.move_status = OPP_MOVE_DONE;
            return;
        }
        else
        {
            part_streak_mid[m_OPP_DEVICE_3 * Dim::x] = ux;
            part_streak_mid[m_OPP_DEVICE_3 * Dim::y] = uy;
            part_streak_mid[m_OPP_DEVICE_3 * Dim::z] = uz;
        }
    }

    OPP_REAL s_dir[3];
    OPP_REAL v0, v1, v2, v3; 
    size_t axis, face;

    OPP_REAL s_midx = part_pos[m_OPP_DEVICE_2 * Dim::x]; // Old positions
    OPP_REAL s_midy = part_pos[m_OPP_DEVICE_2 * Dim::y]; // Old positions
    OPP_REAL s_midz = part_pos[m_OPP_DEVICE_2 * Dim::z]; // Old positions

    OPP_REAL s_dispx = part_streak_mid[m_OPP_DEVICE_3 * Dim::x];   // distance moved during push
    OPP_REAL s_dispy = part_streak_mid[m_OPP_DEVICE_3 * Dim::y];   // distance moved during push
    OPP_REAL s_dispz = part_streak_mid[m_OPP_DEVICE_3 * Dim::z];   // distance moved during push

    s_dir[0] = (s_dispx>0) ? 1 : -1; // Try to get the movement direction
    s_dir[1] = (s_dispy>0) ? 1 : -1;
    s_dir[2] = (s_dispz>0) ? 1 : -1;

    // Compute the twice the fractional distance to each potential
    // streak/cell face intersection.
    v0 = (s_dispx==0) ? 3.4e38 : (s_dir[0]-s_midx)/s_dispx; // 3.4e38 is max of float
    v1 = (s_dispy==0) ? 3.4e38 : (s_dir[1]-s_midy)/s_dispy;
    v2 = (s_dispz==0) ? 3.4e38 : (s_dir[2]-s_midz)/s_dispz;

    // Determine the fractional length and axis of current streak. The
    // streak ends on either the first face intersected by the
    // particle track or at the end of the particle track.
    //
    //   axis 0,1 or 2 ... streak ends on a x,y or z-face respectively
    //   axis 3        ... streak ends at end of the particle track
    /**/      v3=2,  axis=3;
    if(v0<v3) v3=v0, axis=0;
    if(v1<v3) v3=v1, axis=1;
    if(v2<v3) v3=v2, axis=2;
    v3 *= 0.5;

    // Compute the midpoint and the normalized displacement of the streak
    s_dispx *= v3;
    s_dispy *= v3;
    s_dispz *= v3;
    s_midx += s_dispx;
    s_midy += s_dispy;
    s_midz += s_dispz;

    // Accumulate the streak.  Note: accumulator values are 4 times
    // the total physical charge that passed through the appropriate
    // current quadrant in a time-step

    const OPP_REAL q = part_weight[0] * CONST_DEV_qsp;

    dev_weight_current_to_accumulator__kernel(
        cell_acc,
        &q,
        s_dispx, s_dispy, s_dispz,
        s_midx, s_midy, s_midz);

    // Compute the remaining particle displacment
    part_streak_mid[m_OPP_DEVICE_3 * Dim::x] -= s_dispx;
    part_streak_mid[m_OPP_DEVICE_3 * Dim::y] -= s_dispy;
    part_streak_mid[m_OPP_DEVICE_3 * Dim::z] -= s_dispz;

    // Compute the new particle offset
    part_pos[m_OPP_DEVICE_2 * Dim::x] += (s_dispx + s_dispx);
    part_pos[m_OPP_DEVICE_2 * Dim::y] += (s_dispy + s_dispy);
    part_pos[m_OPP_DEVICE_2 * Dim::z] += (s_dispz + s_dispz);

    // If an end streak, return success (should be ~50% of the time)

    if( axis != 3 ) 
    {
        // Determine if the particle crossed into a local cell or if it
        // hit a boundary and convert the coordinate system accordingly.
        // Note: Crossing into a local cell should happen ~50% of the
        // time; hitting a boundary is usually a rare event.  Note: the
        // entry / exit coordinate for the particle is guaranteed to be
        // +/-1 _exactly_ for the particle.

        v0 = s_dir[axis];

        // TODO: this conditional could be better
        if (axis == 0) part_pos[m_OPP_DEVICE_2 * Dim::x] = v0;
        if (axis == 1) part_pos[m_OPP_DEVICE_2 * Dim::y] = v0;
        if (axis == 2) part_pos[m_OPP_DEVICE_2 * Dim::z] = v0;

        // _exactly_ on the boundary.
        face = axis;
        if( v0>0 ) face += 3;

        part_cid[0] = cell_cell_map[m_OPP_DEVICE_7 * face]; 

        // TODO: this conditional/branching could be better
        if (axis == 0) { part_pos[m_OPP_DEVICE_2 * Dim::x] = -v0; /* printf("0\n"); */ }
        if (axis == 1) { part_pos[m_OPP_DEVICE_2 * Dim::y] = -v0; /* printf("1\n"); */ }
        if (axis == 2) { part_pos[m_OPP_DEVICE_2 * Dim::z] = -v0; /* printf("2\n"); */ }

        m.move_status = OPP_NEED_MOVE;
    }
    else
    {
        m.move_status = OPP_MOVE_DONE;
    }
}

//*******************************************************************************
// Returns true only if another hop is required by the current rank
__device__ bool opp_part_check_status_device(opp_move_var& m, int* map0idx, int particle_index, 
        int& remove_count, int *remove_particle_indices, int *move_particle_indices, int *move_cell_indices, int *move_count) 
{
    m.iteration_one = false;

    if (m.move_status == OPP_MOVE_DONE)
    {
        return false;
    }
    else if (m.move_status == OPP_NEED_REMOVE)
    {
        *map0idx = MAX_CELL_INDEX;
        int removeArrayIndex = atomicAdd(&remove_count, 1);
        remove_particle_indices[removeArrayIndex] = particle_index;

        return false;
    }
    else if (*map0idx >= OPP_cells_set_size_d)
    {
        // map0idx cell is not owned by the current mpi rank (it is in the import exec halo region), need to communicate
        int moveArrayIndex = atomicAdd(move_count, 1);
        move_particle_indices[moveArrayIndex] = particle_index;
        move_cell_indices[moveArrayIndex] = *map0idx;

        // Needs to be removed from the current rank, bdw particle packing will be dCONST_DEV_one just prior exchange and removal
        m.move_status = OPP_NEED_REMOVE; 
        int removeArrayIndex = atomicAdd(&remove_count, 1);
        remove_particle_indices[removeArrayIndex] = particle_index;

        return false;
    }

    // map0idx is an own cell and m.move_status == OPP_NEED_MOVE
    return true;
}

// DEVICE kernel function
//*************************************************************************************************
__global__ void opp_device_all_MoveToCells(
    OPP_INT *__restrict d_cell_index,
    OPP_INT *__restrict arg0_dir,           // part_cid           // OPP_RW
    OPP_REAL *__restrict arg1_dir,          // part_vel           // OPP_RW
    OPP_REAL *__restrict arg2_dir,          // part_pos           // OPP_RW
    OPP_REAL *__restrict arg3_dir,          // part_streak_mid    // OPP_RW
    const OPP_REAL *__restrict arg4_dir,    // part_weight        // OPP_READ
    const OPP_REAL *__restrict arg5_ind,    // cell_inter         // OPP_READ
    OPP_REAL *__restrict arg6_ind,          // cell_acc           // OPP_INC
    const OPP_INT *__restrict arg7_ind,     // cell_cell_map      // OPP_READ
    int *__restrict particle_remove_count,
    int *__restrict particle_remove_indices,
    int *__restrict move_particle_indices,
    int *__restrict move_cell_indices,
    int *__restrict move_count,
    int start,
    int end) 
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid + start < end) 
    {
        int n = tid + start;

        opp_move_var m;
        m.iteration_one = (OPP_comm_iteration_d > 0) ? false : true;

        int* map0idx = nullptr; //MAX_CELL_INDEX;

        do
        {
            map0idx = &(arg0_dir[n]);
            // int map0idx = d_cell_index[n]; // TODO : I dont know why this isn't working ??? 
            // arg0_dir and d_cell_index has same pointer values, but this get stuck!

            //user-supplied kernel call
            dev_push_particles__kernel(
                (m),
                (arg0_dir + n),          // part_cid           // OPP_RW
                (arg1_dir + n),          // part_vel           // OPP_RW
                (arg2_dir + n),          // part_pos           // OPP_RW
                (arg3_dir + n),          // part_streak_mid    // OPP_RW
                (arg4_dir + n),          // part_weight        // OPP_READ
                (arg5_ind + *map0idx),   // cell_inter         // OPP_READ
                (arg6_ind + *map0idx),   // cell_acc           // OPP_INC
                (arg7_ind + *map0idx)    // cell_cell_map      // OPP_READ
            );                

        } while (opp_part_check_status_device(m, map0idx, n, 
            *particle_remove_count, particle_remove_indices, move_particle_indices, move_cell_indices, move_count));
    }
}

//*******************************************************************************
void opp_particle_move__move_deposit(
    opp_set set,        // particles_set
    opp_arg arg0,       // part_cid           // OPP_RW
    opp_arg arg1,       // part_vel           // OPP_RW
    opp_arg arg2,       // part_pos           // OPP_RW
    opp_arg arg3,       // part_streak_mid    // OPP_RW
    opp_arg arg4,       // part_weight        // OPP_READ
    opp_arg arg5,       // cell_inter         // OPP_READ
    opp_arg arg6,       // cell_acc           // OPP_INC
    opp_arg arg7        // cell_cell_map      // OPP_READ
)
{

    if (OPP_DBG) 
        opp_printf("CABANA", "opp_particle_move__move_deposit set_size %d diff %d", set->size, set->diff);

    opp_profiler->start("Move");

    int nargs = 8;
    opp_arg args[nargs];

    args[0]  = std::move(arg0);
    args[1]  = std::move(arg1);
    args[2]  = std::move(arg2);
    args[3]  = std::move(arg3);
    args[4]  = std::move(arg4);
    args[5]  = std::move(arg5);
    args[6]  = std::move(arg6);
    args[7]  = std::move(arg7);

    opp_profiler->start("FMv_HaloSend");
    int set_size = opp_mpi_halo_exchanges_grouped(set, nargs, args, Device_GPU);
    opp_profiler->end("FMv_HaloSend");
    opp_profiler->start("FMv_HaloWait");
    opp_mpi_halo_wait_all(nargs, args);
    opp_profiler->end("FMv_HaloWait");

    if (set_size > 0) 
    {
        do {

            m_OPP_HOST_1 = args[1].dat->set->set_capacity;
            m_OPP_HOST_2 = args[2].dat->set->set_capacity;
            m_OPP_HOST_3 = args[3].dat->set->set_capacity; 
            m_OPP_HOST_5 = args[5].dat->set->set_capacity; 
            m_OPP_HOST_6 = args[6].dat->set->set_capacity; 
            m_OPP_HOST_7 = args[7].size; 
            OPP_cells_set_size = set->cells_set->size; 
         
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(m_OPP_DEVICE_1), 
                                                        &m_OPP_HOST_1, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(m_OPP_DEVICE_2), 
                                                        &m_OPP_HOST_2, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(m_OPP_DEVICE_3), 
                                                        &m_OPP_HOST_3, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(m_OPP_DEVICE_5), 
                                                        &m_OPP_HOST_5, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(m_OPP_DEVICE_6), 
                                                        &m_OPP_HOST_6, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(m_OPP_DEVICE_7), 
                                                        &m_OPP_HOST_7, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(OPP_cells_set_size_d), 
                                                        &OPP_cells_set_size, sizeof(int)));
            cutilSafeCall(hipMemcpyToSymbol(HIP_SYMBOL(OPP_comm_iteration_d), 
                                                        &OPP_comm_iteration, sizeof(int)));

            opp_profiler->start("FMv_init_part");
            opp_init_particle_move(set, nargs, args);
            opp_profiler->end("FMv_init_part");

            if (OPP_iter_end - OPP_iter_start > 0) 
            {

                if (OPP_DBG || OPP_comm_iteration > 3)
                    opp_printf("MOVE", "iter %d start %d end %d - COUNT=%d", OPP_comm_iteration, 
                                    OPP_iter_start, OPP_iter_end, (OPP_iter_end - OPP_iter_start));

                int nthread = OPP_gpu_threads_per_block;
                int nblocks = (OPP_iter_end - OPP_iter_start - 1) / nthread + 1;

                cutilSafeCall(hipDeviceSynchronize());
                opp_profiler->start("FMv_OnlyMoveKernel");
                
                opp_device_all_MoveToCells<<<nblocks, nthread>>>(
                    (OPP_INT *)  set->mesh_relation_dat->data_d,
                    (OPP_INT *)  args[0].data_d,                   // part_cid           // OPP_RW
                    (OPP_REAL *) args[1].data_d,                   // part_vel           // OPP_RW
                    (OPP_REAL *) args[2].data_d,                   // part_pos           // OPP_RW
                    (OPP_REAL *) args[3].data_d,                   // part_streak_mid    // OPP_RW
                    (OPP_REAL *) args[4].data_d,                   // part_weight        // OPP_READ
                    (OPP_REAL *) args[5].data_d,                   // cell_inter         // OPP_READ
                    (OPP_REAL *) args[6].data_d,                   // cell_acc           // OPP_INC
                    (OPP_INT *)  args[7].data_d,                   // cell_cell_map      // OPP_READ
                    (int *)      set->particle_remove_count_d,
                    (int *)      OPP_remove_particle_indices_d,
                    (int*)       OPP_move_particle_indices_d,
                    (int*)       OPP_move_cell_indices_d,
                    (int*)       OPP_move_count_d,
                    OPP_iter_start, 
                    OPP_iter_end);

                cutilSafeCall(hipDeviceSynchronize());
                opp_profiler->end("FMv_OnlyMoveKernel");

            }

        } while (opp_finalize_particle_move(set)); 

    }

    opp_set_dirtybit_grouped(nargs, args, Device_GPU);

    opp_profiler->end("Move");
}

//*************************************************************************************************
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


//*************************************************************************************************
inline void move_deposit_kernel_omp(
    char& opp_move_flag, const bool iter_one_flag, const OPP_INT* opp_thr_c2c, OPP_INT* opp_thr_p2c, // added by code-gen
    OPP_REAL* part_vel, 
    OPP_REAL* part_pos, 
    OPP_REAL* part_streak_mid, 
    const OPP_REAL* part_weight, 
    const OPP_REAL* cell_interp, 
    OPP_REAL* cell_acc)
{
    if (iter_one_flag)
    {
        const OPP_REAL& ex       = cell_interp[CellInterp::ex];
        const OPP_REAL& dexdy    = cell_interp[CellInterp::dexdy];
        const OPP_REAL& dexdz    = cell_interp[CellInterp::dexdz];
        const OPP_REAL& d2exdydz = cell_interp[CellInterp::d2exdydz];
        const OPP_REAL& ey       = cell_interp[CellInterp::ey];
        const OPP_REAL& deydz    = cell_interp[CellInterp::deydz];
        const OPP_REAL& deydx    = cell_interp[CellInterp::deydx];
        const OPP_REAL& d2eydzdx = cell_interp[CellInterp::d2eydzdx];
        const OPP_REAL& ez       = cell_interp[CellInterp::ez];
        const OPP_REAL& dezdx    = cell_interp[CellInterp::dezdx];
        const OPP_REAL& dezdy    = cell_interp[CellInterp::dezdy];
        const OPP_REAL& d2ezdxdy = cell_interp[CellInterp::d2ezdxdy];
        OPP_REAL cbx             = cell_interp[CellInterp::cbx];
        const OPP_REAL& dcbxdx   = cell_interp[CellInterp::dcbxdx];
        OPP_REAL cby             = cell_interp[CellInterp::cby];
        const OPP_REAL& dcbydy   = cell_interp[CellInterp::dcbydy];
        OPP_REAL cbz             = cell_interp[CellInterp::cbz];
        const OPP_REAL& dcbzdz   = cell_interp[CellInterp::dcbzdz];

        const OPP_REAL& dx = part_pos[Dim::x];             // Load position
        const OPP_REAL& dy = part_pos[Dim::y];             // Load position
        const OPP_REAL& dz = part_pos[Dim::z];             // Load position

        const OPP_REAL hax  = CONST_qdt_2mc * ( ( ex + dy*dexdy ) + dz * ( dexdz + dy*d2exdydz ) );
        const OPP_REAL hay  = CONST_qdt_2mc * ( ( ey + dz*deydz ) + dx * ( deydx + dz*d2eydzdx ) );
        const OPP_REAL haz  = CONST_qdt_2mc * ( ( ez + dx*dezdx ) + dy * ( dezdy + dx*d2ezdxdy ) );

        cbx  = cbx + dx*dcbxdx;                     // Interpolate B
        cby  = cby + dy*dcbydy;
        cbz  = cbz + dz*dcbzdz;

        OPP_REAL ux = part_vel[Dim::x];             // Load velocity
        OPP_REAL uy = part_vel[Dim::y];             // Load velocity
        OPP_REAL uz = part_vel[Dim::z];             // Load velocity

        ux  += hax;                                 // Half advance E
        uy  += hay;
        uz  += haz;

        OPP_REAL v0   = CONST_qdt_2mc/sqrt(one + (ux*ux + (uy*uy + uz*uz)));
                                                    // Boris - scalars
        OPP_REAL v1   = cbx*cbx + (cby*cby + cbz*cbz);
        OPP_REAL v2   = (v0*v0)*v1;
        OPP_REAL v3   = v0*(one+v2*(one_third+v2*two_fifteenths));
        OPP_REAL v4   = v3/(one+v1*(v3*v3));
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

        part_vel[Dim::x] = ux;                      // save new velocity
        part_vel[Dim::y] = uy;                      // save new velocity
        part_vel[Dim::z] = uz;                      // save new velocity

        /**/                                        // Get norm displacement
        v0   = one/sqrt(one + (ux*ux+ (uy*uy + uz*uz)));
        ux  *= CONST_cdt_d[Dim::x];
        uy  *= CONST_cdt_d[Dim::y];
        uz  *= CONST_cdt_d[Dim::z];
    
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
        if (  v3<=one &&  v4<=one &&  v5<=one && -v3<=one && -v4<=one && -v5<=one ) 
        {
            part_pos[Dim::x] = v3;            // save new position
            part_pos[Dim::y] = v4;            // save new position
            part_pos[Dim::z] = v5;            // save new position

            const OPP_REAL q = part_weight[0] * CONST_qsp;

            weight_current_to_accumulator_kernel(
                cell_acc,
                &q,
                ux, uy, uz,
                v0, v1, v2);

            opp_move_flag = OPP_MOVE_DONE;
            return;
        }
        else
        {
            part_streak_mid[Dim::x] = ux;
            part_streak_mid[Dim::y] = uy;
            part_streak_mid[Dim::z] = uz;
        }
    }

    OPP_REAL s_dir[3];
    OPP_REAL v0, v1, v2, v3; 
    size_t axis, face;

    OPP_REAL s_midx = part_pos[Dim::x]; // Old positions
    OPP_REAL s_midy = part_pos[Dim::y]; // Old positions
    OPP_REAL s_midz = part_pos[Dim::z]; // Old positions

    OPP_REAL s_dispx = part_streak_mid[Dim::x];   // distance moved during push
    OPP_REAL s_dispy = part_streak_mid[Dim::y];   // distance moved during push
    OPP_REAL s_dispz = part_streak_mid[Dim::z];   // distance moved during push

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

    const OPP_REAL q = part_weight[0] * CONST_qsp;

    weight_current_to_accumulator_kernel(
        cell_acc,
        &q,
        s_dispx, s_dispy, s_dispz,
        s_midx, s_midy, s_midz);

    // Compute the remaining particle displacment
    part_streak_mid[Dim::x] -= s_dispx;
    part_streak_mid[Dim::y] -= s_dispy;
    part_streak_mid[Dim::z] -= s_dispz;

    // Compute the new particle offset
    part_pos[Dim::x] += (s_dispx + s_dispx);
    part_pos[Dim::y] += (s_dispy + s_dispy);
    part_pos[Dim::z] += (s_dispz + s_dispz);

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
        if (axis == 0) part_pos[Dim::x] = v0;
        if (axis == 1) part_pos[Dim::y] = v0;
        if (axis == 2) part_pos[Dim::z] = v0;

        // _exactly_ on the boundary.
        face = axis;
        if( v0>0 ) face += 3;
        
        (*opp_thr_p2c) =  opp_thr_c2c[face];

        // TODO: this conditional/branching could be better
        if (axis == 0) { part_pos[Dim::x] = -v0; /* printf("0\n"); */ }
        if (axis == 1) { part_pos[Dim::y] = -v0; /* printf("1\n"); */ }
        if (axis == 2) { part_pos[Dim::z] = -v0; /* printf("2\n"); */ }

        opp_move_flag = OPP_NEED_MOVE;
    }
    else
    {
        opp_move_flag = OPP_MOVE_DONE;
    }
}
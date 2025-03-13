#pragma once

#include "cabana_defs.h"
#include "math.h"

const OPP_REAL fourth         = (1.0 / 4.0);
const OPP_REAL half           = (1.0 / 2.0);
const OPP_REAL one            = 1.0;
const OPP_REAL one_third      = (1.0 / 3.0);
const OPP_REAL two_fifteenths = (2.0 / 15.0);

//*************************************************************************************************
inline void interpolate_mesh_fields_kernel(
    const OPP_REAL* cell0_e,   
    const OPP_REAL* cell0_b,   
    const OPP_REAL* cell_x_e,  
    const OPP_REAL* cell_y_e,  
    const OPP_REAL* cell_z_e,  
    const OPP_REAL* cell_xy_e, 
    const OPP_REAL* cell_yz_e, 
    const OPP_REAL* cell_xz_e, 
    const OPP_REAL* cell_x_b,  
    const OPP_REAL* cell_y_b,  
    const OPP_REAL* cell_z_b,  
    OPP_REAL* cell0_interp,
    const OPP_INT* cell0_ghost)    
{
    if (cell0_ghost[0] == 0)
    {
        OPP_REAL w0 = 0.0, w1 = 0.0, w2 = 0.0, w3 = 0.0;

        // ex interpolation coefficients
        w0 = cell0_e[Dim::x];                       // pf0->ex;
        w1 = cell_y_e[Dim::x];                      // pfy->ex;
        w2 = cell_z_e[Dim::x];                      // pfz->ex;
        w3 = cell_yz_e[Dim::x];                     // pfyz->ex;

        cell0_interp[CellInterp::ex]       = fourth*( (w3 + w0) + (w1 + w2) );
        cell0_interp[CellInterp::dexdy]    = fourth*( (w3 - w0) + (w1 - w2) );
        cell0_interp[CellInterp::dexdz]    = fourth*( (w3 - w0) - (w1 - w2) );
        cell0_interp[CellInterp::d2exdydz] = fourth*( (w3 + w0) - (w1 + w2) );

        // ey interpolation coefficients
        w0 = cell0_e[Dim::y];
        w1 = cell_z_e[Dim::y];                       // pfz->ey;
        w2 = cell_x_e[Dim::y];                       // pfx->ey;
        w3 = cell_xz_e[Dim::y];                      // pfzx->ey;

        cell0_interp[CellInterp::ey]       = fourth*( (w3 + w0) + (w1 + w2) );
        cell0_interp[CellInterp::deydz]    = fourth*( (w3 - w0) + (w1 - w2) );
        cell0_interp[CellInterp::deydx]    = fourth*( (w3 - w0) - (w1 - w2) );
        cell0_interp[CellInterp::d2eydzdx] = fourth*( (w3 + w0) - (w1 + w2) );

        // ez interpolation coefficients
        w0 = cell0_e[Dim::z];                       // pf0->ez;
        w1 = cell_x_e[Dim::z];                      // pfx->ez;
        w2 = cell_y_e[Dim::z];                      // pfy->ez;
        w3 = cell_xy_e[Dim::z];                     // pfxy->ez;

        cell0_interp[CellInterp::ez]       = fourth*( (w3 + w0) + (w1 + w2) );
        cell0_interp[CellInterp::dezdx]    = fourth*( (w3 - w0) + (w1 - w2) );
        cell0_interp[CellInterp::dezdy]    = fourth*( (w3 - w0) - (w1 - w2) );
        cell0_interp[CellInterp::d2ezdxdy] = fourth*( (w3 + w0) - (w1 + w2) );

        // bx interpolation coefficients
        w0 = cell0_b[Dim::x];                      // pf0->cbx;
        w1 = cell_x_b[Dim::x];                     // pfx->cbx;
        cell0_interp[CellInterp::cbx]    = half*( w1 + w0 );
        cell0_interp[CellInterp::dcbxdx] = half*( w1 - w0 );

        // by interpolation coefficients
        w0 = cell0_b[Dim::y];                      // pf0->cby;
        w1 = cell_y_b[Dim::y];                     // pfy->cby;
        cell0_interp[CellInterp::cby]    = half*( w1 + w0 );
        cell0_interp[CellInterp::dcbydy] = half*( w1 - w0 );

        // bz interpolation coefficients
        w0 = cell0_b[Dim::z];                      // pf0->cbz;
        w1 = cell_z_b[Dim::z];                     // pfz->cbz;
        cell0_interp[CellInterp::cbz]    = half*( w1 + w0 );
        cell0_interp[CellInterp::dcbzdz] = half*( w1 - w0 );
    }
}

//*************************************************************************************************
inline void clear_accumulator_fields_kernel(
    OPP_REAL* acc)
{
    for (int i = 0; i < ACCUMULATOR_ARRAY_LENGTH * DIM; i++)
        acc[i] = 0.0;
}

//*************************************************************************************************
inline void weight_current_to_accumulator_kernel(
        OPP_REAL* cell0_acc,
        const OPP_REAL* q,
        const OPP_REAL& ux, const OPP_REAL& uy, const OPP_REAL& uz,
        const OPP_REAL& dx, const OPP_REAL& dy, const OPP_REAL& dz)
{
    OPP_REAL v0 = 0.0f, v1 = 0.0f, v2 = 0.0f, v3 = 0.0f, v4 = 0.0f, v5 = 0.0f;

    v5 = (*q) * ux * uy * uz * one_third;              // Compute correction
 
    #define CALC_J(X,Y,Z)                                           \
    v4  = (*q)*u##X;   /* v2 = q ux                            */   \
    v1  = v4*d##Y;     /* v1 = q ux dy                         */   \
    v0  = v4-v1;       /* v0 = q ux (1-dy)                     */   \
    v1 += v4;          /* v1 = q ux (1+dy)                     */   \
    v4  = one+d##Z;    /* v4 = 1+dz                            */   \
    v2  = v0*v4;       /* v2 = q ux (1-dy)(1+dz)               */   \
    v3  = v1*v4;       /* v3 = q ux (1+dy)(1+dz)               */   \
    v4  = one-d##Z;    /* v4 = 1-dz                            */   \
    v0 *= v4;          /* v0 = q ux (1-dy)(1-dz)               */   \
    v1 *= v4;          /* v1 = q ux (1+dy)(1-dz)               */   \
    v0 += v5;          /* v0 = q ux [ (1-dy)(1-dz) + uy*uz/3 ] */   \
    v1 -= v5;          /* v1 = q ux [ (1+dy)(1-dz) - uy*uz/3 ] */   \
    v2 -= v5;          /* v2 = q ux [ (1-dy)(1+dz) - uy*uz/3 ] */   \
    v3 += v5;          /* v3 = q ux [ (1+dy)(1+dz) + uy*uz/3 ] */

    CALC_J( x,y,z );
    cell0_acc[CellAcc::jfx + 0] += v0; 
    cell0_acc[CellAcc::jfx + 1] += v1; 
    cell0_acc[CellAcc::jfx + 2] += v2; 
    cell0_acc[CellAcc::jfx + 3] += v3; 

    CALC_J( y,z,x );
    cell0_acc[CellAcc::jfy + 0] += v0; 
    cell0_acc[CellAcc::jfy + 1] += v1; 
    cell0_acc[CellAcc::jfy + 2] += v2; 
    cell0_acc[CellAcc::jfy + 3] += v3; 

    CALC_J( z,x,y );
    cell0_acc[CellAcc::jfz + 0] += v0; 
    cell0_acc[CellAcc::jfz + 1] += v1; 
    cell0_acc[CellAcc::jfz + 2] += v2; 
    cell0_acc[CellAcc::jfz + 3] += v3; 

    #undef CALC_J
}

//*************************************************************************************************
inline void move_deposit_kernel(
    OPP_REAL* part_vel, 
    OPP_REAL* part_pos, 
    OPP_REAL* part_streak_mid, 
    const OPP_REAL* part_weight, 
    const OPP_REAL* cell_interp, 
    OPP_REAL* cell_acc)
{
    if (OPP_DO_ONCE)
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

            OPP_PARTICLE_MOVE_DONE;
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
        
        (*opp_p2c) =  opp_c2c[face];

        // TODO: this conditional/branching could be better
        if (axis == 0) { part_pos[Dim::x] = -v0; /* printf("0\n"); */ }
        if (axis == 1) { part_pos[Dim::y] = -v0; /* printf("1\n"); */ }
        if (axis == 2) { part_pos[Dim::z] = -v0; /* printf("2\n"); */ }

        OPP_PARTICLE_NEED_MOVE;
    }
    else
    {
        OPP_PARTICLE_MOVE_DONE;
    }
}

//*************************************************************************************************
inline void accumulate_current_to_cells_kernel(
        OPP_REAL* cell0_j, 
        const OPP_REAL* cell0_acc,
        const OPP_REAL* cell_xd_acc, 
        const OPP_REAL* cell_yd_acc, 
        const OPP_REAL* cell_zd_acc, 
        const OPP_REAL* cell_xyd_acc, 
        const OPP_REAL* cell_yzd_acc, 
        const OPP_REAL* cell_xzd_acc,
        const OPP_INT* iter_acc)
{
    if (iter_acc[0] == 1)
    {
        cell0_j[Dim::x] = CONST_acc_coef[Dim::x] * (cell0_acc[CellAcc::jfx + 0] +
                                                    cell_yd_acc[CellAcc::jfx + 1] +
                                                    cell_zd_acc[CellAcc::jfx + 2] +
                                                    cell_yzd_acc[CellAcc::jfx + 3]);

        cell0_j[Dim::y] = CONST_acc_coef[Dim::y] * (cell0_acc[CellAcc::jfy + 0] +
                                                    cell_zd_acc[CellAcc::jfy + 1] +
                                                    cell_xd_acc[CellAcc::jfy + 2] +
                                                    cell_xzd_acc[CellAcc::jfy + 3]);

        cell0_j[Dim::z] = CONST_acc_coef[Dim::z] * (cell0_acc[CellAcc::jfz + 0] +
                                                    cell_xd_acc[CellAcc::jfz + 1] +
                                                    cell_yd_acc[CellAcc::jfz + 2] +
                                                    cell_xyd_acc[CellAcc::jfz + 3]);
    }
}

//*************************************************************************************************
inline void half_advance_b_kernel (
    const OPP_REAL* cell_x_e, 
    const OPP_REAL* cell_y_e, 
    const OPP_REAL* cell_z_e, 
    const OPP_REAL* cell0_e, 
    OPP_REAL* cell0_b,
    const OPP_INT* cell0_ghost)
{
    if (cell0_ghost[0] == 0)
    {
        cell0_b[Dim::x] -= (0.5 * CONST_p[Dim::y] * (cell_y_e[Dim::z] - cell0_e[Dim::z]) 
                            - 0.5 * CONST_p[Dim::z] * (cell_z_e[Dim::y] - cell0_e[Dim::y]));

        cell0_b[Dim::y] -= (0.5 * CONST_p[Dim::z] * (cell_z_e[Dim::x] - cell0_e[Dim::x]) 
                            - 0.5 * CONST_p[Dim::x] * (cell_x_e[Dim::z] - cell0_e[Dim::z]));

        cell0_b[Dim::z] -= (0.5 * CONST_p[Dim::x] * (cell_x_e[Dim::y] - cell0_e[Dim::y]) 
                            - 0.5 * CONST_p[Dim::y] * (cell_y_e[Dim::x] - cell0_e[Dim::x]));
    }
}

//*************************************************************************************************
inline void advance_e_kernel (
    const OPP_REAL* cell_x_b, 
    const OPP_REAL* cell_y_b, 
    const OPP_REAL* cell_z_b, 
    const OPP_REAL* cell0_b, 
    const OPP_REAL* cell0_j, 
    OPP_REAL* cell0_e,
    const OPP_INT* iter_adv_e)
{
    if (iter_adv_e[0] == 1)
    {
        cell0_e[Dim::x] += ( - CONST_dt_eps0 * cell0_j[Dim::x] ) + 
            ( CONST_p[Dim::y] * (cell0_b[Dim::z] - cell_y_b[Dim::z]) - 
            CONST_p[Dim::z] * (cell0_b[Dim::y] - cell_z_b[Dim::y]) );

        cell0_e[Dim::y] += ( - CONST_dt_eps0 * cell0_j[Dim::y] ) +            
            ( CONST_p[Dim::z] * (cell0_b[Dim::x] - cell_z_b[Dim::x]) - 
            CONST_p[Dim::x] * (cell0_b[Dim::z] - cell_x_b[Dim::z]) );

        cell0_e[Dim::z] += ( - CONST_dt_eps0 * cell0_j[Dim::z] ) +           
            ( CONST_p[Dim::x] * (cell0_b[Dim::y] - cell_x_b[Dim::y]) - 
            CONST_p[Dim::y] * (cell0_b[Dim::x] - cell_y_b[Dim::x]) );
    }
}

//*******************************************************************************
inline void get_max_x_values_kernel(
    const OPP_REAL* cell_j,
    OPP_REAL* max_j,
    const OPP_REAL* cell_e,
    OPP_REAL* max_e,
    const OPP_REAL* cell_b,
    OPP_REAL* max_b)
{
    *max_j = MAX(*cell_j, *max_j);
    
    *max_e = MAX(*cell_e, *max_e);

    *max_b = MAX(*cell_b, *max_b);
}

//*******************************************************************************
inline void update_ghosts_B_kernel(
    const OPP_INT* c_mask_ugb,
    const OPP_REAL* from_cell,
    OPP_REAL* to_cell,
    const OPP_INT* m_idx)
{
    if (c_mask_ugb[*m_idx] == 1)
    {
        to_cell[Dim::x] = from_cell[Dim::x];
        to_cell[Dim::y] = from_cell[Dim::y];
        to_cell[Dim::z] = from_cell[Dim::z];
    }  
}
//*******************************************************************************
inline void update_ghosts_kernel(
    const OPP_INT* c_mask_ug,
    const OPP_REAL* from_cell,
    OPP_REAL* to_cell,
    const OPP_INT* m_idx, 
    const OPP_INT* dim)
{
    if (c_mask_ug[*m_idx] == 1)
        to_cell[*dim] += from_cell[*dim];
}


//*******************************************************************************
inline void compute_energy_kernel(
    const OPP_INT* cell0_ghost, 
    const OPP_REAL* cell_field, 
    OPP_REAL* energy)
{
    if (cell0_ghost[0] == 0) 
    {
        energy[0] += cell_field[Dim::x] * cell_field[Dim::x] +
                    cell_field[Dim::y] * cell_field[Dim::y] +
                    cell_field[Dim::z] * cell_field[Dim::z];
    }
};

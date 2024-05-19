#pragma once

//*********************************************
// USER WRITTEN CODE
//*********************************************


// #include "simpic.h"

//*************************************************************************************************
/*Interpolates values for particle between two grid points.
* Used to compute E field in particle position*/
void weight_f2p_kernel(
        const double* node0_field_E,  //LHS
        const double* node1_field_E,  //RHS
        const double* particle0_position_x,
        double* particle0_field_E
    )
{
    double xx = ((particle0_position_x[0] - CONST_xl[0]) / CONST_dx[0]); // Makes Global position to local position comapared to the cell
    int n = int(xx);
    double frac = (xx - n);

    particle0_field_E[0] = ((frac * node1_field_E[0]) + ((1.0 - frac) * node0_field_E[0]));
};


//*************************************************************************************************
void move_kernel(
        const double* part_field_E,
        double* part_velocity_x,
        double* part_position_x
    )
{
    if (OPP_DO_ONCE)
    {
        part_velocity_x[0] += (CONST_qm[0] * part_field_E[0]);
        part_position_x[0] += (part_velocity_x[0] * CONST_dt[0]);
    }

    // since particle cell index can be easily calculated with global positioning, no need to search by iterating
    if ((part_position_x[0] > CONST_xl[0]) && (part_position_x[0] < CONST_xr[0]))
    {    
        double xx = ((part_position_x[0] - CONST_xl[0]) / CONST_dx[0]); // Makes Global position to local position comapared to the cell   
        opp_p2c[0] = int(xx);
        OPP_PARTICLE_MOVE_DONE;
    }
    else if ((part_position_x[0] >= CONST_xl[0]) && (CONST_rank[0] == 0) ||
             (part_position_x[0] <= CONST_xr[0]) && (CONST_rank[0] == (CONST_comm_size[0]-1)))
    {
        opp_p2c[0] = MAX_CELL_INDEX;
        OPP_PARTICLE_NEED_REMOVE;
    }
    else
    {
        opp_p2c[0] = (part_position_x[0] < CONST_xl[0]) ? 
                        CONST_neighbour_cell[Dir::Left] : CONST_neighbour_cell[Dir::Right];
        OPP_PARTICLE_NEED_MOVE;
    }
};


//*************************************************************************************************
/*Extrapolates values from particles to grid points.
 * Used to calculate density at the grid points*/
void weight_p2f_kernel(
        double* node0_field_J, 
        double* node1_field_J,
        const double* particle0_position_x
    )
{
    double xx = ((particle0_position_x[0] - CONST_xl[0]) / CONST_dx[0]); // Makes Global position to local position comapared to the cell
    int n = int(xx);
    double frac = (xx - n);

    (*node0_field_J) += (CONST_qscale[0] * (1.0 - frac));  // Can change qscale to be calculated from particle data
    (*node1_field_J) += (CONST_qscale[0] * frac);
}


//*************************************************************************************************
void sum_laplace_kernel(
        const double* node0_xlocal,
        double* node0_field_P
    )
{
    double rv = 0.0;
    double lv = CONST_lhs_voltage[0];

    double frac = ((*node0_xlocal) / CONST_L[0]);
    (*node0_field_P) += (frac * rv + (1. - frac) * lv);
}


//*************************************************************************************************

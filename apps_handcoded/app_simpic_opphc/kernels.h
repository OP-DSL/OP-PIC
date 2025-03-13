#pragma once

//*********************************************
// USER WRITTEN CODE
//*********************************************


#include "simpic.h"

//*************************************************************************************************
/*Interpolates values for particle between two grid points.
* Used to compute E field in particle position*/
void weight_fields_to_particles__kernel(
        const double* node0_field_E,  //LHS
        const double* node1_field_E,  //RHS
        const double* particle0_position_x,
        double* particle0_field_E
    )
{
    double xx = (((*particle0_position_x) - xl)/dx); // Makes Global position to local position comapared to the cell
    int n = int(xx);
    double frac = (xx - n);

    (*particle0_field_E) = ((frac * (*node1_field_E)) + ((1.0 - frac) * (*node0_field_E)));
};


//*************************************************************************************************
void particle_push__kernel(
        const double* part_field_E,
        double* part_velocity_x,
        double* part_position_x,
        int* part_cell_index
    )
{
    (*part_velocity_x) += (qm * (*part_field_E));
    (*part_position_x) += ((*part_velocity_x) * dt);

    double xx = (((*part_position_x) - xl)/dx); // Makes Global position to local position comapared to the cell
    (*part_cell_index) = int(xx);

    // since particle cell index can be easily calculated with global positioning, no need to search by iterating
    OPP_PARTICLE_MOVE_DONE;
};


//*************************************************************************************************
/*Extrapolates values from particles to grid points.
 * Used to calculate density at the grid points*/
void weight_particles_to_fields__kernel(
        double* node0_field_J, 
        double* node1_field_J,
        const double* particle0_position_x
    )
{
    double xx = (((*particle0_position_x) - xl)/dx); // Makes Global position to local position comapared to the cell
    int n = int(xx);
    double frac = (xx - n);

    (*node0_field_J) += (qscale * (1.0 - frac));  // Can change qscale to be calculated from particle data
    (*node1_field_J) += (qscale * frac);
}


//*************************************************************************************************
void field_solve_sum_laplace__kernel(
        const double* node0_xlocal,
        double* node0_field_P
    )
{
    double max_time = t;

    double rv = rhsV(max_time);
    double lv = lhsV(max_time);

    double frac = ((*node0_xlocal) / L);
    (*node0_field_P) += (frac * rv + (1. - frac) * lv);
}


//*************************************************************************************************

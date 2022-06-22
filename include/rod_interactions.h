//
//  This file is part of the FFEA simulation package
//
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file.
//
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
//
//  To help us fund FFEA development, we humbly ask that you cite
//  the research papers on the package.
//

/*
 *      rod_interactions.h
 *	Author: Ryan Cocking, University of Leeds
 *	Email: bsrctb@leeds.ac.uk
 */

#ifndef ROD_INTERACTIONS
#define ROD_INTERACTIONS

#include "rod_math_v9.h"

namespace rod
{

struct InteractionData
{
    int rod_id_self;
    int rod_id_nbr;
    int element_id_self;
    int element_id_nbr;
    float radius_self;
    float radius_nbr;
    float contact_self[3];
    float contact_nbr[3];

    InteractionData(
        int rod_id_a,
        int rod_id_b,
        int elem_id_a,
        int elem_id_b,
        float rad_a,
        float rad_b,
        float c_a[3],
        float c_b[3]);
};

void rod_distance_correction(
    float c_a[3],
    float c_b[3],
    float r_a[3],
    float r_b[3],
    float p_a[3],
    float p_b[3],
    OUT float c_a_out[3],
    float c_b_out[3]);

void get_shortest_distance_to_rod(
    float p_a[3],
    float p_b[3],
    float r_a[3],
    float r_b[3],
    OUT float c_a[3],
    float c_b[3]);

void set_element_neighbours(
    int rod_id_a,
    int rod_id_b,
    int elem_id_a,
    int elem_id_b,
    float p_a[3],
    float p_b[3],
    float r_a[3],
    float r_b[3],
    float radius_a,
    float radius_b,
    std::vector<InteractionData> &neighbours_a,
    std::vector<InteractionData> &neighbours_b);

float steric_energy_linear(
    float force_scaling_factor,
    float intersect_distance,
    float radius_sum);

float perturbed_intersection_distance(
    int perturb_dim,
    float perturb_delta,
    float contact_a[3],
    float contact_b[3],
    float radius_sum);

std::vector<float> element_steric_force(
    float delta,
    float force_strength,
    float radius_sum,
    float contact_self[3],
    float contact_neighb[3]);

std::vector<float> node_force_interpolation(
    float contact[3],
    float node_1[3],
    float element_length,
    const std::vector<float> &element_force);
}
#endif

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
        int elem_id_self;
        int elem_id_nbr;
        float radius_self;
        float radius_nbr;
        std::array<float, 3> contact_self;
        std::array<float, 3> contact_nbr;
        std::array<float, 3> c_ab;
        std::array<float, 3> img_shift;
        std::array<float, 3> r_self;
        std::array<float, 3> r_nbr;
        float epsilon;
        float sigma;
        float r_min;
        float r_min_inv;

        InteractionData(
            int rod_id_a,
            int rod_id_b,
            int elem_id_a,
            int elem_id_b,
            float radius_a,
            float radius_b,
            const std::array<float, 3> &c_a,
            const std::array<float, 3> &c_b,
            const std::array<float, 3> &shift,
            const std::array<float, 3> &r_a,
            const std::array<float, 3> &r_b);

        InteractionData(
            int rod_id_a,
            int rod_id_b,
            int elem_id_a,
            int elem_id_b,
            float radius_a,
            float radius_b,
            const float3 &c_a,
            const float3&c_b,
            const float3&shift,
            const float3&r_a,
            const float3&r_b,
            float eps,
            float sig);

        bool elements_intersect() const;
    };


    // previously: snap_to_nodes
    void finite_length_correction(
        float3 &c,
        const float3 &r,
        const float3 &p);

    void finite_length_correction(
        const float3 &c,
        const float3 &r,
        const float3 &p,
        OUT float3 &c_out);

    void nearest_node_correction(
        float3 &c_a,
        const float3 &c_b,
        const float3 &r_a,
        const float3 &r_b,
        const float3 &p_a,
        const float3 &p_b);

    void element_minimum_displacement(
        const float3 &p_a,
        const float3 &p_b,
        const float3 &r_a,
        const float3 &r_b,
        float3 &c_a,
        float3 &c_b);

    // Not implemented?
    //std::vector<int> nearest_periodic_image(
    //    float3 &p_a,
    //    float3 &p_b,
    //    float3 &r_a,
    //    float3 &r_b,
    //    const std::vector<float> &box_dim);
    // Not implemented?
    //std::vector<int> nearest_periodic_image(
    //    float3 &displacement,
    //    const std::vector<float> &box_dim);

    std::vector<int> nearest_periodic_image(
        const float3 &a,
        const float3 &b,
        const std::vector<float> &box_dim);

    void set_steric_nbrs(
        int rod_id_a,
        int rod_id_b,
        int elem_id_a,
        int elem_id_b,
        const float3 &p_a,
        const float3 &p_b,
        float3 &r_a,
        float3 &r_b,
        float radius_a,
        float radius_b,
        std::vector<InteractionData> &neighbours_a,
        std::vector<InteractionData> &neighbours_b,
        bool periodic,
        const std::vector<float> &box_dim);

    float steric_energy_linear(
        float force_scaling_factor,
        float intersect_distance);

    float steric_energy_squared(
        float force_scaling_factor,
        float intersect_distance);

    // Perturbation in +x/-x/+y/-y/+z/-z
    float intersection_distance(
        int dim,
        float delta,
        const float3 &c_a,
        const float3 &c_b,
        float radius_sum);

    // Perturbation in +c_ab/-c_ab
    float intersection_distance(
        float delta,
        const float3 &c_a,
        const float3 &c_b,
        float radius_sum);

    // No perturbation
    float intersection_distance(
        const float3 &c_a,
        const float3 &c_b,
        float radius_sum);

    std::vector<float> element_steric_energy(
        float delta,
        float force_strength,
        float radius_sum,
        const float3 &c_a,
        const float3 &c_b);

    std::vector<float> node_force_interpolation(
        const float3 & contact,
        const float3 & node_1,
        float element_length,
        const std::vector<float> &element_force);

    // ! - should ideally be in their own file
    /*
    ================================================================================
        VAN DER WAALS INTERACTIONS
    ================================================================================
    */

    float vdw_energy_6_12(float r_inv, float eps, float sig);
    float vdw_force_6_12(float r_inv, float eps, float sig);
    float vdw_energy_interp(float r, float eps, float r_min_inv);
    float vdw_force_interp(float r, float eps, float r_min_inv);

    struct VDWSite
    {
        int rod_id;     // Parent rod index
        int elem_id;
        int site_id;    // Site index within parent rod
        int vdw_type;
        float L_rod;    // Fraction of rod length at which site is located, 0 <= L <= 1
        float L_elem;
        float3 pos;

        VDWSite(const int rodid, const int siteid, const int vdwtype, const float lrod,  const std::vector<float> &r_rod,
            const float contour_length, const int num_nodes);
        void print_info() const;
        void get_parent_element(
            const std::vector<float> &p_rod,
            const float contour_length,
            const float norm_length_along_rod,
            const int num_nodes);
        void update_position(const std::vector<float> &r_rod);

    };

    void set_vdw_nbrs(
        VDWSite site_a,
        VDWSite site_b,
        const float3 &p_a,
        const float3 &p_b,
        const float3 &r_a,
        const float3 &r_b,
        float radius_a,
        float radius_b,
        std::vector<InteractionData> &nbr_a,
        std::vector<InteractionData> &nbr_b,
        bool periodic,
        std::vector<float> box_dim,
        float vdw_cutoff,
        float epsilon,
        float sigma);

}
#endif

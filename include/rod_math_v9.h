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
 *      rod_math_v9.h
 *	Author: Rob Welch, University of Leeds
 *	Email: py12rw@leeds.ac.uk

 *	Author: Ryan Cocking, University of Leeds
 *	Email: bsrctb@leeds.ac.uk
 */
#ifndef ROD_MATH
#define ROD_MATH

#define OUT               ///< This is used to denote when a function modifies one of its parameters
#define _USE_MATH_DEFINES ///<  This has to come before including cmath

#include <cmath>
#include <iostream>
#include <assert.h>
#include "dimensions.h"
#include <stdlib.h>
#include <boost/math/special_functions/fpclassify.hpp>
#include <vector>
#include <string>
//#include <fenv.h>

namespace rod
{
    typedef std::array<float, 2> float2;
    typedef std::array<float, 3> float3;
    typedef std::array<float, 4> float4;
    typedef std::array<float, 6> float6;
    typedef std::array<float, 9> float9;
    typedef std::array<float3, 2> float2x3;
    typedef std::array<float4, 2> float2x4;
    typedef std::array<float3, 3> float3x3;
    typedef std::array<float3, 4> float4x3;
    typedef std::array<float4, 4> float4x4;
    typedef std::array<int, 3> int3;
    extern bool dbg_print;

    const static bool debug_nan = false;

    const double boltzmann_constant = 1.3806503e-23 / mesoDimensions::Energy;

    // A less weird way to access the contents of our arrays representing our vectors
    static const int x = 0;
    static const int y = 1;
    static const int z = 2;

    // For clarity, anytime you see an index [x], we're referring to the x
    // dimension.

    // Computing the energy for a perturbation requires four segments, defined
    // by 5 nodes. They are i-2 to i+1, listed here.
    static const int im2 = 0; ///< index of i-2nd thing
    static const int im1 = 1; ///< index of i-1st thing
    static const int i = 2;   ///< index of ith thing
    static const int ip1 = 3; ///< index of i+1th thing

#define OMP_SIMD_INTERNAL _Pragma("omp simd")

#define vec3d(x) for (int x = 0; x < 3; ++x) ///< Shorthand to loop over elements of our 1d arrays representing 3d vectors

    static const float rod_software_version = 1.6;

    void rod_abort(std::string message);

    // These are just generic vector functions that will be replaced by mat_vec_fns at some point

    template<typename T, size_t N>
    void print_array(std::string array_name, const std::array<T, N>& arr);
    template<typename T>
    void print_array(std::string array_name, const std::vector<T>& vec);
    void print_array(std::string array_name, const std::vector<float> &vec, int start, int end);
    void write_vector(FILE *file_ptr, const std::vector<int> &vec, bool new_line);
    void write_vector(FILE *file_ptr, const std::vector<float> &vec, float unit_scale_factor, bool new_line);
    void print_vector(std::string vector_name, const std::vector<float> &vec);
    void print_vector(std::string vector_name, const std::vector<int> &vec);
    void print_vector(std::string vector_name, std::vector<float>::iterator start, std::vector<float>::iterator end);
    void print_vector(std::string vector_name, const std::vector<float> &vec, int start_ind, int end_ind);
    void print_vector(std::string vector_name, const std::vector<int> &vec, int start_ind, int end_ind);
    std::vector<float> slice_vector(std::vector<float> vec, int start_index, int end_index);
    void normalize(const float3 &in, OUT float3 &out);
    void normalize(std::vector<float> in, OUT std::vector<float> out);
    void normalize_unsafe(const float3 &in, OUT float3 &out);
    float absolute(const float3 &in);
    float absolute(const std::vector<float> &in);
    void cross_product(const float3 &a, const float3 &b, float3 &out);
    void cross_product_unsafe(const float3 &a, const float3 &b, float3 &out);
    void get_rotation_matrix(const float3 &a, const float3 &b, float9 &rotation_matrix);
    void get_cartesian_rotation_matrix(int dim, float angle, float9 &rotation_matrix);
    void apply_rotation_matrix(arr3_view<float, float3> vec, const float9 &matrix, OUT float3 &rotated_vec);
    void apply_rotation_matrix_row(arr3_view<float, float3> vec, const float9 &matrix, OUT float3 &rotated_vec);
    void matmul_3x3_3x3(const float9 &a, const float9 &b, OUT float9 &out);
    float dot_product_3x1(const float3 &a, const float3 &b);

    // These are utility functions specific to the math for the rods
    void get_p_i(const float3 &curr_r, const float3 &next_r, OUT float3 &p_i);
    void get_element_midpoint(const float3 &p_i, const float3 &r_i, OUT float3 &r_mid);
    void rodrigues_rotation(const float3 &v, const float3 &k, float theta, OUT float3 &v_rot);
    float safe_cos(float in);
    float get_l_i(const float3 &p_i, const float3 &p_im1);
    float get_signed_angle(const float3 &m1, const float3 &m2, const float3 &l);

    /*-----------------------*/
    /* Update Material Frame */
    /*-----------------------*/

    void perpendicularize(const float3 &m_i, const float3 &p_i, OUT float3 &m_i_prime);
    void update_m1_matrix(float3 &m_i, const float3 &p_i, const float3 &p_i_prime, float3 &m_i_prime);

    /*------------------*/
    /* Compute Energies */
    /*------------------*/

    float get_stretch_energy(float k, float3 &p_i, float3 &p_i_equil);
    void parallel_transport(float3 &m, float3 &m_prime, const float3 &p_im1, const float3 &p_i);
    float get_twist_energy(float beta, float3 &m_i, float3 &m_im1, float3 &m_i_equil, float3 &m_im1_equil, float3 &p_im1, float3 &p_i, float3 &p_im1_equil, float3 &p_i_equil);
    void get_kb_i(const float3 &p_im1, const float3 &p_i, OUT float3 &kb_i);
    void get_omega_j_i(const float3 &kb_i, const float3 &n_j, const float3 &m_j, OUT float2 &omega_j_i);
    float get_bend_energy(const float2 &omega_i_im1, const float2 &omega_i_im1_equil, const float4 &B_equil);

    float get_bend_energy_from_p(
        const float3 &p_im1,
        const float3 &p_i,
        const float3 &p_im1_equil,
        const float3 &p_i_equil,
        const float3 &n_im1,
        const float3 &m_im1,
        const float3 &n_im1_equil,
        const float3 &m_im1_equil,
        const float3 &n_i,
        const float3 &m_i,
        const float3 &n_i_equil,
        const float3 &m_i_equil,
        const float4 &B_i_equil,
        const float4 &B_im1_equil);

    float get_weights(const float3 &a, const float3 &b);
    void get_mutual_element_inverse(const float3 &pim1, const float3 &pi, float weight, OUT float3 &mutual_element);
    void get_mutual_axes_inverse(const float3 &mim1, const float3 &mi, float weight, OUT float3 &m_mutual);

    float get_bend_energy_mutual_parallel_transport(
        const float3 &p_im1,
        const float3 &p_i,
        const float3 &p_im1_equil,
        const float3 &p_i_equil,
        const float3 &n_im1,
        float3 &m_im1,
        const float3 &n_im1_equil,
        float3 &m_im1_equil,
        const float3 &n_i,
        float3 &m_i,
        const float3 &n_i_equil,
        float3 &m_i_equil,
        const float4 &B_i_equil,
        const float4 &B_im1_equil);

    /*----------*/
    /* Dynamics */
    /*----------*/

    float get_translational_friction(float viscosity, float radius, bool rotational);
    float get_rotational_friction(float viscosity, float radius, float length, bool safe);
    float get_force(float bend_energy, float stretch_energy, float delta_x);
    float get_torque(float twist_energy, float delta_theta);
    float get_delta_r(float friction, float timestep, float force, float noise, float external_force);
    float get_delta_r(float friction, float timestep, const std::vector<float> &forces);
    float get_delta_r(float friction, float timestep, const std::vector<float> &forces, float background_flow);
    float get_noise(float timestep, float kT, float friction, float random_number);

    /*------------*/
    /* Shorthands */
    /*------------*/

    void load_p(float4x3 &p, const std::vector<float> &r, int offset);
    void load_m(float4x3 &m_loaded, const std::vector<float> &m, int offset);
    void normalize_all(float4x3 &p);
    void absolute_all(const float4x3 &p, float4 &absolutes);
    void cross_all(const float4x3 &p, const float4x3 &m, OUT float4x3 &n);
    void delta_e_all(const float4x3 &e, const float4x3 &new_e, OUT float4x3 &delta_e);
    // void update_m1_all(float4x3 &m, float4 &absolutes, float4x3 &t, float4x3 &delta_e, OUT float4x3 &m_out); // This doesn't exist?
    void update_m1_matrix_all(float4x3 &m, const float4x3 &p, const float4x3 &p_prime, OUT float4x3 &m_prime, int start_cutoff, int end_cutoff);
    void fix_m1_all(float4x3 &m, float4x3 &new_t);
    void update_and_fix_m1_all(float4x3 &old_e, float4x3 &new_e, float4x3 &m);
    void set_cutoff_values(int e_i_node_no, int length, OUT int start_cutoff, int end_cutoff);

    /*----------*/
    /* Utility  */
    /*----------*/

    bool not_simulation_destroying(float x);
    bool not_simulation_destroying(const float3 &x);
    void load_B_all(float4x4 &B, const std::vector<float> &B_matrix, int offset);
    void make_diagonal_B_matrix(float B, OUT float4 &B_matrix);
    void set_cutoff_values(int e_i_node_no, int length, OUT int *start_cutoff, int *end_cutoff);
    float get_absolute_length_from_array(const std::vector<float> &array, int node_no, int length);
    void get_centroid_generic(const std::vector<float> &r, OUT float3 &centroid);

    /*-------------------------------*/
    /* Move the node, get the energy */
    /*-------------------------------*/

    void get_perturbation_energy(
        float perturbation_amount,
        int perturbation_dimension,
        std::vector<float> &B_matrix, //@todo This was previously named B_equil and len 4 in the header
        std::vector<float> &material_params,
        int start_cutoff,
        int end_cutoff,
        int p_i_node_no,
        std::vector<float> &r_all,
        std::vector<float> &r_all_equil,
        std::vector<float> &m_all,
        std::vector<float> &m_all_equil,
        float3 &energies);

}
#endif

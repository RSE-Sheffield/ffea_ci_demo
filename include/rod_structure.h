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
 *      rod_structure.h
 *	Author: Rob Welch, University of Leeds
 *	Email: py12rw@leeds.ac.uk

 *	Author: Ryan Cocking, University of Leeds
 *	Email: bsrctb@leeds.ac.uk
 */
#ifndef ROD_STRUCTURE
#define ROD_STRUCTURE

#include <iostream>
#include <fstream>
#include <cassert>
#include <boost/algorithm/string.hpp>
#include <string>
#include <vector>
#include <omp.h>
#include "RngStream.h"
#include <stdio.h>
#include "rod_interactions.h"

namespace rod
{

    std::vector<float> stof_vec(std::vector<std::string> vec_in, int length);
    InteractionData get_interaction_data(int elem_id_self, int elem_id_nbr, const std::vector<std::vector<InteractionData>> &nbr_list);
    struct Rod
    {
        /** Rod metadata **/
        int length;                   /** The length, L, of an array in the rod (x3 the number of nodes, N, in x, y and z) */
        int num_nodes;                /** The number of nodes in the rod */
        int num_rods;                 /** When the system is set up in ffeatools, this will be set */
        int rod_no;                   /** Each rod will be created with a unique ID */
        int line_start = 0;           /** Keeps track of whether the header file has been read */
        double rod_version = 999.999; /** If version number unknown, assume latest **/
        bool computed_rest_energy = false;
        int num_vdw_sites = 0;        /** The number of attractive van der Waals sites on the rod. Arranged by site index */

        /** Global simulation parameters - eventually read in from the .ffea file **/
        float viscosity = 0.6913 * pow(10, -3) / (mesoDimensions::pressure * mesoDimensions::time);  // denominator: poiseuille
        float timestep = 1e-12 / mesoDimensions::time;
        float kT = 0;
        float perturbation_amount = 0.001 * pow(10, -9) / mesoDimensions::length; /** Amount by which nodes are perturbed during numerical differentiation. May want to override with a local value depending on the scale of the simulation. **/
        int calc_noise = 0;
        int calc_steric = 0;  // Repulsive overlap of elements
        int calc_vdw = 0;   // Attractive protein-protein interactions
        int pbc = 0;
        float max_steric_energy = 50;    /** Potential energy when two rod elements are fully overlapped [energy units] **/
        std::string flow_profile;       // the type of background flow experienced by the rod (set by .ffea file)
        float flow_velocity[3] = {0};   // background flow imposed on the rod (set by .ffea file)
        float shear_rate = 0;           //

        float translational_friction;
        float rotational_friction; /** these will have to be changed when I end up implementing per-element radius, to be computed on-the-fly instead most likely **/

        /** Each set of rod data is stored in a single, c-style array, most of which go as {x,y,z,x,y,z...} */

        float *equil_r;                              /** Equilibrium configuration of the rod nodes. */
        float *equil_m;                              /** Equilibrium configuration of the material frame. */
        float *current_r;                            /** Current configuration of the rod nodes. */
        float *current_m;                            /** Current configuration of the material frame. */
        float *internal_perturbed_x_energy_positive; /** Energies associated with the perturbations we do in order to get dE. Given as [stretch, bend, twist, stretch, bend, twist...]**/
        float *internal_perturbed_y_energy_positive;
        float *internal_perturbed_z_energy_positive;
        float *internal_twisted_energy_positive;
        float *internal_perturbed_x_energy_negative; /** Stretch, bend, twist, stretch, bend, twist... */
        float *internal_perturbed_y_energy_negative;
        float *internal_perturbed_z_energy_negative;
        float *internal_twisted_energy_negative;
        float *material_params;                  /** Stretch, twist, radius, stretch, twist, radius... **/
        float *B_matrix;                         /** Contents of the bending modulus matrix for each node, as a 1-d array. Given as [a_1_1, a_1_2, a_2,1, a_2_2, a_1_1...]. **/
        float *steric_perturbed_energy_positive; // Length 2L array: energies from steric interactions at the start (0) and end (1) nodes of each rod element, i. Given as [xi0 yi0 zi0, xi1 yi1 zi1, ...]
        float *steric_perturbed_energy_negative;
        float *steric_energy;                    // Length L array : steric repulsion energy of elements (should be L/3 array, but previous decisions make it easier this way; xi=yi=zi, and so on)
        float *steric_force;                     // Length L array: steric repulsive force interpolated onto nodes from elements [x, y, z, ...]
        int *num_steric_nbrs; // Length L/3 array. Keeps track of how many neighbours each rod element has.
        float *vdw_energy;
        float *vdw_force;
        int *num_vdw_nbrs;
        std::vector<float> vdw_site_pos;     // Length = num_vdw_sites * 3. Must be a vector since it has to allow for there being no VDW sites on the rod.
        float *applied_forces;              /** Another [x,y,z,x,y,z...] array, this one containing the force vectors acting on each node in the rod. **/
        bool *pinned_nodes;                 /** This array is the length of the number of nodes in the rod, and it contains a boolean stating whether that node is pinned (true) or not (false). **/

        bool interface_at_start = false;    /** if this is true, the positioning of the start node is being handled by a rod-blob interface, so its energies will be ignored. **/
        bool interface_at_end = false;      /** if this is true, the positioning of the end node is being handled by a rod-blob interface, so its energies will be ignored. **/
        bool restarting = false;            /** If this is true, the rod will skip writing a frame of the trajectory (this is normally done so that the trajectory starts with correct box positioning) **/

        std::vector<std::vector<InteractionData>> steric_nbrs; /** Steric interaction neighbour list **/
        std::vector<std::vector<InteractionData>> vdw_nbrs;
        std::vector<VDWSite> vdw_sites;  // Attractive van der Waals binding sites that lie along the rod

        /** Unit conversion factors - the input\output files are in SI, but internally it uses FFEA's units as determined in dimensions.h **/
        float bending_response_factor;
        float spring_constant_factor;
        float twist_constant_factor;
        float viscosity_constant_factor;

        std::string rod_filename;
        std::string vdw_filename;
        FILE *file_ptr;
        int frame_no = 0;
        int step_no = 0;  // ! - redundant?
        Rod(int length, int set_rod_no);
        Rod(std::string path, int set_rod_no);
        Rod set_units();
        Rod compute_rest_energy();
        Rod do_timestep(std::shared_ptr<std::vector<RngStream>> &rng);
        Rod add_force(float force[4], int node_index);
        Rod pin_node(bool pin_state, int node_index);
        Rod load_header(std::string filename);
        Rod load_contents(std::string filename);
        Rod load_vdw(const std::string filename);
        Rod write_frame_to_file();
        Rod write_mat_params_array(float *array_ptr, int array_len, float stretch_scale_factor, float twist_scale_factor, float length_scale_factor);
        Rod change_filename(std::string new_filename);
        Rod equilibrate_rod(std::shared_ptr<std::vector<RngStream>> &rng);
        Rod translate_rod(float *r, float translation_vec[3]);
        Rod rotate_rod(float euler_angles[3]);
        Rod scale_rod(float scale);
        Rod get_centroid(float *r, float centroid[3]);
        Rod get_min_max(float *r, OUT float min[3], float max[3]);
        Rod get_p(int index, OUT float p[3], bool equil);
        Rod get_r(int node_index, OUT float r[3], bool equil);
        float get_radius(int node_index);
        float contour_length();
        float end_to_end_length();
        int get_num_nbrs(int element_index, const std::vector<std::vector<InteractionData>> &nbr_list);
        int get_num_vdw_sites();
        int get_num_nodes();
        Rod check_nbr_list_dim(std::vector<std::vector<InteractionData>> &nbr_list);
        void reset_nbr_list(std::vector<std::vector<InteractionData>> &nbr_list);
        Rod print_node_positions();
        std::vector<float> net_steric_force_nbrs(int elem_id);
        std::vector<float> net_vdw_force_nbrs(int elem_id);
        void do_steric();
        void do_vdw();
    };

} //end namespace
#endif

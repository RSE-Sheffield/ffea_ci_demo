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
 *      rod_structure.cpp
 *	Author: Rob Welch, University of Leeds
 *	Email: py12rw@leeds.ac.uk

 *	Author: Ryan Cocking, University of Leeds
 *	Email: bsrctb@leeds.ac.uk
 */

#include "rod_structure.h"

namespace rod
{

/** Easy access to 1-d arrays */
#define odx(x) x * 3
#define ody(y) (y * 3) + 1
#define odz(z) (z * 3) + 2
#define bend_index 0
#define stretch_index 1
#define twist_index 2

    /**---------**/
    /** Utility **/
    /**---------**/

    /**
    Convert a vector containing strings to a vector of floats.
    */
    std::vector<float> stof_vec(std::vector<std::string> vec_in)
    {
        std::vector<float> vec_out(vec_in.size());
        for (unsigned int i = 0; i < vec_in.size(); i++)
        {
            vec_out[i] = std::stof(vec_in[i]);
        }
        return vec_out;
    }

    /**
    Generate a random number between A and B, given an array of RngStream
    objects, and the id of the RngStream objects to be used.
    */
    float random_number(float A, float B, RngStream rng[], int thread_id)
    {
        return ((A) + ((B) - (A)) * (rng[thread_id].RandU01()));
    }

    /**-----**/
    /** Rod **/
    /**-----**/

    /**
    Create a rod of a known length. The rod_no is an arbitrary identifier.
    Note that this creates a rod without initialising the contents of the
    arrays.
    */
    Rod::Rod(int length, int set_rod_no) : equil_r(new float[length]),
                                           equil_m(new float[length]),
                                           current_r(new float[length]),
                                           current_m(new float[length]),
                                           internal_perturbed_x_energy_positive(new float[length]),
                                           internal_perturbed_y_energy_positive(new float[length]),
                                           internal_perturbed_z_energy_positive(new float[length]),
                                           internal_twisted_energy_positive(new float[length]),
                                           internal_perturbed_x_energy_negative(new float[length]),
                                           internal_perturbed_y_energy_negative(new float[length]),
                                           internal_perturbed_z_energy_negative(new float[length]),
                                           internal_twisted_energy_negative(new float[length]),
                                           material_params(new float[length]),
                                           B_matrix(new float[length + (length / 3)]),
                                           steric_perturbed_energy_positive(new float[2 * length]),
                                           steric_perturbed_energy_negative(new float[2 * length]),
                                           steric_force(new float[length]),
                                           applied_forces(new float[length + (length / 3)]),
                                           pinned_nodes(new bool[length / 3]){};

    /**
    Create a rod from a file. This won't do anything other than create the
    object - it won't even set up any arrays, because it doesn't know how
    long they'll have to be. After this, you have to call load_header and
    load_contents, which actually do the dirty work.
    */
    Rod::Rod(std::string path, int set_rod_no) : /** When we initialize from a file, we don't allocate arrays until we've loaded the file. **/
                                                 line_start(0)
    {
        rod_no = set_rod_no;
    };

    /**
    The contents of the rod, by default, are specified in SI units. Although
    it's possible to do everything in SI, you'll get more precision out of
    FFEA units. This function will convert all the units into FFEA units.
    When the file is written, the units are converted back automagically.
    The units are specified in mesoDimensions.h.
    */
    Rod Rod::set_units()
    {
        /** Translate our units into the units specified in FFEA's mesoDimensions header file **/
        bending_response_factor = pow(mesoDimensions::length, 4) * mesoDimensions::pressure;
        spring_constant_factor = mesoDimensions::force;
        twist_constant_factor = mesoDimensions::force * mesoDimensions::length * mesoDimensions::length;

        /** And now the rod itself **/
        for (int i = 0; i < length; i++)
        {
            equil_r[i] /= mesoDimensions::length;
            equil_m[i] /= mesoDimensions::length;
            current_m[i] /= mesoDimensions::length;
            current_r[i] /= mesoDimensions::length;
            internal_perturbed_x_energy_positive[i] /= mesoDimensions::Energy;
            internal_perturbed_y_energy_positive[i] /= mesoDimensions::Energy;
            internal_perturbed_z_energy_positive[i] /= mesoDimensions::Energy;
            internal_twisted_energy_positive[i] /= mesoDimensions::Energy;
            internal_perturbed_x_energy_negative[i] /= mesoDimensions::Energy;
            internal_perturbed_y_energy_negative[i] /= mesoDimensions::Energy;
            internal_perturbed_z_energy_negative[i] /= mesoDimensions::Energy;
            internal_twisted_energy_negative[i] /= mesoDimensions::Energy;

            if (i % 3 == 0)
            {
                material_params[i] /= spring_constant_factor;
            }

            if (i % 3 == 1)
            {
                material_params[i] /= twist_constant_factor;
            }

            if (i % 3 == 2)
            {
                material_params[i] /= mesoDimensions::length;
            }
        }

        for (int i = 0; i < length + length / 3; i++)
        {
            B_matrix[i] /= bending_response_factor;
        }

        for (int i = 0; i < 2 * length; i++)
        {
            steric_perturbed_energy_positive[i] /= mesoDimensions::Energy;
            steric_perturbed_energy_negative[i] /= mesoDimensions::Energy;
        }

        std::cout << "Set units of rod" << this->rod_no << "\n";

        return *this; /** Return a pointer to the object itself instead of void. Allows for method chaining! **/
    }

    /**---------**/
    /** Updates **/
    /**---------**/

    /**
    Do a timestep.
    This function contains three loops. Two over nodes and one over elements. The
    first loop (nodes) populates the contents of the energy arrays, which we use to
    work out delta E. The second one (elements) works out the energy from all steric
    interactions between neighbouring elements. The third loop (nodes) uses energies
    to compute dynamics and applies those dynamics to the position arrays.
    */
    Rod Rod::do_timestep(RngStream rng[])
    { // Most exciting method

        // if there is a rod-blob interface, this will avoid doing dynamics
        // on nodes which are attached to blobs
        int end_node = this->get_num_nodes();

        int node_min = 0;
        //if (this->interface_at_start){
        //    node_min = 1;
        //}

        //if (this->interface_at_end){
        //    end_node = get_num_nodes()-1;
        //}

        if (rod::dbg_print)
        {
            std::cout << "Rod: " << this->rod_no << "\n";
            std::cout << "Num nodes: " << this->get_num_nodes() << "\n"; //temp
            std::cout << "End node: " << end_node << "\n";               //temp
            std::cout << "Node min: " << node_min << "\n";
        } //temp

//The first loop is over all the nodes, and it computes all the energies for each one
#pragma omp parallel for schedule(dynamic) //most of the execution time is spent in this first loop
        for (int node_no = 0; node_no < end_node; node_no++)
        {

            //std::cout << "node " << node_no << "\n";

            if (rod::dbg_print)
            {
                std::cout << "IT BREAKS DURING NODE " << node_no << " ENERGY\n";
            }

            // if the node is pinned, we go to the next iteration of the loop (e.g. the next node)
            if (pinned_nodes[node_no] == true)
            {
                continue;
            }

            if (node_no < node_min)
            {
                continue;
            }

            // the cutoff values tell us how many nodes in our 5 node slice 'don't exist'
            // e.g. if we are at node i = n (at the end of the rod) the end cutoff will be 2
            int start_cutoff_val;
            int end_cutoff_val;
            int *start_cutoff = &start_cutoff_val;
            int *end_cutoff = &end_cutoff_val; // for the multiple return values
            set_cutoff_values(node_no, this->get_num_nodes(), start_cutoff, end_cutoff);

            // We need this e now because we need the previous value of e to do a material frame update
            // If you're curious about the [4][3] check out the get_perturbation_energy docs
            float p[4][3];
            load_p(p, current_r, node_no);

            float energies[3]; //bend, stretch, twist (temporary variable).

            // todo: these can have fewer arguments - maybe perturbation amount and cutoff values

            // We move the node backwards and forwards in each degree of freedom, so we end up calling get_perturbation_energy eight whole times
            // Fill the temporary variable with energies ( we basically pass the entire state of the rod to get_perturbation_energy)
            if (rod::dbg_print)
            {
                std::cout << "getting perturbation energy 1...\n";
            } //temp
            get_perturbation_energy(
                perturbation_amount * 0.5, //half one way, half the other
                x,                         // dimension (x, y, z are array indices, defined to be 1, 2, 3 at the top of this file, twist is = 4)
                B_matrix,
                material_params,
                start_cutoff_val,
                end_cutoff_val,
                node_no,
                current_r,
                equil_r,
                current_m,
                equil_m,
                energies); // this is a void function with energies as the output

            // transfer them from the temporary variable to the real thing
            internal_perturbed_x_energy_positive[node_no * 3] = energies[stretch_index];
            internal_perturbed_x_energy_positive[(node_no * 3) + 1] = energies[bend_index];
            internal_perturbed_x_energy_positive[(node_no * 3) + 2] = energies[twist_index];

            if (rod::dbg_print)
            {
                std::cout << "getting perturbation energy 2...\n";
            }                        //temp
            get_perturbation_energy( //from rod_math
                perturbation_amount * 0.5,
                y, // dimension
                B_matrix,
                material_params,
                start_cutoff_val,
                end_cutoff_val,
                node_no,
                current_r,
                equil_r,
                current_m,
                equil_m,
                energies);

            internal_perturbed_y_energy_positive[node_no * 3] = energies[stretch_index];
            internal_perturbed_y_energy_positive[(node_no * 3) + 1] = energies[bend_index];
            internal_perturbed_y_energy_positive[(node_no * 3) + 2] = energies[twist_index];

            get_perturbation_energy( //from rod_math
                perturbation_amount * 0.5,
                z, // dimension
                B_matrix,
                material_params,
                start_cutoff_val,
                end_cutoff_val,
                node_no,
                current_r,
                equil_r,
                current_m,
                equil_m,
                energies);

            internal_perturbed_z_energy_positive[node_no * 3] = energies[stretch_index];
            internal_perturbed_z_energy_positive[(node_no * 3) + 1] = energies[bend_index];
            internal_perturbed_z_energy_positive[(node_no * 3) + 2] = energies[twist_index];

            float twist_perturbation = 0.006283185; // 2pi/1000
            get_perturbation_energy(                //from rod_math
                twist_perturbation * 0.5,
                4, // twist dimension = 4
                B_matrix,
                material_params,
                start_cutoff_val,
                end_cutoff_val,
                node_no,
                current_r,
                equil_r,
                current_m,
                equil_m,
                energies);

            internal_twisted_energy_positive[node_no * 3] = energies[stretch_index];
            internal_twisted_energy_positive[(node_no * 3) + 1] = energies[bend_index];
            internal_twisted_energy_positive[(node_no * 3) + 2] = energies[twist_index];

            get_perturbation_energy(
                perturbation_amount * -0.5,
                x, // dimension
                B_matrix,
                material_params,
                start_cutoff_val,
                end_cutoff_val,
                node_no,
                current_r,
                equil_r,
                current_m,
                equil_m,
                energies);

            internal_perturbed_x_energy_negative[node_no * 3] = energies[stretch_index];
            internal_perturbed_x_energy_negative[(node_no * 3) + 1] = energies[bend_index];
            internal_perturbed_x_energy_negative[(node_no * 3) + 2] = energies[twist_index];

            get_perturbation_energy( //from rod_math
                perturbation_amount * -0.5,
                y, // dimension
                B_matrix,
                material_params,
                start_cutoff_val,
                end_cutoff_val,
                node_no,
                current_r,
                equil_r,
                current_m,
                equil_m,
                energies);

            internal_perturbed_y_energy_negative[node_no * 3] = energies[stretch_index];
            internal_perturbed_y_energy_negative[(node_no * 3) + 1] = energies[bend_index];
            internal_perturbed_y_energy_negative[(node_no * 3) + 2] = energies[twist_index];

            get_perturbation_energy( //from rod_math
                perturbation_amount * -0.5,
                z, // dimension
                B_matrix,
                material_params,
                start_cutoff_val,
                end_cutoff_val,
                node_no,
                current_r,
                equil_r,
                current_m,
                equil_m,
                energies);

            internal_perturbed_z_energy_negative[node_no * 3] = energies[stretch_index];
            internal_perturbed_z_energy_negative[(node_no * 3) + 1] = energies[bend_index];
            internal_perturbed_z_energy_negative[(node_no * 3) + 2] = energies[twist_index];

            get_perturbation_energy( //from rod_math
                twist_perturbation * -0.5,
                4, // twist dimension = 4
                B_matrix,
                material_params,
                start_cutoff_val,
                end_cutoff_val,
                node_no,
                current_r,
                equil_r,
                current_m,
                equil_m,
                energies);

            internal_twisted_energy_negative[node_no * 3] = energies[stretch_index];
            internal_twisted_energy_negative[(node_no * 3) + 1] = energies[bend_index];
            internal_twisted_energy_negative[(node_no * 3) + 2] = energies[twist_index];

        } // exit internal energy loop

        // Rod-rod interactions: loop over all elements, interpolate energies onto nodes
        if (rod::dbg_print)
        {
            std::cout << "rod.calc_steric_rod: " << this->calc_steric_rod << std::endl;
        }
        int end_elem = this->get_num_nodes() - 1;
        if (this->calc_steric_rod == 1)
        {
#pragma omp parallel for schedule(dynamic)
            for (int node_no = 0; node_no < end_elem; node_no++)
            {
                int num_neighbours = get_num_steric_neighbours(node_no);

                if (rod::dbg_print)
                {
                    std::cout << "element: " << node_no << std::endl;
                }

                for (int neighbour_no = 0; neighbour_no < num_neighbours; neighbour_no++)
                {
                    float r_a[3] = {0, 0, 0};
                    float p_a[3] = {0, 0, 0};
                    float radius_a = 0;
                    float radius_b = 0;
                    float c_a[3] = {0, 0, 0};
                    float c_b[3] = {0, 0, 0};
                    float energies[2] = {0, 0};

                    if (rod::dbg_print)
                    {
                        std::cout << "neighbour_no: " << neighbour_no << " of " << num_neighbours << std::endl;
                    }
                    get_r(node_no, r_a, false);
                    get_p(node_no, p_a, false);
                    radius_a = get_radius(node_no);
                    // ! - implement new neighbour style!

                    // positive x
                    rod::get_steric_perturbation_energy(
                        perturbation_amount * 0.5,
                        x,
                        steric_force_factor,
                        r_a,
                        p_a,
                        c_a,
                        c_b,
                        radius_a,
                        radius_b,
                        energies);

                    steric_perturbed_energy_positive[node_no * 3] += energies[0];
                    steric_perturbed_energy_positive[(node_no * 3) + 3] += energies[1];

                    // positive y
                    rod::get_steric_perturbation_energy(
                        perturbation_amount * 0.5,
                        y,
                        steric_force_factor,
                        r_a,
                        p_a,
                        c_a,
                        c_b,
                        radius_a,
                        radius_b,
                        energies);

                    steric_perturbed_energy_positive[(node_no * 3) + 1] += energies[0];
                    steric_perturbed_energy_positive[(node_no * 3) + 1 + 3] += energies[1];

                    // positive z
                    rod::get_steric_perturbation_energy(
                        perturbation_amount * 0.5,
                        z,
                        steric_force_factor,
                        r_a,
                        p_a,
                        c_a,
                        c_b,
                        radius_a,
                        radius_b,
                        energies);

                    steric_perturbed_energy_positive[(node_no * 3) + 2] += energies[0];
                    steric_perturbed_energy_positive[(node_no * 3) + 2 + 3] += energies[1];

                    // negative x
                    rod::get_steric_perturbation_energy(
                        perturbation_amount * -0.5,
                        x,
                        steric_force_factor,
                        r_a,
                        p_a,
                        c_a,
                        c_b,
                        radius_a,
                        radius_b,
                        energies);

                    steric_perturbed_energy_negative[(node_no * 3)] += energies[0];
                    steric_perturbed_energy_negative[(node_no * 3) + 3] += energies[1];

                    // negative y
                    rod::get_steric_perturbation_energy(
                        perturbation_amount * -0.5,
                        y,
                        steric_force_factor,
                        r_a,
                        p_a,
                        c_a,
                        c_b,
                        radius_a,
                        radius_b,
                        energies);

                    steric_perturbed_energy_negative[(node_no * 3) + 1] += energies[0];
                    steric_perturbed_energy_negative[(node_no * 3) + 1 + 3] += energies[1];

                    // negative z
                    rod::get_steric_perturbation_energy(
                        perturbation_amount * -0.5,
                        z,
                        steric_force_factor,
                        r_a,
                        p_a,
                        c_a,
                        c_b,
                        radius_a,
                        radius_b,
                        energies);

                    steric_perturbed_energy_negative[(node_no * 3) + 2] += energies[0];
                    steric_perturbed_energy_negative[(node_no * 3) + 2 + 3] += energies[1];

                } // exit steric loop

            } //exit interaction energy loop

        } //exit steric if-statement

        //This loop is for the dynamics
        for (int node_no = 0; node_no < end_node; node_no++)
        {

            // If the node is pinned, there's nothing to do
            if (pinned_nodes[node_no] == true)
            {
                continue;
            }

            if (node_no < node_min)
            {
                continue;
            }

// Grab thread ID from openMP (needed for RNG)
#ifdef USE_OPENMP
            int thread_id = omp_get_thread_num();
#else
            int thread_id = 0;
#endif

            // Get friction, needed for delta r and delta theta
            float translational_friction = get_translational_friction(this->viscosity, get_radius(node_no), false);
            //float length_for_friction = (get_absolute_length_from_array(equil_r, node_no, this->length) + get_absolute_length_from_array(equil_r, node_no-1, this->length))/2;
            float length_for_friction = 1e-8 / mesoDimensions::length;

            float rotational_friction = get_rotational_friction(this->viscosity, get_radius(node_no), length_for_friction, true);

            // Need these again
            float twist_perturbation = 0.006283185; // 2pi/1000

            // The material frame update requires that we grab the ith segment as it was before we did any dynamics
            float previous_p_i[3];
            float r_i[3];
            float r_ip1[3];
            r_i[0] = current_r[node_no * 3];
            r_i[1] = current_r[(node_no * 3) + 1];
            r_i[2] = current_r[(node_no * 3) + 2];
            r_ip1[0] = current_r[(node_no * 3) + 3];
            r_ip1[1] = current_r[(node_no * 3) + 4];
            r_ip1[2] = current_r[(node_no * 3) + 5];
            get_p_i(r_i, r_ip1, previous_p_i);

            // Get fluctuating force
            float x_noise = 0;
            float y_noise = 0;
            float z_noise = 0;
            float twist_noise = 0;
            if (this->calc_noise == 1)
            {
                x_noise = get_noise(timestep, kT, translational_friction, random_number(-0.5, 0.5, rng, thread_id));
                y_noise = get_noise(timestep, kT, translational_friction, random_number(-0.5, 0.5, rng, thread_id));
                z_noise = get_noise(timestep, kT, translational_friction, random_number(-0.5, 0.5, rng, thread_id));
                twist_noise = get_noise(timestep, kT, rotational_friction, random_number(-0.5, 0.5, rng, thread_id));
            }

            // Sum our energies and use them to compute the force
            // ! Not sure why we have negative - positive here, I'd have assumed it was the other way around
            float x_force = (internal_perturbed_x_energy_negative[node_no * 3] + internal_perturbed_x_energy_negative[(node_no * 3) + 1] + internal_perturbed_x_energy_negative[(node_no * 3) + 2] - (internal_perturbed_x_energy_positive[node_no * 3] + internal_perturbed_x_energy_positive[(node_no * 3) + 1] + internal_perturbed_x_energy_positive[(node_no * 3) + 2])) / perturbation_amount;
            float y_force = (internal_perturbed_y_energy_negative[node_no * 3] + internal_perturbed_y_energy_negative[(node_no * 3) + 1] + internal_perturbed_y_energy_negative[(node_no * 3) + 2] - (internal_perturbed_y_energy_positive[node_no * 3] + internal_perturbed_y_energy_positive[(node_no * 3) + 1] + internal_perturbed_y_energy_positive[(node_no * 3) + 2])) / perturbation_amount;
            float z_force = (internal_perturbed_z_energy_negative[node_no * 3] + internal_perturbed_z_energy_negative[(node_no * 3) + 1] + internal_perturbed_z_energy_negative[(node_no * 3) + 2] - (internal_perturbed_z_energy_positive[node_no * 3] + internal_perturbed_z_energy_positive[(node_no * 3) + 1] + internal_perturbed_z_energy_positive[(node_no * 3) + 2])) / perturbation_amount;
            float twist_force = (internal_twisted_energy_negative[node_no * 3] + internal_twisted_energy_negative[(node_no * 3) + 1] + internal_twisted_energy_negative[(node_no * 3) + 2] - (internal_twisted_energy_positive[node_no * 3] + internal_twisted_energy_positive[(node_no * 3) + 1] + internal_twisted_energy_positive[(node_no * 3) + 2])) / twist_perturbation;

            // Rod-rod steric interactions
            if (this->calc_steric_rod == 1)
            {

                float steric_positive[3] = {0};
                float steric_negative[3] = {0};

                float check = 1000;
                if (std::abs(steric_positive[0]) >= check or std::abs(steric_positive[1]) >= check or std::abs(steric_positive[2]) >= check)
                {
                    rod::print_array("steric_energy_positive", steric_positive, 3);
                    throw std::runtime_error("Steric energy exceeds 1000 kT");
                }
                if (std::abs(steric_negative[0]) >= check or std::abs(steric_negative[1]) >= check or std::abs(steric_negative[2]) >= check)
                {
                    rod::print_array("steric_energy_negative", steric_negative, 3);
                    throw std::runtime_error("Steric energy exceeds 1000 kT");
                }
            }

            // Get applied force, if any
            float applied_force_x = applied_forces[node_no * 4];
            float applied_force_y = applied_forces[(node_no * 4) + 1];
            float applied_force_z = applied_forces[(node_no * 4) + 2];
            float applied_force_twist = applied_forces[(node_no * 4) + 3];

            // Get delta r and delta twist
            float delta_r_x = get_delta_r(translational_friction, timestep, x_force, x_noise, applied_force_x); //from rod_math
            float delta_r_y = get_delta_r(translational_friction, timestep, y_force, y_noise, applied_force_y);
            float delta_r_z = get_delta_r(translational_friction, timestep, z_force, z_noise, applied_force_z);
            float delta_twist = get_delta_r(rotational_friction, timestep, twist_force, twist_noise, applied_force_twist);
            //std::cout << "delta_twist: " << delta_twist << "\n";

            // Apply our delta x
            current_r[node_no * 3] += delta_r_x;
            current_r[(node_no * 3) + 1] += delta_r_y;
            current_r[(node_no * 3) + 2] += delta_r_z;
            if (rod::dbg_print)
            {
                std::cout << "delta_r: [" << delta_r_x << ", " << delta_r_y << ", " << delta_r_z << "]" << std::endl;
            }

            // A wee sanity check to stop your simulations from exploding horribly
            float p[3] = {0, 0, 0};
            this->get_p(node_no, p, false);
            float check = rod::absolute(p) * int(1e3);
            if (std::abs(delta_r_x) >= check or std::abs(delta_r_y) >= check or std::abs(delta_r_z) >= check)
            {
                std::cout << "node " << node_no << "\n";
                std::cout << "dynamics: " << delta_r_x << ", " << delta_r_y << ", " << delta_r_z << "\n";
                rod_abort("Rod dynamics explosion. Bring your debugger.");
            }

            // If we're applying delta twist, we must load our new p_i back in
            if (node_no != this->get_num_nodes() - 1)
            { // The last node has no p_i, so it can't rotate
                float m_to_rotate[3];
                m_to_rotate[0] = current_m[node_no * 3];
                m_to_rotate[1] = current_m[(node_no * 3) + 1];
                m_to_rotate[2] = current_m[(node_no * 3) + 2]; // take the relevant info out of the data structure
                float p_i[3];
                float r_i[3];
                float r_ip1[3]; // we need all of these from the rod again
                r_i[0] = current_r[node_no * 3];
                r_i[1] = current_r[(node_no * 3) + 1];
                r_i[2] = current_r[(node_no * 3) + 2];
                r_ip1[0] = current_r[(node_no * 3) + 3];
                r_ip1[1] = current_r[(node_no * 3) + 1 + 3];
                r_ip1[2] = current_r[(node_no * 3) + 2 + 3];
                get_p_i(r_i, r_ip1, p_i);                                       //from rod_math
                rodrigues_rotation(m_to_rotate, p_i, delta_twist, m_to_rotate); // work out the actual rotated value
                current_m[node_no * 3] = m_to_rotate[0];
                current_m[(node_no * 3) + 1] = m_to_rotate[1];
                current_m[(node_no * 3) + 2] = m_to_rotate[2]; //put it back in the data structure
            }

            // If the element has moved, we need to update the material frame to have moved accordingly
            if (node_no != this->get_num_nodes() - 1)
            { // for last node index, material frame doesn't exist!
                float current_p_i[3];
                float m_to_fix[3];
                float m_i_prime[3];
                float r_i[3];
                float r_ip1[3]; //now: grab the quantities we need out of the data structure
                r_i[0] = current_r[node_no * 3];
                r_i[1] = current_r[(node_no * 3) + 1];
                r_i[2] = current_r[(node_no * 3) + 2];
                r_ip1[0] = current_r[(node_no * 3) + 3];
                r_ip1[1] = current_r[(node_no * 3) + 1 + 3];
                r_ip1[2] = current_r[(node_no * 3) + 2 + 3];
                for (int i = 0; i < 3; i++)
                {
                    current_p_i[i] = r_ip1[i] - r_i[i];
                }
                m_to_fix[0] = current_m[node_no * 3];
                m_to_fix[1] = current_m[(node_no * 3) + 1];
                m_to_fix[2] = current_m[(node_no * 3) + 2];
                update_m1_matrix(m_to_fix, previous_p_i, current_p_i, m_i_prime); //from rod_math - notice we're using the previous p_i we grabbed at the start of the function
                current_m[node_no * 3] = m_i_prime[0];
                current_m[(node_no * 3) + 1] = m_i_prime[1];
                current_m[(node_no * 3) + 2] = m_i_prime[2]; // back into the data structure you go
            }

            step_no += 1; //we just did one timestep so increment this
        }

        if (this->calc_steric_rod == 1)
        {
            this->reset_neighbour_list();
        }

        return *this;
    }

    /**----**/
    /** IO **/
    /**----**/

    /**
    Load the header info from a .rodtraj file (that's everything before the
    ---END HEADER--- line). Not all the info is read, some of it is for
    clarity. This populates some rod variables:
        - length - total length of the array (normally 3x the number of nodes)
        - num_nodes - number of nodes in the rod
        - num_rods - number of rods in the simulation. Not used right now.
        - line_start - number of the line at which the trajectory begins. This
        variable is used by load_contents later on, to skip the header.
        - version - which version of the algorithm is this made for?
    This method will also allocate the memory for all the arrays in the
    rod. Descriptions of those array are in rod_structure.h. Finally, it
    sets some default values for global simulation parameters. Eventually,
    these will be overwritten by parameters from the .ffea file.
    */
    Rod Rod::load_header(std::string filename)
    {
        rod_filename = filename;
        file_ptr = fopen(filename.c_str(), "a");

        /** This string denotes where the header info ends */
        const std::string rod_connections = "---END HEADER---";

        /** Check that we can load our input file */
        std::ifstream infile(filename);
        if (!infile)
        {
            std::cout << "Walrus IO error. Does input file exist?" << std::endl;
        }

        /** Iterate through file up to rod connections marker. */
        int n = 0;
        bool length_set = false; /** Use this to prevent rods with bad header info from being initialized **/
        for (std::string line; getline(infile, line);)
        {

            if (n == 0)
            {
                assert(line == "format,ffea_rod");
            } /** Check that format is valid FFEA_rod */
            if (n > 0 and line != rod_connections)
            {
                /** Extract data from lines and deposit it into object */
                std::vector<std::string> line_vec;
                boost::split(line_vec, line, boost::is_any_of(","));

                /** Read in the contents */
                if (line_vec[0] == "version")
                {
                    this->rod_version = std::stod(line_vec[1]);
                }
                if (line_vec[0] == "length")
                {
                    this->length = std::stoi(line_vec[1]);
                    length_set = true;
                }
                if (line_vec[0] == "num_elements")
                {
                    this->num_nodes = std::stoi(line_vec[1]);
                }
                if (line_vec[0] == "num_rods")
                {
                    this->num_rods = std::stoi(line_vec[1]);
                }
            }
            if (line == rod_connections)
            {
                this->line_start = n; /** Set this variable so other methods know we've read the header */
            }
            n++;
        }

        assert(length_set == true && "Length of the rod from file has not been set. File may be missing or formatted incorrectly.");
        assert(this->num_nodes != 0 && "Rod has zero nodes.");

        /** Warn the user if there is a file version mismatch */
        if (fabs(this->rod_version - rod_software_version) > 0.0000001)
        {
            std::cout << "Serious warning: The rod input was generated by a different version of the software to the one you are using. \n";
            std::cout << "Rod file version: " << this->rod_version << ", software version: " << rod_software_version << ".\n";
        }

        /** Now that we know the rod length, we can allocate the memory **/
        equil_r = static_cast<float *>(malloc(sizeof(float) * length));
        equil_m = static_cast<float *>(malloc(sizeof(float) * length));
        current_r = static_cast<float *>(malloc(sizeof(float) * length));
        current_m = static_cast<float *>(malloc(sizeof(float) * length));
        internal_perturbed_x_energy_positive = static_cast<float *>(malloc(sizeof(float) * length));
        internal_perturbed_y_energy_positive = static_cast<float *>(malloc(sizeof(float) * length));
        internal_perturbed_z_energy_positive = static_cast<float *>(malloc(sizeof(float) * length));
        internal_twisted_energy_positive = static_cast<float *>(malloc(sizeof(float) * length));
        internal_perturbed_x_energy_negative = static_cast<float *>(malloc(sizeof(float) * length));
        internal_perturbed_y_energy_negative = static_cast<float *>(malloc(sizeof(float) * length));
        internal_perturbed_z_energy_negative = static_cast<float *>(malloc(sizeof(float) * length));
        internal_twisted_energy_negative = static_cast<float *>(malloc(sizeof(float) * length));
        material_params = static_cast<float *>(malloc(sizeof(float) * length));
        B_matrix = static_cast<float *>(malloc(sizeof(float) * (length + (length / 3))));
        steric_perturbed_energy_positive = static_cast<float *>(malloc(sizeof(float) * 2 * length));
        steric_perturbed_energy_negative = static_cast<float *>(malloc(sizeof(float) * 2 * length));
        steric_force = static_cast<float *>(malloc(sizeof(float) * length));
        applied_forces = static_cast<float *>(malloc(sizeof(float) * (length + (length / 3))));
        pinned_nodes = static_cast<bool *>(malloc(sizeof(bool) * length / 3));
        steric_neighbours = std::vector<std::vector<InteractionData>>((length / 3) - 1);

        for (int i = 0; i < length / 3; i++)
        {
            pinned_nodes[i] = false;
        }

        for (int i = 0; i < length + (length / 3); i++)
        {
            applied_forces[i] = 0;
        }

        // These are hardcoded default values - the real values are loaded in from the .ffea script
        viscosity_constant_factor = mesoDimensions::pressure * mesoDimensions::time; ///poiseuille
        this->viscosity = 0.6913 * pow(10, -3) / viscosity_constant_factor;
        this->timestep = 1e-12 / mesoDimensions::time;
        this->kT = 0;
        this->perturbation_amount = 0.001 * pow(10, -9) / mesoDimensions::length; // todo: set this dynamically, maybe 1/1000 equilibrium length?
        this->calc_noise = 1;
        this->calc_steric_rod = 0;
        this->steric_force_factor = (10 * mesoDimensions::Energy / 1e-9) / mesoDimensions::force;  // maximum at 41 pN (Biotin/streptavidin ~ 160 pN)

        return *this;
    }

    /**
    Add a constant force that acts upon the rod every timestep.
    Parameters: force[4], a 4-element array that specifies a 3-D force
    vector, with the last element being a torque instead.
    Node_index: the index of the node to apply these forces on.
    Note: to remove the force, you must call this function again with 0s,
    or it will continue appyling the force.
    */
    Rod Rod::add_force(float force[4], int node_index)
    {
        this->applied_forces[node_index * 4] = force[0] / mesoDimensions::force;
        this->applied_forces[(node_index * 4) + 1] = force[1] / mesoDimensions::force;
        this->applied_forces[(node_index * 4) + 2] = force[2] / mesoDimensions::force;
        this->applied_forces[(node_index * 4) + 3] = force[3] / (mesoDimensions::force * mesoDimensions::length);
        return *this;
    }

    /**
    Load the current state of the rod. This function expects that load_header
    has already been called. This populates all of the already-initialised
    arrays containing the state of the rod. Note that it only contains the
    current state of the rod - the FFEA_rod python class is the only one
    that loads the rod trajectory.
    */
    Rod Rod::load_contents(std::string filename)
    {

        /** Make sure this method isn't called before loading header info */
        if (line_start == 0)
        {
            std::cout << "Rod file at " << filename << "was not found. \n";
            std::cout << "Rod version: " << this->rod_version << ". Length =  " << this->length << "\n";
            assert(line_start != 0 && "Rod header\rod file not found.");
        }

        std::ifstream infile(filename);
        int n = 0;
        int line_of_last_frame = 0;

        /** Get the line number of the last frame */
        for (std::string line; getline(infile, line);)
        {
            std::vector<std::string> line_vec;
            boost::split(line_vec, line, boost::is_any_of(" "));
            if (line_vec[0] == "FRAME")
            {
                line_of_last_frame = n;
            }
            n++;
        }

        /** Check that we got the last frame */
        assert(line_of_last_frame != 0 && "Could not find last frame of rod.");
        if (rod::dbg_print)
        {
            std::cout << "[rod] Last frame: line " << line_of_last_frame << "\n";
        }

        /** Seek back to the start of the file */
        infile.clear();
        infile.seekg(0, std::ios::beg);
        n = 0;

        for (std::string line; getline(infile, line);)
        {
            /** Find the last frame from the line number we got earlier */
            if (n > line_of_last_frame)
            {
                /** Convert each string into a vector<string> */
                std::vector<std::string> line_vec;
                boost::split(line_vec, line, boost::is_any_of(","));
                /** Then convert that into a vector<float> */
                std::vector<float> line_vec_float;
                line_vec_float = stof_vec(line_vec);

                /** Check we're not going to overflow and ruin someone's life when we write into the array*/
                assert(
                    (unsigned)length == line_vec_float.size() ||
                    (unsigned)length + (length / 3) == line_vec_float.size() ||
                    (unsigned)2 * length == line_vec_float.size()); //it's definitely fine to cast length

                /** Set our rod data arrays to the raw .data() from the vector. */
                if (n == line_of_last_frame + 1)
                {
                    for (int i = 0; i < length; i++)
                        equil_r[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 2)
                {
                    for (int i = 0; i < length; i++)
                        equil_m[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 3)
                {
                    for (int i = 0; i < length; i++)
                        current_r[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 4)
                {
                    for (int i = 0; i < length; i++)
                        current_m[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 5)
                {
                    for (int i = 0; i < length; i++)
                        internal_perturbed_x_energy_positive[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 6)
                {
                    for (int i = 0; i < length; i++)
                        internal_perturbed_y_energy_positive[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 7)
                {
                    for (int i = 0; i < length; i++)
                        internal_perturbed_z_energy_positive[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 8)
                {
                    for (int i = 0; i < length; i++)
                        internal_twisted_energy_positive[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 9)
                {
                    for (int i = 0; i < length; i++)
                        internal_perturbed_x_energy_negative[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 10)
                {
                    for (int i = 0; i < length; i++)
                        internal_perturbed_y_energy_negative[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 11)
                {
                    for (int i = 0; i < length; i++)
                        internal_perturbed_z_energy_negative[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 12)
                {
                    for (int i = 0; i < length; i++)
                        internal_twisted_energy_negative[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 13)
                {
                    for (int i = 0; i < length; i++)
                        material_params[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 14)
                {
                    for (int i = 0; i < length + (length / 3); i++)
                        B_matrix[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 15)
                {
                    for (int i = 0; i < 2 * length; i++)
                        steric_perturbed_energy_positive[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 16)
                {
                    for (int i = 0; i < 2 * length; i++)
                        steric_perturbed_energy_negative[i] = line_vec_float.data()[i];
                }
                if (n == line_of_last_frame + 17)
                {
                    for (int i = 0; i < length; i++)
                        steric_force[i] = line_vec_float.data()[i];
                }
            }
            n++;
        }
        std::cout << "Loaded contents of rod " << this->rod_no << " with " << this->get_num_nodes() << " nodes." << std::endl;

        return *this;
    }

    /**
    Write the current state of the rod to a file specified by the pointer
    *file_ptr. This will convert from MesoDimensions to SI units, so if your
    values are already in SI units, they'll be wrong.
    */
    Rod Rod::write_frame_to_file()
    {
        this->frame_no += 1;
        std::fprintf(file_ptr, "FRAME %i ROD %i\n", frame_no, rod_no);
        write_array(file_ptr, equil_r, length, mesoDimensions::length, true);
        write_array(file_ptr, equil_m, length, mesoDimensions::length, true);
        write_array(file_ptr, current_r, length, mesoDimensions::length, true);
        write_array(file_ptr, current_m, length, mesoDimensions::length, true);
        write_array(file_ptr, internal_perturbed_x_energy_positive, length, mesoDimensions::Energy, true);
        write_array(file_ptr, internal_perturbed_y_energy_positive, length, mesoDimensions::Energy, true);
        write_array(file_ptr, internal_perturbed_z_energy_positive, length, mesoDimensions::Energy, true);
        write_array(file_ptr, internal_twisted_energy_positive, length, mesoDimensions::Energy, true);
        write_array(file_ptr, internal_perturbed_x_energy_negative, length, mesoDimensions::Energy, true);
        write_array(file_ptr, internal_perturbed_y_energy_negative, length, mesoDimensions::Energy, true);
        write_array(file_ptr, internal_perturbed_z_energy_negative, length, mesoDimensions::Energy, true);
        write_array(file_ptr, internal_twisted_energy_negative, length, mesoDimensions::Energy, true);
        write_mat_params_array(material_params, length, spring_constant_factor, twist_constant_factor, mesoDimensions::length);
        write_array(file_ptr, B_matrix, length + (length / 3), bending_response_factor, true);
        write_array(file_ptr, steric_perturbed_energy_positive, 2 * length, mesoDimensions::Energy, true);
        write_array(file_ptr, steric_perturbed_energy_negative, 2 * length, mesoDimensions::Energy, true);
        write_array(file_ptr, steric_force, length, mesoDimensions::force, true);
        fflush(file_ptr);
        return *this;
    }

    /**
    This function is almost identical to the one above, but it appllies
    different scale factors for objects in the array,
    */
    Rod Rod::write_mat_params_array(float *array_ptr, int array_len, float stretch_scale_factor, float twist_scale_factor, float length_scale_factor)
    {
        float scale_factors[3] = {stretch_scale_factor, twist_scale_factor, length_scale_factor};
        for (int i = 0; i < array_len; i++)
        {
            if (i < array_len - 1)
            {
                std::fprintf(file_ptr, "%e,", array_ptr[i] * scale_factors[i % 3]);
            }
            else
            {
                std::fprintf(file_ptr, "%e", array_ptr[i] * scale_factors[i % 3]);
            }
        }
        std::fprintf(file_ptr, "\n");
        return *this;
    }

    /**
    Close the previous file and create a new file, assigning that to the rod
    variable *file_ptr. This will also copy the contents of the previous
    file into this one.
    */
    Rod Rod::change_filename(std::string new_filename)
    {
        /** Check if output file exists */
        std::ifstream out_test(new_filename);
        if (out_test.good())
        {
            this->rod_filename = new_filename;
            return *this;
        }

        /** Get contents of current file */
        std::string current_file_contents = "";
        std::ifstream infile(this->rod_filename);
        for (std::string line; getline(infile, line);)
        {
            if (line == "---END HEADER---")
            {
                break;
            }
            current_file_contents += line + "\n";
        }
        current_file_contents += "---END HEADER---\n";
        infile.close();

        /** Create new file **/
        std::ofstream outfile(new_filename);
        outfile << current_file_contents;
        outfile.close();

        /** Update member variables (filename, pointer etc) */
        this->rod_filename = new_filename;
        fclose(file_ptr);
        file_ptr = fopen(new_filename.c_str(), "a");
        return *this;
    }

    /**
    Run the simulation for an arbitrary amount of time. If you start a
    rod exactly in its equilibrium state, chances are it's not going to be
    equilibrated, which can throw off some tests. It runs for a totally
    arbitrary 1e-7 seconds and does not save the trajectory from the
    equilibration.
    */
    Rod Rod::equilibrate_rod(RngStream rng[])
    {
        int no_steps = 1e-7 / timestep; // this is arbitrary
        for (int i = 0; i < no_steps; i++)
        {
            this->do_timestep(rng);
        }
        return *this;
    }

    /**
    Translate every node in the rod by a given translation vector,
    translation_vec. The parameter float* r is the pointer to any array of
    node positions, e.g. this->current_r or this->equil_r. No return values,
    it just updates those arrays
    */
    Rod Rod::translate_rod(float *r, float translation_vec[3])
    {
        for (int i = 0; i < this->length; i += 3)
        {
            r[i] += translation_vec[0];
            r[i + 1] += translation_vec[1];
            r[i + 2] += translation_vec[2];
        }
        return *this;
    }

    /**
    Rotates the rod by the euler angles alpha, beta and gamma (or x, y
    and z if you prefer. This will update all the node positions AND the
    material frames will rotate as well. The rotations happen relative to
    each centroid, so if current_r and equil_r have different centroids,
    they will be rotated about different points.
    */
    Rod Rod::rotate_rod(float euler_angles[3])
    {
        /** Put rod centroid on 0,0,0 */
        float equil_centroid[3];
        get_centroid_generic(this->equil_r, this->length, equil_centroid);
        float equil_to_translate[3] = {-equil_centroid[0], -equil_centroid[1], -equil_centroid[2]};
        this->translate_rod(this->equil_r, equil_to_translate);

        float current_centroid[3];
        get_centroid_generic(this->current_r, this->length, current_centroid);
        float current_to_translate[3] = {-current_centroid[0], -current_centroid[1], -current_centroid[2]};
        this->translate_rod(this->current_r, current_to_translate);

        /** Construct rotation matrix from euler angles **/
        float xe = euler_angles[0];
        float ye = euler_angles[1];
        float ze = euler_angles[2];

        float Rx[9] = {1, 0, 0, 0, cosf(xe), -sinf(xe), 0, sinf(xe), cosf(xe)};
        float Ry[9] = {cosf(ye), 0, sinf(ye), 0, 1, 0, -sinf(ye), 0, cosf(ye)};
        float Rz[9] = {cosf(ze), -sinf(ze), 0, sinf(ze), cosf(ze), 0, 0, 0, 1};

        float RyRx[9];
        matmul_3x3_3x3(Ry, Rx, RyRx);
        float RzRyRx[9];
        matmul_3x3_3x3(Rz, RyRx, RzRyRx);

        /** Apply rotation matrix **/
        for (int i = 0; i < length; i += 3)
        {
            float temp_vec[3];
            /** do not try to improve this */
            apply_rotation_matrix(this->equil_r + i, RzRyRx, temp_vec);
            equil_r[i] = temp_vec[0];
            equil_r[i + 1] = temp_vec[1];
            equil_r[i + 2] = temp_vec[2];
            apply_rotation_matrix(this->current_r + i, RzRyRx, temp_vec);
            current_r[i] = temp_vec[0];
            current_r[i + 1] = temp_vec[1];
            current_r[i + 2] = temp_vec[2];
            apply_rotation_matrix(this->equil_m + i, RzRyRx, temp_vec);
            equil_m[i] = temp_vec[0];
            equil_m[i + 1] = temp_vec[1];
            equil_m[i + 2] = temp_vec[2];
            apply_rotation_matrix(this->current_m + i, RzRyRx, temp_vec);
            current_m[i] = temp_vec[0];
            current_m[i + 1] = temp_vec[1];
            current_m[i + 2] = temp_vec[2];
        }

        /** Move centroids back */
        this->translate_rod(this->current_r, current_centroid);
        this->translate_rod(this->equil_r, equil_centroid);

        return *this;
    }

    /**
     * Scale the rod by a float. No return values, it just updates the
     * arrays current_r and equil_r. It doesn't modify m, that'll be
     * normalized away anyway.
    */
    Rod Rod::scale_rod(float scale)
    {
        for (int i = 0; i < this->length; i += 3)
        {
            this->current_r[i] *= scale;
            this->current_r[i + 1] *= scale;
            this->current_r[i + 2] *= scale;
            this->equil_r[i] *= scale;
            this->equil_r[i + 1] *= scale;
            this->equil_r[i + 2] *= scale;
        }
        return *this;
    }

    /**
     * Get a centroid for the current frame of the rod. Note: you must supply
     * an array (either current_r or equil_r).
    */
    Rod Rod::get_centroid(float *r, OUT float centroid[3])
    {
        vec3d(n) { centroid[n] = 0; }
        for (int i = 0; i < this->length; i += 3)
        {
            centroid[0] += r[i];
            centroid[1] += r[i + 1];
            centroid[2] += r[i + 2];
        }
        vec3d(n) { centroid[n] /= this->get_num_nodes(); }
        return *this;
    }

    Rod Rod::get_min_max(float *r, OUT float min[3], float max[3])
    {
        for (int i = 0; i < this->length; i += 3)
        {
            if (r[i + x] > max[x])
            {
                max[x] = r[i + x];
            }
            if (r[i + y] > max[x])
            {
                max[y] = r[i + y];
            }
            if (r[i + z] > max[z])
            {
                max[z] = r[i + z];
            }
            if (r[i + x] < min[x])
            {
                min[x] = r[i + x];
            }
            if (r[i + y] < min[y])
            {
                min[y] = r[i + y];
            }
            if (r[i + z] < min[z])
            {
                min[z] = r[i + z];
            }
        }
        return *this;
    }

    /**
     * Get the rod element for the equilibrium or current structure, given
     * an element index.
     */
    Rod Rod::get_p(int index, OUT float p[3], bool equil)
    {
        if (equil)
        {
            vec3d(n) { p[n] = equil_r[(index * 3) + 3 + n] - equil_r[(index * 3) + n]; }
        }
        else
        {
            if (rod::dbg_print)
            {
                std::cout << "rod " << this->rod_no << ", index = " << index << "\n";
                float r_i[3] = {0, 0, 0};
                float r_ip1[3] = {0, 0, 0};
                vec3d(n) { r_i[n] = current_r[(index * 3) + n]; }
                vec3d(n) { r_ip1[n] = current_r[(index * 3) + 3 + n]; }
                rod::print_array("  r_i", r_i, 3);
                rod::print_array("  r_ip1", r_ip1, 3);
            }
            // current_r is 1D, so coordinates for rod nodes that are adjacent in space
            // are separated in the array by 3 items (x, y, and z)
            vec3d(n) { p[n] = current_r[(index * 3) + 3 + n] - current_r[(index * 3) + n]; }
        }
        if (rod::dbg_print)
        {
            print_array("  p", p, 3);
        }
        assert(std::abs(rod::absolute(p)) > 1e-7 && "Length of rod element is zero");
        return *this;
    }

    /**
     * Get the rod node position for the equilibrium or current structure, given
     * a node index.
     */
    Rod Rod::get_r(int node_index, OUT float r_i[3], bool equil)
    {
        if (equil)
        {
            vec3d(n) { r_i[n] = equil_r[(node_index * 3) + n]; }
        }
        else
        {
            vec3d(n) { r_i[n] = current_r[(node_index * 3) + n]; }
        }
        return *this;
    }

    float Rod::get_radius(int node_index)
    {
        return material_params[(node_index * 3) + 2];
    }

    int Rod::get_num_steric_neighbours(int element_index)
    {
        return steric_neighbours.at(element_index).size();
    }

    int Rod::get_num_nodes()
    {
        if (rod::dbg_print)
        {
            if (this->num_nodes == 0)
            {
                throw std::invalid_argument("Rod cannot have zero nodes.");
            }
        }
        return this->num_nodes;
    }


    Rod Rod::check_neighbour_list_dimensions()
    {

        int num_rows = steric_neighbours.size();
        int num_cols = 0;
        std::string msg;

        msg = "Number of rows in neighbour list (" + std::to_string(num_rows) + ") should be equal to number of rod elements (" + std::to_string(this->get_num_nodes() - 1) + ").";
        if (num_rows != this->get_num_nodes() - 1)
        {
            throw std::logic_error(msg);
        }

        if (rod::dbg_print)
        {
            std::cout << "Neighbour list dimensions of rod " << this->rod_no << " OK" << std::endl;
        }
        return *this;
    }

    InteractionData Rod::get_interaction_data(int rod_id_self, int rod_id_neighb,
        int elem_id_self, int elem_id_neighb)
    {
        InteractionData data = this->steric_neighbours.at(elem_id_self).at(elem_id_neighb);

        if (rod::dbg_print)
        {
            bool rod_ok = (rod_id_self == data.rod_id_self) and (rod_id_neighb = data.rod_id_neighb);
            bool elem_ok = (elem_id_self == data.element_id_self) and (elem_id_neighb = data.element_id_neighb);
            if (rod_ok and elem_ok)
            {
                return data;
            }
            else
            {
                std::cout << rod_id_self << ", " << data.rod_id_self << "; "
                          << rod_id_neighb << ", " << data.rod_id_neighb << "; "
                          << elem_id_self << ", " << data.element_id_self << "; "
                          << elem_id_neighb << ", " << data.element_id_neighb << "\n";
                throw std::runtime_error("Inconsistency between vector indexing and struct contents");
            }
        }
        else
        {
            return data;
        }
    }

    Rod Rod::reset_neighbour_list()
    {
        for (int i = 0; i < this->get_num_nodes() - 1; i++)
        {
            this->steric_neighbours.at(i).clear();
        }
        if (rod::dbg_print)
        {
            std::cout << "Reset neighbour list of rod " << this->rod_no << std::endl;
        }
        return *this;
    }

    // Just a silly debug function that prints all the positional data of the rod
    Rod Rod::print_node_positions()
    {
        float r[3] = {0, 0, 0};
        for (int i = 0; i < this->get_num_nodes(); i++)
        {
            std::string msg = "rod " + std::to_string(this->rod_no) + " r" + std::to_string(i);
            vec3d(n) { r[n] = this->current_r[(i * 3) + n]; }
            rod::print_array(msg, r, 3);
        }
        return *this;
    }

    /**
     * @brief Steric repulsive force on a rod element, a, due to its collision
     * with a neighbouring element, b.
     */
    std::array<float, 3> Rod::element_steric_force(int element_index, float radius_neighbour,
        float contact_a[3], float contact_b[3])
    {
        std::array<float, 6> energy = { 0 };  // +x, -x, +y, -y, +z, -z
        std::array<float, 3> force = { 0 };   // x, y, z
        float delta = this->perturbation_amount;
        float force_strength = this->steric_force_factor;
        float radius_sum = this->get_radius(element_index) + radius_neighbour;

        energy[0] = element_energy_from_perturbation(0,  delta, force_strength, contact_a, contact_b, radius_sum);
        energy[1] = element_energy_from_perturbation(0, -delta, force_strength, contact_a, contact_b, radius_sum);
        energy[2] = element_energy_from_perturbation(1,  delta, force_strength, contact_a, contact_b, radius_sum);
        energy[3] = element_energy_from_perturbation(1, -delta, force_strength, contact_a, contact_b, radius_sum);
        energy[4] = element_energy_from_perturbation(2,  delta, force_strength, contact_a, contact_b, radius_sum);
        energy[5] = element_energy_from_perturbation(2, -delta, force_strength, contact_a, contact_b, radius_sum);

        vec3d(n) { force[n] = -(energy[n] - energy[n + 1]) / delta; }

        return force;
    }


} //end namespace

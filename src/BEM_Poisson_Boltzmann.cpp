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
 *      BEM_Poisson_Boltzmann.h
 *	Author: Robin Richardson, University of Leeds
 *	Email: pyrar@leeds.ac.uk
 */

#define _USE_MATH_DEFINES
#include <cmath>

#include "FFEA_return_codes.h"
#include "BEM_Poisson_Boltzmann.h"
#include "NearestNeighbourLinkedListCube.h"
#include "SparseMatrixUnknownPattern.h"

#include "GaussianQuadrature_tri.h"
#include "GaussianQuadrature_1d.h"

BEM_Poisson_Boltzmann::BEM_Poisson_Boltzmann() {
    num_faces = 0;
    lookup = nullptr;
    kappa = 0;
}

BEM_Poisson_Boltzmann::~BEM_Poisson_Boltzmann() {
    num_faces = 0;
    lookup = nullptr;
    kappa = 0;
    
    mat_C.reset();
    mat_D.reset();
}

void BEM_Poisson_Boltzmann::init(NearestNeighbourLinkedListCube *lookup) {
    this->lookup = lookup;
    this->num_faces = lookup->get_pool_size();

    /* Create and initialise our sparse matrices */
    mat_C = std::make_unique<SparseMatrixUnknownPattern>();
    if (!mat_C) throw FFEAException("Could not allocate C matrix");
    mat_C->init(num_faces, 100);

    mat_D = std::make_unique<SparseMatrixUnknownPattern>();
    if (!mat_D) throw FFEAException("Could not allocate D matrix");
    mat_D->init(num_faces, 100);
}

/*
 * Sets the inverse debye screening length for the system
 */
void BEM_Poisson_Boltzmann::set_kappa(scalar kappa) {
    this->kappa = kappa;
}

/*
 *
 */
void BEM_Poisson_Boltzmann::build_BEM_matrices() {
    std::array<arr3, 4> gqp = {};

    /* Clear the C and D matrices */
    mat_C->zero();
    mat_D->zero();

    /* For each face, calculate the interaction with all other relevant faces and add the contribution to mat_C and mat_D */
    for (int i = 0; i < num_faces; i++) {

        // get the ith face
        LinkedListNode<Face> *l_i = lookup->get_from_pool(i);
        Face *f = l_i->obj;

        // Create matrix C diagonal (self term) for constant element case
        mat_C->set_diagonal_element(i, -.5 * 4.0 * M_PI);

        // Create matrix D diagonal (self term) for constant element case
        mat_D->set_diagonal_element(i, -(self_term(f->centroid, f->n[1]->pos, f->n[2]->pos, 6) +
                self_term(f->centroid, f->n[2]->pos, f->n[0]->pos, 6) +
                self_term(f->centroid, f->n[0]->pos, f->n[1]->pos, 6)));


        // calculate the 4 gauss points for this face
        f->barycentric_calc_point(1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, gqp[0]);
        f->barycentric_calc_point(.6, .2, .2, gqp[1]);
        f->barycentric_calc_point(.2, .6, .2, gqp[2]);
        f->barycentric_calc_point(.2, .2, .6, gqp[3]);

        // Perform the necessary integrals for the face f with all faces in its own cell (in the lookup grid)
        perform_integrals_for_lookup_cell_self(l_i, gqp);

        // Perform integrals for the 9 adjacent cells above and below the plane containg f's cell
        for (int x = -1; x < 2; x++)
            for (int z = -1; z < 2; z++) {
                perform_integrals_for_lookup_cell_relative(l_i, gqp, x, -1, z);
                perform_integrals_for_lookup_cell_relative(l_i, gqp, x, +1, z);
            }

        // Perform integrals for the remaining 8 adjacent cells, lying in the y-plane of f's cell
        for (int x = -1; x < 2; x++) {
            perform_integrals_for_lookup_cell_relative(l_i, gqp, x, 0, -1);
            perform_integrals_for_lookup_cell_relative(l_i, gqp, x, 0, +1);
        }
        perform_integrals_for_lookup_cell_relative(l_i, gqp, -1, 0, 0);
        perform_integrals_for_lookup_cell_relative(l_i, gqp, +1, 0, 0);
    }
}

void BEM_Poisson_Boltzmann::perform_integrals_for_lookup_cell_self(const LinkedListNode<Face> *l_i, std::array<arr3, 4> &gqp) {
    scalar mat_C_contribution, mat_D_contribution;

    // Get the top of the stack in the cell in which face f lies
    LinkedListNode<Face> *l_j = lookup->get_top_of_stack(l_i->x, l_i->y, l_i->z);
    while (l_j != nullptr) {
        if (l_j->index != l_i->index) {
            gauss_quadrature_4_point(gqp,
                    l_j->obj->centroid,
                    mat_D_contribution,
                    mat_C_contribution,
                    l_i->obj);

            mat_C->add_off_diagonal_element(l_i->index, l_j->index, mat_C_contribution);
            mat_D->add_off_diagonal_element(l_i->index, l_j->index, mat_D_contribution);
        }
        l_j = l_j->next;
    }
}

void BEM_Poisson_Boltzmann::perform_integrals_for_lookup_cell_relative(const LinkedListNode<Face> *l_i, std::array<arr3, 4>& gqp, int dx, int dy, int dz) {
    scalar mat_C_contribution, mat_D_contribution;

    LinkedListNode<Face> *l_j = lookup->get_top_of_stack(l_i->x + dx, l_i->y + dy, l_i->z + dz);
    while (l_j != nullptr) {
        gauss_quadrature_4_point(gqp,
                l_j->obj->centroid,
                mat_D_contribution,
                mat_C_contribution,
                l_i->obj);

        mat_C->add_off_diagonal_element(l_i->index, l_j->index, mat_C_contribution);
        mat_D->add_off_diagonal_element(l_i->index, l_j->index, mat_D_contribution);

        l_j = l_j->next;
    }
}

void BEM_Poisson_Boltzmann::print_matrices() {
    printf("\n\nmat_C:\n\n");
    mat_C->print();

    printf("\n\nmat_D:\n\n");
    mat_D->print();
}

std::unique_ptr<SparseMatrixUnknownPattern> &BEM_Poisson_Boltzmann::get_C() {
    return mat_C;
}

std::unique_ptr<SparseMatrixUnknownPattern> &BEM_Poisson_Boltzmann::get_D() {
    return mat_D;
}

/* Returns the value of the fundamental solution u multiplied by 4*pi */
scalar BEM_Poisson_Boltzmann::u_4pi(scalar r) {
    return exp(-kappa * r) / r;
}

/* Returns the radial component of grad of u multiplied by 4*pi (all other components are zero) */
scalar BEM_Poisson_Boltzmann::grad_u_4pi(scalar r, scalar r2) {
    return -exp(-kappa * r) * (1.0 / r2 + kappa / r);
}

/*
                scalar BEM_Poisson_Boltzmann::screened_R_theta(scalar r_perp_mag, scalar half_theta_max, scalar theta_bar, scalar xi)
                {
                        return exp( (-kappa * r_perp_mag) / cos(half_theta_max * (xi + 1) - theta_bar) );
                }
 */

void BEM_Poisson_Boltzmann::gauss_quadrature_4_point(std::array<arr3, 4> &gqp, arr3 &p, scalar &int_u, scalar &int_du, Face *f) {
    std::array<scalar, 4> r_hat_dot_n, r2, r_mag;
    arr3 r;
    
    for (int i = 0; i < 4; i++) {
        sub(gqp[i], p, r);
        r2[i] = magnitude2(r);
        r_mag[i] = sqrt(r2[i]);
        r_hat_dot_n[i] = dot(r, f->normal) / r_mag[i];
    }

    int_u = -f->area * (-.5625 * u_4pi(r_mag[0]) + 0.5208333333 * (u_4pi(r_mag[1]) + u_4pi(r_mag[2]) + u_4pi(r_mag[3])));
    int_du = f->area * (-.5625 * grad_u_4pi(r_mag[0], r2[0]) * r_hat_dot_n[0] + 0.5208333333 * (
            grad_u_4pi(r_mag[1], r2[1]) * r_hat_dot_n[1] +
            grad_u_4pi(r_mag[2], r2[2]) * r_hat_dot_n[2] +
            grad_u_4pi(r_mag[3], r2[3]) * r_hat_dot_n[3]));
}

scalar BEM_Poisson_Boltzmann::self_term(const arr3 &n0, const arr3 &n1, const arr3 &n2, int precision) {
    // Calculate the triangle side vectors, r_01, r_02 and r_12
    arr3 r_01, r_02, r_12;
    sub(n1, n0, r_01);
    sub(n2, n0, r_02);
    sub(n2, n1, r_12);

    // Get the lengths of these vectors
    const scalar r_01_sq = magnitude2(r_01);
    const scalar r_02_sq = magnitude2(r_02);
    const scalar r_12_sq = magnitude2(r_12);

    // Calculate some dot products
    const scalar r_01_dot_r_02 = dot(r_01, r_02);
    const scalar r_01_dot_r_12 = dot(r_01, r_12);

    // Get the size of the angle through which we need to integrate
    const scalar theta_max = acos(r_01_dot_r_02 / sqrt(r_01_sq * r_02_sq));

    // Get the vector between node 0 and the point perpendicularly opposite (on line node 1 to 2)
    arr3 r_0_perp;
    for (int i = 0; i < 3; ++i) {
        r_0_perp[i] = r_01[i] - r_01_dot_r_12 * r_12[i] / r_12_sq;
    }

    // Get the perpendicular distance between node 0 and line 1 to 2
    const scalar L_perp = magnitude(r_0_perp);

    // Get the angle between side 01 and the perpedicular line between 0 and side 12.
    const scalar theta_star = acos(L_perp / sqrt(r_01_sq));

    return integrate_function_1d_tri(theta_max, L_perp, theta_star, precision);
}

scalar BEM_Poisson_Boltzmann::f_1d(scalar r) {
    return exp(-kappa * r);
}

scalar BEM_Poisson_Boltzmann::f_3d(const arr3 &p, const arr3 &q) {
    scalar dx = p[0] - q[0];
    scalar dy = p[1] - q[1];
    scalar dz = p[2] - q[2];
    scalar r = sqrt(dx * dx + dy * dy + dz * dz);
    return exp(-kappa * r) / r;
}


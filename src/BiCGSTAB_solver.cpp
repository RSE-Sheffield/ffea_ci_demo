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

#include "BiCGSTAB_solver.h"

#include "mat_vec_fns_II.h"

void BiCGSTAB_solver::init(const int N, const scalar _tol, const int _max_num_iterations) {
    if (_max_num_iterations == 0) {
        throw FFEAException("BiCGSTAB_solver::max_num_iterations should not be set to zero.");
    }
    this->tol = _tol;
    this->max_num_iterations = _max_num_iterations;
    
    // Allocate all required memory
    try {
        inv_M = std::vector<scalar>(N, 0);
        r = std::vector<scalar>(N, 0);
        r_hat = std::vector<scalar>(N, 0);
        p = std::vector<scalar>(N, 0);
        p_hat = std::vector<scalar>(N, 0);
        q = std::vector<scalar>(N, 0);
        s = std::vector<scalar>(N, 0);
        s_hat = std::vector<scalar>(N, 0);
        t = std::vector<scalar>(N, 0);
    } catch (std::bad_alloc &) {
        throw FFEAException("While initialising BiCGSTAB_solver, could not allocate memory for vectors.");
    }
}

void BiCGSTAB_solver::solve(std::unique_ptr<SparseMatrixUnknownPattern> &A, std::vector<scalar> &x, std::vector<scalar> &b, const int num_iterations) {
    scalar rho_last = 1, alpha = 1, omega = 1, rho, beta;

    // Get the inverse of the diagonal of the matrix to use as a preconditioner
    A->calc_inverse_diagonal(inv_M);

    // Get the initial residual vector using a first guess at the solution x (r_0 = b - A x_0)
    get_residual_vector(r, b, A, x);

    // r_hat_0 <- r_0
    store(r, r_hat);

    // initialise q vector to zero
    initialise(q);

    // If num_iterations==0, we use max_num_iterations mode and check for convergence
    const int do_iterations = num_iterations == 0 ? max_num_iterations : num_iterations;
    
    // loop till convergence (or a threshold maximum number of iterations have elapsed -> convergence failure)
    for (int i = 0; i < do_iterations; i++) {
        // rho_i = dot(r_hat_0, r_i_minus_1);
        rho = dot(r_hat, r);

        if (rho == 0) {
            throw FFEAException("In BiCGSTAB_solver solve(), rho is zero. Solver stopping on iteration %d.", i);
        }

        // beta = (rho_i/rho_i_minus_1) (alpha/omega_i_minus_1)
        beta = (rho / rho_last) * (alpha / omega);

        // p = r + beta * (p - omega * q);
        complicated_machine(p, r, beta, p, -omega, q);

        // apply preconditioner to p (p_hat = inv_M * p_i)
        prod(inv_M, p, p_hat);

        // q = A p_hat
        A->apply(p_hat, q);

        // alpha = rho_i / dot(r_hat_0, q)
        alpha = rho / dot(r_hat, q);

        // s = r - alpha * q;
        scalar_vector_add(s, r, -alpha, q);

        // x = x + alpha * p_hat
        scalar_vector_add(x, x, alpha, p_hat);

        // Check for convergence (if max_num_iterations mode)
        if (num_iterations == 0 && dot(s, s) < tol) {
            return;
        }

        // apply preconditioner to s (s_hat = inv_M * s)
        prod(inv_M, s, s_hat);

        // t = A s_hat
        A->apply(s_hat, t);

        // omega_i = dot(t, s) / dot(t, t)
        omega = dot(t, s) / dot(t, t);

        // x = x + omega * s_hat;
        scalar_vector_add(x, x, omega, s_hat);

        // r = s - omega * t;
        scalar_vector_add(r, s, -omega, t);

        //		printf("BiCGSTAB_solver: %f\n", dot(s,s,N));
    }
    if (num_iterations == 0) {
        throw FFEAException("Bi-Conjugate Gradient Stabilised solver could not converge in %d iterations.", do_iterations);
    }
}

void BiCGSTAB_solver::get_residual_vector(std::vector<scalar> &r, std::vector<scalar> &b, std::unique_ptr<SparseMatrixUnknownPattern> &A, std::vector<scalar> &x) {
    // Ax
    A->apply(x, r);

    // r = b - Ax
    sub(b, r, r);
}

/* Carries out the operation x = y + c*z, where c is a scalar, and x, y and z are vectors of length N. */
void BiCGSTAB_solver::scalar_vector_add(std::vector<scalar> &x, const std::vector<scalar> &y, scalar c, const std::vector<scalar> &z) {
    for (unsigned int i = 0; i < x.size(); i++)
        x[i] = y[i] + c * z[i];
}

/* Carries out the operation w = x + a * (y + b * z) */
// p = r + beta(p - omega * v)
void BiCGSTAB_solver::complicated_machine(std::vector<scalar> &w, const std::vector<scalar> &x, scalar a, const std::vector<scalar> &y, scalar b, const std::vector<scalar> &z) {
    for (unsigned int i = 0; i < w.size(); i++)
        w[i] = x[i] + a * (y[i] + b * z[i]);
}

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

BiCGSTAB_solver::BiCGSTAB_solver() {
    N = 0;
    tol = 0;
    max_num_iterations = 0;
}

BiCGSTAB_solver::~BiCGSTAB_solver() {
    inv_M.clear();
    r.clear();
    r_hat.clear();
    p.clear();
    p_hat.clear();
    q.clear();
    s.clear();
    s_hat.clear();
    t.clear();

    N = 0;
    tol = 0;
    max_num_iterations = 0;
}

void BiCGSTAB_solver::init(int N, scalar tol, int max_num_iterations) {
    this->N = N;
    this->tol = tol;
    this->max_num_iterations = max_num_iterations;

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

void BiCGSTAB_solver::solve(std::unique_ptr<SparseMatrixUnknownPattern> &A, std::vector<scalar> &x, std::vector<scalar> &b) {
    scalar rho_last = 1, alpha = 1, omega = 1, rho, beta;

    // Get the inverse of the diagonal of the matrix to use as a preconditioner
    A->calc_inverse_diagonal(inv_M);

    // Get the initial residual vector using a first guess at the solution x (r_0 = b - A x_0)
    get_residual_vector(r, b, A, x, N);

    // r_hat_0 <- r_0
    copy_vector(r_hat, r, N);

    // initialise q vector to zero
    zero(q, N);

    // loop till convergence (or a threshold maximum number of iterations have elapsed -> convergence failure)
    for (int i = 0; i < max_num_iterations; i++) {

        // rho_i = dot(r_hat_0, r_i_minus_1);
        rho = dot(r_hat, r, N);

        if (rho == 0) {
            throw FFEAException("In BiCGSTAB_solver solve(), rho is zero. Solver stopping on iteration %d.", i);
        }

        // beta = (rho_i/rho_i_minus_1) (alpha/omega_i_minus_1)
        beta = (rho / rho_last) * (alpha / omega);

        // p = r + beta * (p - omega * q);
        complicated_machine(p, r, beta, p, -omega, q, N);

        // apply preconditioner to p (p_hat = inv_M p_i)
        apply_diagonal_matrix(p_hat, inv_M, p, N);

        // q = A p_hat
        A->apply(p_hat, q);

        // alpha = rho_i / dot(r_hat_0, q)
        alpha = rho / dot(r_hat, q, N);

        // s = r - alpha * q;
        scalar_vector_add(s, r, -alpha, q, N);

        // x = x + alpha * p_hat
        scalar_vector_add(x, x, alpha, p_hat, N);

        // Check for convergence
        if (dot(s, s, N) < tol) {
            //			printf("Convergence reached on iteration %d.\n", i);
            return;
        }

        // apply preconditioner to s (s_hat = inv_M s)
        apply_diagonal_matrix(s_hat, inv_M, s, N);

        // t = A s_hat
        A->apply(s_hat, t);

        // omega_i = dot(t, s) / dot(t, t)
        omega = dot(t, s, N) / dot(t, t, N);

        // x = x + omega * s_hat;
        scalar_vector_add(x, x, omega, s_hat, N);

        // r = s - omega * t;
        scalar_vector_add(r, s, -omega, t, N);

        //		printf("BiCGSTAB_solver: %f\n", dot(s,s,N));
    }

    throw FFEAException("Bi-Conjugate Gradient Stabilised solver could not converge in max_num_iterations.");
}

void BiCGSTAB_solver::solve(std::unique_ptr<SparseMatrixUnknownPattern> &A, std::vector<scalar> &x, std::vector<scalar> &b, int num_iterations) {
    scalar rho_last = 1, alpha = 1, omega = 1, rho, beta;

    // Get the inverse of the diagonal of the matrix to use as a preconditioner
    A->calc_inverse_diagonal(inv_M);

    // Get the initial residual vector using a first guess at the solution x (r_0 = b - A x_0)
    get_residual_vector(r, b, A, x, N);

    // r_hat_0 <- r_0
    copy_vector(r_hat, r, N);

    // initialise q vector to zero
    zero(q, N);

    // loop till convergence (or a threshold maximum number of iterations have elapsed -> convergence failure)
    int i;
    for (i = 0; i < num_iterations; i++) {

        // rho_i = dot(r_hat_0, r_i_minus_1);
        rho = dot(r_hat, r, N);

        if (rho == 0) {
            throw FFEAException("In BiCGSTAB_solver solve(), rho is zero. Solver stopping on iteration %d.", i);
        }

        // beta = (rho_i/rho_i_minus_1) (alpha/omega_i_minus_1)
        beta = (rho / rho_last) * (alpha / omega);

        // p = r + beta * (p - omega * q);
        complicated_machine(p, r, beta, p, -omega, q, N);

        // apply preconditioner to p (p_hat = inv_M p_i)
        apply_diagonal_matrix(p_hat, inv_M, p, N);

        // q = A p_hat
        A->apply(p_hat, q);

        // alpha = rho_i / dot(r_hat_0, q)
        alpha = rho / dot(r_hat, q, N);

        // s = r - alpha * q;
        scalar_vector_add(s, r, -alpha, q, N);

        // x = x + alpha * p_hat
        scalar_vector_add(x, x, alpha, p_hat, N);

        // apply preconditioner to s (s_hat = inv_M s)
        apply_diagonal_matrix(s_hat, inv_M, s, N);

        // t = A s_hat
        A->apply(s_hat, t);

        // omega_i = dot(t, s) / dot(t, t)
        omega = dot(t, s, N) / dot(t, t, N);

        // x = x + omega * s_hat;
        scalar_vector_add(x, x, omega, s_hat, N);

        // r = s - omega * t;
        scalar_vector_add(r, s, -omega, t, N);

        //		printf("BiCGSTAB_solver: %f\n", dot(s,s,N));
    }
}

void BiCGSTAB_solver::get_residual_vector(std::vector<scalar> &r, std::vector<scalar> &b, std::unique_ptr<SparseMatrixUnknownPattern> &A, std::vector<scalar> &x, const int N) {
    // Ax
    A->apply(x, r);

    // r = b - Ax
    for (int i = 0; i < N; i++)
        r[i] = b[i] - r[i];
}

/* Copies the contents of vector b into vector a (a <- b) */
void BiCGSTAB_solver::copy_vector(std::vector<scalar> &a, const std::vector<scalar> &b, const int N) {
    for (int i = 0; i < N; i++)
        a[i] = b[i];
}

/* Returns the dot product of vectors a and b, of length N */
scalar BiCGSTAB_solver::dot(const std::vector<scalar> &a, const std::vector<scalar> &b, const int N) {
    scalar result = 0;
    for (int i = 0; i < N; i++)
        result += a[i] * b[i];

    return result;
}

/* Calculates y = Mx for diagonal matrix and vectors of dimension N */
inline void BiCGSTAB_solver::apply_diagonal_matrix(std::vector<scalar> &y, const std::vector<scalar> &M, const std::vector<scalar> &x, const int N) {
    for (int i = 0; i < N; i++)
        y[i] = M[i] * x[i];
}

/* Sets the given vector (length N) to zero */
void BiCGSTAB_solver::zero(std::vector<scalar> &x, const int N) {
    for (int i = 0; i < N; i++)
        x[i] = 0;
}

/* Carries out the operation x = y + c*z, where c is a scalar, and x, y and z are vectors of length N. */
void BiCGSTAB_solver::scalar_vector_add(std::vector<scalar> &x, const std::vector<scalar> &y, scalar c, const std::vector<scalar> &z, const int N) {
    for (int i = 0; i < N; i++)
        x[i] = y[i] + c * z[i];
}

/* Carries out the operation w = x + a * (y + b * z) */
// p = r + beta(p - omega * v)

void BiCGSTAB_solver::complicated_machine(std::vector<scalar> &w, const std::vector<scalar> &x, scalar a, const std::vector<scalar> &y, scalar b, const std::vector<scalar> &z, const int N) {
    for (int i = 0; i < N; i++)
        w[i] = x[i] + a * (y[i] + b * z[i]);
}

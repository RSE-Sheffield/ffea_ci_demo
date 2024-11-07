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

#include "CG_solver.h"

#include "mat_vec_fns_II.h"

CG_solver::CG_solver() {
    N = 0;
    tol = 0;
    max_num_iterations = 0;
    inv_M = {};
    d = {};
    r = {};
    q = {};
    s = {};
}

CG_solver::~CG_solver() {
    N = 0;
    tol = 0;
    max_num_iterations = 0;
    inv_M.clear();
    d.clear();
    r.clear();
    q.clear();
    s.clear();
}

void CG_solver::init(int N, scalar tol, int max_num_iterations) {
    this->N = N;
    this->tol = tol;
    this->max_num_iterations = max_num_iterations;

    // Allocate all required memory
    try {
        inv_M = std::vector<scalar>(N, 0);
        d = std::vector<scalar>(N, 0);
        r = std::vector<scalar>(N, 0);
        q = std::vector<scalar>(N, 0);
        s = std::vector<scalar>(N, 0);
    } catch(std::bad_alloc &) {
        throw FFEAException("While initialising CG_solver, could not allocate memory for vectors.\n");
    }
}

void CG_solver::solve(std::shared_ptr<SparseMatrixFixedPattern> &A, std::vector<scalar> &x, const std::vector<scalar> &b) {
    // Get the preconditioner matrix (inverse of the matrix diagonal)
    A->calc_inverse_diagonal(inv_M);

    scalar delta_new, delta_old, alpha;

    // Get the residual vector for Ax = b with current x
    delta_new = conjugate_gradient_residual(A, x, b);

    for (int i = 0; i < max_num_iterations; i++) {

        // Once convergence is achieved, return
        if (residual2() < tol) {
            std::cout << "CG_solver: Convergence reached on iteration " << i << "\n"; // DEBUGGO
            return;
        }

        // q = A * d
        A->apply(d, q);

        alpha = delta_new / dot(d, q);

        parallel_vector_add_self(x, alpha, d);
        parallel_vector_add_self(r, -alpha, q);

        delta_old = delta_new;
        delta_new = parallel_apply_preconditioner();

        parallel_vector_add(d, (delta_new / delta_old), s);

    }
    // If desired convergence was not reached in the set number of iterations...
    throw FFEAException("CG_solver: Could not converge after %d iterations.\n\tEither epsilon or max_iterations_cg are set too low, or something went wrong with the simulation.\n", max_num_iterations);
}

void CG_solver::solve(std::shared_ptr<SparseMatrixFixedPattern> &A, std::vector<scalar> &x, const std::vector<scalar> &b, int num_iterations) {
    // Get the preconditioner matrix (inverse of the matrix diagonal)
    A->calc_inverse_diagonal(inv_M);

    scalar delta_new, delta_old, alpha;

    // Get the residual vector for Ax = b with current x
    delta_new = conjugate_gradient_residual(A, x, b);

    for (int i = 0; i < num_iterations; i++) {
        // q = A * d
        A->apply(d, q);

        alpha = delta_new / dot(d, q);

        parallel_vector_add_self(x, alpha, d);
        parallel_vector_add_self(r, -alpha, q);

        delta_old = delta_new;
        delta_new = parallel_apply_preconditioner();

        parallel_vector_add(d, (delta_new / delta_old), s);
    }
}

scalar CG_solver::conjugate_gradient_residual(std::shared_ptr<SparseMatrixFixedPattern> &A, const std::vector<scalar> &x, const std::vector<scalar> &b) {
    // Ax
    A->apply(x, r);

    // r = b - Ax
    // d = inv_M * r
    // delta_new = r . d
    scalar delta_new = 0;
    for (int i = 0; i < N; i++) {
        r[i] = b[i] - r[i];
        d[i] = inv_M[i] * r[i];
        delta_new += r[i] * d[i];
    }

    return delta_new;
}

scalar CG_solver::residual2() {
    int i;
    scalar r2 = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
//#pragma omp parallel for default(none) private(i) reduction(+:r2)
#endif
    for (i = 0; i < N; i++) {
        r2 += r[i] * r[i];
    }

    return r2;
}

void CG_solver::parallel_vector_add_self(std::vector<scalar> &v1, scalar a, const std::vector<scalar> &v2) {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) shared(v1, a, v2)
#endif
    for (int i = 0; i < N; i++) {
        v1[i] += a * v2[i];
    }
}

void CG_solver::parallel_vector_add(std::vector<scalar> &v1, scalar a, const std::vector<scalar> &v2) {
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) shared(v1, a, v2)
#endif
    for (int i = 0; i < N; i++) {
        v1[i] = v2[i] + a * v1[i];
    }
}

scalar CG_solver::parallel_apply_preconditioner() {
    scalar delta_new = 0;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
//#pragma omp parallel for default(none) reduction(+:delta_new)
#endif
    for (int i = 0; i < N; i++) {
        s[i] = inv_M[i] * r[i];
        delta_new += r[i] * s[i];
    }

    return delta_new;
}

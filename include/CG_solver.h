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

#ifndef CG_SOLVER_H_INCLUDED
#define CG_SOLVER_H_INCLUDED

#include <stdlib.h>
#include <stdio.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "SparseMatrixFixedPattern.h"

class CG_solver {
public:
    CG_solver();
    ~CG_solver();

    void init(int N, scalar tol, int max_num_iterations);

    void solve(std::shared_ptr<SparseMatrixFixedPattern> &A, std::vector<scalar> &x, const std::vector<scalar> &b);

    void solve(std::shared_ptr<SparseMatrixFixedPattern> &A, std::vector<scalar> &x, const std::vector<scalar> &b, int num_iterations);

private:

    int N; ///< The number of unknowns

    scalar tol; ///< The convergence tolerance threshold

    int max_num_iterations; ///< Maximum number of iterations before giving up 

    std::vector<scalar>inv_M; ///< The preconditioner matrix (inverse of the diagonal)

    std::vector<scalar> d; ///< Vector needed for use by conjugate gradient solver
    std::vector<scalar> r; ///< Vector needed for use by conjugate gradient solver
    std::vector<scalar> q; ///< Vector needed for use by conjugate gradient solver
    std::vector<scalar> s; ///< Vector needed for use by conjugate gradient solver

    scalar conjugate_gradient_residual(std::shared_ptr<SparseMatrixFixedPattern> &, const std::vector<scalar> &x, const std::vector<scalar> &b);

    scalar residual2();

    void parallel_vector_add_self(std::vector<scalar> &v1, scalar a, const std::vector<scalar> &v2);

    void parallel_vector_add(std::vector<scalar> &v1, scalar a, const std::vector<scalar> &v2);

    scalar parallel_apply_preconditioner();

    /** Returns the dot product of vectors a and b, of length N */
    scalar dot(const std::vector<scalar> &a, const std::vector<scalar> &b);
};

#endif

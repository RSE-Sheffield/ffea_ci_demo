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

#ifndef BICGSTAB_SOLVER_H_INCLUDED
#define BICGSTAB_SOLVER_H_INCLUDED

#include <memory>
#include <stdlib.h>
#include <stdio.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "SparseMatrixUnknownPattern.h"

class BiCGSTAB_solver {
public:
    BiCGSTAB_solver() = default;

    ~BiCGSTAB_solver() = default;

    void init(int N, scalar _tol, int _max_num_iterations);
    /**
     * @param num_iterations If 0, max_num_iterations will instead be used and the method will return at convergence
     */
    void solve(unique_ptr<SparseMatrixUnknownPattern> &A, std::vector<scalar> &x, std::vector<scalar> &b, int num_iterations = 0);

private:
    /// The convergence tolerance threshold
    scalar tol = 0;

    /// Maximum number of iterations before giving up
    int max_num_iterations = 0;

    /// The inverse of the preconditioner matrix
    std::vector<scalar> inv_M;

    //@{
    /// The residual vectors
    std::vector<scalar> r, r_hat;
    //@}

    //@{
    /// Other necessary vectors
    std::vector<scalar> p, p_hat, q, s, s_hat, t;
    //@}

    /** Calculates the residual vector, r, for the matrix equation Ax = b. Specifically, r = b - Ax. */
    void get_residual_vector(std::vector<scalar> &r, std::vector<scalar> &b, std::unique_ptr<SparseMatrixUnknownPattern> &A, std::vector<scalar> &x);
    
    /** Carries out the operation x = y + c*z, where c is a scalar, and x, y and z are vectors of length N. */
    void scalar_vector_add(std::vector<scalar> &x, const std::vector<scalar> &y, scalar c, const std::vector<scalar> &z);

    /** Carries out the operation w = x + a * (y + b * z) */
    // p = r + beta(p - omega * v)
    void complicated_machine(std::vector<scalar> &w, const std::vector<scalar> &x, scalar a, const std::vector<scalar> &y, scalar b, const std::vector<scalar> &z);
};

#endif

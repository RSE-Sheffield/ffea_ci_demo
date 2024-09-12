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

#ifndef CONJUGATEGRADIENTSOLVER_HPP_INCLUDED
#define CONJUGATEGRADIENTSOLVER_HPP_INCLUDED


#include <stdio.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "ConjugateGradientSolver.h"
#include "SparseMatrixTypes.h"

class ConjugateGradientSolver : public Solver {
public:
    /** Constructor */
    ConjugateGradientSolver();

    /** Destructor */
    ~ConjugateGradientSolver();

    /** Builds the sparse mass matrix and allocates the various work vectors required for conjugate gradient */
    int init(std::vector<mesh_node> &node, std::vector<tetra_element_linear> &elem, const SimulationParams &params, const std::vector<int> &pinned_nodes_list, const set<int> &bsite_pinned_node_list);

    /** Applies conjugate gradient with a Jacobi preconditioner to solve the system Mx = f */
    int solve(std::vector<arr3> &x);

    /** Applies the mass matrix to the given vector, 'in', putting the result in 'result'*/
    void apply_matrix(scalar *in, scalar *result);

private:

    /** Error tolerance threshold (normalised and squared) to determine when solution has converged */
    scalar epsilon2;

    /** Maximum number of iterations the solver should use before giving up (as solution is not converging) */
    int i_max;

    /** Number of rows in original matrix */
    int num_rows;

    /** Array of all non-zero entries comprising the mass matrix, in the order they appear in the matrix
     * when scanned left to right across rows first (and columns secondary)
     */
    sparse_entry *entry;

    /**
     * Array of size (num_nodes+1) containing the index (in the sparse_entry array above)
     * of the first non-zero entry in each row of the original matrix.
     * The final entry in this vector, key[num_nodes], is the total number of non-zero entries
     * in the entire matrix.
     */
    int *key;

    /** Jacobi preconditioner (inverse of the mass matrix diagonal) */
    scalar *preconditioner;

    /** Work vectors */
    std::vector<arr3> d, r, q, s, f;

    /* */
    scalar conjugate_gradient_residual_assume_x_zero(std::vector<arr3> &b);

    /* */
    scalar parallel_sparse_matrix_apply();

    void parallel_vector_add_self(std::vector<arr3> &v1, scalar a, std::vector<arr3> &v2, int vec_size);

    void parallel_vector_add(std::vector<arr3> &v1, scalar a, std::vector<arr3> &v2, int vec_size);

    /* */
    scalar parallel_apply_preconditioner();

    /* */
    scalar residual2();
};

#endif

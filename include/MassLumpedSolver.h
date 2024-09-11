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

#ifndef MASSLUMPEDSOLVER_HPP_INCLUDED
#define MASSLUMPEDSOLVER_HPP_INCLUDED


#include <stdio.h>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mesh_node.h"
#include "tetra_element_linear.h"
#include "SimulationParams.h"
#include "Solver.h"
#include "SparseMatrixTypes.h"

class MassLumpedSolver : public Solver {
public:

    /** Constructor */
    MassLumpedSolver();

    /** Destructor */
    ~MassLumpedSolver();

    /** Builds the diagonal mass matrix and gets reciprocal of each value */
    int init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list, set<int> bsite_pinned_node_list);

    /** Applies inverse mass matrix (since diagonal: Mx = f => x = f_i/M_i */
    int solve(arr3 *x);

    /** Applies the mass matrix to the given vector, 'in', putting the result in 'result'*/
    void apply_matrix(scalar *in, scalar *result);

private:

    /** Number of rows in original matrix */
    int num_rows;

    /** Mass matrix diagonal inversed */
    scalar *inv_M;
};

#endif

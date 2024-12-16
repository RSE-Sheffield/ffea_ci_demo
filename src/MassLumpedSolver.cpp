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

#include "MassLumpedSolver.h"

/* */
void MassLumpedSolver::init(std::vector<mesh_node> &node, std::vector<tetra_element_linear> &elem, const SimulationParams &params, const std::vector<int> &pinned_nodes_list, const set<int> &bsite_pinned_node_list) {
    // Store the number of rows, error threshold (stopping criterion for solver) and max
    // number of iterations, on this Solver (these quantities will be used a lot)
    this->num_rows = node.size();
    try {
        inv_M = std::vector<scalar>(num_rows);
    } catch(std::bad_alloc &) {
        throw FFEAException("could not allocate inv_M.");
    }

    for (int i = 0; i < num_rows; i++) {
        inv_M[i] = 0;
    }
    for (int n = 0; n < elem.size(); n++) {
        // add mass contribution for this element
        for (int i = 0; i < 10; i++) {
            if (i < 4) {
                int ni = elem[n].n[i]->index;
                inv_M[ni] += .25 * elem[n].rho * elem[n].vol_0;
            } else {
                int ni = elem[n].n[i]->index;
                inv_M[ni] = 1;
            }
        }
    }

    // set elements corresponding to unmovable 'pinned' nodes to 1
    for (int i = 0; i < pinned_nodes_list.size(); ++i) {
        inv_M[pinned_nodes_list[i]] = 1.0;
    }

    // inverse
    for (int i = 0; i < num_rows; i++) {
        inv_M[i] = 1.0 / inv_M[i];
    }
}

/* */
void MassLumpedSolver::solve(std::vector<arr3> &x) {
    int i = 0;
    for (i = 0; i < num_rows; i++) {
        x[i][0] = x[i][0] * inv_M[i];
        x[i][1] = x[i][1] * inv_M[i];
        x[i][2] = x[i][2] * inv_M[i];
    }
}

/* */
void MassLumpedSolver::apply_matrix(const std::vector<scalar> &in, std::vector<scalar> &result) {
    for (int i = 0; i < num_rows; i++) {
        result[i] *= 1.0 / inv_M[i];
    }
}


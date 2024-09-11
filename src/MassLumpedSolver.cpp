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

/**/
MassLumpedSolver::MassLumpedSolver() {
    num_rows = 0;
    inv_M = nullptr;
}

/* */
MassLumpedSolver::~MassLumpedSolver() {
    delete[] inv_M;
    inv_M = nullptr;
    num_rows = 0;
}

/* */
int MassLumpedSolver::init(int num_nodes, int num_elements, mesh_node *node, tetra_element_linear *elem, SimulationParams *params, int num_pinned_nodes, int *pinned_nodes_list, set<int> bsite_pinned_node_list) {
    int n, i, ni;

    // Store the number of rows, error threshold (stopping criterion for solver) and max
    // number of iterations, on this Solver (these quantities will be used a lot)
    this->num_rows = num_nodes;
    inv_M = new scalar[num_rows];

    if (!inv_M) {
        FFEA_ERROR_MESSG("could not allocate inv_M\n");
    }

    for (i = 0; i < num_rows; i++) {
        inv_M[i] = 0;
    }
    for (n = 0; n < num_elements; n++) {
        // add mass contribution for this element
        for (i = 0; i < 10; i++) {
            if (i < 4) {
                ni = elem[n].n[i]->index;
                inv_M[ni] += .25 * elem[n].rho * elem[n].vol_0;
            } else {
                ni = elem[n].n[i]->index;
                inv_M[ni] = 1;
            }
        }
    }

    // set elements corresponding to unmovable 'pinned' nodes to 1
    for (i = 0; i < num_pinned_nodes; i++) {
        inv_M[pinned_nodes_list[i]] = 1.0;
    }

    // inverse
    for (i = 0; i < num_rows; i++) {
        inv_M[i] = 1.0 / inv_M[i];
    }

    return FFEA_OK;
}

/* */
int MassLumpedSolver::solve(arr3 *x) {
    int i = 0;
    for (i = 0; i < num_rows; i++) {
        x[i][0] = x[i][0] * inv_M[i];
        x[i][1] = x[i][1] * inv_M[i];
        x[i][2] = x[i][2] * inv_M[i];
    }
    return FFEA_OK;
}

/* */
void MassLumpedSolver::apply_matrix(scalar *in, scalar *result) {
    for (int i = 0; i < num_rows; i++) {
        result[i] *= 1.0 / inv_M[i];
    }
}


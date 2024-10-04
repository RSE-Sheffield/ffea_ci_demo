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

#ifndef SOLVER_H_INCLUDED
#define SOLVER_H_INCLUDED

#include "tetra_element_linear.h"
#include <set>

class Solver {
public:

    /** Provide an empty constructor (to avoid linker problems) */
    Solver() {
    };

    /** Make the destructor virtual (so that the destructor of derived classes will be called too) */
    virtual ~Solver() {
    };

    /**
     * Initialises the solver (by building whatever representation of the mass matrix it needs)
     * using the given node-element connectivity.
     */
    virtual void init(std::vector<mesh_node> &node, std::vector<tetra_element_linear> &elem, const SimulationParams &params, const std::vector<int> &pinned_nodes_list, const set<int> &bsite_pinned_node_list) = 0;

    /**
     * Solves the linear system Mx = f where f is the force vector (should be 'x' on input), and the mass matrix
     * M has already been constructed (in whatever representation) by the init() function for
     * a particular Blob. The solution is written to x.
     */
    virtual void solve(std::vector<arr3> &x) = 0;

    virtual void apply_matrix(const std::vector<scalar> &in, std::vector<scalar> &result) = 0;
};
#endif

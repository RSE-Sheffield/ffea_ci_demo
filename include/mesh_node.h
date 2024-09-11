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

/*
 *	mesh_node.hpp
 *
 */

#ifndef MESH_NODE_H_INCLUDED
#define MESH_NODE_H_INCLUDED

#include <stddef.h>
#include <stdio.h>
#include "mat_vec_fns.h"

/**
 * Structure for a mesh_node: the points FEM meshes are built from.
 */
class mesh_node {
public:

    mesh_node();
    ~mesh_node();

    void move(int direction, scalar dx);

    void set_pos(scalar x, scalar y, scalar z);

    void print();

    /** Position of node */
    arr3 pos;

    /** Velocity of node */
    arr3 vel;

    /** Electrostatic potential at this node */
    scalar phi;

    int num_element_contributors;

    /** An array of pointers to contributions to the total force on this node. There should be one
     * contribution from each element this node is a part of (so the length will be num_element_contributors).
     */
    arr3 **force_contributions;

    /** Required for some general matrix constructions in which we need to know this node's 'index' in the node vector */
    int index;

    /** Equilibrium position of nodes (for RMSD calculations) */
    arr3 pos_0;

    /** Charge density on this node */
    scalar rho;

    /** Stokes radius of this node */
    scalar stokes_radius;

    /** The drag due to stokes on this node, not including velocity */
    scalar stokes_drag;

    /** Stores whether or not this node is linear (as the order is surface - interior, not linear - secondary) */
    bool linear;

    void set_linear();
    bool am_I_linear();
};

#endif

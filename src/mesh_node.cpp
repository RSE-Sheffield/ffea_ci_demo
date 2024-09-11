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
 *	mesh_node.cpp
 *
 */

#include "mesh_node.h"

#include <cstring>

/*
 * Structure for a mesh_node: the points FEM meshes are built from.
 */
mesh_node::mesh_node() {
    num_element_contributors = 0;
    force_contributions = nullptr;
    memset(&pos, 0, sizeof(arr3));
    memset(&vel, 0, sizeof(arr3));
    phi = 0;
    index = 0;
    memset(&pos_0, 0, sizeof(arr3));
    stokes_radius = 0;
    stokes_drag = 0;
    linear = false;
}

mesh_node::~mesh_node() {
    delete[] force_contributions;
    memset(&pos, 0, sizeof(arr3));
    memset(&vel, 0, sizeof(arr3));
    phi = 0;
    index = 0;
    memset(&pos_0, 0, sizeof(arr3));
    num_element_contributors = 0;
    stokes_radius = 0;
    stokes_drag = 0;
    linear = false;
}

void mesh_node::move(int direction, scalar dx) {
	switch(direction) {
		case(0):
			pos[0] += dx;
			break;
		case(1):
			pos[1] += dx;
			break;
		case(2):
			pos[2] += dx;
			break;
	}
}

void mesh_node::set_pos(scalar x, scalar y, scalar z) {
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
}

void mesh_node::print() {
    printf("pos: %e %e %e\n", pos[0], pos[1], pos[2]);
    printf("vel: %e %e %e\n", vel[0], vel[1], vel[2]);
}

void mesh_node::set_linear() {
    linear = true;
}

bool mesh_node::am_I_linear() {
    return linear;
}

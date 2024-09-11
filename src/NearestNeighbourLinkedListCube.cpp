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

#include "NearestNeighbourLinkedListCube.h"

/* */
int NearestNeighbourLinkedListCube::build_nearest_neighbour_lookup(scalar h) {
    // Clear the grid
    clear();

    // Loop through each Face in the pool, calculate which cell of the
    // grid it belongs in based on its centroid position, and add it to
    // the linked list stack on that cell
    int i;
    for (i = 0; i < num_nodes_in_pool; i++) {

        // calculate which cell the face belongs in
	if (!pool[i].obj->kinetically_active) {
		continue;
	}

	// Do we have the correct centroid?? We do now!
        //pool[i].obj->calc_area_normal_centroid();
        int x = (int) floor(pool[i].obj->centroid[0] / h);
        int y = (int) floor(pool[i].obj->centroid[1] / h);
        int z = (int) floor(pool[i].obj->centroid[2] / h);
        /*
                                        // If face centroid is out of bounds of the box, add its node to the nearest cell on the edge
                                        // of the box. This is necessary in the case of PBC since an entire blob is moved only if its centre
                                        // of mass exceeds the bounds of the box
                                        if(x < 0) {
                                                x = 0;
                                        }
                                        if(x > N - 1) {
                                                x -= N - 1;
                                        }
                                        if(y < 0) {
                                                y = 0;
                                        }
                                        if(y > N - 1) {
                                                y = N - 1;
                                        }
                                        if(z < 0) {
                                                z = 0;
                                        }
                                        if(z > N - 1) {
                                                z = N - 1;
                                        }
         */
        // attempt to add the node to the cell
        if (add_node_to_stack(i, x, y, z) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when trying to add node %d to nearest neighbour stack at (%d %d %d)\n", i, x, y, z);
        }
    }
    return FFEA_OK;
}


int NearestNeighbourLinkedListCube::prebuild_nearest_neighbour_lookup_and_swap(scalar h) {
    // Clear the grid
    clear_layer(shadow_layer);

    // Loop through each Face in the pool, calculate which cell of the
    // grid it belongs in based on its centroid position, and add it to
    // the linked list stack on that cell
    int i;
    for (i = 0; i < num_nodes_in_pool; i++) {

        // calculate which cell the face belongs in
        if (!pool[i].obj->kinetically_active) {
          continue;
        }

       // Do we have the correct centroid?? We do now!
       // pool[i].obj->calc_area_normal_centroid();

        int x = (int) floor(pool[i].obj->centroid[0] / h);
        int y = (int) floor(pool[i].obj->centroid[1] / h);
        int z = (int) floor(pool[i].obj->centroid[2] / h);

        // attempt to add the node to the cell
        if (add_node_to_stack_shadow(i, x, y, z) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when trying to add node %d to nearest neighbour stack at (%d %d %d)\n", i, x, y, z);
        }

    }
    swap_layers(); 
    return FFEA_OK;
}

int NearestNeighbourLinkedListCube::prebuild_nearest_neighbour_lookup(scalar h) {
    // Clear the grid
    clear_layer(shadow_layer);
    can_swap = false;

    // Loop through each Face in the pool, calculate which cell of the
    // grid it belongs in based on its centroid position, and add it to
    // the linked list stack on that cell
    int i;
    for (i = 0; i < num_nodes_in_pool; i++) {

        // calculate which cell the face belongs in
        if (!pool[i].obj->kinetically_active) {
          continue;
        }

       // Do we have the correct centroid?? We do now!
       // pool[i].obj->calc_area_normal_centroid();

        int x = (int) floor(pool[i].obj->centroid[0] / h);
        int y = (int) floor(pool[i].obj->centroid[1] / h);
        int z = (int) floor(pool[i].obj->centroid[2] / h);

        // attempt to add the node to the cell
        if (add_node_to_stack_shadow(i, x, y, z) == FFEA_ERROR) {
            FFEA_ERROR_MESSG("Error when trying to add node %d to nearest neighbour stack at (%d %d %d)\n", i, x, y, z)
        }

    }
    can_swap = true;
    return FFEA_OK;
}

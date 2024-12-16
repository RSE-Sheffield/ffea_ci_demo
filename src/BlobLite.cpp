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

#include "BlobLite.h"

BlobLite::BlobLite() {
  // initialise everything to zero 
  num_nodes = 0;
  num_surface_nodes = 0;
  num_interior_nodes = 0;
  num_elements = 0;
  num_surface_elements = 0;
  num_interior_elements = 0; 
  blob_state = FFEA_BLOB_IS_DYNAMIC; 
} 

BlobLite::~BlobLite() { 
  delete[] coord; 
  delete[] elem; 
} 

void BlobLite::load_topology(const char *topology_filename){
    FILE *in = NULL;
    int i, n1, n2, n3, n4, n5, n6, n7, n8, n9, n10;
    const int max_line_size = 50;
    char line[max_line_size];

    // Now open the topology file
    if ((in = fopen(topology_filename, "r")) == NULL) {
        throw FFEAFileException(topology_filename);
    }
    printf("\t\tReading in topology file: %s\n", topology_filename);

    // first line should be the file type "ffea topology file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        throw FFEAException("Error reading first line of topology file.");
    }
    if (strcmp(line, "walrus topology file\n") != 0 && strcmp(line, "ffea topology file\n") != 0) {
        fclose(in);
        throw FFEAException("This is not a 'ffea topology file' (read '%s').", line);
    }

    // read in the total number of elements in the file
    if (fscanf(in, "num_elements %d\n", &num_elements) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of elements.");
    }
    printf("\t\t\tNumber of elements = %d\n", num_elements);

    // read in the number of surface elements in the file
    if (fscanf(in, "num_surface_elements %d\n", &num_surface_elements) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of surface elements.");
    }
    printf("\t\t\tNumber of surface elements = %d\n", num_surface_elements);

    // read in the number of interior elements in the file
    if (fscanf(in, "num_interior_elements %d\n", &num_interior_elements) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of interior elements.");
    }

    printf("\t\t\tNumber of interior elements = %d\n", num_interior_elements);

    // Allocate the memory for all these elements
    elem = new(std::nothrow) int[NUM_NODES_QUADRATIC_TET*num_elements];
    if (elem == NULL) {
        fclose(in);
        throw FFEAException("Unable to allocate memory for element array.");
    }

    // Check for "surface elements:" line
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        throw FFEAException("Error when looking for 'surface elements:' line.");
    }
    if (strcmp(line, "surface elements:\n") != 0) {
        fclose(in);
        throw FFEAException("Could not find 'surface elements:' line (found '%s' instead).", line);
    }

    // Read in all the elements from file
    for (i = 0; i < num_elements; i++) {
        if (i == num_surface_elements) {
            // Check for "interior elements:" line
            if (fgets(line, max_line_size, in) == NULL) {
                fclose(in);
                throw FFEAException("Error when looking for 'interior elements:' line.");
            }
            if (strcmp(line, "interior elements:\n") != 0) {
                fclose(in);
                throw FFEAException("Could not find 'interior elements:' line (found '%s' instead).", line);
            }
        }

        if (fscanf(in, "%d %d %d %d %d %d %d %d %d %d\n", &n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8, &n9, &n10) != 10) {
            fclose(in);
            throw FFEAException("Error reading from elements file at element %d.", i);
        } else {
            // check that none of these reference nodes outside of the node array
            if (n1 < 0 || n1 >= num_nodes ||
                    n2 < 0 || n2 >= num_nodes ||
                    n3 < 0 || n3 >= num_nodes ||
                    n4 < 0 || n4 >= num_nodes ||
                    n5 < 0 || n5 >= num_nodes ||
                    n6 < 0 || n6 >= num_nodes ||
                    n7 < 0 || n7 >= num_nodes ||
                    n8 < 0 || n8 >= num_nodes ||
                    n9 < 0 || n9 >= num_nodes ||
                    n10 < 0 || n10 >= num_nodes) {
                fclose(in);
                throw FFEAException("Error: Element %d references an out of bounds node index.", i);
            }

       // Link element nodes to actual nodes
            store_index_to_elemnode(n1,i,0); 
            store_index_to_elemnode(n2,i,1); 
            store_index_to_elemnode(n3,i,2); 
            store_index_to_elemnode(n4,i,3); 
            store_index_to_elemnode(n5,i,4); 
            store_index_to_elemnode(n6,i,5); 
            store_index_to_elemnode(n7,i,6); 
            store_index_to_elemnode(n8,i,7); 
            store_index_to_elemnode(n9,i,8); 
            store_index_to_elemnode(n10,i,9); 

        }
    }

    fclose(in);

    if (i == 1)
        printf("\t\t\tRead 1 element from %s\n", topology_filename);
    else
        printf("\t\t\tRead %d elements from %s\n", i, topology_filename);
}

   

void BlobLite::store_index_to_elemnode(int node_index, int elem_number, int node_number){ 
    elem[elem_number*NUM_NODES_QUADRATIC_TET + node_number] = node_index;
}


void BlobLite::load_nodes(const char *node_filename, scalar scale) {
    FILE *in = NULL;
    int i;
    double x, y, z;
    const int max_line_size = 50;
    char line[max_line_size];

    // open the node file
    if ((in = fopen(node_filename, "r")) == NULL) {
        throw FFEAFileException(node_filename);
    }
    printf("\t\tReading in nodes file: %s\n", node_filename);

    // first line should be the file type "ffea node file"
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        throw FFEAException("Error reading first line of node file.");
    }
    if (strcmp(line, "walrus node file\n") != 0 && strcmp(line, "ffea node file\n") != 0) {
        fclose(in);
        throw FFEAException("This is not a 'ffea node file' (read '%s').", line);
    }

    // read in the number of nodes in the file
    if (fscanf(in, "num_nodes %d\n", &num_nodes) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of nodes.");
    }
    printf("\t\t\tNumber of nodes = %d\n", num_nodes);

    // read in the number of surface nodes in the file
    if (fscanf(in, "num_surface_nodes %d\n", &num_surface_nodes) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of surface nodes.");
    }
    printf("\t\t\tNumber of surface nodes = %d\n", num_surface_nodes);

    // read in the number of interior nodes in the file
    if (fscanf(in, "num_interior_nodes %d\n", &num_interior_nodes) != 1) {
        fclose(in);
        throw FFEAException("Error reading number of interior nodes.");
    }
    printf("\t\t\tNumber of interior nodes = %d\n", num_interior_nodes);

    // Allocate the memory for all these nodes
    coord = new(std::nothrow) scalar[3*num_nodes];
    if (coord == NULL) throw FFEAException("Failed to allocate node coordinates.");;

    // Check for "surface nodes:" line
    if (fgets(line, max_line_size, in) == NULL) {
        fclose(in);
        throw FFEAException("Error when looking for 'surface nodes:' line.");
    }
    if (strcmp(line, "surface nodes:\n") != 0) {
        fclose(in);
        throw FFEAException("Could not find 'surface nodes:' line (found '%s' instead).", line);
    }

    // Read in all the surface nodes from file
    for (i = 0; i < num_nodes; i++) {

        if (i == num_surface_nodes) {
            // Check for "interior nodes:" line
            if (fgets(line, max_line_size, in) == NULL) {
                fclose(in);
                throw FFEAException("Error when looking for 'interior nodes:' line.");
            }
            if (strcmp(line, "interior nodes:\n") != 0) {
                fclose(in);
                throw FFEAException("Could not find 'interior nodes:' line (found '%s' instead).", line);
            }
        }

        if (fscanf(in, "%le %le %le\n", &x, &y, &z) != 3) {
            fclose(in);
            throw FFEAException("Error reading from nodes file at node %d.", i);
        } else {
            coord[3*i] = scale * (scalar) x;
            coord[3*i+1] = scale * (scalar) y;
            coord[3*i+2] = scale * (scalar) z;
        }
    }

    fclose(in);
    printf("\t\t\tRead %d nodes from %s\n", i, node_filename);
}



void BlobLite::read_nodes_from_file(FILE *trj) {

    const int line_size = 256;
    char state_str[line_size];
    char *result = NULL;

    // If blob is static, don't read any nodes. Simply read the word "STATIC"
    if (blob_state == FFEA_BLOB_IS_STATIC) {
        result = fgets(state_str, line_size, trj);
        if (result == NULL) {
            throw FFEAException("Problem when reading 'STATIC' (expected) line in trajectory file.");
        }
        if (strcmp(state_str, "STATIC\n") != 0) {
            throw FFEAException("When restarting from trajectory file, expected to read 'STATIC', but instead found '%s...'.", state_str);
        }
        return;
    } else {
        result = fgets(state_str, line_size, trj);
        if (result == NULL) {
            throw FFEAException("Problem when reading state line in trajectory file.");
        }
    }

    scalar x, y, z, u;
    for (int i = 0; i < num_nodes; i++) {
        if (fscanf(trj, "%le %le %le %le %le %le %le %le %le %le\n", &x, &y, &z, &u, &u, &u, &u, &u, &u, &u) != 10) {
            throw FFEAException("(When restarting) Error reading from trajectory file, for node %d.", i);
        } else {
          coord[3*i] = (scalar) x / mesoDimensions::length;
          coord[3*i+1] = (scalar) y / mesoDimensions::length;
          coord[3*i+2] = (scalar) z / mesoDimensions::length;
          /* node[i].pos.x /= mesoDimensions::length;
          node[i].pos.y /= mesoDimensions::length;
          node[i].pos.z /= mesoDimensions::length;
          node[i].vel.x /= mesoDimensions::velocity;
          node[i].vel.y /= mesoDimensions::velocity;
          node[i].vel.z /= mesoDimensions::velocity;
          force[i].x /= mesoDimensions::force;
          force[i].y /= mesoDimensions::force;
          force[i].z /= mesoDimensions::force;
          // mesoDimensions... */ 
        } 

    }
}

int BlobLite::icoord_for_elem_node(int elemi, int nodei){
   return 3*elem[elemi * NUM_NODES_QUADRATIC_TET + nodei ];
}

void BlobLite::center_of_coord(arr3 &cm){
   cm[0] = 0;
   cm[1] = 0;
   cm[2] = 0;
   for (int i=0; i<num_nodes; i++){
     cm[0] = coord[ 3*i ];  
     cm[1] = coord[ 3*i + 1];  
     cm[2] = coord[ 3*i + 2];  
   }
   cm[0] /= num_nodes;
   cm[1] /= num_nodes;
   cm[2] /= num_nodes;
}

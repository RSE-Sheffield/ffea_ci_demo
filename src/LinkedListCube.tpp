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

#include "LinkedListCube.h"

/*  */
template <class T>
LinkedListCube<T>::LinkedListCube() {
    N_x = 0;
    N_y = 0;
    N_z = 0;
    max_num_nodes_in_pool = 0;
    num_nodes_in_pool = 0;
    add_index = 0;
    root = NULL;
    pool = NULL;
    can_swap = true; 
    active_layer = 1;
    shadow_layer = 2;
}

/*  */
template <class T>
LinkedListCube<T>::~LinkedListCube() {
    delete[] pool;
    delete[] root;

    N_x = 0;
    N_y = 0;
    N_z = 0;
    max_num_nodes_in_pool = 0;
    num_nodes_in_pool = 0;
    add_index = 0;
    root = NULL;
    pool = NULL;
}

/* */
template <class T>
int LinkedListCube<T>::alloc(int N_x, int N_y, int N_z, int max_num_nodes_in_pool) {
    this->N_x = N_x;
    this->N_y = N_y;
    this->N_z = N_z;
    this->max_num_nodes_in_pool = max_num_nodes_in_pool;
    num_nodes_in_pool = 0;
    root = new(std::nothrow) LinkedListNode<T> * [N_x * N_y * N_z];
    pool = new(std::nothrow) LinkedListNode<T>[max_num_nodes_in_pool];

    if (root == NULL || pool == NULL) {
        FFEA_ERROR_MESSG("Could not allocate memory (for root and pool arrays) in LinkedListCube\n");
    }

    // Make sure all pointers are initialised to NULL
    clear();

    // Set the current addition index to the beginning of the pool array
    add_index = 0;

    return FFEA_OK;
}

/* */
template <class T>
int LinkedListCube<T>::alloc_dual(int N_x, int N_y, int N_z, int max_num_nodes_in_pool) {
    this->N_x = N_x;
    this->N_y = N_y;
    this->N_z = N_z;
    this->max_num_nodes_in_pool = max_num_nodes_in_pool;
    num_nodes_in_pool = 0;
    root1 = new(std::nothrow) LinkedListNode<T> * [N_x * N_y * N_z];
    pool1 = new(std::nothrow) LinkedListNode<T>[max_num_nodes_in_pool];
    root2 = new(std::nothrow) LinkedListNode<T> * [N_x * N_y * N_z];
    pool2 = new(std::nothrow) LinkedListNode<T>[max_num_nodes_in_pool];

    if (root1 == NULL || pool1 == NULL || root2 == NULL || pool2 == NULL) {
        FFEA_ERROR_MESSG("Could not allocate memory (for root and pool arrays) in LinkedListCube\n");
    }

    // Make sure all pointers are initialised to NULL
    pool = pool1;
    root = root1;
    clear();
    pool_shadow = pool2;
    root_shadow = root2;
    clear_layer(2); 

    // Set the current addition index to the beginning of the pool array
    add_index = 0;

    return FFEA_OK;
}

/* */
template <class T>
int LinkedListCube<T>::add_to_pool(T *t) {
    if (add_index >= max_num_nodes_in_pool) {
        FFEA_ERROR_MESSG("In LinkedListCube, attempt to add more nodes to the pool than space has been allocated for.\n");
    }

    // Give object its own representant LinkedListNode in the pool
    pool[add_index].obj = t;
    pool[add_index].index = add_index;
    add_index++;

    num_nodes_in_pool++;

    return FFEA_OK;
}


/* */
template <class T>
int LinkedListCube<T>::add_to_pool_dual(T *t) {
    if (add_index >= max_num_nodes_in_pool) {
        FFEA_ERROR_MESSG("In LinkedListCube, attempt to add more nodes to the pool than space has been allocated for.\n");
    }

    // Give object its own representant LinkedListNode in the pool
    pool1[add_index].obj = t;
    pool1[add_index].index = add_index;
    pool2[add_index].obj = t;
    pool2[add_index].index = add_index;
    add_index++;

    num_nodes_in_pool++;

    return FFEA_OK;
}

/* */
template <class T>
LinkedListNode<T> * LinkedListCube<T>::get_from_pool(int i) {
    return &pool[i];
}

/* */
template <class T>
void LinkedListCube<T>::clear() {
    // Clear the grid
    #pragma omp simd
    for (int i = 0; i < N_x * N_y * N_z; i++)
        root[i] = NULL;

    // Clear the pool
    for (int i = 0; i < max_num_nodes_in_pool; i++) {
        pool[i].next = NULL;
	//pool[i].obj = NULL;
    }
}

/* */
template <class T>
void LinkedListCube<T>::clear_layer(int l) {
	if (l == 1) {
		// Clear the grid
		for (int i = 0; i < N_x * N_y * N_z; i++) {
			root1[i] = NULL;
		}
		// Clear the pool
		for (int i = 0; i < max_num_nodes_in_pool; i++) {
			pool1[i].next = NULL;
		//	pool1[i].obj = NULL;
		}
	} else if (l == 2) {
		// Clear the grid
		for (int i = 0; i < N_x * N_y * N_z; i++) {
			root2[i] = NULL;
		}
		// Clear the pool
		for (int i = 0; i < max_num_nodes_in_pool; i++) {
			pool2[i].next = NULL;
		//	pool2[i].obj = NULL;
		}
	}
}


/* */
template <class T>
void LinkedListCube<T>::clear_shadow_layer() {
    clear_layer(shadow_layer);
}

/* */
template <class T>
int LinkedListCube<T>::add_node_to_stack(int i, int x, int y, int z) {
    // Apply PBC to the face centroid
    pbc(&x, &y, &z);

    if (x < 0 || x >= N_x || y < 0 || y >= N_y || z < 0 || z >= N_z) {
        printf("Face centroid out of bounds of LinkedListCube (coords are [%d, %d, %d])\n", x, y, z);
        return FFEA_ERROR;
    }

    int abs_index = x * N_y * N_z + y * N_z + z;

    pool[i].x = x;
    pool[i].y = y;
    pool[i].z = z;
    pool[i].next = root[abs_index];
    root[abs_index] = &pool[i];
    return FFEA_OK;
}

/* */
template <class T>
int LinkedListCube<T>::add_node_to_stack_shadow(int i, int x, int y, int z) {
    // Apply PBC to the face centroid
    pbc(&x, &y, &z);

    if (x < 0 || x >= N_x || y < 0 || y >= N_y || z < 0 || z >= N_z) {
        printf("Face centroid out of bounds of LinkedListCube (coords are [%d, %d, %d])\n", x, y, z);
        return FFEA_ERROR;
    }

    int abs_index = x * N_y * N_z + y * N_z + z;

    pool_shadow[i].x = x;
    pool_shadow[i].y = y;
    pool_shadow[i].z = z;
    pool_shadow[i].next = root_shadow[abs_index];
    root_shadow[abs_index] = &pool_shadow[i];
    return FFEA_OK;
}

/* */
template <class T>
LinkedListNode<T> * LinkedListCube<T>::get_top_of_stack(int x, int y, int z) {

    pbc(&x, &y, &z);
    if (x < 0 || x >= N_x || y < 0 || y >= N_y || z < 0 || z >= N_z) {
        printf("Error: Looking for stack in out of bounds cell %d %d %d\n", x, y, z);
        return NULL;
    }

    return root[x * N_y * N_z + y * N_z + z];
}

/* */
template <class T>
int LinkedListCube<T>::get_pool_size() {
    return num_nodes_in_pool;
}

template <class T>
void LinkedListCube<T>::get_dim(int *Nx, int *Ny, int *Nz) {
    *Nx = N_x;
    *Ny = N_y;
    *Nz = N_z;
}

template <class T>
void LinkedListCube<T>::pbc(int *x, int *y, int *z) {
    if (*x < 0) {
        *x += N_x;
    } else if (*x >= N_x) {
        *x -= N_x;
    }
    if (*y < 0) {
        *y += N_y;
    } else if (*y >= N_y) {
        *y -= N_y;
    }
    if (*z < 0) {
        *z += N_z;
    } else if (*z >= N_z) {
        *z -= N_z;
    }
}

template <class T>
void LinkedListCube<T>::swap_layers() {

    if (active_layer == 1){
      active_layer = 2;
      shadow_layer = 1;
      pool = pool2;
      root = root2; 
      pool_shadow = pool1; 
      root_shadow = root1;
    } else {
      active_layer = 1;
      shadow_layer = 2;
      pool = pool1;
      root = root1; 
      pool_shadow = pool2; 
      root_shadow = root2;
    } 

}

template <class T>
int LinkedListCube<T>::safely_swap_layers() {

    if (can_swap == true) {
       swap_layers(); 
       return FFEA_OK; 
    } else return FFEA_ERROR;

}

template <class T>
void LinkedListCube<T>::allow_swapping() {

    can_swap = true;

}

template <class T>
void LinkedListCube<T>::forbid_swapping() {

    can_swap = false;

}

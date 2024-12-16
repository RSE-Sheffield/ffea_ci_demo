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

#ifndef LINKEDLISTCUBE_H_INCLUDED
#define LINKEDLISTCUBE_H_INCLUDED

template <class T>
struct LinkedListNode {
    /** Pointer to the object this LinkedListNode represents */
    T *obj;

    /** Pointer to the next LinkedListNode in the list */
    LinkedListNode *next;

    /** Remember this node's index in the pool (for matrix construction purposes) */
    int index;

    //@{
    /** Remember which cell this node has been placed in to avoid recalculations */
    int x, y, z;
    //@} 
};

/**
 * @note The implementation of this class could do with overhauling and modernisation,
 *       however without a specific test suite it's hard to guarantee consistent behaviour.
 */
template <class T>
class LinkedListCube {
public:

    /** Constructor */
    LinkedListCube();

    /** Destructor */
    ~LinkedListCube();

    /** Builds a LinkedListCube of dimensions N_x x N_y x N_z, and an array of LinkedListNodes of size max_num_nodes_in_pool */
    void alloc(int N_x, int N_y, int N_z, int max_num_nodes_in_pool);

    /** Builds 2 layers of LinkedListCubes of dimensions N_x x N_y x N_z, a
      *   and an dual array of LinkedListNodes of size max_num_nodes_in_pool */ 
    void alloc_dual(int N_x, int N_y, int N_z, int max_num_nodes_in_pool);

    /** Adds the specified T object to the pool (at the index given by add_index) */
    void add_to_pool(T *t);

    /** Adds the specified T object to the pool (at the index given by add_index) */
    void add_to_pool_dual(T *t);

    /** Returns pointer to ith object in the pool */
    LinkedListNode<T> * get_from_pool(int i);

    /**
     * Completely clears entire grid of all linked lists by setting all pointers to nullptr in root,
     * and all 'next' pointers to nullptr in the pool of LinkedListNodes.
     */
    void clear();
    void clear_layer(int l);
    void clear_shadow_layer();

    /**
     * Adds the specified pool node to the stack at index (x, y, z) on the grid.
     * Does this by setting the current top of stack node to be the 'next' in the list
     * after node i, with node i becoming the new top of stack node on that cell;
     */
    void add_node_to_stack(int i, int x, int y, int z);
    void add_node_to_stack_shadow(int i, int x, int y, int z);

    /** Returns whatever node is at the top of the (linked list) stack in grid cell (x, y, z) */
    LinkedListNode<T> * get_top_of_stack(int x, int y, int z);

    /** Returns how many objects are in the 'pool' */
    int get_pool_size();

    void get_dim(int *Nx, int *Ny, int *Nz);

    void safely_swap_layers(); 

    void allow_swapping(); 
    void forbid_swapping(); 
    

protected:

    /** The cube has dimensions N_x x N_y x N_z */
    int N_x, N_y, N_z;

    /** The number of linked list nodes available for stacking */
    int max_num_nodes_in_pool;

    /** The number of nodes in use */
    int num_nodes_in_pool;

    /** A cubic grid of pointers */
    LinkedListNode<T> **root1;
    LinkedListNode<T> **root2;
    LinkedListNode<T> **root;
    LinkedListNode<T> **root_shadow;

    /** A "pool" of LinkedListNodes available for stacking on the various stacks in the cube */
    LinkedListNode<T> *pool1;
    LinkedListNode<T> *pool2;
    LinkedListNode<T> *pool;
    LinkedListNode<T> *pool_shadow;

    int active_layer;
    int shadow_layer;

    void swap_layers(); 
    bool can_swap;

private:
    /** The current index in the pool array at which nodes are to be added */
    int add_index;

    void pbc(int *x, int *y, int *z);

};

#include "../src/LinkedListCube.tpp"

#endif

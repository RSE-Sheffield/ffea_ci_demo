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

#ifndef SPARSEMATRIXFIXEDPATTERN_H_INCLUDED
#define SPARSEMATRIXFIXEDPATTERN_H_INCLUDED

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mat_vec_fns.h"
#include "SparseMatrixTypes.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

class SparseMatrixFixedPattern {
public:
    SparseMatrixFixedPattern();

    ~SparseMatrixFixedPattern();

    int init(int num_rows, int num_nonzero_elements, sparse_entry *entry, int *key, sparse_entry_sources *source_list);
    int init(int num_rows, int num_entries, scalar *entries, int *key, int *col_indices);

    /** Reconstruct the matrix by adding up all the contributions from the sources stored in the source list */
    void build();

    /** Applies this matrix to the given vector 'in', writing the result to 'result' */
    void apply(scalar *in, scalar *result);

    /** Applies this matrix to the given vector 'in', writing the result to 'result'. 'in' is made of 'arr3's */
    void apply(arr3 *in, arr3 *result);

    /** Applies each matrix element to vector i.e. each arr3 acts as a scalar */
    void block_apply(arr3 *in, arr3 *result);
    
    /** Applies this matrix to the give vector 'in', writing the result to 'result'. 'in' is made of '*arr3's */
    void block_apply(arr3 **in, arr3 **result);

    /** Applies matrix to vector in and leaves result in in */
    void block_apply(arr3 **in);

     /** Applies this matrix to the given sparse matrix 'in', and returns a new sparse matrix */
    SparseMatrixFixedPattern * apply(SparseMatrixFixedPattern *in);

    void calc_inverse_diagonal(scalar *inv_D);

    void print();

    void print_dense();

    /** Prints dense matrix out to file for analysis. I suggest only letting this function run once (step = 1?) */
    void print_dense_to_file(arr3 *a);

    void print_row_column();

    void check_symmetry();

    void am_i_diagonally_dominant();

    sparse_entry * get_entries();

    int * get_key();

    int get_num_nonzero_elements();

    int get_num_rows();

    int get_num_columns();

private:

    /** Number of rows in matrix */
    int num_rows;

    /** Number of nonzero elements in matrix */
    int num_nonzero_elements;

    /** The offdiagonal matrix element values and column positions */
    sparse_entry *entry;

    /** The key (array of integers, indices to the start of each row in entry array) */
    int *key;

    /** An array of pointers to the diagonal elements of the matrix (for fast calculation of inverse diagonal) */
    scalar **diagonal;

    /** Lists the source of contributions to each corresponding entry in the entry array */
    sparse_entry_sources *source_list;
};

#endif

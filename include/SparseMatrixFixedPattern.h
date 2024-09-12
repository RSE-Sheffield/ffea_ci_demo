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

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <memory>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mat_vec_fns.h"
#include "SparseMatrixTypes.h"

using namespace std;

class SparseMatrixFixedPattern {
public:
    SparseMatrixFixedPattern();

    ~SparseMatrixFixedPattern();

    /**
     * \brief 
     * \param num_rows 
     * \param num_nonzero_elements 
     * \param entry Pass me using std::move()
     * \param key Pass me using std::move()
     * \param source_list 
     * \return 
     */
    int init(int num_rows, int num_nonzero_elements, std::vector<sparse_entry> &&entry, std::vector<int> &&key, std::vector<sparse_entry_sources> source_list);
    int init(int num_rows, const std::vector<scalar> &entries, std::vector<int> &&key, const std::vector<int>& col_indices);

    /** Reconstruct the matrix by adding up all the contributions from the sources stored in the source list */
    void build();

    /** Applies this matrix to the given vector 'in', writing the result to 'result' */
    void apply(const std::vector<scalar> &in, std::vector<scalar> &result) const;

    /** Applies this matrix to the given vector 'in', writing the result to 'result'. 'in' is made of 'arr3's */
    void apply(const std::vector<arr3> &in, std::vector<arr3> &result) const;

    /** Applies each matrix element to vector i.e. each arr3 acts as a scalar */
    void block_apply(const std::vector<arr3> &in, std::vector<arr3> &result) const;
    
    /** Applies this matrix to the give vector 'in', writing the result to 'result'. 'in' is made of '*arr3's */
    void block_apply(const std::vector<arr3*> &in, std::vector<arr3*> &result) const;

    /** Applies matrix to vector in and leaves result in in */
    void block_apply(vector<arr3*> &in) const;

     /** Applies this matrix to the given sparse matrix 'in', and returns a new sparse matrix */
    std::shared_ptr<SparseMatrixFixedPattern> apply(std::shared_ptr<SparseMatrixFixedPattern> &in) const;

    void calc_inverse_diagonal(std::vector<scalar> &inv_D) const;

    void print() const;

    void print_dense() const;

    /** Prints dense matrix out to file for analysis. I suggest only letting this function run once (step = 1?) */
    void print_dense_to_file(std::vector<arr3> &a) const;

    void print_row_column() const;

    void check_symmetry() const;

    void am_i_diagonally_dominant() const;

    std::vector<sparse_entry> &get_entries();

    std::vector<int> &get_key();

    int get_num_nonzero_elements() const;

    int get_num_rows() const;

    int get_num_columns() const;

private:

    /** Number of rows in matrix */
    int num_rows;

    /** Number of nonzero elements in matrix */
    int num_nonzero_elements;

    /** The offdiagonal matrix element values and column positions */
    std::vector<sparse_entry> entry;

    /** The key (array of integers, indices to the start of each row in entry array) */
    std::vector<int> key;

    /** An array of pointers to the diagonal elements of the matrix (for fast calculation of inverse diagonal) */
    std::vector<scalar *>diagonal;

    /** Lists the source of contributions to each corresponding entry in the entry array */
    std::vector<sparse_entry_sources> source_list;
};

#endif

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

#include "SparseMatrixUnknownPattern.h"

void SparseMatrixUnknownPattern::init(int _num_rows, int suggested_initial_size_for_row_vectors) {

    num_rows = _num_rows;

    // Allocate the array of row vectors
    try {
        row = vector<vector<sparse_entry>>(num_rows);
    } catch(std::bad_alloc &) {
        throw FFEAException("Failed to allocate SparseMatrixUnknownPattern::row vector (len %d).", num_rows);
    }

    // Reserve a certain amount of space (for the inevitable large amount of growth that will come)
    for (int i = 0; i < num_rows; i++) {
        row[i].reserve(suggested_initial_size_for_row_vectors);
    }

    // Allocate the array for the diagonal elements (which are treated differently)
    try {
        diagonal = vector<scalar>(num_rows);
    } catch(std::bad_alloc &) {
        throw FFEAException("Failed to allocate SparseMatrixUnknownPattern::diagonal vector (len %d).", num_rows);
    }

    // Initialise diagonal to zero
    for (int i = 0; i < num_rows; i++) {
        diagonal[i] = 0;
    }
    
}

void SparseMatrixUnknownPattern::add_off_diagonal_element(int row_index, int column_index, scalar val) {
    row[row_index].push_back({ column_index, val });
}

void SparseMatrixUnknownPattern::set_diagonal_element(int row_index, scalar val) {
    diagonal[row_index] = val;
}

void SparseMatrixUnknownPattern::calc_inverse_diagonal(std::vector<scalar> &inv_D) {
    for (int i = 0; i < num_rows; i++) {
        inv_D[i] = 1.0 / diagonal[i];
    }
}

void SparseMatrixUnknownPattern::zero() {
    for (int i = 0; i < num_rows; i++) {
        row[i].clear();
        diagonal[i] = 0;
    }
}

/* Applies this matrix to the given vector 'in', writing the result to 'result' */
void SparseMatrixUnknownPattern::apply(const std::vector<scalar> &in, std::vector<scalar> &result) {
    for (int i = 0; i < num_rows; i++) {
        result[i] = diagonal[i] * in[i];
        int row_size = row[i].size();
        for (int j = 0; j < row_size; j++) {
            result[i] += row[i][j].val * in[row[i][j].column_index];
        }
    }
}

void SparseMatrixUnknownPattern::print() {
    int num_nonzero_elements = 0;
    for (int i = 0; i < num_rows; i++) {
        printf("Row %d:  ", i);
        printf("[%d %e] ", i, diagonal[i]);
        for (vector<sparse_entry>::iterator it = row[i].begin(); it != row[i].end(); ++it) {
            printf("[%d %e] ", (*it).column_index, (*it).val);
        }
        printf("\n");
        num_nonzero_elements += row[i].size();
    }
    num_nonzero_elements += num_rows;

    printf("Num nonzero elements = %d, occupation = %f%%\n", num_nonzero_elements, (100.0 * num_nonzero_elements) / (num_rows * num_rows));
}

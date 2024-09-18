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

#ifndef SPARSEMATRIXUNKNOWNPATTERN_H_INCLUDED
#define SPARSEMATRIXUNKNOWNPATTERN_H_INCLUDED

#include <vector>

#include "mat_vec_types.h"
#include "SparseMatrixTypes.h"
#include "FFEA_return_codes.h"

using namespace std;

class SparseMatrixUnknownPattern {
public:
    SparseMatrixUnknownPattern();

    ~SparseMatrixUnknownPattern();

    int init(int num_rows, int suggested_initial_size_for_row_vectors);

    void add_off_diagonal_element(int row_index, int column_index, scalar val);

    void set_diagonal_element(int row_index, scalar val);

    void calc_inverse_diagonal(std::vector<scalar> &inv_D);

    void zero();

    /** Applies this matrix to the given vector 'in', writing the result to 'result' */
    void apply(const std::vector<scalar>& in, std::vector<scalar>& result);

    void print();

private:
    int num_rows;
    vector<vector<sparse_entry>> row;
    vector<scalar> diagonal;
};

#endif

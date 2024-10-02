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

#ifndef SPARSITYPATTERN_H_INCLUDED
#define SPARSITYPATTERN_H_INCLUDED

#include <list>
#include <vector>
#include <memory>

#include "mat_vec_types.h"
#include "FFEA_return_codes.h"

#include "SparseMatrixTypes.h"
#include "SparseMatrixFixedPattern.h"

using namespace std;

typedef struct {
    int column_index;
    vector<scalar *> source_list;

} sparse_contribution_location;

class SparsityPattern {
public:
    SparsityPattern();

    ~SparsityPattern();

    int init(int num_rows);

    /* * */
    int register_contribution(int i, int j, scalar *contrib_memory_loc);

    bool check_for_contribution(int i, int j);

    /** Factory function for making empty fixed sparsity pattern matrices from this sparsity pattern */
    std::shared_ptr<SparseMatrixFixedPattern> create_sparse_matrix();

    void print();

private:
    /** An array of vectors containing the indices of the occupied sites */
    std::vector<list<sparse_contribution_location*>> row;

    int num_nonzero_elements; ///< Total number of nonzero elements in the sparsity pattern */
};

#endif

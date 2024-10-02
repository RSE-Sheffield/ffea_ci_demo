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

#include "SparsityPattern.h"

SparsityPattern::SparsityPattern() {
    row = {};
    num_nonzero_elements = 0;
}

SparsityPattern::~SparsityPattern() {
    list<sparse_contribution_location*>::iterator it;
    for(int i = 0; i < row.size(); ++i) {
        for (it = row[i].begin(); it != row[i].end(); ++it) {
	    delete (*it);
	}	
    }
    row.clear();
    num_nonzero_elements = 0;
}

int SparsityPattern::init(int num_rows) {
    try {
        row = std::vector<list<sparse_contribution_location*>>(num_rows);
    } catch (std::bad_alloc) {
        FFEA_ERROR_MESSG("Could not allocate memory for 'row' array in SparsityPattern\n");
    }

    return FFEA_OK;
}

/* * */
int SparsityPattern::register_contribution(int i, int j, scalar *contrib_memory_loc) {
    list<sparse_contribution_location*>::iterator it;
    for (it = row[i].begin(); it != row[i].end(); ++it) {

        // If element already has sources, add the source to the list
        if ((*it)->column_index == j) {
            (*it)->source_list.push_back(contrib_memory_loc);
            return FFEA_OK;
        } 

        // If we've passed the point where we'd expect our element to be,
        // create a new contribution list and insert it
        if ((*it)->column_index > j) {
            break;
        }
    }

    num_nonzero_elements++;
    sparse_contribution_location *scl = new(std::nothrow) sparse_contribution_location();
    if (!scl) FFEA_ERROR_MESSG("Failed to allocate memory for 'scl' in SparsityPattern::register_contribution\n"); 
    scl->column_index = j;
    scl->source_list.push_back(contrib_memory_loc);
    row[i].insert(it, scl);
    
    return FFEA_OK;
}

bool SparsityPattern::check_for_contribution(int i, int j) {
    list<sparse_contribution_location*>::iterator it;
    for (it = row[i].begin(); it != row[i].end(); ++it) {
        // If element already has sources, add the source to the list
        if ((*it)->column_index == j) {
            return true;
        }
    }
    return false;
}

/* Factory function for making empty fixed sparsity pattern matrices from this sparsity pattern */
std::shared_ptr<SparseMatrixFixedPattern> SparsityPattern::create_sparse_matrix() {
    std::shared_ptr<SparseMatrixFixedPattern> sm = std::make_shared<SparseMatrixFixedPattern>();

    // Generate the array of sparse matrix elements from this sparsity pattern
    std::vector<sparse_entry> entry = std::vector<sparse_entry>(num_nonzero_elements);

    // Also generate the key (array of pointers that take us to the start of each row)
    std::vector<int> key = std::vector<int>(row.size() + 1, 0);

    // Generate the array of sources for each element in the sparse matrix
    std::vector<sparse_entry_sources> source_list = std::vector<sparse_entry_sources>(num_nonzero_elements);

    int pos = 0;
    for (int i = 0; i < row.size(); i++) {
        // Store the index of the start of each row in the key
        key[i] = pos;

        // Dump the data serially, initialising values to zero and column indices to those
        // given by the sparsity pattern
        for (auto it = row[i].begin(); it != row[i].end(); ++it) {
            entry[pos].val = 0;
            entry[pos].column_index = (*it)->column_index;

            source_list[pos].init((*it)->source_list.size());
            for (int j = 0; j < source_list[pos].get_num_sources(); j++) {
                source_list[pos].set_source(j, (*it)->source_list[j]);
            }

            ++pos;
        }
    }

    // Last index in the key array should contain the total number of nonzero elements in the matrix
    key[row.size()] = pos;

    // Initialise the sparse matrix with the array of entries and the row access key
    // Use of std::move() here moves the local scope entry/key into the method (converts them to rval)
    sm->init(row.size(), num_nonzero_elements, std::move(entry), std::move(key), source_list);

    // Return pointer to the newly allocated and initialised sparse matrix
    return sm;
}

void SparsityPattern::print() {
    for (int i = 0; i < row.size(); i++) {
        printf("= ");
        for (list<sparse_contribution_location*>::iterator it = row[i].begin(); it != row[i].end(); ++it) {
            printf("[%d %d] ", (*it)->column_index, (int) ((*it)->source_list.size()));
        }
        printf("\n");
    }
}


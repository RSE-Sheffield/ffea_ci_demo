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

#include "SparseMatrixTypes.h" 

sparse_entry_sources::sparse_entry_sources() {
    sources = {};
}

sparse_entry_sources::~sparse_entry_sources() {
    sources.clear();
}

void sparse_entry_sources::init(int num_sources) {
    try {
        sources = std::vector<scalar*>(num_sources);
    } catch (std::bad_alloc &) {
        throw FFEAException("Could not allocate memory (for sources array).");
    }
}

void sparse_entry_sources::set_source(int i, scalar *s) {
    sources[i] = s;
}

scalar sparse_entry_sources::sum_all_sources() {
    scalar sum = 0.0;
    for (int i = 0; i < sources.size(); i++) {
        sum += *(sources[i]);
    }
    return sum;
}

int sparse_entry_sources::get_num_sources() {
    return sources.size();
}


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

#ifndef SPARSEMATRIXTYPES_H_INCLUDED
#define SPARSEMATRIXTYPES_H_INCLUDED

#include <vector>

#include "mat_vec_types.h"

struct sparse_entry {
    int column_index; ///< The column index of this entry in the original matrix

    scalar val; ///< The value of this entry
};

class sparse_entry_sources {
public:

    void init(int num_sources);

    void set_source(int i, scalar *s);

    scalar sum_all_sources();

    int get_num_sources();

private:
    std::vector<scalar *> sources;
};

#endif

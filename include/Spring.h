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

#ifndef SPRING_H_INCLUDED
#define SPRING_H_INCLUDED

#include <array>
#include "mat_vec_types.h"

struct Spring {
    /* *  Variables */

    scalar k = 0; ///< Spring constant

    scalar l = 0; ///< Equilibrium length

    std::array<int, 2> blob_index; ///< Blobs connected to

    std::array<int, 2> conformation_index; ///< Conformations connected to

    std::array<int, 2> node_index; ///< Nodes connected to

    bool am_i_active = false; ///< Check if spring is active
};

#endif // SPRING_H_INCLUDED

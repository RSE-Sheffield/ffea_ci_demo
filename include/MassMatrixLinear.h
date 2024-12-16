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

#ifndef MASSMATRIXLINEAR_H_INCLUDED
#define MASSMATRIXLINEAR_H_INCLUDED

#define NUM_ELEMENTS_LOWER_TRIANGULAR_4X4 10

#define NUM_LINEAR_SHAPE_FUNCTIONS 4

#include "mat_vec_types.h"
#include "dimensions.h"
class MassMatrixLinear {
public:
    MassMatrixLinear();

    scalar * get_M_alpha_mem_loc(int i, int j);

    void build(scalar density, scalar vol);

    scalar get_M_alpha_value(int i, int j);
    void print_details();
private:
    std::array<scalar, NUM_ELEMENTS_LOWER_TRIANGULAR_4X4> M_alpha;

    void zero();
};


#endif

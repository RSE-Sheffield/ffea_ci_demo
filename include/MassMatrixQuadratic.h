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

#ifndef MASSMATRIXQUADRATIC_H_INCLUDED
#define MASSMATRIXQUADRATIC_H_INCLUDED

#define NUM_ELEMENTS_LOWER_TRIANGULAR_10X10 55

#define NUM_SHAPE_FUNCTIONS 10

#define NUM_TET_GAUSS_QUAD_POINTS 14

#include "SecondOrderFunctions.h"
#include "dimensions.h"
class MassMatrixQuadratic {
public:
    scalar *get_M_alpha_mem_loc(int i, int j);

    void build(std::array<mesh_node*, NUM_NODES_QUADRATIC_TET> &n);

    scalar get_M_alpha_value(int i, int j);
    void print_details();

private:
    std::array<scalar, NUM_ELEMENTS_LOWER_TRIANGULAR_10X10>  M_alpha = {};

    struct tetrahedron_gauss_point {
        scalar W;
        std::array<scalar, 4> eta;
    };

    void add_psi_dot_products(std::array<scalar, NUM_NODES_QUADRATIC_TET> &psi, scalar det_J, scalar weight);
};


#endif

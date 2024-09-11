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

#ifndef POISSONMATRIXQUADRATIC_H_INCLUDED
#define POISSONMATRIXQUADRATIC_H_INCLUDED

#define NUM_ELEMENTS_LOWER_TRIANGULAR_10X10 55

#define GRAD_PSI_1 0
#define GRAD_PSI_2 1
#define GRAD_PSI_3 2
#define GRAD_PSI_4 3
#define GRAD_PSI_5 4
#define GRAD_PSI_6 5
#define GRAD_PSI_7 6
#define GRAD_PSI_8 7
#define GRAD_PSI_9 8
#define GRAD_PSI_10 9

#define NUM_TET_GAUSS_QUAD_POINTS 14

#include "SecondOrderFunctions.h"

class PoissonMatrixQuadratic {
public:
    PoissonMatrixQuadratic();

    scalar * get_K_alpha_mem_loc(int i, int j);

    void build(mesh_node *n[10], scalar epsilon);

    scalar get_K_alpha_value(int i, int j);

private:
    scalar K_alpha[NUM_ELEMENTS_LOWER_TRIANGULAR_10X10];

    struct tetrahedron_gauss_point {
        scalar W;
        scalar eta[4];
    };

    void add_grad_dot_products(arr3 grad_psi[10], scalar det_J, scalar weight);

    scalar grad_dot(arr3 &grad_psi_i, arr3 &grad_psi_j);

    void zero();

    void scale(scalar factor);
};


#endif

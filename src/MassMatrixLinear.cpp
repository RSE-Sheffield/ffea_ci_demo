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

#include "MassMatrixLinear.h"

MassMatrixLinear::MassMatrixLinear()
    : M_alpha({}) {}
scalar *MassMatrixLinear::get_M_alpha_mem_loc(int i, int j) {
    // Symmetric Matrix
    if (i < j) {
        int temp = i;
        i = j;
        j = temp;
    }

    // Get the index of the element corresponding to the position i,j in our lower triangular K_alpha matrix
    const int c = (i * (i + 1)) / 2 + j;

    if (c < 0 || c > 9) {
        return nullptr;
    }

    // Return a pointer to this memory location
    return &(M_alpha[c]);
}

void MassMatrixLinear::build(scalar density, scalar vol) {
	// Assigning the mass based on the solver (0.1 if i == j, 0.05 otherwise for all i,j < 4)  
	int c=0;
	scalar mult = density * vol;
	for(int i = 0; i < 4; ++i) {
		for(int j = 0; j <= i; ++j) {
			if(i == j) {
				M_alpha[c] = 0.1 * mult;
			} else {
				M_alpha[c] = 0.05 * mult;
			}
			++c;
		}
	}
}

scalar MassMatrixLinear::get_M_alpha_value(int i, int j) {
    // Poisson matrix is symmetric, so convert any request for an upper triangular element into its
    // corresponding (equivalent) lower triangular element
    if (i < j) {
        std::swap(i, j);
    }

    // Get the index of the element corresponding to the position i,j in our lower triangular K_alpha matrix
    const int c = (i * (i + 1)) / 2 + j;

    if (c < 0 || c > 9) {
        return 0.0;
    }

    // Return actual value
    return M_alpha[c];
}

void MassMatrixLinear::print_details() {
    // Printing total mass
    int c = 0;
    scalar tot = 0.0;
    for(int i = 0; i < 4; ++i) {
	    for(int j = i; j < 4; ++j) {
	        if(i != j) {
		        tot += 2 * M_alpha[c];
	        } else {
		        tot += M_alpha[c];
	        }
	        c++;
	    }
    }
    // printf("Total Mass = %e\n", tot * mesoDimensions::mass);
}

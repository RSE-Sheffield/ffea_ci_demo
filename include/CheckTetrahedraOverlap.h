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

#ifndef TETRAHEDRAOVERLAP_H_INCLUDED
#define TETRAHEDRAOVERLAP_H_INCLUDED

#include <stddef.h>
#include "mat_vec_types.h"

class checkVars {
	
	public:

		/* checkVars();
		~checkVars(); */ 
		
		// Member variables
		std::array<arr3, 6> e_v1;            ///< vector edge-oriented 
		std::array<arr3, 6> e_v2;            ///< vectors edge-oriented
		std::array<int, 4> masks;  ///< for each face of the first tetrahedron stores the halfspace each vertex of the second tetrahedron belongs to

		std::array<arr3, 4> P_V1; ///< differences between the vertices of the second (first) tetrahedron and the vertex 0  of the first(second) tetrahedron
		std::array<arr3, 4> P_V2; ///< differences between the vertices of the second (first) tetrahedron and the vertex 0  of the first(second) tetrahedron

		std::array<arr4, 4>  Coord_1; ///< vertices coordinates in the affine space
		std::array<arr4, 4>  Coord_2; ///< vertices coordinates in the affine space
		arr3 n;	  ///< variable to store the normals
};

/** Fast Tetrahedron-Tetrahedron Overlap Algorithm, by Fabio Ganovelli, Frederico Ponchio, Claudio Rocchini. ACM 2002. */
bool tet_a_tetII(const arr3 &V1_0, const arr3 &V1_1, const arr3 &V1_2, const arr3 &V1_3,
	            const arr3 &V2_0, const arr3 &V2_1, const arr3 &V2_2, const arr3 &V2_3);

#endif 

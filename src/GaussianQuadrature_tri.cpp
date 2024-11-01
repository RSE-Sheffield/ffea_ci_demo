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

#include "GaussianQuadrature_tri.h"

/* */
scalar GaussianQuadrature_tri::integrate_point_to_face(gq_func* /*f*/, const arr3 &p, Face *face, int precision) {
    arr3 q;
    scalar result = 0;
    int j = gq_precision[precision].index;
    for (int i = 0; i < gq_precision[precision].num_points; i++) {
        face->barycentric_calc_point(gq_triangle[j].zeta_0, gq_triangle[j].zeta_1, gq_triangle[j].zeta_2, q);
        result += gq_triangle[j].W * f_3d(p, q);
        j++;
    }
    return (result * face->area);
}

/* */
scalar GaussianQuadrature_tri::integrate_face_to_face(gq_func* f, Face *f1, Face *f2, int precision) {
    arr3 q;
    scalar result = 0;
    int j = gq_precision[precision].index;
    for (int i = 0; i < gq_precision[precision].num_points; i++) {
        f1->barycentric_calc_point(gq_triangle[j].zeta_0, gq_triangle[j].zeta_1, gq_triangle[j].zeta_2, q);
        result += gq_triangle[j].W * integrate_point_to_face(f, q, f2, precision);
        j++;
    }
    return (result * f1->area);
}


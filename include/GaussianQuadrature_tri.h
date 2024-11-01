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

#ifndef GAUSSIANQUADRATURE_TRI_H_INCLUDED
#define GAUSSIANQUADRATURE_TRI_H_INCLUDED

#include "mat_vec_types.h"
#include "Face.h"

typedef struct {
    int index, num_points;
} precision_lookup;

typedef struct {
    /** Weight */
    scalar W;

    //@{
    /** Barycentric coords of point on triangle */
    scalar zeta_0, zeta_1, zeta_2;
    //@} 
} barycentric_gq_point;


 static barycentric_gq_point gq_triangle[] = {
    // 1 point - precision 1
    {1.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},

    // 3 point - precision 2 
    {1.0 / 3.0, 2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
    {1.0 / 3.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
    {1.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0},

    // 4 point - precision 3 
    {-0.562500000000000, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},

    {0.520833333333333, .6, .2, .2},
    {0.520833333333333, .2, .6, .2},
    {0.520833333333333, .2, .2, .6},

    // 6 point - precision 4 
    {0.109951743655322, 0.816847572980459, 0.091576213509771, 0.091576213509771},
    {0.109951743655322, 0.091576213509771, 0.816847572980459, 0.091576213509771},
    {0.109951743655322, 0.091576213509771, 0.091576213509771, 0.816847572980459},

    {0.223381589678011, 0.108103018168070, 0.445948490915965, 0.445948490915965},
    {0.223381589678011, 0.445948490915965, 0.108103018168070, 0.445948490915965},
    {0.223381589678011, 0.445948490915965, 0.445948490915965, 0.108103018168070},

    // 7 point - precision 5 
    {0.225000000000000, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},

    {0.125939180544827, 0.797426985353087, 0.101286507323456, 0.101286507323456},
    {0.125939180544827, 0.101286507323456, 0.797426985353087, 0.101286507323456},
    {0.125939180544827, 0.101286507323456, 0.101286507323456, 0.797426985353087},

    {0.132394152788506, 0.059715871789770, 0.470142064105115, 0.470142064105115},
    {0.132394152788506, 0.470142064105115, 0.059715871789770, 0.470142064105115},
    {0.132394152788506, 0.470142064105115, 0.470142064105115, 0.059715871789770},

    // 12 point - precision 6
    {0.050844906370207, 0.873821971016996, 0.063089014491502, 0.063089014491502},
    {0.050844906370207, 0.063089014491502, 0.873821971016996, 0.063089014491502},
    {0.050844906370207, 0.063089014491502, 0.063089014491502, 0.873821971016996},

    {0.116786275726379, 0.501426509658179, 0.249286745170910, 0.249286745170910},
    {0.116786275726379, 0.249286745170910, 0.501426509658179, 0.249286745170910},
    {0.116786275726379, 0.249286745170910, 0.249286745170910, 0.501426509658179},

    {0.082851075618374, 0.636502499121399, 0.310352451033785, 0.053145049844816},
    {0.082851075618374, 0.310352451033785, 0.053145049844816, 0.636502499121399},
    {0.082851075618374, 0.053145049844816, 0.636502499121399, 0.310352451033785},
    {0.082851075618374, 0.636502499121399, 0.053145049844816, 0.310352451033785},
    {0.082851075618374, 0.310352451033785, 0.636502499121399, 0.053145049844816},
    {0.082851075618374, 0.053145049844816, 0.310352451033785, 0.636502499121399},

    // 13 point - precision 7
    {-0.149570044467670, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},

    {0.175615257433204, 0.479308067841923, 0.260345966079038, 0.260345966079038},
    {0.175615257433204, 0.260345966079038, 0.479308067841923, 0.260345966079038},
    {0.175615257433204, 0.260345966079038, 0.260345966079038, 0.479308067841923},

    {0.053347235608839, 0.869739794195568, 0.065130102902216, 0.065130102902216},
    {0.053347235608839, 0.065130102902216, 0.869739794195568, 0.065130102902216},
    {0.053347235608839, 0.065130102902216, 0.065130102902216, 0.869739794195568},

    {0.077113760890257, 0.638444188569809, 0.312865496004875, 0.048690315425316},
    {0.077113760890257, 0.048690315425316, 0.638444188569809, 0.312865496004875},
    {0.077113760890257, 0.312865496004875, 0.048690315425316, 0.638444188569809},
    {0.077113760890257, 0.638444188569809, 0.048690315425316, 0.312865496004875},
    {0.077113760890257, 0.312865496004875, 0.638444188569809, 0.048690315425316},
    {0.077113760890257, 0.048690315425316, 0.312865496004875, 0.638444188569809}
};

 static precision_lookup gq_precision[] = {
    {0, 0}, // Ignore this line (no precision 0)
    {0, 1}, // Precision 1: Starts at index 0, 1 point
    {1, 3}, // Precision 2: Starts at index 1, 3 points
    {4, 4}, // Precision 3: Starts at index 4, 4 points
    {8, 6}, // Precision 4: Starts at index 8, 6 points
    {14, 7}, // Precision 5: Starts at index 14, 7 points
    {21, 12}, // Precision 6: Starts at index 21, 12 points
    {33, 13} // Precision 6: Starts at index 33, 13 points
};

class GaussianQuadrature_tri {
public:
    // Typedef the function prototype used by this classes methods for readability
    typedef scalar gq_func(arr3&, arr3&);

    /**
     * Virtual destructor is required for correct inheritance behaviour
     */
    virtual ~GaussianQuadrature_tri() = default;

    /** Integrates function f(p,q), for fixed point p and face coordinate q, at the given precision */
    scalar integrate_point_to_face(gq_func *f, const arr3 &p, Face *face, int precision);

    /** Integrates function f(p, q) between two faces, at the given precision */
    scalar integrate_face_to_face(gq_func *f, Face *f1, Face *f2, int precision);

protected:
    virtual scalar f_3d(const arr3 &p, const arr3 &q) = 0;
};
#endif

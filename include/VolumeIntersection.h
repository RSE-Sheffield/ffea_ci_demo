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

#ifndef VOLUME_INTERSECTION_H_INCLUDED
#define VOLUME_INTERSECTION_H_INCLUDED

#include "mat_vec_fns_II.h"

// scalar volumeIntersection(arr3 (&tetA)[4], arr3 (&tetB)[4]);
/** return the the overlapping volume between tetrahedra tetA and tetB */ 
template <typename T>
T volumeIntersection(const std::array<std::array<T,3>, 4> &tetA, const std::array<std::array<T,3>, 4> &tetB, bool calcCM, std::array<T,3> &cm);
template <typename T>
T volumeIntersectionII(const std::array<T,3> &tetA0, const std::array<T,3> &tetA1, const std::array<T,3> &tetA2, const std::array<T,3> &tetA3, const std::array<T,3> &tetB0, const std::array<T,3> &tetB1, const std::array<T,3> &tetB2, const std::array<T,3> &tetB3, bool calcCM, std::array<T,3> &cm);
/** return the the overlapping volume between tetrahedra tetA and tetB, and the area enclosing this volume */ 
template <typename T>
void volumeAndAreaIntersection(const std::array<std::array<T,3>,4> &tetA, const std::array<std::array<T,3>,4> &tetB, T &vol, T &area);

#endif

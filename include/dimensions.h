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

#ifndef DIMENSIONS_H_INCLUDED
#define DIMENSIONS_H_INCLUDED

#include "mat_vec_types.h"


namespace mesoDimensions {
    constexpr scalar length = 1.7e-10; ///< length = C atom VdW radius
    constexpr scalar Energy = 4.141945559999999e-21; ///< Energy = KbT, T=300K
    constexpr scalar mass = 1.994307387553024e-26; ///< mass = C atom mass
    constexpr scalar charge = 1.602176565e-19;  ///< charge = electron charge
    constexpr scalar area = 2.89e-20; ///< area = l**2
    constexpr scalar volume = 4.913e-30; ///< volume = l**3
    constexpr scalar force = 2.4364705882352941e-11; ///< force = E/l
    constexpr scalar time = 3.7302915560886126e-13; ///< time = l*sqrt(m/E)
    constexpr scalar pressure = 8.4305832688784826e8; ///< pressure = E/l**3
    constexpr scalar velocity = 4.5572845297447225e2; ///< velocity = sqrt(E/m)
}

/// atomic units
namespace atomicDimensions {
    constexpr scalar length = 5.2917721092e-11; /* length */
    constexpr scalar Energy = 4.35974417e-18; /* Energy */
    constexpr scalar mass = 9.10938291e-31; /* mass */
    constexpr scalar charge = 1.602176565e-19; /* charge */
    constexpr scalar area = 2.8002852055707021e-21; /* area */
    constexpr scalar volume = 1.4818471148644432e-31; /* volume */
    constexpr scalar force = 8.2387224544692213e-08; /* force */
    constexpr scalar time = 2.418884326505e-17; /* time */
    constexpr scalar pressure = 2.9421912e13; /* pressure */
    constexpr scalar velocity = 2.1876912633e6; /* velocity */
}

#endif

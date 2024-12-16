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

#include "Steric_solver.h"


/** do_volumeExclusion calculates the force (and not the energy, yet) of two tetrahedra */
void Steric_solver::do_interaction(Face *f1, Face *f2, std::vector<scalar> &blob_corr){
    if (blob_corr.empty()) {
       if (! f1->checkTetraIntersection(f2))
           return;
    } else {
       if (! f1->checkTetraIntersection(f2,
            blob_corr,f1->daddy_blob->blob_index, f2->daddy_blob->blob_index))
           return;
    } 
    do_steric_interaction(f1, f2, blob_corr); 
}

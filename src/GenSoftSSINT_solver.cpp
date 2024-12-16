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

#include "GenSoftSSINT_solver.h"


    /**Calculates GenSoftSSINT forces modified with periodic boundary correction in distance calculation*/
void GenSoftSSINT_solver::do_interaction(Face *f1, Face *f2, std::vector<scalar> &blob_corr){
    bool gensoft = true;
    bool intersection; 
    if (blob_corr.empty()) {
       intersection = f1->checkTetraIntersection(f2);
    } else {
       intersection = f1->checkTetraIntersection(f2,
            blob_corr,f1->daddy_blob->blob_index, f2->daddy_blob->blob_index);
    }

    if (intersection) {
      if (do_steric_interaction(f1, f2, blob_corr))
          gensoft = false; 
    } 

    if (gensoft) {
        do_gensoft_interaction(f1, f2, blob_corr);
    }
}

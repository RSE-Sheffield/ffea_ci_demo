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

#ifndef LJMATRIXHINCLUDED
#define LJMATRIXHINCLUDED

#include <stdio.h>
#include <cstring> 
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "dimensions.h"
#include "SimulationParams.h"

#define LJI(A,B)	((A) * num_ssint_face_types + (B))

#define SSINT_TYPE_UNDEFINED 0
#define SSINT_TYPE_STERIC 1
#define SSINT_TYPE_LJSTERIC 2
#define SSINT_TYPE_LJ 3
#define SSINT_TYPE_GENSOFT 4

using namespace std;

class LJ_pair {
public:
    LJ_pair();
    ~LJ_pair();
    scalar Emin;
    scalar Rmin;
};

class SSINT_matrix {
public:
    SSINT_matrix();
    ~SSINT_matrix(); 

    void init(string ssint_params_fname, string ssint_type, int calc_ssint, scalar ssint_cutoff);

   // void get_SSINT_params(int type1, int type2, map<string, scalar> *parmap);
    map<string, scalar> get_SSINT_params(int type1, int type2);   
    int get_num_types();

private:
    void init_ssint(string ssint_params_fname, string ssint_type, scalar ssint_cutoff);
    void init_steric();
    //LJ_pair *params;
    std::vector<map<string, scalar>> params;
    int num_ssint_face_types;
};

#endif

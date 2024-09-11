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

#ifndef BINDINGMATRIXHINCLUDED
#define BINDINGMATRIXHINCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include "FFEA_return_codes.h"
#include "mat_vec_types.h"
#include "mat_vec_fns.h"
#include "Face.h"

using namespace std;

class BindingSite{

	public:

		BindingSite();
		~BindingSite();

		// Structure
		int num_faces, site_type;

		/** List of faces making up the site */
		vector<Face*> faces;
		
		void print_to_screen();
		void set_num_faces(int num_faces);
		void set_type(int site_type);
		int get_type();
		void add_face(Face *aface);
		set<int> get_nodes();

		void calculate_centroid();
		void get_centroid(arr3 &v); 
		void calculate_area();
		scalar get_area();
		void calculate_characteristic_length();
		scalar get_characteristic_length();

		static bool sites_in_range(BindingSite a, BindingSite b);

	private:

		/** Centroid of the whole site (needs recalculating if simulation has continued) */
		arr3 centroid;

		/** Area of the whole site (needs recalculating if simulation has continued) */
		scalar area;

		/** Area of the whole site (needs recalculating if area has updated */
		scalar characteristic_length;
};

class BindingSite_matrix{

	public:
		
		BindingSite_matrix();
		~BindingSite_matrix();

		int init(string fname);
		int get_num_interaction_types();
		void print_to_screen();

		/** Number of interaction types allowed */
		int num_interaction_types;
		
		/** 2D matrix defining allowed interactions between site types */
		bool **interaction;
};

#endif

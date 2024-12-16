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

#include "BindingSite.h"

BindingSite_matrix::BindingSite_matrix() {
	num_interaction_types = 0;
}

BindingSite_matrix::~BindingSite_matrix() {
	num_interaction_types = 0;
	interaction.clear();
}

void BindingSite_matrix::init(string fname) {

	// Open file
	ifstream fin;
	fin.open(fname.c_str());
	if(fin.fail()) {
		throw FFEAException("'binding_params_fname' %s not found.", fname.c_str());
	}

	// Check if correct file
	int inter;
	string buf_string;
   getline(fin, buf_string); 
	boost::trim(buf_string);
	if(buf_string != "ffea binding site params file") {
		throw FFEAException("Expected 'ffea binding site params file', got '%s'.", buf_string.c_str());
	}

	// Get num_site_types
	fin >> buf_string >> num_interaction_types;
	if (num_interaction_types <= 0 || num_interaction_types > 10) {
		throw FFEAException("In %s, %d is an invalid value for 'num_interaction_types'. Value should be 0 < num_interaction_types <= 10.", fname.c_str(), num_interaction_types);
	}

	// Get all interactions
	try {
		interaction = std::vector<std::vector<bool>>(num_interaction_types);
	} catch (std::bad_alloc &) {
		throw FFEAException("Could not allocate 2D interaction array.");
	}
	for(int i = 0; i < num_interaction_types; ++i) {
	    try {
			interaction[i] = std::vector<bool>(num_interaction_types);
	    } catch (std::bad_alloc &) {
			throw FFEAException("Could not allocate memory for interaction[%d].", i);
	    }
		for(int j = 0; j < num_interaction_types; ++j) {
			if(fin.eof()) {
				throw FFEAException("EOF reached prematurely. For 'num_interaction types = %d', expected at %d x %d matrix of 0's and 1's.", num_interaction_types, num_interaction_types, num_interaction_types);
			}

			try {
				fin >> inter;
			} catch (...) {
				throw FFEAException("In %s, error reading interaction value at Row %d Column %d.", fname.c_str(), i + 1, j + 1);
			}

			if (inter == 0) {
				interaction[i][j] = false;
			} else if (inter == 1) {
				interaction[i][j] = true;
			} else {
				throw FFEAException("Binding Site Param Row %d Column %d must be either 0 or 1.", i + 1, j + 1);
			}
		}
	}

	// Is it symmetric ? It should be, really
	for(int i = 0; i < num_interaction_types; ++i) {
		for(int j = i; j < num_interaction_types; ++j) {
			if(interaction[i][j] != interaction[j][i]) {
				throw FFEAException("Binding Site Param Row %d Column %d must equal to Binding Site Param Row %d Column %d (symettric).", i + 1, j + 1, j + 1, i + 1);
			} 
		}
	}
	fin.close();
}

int BindingSite_matrix::get_num_interaction_types() {	
	return num_interaction_types;
}

void BindingSite_matrix::print_to_screen() {	
	int i, j;
	cout << "Binding Params Matrix:\n" << endl;
	for(i = 0; i < num_interaction_types; ++i) {
		cout << "\tType " << i;
	}
	cout << endl;
	for(i = 0; i < num_interaction_types; ++i) {
		cout << "Type " << i << "\t";
		for(j = 0; j < num_interaction_types; ++j) {
			cout << interaction[i][j] << "\t";
		}
		cout << endl;
	}
}

BindingSite::BindingSite() {
	num_faces = 0;
	site_type = -1;
	faces.clear();
	initialise(centroid);
	area = 0.0;
	characteristic_length = 0.0;
}

BindingSite::~BindingSite() {
	num_faces = 0;
	site_type = -1;
	faces.clear();
	initialise(centroid);
	area = 0.0;
	characteristic_length = 0.0;
}

void BindingSite::print_to_screen() {
	vector<Face*>::iterator it;
	cout << "Binding Site:";
	cout << "type = " << site_type << ", num_faces = " << num_faces << endl;
	cout << "Faces: ";
	for(it = faces.begin(); it != faces.end(); ++it) {
		cout << (*it)->index << " ";
	}
	cout << endl << endl;
}
void BindingSite::set_num_faces(int num_faces) {
	this->num_faces = num_faces;
}

void BindingSite::set_type(int site_type) {
	this->site_type = site_type;
}

int BindingSite::get_type() {
	return site_type;
}

void BindingSite::add_face(Face *aface) {
	faces.push_back(aface);
}

void BindingSite::get_centroid(arr3 &v) {	
   store(centroid, v);
}

void BindingSite::calculate_centroid() {	
	initialise(centroid);
	vector<Face*>::iterator it;
	for(it = faces.begin(); it != faces.end(); ++it) {
		arr3 &face_centroid = (*it)->get_centroid();
		centroid[0] += face_centroid[0];
		centroid[1] += face_centroid[1];
		centroid[2] += face_centroid[2];
	}
	resize(1.0/num_faces, centroid);
}

set<int> BindingSite::get_nodes() {	
	set<int> nodes;
	vector<Face*>::iterator it;
	for(it = faces.begin(); it != faces.end(); ++it) {
		for(int i = 0; i < 3; ++i) {
			nodes.insert((*it)->n[i]->index);
		}	
	}

	return nodes;
}

void BindingSite::calculate_area() {
	area = 0.0;
	vector<Face*>::iterator it;
	for(it = faces.begin(); it != faces.end(); ++it) {
		area += (*it)->get_area();
	}
}

scalar BindingSite::get_area() {
	return area;
}

void BindingSite::calculate_characteristic_length() {
	calculate_area();
	characteristic_length = sqrt(area);
}

scalar BindingSite::get_characteristic_length() {
	return characteristic_length;
}

bool BindingSite::sites_in_range(BindingSite a, BindingSite b) {
	// Calculate separation and see if less than sum of characteristic length scales / radii
	scalar separation;
	arr3 a_cent, b_cent;
	a.calculate_characteristic_length();
	b.calculate_characteristic_length();

	a.calculate_centroid();
	b.calculate_centroid();
	a.get_centroid(a_cent);
	b.get_centroid(b_cent);
	separation = sqrt(pow(a_cent[0] - b_cent[0], 2) + pow(a_cent[1] - b_cent[1], 2) + pow(a_cent[2] - b_cent[2], 2));
	//cout << "Separation = " << separation << ", Limiting distance = " << a.get_characteristic_length() + b.get_characteristic_length() << endl;
	//cout << "Char length a = " << a.get_characteristic_length() << ", Char length b = " << b.get_characteristic_length() << endl;
	if(separation < a.get_characteristic_length() + b.get_characteristic_length()) {
		return true;
	} else {
		return false;
	}
}

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

#include "KineticState.h"

KineticState::KineticState()
    : conformation_index(0)
    , base(-1)
    , target(-1)
    , bound(false)
    , base_site(nullptr)
    , target_site(nullptr) { }

void KineticState::init() {
	// Set conformation index
	conformation_index = 0;

	// No sites to bind
	base = -1;
	target = -1;
	bound = false;
}

void KineticState::init(int conf_index, int base, int target) {
	// Set conformation index
	conformation_index = conf_index;

	// Set the binding site types (from type 'base' to type 'target')
	this->base = base;
	this->target = target;

	if(base == -1 && target == -1) {
		bound = false;
	} else if (base == -1 && target != -1 || base != -1 && target == -1) {
		throw FFEAException("Either both site types must be -1 (not bound) or neither should be -1 (bound).");
	} else {
		bound = true;
	}
}

int KineticState::get_conformation_index() {	
	return conformation_index;
}

bool KineticState::is_bound() {
	return bound;
}

int KineticState::get_base_bsite_type() {
	return base;
}

int KineticState::get_target_bsite_type() {
	return target;
}

void KineticState::set_sites(BindingSite *base, BindingSite *target) {
	base_site = base;
	target_site = target;
}

BindingSite * KineticState::get_base_site() {	
	return base_site;
}

BindingSite * KineticState::get_target_site() {	
	return target_site;
}

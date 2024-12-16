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

#ifndef KINETICSTATES_H_INCLUDED
#define KINETICSTATES_H_INCLUDED

#define FFEA_KINETIC_STATE_UNBOUND	0
#define FFEA_KINETIC_STATE_BOUND	1

#include "FFEA_return_codes.h"
#include "BindingSite.h"

#include <stdio.h>
#include <iostream>
#include <set>
using namespace std;

class KineticState {
	
	public:
		
		KineticState();

		~KineticState() = default;

		void init();
		void init(int conf_index, int from, int to);

		int get_conformation_index();
		int get_base_bsite_type();
		int get_target_bsite_type();
		bool is_bound();
		void set_sites(BindingSite *base, BindingSite *target);
		BindingSite * get_base_site();
		BindingSite * get_target_site();

	private:

		/** @brief The active conformation in this state */
		int conformation_index;

		/** @brief The active binding site type in this blob in this state */
		int base;

		/** @brief The target binding site type on any blob in this state */
		int target;

		/** @brief Is the base site type bound in this state? */
		bool bound;

      //@{
		/** @brief Pointers to the actual binding sites (after state change has occured) */
		BindingSite *base_site, *target_site;
      //@}
};
#endif

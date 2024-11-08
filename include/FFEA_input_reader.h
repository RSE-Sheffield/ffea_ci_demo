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

#ifndef FFEA_INPUT_READER
#define FFEA_INPUT_READER

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include "FFEA_return_codes.h"
#include "mat_vec_types.h"

using namespace std;

/** Class to read FFEA files (pseudo-xml)*/
class FFEA_input_reader {
	
	public:
		FFEA_input_reader();
		~FFEA_input_reader();


		/** Get all lines from ffea, strip them of whitespace and return as a vector object */ 
		void file_to_lines(const string &script_fname, vector<string> &output);

		/** Extract any block from the current block */ 
		void extract_block(const string& block_title, int block_index, const vector<string> &input, vector<string> &output, bool mandatory=true);

		/** Get rvalue from block */ 
		void parse_tag(const string &input, std::array<string, 2> &output);

		/** Specifically return map data */ 
		void parse_map_tag(const string &input, std::array<int, 2> &map_indices, string &map_fname);

		/** Split string around delim and return as strings */ 
		int split_string(const string &input, string *output, const string &delim, size_t output_length);

		/** Split string around delim and return as strings vector. */ 
		int split_string(const string &input, vector<string> &output, const string &delim);

		/** Split string around delim and return as ints */ 
		int split_string(const string &input, int *output, const string &delim, size_t output_length);

		/** Split string around delim and return as scalars */ 
		int split_string(const string &input, scalar *output, const string &delim, size_t output_length);
        
	private:
		string buf_string;
		int copying;		
};

#endif

# 
#  This file is part of the FFEA simulation package
#  
#  Copyright (c) by the Theory and Development FFEA teams,
#  as they appear in the README.md file. 
# 
#  FFEA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  FFEA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
# 
#  To help us fund FFEA development, we humbly ask that you cite 
#  the research papers on the package.
#

set (NEIGHBOURLISTDIR "${PROJECT_BINARY_DIR}/tests/rods/unit/rod_neighbour_list_construction/")
file (COPY rod_neighbour_list_construction.ffeatest DESTINATION ${NEIGHBOURLISTDIR})
add_test(NAME rod_neighbour_list_construction COMMAND ${PROJECT_BINARY_DIR}/src/ffea rod_neighbour_list_construction.ffeatest)

file (COPY create_rods.py DESTINATION ${NEIGHBOURLISTDIR})
file (COPY collider_0.rod DESTINATION ${NEIGHBOURLISTDIR})
file (COPY collider_1.rod DESTINATION ${NEIGHBOURLISTDIR})
file (COPY collider_2.rod DESTINATION ${NEIGHBOURLISTDIR})
file (COPY collider_3.rod DESTINATION ${NEIGHBOURLISTDIR})

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

set (CONNECTIONPROPDIR "${PROJECT_BINARY_DIR}/tests/rods/integration/connection/")
file (COPY connection_propagation.ffeatest DESTINATION ${CONNECTIONPROPDIR})
file (COPY bend.rod DESTINATION ${CONNECTIONPROPDIR})
file (COPY twist.rod DESTINATION ${CONNECTIONPROPDIR})
file (COPY bend.ffea DESTINATION ${CONNECTIONPROPDIR})
file (COPY twist.ffea DESTINATION ${CONNECTIONPROPDIR})
add_subdirectory(sphere_structure)
add_test(NAME connection_propagation COMMAND ${PROJECT_BINARY_DIR}/src/ffea connection_propagation.ffeatest)
set_tests_properties(connection_propagation PROPERTIES LABELS "SLOW_TEST")


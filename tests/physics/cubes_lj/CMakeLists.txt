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

if(Python3_EXECUTABLE)
    add_test(NAME lennardJones COMMAND ${PROJECT_BINARY_DIR}/src/ffea cube_veryFine.ffea)
    set_tests_properties(lennardJones PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)

    add_test(NAME lennardJones_check COMMAND ${Python3_EXECUTABLE} checkE.py)
    set_tests_properties(lennardJones_check PROPERTIES DEPENDS lennardJones ENVIRONMENT_MODIFICATION PYTHONPATH=unset: ENVIRONMENT_MODIFICATION PYTHONHOME=unset:)
endif()

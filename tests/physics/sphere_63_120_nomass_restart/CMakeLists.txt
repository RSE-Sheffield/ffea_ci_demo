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
    add_test(NAME restart_run_I COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_63_120_nomass_10steps.ffea)
    set_tests_properties(restart_run_I PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)

    add_test(NAME restart_run_II COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_63_120_nomass_1-5steps.ffea)
    set_tests_properties(restart_run_II PROPERTIES ENVIRONMENT OMP_NUM_THREADS=1)

    add_test(NAME restart_run_III COMMAND ${PROJECT_BINARY_DIR}/src/ffea sphere_63_120_nomass_6-10steps.ffea)
    set_tests_properties(restart_run_III PROPERTIES DEPENDS restart_run_II ENVIRONMENT OMP_NUM_THREADS=1)

    add_test(NAME restart_check COMMAND ${Python3_EXECUTABLE} checkRuns.py)
    set_tests_properties(restart_check PROPERTIES DEPENDS "restart_run_I;restart_run_II;restart_run_III" ENVIRONMENT_MODIFICATION PYTHONPATH=unset: ENVIRONMENT_MODIFICATION PYTHONHOME=unset:)
endif()

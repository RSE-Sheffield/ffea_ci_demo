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
    set (SYMMETRYDIR "${PROJECT_BINARY_DIR}/tests/rods/integration/equipartition/")
    file (COPY bend.rod DESTINATION ${SYMMETRYDIR})
    file (COPY stretch.rod DESTINATION ${SYMMETRYDIR})
    file (COPY twist.rod DESTINATION ${SYMMETRYDIR})
    file (COPY equipartition_test.py DESTINATION ${SYMMETRYDIR})
    file (COPY realistic.ffea DESTINATION ${SYMMETRYDIR})
    add_test(NAME equipartition_test COMMAND ${Python3_EXECUTABLE} equipartition_test.py)
    set_tests_properties(equipartition_test PROPERTIES ENVIRONMENT_MODIFICATION PYTHONPATH=unset: ENVIRONMENT_MODIFICATION PYTHONHOME=unset:)
endif()

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
    set (CONNECTIONDIR "${PROJECT_BINARY_DIR}/tests/rods/unit/connection_orientation/")
    file (COPY connection.ffeatest DESTINATION ${CONNECTIONDIR})
    file (COPY tet_ascii.1.face DESTINATION ${CONNECTIONDIR})
    file (COPY tet_ascii.1.lj DESTINATION ${CONNECTIONDIR})
    file (COPY tet_ascii.1.mat DESTINATION ${CONNECTIONDIR})
    file (COPY tet_ascii.1.node DESTINATION ${CONNECTIONDIR})
    file (COPY tet_ascii.1.pin DESTINATION ${CONNECTIONDIR})
    file (COPY tet_ascii.1.stokes DESTINATION ${CONNECTIONDIR})
    file (COPY tet_ascii.1.surf DESTINATION ${CONNECTIONDIR})
    file (COPY tet_ascii.1.top DESTINATION ${CONNECTIONDIR})
    file (COPY tet_ascii.1.vdw DESTINATION ${CONNECTIONDIR})
    file (COPY tet_ascii.1.ffea DESTINATION ${CONNECTIONDIR})
    file (COPY bend.rod DESTINATION ${CONNECTIONDIR})
    file (COPY connection.py DESTINATION ${CONNECTIONDIR})
    add_test(NAME connection_test COMMAND ${Python3_EXECUTABLE} connection_orientation.py)
    set_tests_properties(connection_orientation_test PROPERTIES ENVIRONMENT_MODIFICATION PYTHONPATH=unset: ENVIRONMENT_MODIFICATION PYTHONHOME=unset:)
endif()

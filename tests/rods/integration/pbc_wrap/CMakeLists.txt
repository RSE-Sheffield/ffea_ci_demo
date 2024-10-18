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
    set (WRAP_PBC_DIR          pbc_wrap)
    set (WRAP_PBC_CTEST        int_rod_${WRAP_PBC_DIR} )
    set (WRAP_PBC_DIR_PATH     "${PROJECT_BINARY_DIR}/tests/rods/integration/${WRAP_PBC_DIR}/")

    file (COPY run_test.py                          DESTINATION ${WRAP_PBC_DIR_PATH})
    file (COPY create_straight_rod.py               DESTINATION ${WRAP_PBC_DIR_PATH})
    file (COPY wrap.ffea                            DESTINATION ${WRAP_PBC_DIR_PATH})
    file (COPY cleanup.sh                           DESTINATION ${WRAP_PBC_DIR_PATH})

    add_test(NAME ${WRAP_PBC_CTEST} COMMAND ${Python3_EXECUTABLE} run_test.py "${PROJECT_BINARY_DIR}/src/ffea")
    set_tests_properties(${WRAP_PBC_CTEST} PROPERTIES ENVIRONMENT_MODIFICATION PYTHONPATH=unset: ENVIRONMENT_MODIFICATION PYTHONHOME=unset:)
endif()
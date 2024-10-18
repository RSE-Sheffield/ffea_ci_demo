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
    set (COLLISION_FULL_DIR          stericInteract_fullySeparate)
    set (COLLISION_FULL_CTEST        int_rod_${COLLISION_FULL_DIR})
    set (COLLISION_FULL_DIR_PATH     "${PROJECT_BINARY_DIR}/tests/rods/integration/${COLLISION_FULL_DIR}/")

    file (COPY run_test.py                          DESTINATION ${COLLISION_FULL_DIR_PATH})
    file (COPY ffea_from_template.py                DESTINATION ${COLLISION_FULL_DIR_PATH})
    file (COPY create_straight_rod.py               DESTINATION ${COLLISION_FULL_DIR_PATH})
    file (COPY node_node_distance.py                DESTINATION ${COLLISION_FULL_DIR_PATH})
    file (COPY params.yml                           DESTINATION ${COLLISION_FULL_DIR_PATH})
    file (COPY readme.md                            DESTINATION ${COLLISION_FULL_DIR_PATH})
    file (COPY template.ffea                        DESTINATION ${COLLISION_FULL_DIR_PATH})

    add_test(NAME ${COLLISION_FULL_CTEST} COMMAND ${Python3_EXECUTABLE} run_test.py PROPERTIES ENVIRONMENT_MODIFICATION PYTHONPATH=unset: ENVIRONMENT_MODIFICATION PYTHONHOME=unset:)
    set_tests_properties(${COLLISION_FULL_CTEST} PROPERTIES ENVIRONMENT_MODIFICATION PYTHONPATH=unset: ENVIRONMENT_MODIFICATION PYTHONHOME=unset:)
endif()
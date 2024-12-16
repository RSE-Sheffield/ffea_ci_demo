# 

set (POINTLINEDIR "${PROJECT_BINARY_DIR}/tests/rods/unit/point_lies_within_rod_element/")
file (COPY point_lies_within_rod_element.ffeatest DESTINATION ${POINTLINEDIR})
add_test(NAME point_lies_within_rod_element COMMAND ${PROJECT_BINARY_DIR}/src/ffea point_lies_within_rod_element.ffeatest)

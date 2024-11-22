# This script will create a venv into the install directory
# and install the local copy of FFEATools in editable mode
# It is based on GetFFEATools.cmake

set(VENV_DIR "${CMAKE_INSTALL_PREFIX}/venv")
# Delete the venv if it already exists
if (EXISTS "${VENV_DIR}")
    file(REMOVE_RECURSE "${VENV_DIR}")
endif ()

# Locate base Python
#find_package(Python3 ${PYTHON3_EXACT_VERSION_ARG} REQUIRED COMPONENTS Interpreter)
# create venv
set(VENV_EXECUTABLE @INSTALL_Python3_EXECUTABLE@ -m venv)
if(WIN32)
    set(VENV_PIP "${VENV_DIR}\\Scripts\\pip.exe")
else()
    set(VENV_PIP ${VENV_DIR}/bin/pip)
endif()
message(STATUS "Creating venv")
execute_process(COMMAND ${VENV_EXECUTABLE} "${VENV_DIR}" WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}
    COMMAND_ECHO STDOUT
    OUTPUT_QUIET)
message(STATUS "Installing ffeatools and it's dependencies into venv")
execute_process(COMMAND ${CMAKE_COMMAND} -E env CC=gcc ${VENV_PIP} install --editable "${CMAKE_INSTALL_PREFIX}/ffeatools" -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND_ECHO STDOUT
    OUTPUT_QUIET
    ERROR_QUIET) # This puts out alot of junk to stderr
message(STATUS "Installing additional dependencies for testing")
execute_process(COMMAND ${VENV_PIP} install pandas wrap
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND_ECHO STDOUT
    OUTPUT_QUIET
    ERROR_QUIET) # This puts out alot of junk to stderr

set(FFEATOOLS_DOWNLOAD_VERSION "master")

include(FetchContent)
FetchContent_Declare(
  FFEATools
  GIT_REPOSITORY https://bitbucket.org/FFEA/ffeatools.git
  GIT_TAG ${FFEATOOLS_DOWNLOAD_VERSION}
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)
#set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
FetchContent_MakeAvailable(FFEATools)

# Function to find if python module MODULE_NAME is available, and error if it is not available.
macro(ffea_search_python_module MODULE_NAME)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import ${MODULE_NAME}; print(${MODULE_NAME}.__version__) if hasattr(${MODULE_NAME}, '__version__') else print('Unknown');"
    RESULT_VARIABLE _RESULT
    OUTPUT_VARIABLE _MODULE_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endmacro()

# FFEAtools requires Python, but doesn't have it's own CMakeLists
if(PYTHON3_EXACT_VERSION)
    set(PYTHON3_EXACT_VERSION_ARG ${PYTHON3_EXACT_VERSION} EXACT)
endif()

# Has a venv with ffeatools already been created?
set(VENV_DIR ${CMAKE_BINARY_DIR}/venv)
unset(Python3_EXECUTABLE)
unset(Python3_FOUND)
set(ENV{VIRTUAL_ENV} "${VENV_DIR}")
unset(ENV{CONDA_PREFIX})
#set(Python3_ROOT_DIR "${VENV_DIR}")
set(Python3_FIND_VIRTUALENV ONLY)
find_package(Python3 ${PYTHON3_EXACT_VERSION_ARG} QUIET COMPONENTS Interpreter)
#message(STATUS "first python ${Python3_EXECUTABLE}")
unset(ENV{VIRTUAL_ENV})
unset(Python3_FIND_VIRTUALENV)

if(NOT Python3_FOUND)
    # venv python was not found
    # locate system python
    find_package(Python3 ${PYTHON3_EXACT_VERSION_ARG} REQUIRED COMPONENTS Interpreter)
    #message(STATUS "second python ${Python3_EXECUTABLE}")
    # Check the found python has venv
    ffea_search_python_module(venv)
    if(NOT ${_RESULT} STREQUAL "0")
        message(FATAL_ERROR 
            "  Unable to find required python module 'venv'.\n"
            "  Please install this to your python environment, e.g.\n"
            "    python3 -m pip install venv\n"
            "    or\n"
            "    apt-get install python3-venv")
    endif()
    # create venv
    set(VENV_EXECUTABLE ${Python3_EXECUTABLE} -m venv)
    if(WIN32)
        set(VENV_PIP "${VENV_DIR}\\Scripts\\pip.exe")
    else()
        set(VENV_PIP ${VENV_DIR}/bin/pip)
    endif()
    message(STATUS "Creating venv for ffeatools")
    execute_process(COMMAND ${VENV_EXECUTABLE} "${VENV_DIR}" WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMAND_ECHO STDOUT
        OUTPUT_QUIET)
    message(STATUS "Installing ffeatools and it's dependencies into venv")
    execute_process(COMMAND ${CMAKE_COMMAND} -E env CC=gcc ${VENV_PIP} install --editable "${ffeatools_SOURCE_DIR}" -v
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMAND_ECHO STDOUT
        OUTPUT_QUIET
        ERROR_QUIET) # This puts out alot of junk to stderr

    # locate venv python
    unset(Python3_EXECUTABLE)
    message(STATUS "venv dir ${VENV_DIR}")
    set(ENV{VIRTUAL_ENV} "${VENV_DIR}")
    set(Python3_ROOT_DIR "${VENV_DIR}")
    set(Python3_FIND_VIRTUALENV ONLY)
    find_package(Python3 ${PYTHON3_EXACT_VERSION_ARG} REQUIRED COMPONENTS Interpreter)
    unset(ENV{VIRTUAL_ENV})
    unset(Python3_FIND_VIRTUALENV)
    #message(STATUS "third python ${Python3_EXECUTABLE}")
endif()
# check that ffeatools is available
ffea_search_python_module(ffeatools)
if(${_RESULT} STREQUAL "0")
    message(STATUS "Successfully installed ffeatools (version ${_MODULE_VERSION})\n"
        "Python3 venv: ${VENV_PIP}")
else()    
    message(FATAL_ERROR 
        "Something went wrong!\n"
        "Unable to find required python module ffeatools within\n"
        "Python3 venv: ${VENV_PIP}")
endif()

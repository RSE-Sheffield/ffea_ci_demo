# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /localhome/py12rw/Software/src/ffea/ffeatools/FFEA_analysis/Euler_characteristic_analysis

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /localhome/py12rw/Software/src/ffea/ffeatools_build

# Include any dependencies generated for this target.
include CMakeFiles/build_euler_characteristic_table.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/build_euler_characteristic_table.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/build_euler_characteristic_table.dir/flags.make

CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o: CMakeFiles/build_euler_characteristic_table.dir/flags.make
CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o: /localhome/py12rw/Software/src/ffea/ffeatools/FFEA_analysis/Euler_characteristic_analysis/build_euler_characteristic_table.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /localhome/py12rw/Software/src/ffea/ffeatools_build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o -c /localhome/py12rw/Software/src/ffea/ffeatools/FFEA_analysis/Euler_characteristic_analysis/build_euler_characteristic_table.cpp

CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /localhome/py12rw/Software/src/ffea/ffeatools/FFEA_analysis/Euler_characteristic_analysis/build_euler_characteristic_table.cpp > CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.i

CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /localhome/py12rw/Software/src/ffea/ffeatools/FFEA_analysis/Euler_characteristic_analysis/build_euler_characteristic_table.cpp -o CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.s

CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o.requires:
.PHONY : CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o.requires

CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o.provides: CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o.requires
	$(MAKE) -f CMakeFiles/build_euler_characteristic_table.dir/build.make CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o.provides.build
.PHONY : CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o.provides

CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o.provides.build: CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o

# Object files for target build_euler_characteristic_table
build_euler_characteristic_table_OBJECTS = \
"CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o"

# External object files for target build_euler_characteristic_table
build_euler_characteristic_table_EXTERNAL_OBJECTS =

build_euler_characteristic_table: CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o
build_euler_characteristic_table: CMakeFiles/build_euler_characteristic_table.dir/build.make
build_euler_characteristic_table: CMakeFiles/build_euler_characteristic_table.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable build_euler_characteristic_table"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/build_euler_characteristic_table.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/build_euler_characteristic_table.dir/build: build_euler_characteristic_table
.PHONY : CMakeFiles/build_euler_characteristic_table.dir/build

CMakeFiles/build_euler_characteristic_table.dir/requires: CMakeFiles/build_euler_characteristic_table.dir/build_euler_characteristic_table.o.requires
.PHONY : CMakeFiles/build_euler_characteristic_table.dir/requires

CMakeFiles/build_euler_characteristic_table.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/build_euler_characteristic_table.dir/cmake_clean.cmake
.PHONY : CMakeFiles/build_euler_characteristic_table.dir/clean

CMakeFiles/build_euler_characteristic_table.dir/depend:
	cd /localhome/py12rw/Software/src/ffea/ffeatools_build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /localhome/py12rw/Software/src/ffea/ffeatools/FFEA_analysis/Euler_characteristic_analysis /localhome/py12rw/Software/src/ffea/ffeatools/FFEA_analysis/Euler_characteristic_analysis /localhome/py12rw/Software/src/ffea/ffeatools_build /localhome/py12rw/Software/src/ffea/ffeatools_build /localhome/py12rw/Software/src/ffea/ffeatools_build/CMakeFiles/build_euler_characteristic_table.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/build_euler_characteristic_table.dir/depend


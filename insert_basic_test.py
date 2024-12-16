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

# This script requires Python 3 to run!
#
# It will automatically insert some boilerplate CMake stuff into the FFEA
# test directory, allowing ctest to find your unit test so that you can check
# it runs correctly.
#
# Currently, only rod unit tests will be written!
import os

# TODO: Check for Python 3
# TODO: Choose which test directory to put stuff in
# TODO: Check that FFEA source directory exists
# TODO: Automatically search for $FFEA_SRC bash environment variable
# TODO: Search src/ffea_test.cpp to check that test_name is entered correctly

test_dir = input("Enter name of CMake test directory in allcaps, e.g. TESTNAMEDIR: ")
test_name = input("Enter name of unit test in snake case, e.g. this_is_a_test: ")
prefix = input("Enter path to FFEA source code directory (with trailing /):  ")

path_to_test = f"tests/rods/unit/{test_name}"

cmake_text = """# 

set ({0} "${{PROJECT_BINARY_DIR}}/{1}/")
file (COPY {2}.ffeatest DESTINATION ${{{0}}})
add_test(NAME {2} COMMAND ${{PROJECT_BINARY_DIR}}/src/ffea {2}.ffeatest)
""".format(test_dir, path_to_test, test_name)

go = input(f"Unit test template will be written in {prefix + path_to_test}, proceed? (y/n)")
if go != "y" and go != "Y":
    quit()

try:
    os.mkdir(prefix + path_to_test)
except FileExistsError:
    print("Directory exists. Proceeding")
finally:
    os.chdir(prefix + path_to_test)

print(f"Writing {prefix}{path_to_test}/CMakeLists.txt")
with open(f"{prefix}{path_to_test}/CMakeLists.txt", "w") as f:
    f.write(cmake_text)
    
print(f"Writing {prefix}{path_to_test}/{test_name}.ffeatest")
with open(f"{prefix}{path_to_test}/{test_name}.ffeatest", "w") as f:
    f.write(test_name)
    
print(f"Appending to {prefix}{path_to_test}/../CMakeLists.txt")
with open(f"{prefix}{path_to_test}/../CMakeLists.txt", "a") as f:
    f.write(f"add_subdirectory({test_name})\n")
    
print("Done")

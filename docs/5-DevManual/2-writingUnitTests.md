Writing unit tests for FFEA using CTest
=======================================

Unit tests are small, testable parts of code, e.g. writing an array to a file, 
 computing the result of an equation, etc. An integration test is more
 complex, and tests the entire program against a set of criteria, e.g. running 
 a mesh for some time and ensuring it doesn't break or return energies well
 above kT.

Writing a test for a piece of code before or during its development, rather
 than after, will at first seem like more work than is necessary, but it can 
 save you headaches and time later on when you inevitably need to debug it.

The following is not the _only_ way to write a unit test for FFEA, but I think
 it makes sense:

 1. Write the main body of your test as a `static int` function in 
    `FFEA_SRC/src/ffea_test.cpp`. Let's call it `test_name`. All tests in
    this file return 0 for a pass or 1 for a fail.

 2. In the same file, add these lines to the `do_ffea_test()` function:
         if (buffer.str().find("test_name") != std::string::npos ){
         result = ffea_test::test_name();
         }

 3. Create a directory `test_name` and place it in 
    `FFEA_SRC/testing/path/to/test_name`, e.g. for a unit test relating to
    rods it should go in `FFEA_SRC/tests/rods/unit/`.

 4. Inside `test_name` create two files: `CMakeLists.txt` and 
    `test_name.ffeatest`. 

    `test_name.ffeatest` is a file unique to FFEA and contains only a single
    line: `test_name`. When FFEA is given an .ffeatest file as input it runs
    the in 'test mode'.

    `CMakeLists.txt` tells the compiler where to find all the files that your
    test uses, and hence lets CTest know that we want to be able to run it!
    You'll notice that other directories in `FFEA_SRC/testing` have their own
    `CMakeLists.txt`. It should contain the following lines:

         set (TESTDIR "${PROJECT_BINARY_DIR}/tests/path/to/test_name")
         file (COPY test_name.ffeatest DESTINATION ${TESTDIR})
         add_test(NAME test_name COMMAND ${PROJECT_BINARY_DIR}/src/ffea test_name.ffeatest)

    where the name of `TESTDIR` doesn't really matter, as long as it points to
    the directory containing your test.

 5. Append the line `add_subdirectory(test_name)` to `../CMakeLists.txt`; one 
    directory above `test_name`.

 6. Recompile the FFEA source code with CMake.

 7. Run your specific test:
         `cd $FFEA_BUILD`
         `ctest --verbose -R test_name`

 8. If anything went wrong and standard output doesn't give you much info,
    check `$FFEA_BUILD/Testing/Temporary/LastTest.log`.

The python script `$FFEA_SRC/install_basic_test.py` should do most of the
 above for you. Check out some of the already written FFEA tests for
 inspiration when writing more complicated tests!!!

Further reading: https://softwaretestingfundamentals.com/software-testing-levels

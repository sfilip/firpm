firpm library (double precision version)
=========================================

## Installation instructions
This code has been tested on recent Mac OS X (>= 10.14) and Linux installations. In order to compile and use it, a recent 
C++ compiler with C++11 support is necessary (recent g++ and clang compilers have been shown to work). An active internet 
connection and some external utilities and libraries must also be installed and available on your system search paths:
* CMake version 3.12 or newer
* Eigen version 3 or newer
* Google gtest framework for generating the test executables will be downloaded during CMake configuration
* (optional) doxygen version 1.8.3 or newer in order to generate the accompanying documentation

Assuming these prerequisites are taken care of and you are located at the base folder containing the source code and test files,
building the library can be done using the following commands (with the appropriate user privileges):

        mkdir build
        cd build
        cmake ..
        make all
        make install

This series of steps will install the library on your system. Usually the default location is /usr/local, but this behavior can
be changed when calling CMake by specifying an explicit install prefix:

        cmake -DCMAKE_INSTALL_PREFIX=/your/install/path/here ..


With gtest downloaded in the build directory, the *make all* command should have generated two test executables
* firpmlib_scaling_test : should contain code to generate all the example filters used in the article and tests for subsequent bugs
* firpmlib_extensive_test : contains a more extensive set of over 50 different filters, giving an iteration count comparison for uniform initialization, reference scaling and AFP initialization.

The *make test* command will launch these two test executables on your system.


The *make all* target also generates the documentation if Doxygen was found on your system when running CMake. It can also be
generated individually by running the command

        make doc

after CMake was called.


Executing these files will print out information regarding the final reference error (i.e. final delta value) obtained when
executing the Parks-McClellan exchange algorithm, along with iteration count information. Depending on your machine, these tests
can take some time to finish. For example, on a core Intel Xeon(R) E5-1620, it took several minutes to finish running
the firpmlib_extensive_test file.
An example output for one test case from firpm_scaling_test looks like this:

        [ RUN      ] firpm_scaling_test.combfir
        START Parks-McClellan with uniform initialization
        Final Delta     = 1.6066103765385452407e-07
        Iteration count = 12
        FINISH Parks-McClellan with uniform initialization
        START Parks-McClellan with reference scaling
        Final Delta     = 1.6066983887677986099e-07
        Iteration count = 3
        FINISH Parks-McClellan with reference scaling
        Iteration count reduction for final filter: 0.75
        [       OK ] firpm_scaling_test.combfir (263 ms)

It shows information when using both uniform initialization and reference scaling to design that particular filter. If uniform initialization
fails to converge to a result, then the output can look like this (test case extracted from the firpm_extensive_test file):

        [ RUN      ] firpm_extensive_test.extensive12
        START Parks-McClellan with uniform initialization
        WARNING: The exchange algorithm did not converge.
        TRIGGER: numerical instability
        POSSIBLE CAUSES: poor starting reference and/or a too small value for Nmax.
        Iteration count = NC
        FINISH Parks-McClellan with uniform initialization
        START Parks-McClellan with reference scaling
        Final delta     = 3.0071500558839879242e-13
        Iteration count = 7
        FINISH Parks-McClellan with reference scaling
        [       OK ] firpm_extensive_test.extensive12 (24 ms)

## TODOs
* add CI to the project
* ...
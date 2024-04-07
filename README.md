firpm library
==================================

## Folder organization

Prior to November 2019, there were three different folders containing
implementations with different numerical precision used for the computations:
* double (64-bit)
* long double (80-bit on x86 architectures)
* MPFR-based custom extendend precision

They have now been merged into one template version.

## Installation instructions

This code has been tested on recent Mac OS X (>= 10.14) and Linux installations.
In order to compile and use it, a recent C++ compiler with C++14 and OpenMP support is necessary (recent g++ and clang compilers have been shown to work). An active internet
connection and some external utilities and libraries must also be installed and
available on your system search paths:
* CMake version 3.16 or newer
* Eigen version 3.3 or newer
* (optional) GMP version 6.0 or newer for multiple precision support
* (optional) MPFR version 4.0 or newer for multiple precision support
* (optional) [mpreal](https://github.com/advanpix/mpreal) wrapper for MPFR
* (optional) doxygen to generate the accompanying code documentation
* Google gtest framework for generating the test executables will be downloaded during CMake configuration

Assuming these prerequisites are taken care of and you are located at the base
folder, building the library can be done using the following commands
(with the appropriate user privileges):

        mkdir build
        cd build
        cmake ..
        make all
        make install

This series of steps will install the library on your system. Usually the default
location is /usr/local, but this behavior can be changed when calling CMake by
specifying an explicit install prefix:

        cmake -DCMAKE_INSTALL_PREFIX=/your/install/path/here ..


With gtest downloaded in the build directory, the *make all* command should
have generated two test executables:
* firpmlib_scaling_test : should contain code to generate all the example filters used in [1] and tests for subsequent bugs
* firpmlib_extensive_test : contains a more extensive set of over 50 different filters, giving an iteration count comparison for uniform initialization, reference scaling and AFP initialization, when applicable

The *make test* command will launch these two test executables on your system.

The *make all* target also generates the documentation if Doxygen was found on
your system when running CMake. It can also be generated individually by running
the command

        make doc

after CMake was called.


Running the test executables will print out information regarding the final reference error
(i.e., final delta value) obtained when executing the Parks-McClellan exchange algorithm,
along with iteration count information. Depending on your machine, these tests can take
some time to finish.

An example output for one test case from firpm_scaling_test (using the long double
instantiation) of the routines looks like this:

        [ RUN      ] firpm_scaling_test/1.combfir
        START Parks-McClellan with uniform initialization
        Final Delta     = 1.6039185593030203561e-07
        Iteration count = 13
        FINISH Parks-McClellan with uniform initialization
        START Parks-McClellan with reference scaling
        Final Delta     = 1.6066986792983147255e-07
        Iteration count = 3
        FINISH Parks-McClellan with reference scaling
        START Parks-McClellan with AFP
        Final Delta     = 1.606565353078586782e-07
        Iteration count = 4
        FINISH Parks-McClellan with AFP
        Iteration count reduction for final filter  RS: 0.76923076923076916245
        Iteration count reduction for final filter AFP: 0.69230769230769229061
        [       OK ] firpm_scaling_test/1.combfir (378 ms)

## Use

Examples of how to use the library can be found in the **test** folder.

## Licensing

Up until April 7, 2024 the code was GPLv3+ licensed. It has since switched to a more permissive 3-Clause BSD license.

## References
[1] S.-I. Filip, A robust and scalable implementation of the Parks-McClellan
algorithm for designing FIR filters, ACM Trans. Math. Softw., vol. 43,
no. 1, Aug. 2016, Art. no. 7.

[2] S.-I. Filip, Robust tools for weighted Chebyshev approximation and
applications to digital filter design, Ph.D. dissertation, ENS de Lyon, France, 2016.

## Acknowledgements
* Graeme Smecher [@gsmecher]( https://github.com/gsmecher ) - setting up the template version of the code base
* [@jlwehle]( https://github.com/jlwehle ) - various bug reports and fixes

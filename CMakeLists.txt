cmake_minimum_required (VERSION 3.12)
set(PROJECT_NAME_STR firpmlib)
project(${PROJECT_NAME_STR})
# The project version number.
set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   2   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   1   CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)
set(doxy_main_page ${PROJECT_SOURCE_DIR}/README.md)

#-----------------------------------------------------------------------
# dependencies testing
#-----------------------------------------------------------------------
list (INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)

find_package(OpenMP REQUIRED)
find_package(Eigen3 3.3 REQUIRED)
find_package(GMP 6.0.0)
find_package(MPFR 4.0.0)
find_package(Doxygen)

if( MPFR_FOUND AND GMP_FOUND)
    add_definitions(-DHAVE_MPFR)
endif( MPFR_FOUND AND GMP_FOUND)

#-----------------------------------------------------------------------
# build configuration
#-----------------------------------------------------------------------
set(CMAKE_BUILD_TYPE Release)
add_definitions(-std=c++11 -Wall)

#-----------------------------------------------------------------------
# common include directories
#-----------------------------------------------------------------------
if( MPFR_FOUND AND GMP_FOUND )
    set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include ${MPFR_INCLUDES})
else()
    set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
endif( MPFR_FOUND AND GMP_FOUND )

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
if(DOXYGEN_FOUND)
    add_subdirectory(${PROJECT_SOURCE_DIR}/doc)
endif()

#-----------------------------------------------------------------------
# tests
#-----------------------------------------------------------------------
# use the solution proposed at:
# https://github.com/google/googletest/blob/master/googletest/README.md
# for downloading and compiling the googletest framework inside the 
# project build directory

configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . 
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build . 
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker settings
# on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# add googletest directory to our build. This defines the gtest
# and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build 
                 EXCLUDE_FROM_ALL)


# add the test files and generate the test executables
# they can be run by calling *make test* after the initial *make*                 
enable_testing()
set(PROJECT_TEST_SCALING ${PROJECT_NAME_STR}_scaling_test)
set(PROJECT_TEST_EXTENSIVE ${PROJECT_NAME_STR}_extensive_test)
set(CMAKE_CTEST_COMMAND ctest -V)
include_directories(${GTEST_INCLUDE_DIR} ${COMMON_INCLUDES} ${EIGEN3_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIRS}/unsupported/test/mpreal)

set(TEST_SRC_SCALING ${PROJECT_SOURCE_DIR}/test/scaling_tests.cpp)
set(TEST_SRC_EXTENSIVE ${PROJECT_SOURCE_DIR}/test/extensive_tests.cpp)

add_executable(${PROJECT_TEST_SCALING} ${TEST_SRC_SCALING})
add_executable(${PROJECT_TEST_EXTENSIVE} ${TEST_SRC_EXTENSIVE})

if( MPFR_FOUND AND GMP_FOUND )
    target_link_libraries(${PROJECT_TEST_SCALING}
        gtest
        gtest_main
        pthread
        ${GMP_LIBRARIES}
        ${MPFR_LIBRARIES}
        firpm    
    )
else()
    target_link_libraries(${PROJECT_TEST_SCALING}
        gtest
        gtest_main
        pthread
        firpm    
    )
endif()

if( MPFR_FOUND AND GMP_FOUND )
    target_link_libraries(${PROJECT_TEST_EXTENSIVE}
        gtest
        gtest_main
        pthread
        ${GMP_LIBRARIES}
        ${MPFR_LIBRARIES}
        firpm    
    )
else()
    target_link_libraries(${PROJECT_TEST_EXTENSIVE}
        gtest
        gtest_main
        pthread
        firpm    
    )
endif()

add_test(ScalingTests ${PROJECT_TEST_SCALING})
add_test(ExtensiveTests ${PROJECT_TEST_EXTENSIVE})